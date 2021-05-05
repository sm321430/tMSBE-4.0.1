

#include <stdlib.h>
#include <iostream>
#include "classTwoArmCavity.h"
#include "fileIO.cpp"
#include <complex>
#include <sstream>
#include <fstream>
#include <cstring>

#ifdef USE_OPENMP
	#include <omp.h>
#endif

using namespace std;

void TwoArmCavity::initializeZero(double dt, int number_transverse, double *transverse_points, double R_max, double boundary_guard_ratio)
{
	// Set delay
	setWidth(getPosition1() - getPosition0());
	int delay = floor( getWidth()* getRefInd()/(c0*dt));
	if (delay < 1)
	{
			cout << "TwoArmCavity::initializeZero(): WARNING TOO LARGE DT" << endl;
			cout << "dt = " << dt/fs << " [fs]" << endl;
			cout << "Need dt <= " << getWidth()* getRefInd()/(c0*fs) << " [fs]" << endl;
			exit(-1);
	}

	setDelay(delay);
	setNumberOfTimesteps(getDelay()+2);
	
	// Set delay
	if(getRefInd_extraAxis()>0)
	{	
		delay = floor( getWidth()* getRefInd_extraAxis()/(c0*dt));
		if (delay < 1)
		{
				cout << "TwoArmCavity::initializeZero(): Extra Axis WARNING TOO LARGE DT" << endl;
				cout << "dt = " << dt/fs << " [fs]" << endl;
				cout << "Need dt <= " << getWidth()* getRefInd_extraAxis()/(c0*fs) << " [fs]" << endl;
				exit(-1);
		}

		setDelay_extraAxis(delay);
		setNumberOfTimesteps_extraAxis(getDelay_extraAxis()+2);
	}

	cavity_dt = dt;

	// Transverse points
	cavity_transverse_points_num = number_transverse;

	if (cavity_transverse_points_num < 1)
	{
		// .. error
		cout << "Cavity():: num_transverse < 1, need at LEAST 1 point" << endl;
		exit(-1);
	}
	if (electric_field_fp == NULL)
	{
		electric_field_fp = new Cyclic<std::complex<double> >(getNumberOfTimesteps(), cavity_transverse_points_num);
	}
	if (electric_field_bp == NULL)
	{
		electric_field_bp = new Cyclic<std::complex<double> >(getNumberOfTimesteps(), cavity_transverse_points_num);
	}
	if (getNumberOfTimesteps_extraAxis()>0)
	{
		if (electric_field_fm == NULL)
		{
			electric_field_fm = new Cyclic<std::complex<double> >(getNumberOfTimesteps_extraAxis(), cavity_transverse_points_num);
		}
		if (electric_field_bm == NULL)
		{
			electric_field_bm = new Cyclic<std::complex<double> >(getNumberOfTimesteps_extraAxis(), cavity_transverse_points_num);
		}
	} else
	{
		if (electric_field_fm == NULL)
		{
			electric_field_fm = new Cyclic<std::complex<double> >(getNumberOfTimesteps(), cavity_transverse_points_num);
		}
		if (electric_field_bm == NULL)
		{
			electric_field_bm = new Cyclic<std::complex<double> >(getNumberOfTimesteps(), cavity_transverse_points_num);
		}
	}

	if (temp_transverse_array_fp == NULL)
	{
		temp_transverse_array_fp = new std::complex<double>[cavity_transverse_points_num];
		temp_transverse_array_fm = new std::complex<double>[cavity_transverse_points_num];
		temp_transverse_array_bp = new std::complex<double>[cavity_transverse_points_num];
		temp_transverse_array_bm = new std::complex<double>[cavity_transverse_points_num];
	}

	electric_field_fp->zeroAllFields();
	electric_field_fm->zeroAllFields();
	electric_field_bp->zeroAllFields();
	electric_field_bm->zeroAllFields();

	transfer_matrix_MacPol_fp = new std::complex<double>[cavity_transverse_points_num];
	transfer_matrix_MacPol_fm = new std::complex<double>[cavity_transverse_points_num];
	transfer_matrix_MacPol_bp = new std::complex<double>[cavity_transverse_points_num];
	transfer_matrix_MacPol_bm = new std::complex<double>[cavity_transverse_points_num];
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		transfer_matrix_MacPol_fp[i] = 0.0;
		transfer_matrix_MacPol_fm[i] = 0.0;
		transfer_matrix_MacPol_bp[i] = 0.0;
		transfer_matrix_MacPol_bm[i] = 0.0;
	}

	if (cavity_transverse_points_y==NULL)
	{
		cavity_transverse_points_y = new double[cavity_transverse_points_num];
	}
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		cavity_transverse_points_y[i] = transverse_points[i];
	}

	
	output_E_real = new std::ofstream[cavity_transverse_points_num];
	output_E_imag = new std::ofstream[cavity_transverse_points_num];
	output_E_fp_real = new std::ofstream[cavity_transverse_points_num];
	output_E_fp_imag = new std::ofstream[cavity_transverse_points_num];
	output_E_fm_real = new std::ofstream[cavity_transverse_points_num];
	output_E_fm_imag = new std::ofstream[cavity_transverse_points_num];
	output_E_bp_real = new std::ofstream[cavity_transverse_points_num];
	output_E_bp_imag = new std::ofstream[cavity_transverse_points_num];
	output_E_bm_real = new std::ofstream[cavity_transverse_points_num];
	output_E_bm_imag = new std::ofstream[cavity_transverse_points_num];
	
	
	//=======================================================
	// Precalculate weights for evaluation on the bundaries
	//=======================================================
	double z0,z1;
	
	// E^-(X0)
	double t_z = (getPosition1()-getPosition0())/(c0/getRefInd());
	int cavity_Eminus_ind = floor(t_z/cavity_dt);
	z0 =  getPosition1() - ((double)cavity_Eminus_ind)*cavity_dt*(c0/getRefInd());
	z1 =  getPosition1() - ((double)(cavity_Eminus_ind+1.0))*cavity_dt*(c0/getRefInd());
	cavity_weight_x0 = (getPosition0()-z0)/(z1-z0);
		
	
	// E^+(X1)	
	t_z = (getPosition1()-getPosition0())/(c0/getRefInd());
	int cavity_Epluss_ind = floor(t_z/cavity_dt);
	z0 = getPosition0() + ((double)cavity_Epluss_ind)*cavity_dt*(c0/getRefInd());
	z1 = getPosition0() + ((double)(cavity_Epluss_ind+1.0))*cavity_dt*(c0/getRefInd());
	cavity_weight_x1 = (getPosition1()-z0)/(z1-z0);
	
	// E^-(X0) extra axis
	t_z = (getPosition1()-getPosition0())/(c0/getRefInd_extraAxis());
	cavity_Eminus_ind = floor(t_z/cavity_dt);
	z0 =  getPosition1() - ((double)cavity_Eminus_ind)*cavity_dt*(c0/getRefInd());
	z1 =  getPosition1() - ((double)(cavity_Eminus_ind+1.0))*cavity_dt*(c0/getRefInd());
	cavity_weight_ex_x0 = (getPosition0()-z0)/(z1-z0);
		
	
	// E^+(X1) extra axis	
	t_z = (getPosition1()-getPosition0())/(c0/getRefInd_extraAxis());
	cavity_Epluss_ind = floor(t_z/cavity_dt);
	z0 = getPosition0() + ((double)cavity_Epluss_ind)*cavity_dt*(c0/getRefInd_extraAxis());
	z1 = getPosition0() + ((double)(cavity_Epluss_ind+1.0))*cavity_dt*(c0/getRefInd_extraAxis());
	cavity_weight_ex_x1 = (getPosition1()-z0)/(z1-z0);

	double n1  = getRefInd();
	double nim = getRefInd_im();
	std::complex<double> nc = n1 + I*nim; // k = w0(n1 + i*nim)/c0

	double left_th, right_th;
	getCosTh(&left_th, &right_th);
	
	double width = (getPosition1()-getPosition0());
	Ef_transp_x1 = exp( I*w0*nc*right_th*width/c0);
	Ef_transp_x0 = 1.0;
	Eb_transp_x1 = 1.0;
	Eb_transp_x0 = exp( I*w0*nc*left_th*width/c0);

	double n_ex = getRefInd_extraAxis();
	Ef_ex_transp_x1 = exp( I*w0*n_ex*right_th*width/c0);
	Eb_ex_transp_x0 = exp( I*w0*n_ex*left_th*width/c0);

	// Initialize BPM propagator
	double lambda = 2.0*Pi*(c0/getRefInd())/w0;

	freeSpace_propagator = new BPM(cavity_transverse_points_num, R_max);
	freeSpace_propagator->initialize_freeSpace_FFT_BPM(left_th*width, R_max, cavity_transverse_points_y, lambda, boundary_guard_ratio);
	
}

/*=====================================
 * Move fields forward in time
 * */
void TwoArmCavity::updateStorage_freeSpace()
{
	// Apply free space transformation
	#ifdef FFT_BPM_EXTRA_LAYERS
	freeSpace_propagator->compute_freeSpace_FFT_BPM(electric_field_fp->getFirstContainer());
	freeSpace_propagator->compute_freeSpace_FFT_BPM(electric_field_fm->getFirstContainer());
	freeSpace_propagator->compute_freeSpace_FFT_BPM(electric_field_bp->getFirstContainer());
	freeSpace_propagator->compute_freeSpace_FFT_BPM(electric_field_bm->getFirstContainer());
	#endif

	electric_field_fp->updateStorage();
	electric_field_fm->updateStorage();	
	electric_field_bp->updateStorage();
	electric_field_bm->updateStorage();	
	
}

void TwoArmCavity::file_output_write(int output_level)
{
	int delay = getDelay();
	double ni = getRefInd();
	double x1 = getPosition1();
	std::complex<double> E_fp[cavity_transverse_points_num];
	std::complex<double> E_fm[cavity_transverse_points_num];
	std::complex<double> E_bp[cavity_transverse_points_num];
	std::complex<double> E_bm[cavity_transverse_points_num];
	getEfp_front_wall(E_fp);
	getEfm_front_wall(E_fm);
	getEbp_front_wall(E_bp);
	getEbm_front_wall(E_bm);

	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		std::complex<double> E_prop = E_fp[i] + E_fm[i] + E_bp[i] + E_bm[i];

		double tmp = real(E_prop);	
		output_E_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
		tmp = imag(E_prop);
		output_E_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));

		if (output_level == 1)
		{
			tmp = real(E_fp[i]);	
			output_E_fp_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
			tmp = imag(E_fp[i]);
			output_E_fp_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
	
			tmp = real(E_fm[i]);	
			output_E_fm_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
			tmp = imag(E_fm[i]);
			output_E_fm_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));

			tmp = real(E_bp[i]);	
			output_E_bp_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
			tmp = imag(E_bp[i]);
			output_E_bp_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));

			tmp = real(E_bm[i]);	
			output_E_bm_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
			tmp = imag(E_bm[i]);
			output_E_bm_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
		}
	}
}


void TwoArmCavity::file_output_write_backWall(int output_level)
{
	int delay = getDelay();
	double ni = getRefInd();
	double x1 = getPosition1();
	std::complex<double> E_fp[cavity_transverse_points_num];
	std::complex<double> E_fm[cavity_transverse_points_num];
	std::complex<double> E_bp[cavity_transverse_points_num];
	std::complex<double> E_bm[cavity_transverse_points_num];
	getEfp_back_wall(E_fp);
	getEfm_back_wall(E_fm);
	getEbp_back_wall(E_bp);
	getEbm_back_wall(E_bm);

	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		std::complex<double> E_prop = E_fp[i] + E_fm[i] + E_bp[i] + E_bm[i];

		double tmp = real(E_prop);	
		output_E_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
		tmp = imag(E_prop);
		output_E_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));

		if (output_level==1)
		{
			tmp = real(E_fp[i]);	
			output_E_fp_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
			tmp = imag(E_fp[i]);
			output_E_fp_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
	
			tmp = real(E_fm[i]);	
			output_E_fm_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
			tmp = imag(E_fm[i]);
			output_E_fm_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));

			tmp = real(E_bp[i]);	
			output_E_bp_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
			tmp = imag(E_bp[i]);
			output_E_bp_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));

			tmp = real(E_bm[i]);	
			output_E_bm_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
			tmp = imag(E_bm[i]);
			output_E_bm_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
		}
	}
}

void TwoArmCavity::interpolateEfp(double z, std::complex<double> *tmp)
{
	double t_z = (z-getPosition0())/(c0/getRefInd());
	int ind = floor(t_z/cavity_dt);
	double z0,z1;
	std::complex<double> e0[cavity_transverse_points_num],e1[cavity_transverse_points_num];
	getEfp(ind  ,e0);
	getEfp(ind+1,e1);

	z0 =  getPosition0() + ((double)ind)*cavity_dt*(c0/getRefInd());
	z1 =  getPosition0() + ((double)(ind+1.0))*cavity_dt*(c0/getRefInd());
	
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = e0[i] + (z-z0)*(e1[i]-e0[i])/(z1-z0);
	}

	#ifdef FFT_BPM_EXTRA_LAYERS
	cout << "TwoArmCavity::interpolateEpluss_freeSpace() Cannot use 'THIS' BPM function here.." << endl;
	exit(-1);
	//freeSpace_propagator->compute_freeSpace_FFT_BPM(tmp);
	#else
	std::complex<double> transp = exp( I*w0*t_z);
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i]*transp;
	}
	#endif
}

void TwoArmCavity::interpolateEfm(double z, std::complex<double> *tmp)
{
	double t_z = (z-getPosition0())/(c0/getRefInd());
	int ind = floor(t_z/cavity_dt);
	double z0,z1;
	std::complex<double> e0[cavity_transverse_points_num],e1[cavity_transverse_points_num];
	getEfm(ind  ,e0);
	getEfm(ind+1,e1);

	z0 =  getPosition0() + ((double)ind)*cavity_dt*(c0/getRefInd());
	z1 =  getPosition0() + ((double)(ind+1.0))*cavity_dt*(c0/getRefInd());
	
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = e0[i] + (z-z0)*(e1[i]-e0[i])/(z1-z0);
	}

	#ifdef FFT_BPM_EXTRA_LAYERS
	cout << "TwoArmCavity::interpolateEpluss_freeSpace() Cannot use 'THIS' BPM function here.." << endl;
	exit(-1);
	//freeSpace_propagator->compute_freeSpace_FFT_BPM(tmp);
	#else
	std::complex<double> transp = exp( I*w0*t_z);
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i]*transp;
	}
	#endif
}

void TwoArmCavity::interpolateEbp(double z, std::complex<double> *tmp)
{
	double t_z = (getPosition1()-z)/(c0/getRefInd());
	int ind = floor(t_z/cavity_dt);

	double z0,z1;
	std::complex<double> e0[cavity_transverse_points_num],e1[cavity_transverse_points_num];
	
	getEbp(ind  ,e0);
	getEbp(ind+1,e1);
	
	z0 =  getPosition1() - ((double)ind)*cavity_dt*(c0/getRefInd());
	z1 =  getPosition1() - ((double)(ind+1.0))*cavity_dt*(c0/getRefInd());

	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = e0[i] + (z-z0)*(e1[i]-e0[i])/(z1-z0);
	}

	#ifdef FFT_BPM_EXTRA_LAYERS
	cout << "TwoArmCavity::interpolateEminus_freeSpace() Cannot use 'THIS' BPM function here.." << endl;
	exit(-1);
	//freeSpace_propagator->compute_freeSpace_FFT_BPM(tmp);
	#else
	std::complex<double> transp = exp( I*w0*t_z);
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i]*transp;
	}
	#endif
}

void TwoArmCavity::interpolateEbm(double z, std::complex<double> *tmp)
{
	double t_z = (getPosition1()-z)/(c0/getRefInd());
	int ind = floor(t_z/cavity_dt);

	double z0,z1;
	std::complex<double> e0[cavity_transverse_points_num],e1[cavity_transverse_points_num];
	
	getEbm(ind  ,e0);
	getEbm(ind+1,e1);
	
	z0 =  getPosition1() - ((double)ind)*cavity_dt*(c0/getRefInd());
	z1 =  getPosition1() - ((double)(ind+1.0))*cavity_dt*(c0/getRefInd());

	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = e0[i] + (z-z0)*(e1[i]-e0[i])/(z1-z0);
	}

	#ifdef FFT_BPM_EXTRA_LAYERS
	cout << "TwoArmCavity::interpolateEminus_freeSpace() Cannot use 'THIS' BPM function here.." << endl;
	exit(-1);
	//freeSpace_propagator->compute_freeSpace_FFT_BPM(tmp);
	#else
	std::complex<double> transp = exp( I*w0*t_z);
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i]*transp;
	}
	#endif
}

void TwoArmCavity::interpolateEbp_x0(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_bp->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_bp, electric_field_bp->getLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] =  tmp[i] + cavity_weight_x0*(temp_transverse_array_bp[i]-tmp[i]);
	}
}

void TwoArmCavity::interpolateEbm_x0(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_bm->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_bm, electric_field_bm->getLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] =  tmp[i] + cavity_weight_x0*(temp_transverse_array_bm[i]-tmp[i]);
	}
}

void TwoArmCavity::interpolateEbm_ex_x0(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_bm->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_bm, electric_field_bm->getLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] =  tmp[i] + cavity_weight_ex_x0*(temp_transverse_array_bm[i]-tmp[i]);
	}
}

void TwoArmCavity::interpolateEfp_x1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_fp->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_fp, electric_field_fp->getLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_x1*(temp_transverse_array_fp[i]-tmp[i]);
	}
}

void TwoArmCavity::interpolateEfm_x1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_fm->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_fm, electric_field_fm->getLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_x1*(temp_transverse_array_fm[i]-tmp[i]);
	}
}

void TwoArmCavity::interpolateEfm_ex_x1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_fm->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_fm, electric_field_fm->getLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_ex_x1*(temp_transverse_array_fm[i]-tmp[i]);
	}
}


void TwoArmCavity::interpolateEbp_x0_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_bp->get3rdLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_bp, electric_field_bp->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_x0*(temp_transverse_array_bp[i]-tmp[i]);
	}
}

void TwoArmCavity::interpolateEbm_x0_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_bm->get3rdLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_bm, electric_field_bm->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_x0*(temp_transverse_array_bm[i]-tmp[i]);
	}
}

void TwoArmCavity::interpolateEbm_ex_x0_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_bm->get3rdLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_bm, electric_field_bm->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_ex_x0*(temp_transverse_array_bm[i]-tmp[i]);
	}
}

void TwoArmCavity::interpolateEfp_x1_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_fp->get3rdLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_fp, electric_field_fp->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_x1*(temp_transverse_array_fp[i]-tmp[i]);
	}
}

void TwoArmCavity::interpolateEfm_x1_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_fm->get3rdLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_fm, electric_field_fm->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_x1*(temp_transverse_array_fm[i]-tmp[i]);
	}
}

void TwoArmCavity::interpolateEfm_ex_x1_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_fm->get3rdLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_fm, electric_field_fm->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_ex_x1*(temp_transverse_array_fm[i]-tmp[i]);
	}
}

void TwoArmCavity::file_output_open(int out_count, int output_level)
{
	// Setup output to files
	std::stringstream baseName;
	//baseName.str("");
	baseName << getToFileOutputKey() << out_count;
	
	std::stringstream fileName;
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		fileName.str("");
		fileName << baseName.str() << "_E_re_" << getName() << "_T" << i << ".dat";
		openAppendBinary(&(output_E_real[i]), fileName.str());
		
		fileName.str("");
		fileName << baseName.str() << "_E_im_" << getName() << "_T" << i << ".dat";
		openAppendBinary(&(output_E_imag[i]), fileName.str());

		if ( output_level == 1)
		{
			fileName.str("");
			fileName << baseName.str() << "_E_fp_re_" << getName() << "_T" << i << ".dat";
			openAppendBinary(&(output_E_fp_real[i]), fileName.str());
		
			fileName.str("");
			fileName << baseName.str() << "_E_fp_im_" << getName() << "_T" << i << ".dat";
			openAppendBinary(&(output_E_fp_imag[i]), fileName.str());

			fileName.str("");
			fileName << baseName.str() << "_E_fm_re_" << getName() << "_T" << i << ".dat";
			openAppendBinary(&(output_E_fm_real[i]), fileName.str());
		
			fileName.str("");
			fileName << baseName.str() << "_E_fm_im_" << getName() << "_T" << i << ".dat";
			openAppendBinary(&(output_E_fm_imag[i]), fileName.str());

			fileName.str("");
			fileName << baseName.str() << "_E_bp_re_" << getName() << "_T" << i << ".dat";
			openAppendBinary(&(output_E_bp_real[i]), fileName.str());
		
			fileName.str("");
			fileName << baseName.str() << "_E_bp_im_" << getName() << "_T" << i << ".dat";
			openAppendBinary(&(output_E_bp_imag[i]), fileName.str());
	

			fileName.str("");
			fileName << baseName.str() << "_E_bm_re_" << getName() << "_T" << i << ".dat";
			openAppendBinary(&(output_E_bm_real[i]), fileName.str());
		
			fileName.str("");
			fileName << baseName.str() << "_E_bm_im_" << getName() << "_T" << i << ".dat";
			openAppendBinary(&(output_E_bm_imag[i]), fileName.str());
		}
	}
}

void TwoArmCavity::file_output_close(int output_level)
{
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		output_E_real[i].close();
		output_E_imag[i].close();
		
		if ( output_level == 1)
		{
			output_E_fp_real[i].close();
			output_E_fp_imag[i].close();
			output_E_fm_real[i].close();
			output_E_fm_imag[i].close();
			output_E_bp_real[i].close();
			output_E_bp_imag[i].close();
			output_E_bm_real[i].close();
			output_E_bm_imag[i].close();
		}
	}
}

void TwoArmCavity::file_save_variables(int save_count, int cavity_num)
{
	// Save all dynamic variables
        std::stringstream fileName;
	fileName.str("");
	fileName << "save/save_" << save_count << "_TACAV_" << cavity_num << "_Efp.dat";
	electric_field_fp->saveAllContainers(fileName.str());

	fileName.str("");
	fileName << "save/save_" << save_count << "_TACAV_" << cavity_num << "_Efm.dat";
	electric_field_fm->saveAllContainers(fileName.str());

	fileName.str("");
	fileName << "save/save_" << save_count << "_TACAV_" << cavity_num << "_Ebp.dat";
	electric_field_bp->saveAllContainers(fileName.str());

	fileName.str("");
	fileName << "save/save_" << save_count << "_TACAV_" << cavity_num << "_Ebm.dat";
	electric_field_bm->saveAllContainers(fileName.str());
	
}

void TwoArmCavity::file_load_variables(int load_count, int cavity_num)
{
	// Load all dynamic variables
        std::stringstream fileName;
	fileName.str("");
	fileName << "save/save_" << load_count << "_TACAV_" << cavity_num << "_Efp.dat";
	electric_field_fp->loadAllContainers(fileName.str());

	fileName.str("");
	fileName << "save/save_" << load_count << "_TACAV_" << cavity_num << "_Efm.dat";
	electric_field_fm->loadAllContainers(fileName.str());

	fileName.str("");
	fileName << "save/save_" << load_count << "_TACAV_" << cavity_num << "_Ebp.dat";
	electric_field_bp->loadAllContainers(fileName.str());

	fileName.str("");
	fileName << "save/save_" << load_count << "_TACAV_" << cavity_num << "_Ebm.dat";
	electric_field_bm->loadAllContainers(fileName.str());

}

/*
	Zero out all fields
*/
void TwoArmCavity::clearFields()
{
	electric_field_fp->zeroAllFields();
	electric_field_fm->zeroAllFields();
	electric_field_bp->zeroAllFields();
	electric_field_bm->zeroAllFields();
}





