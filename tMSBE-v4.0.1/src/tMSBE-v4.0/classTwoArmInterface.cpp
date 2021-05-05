

#include <stdlib.h>
#include <iostream>
#include "classTwoArmInterface.h"
#include "fileIO.cpp"
#include <complex>
#include <sstream>
#include <fstream>
#include <cstring>

#ifdef USE_OPENMP
	#include <omp.h>
#endif

using namespace std;

void TwoArmInterface::initializeZero(double dt, int number_transverse, double *transverse_points, double R_max, double boundary_guard_ratio)
{
	// Set delay
	setWidth(getPosition1() - getPosition0());
		
	int delay = floor( getWidth()* getRefInd()/(c0*dt));
	if (delay < 1)
	{
			cout << "TwoArmInterface::initializeZero(): WARNING TOO LARGE DT" << endl;
			cout << "dt = " << dt/fs << " [fs]" << endl;
			cout << "Need dt <= " << getWidth()* getRefInd()/(c0*fs) << " [fs]" << endl;
			exit(-1);
	}
	
	setDelay(delay);
	
	setNumberOfTimesteps(getDelay()+2);
	cavity_dt = dt;
	// Transverse points
	cavity_transverse_points_num = number_transverse;

	if (cavity_transverse_points_num < 1)
	{
		// .. error
		cout << "Cavity():: num_transverse < 1, need at LEAST 1 point" << endl;
		exit(-1);
	}
	
	if (cavity_transverse_points_y==NULL)
	{
		cavity_transverse_points_y  = new double[cavity_transverse_points_num];
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
	//t_z = (getPosition1()-getPosition0())/(c0/getRefInd());
	int cavity_Epluss_ind = floor(t_z/cavity_dt);
	z0 = getPosition0() + ((double)cavity_Epluss_ind)*cavity_dt*(c0/getRefInd());
	z1 = getPosition0() + ((double)(cavity_Epluss_ind+1.0))*cavity_dt*(c0/getRefInd());
	cavity_weight_x1 = (getPosition1()-z0)/(z1-z0);


	double n1  = getRefInd();
	double nim = getRefInd_im();
	std::complex<double> nc = n1 + I*nim; // k = w0(n1 + i*nim)/c0

	
	double left_th, right_th;
	getCosTh(&left_th, &right_th);
	
	double width = (getPosition1()-getPosition0());
	Ef_transp_x1 = 1.0;//exp( I*w0*nc*right_th*width/c0);
	Ef_transp_x0 = 1.0;
	Eb_transp_x1 = 1.0;
	Eb_transp_x0 = 1.0;//exp( I*w0*nc*left_th*width/c0);

	//Compute transverse delay variables
	if (cavity_weights_Efp==NULL)
	{
		cavity_weights_Efp  = new double[cavity_transverse_points_num];
		cavity_weights_Efm  = new double[cavity_transverse_points_num];
		cavity_weights_Ebp  = new double[cavity_transverse_points_num];
		cavity_weights_Ebm  = new double[cavity_transverse_points_num];
		Efp_trans_delay     = new std::complex<double>[cavity_transverse_points_num];
		Efm_trans_delay     = new std::complex<double>[cavity_transverse_points_num];
		Ebp_trans_delay     = new std::complex<double>[cavity_transverse_points_num];
		Ebm_trans_delay     = new std::complex<double>[cavity_transverse_points_num];
		delay_index_Efp     = new int[cavity_transverse_points_num];
		delay_index_Efm     = new int[cavity_transverse_points_num];
		delay_index_Ebp     = new int[cavity_transverse_points_num];
		delay_index_Ebm     = new int[cavity_transverse_points_num];
		
	}
	
	if (cavity_transverse_phase_in==NULL)
	{
		cavity_transverse_phase_in  = new std::complex<double>[cavity_transverse_points_num];
		cavity_transverse_phase_out = new std::complex<double>[cavity_transverse_points_num];
	}
	
	// Initialize BPM propagator
	double lambda = 2.0*Pi*(c0/getRefInd())/w0;

	#ifdef TRANS_DELAY	
		double tan_th=tan(getAngle()); //store tangent of angle of incidence for weights scheme
		double t_x=0.0;
		double t_y=0.0;
		for(int i = 0; i < cavity_transverse_points_num; i++)
		{
			t_x=(cavity_transverse_points_y[cavity_transverse_points_num-1]-cavity_transverse_points_y[i])/(c0/getRefInd());
			t_y=t_x*tan_th;		

			delay_index_Efm[i]    = floor(t_y/cavity_dt);
			z0 =  ((double)delay_index_Efm[i])*cavity_dt;
			z1 =  ((double)(delay_index_Efm[i]+1.0))*cavity_dt;
			cavity_weights_Efm[i] = (t_y-z0)/(z1-z0);
			Efm_trans_delay[i]    = exp(I*w0*nc*t_y/c0);	
		
			delay_index_Ebm[i]    = delay_index_Efm[i];
			cavity_weights_Ebm[i] = cavity_weights_Efm[i];
			Ebm_trans_delay[i]    = Efm_trans_delay[i];
		}
		for(int i = 0; i < cavity_transverse_points_num; i++)
		{	
			//Enforce mirroring on the opposite arm	
			delay_index_Efp[i]    = delay_index_Efm[0] - delay_index_Efm[i];
			cavity_weights_Efp[i] = cavity_weights_Efm[cavity_transverse_points_num-1-i];
			Efp_trans_delay[i]    = Efm_trans_delay[cavity_transverse_points_num-1-i];
		
			delay_index_Ebp[i]    = delay_index_Efp[i];
			cavity_weights_Ebp[i] = cavity_weights_Efp[i];	
			Ebp_trans_delay[i]    = Efp_trans_delay[i];
		}

		for(int i = 0; i < cavity_transverse_points_num; i++)
		{
			cavity_transverse_phase_in[i] = exp(I*2.0*Pi*getAngle()*cavity_transverse_points_y[i]/lambda);
			cavity_transverse_phase_out[i] = exp(-I*2.0*Pi*getAngle()*cavity_transverse_points_y[i]/lambda);
		}
	#else
		for(int i = 0; i < cavity_transverse_points_num; i++)
		{

			//SET ARBITRARILY SHORT DELAY WHEN NOT USING TRANSVERSE DELAY
			delay_index_Efm[i]    = 0;
			delay_index_Ebm[i]    = 0;
			delay_index_Efp[i]    = 0;
			delay_index_Ebp[i]    = 0;
		}
		
		for(int i = 0; i < cavity_transverse_points_num; i++)
		{
			cavity_transverse_phase_in[i] = exp(I*2.0*Pi*getAngle()*cavity_transverse_points_y[i]/lambda);
			cavity_transverse_phase_out[i] = exp(-I*2.0*Pi*getAngle()*cavity_transverse_points_y[i]/lambda);
		}
		
	#endif

		if (electric_field_fp == NULL)
		{
			electric_field_fp = new Cyclic<std::complex<double> > * [cavity_transverse_points_num];
			for (int i = 0; i < cavity_transverse_points_num; i++)
			{
				electric_field_fp[i] = new Cyclic<std::complex<double> >(delay_index_Efp[i]+2, 1);
			}
		}
		if (electric_field_fm == NULL)
		{
			electric_field_fm = new Cyclic<std::complex<double> > * [cavity_transverse_points_num];
			for (int i = 0; i < cavity_transverse_points_num; i++)
			{
				electric_field_fm[i] = new Cyclic<std::complex<double> >(delay_index_Efm[i]+2, 1);
			}
		}
		if (electric_field_bp == NULL)
		{
			electric_field_bp = new Cyclic<std::complex<double> > * [cavity_transverse_points_num];
			for (int i = 0; i < cavity_transverse_points_num; i++)
			{
				electric_field_bp[i] = new Cyclic<std::complex<double> >(delay_index_Ebp[i]+2, 1);
			}
		}
		if (electric_field_bm == NULL)
		{
			electric_field_bm = new Cyclic<std::complex<double> > * [cavity_transverse_points_num];
			for (int i = 0; i < cavity_transverse_points_num; i++)
			{
				electric_field_bm[i] = new Cyclic<std::complex<double> >(delay_index_Ebm[i]+2, 1);
			}
		}

		if (temp_transverse_array == NULL)
		{
			temp_transverse_array = new std::complex<double>[cavity_transverse_points_num];
		}

		for( int i =0; i < cavity_transverse_points_num; i++)
		{
			electric_field_fp[i]->zeroAllFields();
			electric_field_fm[i]->zeroAllFields();
			electric_field_bp[i]->zeroAllFields();
			electric_field_bm[i]->zeroAllFields();
		}
	
	//freeSpace_propagator = new BPM(cavity_transverse_points_num, R_max);
	//freeSpace_propagator->initialize_freeSpace_FFT_BPM(left_th*width, R_max, cavity_transverse_points_y, lambda, boundary_guard_ratio);
}

/*=====================================
 * Move fields forward in time
 * */
void TwoArmInterface::updateStorage_freeSpace()
{
	for( int i = 0; i < cavity_transverse_points_num; i++)
	{
		electric_field_fp[i]->updateStorage();
		electric_field_fm[i]->updateStorage();	
		electric_field_bp[i]->updateStorage();
		electric_field_bm[i]->updateStorage();	
	}
}


void TwoArmInterface::file_output_write(int output_level)
{
	int delay = getDelay();
	double ni = getRefInd();
	double x1 = getPosition1();
	std::complex<double> E_fp[cavity_transverse_points_num];
	std::complex<double> E_fm[cavity_transverse_points_num];
	std::complex<double> E_bp[cavity_transverse_points_num];
	std::complex<double> E_bm[cavity_transverse_points_num];
	getBoundary_Efp(E_fp);
	getBoundary_Efm(E_fm);
	getBoundary_Ebp(E_bp);
	getBoundary_Ebm(E_bm);

	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		std::complex<double> E_prop = E_fp[i] + E_fm[i] + E_bp[i] + E_bm[i];

		double tmp = real(E_prop);	
		output_E_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
		tmp = imag(E_prop);
		output_E_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));

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

/* Temporarily kept for reference for future functions
void TwoArmInterface::interpolateEfp(double z, std::complex<double> *tmp)
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
	cout << "TwoArmInterface::interpolateEpluss_freeSpace() Cannot use 'THIS' BPM function here.." << endl;
	exit(-1);
	//freeSpace_propagator->compute_freeSpace_FFT_BPM(tmp);
	#else
	std::complex<double> transp = exp( I*w0*t_z);
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i]*transp;
	}
	#endif
}*/

void TwoArmInterface::getBoundary_Efp(std::complex<double>*tmp)
{
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = getEfp_singlePoint( i , delay_index_Efp[i], Efp_trans_delay[i], cavity_weights_Efp[i]); 
	}
}

std::complex<double>  TwoArmInterface::getEfp_singlePoint(int cav_point, int ind, std::complex<double> cav_transPhase, double cav_weight)
{
	std::complex<double> tmp;
	
	tmp = *(electric_field_fp[cav_point]->get2ndLastContainer());
	
	#ifdef FFT_BPM_EXTRA_LAYERS
		cout << "TwoArmInterface::getEbm_singlePoint() Cannot use 'THIS' BPM function here.." << endl;
		exit(-1);
	#else
		tmp = tmp*cav_transPhase;
	#endif
	return tmp;
}

void TwoArmInterface::getBoundary_Efm(std::complex<double>*tmp)
{
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = getEfm_singlePoint( i , delay_index_Efm[i], Efm_trans_delay[i], cavity_weights_Efm[i]); 
	}
}

std::complex<double>  TwoArmInterface::getEfm_singlePoint(int cav_point, int ind, std::complex<double> cav_transPhase, double cav_weight)
{
	std::complex<double> tmp;
	
	tmp = *(electric_field_fm[cav_point]->get2ndLastContainer());
	
	#ifdef FFT_BPM_EXTRA_LAYERS
		cout << "TwoArmInterface::getEbm_singlePoint() Cannot use 'THIS' BPM function here.." << endl;
		exit(-1);
	#else
		tmp = tmp*cav_transPhase;
	#endif
	return tmp;
}

void TwoArmInterface::getBoundary_Ebp(std::complex<double>*tmp)
{
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = getEbp_singlePoint( i , delay_index_Ebp[i], Ebp_trans_delay[i], cavity_weights_Ebp[i]); 
	}
}

std::complex<double>  TwoArmInterface::getEbp_singlePoint(int cav_point, int ind, std::complex<double> cav_transPhase, double cav_weight)
{
	std::complex<double> tmp;
	
	tmp = *(electric_field_bp[cav_point]->get2ndLastContainer());
	
	#ifdef FFT_BPM_EXTRA_LAYERS
		cout << "TwoArmInterface::getEbm_singlePoint() Cannot use 'THIS' BPM function here.." << endl;
		exit(-1);
	#else
		tmp = tmp*cav_transPhase;
	#endif
	return tmp;
}

void TwoArmInterface::getBoundary_Ebm(std::complex<double>*tmp)
{
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = getEbm_singlePoint( i , delay_index_Ebm[i], Ebm_trans_delay[i], cavity_weights_Ebm[i]); 
	}
}

std::complex<double>  TwoArmInterface::getEbm_singlePoint(int cav_point, int ind, std::complex<double> cav_transPhase, double cav_weight)
{
	std::complex<double> tmp;
	
	tmp = *(electric_field_bm[cav_point]->get2ndLastContainer());
	
	#ifdef FFT_BPM_EXTRA_LAYERS
		cout << "TwoArmInterface::getEbm_singlePoint() Cannot use 'THIS' BPM function here.." << endl;
		exit(-1);
	#else
		tmp = tmp*cav_transPhase;
	#endif
	return tmp;
}


/*Interpolation functions kept as a reference for future versions

void TwoArmInterface::interpolateEbp_x0(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_bp->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array, electric_field_bp->getLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] =  tmp[i] + cavity_weight_x0*(temp_transverse_array[i]-tmp[i]);
	}
}

void TwoArmInterface::interpolateEbm_x0_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_bm->get3rdLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array, electric_field_bm->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_x0*(temp_transverse_array[i]-tmp[i]);
		tmp[i] = tmp[i] * cavity_transverse_phase_out[i];
	}
}*/


void TwoArmInterface::file_output_open(int out_count, int output_level)
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

void TwoArmInterface::file_output_close(int output_level)
{
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		output_E_real[i].close();
		output_E_imag[i].close();
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

void TwoArmInterface::file_save_variables(int save_count, int cavity_num)
{
	// Save all dynamic variables
        std::stringstream fileName;

	for(int i=0; i < cavity_transverse_points_num; i++)
	{
		fileName.str("");
		fileName << "save/save_" << save_count << "_INT_" << cavity_num << "_Efp"<<i<<".dat";
		electric_field_fp[i]->saveAllContainers(fileName.str());

		fileName.str("");
		fileName << "save/save_" << save_count << "_INT_" << cavity_num << "_Efm"<<i<<".dat";
		electric_field_fm[i]->saveAllContainers(fileName.str());

		fileName.str("");
		fileName << "save/save_" << save_count << "_INT_" << cavity_num << "_Ebp"<<i<<".dat";
		electric_field_bp[i]->saveAllContainers(fileName.str());

		fileName.str("");
		fileName << "save/save_" << save_count << "_INT_" << cavity_num << "_Ebm"<<i<<".dat";
		electric_field_bm[i]->saveAllContainers(fileName.str());
	}
}

void TwoArmInterface::file_load_variables(int load_count, int cavity_num)
{
	// Load all dynamic variables
        std::stringstream fileName;
	
	for(int i=0; i < cavity_transverse_points_num; i++)
	{
		fileName.str("");
		fileName << "save/save_" << load_count << "_INT_" << cavity_num << "_Efp"<<i<<".dat";
		electric_field_fp[i]->loadAllContainers(fileName.str());

		fileName.str("");
		fileName << "save/save_" << load_count << "_INT_" << cavity_num << "_Efm"<<i<<".dat";
		electric_field_fm[i]->loadAllContainers(fileName.str());

		fileName.str("");
		fileName << "save/save_" << load_count << "_INT_" << cavity_num << "_Ebp"<<i<<".dat";
		electric_field_bp[i]->loadAllContainers(fileName.str());

		fileName.str("");
		fileName << "save/save_" << load_count << "_INT_" << cavity_num << "_Ebm"<<i<<".dat";
		electric_field_bm[i]->loadAllContainers(fileName.str());
	}
}

/*
	Zero out all fields
*/
void TwoArmInterface::clearFields()
{
	for( int i = 0; i < cavity_transverse_points_num; i++)
	{
		electric_field_fp[i]->zeroAllFields();
		electric_field_fm[i]->zeroAllFields();
		electric_field_bp[i]->zeroAllFields();
		electric_field_bm[i]->zeroAllFields();
	}
}





