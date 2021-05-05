

#include <stdlib.h>
#include <iostream>
#include "classCavity.h"
#include "fileIO.cpp"
#include <complex>
#include <sstream>
#include <fstream>
#include <cstring>

#ifdef USE_OPENMP
	#include <omp.h>
#endif

using namespace std;

void Cavity::initializeZero(double dt, int number_transverse, double *transverse_points, double R_max, double boundary_guard_ratio)
{
	// Set delay
	setWidth(getPosition1() - getPosition0());
	int delay = floor( getWidth()* getRefInd()/(c0*dt));
	if (delay < 1)
	{
			cout << "Cavity::initializeZero(): WARNING TOO LARGE DT" << endl;
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
	if (electric_field_pluss == NULL)
	{
		electric_field_pluss = new Cyclic<std::complex<double> >(getNumberOfTimesteps(), cavity_transverse_points_num);
		electric_field_minus = new Cyclic<std::complex<double> >(getNumberOfTimesteps(), cavity_transverse_points_num);
	}
	
	if (temp_transverse_array_pluss == NULL)
	{
		temp_transverse_array_pluss = new std::complex<double>[cavity_transverse_points_num];
		temp_transverse_array_minus = new std::complex<double>[cavity_transverse_points_num];
	}

	electric_field_pluss->zeroAllFields();
	electric_field_minus->zeroAllFields();

	transfer_matrix_MacPol = new std::complex<double>[cavity_transverse_points_num];
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		transfer_matrix_MacPol[i] = 0.0;
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
	output_E_pluss_real = new std::ofstream[cavity_transverse_points_num];
	output_E_pluss_imag = new std::ofstream[cavity_transverse_points_num];
	output_E_minus_real = new std::ofstream[cavity_transverse_points_num];
	output_E_minus_imag = new std::ofstream[cavity_transverse_points_num];
	
	
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

	double n1  = getRefInd();
	double nim = getRefInd_im();
	std::complex<double> nc = n1 + I*nim; // k = w0(n1 + i*nim)/c0

	double left_th, right_th;
	getCosTh(&left_th, &right_th);
	
	double width = (getPosition1()-getPosition0());
	Epluss_transp_x1 = exp( I*w0*nc*right_th*width/c0);
	Epluss_transp_x0 = 1.0;
	Eminus_transp_x1 = 1.0;
	Eminus_transp_x0 = exp( I*w0*nc*left_th*width/c0);
	
	// Initialize BPM propagator
	double lambda = 2.0*Pi*(c0/getRefInd())/w0;
	
	if (getName().substr(0,8) == "BPMTRANS")
	{
		lens_propagator_pluss = new BPM(cavity_transverse_points_num, R_max);
		lens_propagator_pluss->initialize_focusing_lens_FFT_BPM(cavity_eq_opt_LENS1_focus, cavity_eq_opt_length, cavity_eq_opt_LENS2_focus, R_max, cavity_transverse_points_y, lambda, boundary_guard_ratio);

		lens_propagator_minus = new BPM(cavity_transverse_points_num, R_max);
		lens_propagator_minus->initialize_focusing_lens_FFT_BPM(cavity_eq_opt_LENS2_focus, cavity_eq_opt_length, cavity_eq_opt_LENS1_focus, R_max, cavity_transverse_points_y, lambda, boundary_guard_ratio);
	} else if (getName().substr(0,5) == "BPMHC")
	{
		lens_propagator_pluss = new BPM(cavity_transverse_points_num, R_max);
		lens_propagator_pluss->initialize_focusing_lens_FFT_BPM(cavity_eq_opt_LENS1_focus, cavity_eq_opt_length, 0.0, R_max, cavity_transverse_points_y, lambda, boundary_guard_ratio);

		lens_propagator_minus = new BPM(cavity_transverse_points_num, R_max);
		lens_propagator_minus->initialize_focusing_lens_FFT_BPM(0.0, cavity_eq_opt_length, cavity_eq_opt_LENS1_focus, R_max, cavity_transverse_points_y, lambda, boundary_guard_ratio);
		
	} else if (getName().substr(0,5) == "BPMFS")
	{
		freeSpace_propagator = new BPM(cavity_transverse_points_num, R_max);
		freeSpace_propagator->initialize_freeSpace_FFT_BPM(cavity_eq_opt_length, R_max, cavity_transverse_points_y, lambda, boundary_guard_ratio);
	}	else
	{
		freeSpace_propagator = new BPM(cavity_transverse_points_num, R_max);
		freeSpace_propagator->initialize_freeSpace_FFT_BPM(left_th*width, R_max, cavity_transverse_points_y, lambda, boundary_guard_ratio);
	} 
}

/*=====================================
 * Move fields forward in time
 * */
void Cavity::updateStorage_freeSpace()
{
	// Apply free space transformation
	#ifdef FFT_BPM_EXTRA_LAYERS
	freeSpace_propagator->compute_freeSpace_FFT_BPM(electric_field_pluss->getFirstContainer());
	freeSpace_propagator->compute_freeSpace_FFT_BPM(electric_field_minus->getFirstContainer());
	#endif

	electric_field_pluss->updateStorage();
	electric_field_minus->updateStorage();
	
	// Reset first values
	// Moved to first step in B.C. calculations, does not have to be here..
//	memset(electric_field_pluss->getFirstContainer(),0,cavity_transverse_points_num*sizeof(std::complex<double>));
//	memset(electric_field_minus->getFirstContainer(),0,cavity_transverse_points_num*sizeof(std::complex<double>));
}

void Cavity::updateStorage_freeSpace_forcedBPM()
{
#ifdef FFT_BPM_LENS
	// Apply free space transformation
	freeSpace_propagator->compute_freeSpace_FFT_BPM(electric_field_pluss->getFirstContainer());
	freeSpace_propagator->compute_freeSpace_FFT_BPM(electric_field_minus->getFirstContainer());
#endif

	electric_field_pluss->updateStorage();
	electric_field_minus->updateStorage();	
}

void Cavity::updateStorage_lens_pluss()
{
	// Appy lens transformation
#ifdef FFT_BPM_LENS
	lens_propagator_pluss->compute_lens_FFT_BPM(electric_field_pluss->getFirstContainer());
#endif
	
	electric_field_pluss->updateStorage();
}

void Cavity::updateStorage_lens_minus()
{
	// Appy lens transformation
#ifdef FFT_BPM_LENS
	lens_propagator_minus->compute_lens_FFT_BPM(electric_field_minus->getFirstContainer());
#endif
	
	electric_field_minus->updateStorage();
}

void Cavity::file_output_write_leftWall(int output_level)
{
	int delay = getDelay();
	double ni = getRefInd();
	double x1 = getPosition1();
	std::complex<double> E_pl[cavity_transverse_points_num];
	std::complex<double> E_mi[cavity_transverse_points_num];
	getEpluss_left_wall(E_pl);
	getEminus_left_wall(E_mi);

	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		std::complex<double> E_prop = E_pl[i] + E_mi[i];

		double tmp = real(E_prop);	
		output_E_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
		tmp = imag(E_prop);
		output_E_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));

		if (output_level==1)
		{
			tmp = real(E_pl[i]);	
			output_E_pluss_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
			tmp = imag(E_pl[i]);
			output_E_pluss_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));

			tmp = real(E_mi[i]);	
			output_E_minus_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
			tmp = imag(E_mi[i]);
			output_E_minus_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
		}
	}
}

void Cavity::file_output_write(int output_level)
{
	int delay = getDelay();
	double ni = getRefInd();
	double x1 = getPosition1();
	std::complex<double> E_pl[cavity_transverse_points_num];
	std::complex<double> E_mi[cavity_transverse_points_num];
	getEpluss_right_wall(E_pl);
	getEminus_right_wall(E_mi);

	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		std::complex<double> E_prop = E_pl[i] + E_mi[i];

		double tmp = real(E_prop);	
		output_E_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
		tmp = imag(E_prop);
		output_E_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));

		if (output_level==1)
		{
			tmp = real(E_pl[i]);	
			output_E_pluss_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
			tmp = imag(E_pl[i]);
			output_E_pluss_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));

			tmp = real(E_mi[i]);	
			output_E_minus_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
			tmp = imag(E_mi[i]);
			output_E_minus_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
		}
	}
}

void Cavity::interpolateEpluss(double z, std::complex<double> *tmp)
{
	double t_z = (z-getPosition0())/(c0/getRefInd());
	int ind = floor(t_z/cavity_dt);
	double z0,z1;
	std::complex<double> e0[cavity_transverse_points_num],e1[cavity_transverse_points_num];
	getEpluss(ind  ,e0);
	getEpluss(ind+1,e1);

	z0 =  getPosition0() + ((double)ind)*cavity_dt*(c0/getRefInd());
	z1 =  getPosition0() + ((double)(ind+1.0))*cavity_dt*(c0/getRefInd());
	
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = e0[i] + (z-z0)*(e1[i]-e0[i])/(z1-z0);
	}

	#ifdef FFT_BPM_EXTRA_LAYERS
	cout << "Cavity::interpolateEpluss_freeSpace() Cannot use 'THIS' BPM function here.." << endl;
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

void Cavity::interpolateEminus(double z, std::complex<double> *tmp)
{
	double t_z = (getPosition1()-z)/(c0/getRefInd());
	int ind = floor(t_z/cavity_dt);

	double z0,z1;
	std::complex<double> e0[cavity_transverse_points_num],e1[cavity_transverse_points_num];
	
	getEminus(ind  ,e0);
	getEminus(ind+1,e1);
	
	z0 =  getPosition1() - ((double)ind)*cavity_dt*(c0/getRefInd());
	z1 =  getPosition1() - ((double)(ind+1.0))*cavity_dt*(c0/getRefInd());

	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = e0[i] + (z-z0)*(e1[i]-e0[i])/(z1-z0);
	}

	#ifdef FFT_BPM_EXTRA_LAYERS
	cout << "Cavity::interpolateEminus_freeSpace() Cannot use 'THIS' BPM function here.." << endl;
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


void Cavity::interpolateEminus_x0(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_minus->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_minus, electric_field_minus->getLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] =  tmp[i] + cavity_weight_x0*(temp_transverse_array_minus[i]-tmp[i]);
	}
}

void Cavity::interpolateEpluss_x1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_pluss->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_pluss, electric_field_pluss->getLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_x1*(temp_transverse_array_pluss[i]-tmp[i]);
	}
}


void Cavity::interpolateEminus_x0_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_minus->get3rdLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_minus, electric_field_minus->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_x0*(temp_transverse_array_minus[i]-tmp[i]);
	}
}

void Cavity::interpolateEpluss_x1_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_pluss->get3rdLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	memcpy(temp_transverse_array_pluss, electric_field_pluss->get2ndLastContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + cavity_weight_x1*(temp_transverse_array_pluss[i]-tmp[i]);
	}
}



void Cavity::evaluateEprop(double z, std::complex<double> *tmp)
{
	if ((z >= getPosition0())&&(z <= getPosition1()))
	{

		cout << "Cavity::evaluateEprop() ERROR" << endl;
		cout << "this is the old version of this code, update for new propagator" << endl;
		exit(-1);
	/*
		double ni = getRefInd();
		double aw = getRefInd_im();
		
		std::complex<double> E_pl_z = interpolateEpluss(z);
		std::complex<double> E_mi_z = interpolateEminus(z);
	
		return E_pl_z*exp(I*w0*(z)*ni/c0)*exp(-(z-getPosition0())*aw) + E_mi_z*exp(-I*w0*(z)*ni/c0)*exp((z-getPosition1())*aw);
	*/	
	} else {
		cout << "evaluateEprop: Requesting positions outside of cavity.." << endl;
		cout << "z0 = " << getPosition0() << endl;
		cout << "z1 = " << getPosition1() << endl;
		cout << "z  = " << z << endl;
 		exit(-1);
	}
}


void Cavity::file_output_open(int out_count, int output_level)
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
		
		if ( output_level == 1 )
		{
			fileName.str("");
			fileName << baseName.str() << "_E_pluss_re_" << getName() << "_T" << i << ".dat";
			openAppendBinary(&(output_E_pluss_real[i]), fileName.str());
		
			fileName.str("");
			fileName << baseName.str() << "_E_pluss_im_" << getName() << "_T" << i << ".dat";
			openAppendBinary(&(output_E_pluss_imag[i]), fileName.str());
		
			fileName.str("");
			fileName << baseName.str() << "_E_minus_re_" << getName() << "_T" << i << ".dat";
			openAppendBinary(&(output_E_minus_real[i]), fileName.str());
		
			fileName.str("");
			fileName << baseName.str() << "_E_minus_im_" << getName() << "_T" << i << ".dat";
			openAppendBinary(&(output_E_minus_imag[i]), fileName.str());
		}
	}
}

void Cavity::file_output_close(int output_level)
{
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		output_E_real[i].close();
		output_E_imag[i].close();
		if ( output_level == 1 )
		{
			output_E_pluss_real[i].close();
			output_E_pluss_imag[i].close();
			output_E_minus_real[i].close();
			output_E_minus_imag[i].close();
		}
	}
}

void Cavity::file_save_variables(int save_count, int cavity_num)
{
	// Save all dynamic variables
        std::stringstream fileName;
	fileName.str("");
	fileName << "save/save_" << save_count << "_CAV_" << cavity_num << "_Ep.dat";
	electric_field_pluss->saveAllContainers(fileName.str());
	
	fileName.str("");
	fileName << "save/save_" << save_count << "_CAV_" << cavity_num << "_Em.dat";
	electric_field_minus->saveAllContainers(fileName.str());
}

void Cavity::file_load_variables(int load_count, int cavity_num)
{
	// Load all dynamic variables
        std::stringstream fileName;
	fileName.str("");
	fileName << "save/save_" << load_count << "_CAV_" << cavity_num << "_Ep.dat";
	electric_field_pluss->loadAllContainers(fileName.str());
	
	fileName.str("");
	fileName << "save/save_" << load_count << "_CAV_" << cavity_num << "_Em.dat";
	electric_field_minus->loadAllContainers(fileName.str());
}

/*
	Zero out all fields
*/
void Cavity::clearFields()
{
	electric_field_pluss->zeroAllFields();
	electric_field_minus->zeroAllFields();
}





