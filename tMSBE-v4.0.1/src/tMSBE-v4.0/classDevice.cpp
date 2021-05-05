
/*
	The DEVICE class is intended to be used on types of QW or ABSORBER.
	These objects contain the variable relevant for SBE solvers and such
	Isak Kilen @ 2016
*/
#include <unistd.h>
#include "classDevice.h"
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <complex>
#include <cstring>
#include <cmath>

#include <sys/time.h> // for getTimeOfDay()
#include <limits> // For limit of double
	
#include "fileIO.cpp"

#ifdef USE_OPENMP
	#include <omp.h>
#endif

using namespace std;

/* Print diagnostic information about device
 * */
void Device::Print() const
{
	cout << "Print device:" << endl;
	cout << " -> name  = " << getName() << endl;
	cout << " -> pos_z = " << getPosition() << endl;
	cout << " -> pos_y = " << device_transverse_position << endl;
	cout << " -> dens  = " << getDensity() << endl;
	cout << " -> Eg    = " << getBandgap() << endl;
	cout << " -> pump_scale = " << device_transverse_pump_scale << endl;
	cout << " -> temp_scale = " << device_transverse_temp_scale << endl;
}



/* Construction function
 * When an empty device is created all the new variables are initialized this way
 * */
Device::Device()
{
	resetAllVariables();
}

/* Set all variables to the initial state
 * Pointers to NULL
 * Values to 0 or -1, or something arbitraty.
 * */
void Device::resetAllVariables()
{
//	cout << "Creating empty Device" << endl;
	setName("DEV");
	setToFileOutputKey("out_");
	setPosition(0.0);
	setDensity(0.0);
	setBandgap(0.0);
	setMacroscopicPolarization(0.0);
	setMacroscopicPolarization_m1(0.0);
	setMacroscopicPolarization_m2(0.0);
	setElectricField(0.0);
	computationalDelay = 0.0;
	
	device_focus_E = 1.0;
	device_effective_qw = 1.0;
	device_deph_scale = 1.0;

	device_me = 0;
	device_mh = 0;
	device_mr = 0;
	device_Eb = 0;
	device_dipolemoment = 0;
	device_dcv_hbar = 0;
	device_chemical_potential_e = 0;
	device_chemical_potential_h = 0;
	device_length_qw = 0;
	device_deph_time = 0;
	device_occ_relax_time = 0;
	device_occ_hole_relax_time = 0;
	device_pol_frequency = 0;	
	
	device_background_temp_e = 0;
	device_background_temp_h = 0;
	device_inst_temp_e = 0;
	device_inst_temp_h = 0;
	device_inst_chem_pot_e = 0;
	device_inst_chem_pot_h = 0;
	
	device_inst_temp_e_prev = 0;
	device_inst_temp_h_prev = 0;
	carrier_scattering_rate_approximation_e = NULL;
	carrier_scattering_rate_approximation_h = NULL;
	
	fe_k_inst = NULL;
	fh_k_inst = NULL;
	
	K_max = -1;
	number_K_points = 0;
	K = NULL;
	dK = NULL;
	KdK = NULL;
	
	p_k = NULL;
	ne_k = NULL;
	nh_k = NULL;
	
	fe_k = NULL;
	fh_k = NULL;
	
	Ek_e = NULL;
	Ek_h = NULL;
	sumE_k = NULL;
	
	number_Qsc_points = -1;
	number_Eta_points = -1;
	Q_sc	= NULL;
	dQ_sc	= NULL;
	Eta		= NULL;
	dEta	= NULL;

	number_Th_points = 0;
	theta_2B = NULL;
	dtheta_2B = NULL;
	number_Z_points = 0;
	number_Q_points = 0;
	number_Kp_points = 0;
	Q = NULL;
	dQ = NULL;
	Kp = NULL;
	dKp = NULL;
	wavefunction_e = NULL;
	wavefunction_e_z = NULL;
	wavefunction_h = NULL;
	wavefunction_h_z = NULL;
	coulomb_potential_confinement_ee = NULL;
	coulomb_potential_confinement_hh = NULL;
	coulomb_potential_confinement_eh = NULL;
	coulomb_potential_epsilon_inv = NULL;
	coulomb_potential_normalized_ee = NULL;
	coulomb_potential_normalized_hh = NULL;
	coulomb_potential_normalized_eh = NULL;
	carrier_scattering_rates_e_in = NULL;
	carrier_scattering_rates_e_out = NULL;
	carrier_scattering_rates_h_in = NULL;
	carrier_scattering_rates_h_out = NULL;
	carrier_scattering_rates_e_total = NULL;
	carrier_scattering_rates_h_total = NULL;
	carrier_scattering_index_ee = NULL;
	carrier_scattering_index_hh = NULL;
	carrier_scattering_index_eh = NULL;
	carrier_scattering_index_he = NULL;
	coulomb_potential_epsilon_index = NULL;

	number_Qph_points = -1;
	number_Phi_points = -1;
	Q_ph = NULL;
	dQ_ph = NULL;
	Phi = NULL;
	dPhi = NULL;
	phonon_form_factor_ee = NULL;
	phonon_form_factor_hh = NULL;
	phonon_form_factor_Q = NULL;
	carrier_scattering_phonon_eps_zero = -1.0;
	carrier_scattering_phonon_eps_inf = -1.0;
	carrier_scattering_phonon_hbar_wLO = -1;
	carrier_scattering_phonon_temp = 0;
	carrier_scattering_phonon_n = 0;
	carrier_scattering_phonon_index_e1a = NULL;
	carrier_scattering_phonon_index_e1b = NULL;
	carrier_scattering_phonon_index_e2a = NULL;
	carrier_scattering_phonon_index_e2b = NULL;
	carrier_scattering_phonon_index_h1a = NULL;
	carrier_scattering_phonon_index_h1b = NULL;
	carrier_scattering_phonon_index_h2a = NULL;
	carrier_scattering_phonon_index_h2b = NULL;
	carrier_scattering_rates_phonon_e_in = NULL;
	carrier_scattering_rates_phonon_e_out = NULL;
	carrier_scattering_rates_phonon_h_in = NULL;
	carrier_scattering_rates_phonon_h_out = NULL;
	carrier_scattering_rates_phonon_e_total = NULL;
	carrier_scattering_rates_phonon_h_total = NULL;

	CARRIER_DEVIATION_TOL = 0.0;
	carrier_scattering_prev_ne = NULL;
	carrier_scattering_prev_nh = NULL;
	renormalize_prev_ne = NULL;
	renormalize_prev_nh = NULL;

	carrier_scattering_table_e_in = NULL;
	carrier_scattering_table_e_out = NULL;
	carrier_scattering_table_h_in = NULL;
	carrier_scattering_table_h_out = NULL;
	carrier_scattering_table_gridpoints = NULL;
	carrier_scattering_table_gridpoints_numbers = 0;
	
	
	p_k_tmp = NULL;
	ne_k_tmp = NULL;
	nh_k_tmp = NULL;
	
	p_k_quad=NULL;
	ne_k_quad=NULL;
	nh_k_quad=NULL;
	
	ab4_precompute_counter = 0;
	ab4_prev_state = false;
	p_k_ab4 = NULL;
	ne_k_ab4 = NULL;
	nh_k_ab4 = NULL;
	
	coulomb_matrix_ee = NULL;
	coulomb_matrix_hh = NULL;
	coulomb_matrix_eh = NULL;
	coulomb_matrix_red_size = 0;
	coulomb_matrix_ee_red = NULL;
	coulomb_matrix_hh_red = NULL;
	coulomb_matrix_eh_red = NULL;
	coulomb_matrix_red_index_i = -1;
	coulomb_matrix_red_index_j = -1;
	coulomb_matrix_red_index = NULL;
	electric_field_renormalized = NULL;
	energy_renormalized_e = NULL;
	energy_renormalized_h = NULL;

	SBE_TIME_SCALE = 0;
	SBE_POL_DEPH = 0;
	SBE_OCC_PUMP = 0;
	SBE_OCC_HOLE = 0;
	
	hyst_tnext = 0;
	hyst_count = 0;
	hyst_dens_initial = 0;
	
	device_transverse_position = 0.0;
	device_transverse_pump_scale = 0.0;
	device_transverse_temp_scale = 1.0;

	#ifdef MPI_BALANCE_WORKLOAD
	MPI_load = NULL;
	#endif
	
	#ifdef SBE_USE_RESONANT_PUMP_MODEL
	device_pump_N0 = 5.0e14;				// Initial pump density
	device_pump_W0 = 2.0*Pi*c0/(808*nm);	// Central frequency
	device_pump_E0	= 0.1e7;				// Field contributing to Carriers
	device_pump_ETA	= 4.0e12;				// Pump width
	device_pump_wk = NULL;
	#endif
	
	#ifdef SBE_USE_SPONTAN_EMISSIONS
	device_spont_emission_wk = NULL;
	device_spont_emission_on = true; // Spont Emissions are on by default	
	#endif
}


/* Initialize all arrays and variables to their default value
 * */
void Device::sbe_initialize(double DT)
{
	std::stringstream fileName;
	std::ofstream initOut;

	// Load density, me, mh, number_k_points, K_max, from file	
	file_readInConfigFromFile();
	// Statistics
	#ifdef USE_DEVICE_TIMERS
	MainStat->newSubTimer(getName(),"iterate_rk4");
//	MainStat->newSubTimer(getName(),"c-c scatt");
//	MainStat->newSubTimer(getName(),"renormalize");
	#endif
	
	setRelativeMass();
	setExitonBindEnergy();
	device_deph_time = (1.0/device_deph_scale)*(2.0*device_Eb/hbar);			// In [1/s]
	device_pol_frequency = device_bandgap/hbar;		// Frequency of p_k in [1/s]
	
	cout << "deph time = " << (1.0/device_deph_time)/fs << " [fs/rad]" << endl;
	cout << "pol freq  = " << (1.0/device_pol_frequency)/fs << " [fs/rad]" << endl;
	
	// Set used dipolemoment;
	device_dcv_hbar = device_dipolemoment/hbar;
	
	// Set wave number and density
	if (K == NULL)
	{
		K = new double[number_K_points];
		dK = new double[number_K_points];
		KdK = new double[number_K_points];
	}	
	// K': Type 1 Gauss-Chebychev Grid (a,b)
/*
	for(int i = 1; i <= number_K_points; i++)
	{
		// Type 1 GC grid
		double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_K_points));
		double wi = Pi/((double)number_K_points);
		double omega_zi = 1.0/sqrt(1.0-zi*zi);
		
		double k_a = 0.0;
		double k_b = K_max;
		K[number_K_points-i]  = k_a +  0.5*(k_b-k_a)*(zi+1.0);
		dK[number_K_points-i] = 0.5*(k_b-k_a)*wi*(1.0/omega_zi);

		KdK[number_K_points-i] = K[number_K_points-i]*dK[number_K_points-i];
	}
*/
	for(unsigned i = 0; i < number_K_points; i++) 
	{
		dK[i] = K_max/((double)number_K_points);	// Constant K density
	}
	K[0] = dK[0]/2.0;
	for(unsigned i = 1; i < number_K_points; i++)
	{
		K[i] = K[i-1] + (dK[i-1] + dK[i])/2.0;		// K points in middle of intervals
	}
	for(unsigned i = 0; i < number_K_points; i++) 
	{
		KdK[i] = K[i]*dK[i];
	}
	fileName << getToFileOutputKey() << "K_" << getName() << ".dat";
	saveBinary(fileName.str(), K, number_K_points);

	double lattice_tmp[3];
	lattice_tmp[0] = device_density;
	lattice_tmp[1] = device_background_temp_e;
	lattice_tmp[2] = device_background_temp_h;

	fileName.str(std::string());
	fileName << getToFileOutputKey() << "lattice_setup_" << getName() << ".dat";
	saveBinary(fileName.str(), lattice_tmp, 3);

	// Screening: Carrier-carrier 
	// GC+CS: (N_Q=50, N_eta=20) ~0.03% average relative error, all errors < 1 % For tests with Nk=100, Fermi, Fermi+hole, and Lorentz
	// Similar results for linear interpolation, but error can grow >10% near q=0 because of bad interpolation points in Lorentz
	number_Qsc_points = floor(number_K_points/2.0); // Screening vector
	number_Eta_points = floor(number_K_points/5.0); // Screening angular grid
	
	// Q: Type 1 Gauss-Chebychev Grid (a,b)
	if (Q_sc == NULL)
	{
		Q_sc	= new double[number_Qsc_points];
		dQ_sc	= new double[number_Qsc_points];
		Eta		= new double[number_Eta_points];
		dEta	= new double[number_Eta_points];
	}
	for(int i = 1; i <= number_Qsc_points; i++)
	{
		// Type 1 GC grid
		double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qsc_points));
		double wi = Pi/((double)number_Qsc_points);
		double omega_zi = 1.0/sqrt(1.0-zi*zi);
		
		double q_a = K[0];
		double q_b = K[number_K_points-1];
		Q_sc[number_Qsc_points-i]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
		dQ_sc[number_Qsc_points-i] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
	}
	
	// Q: Type 1 Gauss-Chebychev Grid (a,b)
	for(int i = 1; i <= number_Eta_points; i++)
	{
		// Type 1 GC grid
		double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Eta_points));
		double wi = Pi/((double)number_Eta_points);
		double omega_zi = 1.0/sqrt(1.0-zi*zi);
		
		double phi_a = 0.0;
		double phi_b = Pi;
		Eta[number_Eta_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
		dEta[number_Eta_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
	}

	// Initialize Ek_e, Ek_h
	if (Ek_e == NULL)
	{
		Ek_e = new double[number_K_points];
	}
	if (Ek_h == NULL)
	{
		Ek_h = new double[number_K_points];
	}
	if (sumE_k == NULL)
	{
		sumE_k = new double[number_K_points];
	}
	for(unsigned i = 0; i < number_K_points; i++)
	{
		Ek_e[i]   = hbar*hbar*K[i]*K[i]/(2.0*getElectronMass()*a0*a0);
		Ek_h[i]   = hbar*hbar*K[i]*K[i]/(2.0*getHoleMass()*a0*a0);
		sumE_k[i] = (device_bandgap + Ek_e[i] + Ek_h[i])/hbar;			// Pre-compute sum in [1/s]
	}
	
	#ifdef SBE_USE_RESONANT_PUMP_MODEL
	cout << "USING REAL PUMP MODEL: INITIAL DENSITY FIXED AT: n = " << device_pump_N0 << endl;
	device_density = device_pump_N0; // Force initial density
	#endif
	
	
	// Set chemical potential
	device_chemical_potential_e = misc_get_chem_potential_analytic_2d(device_density, device_me, device_background_temp_e);
	device_chemical_potential_h = misc_get_chem_potential_analytic_2d(device_density, device_mh, device_background_temp_h);

	// Set instantanious temperature initial values similarly
	device_inst_temp_e     = device_background_temp_e;
	device_inst_temp_h     = device_background_temp_h;
	device_inst_chem_pot_e = device_chemical_potential_e;
	device_inst_chem_pot_h = device_chemical_potential_h;
	
	device_inst_temp_e_prev = device_background_temp_e;
	device_inst_temp_h_prev = device_background_temp_h;

	// Initialize background fe_k, fh_k
	if (fe_k == NULL)
	{
		fe_k = new double[number_K_points];
		fh_k = new double[number_K_points];
	}
	if (fe_k_inst == NULL)
	{
		fe_k_inst = new double[number_K_points];
		fh_k_inst = new double[number_K_points];
	}
	if (carrier_scattering_rate_approximation_e == NULL)
	{
		carrier_scattering_rate_approximation_e = new double[number_K_points];
		carrier_scattering_rate_approximation_h = new double[number_K_points];
	}
	for(unsigned i = 0; i < number_K_points; i++)
	{
		fe_k[i]      = misc_get_fermi_distribution(device_background_temp_e, device_chemical_potential_e, device_me, K[i]);
		fh_k[i]      = misc_get_fermi_distribution(device_background_temp_h, device_chemical_potential_h, device_mh, K[i]);
		fe_k_inst[i] = misc_get_fermi_distribution(device_background_temp_e, device_chemical_potential_e, device_me, K[i]);
		fh_k_inst[i] = misc_get_fermi_distribution(device_background_temp_h, device_chemical_potential_h, device_mh, K[i]);

		carrier_scattering_rate_approximation_e[i] = 0;
		carrier_scattering_rate_approximation_h[i] = 0;
	}

/*	
	fileName.str(std::string());
	fileName << getToFileOutputKey() << "fe_" << getName() << ".dat";
	saveBinary(fileName.str(), fe_k, number_K_points);
	
	fileName.str(std::string());
	fileName << getToFileOutputKey() << "fh_" << getName() << ".dat";
	saveBinary(fileName.str(), fh_k, number_K_points);
*/
	
	// Initialize dynamical variables ne_k, nh_k
	if (ne_k == NULL)
	{
		ne_k = new double[number_K_points];
		nh_k = new double[number_K_points];
	}
	for(unsigned i = 0; i < number_K_points; i++)
	{
		ne_k[i] = fe_k[i];			// Initialize to background
		nh_k[i] = fh_k[i];			// Initialize to background
/*
		// Fermi+Hole
		double inversion = ne_k[i] + nh_k[i]-1;
		if(inversion>0.0)
		{
			double fwhm = 2.0; // [a0]
			double tmp = inversion/((K[i]-1.5)*(K[i]-1.5) + (0.5*fwhm)*(0.5*fwhm));
			ne_k[i] = ne_k[i] - tmp*ne_k[i]/(ne_k[i]+nh_k[i]);
			nh_k[i] = nh_k[i] - tmp*nh_k[i]/(ne_k[i]+nh_k[i]);
		}
*/
		
		// Lorentian
	//	double fwhm = 2.0; // [a0]
	//	ne_k[i] = 0.5/((K[i]-6.0)*(K[i]-6.0) + (0.5*fwhm)*(0.5*fwhm));
	//	nh_k[i] = 0.5/((K[i]-6.0)*(K[i]-6.0) + (0.5*fwhm)*(0.5*fwhm));
	}
	
	// Initialize p_k
	if (p_k == NULL)
	{
		p_k = new std::complex<double>[number_K_points];
	}
	for(unsigned i = 0; i < number_K_points; i++)
	{
		p_k[i] = 0;					// Initialize to zero
	}
	
	// Initialize RK4 TMP STORAGE
	// These arrays are passed into carrier scattering routine.
	if (p_k_tmp == NULL)
	{
		p_k_tmp = new std::complex<double>[number_K_points];
		ne_k_tmp = new double[number_K_points];
		nh_k_tmp = new double[number_K_points];
	}
	for(unsigned i = 0; i < number_K_points; i++)
	{
		p_k_tmp[i] = 0;
		ne_k_tmp[i] = 0;
		nh_k_tmp[i] = 0;
	}
	
	// Initialize temporary storage for RK4 steps
	int NUM_TMP_STORAGE = 4;
	if (p_k_quad == NULL)
	{
		p_k_quad = new std::complex<double>*[NUM_TMP_STORAGE];
		for(unsigned i = 0; i < NUM_TMP_STORAGE; i++)
		{
			p_k_quad[i] = new std::complex<double>[number_K_points];
		}
	}
	if (ne_k_quad == NULL)
	{
		ne_k_quad = new double*[NUM_TMP_STORAGE];
		for(unsigned i = 0; i < NUM_TMP_STORAGE; i++)
		{
			ne_k_quad[i] = new double[number_K_points];
		}
	}
	if (nh_k_quad == NULL)
	{
		nh_k_quad = new double*[NUM_TMP_STORAGE];
		for(unsigned i = 0; i < NUM_TMP_STORAGE; i++)
		{
			nh_k_quad[i] = new double[number_K_points];
		}
	}
	for(unsigned i = 0; i < NUM_TMP_STORAGE; i++)
	{
		for(unsigned j = 0; j < number_K_points; j++)
		{
			p_k_quad[i][j]  = 0;
			ne_k_quad[i][j] = 0;
			nh_k_quad[i][j] = 0;
		}
	}
	
	/* Storage for AB4 routine, not complete..
	if (p_k_ab4 == NULL)
	{
		p_k_ab4 = new Cyclic< std::complex<double> >(4,number_K_points);
	}
	if (ne_k_ab4 == NULL)
	{
		ne_k_ab4 = new Cyclic<double>(4,number_K_points);
	}
	if (nh_k_ab4 == NULL)
	{
		nh_k_ab4 = new Cyclic<double>(4,number_K_points);
	}
	ab4_precompute_counter = 0;
	p_k_ab4->zeroAllFields();
	ne_k_ab4->zeroAllFields();
	nh_k_ab4->zeroAllFields();
	*/
	
	// Renormalized energy and fields
	if (electric_field_renormalized == NULL)
	{
		electric_field_renormalized = new std::complex<double>[number_K_points];
		energy_renormalized_e = new double[number_K_points];
		energy_renormalized_h = new double[number_K_points];
	}
	for(unsigned i = 0; i < number_K_points; i++)
	{
		electric_field_renormalized[i] = 0;
		energy_renormalized_e[i] = 0;
		energy_renormalized_h[i] = 0;
	}
	
	//==========================
	// Create scaled variables
	//==========================
	//SBE_TIME_SCALE = 1.0/abs(sumE_k[number_K_points-1]-device_pol_frequency);
	SBE_TIME_SCALE = 1.0/abs(sumE_k[number_K_points-1]);
	device_dcv_hbar = SBE_TIME_SCALE*device_dcv_hbar;
	SBE_POL_DEPH = SBE_TIME_SCALE*device_deph_time;
	SBE_OCC_PUMP = SBE_TIME_SCALE/device_occ_relax_time;
	SBE_OCC_HOLE = SBE_TIME_SCALE/device_occ_hole_relax_time;
	for(unsigned i = 0; i < number_K_points; i++)
	{
		sumE_k[i] = SBE_TIME_SCALE*(sumE_k[i] - device_pol_frequency);
	}
/*
	fileName.str(std::string());
	fileName << getToFileOutputKey() << "SBE_TIME_SCALE_" << getName() << ".dat";
	saveBinary(fileName.str(), &SBE_TIME_SCALE, 1);
*/	
	//===============================================
	// Compute confinement functions

	if (coulomb_potential_confinement_ee == NULL)
	{
		coulomb_potential_confinement_ee = new double[number_K_points];
		coulomb_potential_confinement_hh = new double[number_K_points];
		coulomb_potential_confinement_eh = new double[number_K_points];
	}
	if (phonon_form_factor_ee == NULL)
	{
		phonon_form_factor_ee = new double[number_K_points];
		phonon_form_factor_hh = new double[number_K_points];
		phonon_form_factor_Q = new double[number_K_points];
	}
	
	double dQ_ph_ff = K_max/((double)number_K_points-1.0);
	phonon_form_factor_Q[0] = 0.0;
	for(unsigned i = 1; i < number_K_points; i++)
	{
		phonon_form_factor_Q[i] = phonon_form_factor_Q[i-1] + dQ_ph_ff;
	}

	#ifndef SBE_CONFINEMENT_USE_FILE
	// Assume perfect confinement
	number_Z_points    = 256; // Integration perpendicular to QW, along z-axis

	if (wavefunction_e == NULL)
	{
		wavefunction_e   = new double[number_Z_points];
		wavefunction_e_z = new double[number_Z_points];
		wavefunction_h   = new double[number_Z_points];
		wavefunction_h_z = new double[number_Z_points];
	}

	double dZ = device_length_qw/((double)number_Z_points);
	double Z0 = -0.5*device_length_qw + 0.5*dZ;
	for(int q = 0; q < number_K_points; q++)
	{
		double sum = 0.0;
		for(int i = 0; i < number_Z_points; i++)
		{
			double Zi = Z0 + ((double)i)*dZ; 
			double cosZi = cos((Pi/device_length_qw)*Zi);

			for(int j = 0; j < number_Z_points; j++)
			{
				double Zj = Z0 + ((double)j)*dZ;
				double cosZj = cos((Pi/device_length_qw)*Zj);

				// Note that 1/q is taken out of this expression
				sum += exp(-(K[q]/a0)*fabs(Zi-Zj))*cosZi*cosZi*cosZj*cosZj;
			}

		}
		sum *= dZ*dZ*(2.0/device_length_qw)*(2.0/device_length_qw);
		coulomb_potential_confinement_ee[q] = sum;
		coulomb_potential_confinement_hh[q] = sum;
		coulomb_potential_confinement_eh[q] = sum;
	}
	
	for(int q = 0; q < number_K_points; q++)
	{
		std::complex<double> sum_ff = 0.0;
		for(int i = 0; i < number_Z_points; i++)
		{
			double Zi = Z0 + ((double)i)*dZ;
			double cosZi = cos((Pi/device_length_qw)*Zi);
			sum_ff += exp(I*(phonon_form_factor_Q[q]/a0)*Zi)*cosZi*cosZi;
		}
		sum_ff *= dZ*(2.0/device_length_qw);
		phonon_form_factor_ee[q] = abs(sum_ff)*abs(sum_ff);
		phonon_form_factor_hh[q] = abs(sum_ff)*abs(sum_ff);
	}
/*	
	fileName.str(std::string());
	fileName << getToFileOutputKey() << "coulomb_confinement_e_z_" << getName() << ".dat";
	std::string tmp = fileName.str();
	FILE *coulomb_conf_function = fopen(tmp.c_str(),"w+");
	for(int j = 0; j < number_Z_points; j++)
	{
		double Zj = Z0 + ((double)j)*dZ;
		double cosZj = sqrt(2.0/device_length_qw)*cos((Pi/device_length_qw)*Zj);
		wavefunction_e[j]   = cosZj;
		wavefunction_e_z[j] = Zj;

		fprintf(coulomb_conf_function,"%.6e %.6e\n",Zj,cosZj);
	}
	fclose(coulomb_conf_function);

	fileName.str(std::string());
	fileName << getToFileOutputKey() << "coulomb_confinement_h_z_" << getName() << ".dat";
	tmp = fileName.str();
	coulomb_conf_function = fopen(tmp.c_str(),"w+");
	for(int j = 0; j < number_Z_points; j++)
	{
		double Zj = Z0 + ((double)j)*dZ;
		double cosZj = sqrt(2.0/device_length_qw)*cos((Pi/device_length_qw)*Zj);
		wavefunction_h[j]   = cosZj;
		wavefunction_h_z[j] = Zj;

		fprintf(coulomb_conf_function,"%.6e %.6e\n",Zj,cosZj);
	}
	fclose(coulomb_conf_function);
*/
	
	#else
	// Import wavefunctions from file
	// The variable number_Z_points should match the number of points in the file...
	number_Z_points    = 120; // Number of points in file

	// Assumes the formatting "z f(z)" for each line where f(z) and z are doubles, z is in untis of nm
	// If you use complex wavefunctions,you will need to redfine the import, storage, and integration below

	if (wavefunction_e == NULL)
	{
		wavefunction_e   = new double[number_Z_points];
		wavefunction_e_z = new double[number_Z_points];
		wavefunction_h   = new double[number_Z_points];
		wavefunction_h_z = new double[number_Z_points];
	}
/*
	fileName.str(std::string());
	fileName << "wavefunction_e_" << number_Z_points << ".dat";
	std::string tmp = fileName.str();
	FILE *fid = fopen(tmp.c_str(),"r+");
	if (fid != NULL)
	{
		for(int j = 0; j < number_Z_points; j++)
		{
			fscanf(fid,"%le %le",&(wavefunction_e_z[j]), &(wavefunction_e[j]));
			wavefunction_e_z[j] *= 1.0e-9; // to [nm]
		}
		fclose(fid);
	} else {
		cout << "ERROR: failed to read file: " << fileName.str() << endl;
		exit(-1);
	}
	
	fileName.str(std::string());
	fileName << "wavefunction_h_" << number_Z_points << ".dat";
	tmp = fileName.str();
	fid = fopen(tmp.c_str(),"r+");
	if (fid != NULL)
	{
		for(int j = 0; j < number_Z_points; j++)
		{
			fscanf(fid,"%le %le",&(wavefunction_h_z[j]), &(wavefunction_h[j]));
			wavefunction_h_z[j] *= 1.0e-9; // to [nm]
		}
		fclose(fid);
	} else {
		cout << "ERROR: failed to read file: " << fileName.str() << endl;
		exit(-1);
	}
	*/
	
	fileName.str(std::string());
	fileName << "wavefunctions_eh_" << number_Z_points << ".dat";
	std::string tmp = fileName.str();
	FILE *fid = fopen(tmp.c_str(),"r+");
	if (fid != NULL)
	{
		for(int j = 0; j < number_Z_points; j++)
		{
			int id = fscanf(fid,"%le %le %le",&(wavefunction_e_z[j]), &(wavefunction_e[j]), &(wavefunction_h[j]));
			wavefunction_h_z[j] = wavefunction_e_z[j];
			//wavefunction_e_z[j] *= 1.0e-9; // to [nm]
		}
		fclose(fid);
	} else {
		cout << "ERROR: failed to read file: " << fileName.str() << endl;
		exit(-1);
	}

	// Normalize wavefunctions in case they arent
	double normalize_e = 0.0;
	double normalize_h = 0.0;
	for(int i = 0; i < number_Z_points; i++)
	{
		double dZe = 0.0;
		if (i==0)
		{
			dZe = 0.5*(wavefunction_e_z[1] - wavefunction_e_z[0]);
		} else if (i==number_Z_points-1) {
			dZe = 0.5*(wavefunction_e_z[number_Z_points-1] - wavefunction_e_z[number_Z_points-2]);
		} else {
			dZe = 0.5*(wavefunction_e_z[i+1] - wavefunction_e_z[i-1]);
		}
		normalize_e += wavefunction_e[i]*wavefunction_e[i]*dZe;

		double dZh = 0.0;
		if (i==0)
		{
			dZh = 0.5*(wavefunction_h_z[1] - wavefunction_h_z[0]);
		} else if (i==number_Z_points-1) {
			dZh = 0.5*(wavefunction_h_z[number_Z_points-1] - wavefunction_h_z[number_Z_points-2]);
		} else {
			dZh = 0.5*(wavefunction_h_z[i+1] - wavefunction_h_z[i-1]);
		}
		normalize_h += wavefunction_h[i]*wavefunction_h[i]*dZh;
	}

	for(int i = 0; i < number_Z_points; i++)
	{
		wavefunction_e[i] /= sqrt(normalize_e);
		wavefunction_h[i] /= sqrt(normalize_h);
	}

	for(int q = 0; q < number_K_points; q++)
	{
		// ee
		double sum = 0.0;
		for(int i = 0; i < number_Z_points; i++)
		{
			double Zi  = wavefunction_e_z[i];
			double fZi = wavefunction_e[i];

			double dZ_i = 0.0;
			if (i==0)
			{
				dZ_i = 0.5*(wavefunction_e_z[1] - wavefunction_e_z[0]);
			} else if (i==number_Z_points-1) {
				dZ_i = 0.5*(wavefunction_e_z[number_Z_points-1] - wavefunction_e_z[number_Z_points-2]);
			} else {
				dZ_i = 0.5*(wavefunction_e_z[i+1] - wavefunction_e_z[i-1]);
			}

			for(int j = 0; j < number_Z_points; j++)
			{
				double Zj  = wavefunction_e_z[j];
				double fZj = wavefunction_e[j];

				double dZ_j = 0.0;
				if (j==0)
				{
					dZ_j = 0.5*(wavefunction_e_z[1] - wavefunction_e_z[0]);
				} else if (j==number_Z_points-1) {
					dZ_j = 0.5*(wavefunction_e_z[number_Z_points-1] - wavefunction_e_z[number_Z_points-2]);
				} else {
					dZ_j = 0.5*(wavefunction_e_z[j+1] - wavefunction_e_z[j-1]);
				}

				// Note that 1/q is taken out of this expression
				sum += exp(-(K[q]/a0)*fabs(Zi-Zj))*fZi*fZi*fZj*fZj*dZ_i*dZ_j;
			}
		}
		coulomb_potential_confinement_ee[q] = sum;
		
		// hh
		sum = 0.0;
		for(int i = 0; i < number_Z_points; i++)
		{
			double Zi  = wavefunction_h_z[i];
			double fZi = wavefunction_h[i];

			double dZ_i = 0.0;
			if (i==0)
			{
				dZ_i = 0.5*(wavefunction_h_z[1] - wavefunction_h_z[0]);
			} else if (i==number_Z_points-1) {
				dZ_i = 0.5*(wavefunction_h_z[number_Z_points-1] - wavefunction_h_z[number_Z_points-2]);
			} else {
				dZ_i = 0.5*(wavefunction_h_z[i+1] - wavefunction_h_z[i-1]);
			}

			for(int j = 0; j < number_Z_points; j++)
			{
				double Zj  = wavefunction_h_z[j];
				double fZj = wavefunction_h[j];

				double dZ_j = 0.0;
				if (j==0)
				{
					dZ_j = 0.5*(wavefunction_h_z[1] - wavefunction_h_z[0]);
				} else if (j==number_Z_points-1) {
					dZ_j = 0.5*(wavefunction_h_z[number_Z_points-1] - wavefunction_h_z[number_Z_points-2]);
				} else {
					dZ_j = 0.5*(wavefunction_h_z[j+1] - wavefunction_h_z[j-1]);
				}

				// Note that 1/q is taken out of this expression
				sum += exp(-(K[q]/a0)*fabs(Zi-Zj))*fZi*fZi*fZj*fZj*dZ_i*dZ_j;
			}
		}
		coulomb_potential_confinement_hh[q] = sum;

		// eh and he
		sum = 0.0;
		for(int i = 0; i < number_Z_points; i++)
		{
			double Zi  = wavefunction_h_z[i];
			double fZi = wavefunction_h[i];

			double dZ_i = 0.0;
			if (i==0)
			{
				dZ_i = 0.5*(wavefunction_h_z[1] - wavefunction_h_z[0]);
			} else if (i==number_Z_points-1) {
				dZ_i = 0.5*(wavefunction_h_z[number_Z_points-1] - wavefunction_h_z[number_Z_points-2]);
			} else {
				dZ_i = 0.5*(wavefunction_h_z[i+1] - wavefunction_h_z[i-1]);
			}

			for(int j = 0; j < number_Z_points; j++)
			{
				double Zj  = wavefunction_e_z[j];
				double fZj = wavefunction_e[j];

				double dZ_j = 0.0;
				if (j==0)
				{
					dZ_j = 0.5*(wavefunction_e_z[1] - wavefunction_e_z[0]);
				} else if (j==number_Z_points-1) {
					dZ_j = 0.5*(wavefunction_e_z[number_Z_points-1] - wavefunction_e_z[number_Z_points-2]);
				} else {
					dZ_j = 0.5*(wavefunction_e_z[j+1] - wavefunction_e_z[j-1]);
				}

				// Note that 1/q is taken out of this expression
				sum += exp(-(K[q]/a0)*fabs(Zi-Zj))*fZi*fZi*fZj*fZj*dZ_i*dZ_j;
			}
		}
		coulomb_potential_confinement_eh[q] = sum;
	}
	
	for(int q = 0; q < number_K_points; q++)
	{
		// Phonon Form factor
		std::complex<double> sum_ff = 0.0;
		for(int i = 0; i < number_Z_points; i++)
		{
			double Zi  = wavefunction_e_z[i];
			double fZi = wavefunction_e[i];

			double dZ_i = 0.0;
			if (i==0)
			{
				dZ_i = 0.5*(wavefunction_e_z[1] - wavefunction_e_z[0]);
			} else if (i==number_Z_points-1) {
				dZ_i = 0.5*(wavefunction_e_z[number_Z_points-1] - wavefunction_e_z[number_Z_points-2]);
			} else {
				dZ_i = 0.5*(wavefunction_e_z[i+1] - wavefunction_e_z[i-1]);
			}

			sum_ff += exp(I*(phonon_form_factor_Q[q]/a0)*Zi)*fZi*fZi*dZ_i;
		}
		phonon_form_factor_ee[q] = abs(sum_ff)*abs(sum_ff);

		sum_ff = 0.0;
		for(int i = 0; i < number_Z_points; i++)
		{
			double Zi  = wavefunction_h_z[i];
			double fZi = wavefunction_h[i];

			double dZ_i = 0.0;
			if (i==0)
			{
				dZ_i = 0.5*(wavefunction_h_z[1] - wavefunction_h_z[0]);
			} else if (i==number_Z_points-1) {
				dZ_i = 0.5*(wavefunction_h_z[number_Z_points-1] - wavefunction_h_z[number_Z_points-2]);
			} else {
				dZ_i = 0.5*(wavefunction_h_z[i+1] - wavefunction_h_z[i-1]);
			}

			sum_ff += exp(I*(phonon_form_factor_Q[q]/a0)*Zi)*fZi*fZi*dZ_i;
		}
		phonon_form_factor_hh[q] = abs(sum_ff)*abs(sum_ff);
	}
	//---
/*	
	fileName.str(std::string());
	fileName << getToFileOutputKey() << "coulomb_confinement_e_z_" << getName() << ".dat";
	tmp = fileName.str();
	FILE *coulomb_conf_function = fopen(tmp.c_str(),"w+");
	for(int j = 0; j < number_Z_points; j++)
	{
		fprintf(coulomb_conf_function,"%.6e %.6e\n",wavefunction_e_z[j], wavefunction_e[j]);
	}
	fclose(coulomb_conf_function);

	fileName.str(std::string());
	fileName << getToFileOutputKey() << "coulomb_confinement_h_z_" << getName() << ".dat";
	tmp = fileName.str();
	coulomb_conf_function = fopen(tmp.c_str(),"w+");
	for(int j = 0; j < number_Z_points; j++)
	{
		fprintf(coulomb_conf_function,"%.6e %.6e\n",wavefunction_h_z[j], wavefunction_h[j]);
	}
	fclose(coulomb_conf_function);
*/

	#endif
	
/*
	fileName.str(std::string());
	fileName << getToFileOutputKey() << "coulomb_confinement_ee_k_" << getName() << ".dat";
	tmp = fileName.str();
	FILE *coulomb_conf_function_k = fopen(tmp.c_str(),"w+");
	for(int q = 0; q < number_K_points; q++)
	{
		fprintf(coulomb_conf_function_k,"%.6e %.6e\n",K[q],coulomb_potential_confinement_ee[q]);
	}
	fclose(coulomb_conf_function_k);

	fileName.str(std::string());
	fileName << getToFileOutputKey() << "coulomb_confinement_hh_k_" << getName() << ".dat";
	tmp = fileName.str();
	coulomb_conf_function_k = fopen(tmp.c_str(),"w+");
	for(int q = 0; q < number_K_points; q++)
	{
		fprintf(coulomb_conf_function_k,"%.6e %.6e\n",K[q],coulomb_potential_confinement_hh[q]);
	}
	fclose(coulomb_conf_function_k);

	fileName.str(std::string());
	fileName << getToFileOutputKey() << "coulomb_confinement_eh_k_" << getName() << ".dat";
	tmp = fileName.str();
	coulomb_conf_function_k = fopen(tmp.c_str(),"w+");
	for(int q = 0; q < number_K_points; q++)
	{
		fprintf(coulomb_conf_function_k,"%.6e %.6e\n",K[q],coulomb_potential_confinement_eh[q]);
	}
	fclose(coulomb_conf_function_k);

	fileName.str(std::string());
	fileName << getToFileOutputKey() << "phonon_form_factor_ee_k_" << getName() << ".dat";
	tmp = fileName.str();
	FILE *phonon_form_factor_k = fopen(tmp.c_str(),"w+");
	for(int q = 0; q < number_K_points; q++)
	{
		fprintf(phonon_form_factor_k,"%.6e %.6e\n",phonon_form_factor_Q[q],phonon_form_factor_ee[q]);
	}
	fclose(phonon_form_factor_k);

	fileName.str(std::string());
	fileName << getToFileOutputKey() << "phonon_form_factor_hh_k_" << getName() << ".dat";
	tmp = fileName.str();
	phonon_form_factor_k = fopen(tmp.c_str(),"w+");
	for(int q = 0; q < number_K_points; q++)
	{
		fprintf(phonon_form_factor_k,"%.6e %.6e\n",phonon_form_factor_Q[q],phonon_form_factor_hh[q]);
	}
	fclose(phonon_form_factor_k);
*/
	
	//==================================
	// Full 2D Coulomb matrix

	if (coulomb_matrix_ee == NULL)
	{
		coulomb_matrix_ee = new double*[number_K_points];
		for(unsigned i = 0; i < number_K_points; i++)
		{
			coulomb_matrix_ee[i] = new double[number_K_points];
		}
	}
	if (coulomb_matrix_hh == NULL)
	{
		coulomb_matrix_hh = new double*[number_K_points];
		for(unsigned i = 0; i < number_K_points; i++)
		{
			coulomb_matrix_hh[i] = new double[number_K_points];
		}
	}
	if (coulomb_matrix_eh == NULL)
	{
		coulomb_matrix_eh = new double*[number_K_points];
		for(unsigned i = 0; i < number_K_points; i++)
		{
			coulomb_matrix_eh[i] = new double[number_K_points];
		}
	}

	for(int i = 0; i < number_K_points; i++)
	{
		for(int j = 0; j < number_K_points; j++)
		{
			coulomb_matrix_ee[i][j] = 0;
			coulomb_matrix_hh[i][j] = 0;
			coulomb_matrix_eh[i][j] = 0;
		}
	}

	//========================================
	// Reduced 2D coulomb matrix
	coulomb_matrix_red_index_i = 3;
	coulomb_matrix_red_index_j = 5;
	
	coulomb_matrix_red_size = 0;
	for(int i=0; i<number_K_points; i+=coulomb_matrix_red_index_i )
	{
		for( int iiq = i+1; iiq < number_K_points ; iiq+=coulomb_matrix_red_index_j )
		{
			coulomb_matrix_red_size++;
		}
		for( int iiq = i-1; iiq >=0 ; iiq-=coulomb_matrix_red_index_j )
		{
			coulomb_matrix_red_size++;
		}
	}
	
	if (coulomb_matrix_ee_red == NULL)
	{
		coulomb_matrix_ee_red =  new double[coulomb_matrix_red_size];
		coulomb_matrix_hh_red =  new double[coulomb_matrix_red_size];
		coulomb_matrix_eh_red =  new double[coulomb_matrix_red_size];
		if( coulomb_matrix_ee_red == NULL ){
			printf("allocation of coulomb_matrix_ee_red failed\n");
			exit(-1);
		}
		if( coulomb_matrix_hh_red == NULL ){
			printf("allocation of coulomb_matrix_hh_red failed\n");
			exit(-1);
		}
		if( coulomb_matrix_eh_red == NULL ){
			printf("allocation of coulomb_matrix_eh_red failed\n");
			exit(-1);
		}
	}

	for(int i = 0; i < coulomb_matrix_red_size; i++)
	{
		coulomb_matrix_ee_red[i] = 0;
		coulomb_matrix_hh_red[i] = 0;
		coulomb_matrix_eh_red[i] = 0;
	}

	//==========================================
	// Calculate screening of coulomb potential
	//
	
	if (coulomb_potential_normalized_ee == NULL)
	{
		coulomb_potential_normalized_ee = new double[number_K_points];
		coulomb_potential_normalized_hh = new double[number_K_points];
		coulomb_potential_normalized_eh = new double[number_K_points];
		coulomb_potential_epsilon_inv = new double[number_K_points];
	}
	
	
	#ifdef SBE_USE_FULL_SCREENING
	// Fast algorithm for screening is about 33 times faster
	// Find number of columns in each row
	int screening_col[number_K_points];
	for(int k = 0; k < number_K_points; k++)
	{
		screening_col[k] = 0;
		for(int q = 0; q < number_Qsc_points; q++)
		{
			for(int th_q = 0; th_q < number_Eta_points; th_q++)
			{
				double angle_th_k_q  = Eta[th_q];
				double k_q           = sqrt(K[k ]*K[k ] + Q_sc[q]*Q_sc[q] - 2.0*K[k ]*Q_sc[q]*cos(angle_th_k_q));
				
				if ((k_q >= K[0])&&(k_q <= K[number_K_points-1]))
				{
					screening_col[k] += 1;
				}
			}
		}
	}

	if (coulomb_potential_epsilon_index == NULL)
	{
		coulomb_potential_epsilon_index = new indexSet2d_screening(number_K_points, screening_col);
	}

	int NUM_SCREEN = 0.0;
	
	// Pre calculation of screening
	// Only required if using fast algorithm 
	for(int k = 0; k < number_K_points; k++)
	{
		int count_screen_k = 0;
		
		// coulomb matrix element
		double col_fact = e*e/(2.0*eps0*eps);
		double Ve_q      = (col_fact/(K[k]/a0))*coulomb_potential_confinement_ee[k];
		double Vh_q      = (col_fact/(K[k]/a0))*coulomb_potential_confinement_hh[k];
		//double Veh_q      = (col_fact/(K[k]/a0))*coulomb_potential_confinement_eh[k];
		
		
		//----
		for(int q = 0; q < number_Qsc_points; q++)
		{
			double fact_e_q_tmp = 0.0;
			double fact_h_q_tmp = 0.0;
			
			for(int th_q = 0; th_q < number_Eta_points; th_q++)
			{
				double angle_th_k_q  = Eta[th_q];
				double k_q           = sqrt(K[k ]*K[k ] + Q_sc[q]*Q_sc[q] - 2.0*K[k ]*Q_sc[q]*cos(angle_th_k_q));
				
				if ((k_q >= K[0])&&(k_q <= K[number_K_points-1]))
				{
					int k_q_IND, q_IND;
					double k_q_b1,k_q_b2,k_q_b3;
					double q_b1,q_b2,q_b3;
					
					// Electrons
					double fe_k_q	= misc_interpolate_K_array_parabolic_index(k_q     , K, ne_k, number_K_points, &k_q_IND, &k_q_b1, &k_q_b2, &k_q_b3);
					double fe_q	= misc_interpolate_K_array_parabolic_index(Q_sc[q] , K, ne_k, number_K_points, &q_IND  , &q_b1  , &q_b2  , &q_b3);
					//double fe_k_q  = misc_interpolate_K_array_index(k_q, K, ne_k, number_K_points,0.0,&k_q_alpha,&k_q_IND,&k_q_beta);
					//double fe_q    = misc_interpolate_K_array_index(Q_sc[q], K, ne_k, number_K_points,0.0,&q_alpha,&q_IND,&q_beta);
					double Ee_k_q  = hbar*hbar*(k_q*k_q)/(2.0*device_me*a0*a0);
					double Ee_q    = hbar*hbar*Q_sc[q]*Q_sc[q]/(2.0*device_me*a0*a0);

					// Holes
					double Eh_k_q  = hbar*hbar*(k_q*k_q)/(2.0*device_mh*a0*a0);
					double Eh_q    = hbar*hbar*Q_sc[q]*Q_sc[q]/(2.0*device_mh*a0*a0);
					
					double spin_const = 2.0; // Multiply by 2 to account for the spin
					double integration_constant = Q_sc[q]*dQ_sc[q]*dEta[th_q]*spin_const/(2.0*Pi*Pi*a0*a0); // Multiply by 2.0 for half angular grid
					double const_e = Ve_q*integration_constant/(Ee_k_q - Ee_q);
					double const_h = Vh_q*integration_constant/(Eh_k_q - Eh_q);
					
					coulomb_potential_epsilon_index->updateWeight(k,count_screen_k, const_e, const_h, q_IND, q_b1, q_b2, q_b3, k_q_IND, k_q_b1, k_q_b2, k_q_b3);
					count_screen_k++;
					
					NUM_SCREEN += 1;
				}
			}
			
		}
	}

	cout << "Screening: Number of indicies = " << NUM_SCREEN << endl;

	// Sort index arrays. This has a small influence on the execution time, probably because of low number of elements
	coulomb_potential_epsilon_index->sortIndices();
	
	// Fill coulomb_potential_epsilon_inv[], coulomb_potential_normalized_ee[], coulomb_potential_normalized_hh[], and coulomb_potential_normalized_eh[]
	//carrier_scattering_screening_full(ne_k, nh_k); // Old screening, does not require precalculation
	//carrier_scattering_screening_full_fast(ne_k, nh_k);
	updateScreening_reducedCoulomb_full(ne_k, nh_k);

	#else
		
	// Fill coulomb_potential_epsilon_inv[], coulomb_potential_normalized_ee[], coulomb_potential_normalized_hh[], and coulomb_potential_normalized_eh[]
	updateScreening_StaticPlasmonApproximation(ne_k, nh_k);
	//updateScreening_reducedCoulomb_StaticPlasmonApproximation(ne_k, nh_k);

	#endif
	
/*
	fileName.str(std::string());
	fileName << getToFileOutputKey() << "coulomb_screening_inv_k_" << getName() << ".dat";
	tmp = fileName.str();
	FILE *coulomb_screening_k = fopen(tmp.c_str(),"w+");
	for(int q = 0; q < number_K_points; q++)
	{
		fprintf(coulomb_screening_k,"%.6e %.6e\n",K[q],coulomb_potential_epsilon_inv[q]);
	}
	fclose(coulomb_screening_k);

	
	fileName.str(std::string());
	fileName << getToFileOutputKey() << "coulomb_potential_normalized_ee_k_" << getName() << ".dat";
	tmp = fileName.str();
	FILE *coulomb_potential_k = fopen(tmp.c_str(),"w+");
	for(int q = 0; q < number_K_points; q++)
	{
		fprintf(coulomb_potential_k,"%.6e %.6e\n",K[q],coulomb_potential_normalized_ee[q]);
	}
	fclose(coulomb_potential_k);

	fileName.str(std::string());
	fileName << getToFileOutputKey() << "coulomb_potential_normalized_hh_k_" << getName() << ".dat";
	tmp = fileName.str();
	coulomb_potential_k = fopen(tmp.c_str(),"w+");
	for(int q = 0; q < number_K_points; q++)
	{
		fprintf(coulomb_potential_k,"%.6e %.6e\n",K[q],coulomb_potential_normalized_hh[q]);
	}
	fclose(coulomb_potential_k);

	fileName.str(std::string());
	fileName << getToFileOutputKey() << "coulomb_potential_normalized_eh_k_" << getName() << ".dat";
	tmp = fileName.str();
	coulomb_potential_k = fopen(tmp.c_str(),"w+");
	for(int q = 0; q < number_K_points; q++)
	{
		fprintf(coulomb_potential_k,"%.6e %.6e\n",K[q],coulomb_potential_normalized_eh[q]);
	}
	fclose(coulomb_potential_k);
*/


	// ==========================================
	// Initialize Carrier scattering

	if (carrier_scattering_rates_e_in ==NULL)
	{
		carrier_scattering_rates_e_in  = new double[number_K_points];
		carrier_scattering_rates_e_out = new double[number_K_points];
		carrier_scattering_rates_h_in  = new double[number_K_points];
		carrier_scattering_rates_h_out = new double[number_K_points];
		carrier_scattering_rates_e_total = new double[number_K_points];
		carrier_scattering_rates_h_total = new double[number_K_points];


		carrier_scattering_rates_phonon_e_in    = new double[number_K_points];
		carrier_scattering_rates_phonon_e_out   = new double[number_K_points];
		carrier_scattering_rates_phonon_h_in    = new double[number_K_points];
		carrier_scattering_rates_phonon_h_out   = new double[number_K_points];
		carrier_scattering_rates_phonon_e_total = new double[number_K_points];
		carrier_scattering_rates_phonon_h_total = new double[number_K_points];
	}
	for(int i = 0; i < number_K_points; i++)
	{
		carrier_scattering_rates_e_in[i]  = 0.0;
		carrier_scattering_rates_e_out[i] = 0.0;
		carrier_scattering_rates_h_in[i]  = 0.0;
		carrier_scattering_rates_h_out[i] = 0.0;
		carrier_scattering_rates_e_total[i] = 0.0;
		carrier_scattering_rates_h_total[i] = 0.0;


		carrier_scattering_rates_phonon_e_in[i] = 0.0;
		carrier_scattering_rates_phonon_e_out[i] = 0.0;
		carrier_scattering_rates_phonon_e_total[i] = 0.0;
		carrier_scattering_rates_phonon_h_total[i] = 0.0;
	}
	
	// When to skip re-calculation of carrier scattering rates
	// Tested with spont. emissions, carrier relaxation, on a strong 100fs pulse
	// all measurements 8ps after pulse relaxation
	// For TOL <= 0.05 the Rel. Carrier error is <=1% in the k-range [0-5]
	//    and can get up to 30% at k=15
	//    and the influence on the pulse is < 0.01% in amplitude and < 1e-6 in phase
	//CARRIER_DEVIATION_TOL = 0.05; // [0,5]
	CARRIER_DEVIATION_TOL = 0.05; // [0,5]
	if (carrier_scattering_prev_ne == NULL)
	{
		carrier_scattering_prev_ne = new double[number_K_points];
		carrier_scattering_prev_nh = new double[number_K_points];
	}
	for(int i = 0; i < number_K_points; i++)
	{
		carrier_scattering_prev_ne[i] = 0.0;
		carrier_scattering_prev_nh[i] = 0.0;
	}
	if (renormalize_prev_ne == NULL)
	{
		renormalize_prev_ne = new double[number_K_points];
		renormalize_prev_nh = new double[number_K_points];
	}
	for(int i = 0; i < number_K_points; i++)
	{
		renormalize_prev_ne[i] = 0.0;
		renormalize_prev_nh[i] = 0.0;
	}

	#if defined(USE_ISAK_HOLE_FILLING_TABLE) || defined(USE_ISAK_HOLE_FILLING)


	// Initialize scattering rate
	carrier_scattering_init(ne_k, nh_k); 
	
	// Calculation background scattering
	carrier_scattering_method2_parallel(ne_k, nh_k);
	#endif

	//===================================
	// For hysteris experiment
	hyst_tnext = 0.0;
	hyst_count = 0;
	hyst_dens_initial = device_density;
	
	//================================
	// Update 2D Coulomb matrix
	updateCoulombMatrix();

	//=================================
	// Fill reduced Coulomb matrix
	#ifdef SBE_USE_FULL_SCREENING
	updateScreening_reducedCoulomb_full(ne_k, nh_k);
	#else
	updateScreening_reducedCoulomb_StaticPlasmonApproximation(ne_k, nh_k);
	#endif
	
	//==============================================
	// Calculate Energy and Field renormalization
	renormalize_reducedCoulomb(p_k, ne_k, nh_k);
	renormalize_reducedCoulomb_Energy(ne_k, nh_k); // Compute this once, in case it is not updated
	
	#ifdef SBE_USE_SPONTAN_EMISSIONS
		
	
		srand(time(NULL));
		unsigned seed2 = (unsigned) rand();
		
		#if defined(__ICC) || defined(__INTEL_COMPILER)
			
			vslNewStream( &stream, VSL_BRNG_SFMT19937, seed2); // ICC or ICPC
		#elif defined(__GNUC__) || defined(__GNUG__)
			// Generate random seed
			std::hash<std::string> hash_fn;
			seed2 += (unsigned) hash_fn(getName());
			
			device_generator.seed(seed2);
		#endif
	
		// Spontanious emissions variables
		if (device_spont_emission_wk == NULL)
		{
			device_spont_emission_wk = new double[number_K_points];
		
			double nb = n_vcsel_bgr;
			double dcv = device_dcv_hbar/SBE_TIME_SCALE;
			double gamma_spont = SBE_TIME_SCALE*(dcv*dcv)*hbar*nb*nb*nb/(Pi*Pi*eps*eps0*c0*c0*c0);
		
			for (int i = 0; i < number_K_points; i++)
			{
				double wk = device_pol_frequency + sumE_k[i]/SBE_TIME_SCALE;
				device_spont_emission_wk[i] = -gamma_spont*wk*wk*wk;

				// In case the emissions arent used in this QW, set during initialization
				if (device_spont_emission_on == false )
				{
					device_spont_emission_wk[i] = 0.0;
				}
			}
		}
	#endif

	//================
	// Pump variables
	//================
	#ifdef SBE_USE_RESONANT_PUMP_MODEL
		sbe_set_real_pump_model(device_pump_W0, device_pump_E0, device_pump_ETA, DT);
	#endif


	//===================================================
	// Interpolation: Calculate carrier scattering table
	
	#ifdef USE_ISAK_HOLE_FILLING_TABLE
	double N_start = 1.0; // *1.0e16;
	double N_stop = 5.5; // *1.0e16;
	double N_step = 0.05; // *1.0e16;
	
	double Te_start = 250;
	double Te_stop = 3000;
	double Te_step = 25;
	
	double Th_start = 250;
	double Th_stop = 3000;
	double Th_step = 25;
	
	fileName.str(std::string());
	fileName << "dynamicInterpolationTable_e_in_" << getName() << ".dat";
	carrier_scattering_table_e_in 	= new DynamicInterpolation(fileName.str(),number_K_points,Te_start,Te_stop,Te_step,Th_start,Th_stop,Th_step,N_start,N_stop,N_step);
	fileName.str(std::string());
	fileName << "dynamicInterpolationTable_e_out_" << getName() << ".dat";
	carrier_scattering_table_e_out 	= new DynamicInterpolation(fileName.str(),number_K_points,Te_start,Te_stop,Te_step,Th_start,Th_stop,Th_step,N_start,N_stop,N_step);
	fileName.str(std::string());
	fileName << "dynamicInterpolationTable_h_in_" << getName() << ".dat";
	carrier_scattering_table_h_in 	= new DynamicInterpolation(fileName.str(),number_K_points,Te_start,Te_stop,Te_step,Th_start,Th_stop,Th_step,N_start,N_stop,N_step);
	fileName.str(std::string());
	fileName << "dynamicInterpolationTable_h_out_" << getName() << ".dat";
	carrier_scattering_table_h_out 	= new DynamicInterpolation(fileName.str(),number_K_points,Te_start,Te_stop,Te_step,Th_start,Th_stop,Th_step,N_start,N_stop,N_step);
	carrier_scattering_table_gridpoints_numbers = 0;
	carrier_scattering_table_gridpoints = new double*[8];
	for(int i = 0; i < 8; i++)
	{
		carrier_scattering_table_gridpoints[i] = new double[6];
		for(int j = 0; j < 6; j++)
		{
			carrier_scattering_table_gridpoints[i][j] = 0;
		}
	}

	if (device_density > 1.0e16)
	{
		// Run once
		carrier_scattering_method2_table(ne_k, nh_k);
	}
	#endif

	//=============================
	// To help compute load for MPI
	#ifdef MPI_BALANCE_WORKLOAD
	MPI_load = new myTimer(getName());
	#endif
}


/* Set a delay for when computations are supposed to start
 * Used on in special situations
 * */
void Device::sbe_delay_computations(double delay)
{
	computationalDelay = delay;
}

/* Evaluate RHS of SBE: Polarization
 * Fills the array p_k_new
 * */
void Device::sbe_RHS_pk(double DT, std::complex<double> *p_k_new, std::complex<double> *p_k_old, double *ne_k_old, double *nh_k_old, std::complex<double> E_tmp)
{

	double scale_factor = SBE_TIME_SCALE/hbar;
	
	if (device_density > 1.0e16)
	{
		for(unsigned i = 0; i < number_K_points; i++)
		{
			// Scaled equations
			p_k_new[i] = -I*(sumE_k[i] - energy_renormalized_e[i] - energy_renormalized_h[i])*p_k_old[i]
						 -SBE_POL_DEPH*p_k_old[i]
						 -I*(ne_k_old[i] + nh_k_old[i] - 1.0)*(E_tmp + electric_field_renormalized[i]);
			
			#ifdef SBE_USE_SPONTAN_EMISSIONS
				p_k_new[i] += -device_spont_emission_phasor*device_spont_emission_wk[i]*ne_k_old[i]*nh_k_old[i];
			#endif
		}
	} else {
		for(unsigned i = 0; i < number_K_points; i++)
		{
			// Scaled equations
			p_k_new[i] = -I*(sumE_k[i] - energy_renormalized_e[i] - energy_renormalized_h[i])*p_k_old[i]
						 -SBE_POL_DEPH*p_k_old[i]
						 -I*(ne_k_old[i] + nh_k_old[i] - 1.0)*(E_tmp + electric_field_renormalized[i]);
		}

	}
}

/* Evaluate RHS of SBE: Electrons
 * Fills the array ne_k_new
 * */
void Device::sbe_RHS_ne(double DT, double *ne_k_new, std::complex<double> *p_k_old, double *ne_k_old, double *nh_k_old, std::complex<double> E_tmp)
{
	#ifndef COMPILE_REFLECTION_CALCULATION
	if (device_density > 1.0e16)
	{
		for(unsigned i = 0; i < number_K_points; i++)
		{
			ne_k_new[i] = -2.0*imag((E_tmp + electric_field_renormalized[i])*conj(p_k_old[i]));

			// Hole filling
		#ifdef USE_HOLE_FILLING
			#if defined(USE_ISAK_HOLE_FILLING) || defined(USE_ISAK_HOLE_FILLING_TABLE)
				ne_k_new[i] += SBE_TIME_SCALE*carrier_scattering_rates_e_total[i] + SBE_TIME_SCALE*carrier_scattering_rates_phonon_e_total[i];		
			#else
				ne_k_new[i] += carrier_scattering_rate_approximation_e[i];
			#endif
		#endif
			
			#ifdef SBE_USE_SPONTAN_EMISSIONS
				ne_k_new[i] += device_spont_emission_wk[i]*ne_k_old[i]*nh_k_old[i];
			#endif
			
			#ifdef SBE_USE_RESONANT_PUMP_MODEL
				ne_k_new[i] += (ne_k_old[i] + nh_k_old[i] -1.0)*device_pump_wk[i];
			#else			
				// Long recovery time
				ne_k_new[i] += -(ne_k_old[i]-fe_k[i])*SBE_OCC_PUMP;
			#endif
		}
	} else {
		
		for(unsigned i = 0; i < number_K_points; i++)
		{
			ne_k_new[i] = -2.0*imag((E_tmp + electric_field_renormalized[i])*conj(p_k_old[i]));
//		#if defined(USE_ISAK_HOLE_FILLING) || defined(USE_ISAK_HOLE_FILLING_TABLE)
//			ne_k_new[i] += SBE_TIME_SCALE*carrier_scattering_rates_e_total[i];
//		#elif defined(USE_HOLE_FILLING)
			ne_k_new[i] += carrier_scattering_rate_approximation_e[i];
//		#endif
			ne_k_new[i] += -(ne_k_old[i]-fe_k[i])*SBE_OCC_PUMP;
			
//			#ifdef SBE_USE_SPONTAN_EMISSIONS
//				ne_k_new[i] += device_spont_emission_wk[i]*ne_k_old[i]*nh_k_old[i];
//			#endif
		}
	}
	#endif
}

/* Evaluate RHS of SBE: Holes
 * Fills the array nh_k_new
 * */
void Device::sbe_RHS_nh(double DT, double *nh_k_new, std::complex<double> *p_k_old, double *ne_k_old, double *nh_k_old, std::complex<double> E_tmp)
{
	#ifndef COMPILE_REFLECTION_CALCULATION
	if (device_density > 1.0e16)
	{
		for(unsigned i = 0; i < number_K_points; i++)
		{
			nh_k_new[i] = -2.0*imag((E_tmp + electric_field_renormalized[i])*conj(p_k_old[i]));

			// Hole filling
		#ifdef USE_HOLE_FILLING
			#if defined(USE_ISAK_HOLE_FILLING) || defined(USE_ISAK_HOLE_FILLING_TABLE)
				nh_k_new[i] += SBE_TIME_SCALE*carrier_scattering_rates_h_total[i] + SBE_TIME_SCALE*carrier_scattering_rates_phonon_h_total[i];
			#else
				nh_k_new[i] += carrier_scattering_rate_approximation_h[i];
			#endif
		#endif

			#ifdef SBE_USE_SPONTAN_EMISSIONS
				nh_k_new[i] += device_spont_emission_wk[i]*ne_k_old[i]*nh_k_old[i];
			#endif
			
			#ifdef SBE_USE_RESONANT_PUMP_MODEL
				nh_k_new[i] += (ne_k_old[i] + nh_k_old[i] -1.0)*device_pump_wk[i];
			#else
				// Long recovery time
				nh_k_new[i] += -(nh_k_old[i]-fh_k[i])*SBE_OCC_PUMP;
			#endif
		}
	} else {
		for(unsigned i = 0; i < number_K_points; i++)
		{
			nh_k_new[i] = -2.0*imag((E_tmp + electric_field_renormalized[i])*conj(p_k_old[i]));
//		#if defined(USE_ISAK_HOLE_FILLING) || defined(USE_ISAK_HOLE_FILLING_TABLE)
//			nh_k_new[i] += SBE_TIME_SCALE*carrier_scattering_rates_h_total[i];
//		#elif defined(USE_HOLE_FILLING)
			nh_k_new[i] += carrier_scattering_rate_approximation_h[i];
//		#endif
			nh_k_new[i] += -(nh_k_old[i]-fh_k[i])*SBE_OCC_PUMP;
			
//			#ifdef SBE_USE_SPONTAN_EMISSIONS
//				nh_k_new[i] += device_spont_emission_wk[i]*ne_k_old[i]*nh_k_old[i];
//			#endif
		}
	}
	#endif

}

/* Evaluate RHS of SBE: Polarization
 * Fills the array p_k_new. Used for debugging purposes
 * Does NOT use the renormalized values
 * */
void Device::sbe_RHS_pk_debug(double t, double DT, std::complex<double> *p_k_new, std::complex<double> *p_k_old, double *ne_k_old, double *nh_k_old, std::complex<double> E_tmp)
{
	double scale_factor = SBE_TIME_SCALE/hbar;

	if (device_density > 1.0e16)
	{
		for(unsigned i = 0; i < number_K_points; i++)
		{
			// Scaled equations
			p_k_new[i] =    -I*(sumE_k[i])*p_k_old[i]
						 -SBE_POL_DEPH*p_k_old[i]
						 -I*(ne_k_old[i] + nh_k_old[i] - 1.0)*(E_tmp);
						 
			#ifdef SBE_USE_SPONTAN_EMISSIONS
				p_k_new[i] += -device_spont_emission_phasor*device_spont_emission_wk[i]*ne_k_old[i]*nh_k_old[i];
			#endif
		}
	} else {
		for(unsigned i = 0; i < number_K_points; i++)
		{
			// Scaled equations
			p_k_new[i] = -I*(sumE_k[i])*p_k_old[i]
						 -SBE_POL_DEPH*p_k_old[i]
						 -I*(ne_k_old[i] + nh_k_old[i] - 1.0)*(E_tmp);
		}

	}
}

/* Evaluate RHS of SBE: Polarization
 * Fills the array p_k_new
 * Does NOT use the renormalized values or field or scattering
 * */
void Device::sbe_RHS_pk_FieldFree(double DT, std::complex<double> *p_k_new, std::complex<double> *p_k_old, double *ne_k_old, double *nh_k_old)
{
	if (device_density > 1.0e16)
	{
		for(unsigned i = 0; i < number_K_points; i++)
		{
			// Scaled equations
			p_k_new[i] = -I*(sumE_k[i])*p_k_old[i]
					-SBE_POL_DEPH*p_k_old[i];
						 
			#ifdef SBE_USE_SPONTAN_EMISSIONS
				p_k_new[i] += -device_spont_emission_phasor*device_spont_emission_wk[i]*ne_k_old[i]*nh_k_old[i];
			#endif
		}
	} else {
		for(unsigned i = 0; i < number_K_points; i++)
		{
			// Scaled equations
			p_k_new[i] = -I*(sumE_k[i])*p_k_old[i]
						 -SBE_POL_DEPH*p_k_old[i];
		}

	}
}


/* Evaluate RHS of SBE: Polarization
 * Fills the array p_k_new
 * Does NOT use the renormalized values or scattering
 * */
void Device::sbe_RHS_pk_SmallAmp(double DT, std::complex<double> *p_k_new, std::complex<double> *p_k_old, double *ne_k_old, double *nh_k_old, std::complex<double> E_tmp)
{
	double scale_factor = SBE_TIME_SCALE/hbar;
	
	if (device_density > 1.0e16)
	{
		for(unsigned i = 0; i < number_K_points; i++)
		{
			// Scaled equations
			p_k_new[i] = -I*(sumE_k[i])*p_k_old[i]
					-SBE_POL_DEPH*p_k_old[i]
					-I*(ne_k_old[i] + nh_k_old[i] - 1.0)*(E_tmp);
						 
			#ifdef SBE_USE_SPONTAN_EMISSIONS
				p_k_new[i] += -device_spont_emission_phasor*device_spont_emission_wk[i]*ne_k_old[i]*nh_k_old[i];
			#endif
		}
	} else {
		for(unsigned i = 0; i < number_K_points; i++)
		{
			// Scaled equations
			p_k_new[i] = -I*(sumE_k[i])*p_k_old[i]
						 -SBE_POL_DEPH*p_k_old[i]
						 -I*(ne_k_old[i] + nh_k_old[i] - 1.0)*(E_tmp);
		}

	}
}

/* Evaluate RHS of SBE: Electrons
 * Fills the array ne_k_new
 * Does NOT use the renormalized values or field or hole filling
 * */
void Device::sbe_RHS_ne_FieldFree(double DT, double *ne_k_new, std::complex<double> *p_k_old, double *ne_k_old, double *nh_k_old)
{
	#ifndef COMPILE_REFLECTION_CALCULATION
	if (device_density > 1.0e16)
	{
		for(unsigned i = 0; i < number_K_points; i++)
		{
			ne_k_new[i] = 0.0;
			
			#ifdef SBE_USE_SPONTAN_EMISSIONS
				ne_k_new[i] += device_spont_emission_wk[i]*ne_k_old[i]*nh_k_old[i];
			#endif
			
			#ifdef SBE_USE_RESONANT_PUMP_MODEL
				ne_k_new[i] += (ne_k_old[i] + nh_k_old[i] -1.0)*device_pump_wk[i];
			#else			
				// Long recovery time
				ne_k_new[i] += -(ne_k_old[i]-fe_k[i])*SBE_OCC_PUMP;
			#endif
		}
	} else {
		
		for(unsigned i = 0; i < number_K_points; i++)
		{
			ne_k_new[i] = -(ne_k_old[i]-fe_k[i])*SBE_OCC_PUMP;
		}
		
	}
	#endif
}


/* Evaluate RHS of SBE: Electrons
 * Fills the array ne_k_new
 * Does NOT use the renormalized values or hole filling
 * */
void Device::sbe_RHS_ne_SmallAmp(double DT, double *ne_k_new, std::complex<double> *p_k_old, double *ne_k_old, double *nh_k_old, std::complex<double> E_tmp)
{
	#ifndef COMPILE_REFLECTION_CALCULATION
	if (device_density > 1.0e16)
	{
		for(unsigned i = 0; i < number_K_points; i++)
		{
			ne_k_new[i] = -2.0*imag((E_tmp)*conj(p_k_old[i]));
			
			#ifdef SBE_USE_SPONTAN_EMISSIONS
				ne_k_new[i] += device_spont_emission_wk[i]*ne_k_old[i]*nh_k_old[i];
			#endif
			
			#ifdef SBE_USE_RESONANT_PUMP_MODEL
				ne_k_new[i] += (ne_k_old[i] + nh_k_old[i] -1.0)*device_pump_wk[i];
			#else			
				// Long recovery time
				ne_k_new[i] += -(ne_k_old[i]-fe_k[i])*SBE_OCC_PUMP;
			#endif
		}
	} else {
		
		for(unsigned i = 0; i < number_K_points; i++)
		{
			ne_k_new[i] = -2.0*imag((E_tmp)*conj(p_k_old[i]));
			ne_k_new[i] += -(ne_k_old[i]-fe_k[i])*SBE_OCC_PUMP;
		}
		
	}
	#endif
}

/* Evaluate RHS of SBE: Holes
 * Fills the array nh_k_new
 * Does NOT use the renormalized values or field or hole filling
 * */
void Device::sbe_RHS_nh_FieldFree(double DT, double *nh_k_new, std::complex<double> *p_k_old, double *ne_k_old, double *nh_k_old)
{
	#ifndef COMPILE_REFLECTION_CALCULATION
	if (device_density > 1.0e16)
	{
		for(unsigned i = 0; i < number_K_points; i++)
		{
			nh_k_new[i] = 0.0;

			#ifdef SBE_USE_SPONTAN_EMISSIONS
				nh_k_new[i] += device_spont_emission_wk[i]*ne_k_old[i]*nh_k_old[i];
			#endif
			
			#ifdef SBE_USE_RESONANT_PUMP_MODEL
				nh_k_new[i] += (ne_k_old[i] + nh_k_old[i] -1.0)*device_pump_wk[i];
			#else
				// Long recovery time
				nh_k_new[i] += -(nh_k_old[i]-fh_k[i])*SBE_OCC_PUMP;
			#endif
		}
	} else {
		for(unsigned i = 0; i < number_K_points; i++)
		{
			nh_k_new[i] = -(nh_k_old[i]-fh_k[i])*SBE_OCC_PUMP;
		}
	}
	#endif
}

/* Evaluate RHS of SBE: Holes
 * Fills the array nh_k_new
 * Does NOT use the renormalized values or hole filling
 * */
void Device::sbe_RHS_nh_SmallAmp(double DT, double *nh_k_new, std::complex<double> *p_k_old, double *ne_k_old, double *nh_k_old, std::complex<double> E_tmp)
{
	#ifndef COMPILE_REFLECTION_CALCULATION
	if (device_density > 1.0e16)
	{
		for(unsigned i = 0; i < number_K_points; i++)
		{
			nh_k_new[i] = -2.0*imag((E_tmp)*conj(p_k_old[i]));

			#ifdef SBE_USE_SPONTAN_EMISSIONS
				nh_k_new[i] += device_spont_emission_wk[i]*ne_k_old[i]*nh_k_old[i];
			#endif
			
			#ifdef SBE_USE_RESONANT_PUMP_MODEL
				nh_k_new[i] += (ne_k_old[i] + nh_k_old[i] -1.0)*device_pump_wk[i];
			#else
				// Long recovery time
				nh_k_new[i] += -(nh_k_old[i]-fh_k[i])*SBE_OCC_PUMP;
			#endif
		}
	} else {
		for(unsigned i = 0; i < number_K_points; i++)
		{
			nh_k_new[i] = -2.0*imag((E_tmp)*conj(p_k_old[i]));
			nh_k_new[i] += -(nh_k_old[i]-fh_k[i])*SBE_OCC_PUMP;
		}
	}
	#endif
}




void Device::updateCoulombMatrix()
{
/*
	FILE *fid2 = fopen("matrix_coulomb_screening.dat","w+");
	for(int k = 0; k < number_K_points; k++)
	{
		fprintf(fid2,"%.3e ",coulomb_potential_epsilon_inv[k]);
	}
	fclose(fid2);
*/
	
	for(int k = 0; k < number_K_points; k++)
	{
		for(int q = 0; q < number_K_points; q++)
		{
			// Integrate as normal, no singularity
			double sum_ee_kq = 0.0;
			double sum_hh_kq = 0.0;
			double sum_eh_kq = 0.0;

			for(int phi_kq = 0; phi_kq < number_Eta_points; phi_kq++)
			{
				double k_q = sqrt(K[k]*K[k] + K[q]*K[q] - 2.0*K[k]*K[q]*cos(Eta[phi_kq]));

				if ((k_q >= K[0])&&(k_q <= K[number_K_points-1]))
				{
					double Cee  = misc_interpolate_K_array_parabolic(k_q, K, coulomb_potential_confinement_ee, number_K_points);
					double Chh  = misc_interpolate_K_array_parabolic(k_q, K, coulomb_potential_confinement_hh, number_K_points);
					double Ceh  = misc_interpolate_K_array_parabolic(k_q, K, coulomb_potential_confinement_eh, number_K_points);
					
					double eps_inv = 1.0; // Simple approximation for speed
		//			double eps_inv  = misc_interpolate_K_array_parabolic(k_q, K, coulomb_potential_epsilon_inv, number_K_points);

					sum_ee_kq += Cee*2.0*dEta[phi_kq]/(k_q*eps_inv); // Multiply by because of half angular grid
					sum_hh_kq += Chh*2.0*dEta[phi_kq]/(k_q*eps_inv);
					sum_eh_kq += Ceh*2.0*dEta[phi_kq]/(k_q*eps_inv);
				}
			}



			coulomb_matrix_ee[k][q] = sum_ee_kq;
			coulomb_matrix_hh[k][q] = sum_hh_kq;
			coulomb_matrix_eh[k][q] = sum_eh_kq;
				

			// Constants multiply by kdK from the integrals
			coulomb_matrix_ee[k][q] *= KdK[q]*e*e/(8.0*Pi*Pi*eps*eps0*a0);
			coulomb_matrix_hh[k][q] *= KdK[q]*e*e/(8.0*Pi*Pi*eps*eps0*a0);
			coulomb_matrix_eh[k][q] *= KdK[q]*e*e/(8.0*Pi*Pi*eps*eps0*a0);
			
			// Time scale
			// SBE_TIME_SCALE is from the SBE integration, hbar from the renomarlizations
			coulomb_matrix_ee[k][q] *= SBE_TIME_SCALE/hbar;
			coulomb_matrix_hh[k][q] *= SBE_TIME_SCALE/hbar;
			coulomb_matrix_eh[k][q] *= SBE_TIME_SCALE/hbar;
		}
		
		
		coulomb_matrix_ee[k][k] = 0.0;
		coulomb_matrix_hh[k][k] = 0.0;
		coulomb_matrix_eh[k][k] = 0.0;
	}
/*
	FILE *fid = fopen("matrix_coulomb_ee.dat","w+");
	for(int k = 0; k < number_K_points; k++)
	{
		for(int q = 0; q < number_K_points; q++)
		{
			fprintf(fid,"%.3e ",coulomb_matrix_ee[k][q]);
		}
		fprintf(fid, "\n");
	}
	fclose(fid);
	
	fid = fopen("matrix_coulomb_hh.dat","w+");
	for(int k = 0; k < number_K_points; k++)
	{
		for(int q = 0; q < number_K_points; q++)
		{
			fprintf(fid,"%.3e ",coulomb_matrix_hh[k][q]);
		}
		fprintf(fid, "\n");
	}
	fclose(fid);


	fid = fopen("matrix_coulomb_eh.dat","w+");
	for(int k = 0; k < number_K_points; k++)
	{
		for(int q = 0; q < number_K_points; q++)
		{
			fprintf(fid,"%.3e ",coulomb_matrix_eh[k][q]);
		}
		fprintf(fid, "\n");
	}
	fclose(fid);
*/
}

void Device::updateCoulombMatrixScreened()
{
/*
	FILE *fid2 = fopen("matrix_coulomb_screening.dat","w+");
	for(int k = 0; k < number_K_points; k++)
	{
		fprintf(fid2,"%.3e ",coulomb_potential_epsilon_inv[k]);
	}
	fclose(fid2);
*/
	for(int k = 0; k < number_K_points; k++)
	{
		for(int q = 0; q < number_K_points; q++)
		{
			// Integrate as normal, no singularity
			double sum_ee_kq = 0.0;
			double sum_hh_kq = 0.0;
			double sum_eh_kq = 0.0;

			for(int phi_kq = 0; phi_kq < number_Eta_points; phi_kq++)
			{
				double k_q = sqrt(K[k]*K[k] + K[q]*K[q] - 2.0*K[k]*K[q]*cos(Eta[phi_kq]));

				if ((k_q >= K[0])&&(k_q <= K[number_K_points-1]))
				{
					double Cee  = misc_interpolate_K_array_parabolic(k_q, K, coulomb_potential_confinement_ee, number_K_points);
					double Chh  = misc_interpolate_K_array_parabolic(k_q, K, coulomb_potential_confinement_hh, number_K_points);
					double Ceh  = misc_interpolate_K_array_parabolic(k_q, K, coulomb_potential_confinement_eh, number_K_points);
					
					//double eps_inv = 1.0; // Simple approximation for speed
					double eps_inv  = misc_interpolate_K_array(k_q, K, coulomb_potential_epsilon_inv, number_K_points,1.0);

					sum_ee_kq += Cee*2.0*dEta[phi_kq]/(k_q*eps_inv); // Multiply by because of half angular grid
					sum_hh_kq += Chh*2.0*dEta[phi_kq]/(k_q*eps_inv);
					sum_eh_kq += Ceh*2.0*dEta[phi_kq]/(k_q*eps_inv);
				}
			}



			coulomb_matrix_ee[k][q] = sum_ee_kq;
			coulomb_matrix_hh[k][q] = sum_hh_kq;
			coulomb_matrix_eh[k][q] = sum_eh_kq;
				

			// Constants multiply by kdK from the integrals
			coulomb_matrix_ee[k][q] *= KdK[q]*e*e/(8.0*Pi*Pi*eps*eps0*a0);
			coulomb_matrix_hh[k][q] *= KdK[q]*e*e/(8.0*Pi*Pi*eps*eps0*a0);
			coulomb_matrix_eh[k][q] *= KdK[q]*e*e/(8.0*Pi*Pi*eps*eps0*a0);
			
			// Time scale
			// SBE_TIME_SCALE is from the SBE integration, hbar from the renomarlizations
			coulomb_matrix_ee[k][q] *= SBE_TIME_SCALE/hbar;
			coulomb_matrix_hh[k][q] *= SBE_TIME_SCALE/hbar;
			coulomb_matrix_eh[k][q] *= SBE_TIME_SCALE/hbar;
		}
		
		
		coulomb_matrix_ee[k][k] = 0.0;
		coulomb_matrix_hh[k][k] = 0.0;
		coulomb_matrix_eh[k][k] = 0.0;
	}
/*
	FILE *fid = fopen("matrix_coulomb_screened_ee.dat","w+");
	for(int k = 0; k < number_K_points; k++)
	{
		for(int q = 0; q < number_K_points; q++)
		{
			fprintf(fid,"%.3e ",coulomb_matrix_ee[k][q]);
		}
		fprintf(fid, "\n");
	}
	fclose(fid);
	
	fid = fopen("matrix_coulomb_screened_hh.dat","w+");
	for(int k = 0; k < number_K_points; k++)
	{
		for(int q = 0; q < number_K_points; q++)
		{
			fprintf(fid,"%.3e ",coulomb_matrix_hh[k][q]);
		}
		fprintf(fid, "\n");
	}
	fclose(fid);


	fid = fopen("matrix_coulomb_screened_eh.dat","w+");
	for(int k = 0; k < number_K_points; k++)
	{
		for(int q = 0; q < number_K_points; q++)
		{
			fprintf(fid,"%.3e ",coulomb_matrix_eh[k][q]);
		}
		fprintf(fid, "\n");
	}
	fclose(fid);


	cout << "output screened Vc done" << endl;
	exit(-1);
*/
}

/* Update screening of reduced Coulomb matrix based on carrier numbers
 * */
void Device::updateScreening_reducedCoulomb_StaticPlasmonApproximation(double *ne_k_old, double *nh_k_old)
{
	// Fill coulomb_potential_epsilon_inv[]
	updateScreening_StaticPlasmonApproximation(ne_k_old, nh_k_old);
	
	double red_jh=((double)coulomb_matrix_red_index_j);
	double inv_red_ih=1.0/((double)coulomb_matrix_red_index_i);
	int n = 0;

	// Update reduced 2D Coulomb matrix
	for(int i=0; i<number_K_points; i+=coulomb_matrix_red_index_i ){
		for( int iiq = i+1; iiq < number_K_points; iiq+=coulomb_matrix_red_index_j ){
			double eps_iiq = red_jh/coulomb_potential_epsilon_inv[iiq-i];
			coulomb_matrix_ee_red[n] = coulomb_matrix_ee[i][iiq]*eps_iiq;
			coulomb_matrix_hh_red[n] = coulomb_matrix_hh[i][iiq]*eps_iiq;
			coulomb_matrix_eh_red[n] = coulomb_matrix_eh[i][iiq]*eps_iiq;
			n++;
		}
		for( int iiq = i-1; iiq >=0 ; iiq-=coulomb_matrix_red_index_j ){
			double eps_iiq = red_jh/coulomb_potential_epsilon_inv[i-iiq];
			coulomb_matrix_ee_red[n] = coulomb_matrix_ee[i][iiq]*eps_iiq;
			coulomb_matrix_hh_red[n] = coulomb_matrix_hh[i][iiq]*eps_iiq;
			coulomb_matrix_eh_red[n] = coulomb_matrix_eh[i][iiq]*eps_iiq;
			n++;
		}
	}
}

// Update screening ONLY using static plasmon approx
// Fill array coulomb_potential_epsilon_inv[]
void Device::updateScreening_StaticPlasmonApproximation(double *ne_k_old, double *nh_k_old)
{
	//=============================================================================
	//	calculate screening epsilon_q using static plasmon pole approximation
	//	const_c: numerical factor; fitted to full screening calculation result 
	//=============================================================================

	long int n = 0l;
	double const_c = (8.0e-15l);

	//calculate density; search k with maximum occupation; standard: k=0; due to deviations from Fermi: search
	double h1= 0.0;
	double density_eh= 0.0;
	int q_max;
	for(int i=0; i<number_K_points; i++ ){

		density_eh += ne_k_old[i]*KdK[i];
		double h2 = (ne_k_old[i]+nh_k_old[i]);
		if( h2>h1 ){
			q_max = i;
			h1 = h2;
		}
	}

	density_eh /= (Pi*a0*a0);

	//Bohr radius
	double a_bohr = 4.0*Pi*hbar*hbar*eps*eps0/device_mr/e/e;

	//screening length
	double kappa = 2.0/a_bohr*(device_me/device_mr*ne_k_old[q_max]+device_mh/device_mr*nh_k_old[q_max]);
	
	/*
	//double density_eh = ne_k_old[0]*K[0]*dK[0];
	// Calculate kappa from 
	double kappa_exact = 0.0;
	{
		double fact_e_q = 0.0;
		double fact_h_q = 0.0;
		
		for(int k = 1; k < number_K_points; k++)
		{
			// Use density below
			//density_eh += ne_k_old[k]*K[k]*dK[k];
			
			
			// Electrons
			//double dfe_k    = a0*interpolate_cubic_spline_ne->evaluate_df(K[k]); // sum from i=0
			double dEe_k  = hbar*hbar*K[k]/(device_me*a0);
			double dfe_k = a0*(ne_k_old[k]-ne_k_old[k-1])/(K[k]-K[k-1]); // sum from i=1
			//double dEe_k = a0*(Ee_k[k]-Ee_k[k-1])/(K[k]-K[k-1]);
			double fact_e_q_tmp = dfe_k/dEe_k;
			
			// Holes
			//double dfh_k    = a0*interpolate_cubic_spline_nh->evaluate_df(K[k]); // sum from i=0
			double dEh_k  = hbar*hbar*K[k]/(device_mh*a0);
			double dfh_k = a0*(nh_k_old[k]-nh_k_old[k-1])/(K[k]-K[k-1]); // sum from i=1
			//double dEh_k = a0*(Eh_k[k]-Eh_k[k-1])/(K[k]-K[k-1]);
			double fact_h_q_tmp = dfh_k/dEh_k;
			
			fact_e_q += K[k]*dK[k]*fact_e_q_tmp;
			fact_h_q += K[k]*dK[k]*fact_h_q_tmp;
		}
		
		double integration_volume = 1.0/(Pi*a0*a0);
		fact_e_q *= integration_volume;
		fact_h_q *= integration_volume;
		
		kappa_exact = -(fact_e_q + fact_h_q)*(e*e/(2.0*eps*eps0));
	}
	//density_eh /= (Pi*a0*a0);
	
	cout << "kappa(est)	= " << kappa << endl;
	cout << "kappa(exact)	= " << kappa_exact << endl;
	cout << "ratio	= " << kappa/kappa_exact << endl;
	kappa = kappa_exact;
	*/

	//calculate prefactor for k-dependent plasma frequency
	double omega_pl_fact = density_eh*e*e*device_length_qw/eps/eps0/device_mr*0.5;
	double omega_q_fact = const_c*kappa*kappa*kappa/omega_pl_fact;

	//calculate epsilon_q(k)=(1-omega_pl^2/omega_q^2)
	for(int i=0; i < number_K_points; i++ )
	{
		double dq = K[i]/a0;

		//omega_pl^2
		double omega_pl = omega_pl_fact*dq;
		double dq_kappa = dq/kappa;

		//omega_q^2
		double omega_q = 1.0 + dq_kappa + omega_q_fact*dq_kappa*dq_kappa*dq_kappa;
		coulomb_potential_epsilon_inv[i] = omega_q/(omega_q-1.0);
	}

	// Calculate screened 1D Coulomb matrix
	for(int k = 0; k < number_K_points; k++)
	{
		// Factor a0*e*e/(2.0*eps0*eps) has been divided out
		double eps_inv = 1.0/(K[k]*coulomb_potential_epsilon_inv[k]);
		coulomb_potential_normalized_ee[k] = coulomb_potential_confinement_ee[k]*eps_inv; // screened matrix Vq
		coulomb_potential_normalized_hh[k] = coulomb_potential_confinement_hh[k]*eps_inv; // screened matrix Vq
		coulomb_potential_normalized_eh[k] = coulomb_potential_confinement_eh[k]*eps_inv; // screened matrix Vq
	}
}

/* Update screening of reduced Coulomb matrix based on carrier numbers
 * uses full screening calculation
 * Requires that the full screening method has been initializied
 * */
void Device::updateScreening_reducedCoulomb_full(double *ne_k_old, double *nh_k_old)
{
	// Calculate full screening based on carrier numbers
	// This will also fill 
	// Determine if there is a significant change in the CARRIRES


	carrier_scattering_screening_full_fast(ne_k_old, nh_k_old);

	double red_jh=((double)coulomb_matrix_red_index_j);
	double inv_red_ih=1.0/((double)coulomb_matrix_red_index_i);
	int n = 0;

	// Update reduced 2D Coulomb matrix
	for(int i=0; i<number_K_points; i+=coulomb_matrix_red_index_i ){
		for( int iiq = i+1; iiq < number_K_points; iiq+=coulomb_matrix_red_index_j ){
			double eps_inv = red_jh/coulomb_potential_epsilon_inv[iiq-i];
			coulomb_matrix_ee_red[n] = coulomb_matrix_ee[i][iiq]*eps_inv;
			coulomb_matrix_hh_red[n] = coulomb_matrix_hh[i][iiq]*eps_inv;
			coulomb_matrix_eh_red[n] = coulomb_matrix_eh[i][iiq]*eps_inv;
			n++;
		}
		for( int iiq = i-1; iiq >=0 ; iiq-=coulomb_matrix_red_index_j ){
			double eps_inv = red_jh/coulomb_potential_epsilon_inv[i-iiq];
			coulomb_matrix_ee_red[n] = coulomb_matrix_ee[i][iiq]*eps_inv;
			coulomb_matrix_hh_red[n] = coulomb_matrix_hh[i][iiq]*eps_inv;
			coulomb_matrix_eh_red[n] = coulomb_matrix_eh[i][iiq]*eps_inv;
			n++;
		}
	}
}


/* LEGACY: Renormalize energy and field without screening or reduced Coulomb matrix
 *	Renormalize energy and field
 *	 - Divide by hbar in the coulomb matrix to get correct units
 *	 - energy_renormalized in units [1/s]
 *	 - electric_field_renormalized in units [hbar]
 * */
void Device::renormalize_Unscreened(std::complex<double> *p_k_new, double *ne_k_new, double *nh_k_new)
{
	for (unsigned i = 0; i < number_K_points; i++)
	{
		electric_field_renormalized[i] = 0;
		energy_renormalized_e[i] = 0;
		energy_renormalized_h[i] = 0;
		
		for(unsigned j = 0; j < number_K_points; j++)
		{
			electric_field_renormalized[i] += coulomb_matrix_eh[i][j]*p_k_new[j];
			energy_renormalized_e[i] += coulomb_matrix_ee[i][j]*ne_k_new[j];
			energy_renormalized_h[i] += coulomb_matrix_hh[i][j]*nh_k_new[j];
		}
		
//		electric_field_renormalized[i] /= hbar;
//		energy_renormalized[i] /= hbar;
	}
}
/* Zero energy and field renormalizations */
void Device::zero_reducedCoulomb()
{
	for (unsigned i = 0; i < number_K_points; i++)
	{
		electric_field_renormalized[i]  = 0.0;
		energy_renormalized_e[i]	= 0.0;
		energy_renormalized_h[i]	= 0.0;
	}
}

/* Renormalize energy and field with screening and reduced Coulomb matrix
 *	Renormalize energy and field
 *	 - Divide by hbar in the coulomb matrix to get correct units
 *	 - energy_renormalized in units [1/s]
 *	 - electric_field_renormalized in units [hbar]
 * */
void Device::renormalize_reducedCoulomb(std::complex<double> *p_k_old, double *ne_k_old, double *nh_k_old)
{
//	#ifdef USE_DEVICE_TIMERS
//	MainStat->start(getName(),"renormalize");
//	#endif
	long int n = 0l;
	double red_jh=((double)coulomb_matrix_red_index_j);
	double inv_red_ih=1.0/((double)coulomb_matrix_red_index_i);

	bool calc = false;

	if (device_density > 1.0e16)
	{

		// Check if distributions have changed
		double err_e_max = 0.0;
		double err_h_max = 0.0;
		double err_e,err_h;
		if (CARRIER_DEVIATION_TOL > 0)
		{
			for(int i = 0; i < number_K_points; i++)
			{
				if (fabs(renormalize_prev_ne[i] + ne_k_old[i]) > 0.0)
				{
					err_e = fabs(renormalize_prev_ne[i] - ne_k_old[i])/fabs(renormalize_prev_ne[i] + ne_k_old[i]);
					if (err_e > err_e_max)
					{
						err_e_max = err_e;
					}
				}

				if (fabs(renormalize_prev_nh[i] + nh_k_old[i]) > 0.0)
				{
					err_h = fabs(renormalize_prev_nh[i] - nh_k_old[i])/fabs(renormalize_prev_nh[i] + nh_k_old[i]);
					if (err_h > err_h_max)
					{
						err_h_max = err_h;
					}
				}


			}

			if ((err_e_max > CARRIER_DEVIATION_TOL)||(err_h_max > CARRIER_DEVIATION_TOL))
			{
				calc = true;
			}
		} else {
			calc = true;
		}
	} else {
		calc = true;
	}

	if (calc)
	{
		// Store previous carriers for error evaluation
		memcpy(renormalize_prev_ne, ne_k_old, number_K_points*sizeof(double));
		memcpy(renormalize_prev_nh, nh_k_old, number_K_points*sizeof(double));

	#ifdef SBE_UPDATE_SCREENING
		// Update screening at each timestep, otherwise the screening is only calculated during initialization
		#ifdef SBE_USE_FULL_SCREENING

		// Use full screening
		updateScreening_reducedCoulomb_full(ne_k_old, nh_k_old);

		#else
		
		// Basic screening: static plasmon approximation
		updateScreening_reducedCoulomb_StaticPlasmonApproximation(ne_k_old, nh_k_old);

		#endif
	#endif
		
		
		// Fill arrays energy_renormalized[]
	#ifdef SBE_UPDATE_BG_RENORM
		renormalize_reducedCoulomb_Energy(ne_k_old, nh_k_old);
	#endif
	
	}

	// Fill arrays electric_field_renormalized[]
	renormalize_reducedCoulomb_Field(p_k_old);
	
//	#ifdef USE_DEVICE_TIMERS
//	MainStat->stop(getName(),"renormalize");
//	#endif
}

/* Fill array energy_renormalized[] 
 * Requires the matrix coulomb_matrix_red be computed ahead of time
 * */
void Device::renormalize_reducedCoulomb_Energy(double *ne_k_old, double *nh_k_old)
{
	long int n = 0l;
	double red_jh=((double)coulomb_matrix_red_index_j);
	double inv_red_ih=1.0/((double)coulomb_matrix_red_index_i);

	// Calculate renormalization
	n=0l;
	for(int i=0; i<number_K_points; i+=coulomb_matrix_red_index_i ){
		energy_renormalized_e[i] = 0.0;
		energy_renormalized_h[i] = 0.0;
		
		for( int iiq = i+1; iiq < number_K_points ; iiq+=coulomb_matrix_red_index_j )
		{
			energy_renormalized_e[i] +=  coulomb_matrix_ee_red[n] * ne_k_old[iiq];
			energy_renormalized_h[i] +=  coulomb_matrix_hh_red[n] * nh_k_old[iiq];
			n++;
		}
		for( int iiq = i-1; iiq >=0 ; iiq-=coulomb_matrix_red_index_j )
		{
			energy_renormalized_e[i] +=  coulomb_matrix_ee_red[n] * ne_k_old[iiq];
			energy_renormalized_h[i] +=  coulomb_matrix_hh_red[n] * nh_k_old[iiq];
			n++;
		}
	}
	if (coulomb_matrix_red_index_i > 1)
	{
		int i = 0;
		for(i=0; i<(number_K_points-coulomb_matrix_red_index_i); i+=coulomb_matrix_red_index_i )
		{
/*
			// Set parabolic coefficients
			int i0 = i;
			int i1 = i+coulomb_matrix_red_index_i;
			int i2 = i+2*coulomb_matrix_red_index_i;
			while (i2 >= number_K_points)
			{
				// On left edge
				i0 -= coulomb_matrix_red_index_i;
				i1 -= coulomb_matrix_red_index_i;
				i2 -= coulomb_matrix_red_index_i;
			}
			
			double t0 = K[i0];
			double t1 = K[i1];
			double t2 = K[i2];
			double delta0 = (t0-t1)*(t0-t2);
			double delta1 = (t1-t0)*(t1-t2);
			double delta2 = (t2-t0)*(t2-t1);
*/


			for( int j=1; j<coulomb_matrix_red_index_i; j++ )
			{
				// Linear
				double dL = (K[i+j]-K[i])/(K[i+coulomb_matrix_red_index_i]-K[i]);
				energy_renormalized_e[i+j] = energy_renormalized_e[i]+dL*(energy_renormalized_e[i+coulomb_matrix_red_index_i]-energy_renormalized_e[i]);
				energy_renormalized_h[i+j] = energy_renormalized_h[i]+dL*(energy_renormalized_h[i+coulomb_matrix_red_index_i]-energy_renormalized_h[i]);
/*
				// Parabolic
				double L0 = (K[i+j]-t1)*(K[i+j]-t2)/delta0;
				double L1 = (K[i+j]-t0)*(K[i+j]-t2)/delta1;
				double L2 = (K[i+j]-t0)*(K[i+j]-t1)/delta2;
				energy_renormalized_e[i+j] = energy_renormalized_e[i0]*L0 + energy_renormalized_e[i1]*L1 + energy_renormalized_e[i2]*L2;
				energy_renormalized_h[i+j] = energy_renormalized_h[i0]*L0 + energy_renormalized_h[i1]*L1 + energy_renormalized_h[i2]*L2;
*/
			}
		}

		// Interpolate the last few steps of the grid
		for(int j=(i-coulomb_matrix_red_index_i+1); j<number_K_points; j++)
		{
			double dL = (K[j]-K[i-coulomb_matrix_red_index_i])/(K[i]-K[i-coulomb_matrix_red_index_i]);
			energy_renormalized_e[j] = energy_renormalized_e[i-coulomb_matrix_red_index_i]+dL*(energy_renormalized_e[i]-energy_renormalized_e[i-coulomb_matrix_red_index_i]);
			energy_renormalized_h[j] = energy_renormalized_h[i-coulomb_matrix_red_index_i]+dL*(energy_renormalized_h[i]-energy_renormalized_h[i-coulomb_matrix_red_index_i]);
		}
	}
}

/* Fill array electric_field_renormalized[]
 * Requires the matrix coulomb_matrix_red be computed ahead of time
 * */
void Device::renormalize_reducedCoulomb_Field(std::complex<double> *p_k_old)
{
	long int n = 0l;
	double red_jh=((double)coulomb_matrix_red_index_j);
	double inv_red_ih=1.0/((double)coulomb_matrix_red_index_i);

	// Calculate renormalization
	n=0l;
	for( int i=0; i<number_K_points; i+=coulomb_matrix_red_index_i )
	{
		electric_field_renormalized[i]  = 0.0;
		
		for( int iiq = i+1; iiq < number_K_points ; iiq+=coulomb_matrix_red_index_j )
		{
			electric_field_renormalized[i]  +=  coulomb_matrix_eh_red[n] * p_k_old[iiq];
			n++;
		}
		for( int iiq = i-1; iiq >=0 ; iiq-=coulomb_matrix_red_index_j )
		{
			electric_field_renormalized[i]  +=  coulomb_matrix_eh_red[n] * p_k_old[iiq];
			n++;
		}
	}
	if (coulomb_matrix_red_index_i > 1)
	{
		int i = 0;
		for(i=0; i<(number_K_points-coulomb_matrix_red_index_i); i+=coulomb_matrix_red_index_i )
		{
/*
			// Set parabolic coefficients
			int i0 = i;
			int i1 = i+coulomb_matrix_red_index_i;
			int i2 = i+2*coulomb_matrix_red_index_i;
			while (i2 >= number_K_points)
			{
				// On left edge
				i0 -= coulomb_matrix_red_index_i;
				i1 -= coulomb_matrix_red_index_i;
				i2 -= coulomb_matrix_red_index_i;
			}
			
			double t0 = K[i0];
			double t1 = K[i1];
			double t2 = K[i2];
			double delta0 = (t0-t1)*(t0-t2);
			double delta1 = (t1-t0)*(t1-t2);
			double delta2 = (t2-t0)*(t2-t1);
*/

			for( int j=1; j<coulomb_matrix_red_index_i; j++ )
			{
				// Linear
				double dL = (K[i+j]-K[i])/(K[i+coulomb_matrix_red_index_i]-K[i]);
				electric_field_renormalized[i+j] = electric_field_renormalized[i]+dL*(electric_field_renormalized[i+coulomb_matrix_red_index_i]-electric_field_renormalized[i]);
/*				
				// Parabolic
				double L0 = (K[i+j]-t1)*(K[i+j]-t2)/delta0;
				double L1 = (K[i+j]-t0)*(K[i+j]-t2)/delta1;
				double L2 = (K[i+j]-t0)*(K[i+j]-t1)/delta2;
				electric_field_renormalized[i+j] = electric_field_renormalized[i0]*L0 + electric_field_renormalized[i1]*L1 + electric_field_renormalized[i2]*L2;
*/
			}
		}
		
		// Interpolate the last few steps of the grid
		for( int j=(i-coulomb_matrix_red_index_i+1); j<number_K_points; j++)
		{
			double dL = (K[j]-K[i-coulomb_matrix_red_index_i])/(K[i]-K[i-coulomb_matrix_red_index_i]);
			electric_field_renormalized[j] = electric_field_renormalized[i-coulomb_matrix_red_index_i]+dL*(electric_field_renormalized[i]-electric_field_renormalized[i-coulomb_matrix_red_index_i]);
		}
	}
}


void Device::sbe_iterate(double t_sim, double DT)
{
	#ifdef MPI_BALANCE_WORKLOAD
	MPI_load->start();
	#endif

	#ifdef USE_DEVICE_TIMERS
	MainStat->start(getName());
	#endif
	// FOR HYSTERESIS PROGRAM
	//sbe_hysteresis_update_background_density(t_sim);

	//========================
	// Check thresholding
	//========================
	
	double normE = abs(getElectricField())/pulse_amplitude_threshold;
	double normE_field = normE/pulse_amplitude_threshold_field;
	double normE_analytic = normE/pulse_amplitude_threshold_analytic;
	bool SKIP_EQ          = false;
	bool SKIP_EQ_FIELD    = false;
	bool SKIP_EQ_ANALYTIC = false;

	if ((normE_analytic<=1)&&(SBE_ANALYTIC_THRESHOLDING == 1))
	{
		//Analytic update
		analyticDelay++;
		thresholdZeroFlag++;
		SKIP_EQ_ANALYTIC = true;

	} else if ((normE_field<=1)&&(SBE_FIELD_THRESHOLDING == 1))
	{
		//Field free update
		thresholdZeroFlag++;	
		SKIP_EQ_FIELD = true;
		
	} else if ((normE<=1)&&(SBE_COULOMB_THRESHOLDING == 1))
	{
		//Coulomb free update
		thresholdZeroFlag++;	
		SKIP_EQ = true;
	} else {
		
		//Full SBE
		thresholdZeroFlag=0; //Reset zeroing flag
	}

		
	//======================
	// Spontanious emission
	//======================
	#ifdef SBE_USE_SPONTAN_EMISSIONS
		double randnr;
		#if defined(__ICC) || defined(__INTEL_COMPILER)
			vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &randnr, -1.0, 1.0 ); // ICC or ICPC
		#elif defined(__GNUC__) || defined(__GNUG__)
			std::uniform_real_distribution<double> device_distribution(-1.0,1.0); // G++
			randnr = device_distribution(device_generator); //G++
		#endif

		double b0 = 1.0;

		device_spont_emission_phasor = b0*exp(I*Pi*randnr);
	#endif
	
	//When using field, analytic or coulomb threshlding, zero terms on first iteration
	if (thresholdZeroFlag==1)
	{
		// Remove any influence of the 2nBorn carrier scattering in this scenario
		#ifdef USE_HOLE_FILLING
			#if defined(USE_ISAK_HOLE_FILLING) || defined(USE_ISAK_HOLE_FILLING_TABLE)
				memset(carrier_scattering_rates_e_total,0,number_K_points*sizeof(double));
				memset(carrier_scattering_rates_h_total,0,number_K_points*sizeof(double));
			#else
				memset(carrier_scattering_rate_approximation_e,0,number_K_points*sizeof(double));
				memset(carrier_scattering_rate_approximation_h,0,number_K_points*sizeof(double));
			#endif
		#endif

	}
	
	//=================================
	// Integrator selection
	//=================================

	if (SKIP_EQ_ANALYTIC)
	{
		//First three steps skipped zero out macroscopic polarization
		if(analyticDelay<4)	
		{	
			std::complex<double> tmp_macpol = 0.0;
			macroscopic_polarization_m2 = getMacroscopicPolarization_m1();
			macroscopic_polarization_m1 = getMacroscopicPolarization();
			setMacroscopicPolarization(tmp_macpol);
		}
	} else{
		if (analyticDelay>0)
		{
			sbe_updateAnalytic(DT); //Update analytic form if steps were skipped
			analyticDelay=0;
		}


		#ifdef USE_DEVICE_TIMERS
		MainStat->start(getName(), "iterate_rk4");
		#endif

		// Integrate using RK4 only
		sbe_iterate_rk4(t_sim, DT,SKIP_EQ, SKIP_EQ_FIELD);
		
		// Integrate using AB4
		//sbe_iterate_ab4(t_sim, DT,SKIP_EQ, SKIP_EQ_FIELD);
		
		// Integrate using ABM4
		//sbe_iterate_abm4(t_sim, DT,SKIP_EQ, SKIP_EQ_FIELD);
		
		#ifdef USE_DEVICE_TIMERS
		MainStat->stop(getName(), "iterate_rk4");
		#endif

		//=================================
		// Macroscopic Polarization
		//=================================
		
		// Compute macroscopic polarization
		std::complex<double> tmp_macpol = 0.0;
		for(unsigned i = 0; i < number_K_points; i++)
		{
			tmp_macpol += KdK[i]*p_k[i];
		}
		tmp_macpol *= device_dipolemoment/(Pi*a0*a0);	//In rotating phase approximation
		tmp_macpol *= exp(-I*(device_pol_frequency - w0)*(t_sim + DT));
		macroscopic_polarization_m2 = getMacroscopicPolarization_m1();
		macroscopic_polarization_m1 = getMacroscopicPolarization();
		setMacroscopicPolarization(tmp_macpol);
		
		#ifdef USE_DEVICE_TIMERS
		MainStat->stop(getName());
		#endif
		#ifdef MPI_BALANCE_WORKLOAD
		MPI_load->stop();
		#endif
	}
}

void Device::sbe_updateAnalytic( double DT )
{
	for (unsigned i = 0; i < number_K_points; i++)
	{
		p_k[i]  = 0.0;
		ne_k[i] = fe_k[i] - (fe_k[i]-ne_k[i])*exp(-DT*analyticDelay*SBE_OCC_PUMP);
		nh_k[i] = fh_k[i] - (fh_k[i]-nh_k[i])*exp(-DT*analyticDelay*SBE_OCC_PUMP);
	}
}

void Device::sbe_iterate_abm4(double t_sim, double DT, bool SKIP_EQ, bool SKIP_EQ_FIELD)
{
	std::complex<double> E_tmp;
	
	if (ab4_prev_state != SKIP_EQ)
	{
		ab4_precompute_counter = 0;
		ab4_prev_state = SKIP_EQ;
	}
	
	double DT_T = DT/SBE_TIME_SCALE;
	
	if (ab4_precompute_counter < 3)
	{
		//===============
		// Store results
		//===============
		 p_k_ab4->updateStorage();
		ne_k_ab4->updateStorage();
		nh_k_ab4->updateStorage();
		E_tmp = device_dcv_hbar*getElectricField()*exp(I*(device_pol_frequency - w0)*t_sim);
		
		if (SKIP_EQ)
		{
			sbe_RHS_pk_SmallAmp(DT_T,  p_k_ab4->getContainerNr(0), p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_ne_SmallAmp(DT_T, ne_k_ab4->getContainerNr(0), p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_nh_SmallAmp(DT_T, nh_k_ab4->getContainerNr(0), p_k, ne_k, nh_k, E_tmp);
		} else {
		
			renormalize_reducedCoulomb(p_k, ne_k, nh_k);
	#ifdef USE_HOLE_FILLING
		carrier_scattering_calculate(2l, 4l,ne_k, nh_k, SKIP_EQ);
	#endif
			sbe_RHS_pk(DT_T,  p_k_ab4->getContainerNr(0), p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_ne(DT_T, ne_k_ab4->getContainerNr(0), p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_nh(DT_T, nh_k_ab4->getContainerNr(0), p_k, ne_k, nh_k, E_tmp);
		}
		ab4_precompute_counter++;
		
		sbe_iterate_rk4(t_sim, DT,SKIP_EQ, SKIP_EQ_FIELD);
		
	} else {
	
		// Store previous state
		ab4_prev_state = SKIP_EQ;
	
		// AB4
		p_k_ab4->updateStorage();
		ne_k_ab4->updateStorage();
		nh_k_ab4->updateStorage();
	
		std::complex<double> *p_k_f0 = p_k_ab4->getContainerNr(0);
		std::complex<double> *p_k_f1 = p_k_ab4->getContainerNr(1);
		std::complex<double> *p_k_f2 = p_k_ab4->getContainerNr(2);
		std::complex<double> *p_k_f3 = p_k_ab4->getContainerNr(3);
	
		double *ne_k_f0 = ne_k_ab4->getContainerNr(0);
		double *ne_k_f1 = ne_k_ab4->getContainerNr(1);
		double *ne_k_f2 = ne_k_ab4->getContainerNr(2);
		double *ne_k_f3 = ne_k_ab4->getContainerNr(3);
	
		double *nh_k_f0 = nh_k_ab4->getContainerNr(0);
		double *nh_k_f1 = nh_k_ab4->getContainerNr(1);
		double *nh_k_f2 = nh_k_ab4->getContainerNr(2);
		double *nh_k_f3 = nh_k_ab4->getContainerNr(3);
		
		E_tmp = device_dcv_hbar*getElectricField()*exp(I*(device_pol_frequency - w0)*t_sim);
		
		if (SKIP_EQ)
		{
			sbe_RHS_pk_SmallAmp(DT_T,  p_k_f0, p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_ne_SmallAmp(DT_T, ne_k_f0, p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_nh_SmallAmp(DT_T, nh_k_f0, p_k, ne_k, nh_k, E_tmp);
		} else {
			renormalize_reducedCoulomb(p_k, ne_k, nh_k);
	
		#ifdef USE_HOLE_FILLING
			carrier_scattering_calculate(2l, 4l,ne_k, nh_k, SKIP_EQ);
		#endif
	
			sbe_RHS_pk(DT_T,  p_k_f0, p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_ne(DT_T, ne_k_f0, p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_nh(DT_T, nh_k_f0, p_k, ne_k, nh_k, E_tmp);
		}
	
		// Predictor: AB4
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k_quad[0][i] =  p_k[i] + DT_T*(55.0* p_k_f0[i] - 59.0* p_k_f1[i] + 37.0* p_k_f2[i] - 9.0* p_k_f3[i])/24.0;
			ne_k_quad[0][i] = ne_k[i] + DT_T*(55.0*ne_k_f0[i] - 59.0*ne_k_f1[i] + 37.0*ne_k_f2[i] - 9.0*ne_k_f3[i])/24.0;
			nh_k_quad[0][i] = nh_k[i] + DT_T*(55.0*nh_k_f0[i] - 59.0*nh_k_f1[i] + 37.0*nh_k_f2[i] - 9.0*nh_k_f3[i])/24.0;
		}
		
		
		// Corrector: AM4
		E_tmp = device_dcv_hbar*getElectricField()*exp(I*(device_pol_frequency - w0)*(t_sim + DT));
		if (SKIP_EQ)
		{
			sbe_RHS_pk_SmallAmp(DT_T,  p_k_quad[1], p_k_quad[0], ne_k_quad[0], nh_k_quad[0], E_tmp);
			sbe_RHS_ne_SmallAmp(DT_T, ne_k_quad[1], p_k_quad[0], ne_k_quad[0], nh_k_quad[0], E_tmp);
			sbe_RHS_nh_SmallAmp(DT_T, nh_k_quad[1], p_k_quad[0], ne_k_quad[0], nh_k_quad[0], E_tmp);
		} else {
			renormalize_reducedCoulomb(p_k_quad[0], ne_k_quad[0], nh_k_quad[0]);
	
		#ifdef USE_HOLE_FILLING
			carrier_scattering_calculate(2l, 4l, ne_k_quad[0], nh_k_quad[0], SKIP_EQ);
		#endif
		
			sbe_RHS_pk(DT_T,  p_k_quad[1], p_k_quad[0], ne_k_quad[0], nh_k_quad[0], E_tmp);
			sbe_RHS_ne(DT_T, ne_k_quad[1], p_k_quad[0], ne_k_quad[0], nh_k_quad[0], E_tmp);
			sbe_RHS_nh(DT_T, nh_k_quad[1], p_k_quad[0], ne_k_quad[0], nh_k_quad[0], E_tmp);
		}		
		
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k[i] =  p_k[i] + DT_T*(9.0* p_k_quad[1][i] + 19.0* p_k_f0[i] - 5.0* p_k_f1[i] +  p_k_f2[i])/24.0;
			ne_k[i] = ne_k[i] + DT_T*(9.0*ne_k_quad[1][i] + 19.0*ne_k_f0[i] - 5.0*ne_k_f1[i] + ne_k_f2[i])/24.0;
			nh_k[i] = nh_k[i] + DT_T*(9.0*nh_k_quad[1][i] + 19.0*nh_k_f0[i] - 5.0*nh_k_f1[i] + nh_k_f2[i])/24.0;
		}
	}
}

void Device::sbe_iterate_ab4(double t_sim, double DT, bool SKIP_EQ, bool SKIP_EQ_FIELD)
{

	std::complex<double> E_tmp;
	
	if (ab4_prev_state != SKIP_EQ)
	{
		ab4_precompute_counter = 0;
		ab4_prev_state = SKIP_EQ;
	}

	double DT_T = DT/SBE_TIME_SCALE;
	
	if (ab4_precompute_counter < 3)
	{
		//======================
		// Store f(t,y) for AB
		//=====================
		 p_k_ab4->updateStorage();
		ne_k_ab4->updateStorage();
		nh_k_ab4->updateStorage();
		
		E_tmp = device_dcv_hbar*getElectricField()*exp(I*(device_pol_frequency - w0)*t_sim);
		
		if (SKIP_EQ)
		{
			
			sbe_RHS_pk_SmallAmp(DT_T,  p_k_ab4->getContainerNr(0), p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_ne_SmallAmp(DT_T, ne_k_ab4->getContainerNr(0), p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_nh_SmallAmp(DT_T, nh_k_ab4->getContainerNr(0), p_k, ne_k, nh_k, E_tmp);
		} else {
		
			renormalize_reducedCoulomb(p_k, ne_k, nh_k);
	#ifdef USE_HOLE_FILLING
		carrier_scattering_calculate(2l, 4l,ne_k, nh_k, SKIP_EQ);
	#endif
		
			sbe_RHS_pk(DT_T,  p_k_ab4->getContainerNr(0), p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_ne(DT_T, ne_k_ab4->getContainerNr(0), p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_nh(DT_T, nh_k_ab4->getContainerNr(0), p_k, ne_k, nh_k, E_tmp);
		}
		
		ab4_precompute_counter++;
		
		// DO RK4 step
		
		sbe_iterate_rk4(t_sim, DT,SKIP_EQ, SKIP_EQ_FIELD);
		
	} else {
	
		// Store previous state
		ab4_prev_state = SKIP_EQ;
	
		// AB4
		p_k_ab4->updateStorage();
		ne_k_ab4->updateStorage();
		nh_k_ab4->updateStorage();
	
		std::complex<double> *p_k_f0 = p_k_ab4->getContainerNr(0);
		std::complex<double> *p_k_f1 = p_k_ab4->getContainerNr(1);
		std::complex<double> *p_k_f2 = p_k_ab4->getContainerNr(2);
		std::complex<double> *p_k_f3 = p_k_ab4->getContainerNr(3);
	
		double *ne_k_f0 = ne_k_ab4->getContainerNr(0);
		double *ne_k_f1 = ne_k_ab4->getContainerNr(1);
		double *ne_k_f2 = ne_k_ab4->getContainerNr(2);
		double *ne_k_f3 = ne_k_ab4->getContainerNr(3);
	
		double *nh_k_f0 = nh_k_ab4->getContainerNr(0);
		double *nh_k_f1 = nh_k_ab4->getContainerNr(1);
		double *nh_k_f2 = nh_k_ab4->getContainerNr(2);
		double *nh_k_f3 = nh_k_ab4->getContainerNr(3);	
		
		E_tmp = device_dcv_hbar*getElectricField()*exp(I*(device_pol_frequency - w0)*t_sim);
		
		if (SKIP_EQ)
		{
			sbe_RHS_pk_SmallAmp(DT_T,  p_k_f0, p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_ne_SmallAmp(DT_T, ne_k_f0, p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_nh_SmallAmp(DT_T, nh_k_f0, p_k, ne_k, nh_k, E_tmp);
		} else {
		
			renormalize_reducedCoulomb(p_k, ne_k, nh_k);
		#ifdef USE_HOLE_FILLING
			carrier_scattering_calculate(2l, 4l,ne_k, nh_k, SKIP_EQ);
		#endif

			sbe_RHS_pk(DT_T,  p_k_f0, p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_ne(DT_T, ne_k_f0, p_k, ne_k, nh_k, E_tmp);
			sbe_RHS_nh(DT_T, nh_k_f0, p_k, ne_k, nh_k, E_tmp);
		}
	
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k[i] =  p_k[i] + DT_T*(55.0* p_k_f0[i] - 59.0* p_k_f1[i] + 37.0* p_k_f2[i] - 9.0* p_k_f3[i])/24.0;
			ne_k[i] = ne_k[i] + DT_T*(55.0*ne_k_f0[i] - 59.0*ne_k_f1[i] + 37.0*ne_k_f2[i] - 9.0*ne_k_f3[i])/24.0;
			nh_k[i] = nh_k[i] + DT_T*(55.0*nh_k_f0[i] - 59.0*nh_k_f1[i] + 37.0*nh_k_f2[i] - 9.0*nh_k_f3[i])/24.0;
		}
		
	}
}

void Device::sbe_iterate_rk4(double t_sim, double DT, bool SKIP_EQ, bool SKIP_EQ_FIELD)
{
	std::complex<double> E_tmp;
	
	double DT_T = DT/SBE_TIME_SCALE;
	
	if (SKIP_EQ_FIELD)
	{
		// RK4
		// Step 1, t = t_sim
		sbe_RHS_pk_FieldFree(DT_T,  p_k_quad[0], p_k, ne_k, nh_k);
		sbe_RHS_ne_FieldFree(DT_T, ne_k_quad[0], p_k, ne_k, nh_k);
		sbe_RHS_nh_FieldFree(DT_T, nh_k_quad[0], p_k, ne_k, nh_k);
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k_tmp[i] =  p_k[i] +  p_k_quad[0][i]*DT_T/2.0;
			ne_k_tmp[i] = ne_k[i] + ne_k_quad[0][i]*DT_T/2.0;
			nh_k_tmp[i] = nh_k[i] + nh_k_quad[0][i]*DT_T/2.0;
		}
		
		// Step 2, t = t_sim + DT_T/2
		
		sbe_RHS_pk_FieldFree(DT_T,  p_k_quad[1], p_k_tmp, ne_k_tmp, nh_k_tmp);
		sbe_RHS_ne_FieldFree(DT_T, ne_k_quad[1], p_k_tmp, ne_k_tmp, nh_k_tmp);
		sbe_RHS_nh_FieldFree(DT_T, nh_k_quad[1], p_k_tmp, ne_k_tmp, nh_k_tmp);
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k_tmp[i] =  p_k[i] +  p_k_quad[1][i]*DT_T/2.0;
			ne_k_tmp[i] = ne_k[i] + ne_k_quad[1][i]*DT_T/2.0;
			nh_k_tmp[i] = nh_k[i] + nh_k_quad[1][i]*DT_T/2.0;
		}
		
		// Step 3, t = t_sim + DT/2
		sbe_RHS_pk_FieldFree(DT_T,  p_k_quad[2], p_k_tmp, ne_k_tmp, nh_k_tmp);
		sbe_RHS_ne_FieldFree(DT_T, ne_k_quad[2], p_k_tmp, ne_k_tmp, nh_k_tmp);
		sbe_RHS_nh_FieldFree(DT_T, nh_k_quad[2], p_k_tmp, ne_k_tmp, nh_k_tmp);
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k_tmp[i] =  p_k[i] +  p_k_quad[2][i]*DT_T;
			ne_k_tmp[i] = ne_k[i] + ne_k_quad[2][i]*DT_T;
			nh_k_tmp[i] = nh_k[i] + nh_k_quad[2][i]*DT_T;
		}
		
		// Step 4, t = t_sim + DT
		sbe_RHS_pk_FieldFree(DT_T,  p_k_quad[3], p_k_tmp, ne_k_tmp, nh_k_tmp);
		sbe_RHS_ne_FieldFree(DT_T, ne_k_quad[3], p_k_tmp, ne_k_tmp, nh_k_tmp);
		sbe_RHS_nh_FieldFree(DT_T, nh_k_quad[3], p_k_tmp, ne_k_tmp, nh_k_tmp);

		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k[i] =  p_k[i]  + DT_T*( p_k_quad[0][i] + 2.0*( p_k_quad[1][i] +  p_k_quad[2][i]) +  p_k_quad[3][i])/6.0;
			ne_k[i] = ne_k[i] + DT_T*(ne_k_quad[0][i] + 2.0*(ne_k_quad[1][i] + ne_k_quad[2][i]) + ne_k_quad[3][i])/6.0;
			nh_k[i] = nh_k[i] + DT_T*(nh_k_quad[0][i] + 2.0*(nh_k_quad[1][i] + nh_k_quad[2][i]) + nh_k_quad[3][i])/6.0;
		}	

	}else if (SKIP_EQ)
	{
		// RK4
		// Step 1, t = t_sim
		E_tmp = device_dcv_hbar*getElectricField()*exp(I*(device_pol_frequency - w0)*t_sim);

		sbe_RHS_pk_SmallAmp(DT_T,  p_k_quad[0], p_k, ne_k, nh_k, E_tmp);
		sbe_RHS_ne_SmallAmp(DT_T, ne_k_quad[0], p_k, ne_k, nh_k, E_tmp);
		sbe_RHS_nh_SmallAmp(DT_T, nh_k_quad[0], p_k, ne_k, nh_k, E_tmp);
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k_tmp[i] =  p_k[i] +  p_k_quad[0][i]*DT_T/2.0;
			ne_k_tmp[i] = ne_k[i] + ne_k_quad[0][i]*DT_T/2.0;
			nh_k_tmp[i] = nh_k[i] + nh_k_quad[0][i]*DT_T/2.0;
		}
		
		// Step 2, t = t_sim + DT_T/2
		E_tmp = device_dcv_hbar*getElectricField_tp05()*exp(I*(device_pol_frequency - w0)*(t_sim + DT/2.0));
		
		sbe_RHS_pk_SmallAmp(DT_T,  p_k_quad[1], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		sbe_RHS_ne_SmallAmp(DT_T, ne_k_quad[1], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		sbe_RHS_nh_SmallAmp(DT_T, nh_k_quad[1], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k_tmp[i] =  p_k[i] +  p_k_quad[1][i]*DT_T/2.0;
			ne_k_tmp[i] = ne_k[i] + ne_k_quad[1][i]*DT_T/2.0;
			nh_k_tmp[i] = nh_k[i] + nh_k_quad[1][i]*DT_T/2.0;
		}
		
		// Step 3, t = t_sim + DT/2
		sbe_RHS_pk_SmallAmp(DT_T,  p_k_quad[2], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		sbe_RHS_ne_SmallAmp(DT_T, ne_k_quad[2], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		sbe_RHS_nh_SmallAmp(DT_T, nh_k_quad[2], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k_tmp[i] =  p_k[i] +  p_k_quad[2][i]*DT_T;
			ne_k_tmp[i] = ne_k[i] + ne_k_quad[2][i]*DT_T;
			nh_k_tmp[i] = nh_k[i] + nh_k_quad[2][i]*DT_T;
		}
		
		// Step 4, t = t_sim + DT
		E_tmp = device_dcv_hbar*getElectricField_tp1()*exp(I*(device_pol_frequency - w0)*(t_sim + DT));
	
		sbe_RHS_pk_SmallAmp(DT_T,  p_k_quad[3], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		sbe_RHS_ne_SmallAmp(DT_T, ne_k_quad[3], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		sbe_RHS_nh_SmallAmp(DT_T, nh_k_quad[3], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);

		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k[i] =  p_k[i]  + DT_T*( p_k_quad[0][i] + 2.0*( p_k_quad[1][i] +  p_k_quad[2][i]) +  p_k_quad[3][i])/6.0;
			ne_k[i] = ne_k[i] + DT_T*(ne_k_quad[0][i] + 2.0*(ne_k_quad[1][i] + ne_k_quad[2][i]) + ne_k_quad[3][i])/6.0;
			nh_k[i] = nh_k[i] + DT_T*(nh_k_quad[0][i] + 2.0*(nh_k_quad[1][i] + nh_k_quad[2][i]) + nh_k_quad[3][i])/6.0;
		}	
	} else {
		
		// RK4
		// Step 1, t = t_sim
		E_tmp = device_dcv_hbar*getElectricField()*exp(I*(device_pol_frequency - w0)*t_sim);
		renormalize_reducedCoulomb(p_k, ne_k, nh_k);
	
		#ifdef USE_HOLE_FILLING
			carrier_scattering_calculate(2l,1l,ne_k, nh_k, SKIP_EQ);
		#endif
		
		sbe_RHS_pk(DT_T,  p_k_quad[0], p_k, ne_k, nh_k, E_tmp);
		sbe_RHS_ne(DT_T, ne_k_quad[0], p_k, ne_k, nh_k, E_tmp);
		sbe_RHS_nh(DT_T, nh_k_quad[0], p_k, ne_k, nh_k, E_tmp);
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k_tmp[i] =  p_k[i] +  p_k_quad[0][i]*DT_T/2.0;
			ne_k_tmp[i] = ne_k[i] + ne_k_quad[0][i]*DT_T/2.0;
			nh_k_tmp[i] = nh_k[i] + nh_k_quad[0][i]*DT_T/2.0;
		}
		
		// Step 2, t = t_sim + DT/2
		E_tmp = device_dcv_hbar*getElectricField_tp05()*exp(I*(device_pol_frequency - w0)*(t_sim + DT/2.0));
		renormalize_reducedCoulomb(p_k, ne_k, nh_k);
		
		#ifdef USE_HOLE_FILLING
			carrier_scattering_calculate(2l,2l,ne_k_tmp, nh_k_tmp, SKIP_EQ);
		#endif
		
		sbe_RHS_pk(DT_T,  p_k_quad[1], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		sbe_RHS_ne(DT_T, ne_k_quad[1], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		sbe_RHS_nh(DT_T, nh_k_quad[1], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k_tmp[i] =  p_k[i] +  p_k_quad[1][i]*DT_T/2.0;
			ne_k_tmp[i] = ne_k[i] + ne_k_quad[1][i]*DT_T/2.0;
			nh_k_tmp[i] = nh_k[i] + nh_k_quad[1][i]*DT_T/2.0;
		}
		
		// Step 3, t = t_sim + DT/2
		//E_tmp = device_dcv_hbar*getElectricField()*exp(I*(device_pol_frequency - w0)*(t_sim + DT/2.0));
		renormalize_reducedCoulomb(p_k, ne_k, nh_k);
		
		#ifdef USE_HOLE_FILLING
			carrier_scattering_calculate(2l, 3l,ne_k_tmp, nh_k_tmp, SKIP_EQ);
		#endif
		
		sbe_RHS_pk(DT_T,  p_k_quad[2], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		sbe_RHS_ne(DT_T, ne_k_quad[2], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		sbe_RHS_nh(DT_T, nh_k_quad[2], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k_tmp[i] =  p_k[i] +  p_k_quad[2][i]*DT_T;
			ne_k_tmp[i] = ne_k[i] + ne_k_quad[2][i]*DT_T;
			nh_k_tmp[i] = nh_k[i] + nh_k_quad[2][i]*DT_T;
		}
		
		// Step 4, t = t_sim + DT
		E_tmp = device_dcv_hbar*getElectricField_tp1()*exp(I*(device_pol_frequency - w0)*(t_sim + DT));
		renormalize_reducedCoulomb(p_k, ne_k, nh_k);
		
		#ifdef USE_HOLE_FILLING
			carrier_scattering_calculate(2l, 4l,ne_k_tmp, nh_k_tmp, SKIP_EQ);
		#endif
		
		sbe_RHS_pk(DT_T,  p_k_quad[3], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		sbe_RHS_ne(DT_T, ne_k_quad[3], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		sbe_RHS_nh(DT_T, nh_k_quad[3], p_k_tmp, ne_k_tmp, nh_k_tmp, E_tmp);
		
		for(unsigned i = 0; i < number_K_points; i++)
		{
			 p_k[i] =  p_k[i] + DT_T*( p_k_quad[0][i] + 2.0*( p_k_quad[1][i] +  p_k_quad[2][i]) +  p_k_quad[3][i])/6.0;
			ne_k[i] = ne_k[i] + DT_T*(ne_k_quad[0][i] + 2.0*(ne_k_quad[1][i] + ne_k_quad[2][i]) + ne_k_quad[3][i])/6.0;
			nh_k[i] = nh_k[i] + DT_T*(nh_k_quad[0][i] + 2.0*(nh_k_quad[1][i] + nh_k_quad[2][i]) + nh_k_quad[3][i])/6.0;
		}	
	}
}


/* Hysteresis testing in QWs
 * This function will GRADUALLY modify the background carrier density in the QW
 * between until either of TFINAL or DENS_STOP is reached.
 * TSTART -> When to start changing
 * TFINAL -> When to stop changing
 * DENS_S -> Starting density, SHOULD BE THE starting background denisty in the QW
 * DENS_F -> Final density at TFINAL.
 * DENS_STOP -> Stopping condition 2: Stops when this density is reached
 * */
void Device::sbe_hysteresis_update_background_density(double t_sim)
{
	if (device_name.substr(0,2) == "QW")
	{
		// START AND END OF TUNING
		double TSTART = 100.0*ns;
		double TFINAL = 150.0*ns; // Time when done tuning

		// DENSITY INITIAL AND FINAL
		double DENS_S = 5.0e16; // Initial density, should be pump density in QW
		double DENS_F = 2.0e16; // Final density


		double DENS_STOP = 2.3e16; // if This density is reached, then stop reducing density

		// Initialize time
		if ((t_sim > TSTART)&&(t_sim <= TFINAL))
		{
				// If DENS_S > DENS_F
				if (device_density > DENS_STOP)
				{
						// Tune from start to finish smoothly
						device_density = DENS_S + (DENS_F - DENS_S)*(t_sim - TSTART)/(TFINAL-TSTART);
				
						// Set chemical potential
						device_chemical_potential_e = misc_get_chem_potential_analytic_2d(device_density, device_me, device_background_temp_e);
						device_chemical_potential_h = misc_get_chem_potential_analytic_2d(device_density, device_mh, device_background_temp_h);

						for(unsigned i = 0; i < number_K_points; i++)
						{
								fe_k[i] = misc_get_fermi_distribution(device_background_temp_e, device_chemical_potential_e, device_me, K[i]);
								fh_k[i] = misc_get_fermi_distribution(device_background_temp_h, device_chemical_potential_h, device_mh, K[i]);
						}
				} else {
						// Force correct
						device_density = DENS_STOP;
				
						// Set chemical potential
						device_chemical_potential_e = misc_get_chem_potential_analytic_2d(device_density, device_me, device_background_temp_e);
						device_chemical_potential_h = misc_get_chem_potential_analytic_2d(device_density, device_mh, device_background_temp_h);

						for(unsigned i = 0; i < number_K_points; i++)
						{
							fe_k[i] = misc_get_fermi_distribution(device_background_temp_e, device_chemical_potential_e, device_me, K[i]);
							fh_k[i] = misc_get_fermi_distribution(device_background_temp_h, device_chemical_potential_h, device_mh, K[i]);
						}
				}
		}
		// Does nothing otherwise
	}
}

void Device::sbe_set_real_pump_model(double W0, double E0, double ETA, double DT)
{
	device_pump_W0 = W0;	// Central frequency
	device_pump_E0 = E0; 	// Field contributing to Carriers
	device_pump_ETA = ETA; 	// Pump width

	
	double dcv = device_dcv_hbar/SBE_TIME_SCALE;
	double const1 = -SBE_TIME_SCALE*(dcv*dcv)*device_pump_E0*device_pump_E0;
	
	if (device_pump_wk == NULL)
	{
		device_pump_wk = new double[number_K_points];
	} 
	
	for(unsigned i = 0; i < number_K_points; i++)
	{
		double wk = device_pol_frequency + sumE_k[i]/SBE_TIME_SCALE;
		device_pump_wk[i] = const1*(2.0*device_pump_ETA/((wk - device_pump_W0)*(wk - device_pump_W0) + device_pump_ETA*device_pump_ETA));
	}
	
	if (DT > 0)
	{
	
		// sbe_iterate
		double Nsum_old = 0, Nsum;
		for(int i = 0; i < number_K_points; i++)
		{
			Nsum_old += ne_k[i]*KdK[i];
		}
		
		cout << "Iterating qw distributions into equilibirum... " << endl;
/*
		double err = 1;
		double t=0;
		while (err > 1.0e-4)
		{
			this->sbe_iterate(0, DT);
			
			
			Nsum = 0;
			for(int i =0; i < number_K_points; i++)
			{
				Nsum += ne_k[i]*K[i]*dK[i];
			}
			
			err = abs(Nsum-Nsum_old)/Nsum_old;
		
			Nsum_old = Nsum;
			
			t += DT;
		}
		Nsum_old /= (Pi*a0*a0);
*/
		//cout << " in t = " << t/ps << " [ps], with Ne = " << Nsum_old << " [1/m^2]" << endl;
	}
	
	double Nsum = 0;
	for(int i = 0; i < number_K_points; i++)
	{
		Nsum += ne_k[i]*KdK[i];
	}
	Nsum /= (Pi*a0*a0);
	
	std::stringstream fileName;
	fileName << getToFileOutputKey() << "equilibrium_carrier_density_" << getName() << ".dat";
	saveBinary(fileName.str(), &Nsum, 1);
}

void Device::sbe_set_background_carrier_density(double new_dens)
{
	cout << "Device: Changing density to = " << new_dens << endl;
	device_density = new_dens;
	
	// Set chemical potential
	device_chemical_potential_e = misc_get_chem_potential_analytic_2d(device_density, device_me, device_background_temp_e);
	device_chemical_potential_h = misc_get_chem_potential_analytic_2d(device_density, device_mh, device_background_temp_h);

	// Initialize background fe_k, fh_k
	if (fe_k == NULL)
	{
		fe_k = new double[number_K_points];
	}
	if (fh_k == NULL)
	{
		fh_k = new double[number_K_points];
	}
	for(unsigned i = 0; i < number_K_points; i++)
	{
		fe_k[i]      = misc_get_fermi_distribution(device_background_temp_e, device_chemical_potential_e, device_me, K[i]);
		fh_k[i]      = misc_get_fermi_distribution(device_background_temp_h, device_chemical_potential_h, device_mh, K[i]);
	}
}


/* Calculate full screening and screened coulomb potential
 * Screening formula found on page 140 in: 
 * "Quantum Theory of Optical and Electronic Properties of Semiconductors" 
 * 5th Edition by Hartmut Haug and Stephan W. Koch
 * */
void Device::carrier_scattering_screening_full(double *ne_k_old, double *nh_k_old)
{
	
	
	// Screened V_q at q=0
	// V_0 = V_scale/kappa
	/*
	coulomb_potential_normalized_q0 = 0.0;
	{
		double fact_e_q = 0.0;
		double fact_h_q = 0.0;
		
		for(int k = 0; k < number_K_points; k++)
		{
			// Electrons
			double dfe_k    = a0*interpolate_cubic_spline_ne->evaluate_df(K[k]);
			//double dfe_k = a0*(ne_k_old[k]-ne_k_old[k-1])/(K[k]-K[k-1]);
			double dEe_k  = hbar*hbar*K[k]/(device_me*a0);
			double fact_e_q_tmp = dfe_k/dEe_k;
			
			// Holes
			double dfh_k    = a0*interpolate_cubic_spline_nh->evaluate_df(K[k]);
			//double dfh_k = a0*(nh_k_old[k]-nh_k_old[k-1])/(K[k]-K[k-1]);
			double dEh_k  = hbar*hbar*K[k]/(device_mh*a0);
			double fact_h_q_tmp = dfh_k/dEh_k;
			
			fact_e_q += K[k]*dK[k]*fact_e_q_tmp;
			fact_h_q += K[k]*dK[k]*fact_h_q_tmp;
		}
		
		double integration_volume = 1.0/(Pi*a0*a0);
		fact_e_q *= integration_volume;
		fact_h_q *= integration_volume;
		
		coulomb_potential_normalized_q0 = 1.0/(-(fact_e_q + fact_h_q));
	}
	coulomb_potential_normalized_q0 = coulomb_potential_normalized_q0/(e*e*a0/(2.0*eps*eps0)); // Normalized
	*/
	
	
	for(int k = 0; k < number_K_points; k++)
	{
		double fact_e_q = 0.0;
		double fact_h_q = 0.0;
		for(int q = 0; q < number_Qsc_points; q++)
		{
			double fact_e_q_tmp = 0.0;
			double fact_h_q_tmp = 0.0;
			
			for(int th_q = 0; th_q < number_Eta_points; th_q++)
			{
				double angle_th_k_q  = Eta[th_q];
				double k_q           = sqrt(K[k ]*K[k ] + Q_sc[q]*Q_sc[q] - 2.0*K[k ]*Q_sc[q]*cos(angle_th_k_q));
				
				if ((k_q >= K[0])&&(k_q <= K[number_K_points-1]))
				{
					// Electrons
					//double fe_k_q  = interpolate_cubic_spline_ne->evaluate(k_q);
					//double fe_q    = interpolate_cubic_spline_ne->evaluate(Q_sc[q]);
					double fe_k_q	= misc_interpolate_K_array_parabolic(k_q, K, ne_k_old, number_K_points);
					double fe_q	= misc_interpolate_K_array_parabolic(Q_sc[q], K, ne_k_old, number_K_points);
					//double fe_k_q  = misc_interpolate_K_array(k_q, K, ne_k_old, number_K_points,0.0);
					//double fe_q    = misc_interpolate_K_array(Q_sc[q], K, ne_k_old, number_K_points,0.0);
					double Ee_k_q  = hbar*hbar*(k_q*k_q)/(2.0*device_me*a0*a0);
					double Ee_q    = hbar*hbar*Q_sc[q]*Q_sc[q]/(2.0*device_me*a0*a0);
					fact_e_q_tmp   += ((fe_k_q-fe_q)/(Ee_k_q - Ee_q))*dEta[th_q];

					// Holes
					//double fh_k_q  = interpolate_cubic_spline_nh->evaluate(k_q);
					//double fh_q    = interpolate_cubic_spline_nh->evaluate(Q_sc[q]);
					double fh_k_q	= misc_interpolate_K_array_parabolic(k_q, K, nh_k_old, number_K_points);
					double fh_q	= misc_interpolate_K_array_parabolic(Q_sc[q], K, nh_k_old, number_K_points);
					//double fh_k_q  = misc_interpolate_K_array(k_q   , K, nh_k_old, number_K_points, 0.0);
					//double fh_q    = misc_interpolate_K_array(Q_sc[q], K, nh_k_old, number_K_points, 0.0);
					double Eh_k_q  = hbar*hbar*(k_q*k_q)/(2.0*device_mh*a0*a0);
					double Eh_q    = hbar*hbar*Q_sc[q]*Q_sc[q]/(2.0*device_mh*a0*a0);
					fact_h_q_tmp   += ((fh_k_q-fh_q)/(Eh_k_q - Eh_q))*dEta[th_q];
				}
			}
			
			fact_e_q += Q_sc[q]*dQ_sc[q]*fact_e_q_tmp;
			fact_h_q += Q_sc[q]*dQ_sc[q]*fact_h_q_tmp;
		}
		
		double integration_volume = 1.0/(2.0*Pi*Pi*a0*a0); // Multiply by 2.0 for the half interval of angles
		fact_e_q *= integration_volume;
		fact_h_q *= integration_volume;
		
		// coulomb matrix element
		double col_fact = e*e/(2.0*eps0*eps);
		//double Ve_q      = (col_fact/(K[k]/a0))*coulomb_potential_confinement_ee[k];
		//double Vh_q      = (col_fact/(K[k]/a0))*coulomb_potential_confinement_hh[k];
		double Veh_q     = (col_fact/(K[k]/a0))*coulomb_potential_confinement_eh[k];
		
		coulomb_potential_epsilon_inv[k] = (1.0 - 2.0*Veh_q*(fact_e_q + fact_h_q)); // Veh instead of Ve or Vh. Jorg change April 10 2018.
	}
	
	// Calculate screened coulomb matrix
	for(int k = 0; k < number_K_points; k++)
	{
		//double col_fact = e*e*a0/(2.0*eps0*eps); // factor scaled out below
		double inv_k_eps = 1.0/(K[k]*coulomb_potential_epsilon_inv[k]);
		coulomb_potential_normalized_ee[k] = inv_k_eps*coulomb_potential_confinement_ee[k];; // screened matrix Vq
		coulomb_potential_normalized_hh[k] = inv_k_eps*coulomb_potential_confinement_hh[k]; // screened matrix Vq
		coulomb_potential_normalized_eh[k] = inv_k_eps*coulomb_potential_confinement_eh[k]; // screened matrix Vq
	}
}

/* Calculate FAST screening and screened coulomb potential
 * from the pre-calculated index set found in
 * coulomb_potential_epsilon_index
 * */
void Device::carrier_scattering_screening_full_fast(double *ne_k_old, double *nh_k_old)
{
	//double col_fact = e*e/(2.0*eps0*eps);

	#ifdef USE_OPENMP
	#pragma omp parallel for schedule(static,OMP_THREADS_LEVEL_2) num_threads(OMP_THREADS_LEVEL_2)
	#endif
	for(int k = 0; k < number_K_points; k++)
	{
		double fact_q = 0.0;
		int index_ee_max = coulomb_potential_epsilon_index->get_num_cols(k);
		for(int index_ee = 0; index_ee < index_ee_max; index_ee++)
		{
			int i_q 	= coulomb_potential_epsilon_index->I_q(k,index_ee);
			int i_kq 	= coulomb_potential_epsilon_index->I_kq(k,index_ee);
			
			double q_b0, q_b1, q_b2, kq_b0, kq_b1, kq_b2;
			coulomb_potential_epsilon_index->beta_q(k,index_ee , &q_b0 , &q_b1 , &q_b2);
			coulomb_potential_epsilon_index->beta_kq(k,index_ee, &kq_b0, &kq_b1, &kq_b2);

			double fe_q 	= ne_k_old[i_q]*q_b0   + ne_k_old[i_q+1]*q_b1   + ne_k_old[i_q+2]*q_b2;
			double fe_k_q 	= ne_k_old[i_kq]*kq_b0 + ne_k_old[i_kq+1]*kq_b1 + ne_k_old[i_kq+2]*kq_b2;
			
			double fh_q 	= nh_k_old[i_q]*q_b0   + nh_k_old[i_q+1]*q_b1   + nh_k_old[i_q+2]*q_b2;
			double fh_k_q 	= nh_k_old[i_kq]*kq_b0 + nh_k_old[i_kq+1]*kq_b1 + nh_k_old[i_kq+2]*kq_b2;

			double weight_e = coulomb_potential_epsilon_index->F_e(k,index_ee);
			double weight_h = coulomb_potential_epsilon_index->F_h(k,index_ee);

			fact_q +=	weight_e*(fe_k_q-fe_q) + weight_h*(fh_k_q-fh_q);
		}

		coulomb_potential_epsilon_inv[k] = (1.0 - fact_q); // Calculate eps_inv = 1/eps, such that V_screened = V_q*eps = V_q/eps_inv 
	}	
	
	// Calculate screened 1D coulomb matrix (The factor col_fact*a0 has already been divided away)
	for(int k = 0; k < number_K_points; k++)
	{
		double eps_inv = 1.0/(K[k]*coulomb_potential_epsilon_inv[k]);
		coulomb_potential_normalized_ee[k] = coulomb_potential_confinement_ee[k]*eps_inv; // screened matrix Vq
		coulomb_potential_normalized_hh[k] = coulomb_potential_confinement_hh[k]*eps_inv; // screened matrix Vq
		coulomb_potential_normalized_eh[k] = coulomb_potential_confinement_eh[k]*eps_inv; // screened matrix Vq
	}
}

void Device::carrier_scattering_rate_approx(double * ne_k_current, double * nh_k_current)
{
	//====================================
	// Find instantanious temp
	//====================================
	carrier_scattering_rate_set_instant_temperature(ne_k_current, nh_k_current);
	
	device_inst_temp_e_prev = device_inst_temp_e;
	device_inst_temp_h_prev = device_inst_temp_h;
	
	carrier_scattering_rate_set_background_fermi();

	for(int i  =0; i < number_K_points; i++)
	{
		carrier_scattering_rate_approximation_e[i] = -(ne_k_current[i]-fe_k_inst[i])*SBE_OCC_HOLE;
		carrier_scattering_rate_approximation_h[i] = -(nh_k_current[i]-fh_k_inst[i])*SBE_OCC_HOLE;
	}

	// Particle number balance
	double error_ne = 0.0;
	double error_ne_p = 0.0;
	double error_ne_m = 0.0;
	for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
	{
		double tmp = carrier_scattering_rate_approximation_e[k]*KdK[k];
		error_ne += tmp;
		if (carrier_scattering_rate_approximation_e[k] < 0)
		{
			error_ne_m += tmp;
		} else {
			error_ne_p += tmp;
		}
	}
	
	if (error_ne < 0.0)
	{
		if (error_ne_m < 0.0)
		{
			double tmp = 1.0-error_ne/error_ne_m;
			for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
			{
				if (carrier_scattering_rate_approximation_e[k] < 0.0)
				{
					carrier_scattering_rate_approximation_e[k] *= tmp;
				}
			}
		}
	} else {
		if (error_ne_p > 0.0)
		{
			double tmp = 1.0-error_ne/error_ne_p;
			for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
			{
				if (carrier_scattering_rate_approximation_e[k] > 0)
				{
					carrier_scattering_rate_approximation_e[k] *= tmp;
				}
			}
		}
	}
	
	double error_nh = 0.0;
	double error_nh_p = 0.0;
	double error_nh_m = 0.0;
	for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
	{
		double tmp = carrier_scattering_rate_approximation_h[k]*KdK[k];
		error_nh += tmp;
		if (carrier_scattering_rate_approximation_h[k] < 0)
		{
			error_nh_m += tmp;
		} else {
			error_nh_p += tmp;
		}
	}
	
	if (error_nh < 0.0)
	{
		if (error_nh_m < 0.0)
		{
			double tmp = 1.0-error_nh/error_nh_m;
			for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
			{
				if (carrier_scattering_rate_approximation_h[k] < 0.0)
				{
					carrier_scattering_rate_approximation_h[k] *= tmp;
				}
			}
		}
	} else {
		if (error_nh_p > 0.0)
		{
			double tmp = 1.0-error_nh/error_nh_p;
			for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
			{
				if (carrier_scattering_rate_approximation_h[k] > 0)
				{
					carrier_scattering_rate_approximation_h[k] *= tmp;
				}
			}
		}
	}
	// carrier-carrier scattering: Conservation of Energy
	{
		double sumE_e = 0.0;
		double sumE_h = 0.0;
		for(int k = 0; k < number_K_points; k++) 
		{
			sumE_e += Ek_e[k]*carrier_scattering_rate_approximation_e[k]*KdK[k];
			sumE_h += Ek_h[k]*carrier_scattering_rate_approximation_h[k]*KdK[k];
		}

		if ((fabs(sumE_e) >0.0)&&(fabs(sumE_h) > 0.0))
		{
			// Minimum energy
			if (fabs(sumE_e) < fabs(sumE_h))
			{
				double const1 = fabs(sumE_e/sumE_h);
				for(int k = 0; k < number_K_points; k++) 
				{
					carrier_scattering_rate_approximation_h[k] *= const1;
				}
			} else {
				double const1 = fabs(sumE_h/sumE_e);
				for(int k = 0; k < number_K_points; k++) 
				{
					carrier_scattering_rate_approximation_e[k] *= const1;
				}

			}
			
		} else if ((fabs(sumE_e) == 0.0)&&(fabs(sumE_h) > 0.0))
		{
			for(int k = 0; k < number_K_points; k++) 
			{
				carrier_scattering_rate_approximation_h[k] = 0.0; // set min value
			}
		} else if ((fabs(sumE_e) > 0.0)&&(fabs(sumE_h) == 0.0))
		{
			for(int k = 0; k < number_K_points; k++) 
			{
				carrier_scattering_rate_approximation_e[k] = 0.0; // set min value
			}
		}
	}
}


void Device::carrier_scattering_method2_table(double *ne_k_current, double *nh_k_current)
{
	// Initial values
	double current_Te = device_inst_temp_e;
	double current_Th = device_inst_temp_h;
	double Nsum_e = 0, Nsum_h = 0;

	double current_density_e = 0.0;
	double current_density_h = 0.0;
	for(unsigned i = 0; i < number_K_points; i++)
	{
		current_density_e += ne_k_current[i]*KdK[i];
		current_density_h += nh_k_current[i]*KdK[i];
	}
	current_density_e /= (Pi*a0*a0);
	current_density_h /= (Pi*a0*a0);

//	double current_density_e = misc_simpsons_quadrature_kdk(ne_k_current);
//	double current_density_h = misc_simpsons_quadrature_kdk(nh_k_current);
	
	//cout << "Trying to evaluate the point = (" << current_Te <<  ", " << current_Th << ", " << current_density_e << ")" << endl;

	
	double current_density = current_density_e/1.0e16; // e and h density is the same


	// Fill carrier_scattering_table_gridpoints, assuming SAME GRIDS FOR e-in, e-out, h-in, h-out
	carrier_scattering_table_e_in->prepare_interpolation(current_Te, current_Th, current_density, carrier_scattering_table_gridpoints, &carrier_scattering_table_gridpoints_numbers);
	//carrier_scattering_table_e_out->prepare_interpolation(current_Te, current_Th, current_density, carrier_scattering_table_gridpoints, &carrier_scattering_table_gridpoints_numbers);
	//carrier_scattering_table_h_in->prepare_interpolation(current_Te, current_Th, current_density, carrier_scattering_table_gridpoints, &carrier_scattering_table_gridpoints_numbers);
	//carrier_scattering_table_h_out->prepare_interpolation(current_Te, current_Th, current_density, carrier_scattering_table_gridpoints, &carrier_scattering_table_gridpoints_numbers);
	
	
	/*
	if (carrier_scattering_table_gridpoints_numbers == 0)
	{
		//cout << "No new points needed" << endl;
	} else {
		cout << "Number of new points requested = " << carrier_scattering_table_gridpoints_numbers << endl;
		for(int i =0; i < 8; i++)
		{
			if (i < carrier_scattering_table_gridpoints_numbers)
			{
				cout << "[new]: ";
			} else {
				cout << "[---]: ";
			}
			for(int j = 0; j < 6; j++)
			{
				cout << carrier_scattering_table_gridpoints[i][j] << " ";
			}
			cout << endl;
		}
	}
	*/
	
	double ne_k_dens[number_K_points]; // Tmp array
	double nh_k_dens[number_K_points]; // Tmp array
	
	for(int i = 0; i < 8; i++)
	{
		if (i < carrier_scattering_table_gridpoints_numbers)
		{
			//cout << "Calculating gridpoint: (Te,Th,dens) = (" << carrier_scattering_table_gridpoints[i][0] << ", " << carrier_scattering_table_gridpoints[i][1] << ", " << carrier_scattering_table_gridpoints[i][2]*1.0e16 << ") ...." << endl;
			// Do calculation
			// Update densities and chem potential for new density
			double chem_pot_e = misc_get_chem_potential_analytic_2d(carrier_scattering_table_gridpoints[i][2]*1.0e16, device_me, carrier_scattering_table_gridpoints[i][0]);
			double chem_pot_h = misc_get_chem_potential_analytic_2d(carrier_scattering_table_gridpoints[i][2]*1.0e16, device_mh, carrier_scattering_table_gridpoints[i][1]);
			for(unsigned j = 0; j < number_K_points; j++)
			{
				ne_k_dens[j]      = misc_get_fermi_distribution(carrier_scattering_table_gridpoints[i][0], chem_pot_e, device_me, K[j]);
				nh_k_dens[j]      = misc_get_fermi_distribution(carrier_scattering_table_gridpoints[i][1], chem_pot_h, device_mh, K[j]);
			}
			
			// Calculate screening and scattering rates
			//updateScreening_reducedCoulomb_full(ne_k_dens, nh_k_dens);
			carrier_scattering_method2_parallel(ne_k_dens, nh_k_dens);
			

			// Store gridpoints
			carrier_scattering_table_e_in->update_interpolation_grid(carrier_scattering_table_gridpoints[i][3], carrier_scattering_table_gridpoints[i][4], carrier_scattering_table_gridpoints[i][5], carrier_scattering_rates_e_in);
			carrier_scattering_table_e_out->update_interpolation_grid(carrier_scattering_table_gridpoints[i][3], carrier_scattering_table_gridpoints[i][4], carrier_scattering_table_gridpoints[i][5], carrier_scattering_rates_e_out);
			carrier_scattering_table_h_in->update_interpolation_grid(carrier_scattering_table_gridpoints[i][3], carrier_scattering_table_gridpoints[i][4], carrier_scattering_table_gridpoints[i][5], carrier_scattering_rates_h_in);
			carrier_scattering_table_h_out->update_interpolation_grid(carrier_scattering_table_gridpoints[i][3], carrier_scattering_table_gridpoints[i][4], carrier_scattering_table_gridpoints[i][5], carrier_scattering_rates_h_out);
		}
	}
	
	// Calculate the requested value
	carrier_scattering_table_e_in-> evalF(carrier_scattering_rates_e_in, current_Te, current_Th, current_density);
	carrier_scattering_table_e_out->evalF(carrier_scattering_rates_e_out, current_Te, current_Th, current_density);
	carrier_scattering_table_h_in-> evalF(carrier_scattering_rates_h_in, current_Te, current_Th, current_density);
	carrier_scattering_table_h_out->evalF(carrier_scattering_rates_h_out, current_Te, current_Th, current_density);
	
	// Calculate total scattering
	for(int j = 0; j < number_K_points; j++)
	{
		// Method 2
		carrier_scattering_rates_e_total[j] = carrier_scattering_rates_e_in[j] - ne_k_current[j]*carrier_scattering_rates_e_out[j];
		carrier_scattering_rates_h_total[j] = carrier_scattering_rates_h_in[j] - nh_k_current[j]*carrier_scattering_rates_h_out[j];
	}
	
	// Save data to file: Internally this ensures that a save only happens when 80 new points are computed
	//carrier_scattering_table_e_in->file_save(80);
	//carrier_scattering_table_e_out->file_save(80);
	//carrier_scattering_table_h_in->file_save(80);
	//carrier_scattering_table_h_out->file_save(80);

	
	//========================================
	// Conservation of particle numbers
	// Jorg particle number balance
	double error_ne = 0.0;
	double error_ne_p = 0.0;
	double error_ne_m = 0.0;
	for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
	{
		double tmp = carrier_scattering_rates_e_total[k]*KdK[k];
		error_ne += tmp;
		if (carrier_scattering_rates_e_total[k] < 0)
		{
			error_ne_m += tmp;
		} else {
			error_ne_p += tmp;
		}
	}
	
	if (error_ne < 0.0)
	{
		if (error_ne_m < 0.0)
		{
			double tmp = 1.0-error_ne/error_ne_m;
			for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
			{
				if (carrier_scattering_rates_e_total[k] < 0.0)
				{
					carrier_scattering_rates_e_total[k] *= tmp;
				}
			}
		}
	} else {
		if (error_ne_p > 0.0)
		{
			double tmp = 1.0-error_ne/error_ne_p;
			for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
			{
				if (carrier_scattering_rates_e_total[k] > 0)
				{
					carrier_scattering_rates_e_total[k] *= tmp;
				}
			}
		}
	}
	
	double error_nh = 0.0;
	double error_nh_p = 0.0;
	double error_nh_m = 0.0;
	for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
	{
		double tmp = carrier_scattering_rates_h_total[k]*KdK[k];
		error_nh += tmp;
		if (carrier_scattering_rates_h_total[k] < 0)
		{
			error_nh_m += tmp;
		} else {
			error_nh_p += tmp;
		}
	}
	
	if (error_nh < 0.0)
	{
		if (error_nh_m < 0.0)
		{
			double tmp = 1.0-error_nh/error_nh_m;
			for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
			{
				if (carrier_scattering_rates_h_total[k] < 0.0)
				{
					carrier_scattering_rates_h_total[k] *= tmp;
				}
			}
		}
	} else {
		if (error_nh_p > 0.0)
		{
			double tmp = 1.0-error_nh/error_nh_p;
			for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
			{
				if (carrier_scattering_rates_h_total[k] > 0)
				{
					carrier_scattering_rates_h_total[k] *= tmp;
				}
			}
		}
	}
}

/* Carrier scattering calculation wrapper
 * This function is called by RK4 to run whatever carrier scattering function is selected
 * The two arguments time_counter and rk_step is for Jorg's 2ndBorn calculation to distinguish init.
 * Scattering functions should be init. before running this function.
 * */
void Device::carrier_scattering_calculate(long int time_counter, long int rk_step, double *ne_k_current, double *nh_k_current, bool SKIP_EQ)
{
//	#ifdef USE_DEVICE_TIMERS
//	MainStat->start(getName(),"c-c scatt");
//	#endif

	if (SKIP_EQ)
	{
		carrier_scattering_rate_approx(ne_k_current, nh_k_current);

		cout << "..." << endl;
		cout << "maybe problem.." << endl;
		exit(-1);

	} else {

		// Check if distributions have changed
		bool calc = false;
		double err_e_max = 0.0;
		double err_h_max = 0.0;
		double err_e,err_h;
		if (CARRIER_DEVIATION_TOL > 0)
		{
			for(int i = 0; i < number_K_points; i++)
			{
				//if ((abs(ne_k_current[i] -carrier_scattering_prev_ne[i]) > 1e-3)||(abs(carrier_scattering_prev_nh[i]-nh_k_current[i]) > 1e-3))
				{
					if (fabs(carrier_scattering_prev_ne[i]+ne_k_current[i]) > 0.0)
					{
						double tmp_err_e = fabs(carrier_scattering_prev_ne[i] - ne_k_current[i])/fabs(carrier_scattering_prev_ne[i]+ ne_k_current[i]);
						if (tmp_err_e > CARRIER_DEVIATION_TOL)
						{
							err_e_max = tmp_err_e;
						}
					} else {
						err_e_max = CARRIER_DEVIATION_TOL*2.0;

					}
					if (fabs(carrier_scattering_prev_nh[i]+nh_k_current[i]) > 0.0)
					{
						double tmp_err_h = fabs(carrier_scattering_prev_nh[i] - nh_k_current[i])/fabs(carrier_scattering_prev_nh[i]+nh_k_current[i]);
						if (tmp_err_h > CARRIER_DEVIATION_TOL)
						{
							err_h_max = tmp_err_h;
						}
					} else {
						err_h_max = CARRIER_DEVIATION_TOL*2.0;
					}
				}
			}

			if ((err_e_max > CARRIER_DEVIATION_TOL)||(err_h_max > CARRIER_DEVIATION_TOL))
			{
				calc = true;

				// Store previous carriers for error evaluation
				memcpy(carrier_scattering_prev_ne, ne_k_current, number_K_points*sizeof(double));
				memcpy(carrier_scattering_prev_nh, nh_k_current, number_K_points*sizeof(double));
			}
		} else {
			calc = true;

			// Store previous carriers for error evaluation
			memcpy(carrier_scattering_prev_ne, ne_k_current, number_K_points*sizeof(double));
			memcpy(carrier_scattering_prev_nh, nh_k_current, number_K_points*sizeof(double));
		}

		if (device_density > 1.0e16)
		{
		#if defined(USE_ISAK_HOLE_FILLING)
		
			if (calc == true)
			{
				carrier_scattering_method2_parallel(ne_k_current, nh_k_current);
			} else {
				// Approximate behavior based on previous timestep 
				carrier_scattering_method2_parallel_constant_rate(ne_k_current, nh_k_current);
			}
			
		
		#elif defined(USE_ISAK_HOLE_FILLING_TABLE)
		
			if (calc == true)
			{
				// Find instantanious temp
				carrier_scattering_rate_set_instant_temperature(ne_k_current, nh_k_current);

				// Only updates interpolation if calc==true
				carrier_scattering_method2_table(ne_k_current, nh_k_current);
			}
		
		#else
			if (calc == true)
			{
				carrier_scattering_rate_approx(ne_k_current, nh_k_current);
			}
		
		#endif
		} else {
			carrier_scattering_rate_approx(ne_k_current, nh_k_current);
		}
	}

//	#ifdef USE_DEVICE_TIMERS
//	MainStat->stop(getName(),"c-c scatt");
//	#endif

}

/* Initialization of higher order carrier scattering terms
 * Initialize :
 * 1. Angular grid points
 * 2. Momentum grids Q, K'
 * 3. Coulomb Screening function and index set for fast screening
 * 4. Confinement function
 * 5. Screened Coulomb potential
 * 6. Index sets for fast calculation of carrier scattering.
 * 
 * Using: V_col(q) = (e*e/2*eps0*epsRel)*1/(q)
 * Using: F_conf = sqrt(2/Lqw)cos(z*pi/Lqw)
 * 
 * Assumptions are that ee,eh,he,hh have identical confinement functions
 * and thus equal coulomb potentials
 * */
void Device::carrier_scattering_init(double *ne_k_old, double *nh_k_old)
{
	//===============================
	// Initialize the 2nd Born arrays	

	// Carrier-Carrier scattering
	// GC+CS: (N_Q=Nk/4, N_k'=Nk/4, N_th = 20) ~3 % average relative error, all errors < 10% For tests with Nk=100, Fermi+hole, Loretizan at k=6
	// Linear interpolation gives similar result but ~10% average relative error for test Fermi+hole, Loretizan
	// GC+CS: (N_Q=Nk/8, N_k'=Nk/8, N_th = 10) ~16 % average relative error, all errors < 100% For tests with Nk=100, Fermi+hole, Loretizan at k=6
	// Linear interpolation gives similar result but ~22% average relative error for test Fermi+hole, Loretizan
	number_Q_points		= 10;
	number_Kp_points	= 30;
	number_Th_points	= 6; // Angular integration, carrier-carrier scattering

	// Phonon-carrier scattering
	// GC+CS: (N_Q=10, N_Phi=20) ~0.14% average relative error, all errors < 1 % For tests with Nk=100, Fermi, Lorentz, 3000K Fermi -> 300K lattice
	// Similar results for linear interpolation, but error grows exponentially (~10^6) when interpolation is bad such as for Lorenzian or Fermi+hole
	number_Qph_points	= floor(number_K_points/20.0); // Phonon integral resolution
	number_Phi_points	= floor(number_K_points/20.0); // Angular integration, phonons

	carrier_scattering_phonon_eps_zero = 13.46979; // dielectric constant at low frequency
	carrier_scattering_phonon_eps_inf  = 11.23738; // dielectric constant at high frequency
	carrier_scattering_phonon_hbar_wLO = (33.95/1000.0)*(e); // [LO-Phonon-energy ]


	//================== CC scattering =================================
	dQ = new double[number_Q_points];
	Q = new double[number_Q_points];
	dKp = new double[number_Kp_points];
	Kp = new double[number_Kp_points];
	
	
	// Initialize GC grids
	// Q: Type 1 Gauss-Chebychev Grid (a,b)
	for(int i = 1; i <= number_Q_points; i++)
	{
		// Type 1 GC grid
		double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Q_points));
		double wi = Pi/((double)number_Q_points);
		double omega_zi = 1.0/sqrt(1.0-zi*zi);
		
		double q_a = K[0];
		double q_b = K[number_K_points-1];
		Q[number_Q_points-i]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
		dQ[number_Q_points-i] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
	}

	// Calculate Kp GRID (0,K)
	// K': Type 1 Gauss-Chebychev Grid (a,b)
	for(int i = 1; i <= number_Kp_points; i++)
	{
		// Type 1 GC grid
		double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Kp_points));
		double wi = Pi/((double)number_Kp_points);
		double omega_zi = 1.0/sqrt(1.0-zi*zi);
		
		double kp_a = 0.0;
		double kp_b = K_max;
		Kp[number_Kp_points-i]  = kp_a +  0.5*(kp_b-kp_a)*(zi+1.0);
		dKp[number_Kp_points-i] = 0.5*(kp_b-kp_a)*wi*(1.0/omega_zi);
	}
	if (number_Th_points % 2 != 0)
	{
		number_Th_points += 1; // make even
		//cout << "2nd Born initialization: angular grid should have even number of points.." << endl; 
		//cout << "number_Th_points = " << number_Th_points << endl;
		//exit(-1);
	}
	theta_2B  = new double[number_Th_points];
	dtheta_2B = new double[number_Th_points];

	for(int i = 1; i <= number_Th_points; i++)
	{
		// Type 1 GC grid
		double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Th_points));
		double wi = Pi/((double)number_Th_points);
		double omega_zi = 1.0/sqrt(1.0-zi*zi);
		
		double phi_a = 0.0; 
		double phi_b = Pi;
		theta_2B[number_Th_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
		dtheta_2B[number_Th_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
	}
	
	//================== CP scattering =================================
	
	Q_ph  = new double[number_Qph_points];
	dQ_ph = new double[number_Qph_points];
	// Initialize GC grids
	// Q: Type 1 Gauss-Chebychev Grid (a,b)
	for(int i = 1; i <= number_Qph_points; i++)
	{
		// Type 1 GC grid
		double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
		double wi = Pi/((double)number_Qph_points);
		double omega_zi = 1.0/sqrt(1.0-zi*zi);
		
		double q_a = 0.0;
		double q_b = K_max;
		Q_ph[number_Qph_points-i]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
		dQ_ph[number_Qph_points-i] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
	}
	
	if (number_Phi_points % 2 != 0)
	{
		number_Phi_points += 1; // make even
		//cout << "2nd Born initialization: angular grid should have even number of points.." << endl; 
		//cout << "number_Phi_points = " << number_Phi_points << endl;
		//exit(-1);
	}
	Phi  = new double[number_Phi_points];
	dPhi = new double[number_Phi_points];
	
	double phi_a = 0.0; 
	double phi_b = Pi/2.0;
	for(int i = 1; i <= number_Phi_points; i++)
	{
		// Type 1 GC grid
		double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Phi_points));
		double wi = Pi/((double)number_Phi_points);
		double omega_zi = 1.0/sqrt(1.0-zi*zi);
		
		Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
		dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
	}

	//=========== Start carrier-carrier scattering ========================
	double mass_hole_electron = device_mh/device_me; // Rato of the mass of hole / electron
	double mass_electron_hole = device_me/device_mh; // Rato of the mass of electron / hole
	
	double delta_scale_e = device_me*a0*a0/(hbar*hbar); // Factors scaled out of delta(E)
	double delta_scale_h = device_mh*a0*a0/(hbar*hbar); // Factors scaled out of delta(E)
	double col_scale = a0*e*e/(2.0*eps0*eps); // Factor scaled out of Vq
	double volume_scale = (1.0/(16.0*Pi*Pi*Pi*Pi*a0*a0*a0*a0)); // integration k'dk' qdq dth
	double common_constants = (Pi/hbar)*col_scale*col_scale*volume_scale;


	double theta_zi_array[number_Th_points];
	double theta_wi_array[number_Th_points];
	double theta_omega_zi_array[number_Th_points];
	for(int i = 1; i <= number_Th_points; i++)
	{
		// Type 1 GC grid
		theta_zi_array[i-1]		= cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Th_points));
		theta_wi_array[i-1]		= Pi/((double)number_Th_points);
		theta_omega_zi_array[i-1]	= 1.0/sqrt(1.0-theta_zi_array[i-1]*theta_zi_array[i-1]);
	}
	
	//=====================================================================
	// Precalculate carrier scattering index set for faster 2B evaluation
	 
	 // Initialize counters
	int **MAX_LENGTH_ee_k_kp = new int*[number_K_points];
	for(int k = 0; k < number_K_points; k++)
	{
		MAX_LENGTH_ee_k_kp[k] = new int[number_Kp_points];
		for(int kp = 0; kp < number_Kp_points; kp++)
		{
			MAX_LENGTH_ee_k_kp[k][kp] = 0;
		}
	}
	
	int **MAX_LENGTH_eh_k_kp = new int*[number_K_points];
	for(int k = 0; k < number_K_points; k++)
	{
		MAX_LENGTH_eh_k_kp[k] = new int[number_Kp_points];
		for(int kp = 0; kp < number_Kp_points; kp++)
		{
			MAX_LENGTH_eh_k_kp[k][kp] = 0;
		}
	}
	
	int **MAX_LENGTH_he_k_kp = new int*[number_K_points];
	for(int k = 0; k < number_K_points; k++)
	{
		MAX_LENGTH_he_k_kp[k] = new int[number_Kp_points];
		for(int kp = 0; kp < number_Kp_points; kp++)
		{
			MAX_LENGTH_he_k_kp[k][kp] = 0;
		}
	}
	 
	int MAX_LENGTH_ee = 0;
	int MAX_LENGTH_eh = 0;
	int MAX_LENGTH_he = 0;

	for(int k = 0; k < number_K_points; k++)
	{
		
		// Constants used in boundaries of eh/he integration
		double eh_const1 = 2.0*K[k]*mass_hole_electron/(1.0+mass_hole_electron);
		double eh_const2 = 2.0/(1.0+mass_hole_electron);
		double he_const1 = 2.0*K[k]/(1.0+mass_hole_electron);
		double he_const2 = 2.0*mass_hole_electron/(1.0+mass_hole_electron);
		
		
		// ========================= ee and hh scattering ==================================
		for(int kp = 0; kp < number_Kp_points; kp++) // Momentum magnitude kp
		{
			// Calculate Angular grid
			double th_a = 0.0; // Using half the angular interval. Multiply solution by 2.0
			double th_b = Pi;
			if (Kp[kp] > K[k])
			{
				// whole range
				th_a = 0.0;
				th_b = Pi;
			} else {
				
				th_a = acos(Kp[kp]/K[k]);
				th_b = Pi-th_a;
			}
			
			for(int i = 1; i <= number_Th_points; i++)
			{
				// Type 1 GC grid
				double zi = theta_zi_array[i-1];
				double wi = theta_wi_array[i-1];
				double omega_zi = theta_omega_zi_array[i-1];
				
				theta_2B[number_Th_points-i]  = th_a +  0.5*(th_b-th_a)*(zi+1.0);
				dtheta_2B[number_Th_points-i] = 0.5*(th_b-th_a)*wi*(1.0/omega_zi);
			}
			
			for(int q = 0; q < number_Q_points; q++) // Momentum magnitude q
			{
				// e-e and h-h scattering
				for(int th_k = 0; th_k < number_Th_points; th_k++) // Angle between k and q
				{
					double angle_th_k_q  	= theta_2B[th_k];
					double cos_k_q  		= cos(theta_2B[th_k]);

					double acos_arg = cos_k_q*K[k]/Kp[kp];
					if (fabs(acos_arg) < 1.0)
					{
						// Removed branch by multiplication of 2.0 below
						double angle_th_kp_q    = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]
						double angle_th_k_kp    = angle_th_kp_q-angle_th_k_q; // Angle is between [-pi,2*pi]
						if (angle_th_k_kp > Pi)
						{
							angle_th_k_kp = 2.0*Pi - angle_th_k_kp; // Angle between [0,pi]
						}
						// For speed
						double cos_kp_q         = cos(angle_th_kp_q);
						double cos_k_kp    		= cos(angle_th_k_kp);
						double k_q              = sqrt(K[k ]*K[k ]   + Q[q]*Q[q] -2.0*K[k ] *Q[q]*cos_k_q);
						double kp_q             = sqrt(Kp[kp]*Kp[kp] + Q[q]*Q[q] -2.0*Kp[kp]*Q[q]*cos_kp_q);
						double k_kp				= sqrt(K[k]*K[k]+Kp[kp]*Kp[kp] - 2.0*K[k]*Kp[kp]*cos_k_kp); // Units [-]
						
						if (((((k_q <= K[number_K_points-1])&&(kp_q <= K[number_K_points-1]))&&(k_kp <= K[number_K_points-1]))&&(((k_q >= K[0])&&(kp_q >= K[0]))&&(k_kp >= K[0])))&&((Kp[kp]>=K[0])&&(Kp[kp]<=K[number_K_points-1])))
						{
							MAX_LENGTH_ee++; // Total number of indices
							MAX_LENGTH_ee_k_kp[k][kp]++; // Number of indices sorted by kp
						}
					} else {
						cout << "problem with e-e grid.." << endl;
						exit(-1);
					}
				}
			}
			
			// =============== eh scattering ===================
			for(int q = 0; q < number_Q_points; q++) // Momentum magnitude q
			{
				// Calculate Angular grid
				if ((eh_const1 < Q[q])&&( Q[q]-eh_const1 >= Kp[kp]*eh_const2)) // Check if there are ANY possible angles
				{
					// No angles available
				} else {
				
					double th_a = 0.0; // If Using half the angular interval. Multiply solution by 2.0
					double th_b = Pi;
					if (Kp[kp]*eh_const2 > Q[q]+eh_const1)
					{
						// Use whole range
						th_a = 0.0;
						th_b = Pi;
					} else {
						
						// There is at least one maximal angle
						th_b = acos((Q[q]-Kp[kp]*eh_const2)/eh_const1);
						
						// There can also be a minimum angle > 0.0
						if (Q[q]+Kp[kp]*eh_const2 < eh_const1)
						{
							th_a = acos((Q[q]+Kp[kp]*eh_const2)/eh_const1);
						}
					}
					
					for(int i = 1; i <= number_Th_points; i++)
					{
						// Type 1 GC grid
						double zi = theta_zi_array[i-1];
						double wi = theta_wi_array[i-1];
						double omega_zi = theta_omega_zi_array[i-1];
						
						theta_2B[number_Th_points-i]  = th_a +  0.5*(th_b-th_a)*(zi+1.0);
						dtheta_2B[number_Th_points-i] = 0.5*(th_b-th_a)*wi*(1.0/omega_zi);
					}
					
					
					// e-h and h-e scattering
					for(int th_k = 0; th_k < number_Th_points; th_k++) // Angle between k and q
					{
						double angle_th_k_q = theta_2B[th_k];
						double cos_k_q      = cos(theta_2B[th_k]);
						double cos_arg_eh   = -0.5*(Q[q]/Kp[kp])*(1.0+mass_hole_electron)+(mass_hole_electron)*(K[k]/Kp[kp])*cos_k_q;
						
						if (fabs(cos_arg_eh)<1.0)
						{
							// Removed branch by multiplication of 2.0 below
							double angle_th_kp_q  = acos(cos_arg_eh);
							double angle_th_k_kp  = angle_th_kp_q-angle_th_k_q; // Cos(th) is symetric
							
							// For speed
							double cos_kp_q = cos(angle_th_kp_q);
							double k_q  = sqrt(K[k ]*K[k ]   + Q[q]*Q[q] - 2.0*K[k ]*Q[q]*cos_k_q);
							double kp_q = sqrt(Kp[kp]*Kp[kp] + Q[q]*Q[q] + 2.0*Kp[kp]*Q[q]*cos_kp_q);
							
							if ((((k_q <= K[number_K_points-1])&&(kp_q <= K[number_K_points-1]))&&((k_q >= K[0])&&(kp_q >= K[0])))&&((Kp[kp]>=K[0])&&(Kp[kp]<=K[number_K_points-1])))
							{
								MAX_LENGTH_eh++;
								MAX_LENGTH_eh_k_kp[k][kp]++; // Number of indices sorted by kp
							}
						} else {
							cout << "problem with e-h grid.." << endl;
							exit(-1);
						}
					}
				} // End eh-angles
			}
			
			// =============== he scattering ===================
			for(int q = 0; q < number_Q_points; q++) // Momentum magnitude q
			{
				// Calculate Angular grid
				if ((he_const1 < Q[q])&&( Q[q]-he_const1 >= Kp[kp]*he_const2)) // Check if there are ANY possible angles
				{
					// No angles available
				} else {
				
					double th_a = 0.0; // If Using half the angular interval. Multiply solution by 2.0
					double th_b = Pi;
					if (Kp[kp]*he_const2 > Q[q]+he_const1)
					{
						// Use whole range
						th_a = 0.0;
						th_b = Pi;
					} else {
						
						// There is at least one maximal angle
						th_b = acos((Q[q]-Kp[kp]*he_const2)/he_const1);
						
						// There can also be a minimum angle > 0.0
						if (Q[q]+Kp[kp]*he_const2 < he_const1)
						{
							th_a = acos((Q[q]+Kp[kp]*he_const2)/he_const1);
						}
					}
					
					for(int i = 1; i <= number_Th_points; i++)
					{
						// Type 1 GC grid
						double zi = theta_zi_array[i-1];
						double wi = theta_wi_array[i-1];
						double omega_zi = theta_omega_zi_array[i-1];
						
						theta_2B[number_Th_points-i]  = th_a +  0.5*(th_b-th_a)*(zi+1.0);
						dtheta_2B[number_Th_points-i] = 0.5*(th_b-th_a)*wi*(1.0/omega_zi);
					}
					
					// e-h and h-e scattering
					for(int th_k = 0; th_k < number_Th_points; th_k++) // Angle between k and q
					{
						double angle_th_k_q = theta_2B[th_k];
						double cos_k_q      = cos(theta_2B[th_k]);
						double cos_arg_he = -0.5*(Q[q]/Kp[kp])*(1.0+mass_electron_hole)+(mass_electron_hole)*(K[k]/Kp[kp])*cos_k_q;
						
						if (fabs(cos_arg_he)<1.0)
						{

							double angle_th_kp_q  = acos(cos_arg_he);
							double angle_th_k_kp  = angle_th_kp_q-angle_th_k_q; // Cos(th) is symetric

							// For speed
							double cos_kp_q = cos(angle_th_kp_q);
							double k_q  = sqrt(K[k ]*K[k ]   + Q[q]*Q[q] - 2.0*K[k ]*Q[q]*cos_k_q);
							double kp_q = sqrt(Kp[kp]*Kp[kp] + Q[q]*Q[q] + 2.0*Kp[kp]*Q[q]*cos_kp_q);
							
							if ((((k_q <= K[number_K_points-1])&&(kp_q <= K[number_K_points-1]))&&((k_q >= K[0])&&(kp_q >= K[0])))&&((Kp[kp]>=K[0])&&(Kp[kp]<=K[number_K_points-1])))
							{
								MAX_LENGTH_he++;
								MAX_LENGTH_he_k_kp[k][kp]++; // Number of indices sorted by kp
							}
						} else {
							cout << "problem with h-e grid.." << endl;
							exit(-1);
						}
					}
				} // End he-angles
			}
		}
	}
	
	cout << "carrier scattering # ee index = " << MAX_LENGTH_ee << endl;
	cout << "carrier scattering # eh index = " << MAX_LENGTH_eh << endl;
	cout << "carrier scattering # he index = " << MAX_LENGTH_he << endl;


	// Prepare each index set
	if (carrier_scattering_index_ee == NULL)
	{
		carrier_scattering_index_ee = new indexSet3d_cc(number_K_points, number_Kp_points, MAX_LENGTH_ee_k_kp);
	}
	if (carrier_scattering_index_hh == NULL)
	{
		carrier_scattering_index_hh = new indexSet3d_cc(number_K_points, number_Kp_points, MAX_LENGTH_ee_k_kp);
	}
	if (carrier_scattering_index_eh == NULL)
	{
		carrier_scattering_index_eh = new indexSet3d_cc(number_K_points, number_Kp_points, MAX_LENGTH_eh_k_kp);
	}
	if (carrier_scattering_index_he == NULL)
	{
		carrier_scattering_index_he = new indexSet3d_cc(number_K_points, number_Kp_points, MAX_LENGTH_he_k_kp);
	}


	for(int k = 0; k < number_K_points; k++)
	{
		// Constants used in boundaries of eh/he integration
		double eh_const1 = 2.0*K[k]*mass_hole_electron/(1.0+mass_hole_electron);
		double eh_const2 = 2.0/(1.0+mass_hole_electron);
		double he_const1 = 2.0*K[k]/(1.0+mass_hole_electron);
		double he_const2 = 2.0*mass_hole_electron/(1.0+mass_hole_electron);
		
		
		// ========================= ee and hh scattering ==================================
		for(int kp = 0; kp < number_Kp_points; kp++) // Momentum magnitude kp
		{
			int count_k_kp_ee = 0;
			int count_k_kp_eh = 0;
			int count_k_kp_he = 0;

			// Calculate Angular grid
			double th_a = 0.0; // Using half the angular interval. Multiply solution by 2.0
			double th_b = Pi;
			if (Kp[kp] > K[k])
			{
				// whole range
				th_a = 0.0;
				th_b = Pi;
			} else {
				
				th_a = acos(Kp[kp]/K[k]);
				th_b = Pi-th_a;
			}
			
			for(int i = 1; i <= number_Th_points; i++)
			{
				// Type 1 GC grid
				double zi = theta_zi_array[i-1];
				double wi = theta_wi_array[i-1];
				double omega_zi = theta_omega_zi_array[i-1];
				
				theta_2B[number_Th_points-1-(i-1)]  = th_a +  0.5*(th_b-th_a)*(zi+1.0);
				dtheta_2B[number_Th_points-1-(i-1)] = 0.5*(th_b-th_a)*wi*(1.0/omega_zi);
			}
			
			for(int q = 0; q < number_Q_points; q++) // Momentum magnitude q
			{
				// V[q] calculation
				double Ve_q  = misc_interpolate_K_array_parabolic(Q[q], K, coulomb_potential_normalized_ee, number_K_points);
				double Vh_q  = misc_interpolate_K_array_parabolic(Q[q], K, coulomb_potential_normalized_hh, number_K_points);
				double Veh_q = misc_interpolate_K_array_parabolic(Q[q], K, coulomb_potential_normalized_eh, number_K_points);
				double Veh_qVeh_q = Veh_q*Veh_q;
				
				
				// e-e and h-h scattering
				for(int th_k = 0; th_k < number_Th_points; th_k++) // Angle between k and q
				{
					double angle_th_k_q  	= theta_2B[th_k];
					double cos_k_q  		= cos(theta_2B[th_k]);

					double acos_arg = cos_k_q*K[k]/Kp[kp];
					if (fabs(acos_arg) < 1.0)
					{
						// Removed branch by multiplication of 2.0 below
						double angle_th_kp_q    = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]
						double angle_th_k_kp    = angle_th_kp_q-angle_th_k_q; // Angle is between [-pi,2*pi]
						if (angle_th_k_kp > Pi)
						{
							angle_th_k_kp = 2.0*Pi - angle_th_k_kp; // Angle between [0,pi]
						//} else if (angle_th_k_kp < -Pi)
						//{
						//	angle_th_k_kp = -2.0*Pi - angle_th_k_kp; // Angle between [-pi,0]
						//	cout << "Strange error.." << endl;
						//	exit(-1);
						}
						// For speed
						double cos_kp_q         = acos_arg;
						double cos_k_kp    	= cos(angle_th_k_kp);
						double k_q              = sqrt(K[k ]*K[k ]   + Q[q]*Q[q] -2.0*K[k ] *Q[q]*cos_k_q);
						double kp_q             = sqrt(Kp[kp]*Kp[kp] + Q[q]*Q[q] -2.0*Kp[kp]*Q[q]*cos_kp_q);
						double k_kp		= sqrt(K[k]*K[k]+Kp[kp]*Kp[kp] - 2.0*K[k]*Kp[kp]*cos_k_kp);
						
						if (((((k_q <= K[number_K_points-1])&&(kp_q <= K[number_K_points-1]))&&(k_kp <= K[number_K_points-1]))&&(((k_q >= K[0])&&(kp_q >= K[0]))&&(k_kp >= K[0])))&&((Kp[kp]>=K[0])&&(Kp[kp]<=K[number_K_points-1])))
						{
							// V[k-kp] calculation
							double Ve_k_kp = misc_interpolate_K_array_parabolic(k_kp, K, coulomb_potential_normalized_ee, number_K_points);
							double Vh_k_kp = misc_interpolate_K_array_parabolic(k_kp, K, coulomb_potential_normalized_hh, number_K_points);
							//double Ve_k_kp = 0.0;
							//double Vh_k_kp = 0.0;
							
							//==========================
							//      e-e and h-h scattering 
							//==========================
							// Energy difference
							double DF_ee   =  1.0/(Kp[kp]*Q[q]*fabs(sin(angle_th_kp_q)));
							double DF_hh   =  1.0/(Kp[kp]*Q[q]*fabs(sin(angle_th_kp_q)));
							
							double beta_kp, beta_k_q ,beta_kp_q;
							int ind_kp, ind_k_q,ind_kp_q;

							// Occupation number calculation
							double fe_kp    = misc_interpolate_K_array_linear_index(Kp[kp], K, ne_k_old, number_K_points, &ind_kp  , &beta_kp);
							double fe_k_q   = misc_interpolate_K_array_linear_index(k_q   , K, ne_k_old, number_K_points, &ind_k_q , &beta_k_q);
							double fe_kp_q  = misc_interpolate_K_array_linear_index(kp_q  , K, ne_k_old, number_K_points, &ind_kp_q, &beta_kp_q);
							
							// Integral volume constants
							double integral_volume = Kp[kp]*dKp[kp]*Q[q]*dQ[q]*dtheta_2B[th_k]; // Constants from q*dq*kp*dkp*dt1*dt2
							
							double Fe  = 3.5*integral_volume*4.0*Ve_q*(Ve_q - Ve_k_kp)*DF_ee*delta_scale_e*common_constants;
							double Fh  = 3.5*integral_volume*4.0*Vh_q*(Vh_q - Vh_k_kp)*DF_hh*delta_scale_h*common_constants;

							carrier_scattering_index_ee->updateWeight(k, kp, count_k_kp_ee, Fe, ind_kp, beta_kp, ind_k_q, beta_k_q, ind_kp_q, beta_kp_q);
							carrier_scattering_index_hh->updateWeight(k, kp, count_k_kp_ee, Fh, ind_kp, beta_kp, ind_k_q, beta_k_q, ind_kp_q, beta_kp_q);

							count_k_kp_ee += 1;
						}
					} else {
						cout << "problem with e-e grid.." << endl;
						exit(-1);
					}
				}
			}
			
			// =============== eh scattering ===================
			for(int q = 0; q < number_Q_points; q++) // Momentum magnitude q
			{
				// V[q] calculation
				double Veh_q = misc_interpolate_K_array_parabolic(Q[q], K, coulomb_potential_normalized_eh, number_K_points);
				double Veh_qVeh_q = Veh_q*Veh_q;
				
				// Calculate Angular grid
				if ((eh_const1 < Q[q])&&( Q[q]-eh_const1 >= Kp[kp]*eh_const2)) // Check if there are ANY possible angles
				{
					// No angles available
				} else {
				
					double th_a = 0.0; // If Using half the angular interval. Multiply solution by 2.0
					double th_b = Pi;
					if (Kp[kp]*eh_const2 > Q[q]+eh_const1)
					{
						// Use whole range
						th_a = 0.0;
						th_b = Pi;
					} else {
						
						// There is at least one maximal angle
						th_b = acos((Q[q]-Kp[kp]*eh_const2)/eh_const1);
						
						// There can also be a minimum angle > 0.0
						if (Q[q]+Kp[kp]*eh_const2 < eh_const1)
						{
							th_a = acos((Q[q]+Kp[kp]*eh_const2)/eh_const1);
						}
					}
					
					for(int i = 1; i <= number_Th_points; i++)
					{
						// Type 1 GC grid
						double zi = theta_zi_array[i-1];
						double wi = theta_wi_array[i-1];
						double omega_zi = theta_omega_zi_array[i-1];
						
						theta_2B[number_Th_points-1-(i-1)]  = th_a +  0.5*(th_b-th_a)*(zi+1.0);
						dtheta_2B[number_Th_points-1-(i-1)] = 0.5*(th_b-th_a)*wi*(1.0/omega_zi);
					}
					
					
					// e-h scattering
					for(int th_k = 0; th_k < number_Th_points; th_k++) // Angle between k and q
					{
						double angle_th_k_q = theta_2B[th_k];
						double cos_k_q      = cos(theta_2B[th_k]);
						double cos_arg_eh   = -0.5*(Q[q]/Kp[kp])*(1.0+mass_hole_electron)+(mass_hole_electron)*(K[k]/Kp[kp])*cos_k_q;
						
						if (fabs(cos_arg_eh)<1.0)
						{
							// Removed branch by multiplication of 2.0 below
							double angle_th_kp_q  = acos(cos_arg_eh);
							double angle_th_k_kp  = angle_th_kp_q-angle_th_k_q; // Cos(th) is symetric
							
							double cos_kp_q = cos_arg_eh;
							double k_q  = sqrt(K[k ]*K[k ]   + Q[q]*Q[q] - 2.0*K[k ]*Q[q]*cos_k_q);
							double kp_q = sqrt(Kp[kp]*Kp[kp] + Q[q]*Q[q] + 2.0*Kp[kp]*Q[q]*cos_kp_q);
							
							if ((((k_q <= K[number_K_points-1])&&(kp_q <= K[number_K_points-1]))&&((k_q >= K[0])&&(kp_q >= K[0])))&&((Kp[kp]>=K[0])&&(Kp[kp]<=K[number_K_points-1])))
							{
								//===========================
								//      e-h scattering 
								//===========================
								// Energy difference
								double DF_eh   =  1.0/(Kp[kp]*Q[q]*fabs(sin(angle_th_kp_q)));

								double beta_kp, beta_k_q ,beta_kp_q;
								int ind_kp, ind_k_q,ind_kp_q;

								// Occupation number calculation
								double fe_k_q  = misc_interpolate_K_array_linear_index(k_q   , K, ne_k_old, number_K_points, &ind_k_q , &beta_k_q);
								double fh_kp   = misc_interpolate_K_array_linear_index(Kp[kp], K, nh_k_old, number_K_points, &ind_kp  , &beta_kp);
								double fh_kp_q = misc_interpolate_K_array_linear_index(kp_q  , K, nh_k_old, number_K_points, &ind_kp_q, &beta_kp_q);

								// Integral volume constants
								double integral_volume = Kp[kp]*dKp[kp]*Q[q]*dQ[q]*dtheta_2B[th_k]; // Constants from q*dq*kp*dkp*dt1*dt2
								double F  = 3.5*integral_volume*4.0*Veh_qVeh_q*DF_eh*delta_scale_h*common_constants;

								carrier_scattering_index_eh->updateWeight(k, kp, count_k_kp_eh, F, ind_kp, beta_kp, ind_k_q, beta_k_q, ind_kp_q, beta_kp_q);
								count_k_kp_eh += 1;
							}
						} else {
							cout << "problem with e-h grid.." << endl;
							exit(-1);
						}
					}
				} // End eh-angles
			}
			
			// =============== he scattering ===================
			for(int q = 0; q < number_Q_points; q++) // Momentum magnitude q
			{
				// V[q] calculation
				double Veh_q = misc_interpolate_K_array_parabolic(Q[q], K, coulomb_potential_normalized_eh, number_K_points);
				double Veh_qVeh_q = Veh_q*Veh_q;
				
				// Calculate Angular grid
				if ((he_const1 < Q[q])&&( Q[q]-he_const1 >= Kp[kp]*he_const2)) // Check if there are ANY possible angles
				{
					// No angles available
				} else {
				
					double th_a = 0.0; // If Using half the angular interval. Multiply solution by 2.0
					double th_b = Pi;
					if (Kp[kp]*he_const2 > Q[q]+he_const1)
					{
						// Use whole range
						th_a = 0.0;
						th_b = Pi;
					} else {
						
						// There is at least one maximal angle
						th_b = acos((Q[q]-Kp[kp]*he_const2)/he_const1);
						
						// There can also be a minimum angle > 0.0
						if (Q[q]+Kp[kp]*he_const2 < he_const1)
						{
							th_a = acos((Q[q]+Kp[kp]*he_const2)/he_const1);
						}
					}
					
					for(int i = 1; i <= number_Th_points; i++)
					{
						// Type 1 GC grid
						double zi = theta_zi_array[i-1];
						double wi = theta_wi_array[i-1];
						double omega_zi = theta_omega_zi_array[i-1];
						
						theta_2B[number_Th_points-1-(i-1)]  = th_a +  0.5*(th_b-th_a)*(zi+1.0);
						dtheta_2B[number_Th_points-1-(i-1)] = 0.5*(th_b-th_a)*wi*(1.0/omega_zi);
					}
					
					// h-e scattering
					for(int th_k = 0; th_k < number_Th_points; th_k++) // Angle between k and q
					{
						double angle_th_k_q = theta_2B[th_k];
						double cos_k_q      = cos(theta_2B[th_k]);
						double cos_arg_he = -0.5*(Q[q]/Kp[kp])*(1.0+mass_electron_hole)+(mass_electron_hole)*(K[k]/Kp[kp])*cos_k_q;
						
						if (fabs(cos_arg_he)<1.0)
						{

							double angle_th_kp_q  = acos(cos_arg_he);
							double angle_th_k_kp  = angle_th_kp_q-angle_th_k_q; // Cos(th) is symetric

							double cos_kp_q = cos_arg_he;
							double k_q  = sqrt(K[k ]*K[k ]   + Q[q]*Q[q] - 2.0*K[k ]*Q[q]*cos_k_q);
							double kp_q = sqrt(Kp[kp]*Kp[kp] + Q[q]*Q[q] + 2.0*Kp[kp]*Q[q]*cos_kp_q);
							
							if ((((k_q <= K[number_K_points-1])&&(kp_q <= K[number_K_points-1]))&&((k_q >= K[0])&&(kp_q >= K[0])))&&((Kp[kp]>=K[0])&&(Kp[kp]<=K[number_K_points-1])))
							{
								//===========================
								//      h-e scattering 
								//===========================
								// Energy difference
								double DF_he   =  1.0/(Kp[kp]*Q[q]*fabs(sin(angle_th_kp_q)));

								double beta_kp, beta_k_q ,beta_kp_q;
								int ind_kp, ind_k_q,ind_kp_q;

								// Occupation number calculation
								double fh_k_q  = misc_interpolate_K_array_linear_index(k_q   , K, nh_k_old, number_K_points, &ind_k_q , &beta_k_q);
								double fe_kp   = misc_interpolate_K_array_linear_index(Kp[kp], K, ne_k_old, number_K_points, &ind_kp  , &beta_kp);
								double fe_kp_q = misc_interpolate_K_array_linear_index(kp_q  , K, ne_k_old, number_K_points, &ind_kp_q, &beta_kp_q);
								
								// Integral volume constants
								double integral_volume = Kp[kp]*dKp[kp]*Q[q]*dQ[q]*dtheta_2B[th_k]; // Constants from q*dq*kp*dkp*dt1*dt2
								double F  = 2.5*integral_volume*4.0*Veh_qVeh_q*DF_he*delta_scale_e*common_constants;
								
								carrier_scattering_index_he->updateWeight(k, kp, count_k_kp_he, F, ind_kp, beta_kp, ind_k_q, beta_k_q, ind_kp_q, beta_kp_q);
								count_k_kp_he += 1;
							}
						} else {
							cout << "problem with h-e grid.." << endl;
							exit(-1);
						}
					}
				} // End he-angles
			}
		}
	}
	// Sort arrays
	carrier_scattering_index_ee->sortIndices();
	carrier_scattering_index_hh->sortIndices();
	carrier_scattering_index_eh->sortIndices();
	carrier_scattering_index_he->sortIndices();
	

	// ============================
	// Phonon scattering
	double hbar_wLO = carrier_scattering_phonon_hbar_wLO;
	carrier_scattering_phonon_temp = device_background_temp_e; // Lattice temperature [K]
	carrier_scattering_phonon_n = 1.0/(exp(hbar_wLO/(kB*carrier_scattering_phonon_temp))-1.0);

	cout << "phonon: eps_zero  = " << carrier_scattering_phonon_eps_zero << endl;
	cout << "phonon: eps_inf   = " << carrier_scattering_phonon_eps_inf << endl;
	cout << "phonon: hbar*w_LO = " << 1000.0*carrier_scattering_phonon_hbar_wLO/e << " [meV]" << endl;
	cout << "phonon: temp      = " << carrier_scattering_phonon_temp << " [K]" << endl;
	cout << "phonon: n         = " << carrier_scattering_phonon_n << endl;


	double phi_zi_array[number_Phi_points];
	double phi_wi_array[number_Phi_points];
	double phi_omega_zi_array[number_Phi_points];
	for(int i = 1; i <= number_Phi_points; i++)
	{
		// Type 1 GC grid
		phi_zi_array[i-1]			= cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Phi_points));
		phi_wi_array[i-1]			= Pi/((double)number_Phi_points);
		phi_omega_zi_array[i-1]	= 1.0/sqrt(1.0-phi_zi_array[i-1]*phi_zi_array[i-1]);
	}
	
	double phonon_matrix_element_constants = 0.5*carrier_scattering_phonon_hbar_wLO*(e*e*a0*a0/eps0)*(1.0/carrier_scattering_phonon_eps_inf - 1.0/carrier_scattering_phonon_eps_zero);
	double common_constants_phonon = (Pi/hbar)*(1.0/(8.0*Pi*Pi*Pi*a0*a0*a0));
	double phonon_energy_level_e = 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar);
	double phonon_energy_level_h = 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar);

	
	 // Initialize counters
	int *MAX_LENGTH_e1 = new int[number_K_points];
	for(int k = 0; k < number_K_points; k++)
	{
		MAX_LENGTH_e1[k] = 0;
	}
	int *MAX_LENGTH_e2 = new int[number_K_points];
	for(int k = 0; k < number_K_points; k++)
	{
		MAX_LENGTH_e2[k] = 0;
	}
	int *MAX_LENGTH_h1 = new int[number_K_points];
	for(int k = 0; k < number_K_points; k++)
	{
		MAX_LENGTH_h1[k] = 0;
	}
	int *MAX_LENGTH_h2 = new int[number_K_points];
	for(int k = 0; k < number_K_points; k++)
	{
		MAX_LENGTH_h2[k] = 0;
	}

	for(int k = 0; k < number_K_points; k++)
	{
		//================== e dynamics ======================
		
		// Expression 1: (-Q_plane*Q_plane + 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane) < 1
		double phi_minus = asin((-K[k] + sqrt(K[k]*K[k] + phonon_energy_level_e))/K_max);
		double phi_pluss = Pi/2.0;
		
		double asin_arg_pluss = ( K[k] + sqrt(K[k]*K[k] + phonon_energy_level_e))/K_max;
		if (asin_arg_pluss < 1.0)
		{
			phi_pluss = asin(asin_arg_pluss);
		}
		
		// GRID 1
		double phi_a = phi_minus; // Using half the angular interval. Multiply solution by 2.0
		double phi_b = phi_pluss;
		
		for(int i = 1; i <= number_Phi_points; i++)
		{
			// Type 1 GC grid
			double zi = phi_zi_array[i-1];
			double wi = phi_wi_array[i-1];
			double omega_zi = phi_omega_zi_array[i-1];
			
			Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
			dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
		}
		
		for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
		{
			double q_minus = (-K[k] + sqrt(K[k]*K[k] + phonon_energy_level_e))/fabs(sin(Phi[phi_q]));
			//double q_pluss = ( K[k] + sqrt(K[k]*K[k] + phonon_energy_level_e))/fabs(sin(Phi[phi_q]));
			
			// Q: Type 1 Gauss-Chebychev Grid (a,b)
			double q_a = q_minus;
			double q_b = K_max;
			for(int i = 1; i <= number_Qph_points; i++)
			{
				// Type 1 GC grid
				double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
				double wi = Pi/((double)number_Qph_points);
				double omega_zi = 1.0/sqrt(1.0-zi*zi);
				
				
				Q_ph[number_Qph_points-1-(i-1)]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
				dQ_ph[number_Qph_points-1-(i-1)] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
			}
			
			for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
			{
				double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
				double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
				
				// First ee-phonon scattering
				double acos_arg = (-Q_plane*Q_plane + 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
				double f_acos_arg = fabs(acos_arg);
				if (f_acos_arg < 1.0)
				{
					// Removed branch by multiplication of 2.0 below
					double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]

					double cos_k_q  = acos_arg;
					double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane + 2.0*K[k ] *Q_plane*cos_k_q);

					if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
					{
						//=============================
						//      e-P scattering part 1
						//=============================
						MAX_LENGTH_e1[k] += 1;
					} 
				} else { // end acos_arg < 1
				
					for(int i = 0; i < number_Qph_points; i++)
					{
						cout << "q = " << Q_ph[q] << endl;
					}
					cout << "problems with limits 1e..." << endl;
					exit(-1);
				} 
			}
		}
		
		// GRID 2
		phi_a = phi_pluss; // Using half the angular interval. Multiply solution by 2.0
		phi_b = Pi/2.0;
		
		for(int i = 1; i <= number_Phi_points; i++)
		{
			// Type 1 GC grid
			double zi = phi_zi_array[i-1];
			double wi = phi_wi_array[i-1];
			double omega_zi = phi_omega_zi_array[i-1];
			
			Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
			dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
		}
		
		for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
		{
			double q_minus = (-K[k] + sqrt(K[k]*K[k] + phonon_energy_level_e))/fabs(sin(Phi[phi_q]));
			double q_pluss = ( K[k] + sqrt(K[k]*K[k] + phonon_energy_level_e))/fabs(sin(Phi[phi_q]));
			
			// Q: Type 1 Gauss-Chebychev Grid (a,b)
			double q_a = q_minus;
			double q_b = q_pluss;
			for(int i = 1; i <= number_Qph_points; i++)
			{
				// Type 1 GC grid
				double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
				double wi = Pi/((double)number_Qph_points);
				double omega_zi = 1.0/sqrt(1.0-zi*zi);
				
				
				Q_ph[number_Qph_points-i]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
				dQ_ph[number_Qph_points-i] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
			}
			
			
			for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
			{
				double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
				double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
				
				// First ee-phonon scattering
				double acos_arg = (-Q_plane*Q_plane + 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
				double f_acos_arg = fabs(acos_arg);
				if (f_acos_arg < 1.0)
				{
					// Removed branch by multiplication of 2.0 below
					double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]

					double cos_k_q  = acos_arg;
					double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane + 2.0*K[k ] *Q_plane*cos_k_q);

					if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
					{
						//=============================
						//      e-P scattering part 1
						//=============================
						MAX_LENGTH_e1[k] += 1;
					} 
				} else { // end acos_arg < 1
				
					for(int i = 0; i < number_Qph_points; i++)
					{
						cout << "q = " << Q_ph[q] << endl;
					}
					cout << "problems with limits 1e..." << endl;
					exit(-1);
				} 
			}
		}
		
		
		// Expression 2: (Q_plane*Q_plane + 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane) < 1
		if (K[k] > sqrt(phonon_energy_level_e))
		{
			double phi_minus = asin((K[k] - sqrt(K[k]*K[k] - phonon_energy_level_e))/K_max);
			double phi_pluss = Pi/2.0;
			
			double asin_arg_pluss = (K[k] + sqrt(K[k]*K[k] - phonon_energy_level_e))/K_max;
			if (asin_arg_pluss < 1.0)
			{
				phi_pluss = asin(asin_arg_pluss);
			}
			
			// GRID 1
			double phi_a = phi_minus; // Using half the angular interval. Multiply solution by 2.0
			double phi_b = phi_pluss;
			
			for(int i = 1; i <= number_Phi_points; i++)
			{
				// Type 1 GC grid
				double zi = phi_zi_array[i-1];
				double wi = phi_wi_array[i-1];
				double omega_zi = phi_omega_zi_array[i-1];
				
				Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
				dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
			}
			
			for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
			{
				double q_minus = (K[k] - sqrt(K[k]*K[k] - phonon_energy_level_e))/sin(Phi[phi_q]);
				//double q_pluss = (K[k] + sqrt(K[k]*K[k] - phonon_energy_level_e))/sin(Phi[phi_q]);
				
				// Q: Type 1 Gauss-Chebychev Grid (a,b)
				double q_a = q_minus;
				double q_b = K_max;
				for(int i = 1; i <= number_Qph_points; i++)
				{
					// Type 1 GC grid
					double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
					double wi = Pi/((double)number_Qph_points);
					double omega_zi = 1.0/sqrt(1.0-zi*zi);
					
					Q_ph[number_Qph_points-1-(i-1)]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
					dQ_ph[number_Qph_points-1-(i-1)] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
				}
				
				for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
				{
					double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
					double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
					
					// Second ee-phonon scattering
					double acos_arg = (Q_plane*Q_plane + 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
					double f_acos_arg = fabs(acos_arg);
					if (f_acos_arg < 1.0)
					{
						// Removed branch by multiplication of 2.0 below
						double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]

						double cos_k_q  = acos_arg;
						double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane - 2.0*K[k ] *Q_plane*cos_k_q);

						if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
						{
							//=============================
							//      e-P scattering part 2
							//=============================
							MAX_LENGTH_e2[k] += 1;
						} 
					} else { // end acos_arg < 1
						
						cout << "problems with limits 2e..." << endl;
						exit(-1);
					}
				}
			}
			
			// GRID 2
			phi_a = phi_pluss; // Using half the angular interval. Multiply solution by 2.0
			phi_b = Pi/2.0;
			
			for(int i = 1; i <= number_Phi_points; i++)
			{
				// Type 1 GC grid
				double zi = phi_zi_array[i-1];
				double wi = phi_wi_array[i-1];
				double omega_zi = phi_omega_zi_array[i-1];
				
				Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
				dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
			}
			
			for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
			{
				double q_minus = (K[k] - sqrt(K[k]*K[k] - phonon_energy_level_e))/sin(Phi[phi_q]);
				double q_pluss = (K[k] + sqrt(K[k]*K[k] - phonon_energy_level_e))/sin(Phi[phi_q]);
				
				// Q: Type 1 Gauss-Chebychev Grid (a,b)
				double q_a = q_minus;
				double q_b = q_pluss;
				for(int i = 1; i <= number_Qph_points; i++)
				{
					// Type 1 GC grid
					double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
					double wi = Pi/((double)number_Qph_points);
					double omega_zi = 1.0/sqrt(1.0-zi*zi);
					
					Q_ph[number_Qph_points-1-(i-1)]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
					dQ_ph[number_Qph_points-1-(i-1)] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
				}
				
				for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
				{
					double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
					double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
					
					// Second ee-phonon scattering
					double acos_arg = (Q_plane*Q_plane + 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
					double f_acos_arg = fabs(acos_arg);
					if (f_acos_arg < 1.0)
					{
						// Removed branch by multiplication of 2.0 below
						double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]

						double cos_k_q  = acos_arg;
						double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane - 2.0*K[k ] *Q_plane*cos_k_q);

						if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
						{
							//=============================
							//      e-P scattering part 2
							//=============================
							MAX_LENGTH_e2[k] += 1;
						} 
					} else { // end acos_arg < 1
						
						cout << "problems with limits 2e..." << endl;
						exit(-1);
					}
				}
			}
		}
		
		
		//=========== h dynamics =========================
		// Expression 1: (-Q_plane*Q_plane + 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane) < 1
		phi_minus = asin((-K[k] + sqrt(K[k]*K[k] + phonon_energy_level_h))/K_max);
		phi_pluss = Pi/2.0;
		
		asin_arg_pluss = ( K[k] + sqrt(K[k]*K[k] + phonon_energy_level_h))/K_max;
		if (asin_arg_pluss < 1.0)
		{
			phi_pluss = asin(asin_arg_pluss);
		}
		
		// GRID 1
		phi_a = phi_minus; // Using half the angular interval. Multiply solution by 2.0
		phi_b = phi_pluss;
		
		for(int i = 1; i <= number_Phi_points; i++)
		{
			// Type 1 GC grid
			double zi = phi_zi_array[i-1];
			double wi = phi_wi_array[i-1];
			double omega_zi = phi_omega_zi_array[i-1];
			
			Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
			dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
		}
		
		for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
		{
			double q_minus = (-K[k] + sqrt(K[k]*K[k] + phonon_energy_level_h))/fabs(sin(Phi[phi_q]));
			//double q_pluss = ( K[k] + sqrt(K[k]*K[k] + phonon_energy_level_h))/fabs(sin(Phi[phi_q]));
			
			// Q: Type 1 Gauss-Chebychev Grid (a,b)
			double q_a = q_minus;
			double q_b = K_max;
			for(int i = 1; i <= number_Qph_points; i++)
			{
				// Type 1 GC grid
				double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
				double wi = Pi/((double)number_Qph_points);
				double omega_zi = 1.0/sqrt(1.0-zi*zi);
				
				
				Q_ph[number_Qph_points-1-(i-1)]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
				dQ_ph[number_Qph_points-1-(i-1)] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
			}
			
			for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
			{
				double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
				double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
				
				double acos_arg = (-Q_plane*Q_plane + 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
				double f_acos_arg = fabs(acos_arg);
				if (f_acos_arg < 1.0)
				{
					// Removed branch by multiplication of 2.0 below
					double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]

					double cos_k_q  = acos_arg;
					double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane + 2.0*K[k ] *Q_plane*cos_k_q);

					if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
					{
						//==============================
						//      h-P scattering  part 1
						//==============================
						MAX_LENGTH_h1[k] += 1;
					}
				} else { // end acos_arg < 1
				
					for(int i = 0; i < number_Qph_points; i++)
					{
						cout << "q = " << Q_ph[q] << endl;
					}
					cout << "problems with limits 1h..." << endl;
					exit(-1);
				} 
			}
		}
		
		// GRID 2
		phi_a = phi_pluss; // Using half the angular interval. Multiply solution by 2.0
		phi_b = Pi/2.0;
		
		for(int i = 1; i <= number_Phi_points; i++)
		{
			// Type 1 GC grid
			double zi = phi_zi_array[i-1];
			double wi = phi_wi_array[i-1];
			double omega_zi = phi_omega_zi_array[i-1];
			
			Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
			dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
		}
		
		for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
		{
			double q_minus = (-K[k] + sqrt(K[k]*K[k] + phonon_energy_level_h))/fabs(sin(Phi[phi_q]));
			double q_pluss = ( K[k] + sqrt(K[k]*K[k] + phonon_energy_level_h))/fabs(sin(Phi[phi_q]));
			
			// Q: Type 1 Gauss-Chebychev Grid (a,b)
			double q_a = q_minus;
			double q_b = q_pluss;
			for(int i = 1; i <= number_Qph_points; i++)
			{
				// Type 1 GC grid
				double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
				double wi = Pi/((double)number_Qph_points);
				double omega_zi = 1.0/sqrt(1.0-zi*zi);
				
				
				Q_ph[number_Qph_points-i]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
				dQ_ph[number_Qph_points-i] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
			}
			
			
			for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
			{
				double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
				double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
				
				
				double acos_arg = (-Q_plane*Q_plane + 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
				double f_acos_arg = fabs(acos_arg);
				if (f_acos_arg < 1.0)
				{
					// Removed branch by multiplication of 2.0 below
					double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]

					double cos_k_q  = acos_arg;
					double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane + 2.0*K[k ] *Q_plane*cos_k_q);

					if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
					{
						//==============================
						//      h-P scattering  part 1
						//==============================
						MAX_LENGTH_h1[k] += 1;
					}
				} else { // end acos_arg < 1
				
					for(int i = 0; i < number_Qph_points; i++)
					{
						cout << "q = " << Q_ph[q] << endl;
					}
					cout << "problems with limits 1h..." << endl;
					exit(-1);
				} 
			}
		}
		
		
		// Expression 2: (Q_plane*Q_plane + 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane) < 1
		if (K[k] > sqrt(phonon_energy_level_h))
		{
			double phi_minus = asin((K[k] - sqrt(K[k]*K[k] - phonon_energy_level_h))/K_max);
			double phi_pluss = Pi/2.0;
			
			double asin_arg_pluss = (K[k] + sqrt(K[k]*K[k] - phonon_energy_level_h))/K_max;
			if (asin_arg_pluss < 1.0)
			{
				phi_pluss = asin(asin_arg_pluss);
			}
			
			// GRID 1
			double phi_a = phi_minus; // Using half the angular interval. Multiply solution by 2.0
			double phi_b = phi_pluss;
			
			for(int i = 1; i <= number_Phi_points; i++)
			{
				// Type 1 GC grid
				double zi = phi_zi_array[i-1];
				double wi = phi_wi_array[i-1];
				double omega_zi = phi_omega_zi_array[i-1];
				
				Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
				dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
			}
			
			for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
			{
				double q_minus = (K[k] - sqrt(K[k]*K[k] - phonon_energy_level_h))/sin(Phi[phi_q]);
				//double q_pluss = (K[k] + sqrt(K[k]*K[k] - phonon_energy_level_h))/sin(Phi[phi_q]);
				
				// Q: Type 1 Gauss-Chebychev Grid (a,b)
				double q_a = q_minus;
				double q_b = K_max;
				for(int i = 1; i <= number_Qph_points; i++)
				{
					// Type 1 GC grid
					double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
					double wi = Pi/((double)number_Qph_points);
					double omega_zi = 1.0/sqrt(1.0-zi*zi);
					
					Q_ph[number_Qph_points-1-(i-1)]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
					dQ_ph[number_Qph_points-1-(i-1)] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
				}
				
				for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
				{
					double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
					double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
					
					// Second hh-phonon scattering
					double acos_arg = (Q_plane*Q_plane + 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
					double f_acos_arg = fabs(acos_arg);
					if (f_acos_arg < 1.0)
					{
						// Removed branch by multiplication of 2.0 below
						double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]

						double cos_k_q  = acos_arg;
						double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane - 2.0*K[k ] *Q_plane*cos_k_q);

						if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
						{
							//==============================
							//      h-P scattering  part 2
							//==============================
							MAX_LENGTH_h2[k] += 1;
						}
					} else { // end acos_arg < 1
					
						for(int i = 0; i < number_Qph_points; i++)
						{
							cout << "q = " << Q_ph[q] << endl;
						}
						cout << "problems with limits 1h..." << endl;
						exit(-1);
					} 
				}
			}
			
			// GRID 2
			phi_a = phi_pluss; // Using half the angular interval. Multiply solution by 2.0
			phi_b = Pi/2.0;
			
			for(int i = 1; i <= number_Phi_points; i++)
			{
				// Type 1 GC grid
				double zi = phi_zi_array[i-1];
				double wi = phi_wi_array[i-1];
				double omega_zi = phi_omega_zi_array[i-1];
				
				Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
				dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
			}
			
			for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
			{
				double q_minus = (K[k] - sqrt(K[k]*K[k] - phonon_energy_level_h))/sin(Phi[phi_q]);
				double q_pluss = (K[k] + sqrt(K[k]*K[k] - phonon_energy_level_h))/sin(Phi[phi_q]);
				
				// Q: Type 1 Gauss-Chebychev Grid (a,b)
				double q_a = q_minus;
				double q_b = q_pluss;
				for(int i = 1; i <= number_Qph_points; i++)
				{
					// Type 1 GC grid
					double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
					double wi = Pi/((double)number_Qph_points);
					double omega_zi = 1.0/sqrt(1.0-zi*zi);
					
					Q_ph[number_Qph_points-1-(i-1)]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
					dQ_ph[number_Qph_points-1-(i-1)] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
				}
				
				for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
				{
					double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
					double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
					
					// Second hh-phonon scattering
					double acos_arg = (Q_plane*Q_plane + 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
					double f_acos_arg = fabs(acos_arg);
					if (f_acos_arg < 1.0)
					{
						// Removed branch by multiplication of 2.0 below
						double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]

						double cos_k_q  = acos_arg;
						double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane - 2.0*K[k ] *Q_plane*cos_k_q);

						if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
						{
							//==============================
							//      h-P scattering  part 2
							//==============================
							MAX_LENGTH_h2[k] += 1;
						}
					} else { // end acos_arg < 1
					
						for(int i = 0; i < number_Qph_points; i++)
						{
							cout << "q = " << Q_ph[q] << endl;
						}
						cout << "problems with limits 1h..." << endl;
						exit(-1);
					}
				}
			}
		}
	}

	// Prepare each index set
	if (carrier_scattering_phonon_index_e1a == NULL)
	{
		carrier_scattering_phonon_index_e1a = new indexSet2d_phonon(number_K_points, MAX_LENGTH_e1);
		carrier_scattering_phonon_index_e1b = new indexSet2d_phonon(number_K_points, MAX_LENGTH_e1);
		carrier_scattering_phonon_index_e2a = new indexSet2d_phonon(number_K_points, MAX_LENGTH_e2);
		carrier_scattering_phonon_index_e2b = new indexSet2d_phonon(number_K_points, MAX_LENGTH_e2);
		carrier_scattering_phonon_index_h1a = new indexSet2d_phonon(number_K_points, MAX_LENGTH_h1);
		carrier_scattering_phonon_index_h1b = new indexSet2d_phonon(number_K_points, MAX_LENGTH_h1);
		carrier_scattering_phonon_index_h2a = new indexSet2d_phonon(number_K_points, MAX_LENGTH_h2);
		carrier_scattering_phonon_index_h2b = new indexSet2d_phonon(number_K_points, MAX_LENGTH_h2);
	}


	for(int k = 0; k < number_K_points; k++)
	{
		int num_k_count_e1 = 0;
		int num_k_count_e2 = 0;
		int num_k_count_h1 = 0;
		int num_k_count_h2 = 0;
		
		//================== e dynamics ======================
		
		// Expression 1: (-Q_plane*Q_plane + 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane) < 1
		double phi_minus = asin((-K[k] + sqrt(K[k]*K[k] + phonon_energy_level_e))/K_max);
		double phi_pluss = Pi/2.0;
		
		double asin_arg_pluss = ( K[k] + sqrt(K[k]*K[k] + phonon_energy_level_e))/K_max;
		if (asin_arg_pluss < 1.0)
		{
			phi_pluss = asin(asin_arg_pluss);
		}
		
		// GRID 1
		double phi_a = phi_minus; // Using half the angular interval. Multiply solution by 2.0
		double phi_b = phi_pluss;
		
		for(int i = 1; i <= number_Phi_points; i++)
		{
			// Type 1 GC grid
			double zi = phi_zi_array[i-1];
			double wi = phi_wi_array[i-1];
			double omega_zi = phi_omega_zi_array[i-1];
			
			Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
			dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
		}
		
		for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
		{
			double q_minus = (-K[k] + sqrt(K[k]*K[k] + phonon_energy_level_e))/fabs(sin(Phi[phi_q]));
			//double q_pluss = ( K[k] + sqrt(K[k]*K[k] + phonon_energy_level_e))/fabs(sin(Phi[phi_q]));
			
			// Q: Type 1 Gauss-Chebychev Grid (a,b)
			double q_a = q_minus;
			double q_b = K_max;
			for(int i = 1; i <= number_Qph_points; i++)
			{
				// Type 1 GC grid
				double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
				double wi = Pi/((double)number_Qph_points);
				double omega_zi = 1.0/sqrt(1.0-zi*zi);
				
				
				Q_ph[number_Qph_points-1-(i-1)]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
				dQ_ph[number_Qph_points-1-(i-1)] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
			}
			
			for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
			{
				double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
				double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
				
				// First ee-phonon scattering
				double acos_arg = (-Q_plane*Q_plane + 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
				double f_acos_arg = fabs(acos_arg);
				if (f_acos_arg < 1.0)
				{
					// Removed branch by multiplication of 2.0 below
					double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]
					
					double cos_k_q  = acos_arg;
					double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane + 2.0*K[k ] *Q_plane*cos_k_q);

					if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
					{
						// Calculate g: Tilmann Kuhn "Density matrix theory of coherent ultrafast dynamics" 1998 p.173
						// Calculate g: Inès Waldmüller "Intersubband Dynamics in Semiconductor Quantum Wells: Linear and Nonlinear Response of Quantum Confined Electrons" 2005 p.52
						// convert from gauss units by 1.0/(4*pi*eps0)
						double form_factor_e	= misc_interpolate_K_array_parabolic(Q_perp  , phonon_form_factor_Q, phonon_form_factor_ee, number_K_points);
						double We_q		= form_factor_e/(Q_ph[q]*Q_ph[q] + Q_perp*Q_perp);
						double g_e2		= phonon_matrix_element_constants*We_q;
						
						//==========================
						// e-P scattering part. 1
						//==========================
						// Energy difference
						double DF_ee   =  1.0/(K[k]*Q_plane*fabs(sin(angle_th_k_q)));

						int index_k_q;
						double beta_k_q;

						// Occupation number calculation
						//double fe_k_q	= interpolate_cubic_spline_ne->evaluate(k_q);
						double fe_k_q	= misc_interpolate_K_array_linear_index(k_q  , K, ne_k_old, number_K_points,&index_k_q,&beta_k_q);

						// Integral volume constants
						double integral_volume_e = delta_scale_e*Q_ph[q]*Q_ph[q]*dQ_ph[q]*sin(Phi[phi_q])*dPhi[phi_q]*common_constants_phonon; // Constants from q*q*dq*sin(th1)*dt1

						double Fe_a      = integral_volume_e*4.0*g_e2*DF_ee*carrier_scattering_phonon_n;
						double Fe_b      = integral_volume_e*4.0*g_e2*DF_ee*(1.0+carrier_scattering_phonon_n);
						
						carrier_scattering_phonon_index_e1a->updateWeight(num_k_count_e1, Fe_a, k, q, index_k_q, beta_k_q);
						carrier_scattering_phonon_index_e1b->updateWeight(num_k_count_e1, Fe_b, k, q, index_k_q, beta_k_q);
						num_k_count_e1 += 1;
						
					} 
				} else { // end acos_arg < 1
				
					for(int i = 0; i < number_Qph_points; i++)
					{
						cout << "q = " << Q_ph[q] << endl;
					}
					cout << "problems with limits 1e..." << endl;
					exit(-1);
				} 
			}
		}
		
		// GRID 2
		phi_a = phi_pluss; // Using half the angular interval. Multiply solution by 2.0
		phi_b = Pi/2.0;
		
		for(int i = 1; i <= number_Phi_points; i++)
		{
			// Type 1 GC grid
			double zi = phi_zi_array[i-1];
			double wi = phi_wi_array[i-1];
			double omega_zi = phi_omega_zi_array[i-1];
			
			Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
			dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
		}
		
		for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
		{
			double q_minus = (-K[k] + sqrt(K[k]*K[k] + phonon_energy_level_e))/fabs(sin(Phi[phi_q]));
			double q_pluss = ( K[k] + sqrt(K[k]*K[k] + phonon_energy_level_e))/fabs(sin(Phi[phi_q]));
			
			// Q: Type 1 Gauss-Chebychev Grid (a,b)
			double q_a = q_minus;
			double q_b = q_pluss;
			for(int i = 1; i <= number_Qph_points; i++)
			{
				// Type 1 GC grid
				double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
				double wi = Pi/((double)number_Qph_points);
				double omega_zi = 1.0/sqrt(1.0-zi*zi);
				
				
				Q_ph[number_Qph_points-i]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
				dQ_ph[number_Qph_points-i] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
			}
			
			
			for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
			{
				double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
				double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
				
				// First ee-phonon scattering
				double acos_arg = (-Q_plane*Q_plane + 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
				double f_acos_arg = fabs(acos_arg);
				if (f_acos_arg < 1.0)
				{
					// Removed branch by multiplication of 2.0 below
					double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]
					
					double cos_k_q  = acos_arg;
					double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane + 2.0*K[k ] *Q_plane*cos_k_q);

					if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
					{
						// Calculate g: Tilmann Kuhn "Density matrix theory of coherent ultrafast dynamics" 1998 p.173
						// Calculate g: Inès Waldmüller "Intersubband Dynamics in Semiconductor Quantum Wells: Linear and Nonlinear Response of Quantum Confined Electrons" 2005 p.52
						// convert from gauss units by 1.0/(4*pi*eps0)
						double form_factor_e	= misc_interpolate_K_array_parabolic(Q_perp  , phonon_form_factor_Q, phonon_form_factor_ee, number_K_points);
						double We_q		= form_factor_e/(Q_ph[q]*Q_ph[q] + Q_perp*Q_perp);
						double g_e2		= phonon_matrix_element_constants*We_q;
						
						//==========================
						// e-P scattering part. 1
						//==========================
						// Energy difference
						double DF_ee   =  1.0/(K[k]*Q_plane*fabs(sin(angle_th_k_q)));

						int index_k_q;
						double beta_k_q;

						// Occupation number calculation
						//double fe_k_q	= interpolate_cubic_spline_ne->evaluate(k_q);
						double fe_k_q	= misc_interpolate_K_array_linear_index(k_q  , K, ne_k_old, number_K_points, &index_k_q, &beta_k_q);

						// Integral volume constants
						double integral_volume_e = delta_scale_e*Q_ph[q]*Q_ph[q]*dQ_ph[q]*sin(Phi[phi_q])*dPhi[phi_q]*common_constants_phonon; // Constants from q*q*dq*sin(th1)*dt1

						double Fe_a      = integral_volume_e*4.0*g_e2*DF_ee*carrier_scattering_phonon_n;
						double Fe_b      = integral_volume_e*4.0*g_e2*DF_ee*(1.0+carrier_scattering_phonon_n);
						
						carrier_scattering_phonon_index_e1a->updateWeight(num_k_count_e1, Fe_a, k, q, index_k_q, beta_k_q);
						carrier_scattering_phonon_index_e1b->updateWeight(num_k_count_e1, Fe_b, k, q, index_k_q, beta_k_q);
						num_k_count_e1 += 1;
					} 
				} else { // end acos_arg < 1
				
					for(int i = 0; i < number_Qph_points; i++)
					{
						cout << "q = " << Q_ph[q] << endl;
					}
					cout << "problems with limits 1e..." << endl;
					exit(-1);
				} 
			}
		}
		
		
		// Expression 2: (Q_plane*Q_plane + 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane) < 1
		if (K[k] > sqrt(phonon_energy_level_e))
		{
			double phi_minus = asin((K[k] - sqrt(K[k]*K[k] - phonon_energy_level_e))/K_max);
			double phi_pluss = Pi/2.0;
			
			double asin_arg_pluss = (K[k] + sqrt(K[k]*K[k] - phonon_energy_level_e))/K_max;
			if (asin_arg_pluss < 1.0)
			{
				phi_pluss = asin(asin_arg_pluss);
			}
			
			// GRID 1
			double phi_a = phi_minus; // Using half the angular interval. Multiply solution by 2.0
			double phi_b = phi_pluss;
			
			for(int i = 1; i <= number_Phi_points; i++)
			{
				// Type 1 GC grid
				double zi = phi_zi_array[i-1];
				double wi = phi_wi_array[i-1];
				double omega_zi = phi_omega_zi_array[i-1];
				
				Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
				dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
			}
			
			for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
			{
				double q_minus = (K[k] - sqrt(K[k]*K[k] - phonon_energy_level_e))/sin(Phi[phi_q]);
				//double q_pluss = (K[k] + sqrt(K[k]*K[k] - phonon_energy_level_e))/sin(Phi[phi_q]);
				
				// Q: Type 1 Gauss-Chebychev Grid (a,b)
				double q_a = q_minus;
				double q_b = K_max;
				for(int i = 1; i <= number_Qph_points; i++)
				{
					// Type 1 GC grid
					double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
					double wi = Pi/((double)number_Qph_points);
					double omega_zi = 1.0/sqrt(1.0-zi*zi);
					
					Q_ph[number_Qph_points-1-(i-1)]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
					dQ_ph[number_Qph_points-1-(i-1)] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
				}
				
				for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
				{
					double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
					double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
					
					// Second ee-phonon scattering
					double acos_arg = (Q_plane*Q_plane + 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
					double f_acos_arg = fabs(acos_arg);
					if (f_acos_arg < 1.0)
					{
						// Removed branch by multiplication of 2.0 below
						double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]
						
						double cos_k_q  = acos_arg;
						double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane - 2.0*K[k ] *Q_plane*cos_k_q);

						if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
						{
							// Calculate g: Tilmann Kuhn "Density matrix theory of coherent ultrafast dynamics" 1998 p.173
							// Calculate g: Inès Waldmüller "Intersubband Dynamics in Semiconductor Quantum Wells: Linear and Nonlinear Response of Quantum Confined Electrons" 2005 p.52
							// convert from gauss units by 1.0/(4*pi*eps0)
							double form_factor_e	= misc_interpolate_K_array_parabolic(Q_perp  , phonon_form_factor_Q, phonon_form_factor_ee, number_K_points);
							double We_q		= form_factor_e/(Q_ph[q]*Q_ph[q] + Q_perp*Q_perp);
							double g_e2		= phonon_matrix_element_constants*We_q;

							//==========================
							// e-P scattering part. 2
							//==========================
							// Energy difference
							double DF_ee   =  1.0/(K[k]*Q_plane*fabs(sin(angle_th_k_q)));
							
							int index_k_q;
							double beta_k_q;

							// Occupation number calculation
							//double fe_k_q	= interpolate_cubic_spline_ne->evaluate(k_q);
							double fe_k_q	= misc_interpolate_K_array_linear_index(k_q  , K, ne_k_old, number_K_points, &index_k_q, &beta_k_q);

							// Integral volume constants
							double integral_volume_e = delta_scale_e*Q_ph[q]*Q_ph[q]*dQ_ph[q]*sin(Phi[phi_q])*dPhi[phi_q]*common_constants_phonon; // Constants from q*q*dq*sin(th1)*dt1

							double Fe_a = integral_volume_e*4.0*g_e2*DF_ee*carrier_scattering_phonon_n;
							double Fe_b = integral_volume_e*4.0*g_e2*DF_ee*(1.0+carrier_scattering_phonon_n);
							
							carrier_scattering_phonon_index_e2a->updateWeight(num_k_count_e2, Fe_a, k, q, index_k_q, beta_k_q);
							carrier_scattering_phonon_index_e2b->updateWeight(num_k_count_e2, Fe_b, k, q, index_k_q, beta_k_q);
							num_k_count_e2 += 1;
						} 
					} else { // end acos_arg < 1
						
						cout << "problems with limits 2e..." << endl;
						exit(-1);
					}
				}
			}
			
			// GRID 2
			phi_a = phi_pluss; // Using half the angular interval. Multiply solution by 2.0
			phi_b = Pi/2.0;
			
			for(int i = 1; i <= number_Phi_points; i++)
			{
				// Type 1 GC grid
				double zi = phi_zi_array[i-1];
				double wi = phi_wi_array[i-1];
				double omega_zi = phi_omega_zi_array[i-1];
				
				Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
				dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
			}
			
			for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
			{
				double q_minus = (K[k] - sqrt(K[k]*K[k] - phonon_energy_level_e))/sin(Phi[phi_q]);
				double q_pluss = (K[k] + sqrt(K[k]*K[k] - phonon_energy_level_e))/sin(Phi[phi_q]);
				
				// Q: Type 1 Gauss-Chebychev Grid (a,b)
				double q_a = q_minus;
				double q_b = q_pluss;
				for(int i = 1; i <= number_Qph_points; i++)
				{
					// Type 1 GC grid
					double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
					double wi = Pi/((double)number_Qph_points);
					double omega_zi = 1.0/sqrt(1.0-zi*zi);
					
					Q_ph[number_Qph_points-1-(i-1)]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
					dQ_ph[number_Qph_points-1-(i-1)] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
				}
				
				for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
				{
					double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
					double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
					
					// Second ee-phonon scattering
					double acos_arg = (Q_plane*Q_plane + 2.0*device_me*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
					double f_acos_arg = fabs(acos_arg);
					if (f_acos_arg < 1.0)
					{
						// Removed branch by multiplication of 2.0 below
						double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]
						
						double cos_k_q  = acos_arg;
						double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane - 2.0*K[k ] *Q_plane*cos_k_q);

						if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
						{
							// Calculate g: Tilmann Kuhn "Density matrix theory of coherent ultrafast dynamics" 1998 p.173
							// Calculate g: Inès Waldmüller "Intersubband Dynamics in Semiconductor Quantum Wells: Linear and Nonlinear Response of Quantum Confined Electrons" 2005 p.52
							// convert from gauss units by 1.0/(4*pi*eps0)
							double form_factor_e	= misc_interpolate_K_array_parabolic(Q_perp  , phonon_form_factor_Q, phonon_form_factor_ee, number_K_points);
							double We_q		= form_factor_e/(Q_ph[q]*Q_ph[q] + Q_perp*Q_perp);
							double g_e2		= phonon_matrix_element_constants*We_q;

							//==========================
							// e-P scattering part. 2
							//==========================
							// Energy difference
							double DF_ee   =  1.0/(K[k]*Q_plane*fabs(sin(angle_th_k_q)));

							int index_k_q;
							double beta_k_q;

							// Occupation number calculation
							//double fe_k_q	= interpolate_cubic_spline_ne->evaluate(k_q);
							double fe_k_q	= misc_interpolate_K_array_linear_index(k_q  , K, ne_k_old, number_K_points, &index_k_q, &beta_k_q);

							// Integral volume constants
							double integral_volume_e = delta_scale_e*Q_ph[q]*Q_ph[q]*dQ_ph[q]*sin(Phi[phi_q])*dPhi[phi_q]*common_constants_phonon; // Constants from q*q*dq*sin(th1)*dt1

							double Fe_a = integral_volume_e*4.0*g_e2*DF_ee*carrier_scattering_phonon_n;
							double Fe_b = integral_volume_e*4.0*g_e2*DF_ee*(1.0+carrier_scattering_phonon_n);
							
							carrier_scattering_phonon_index_e2a->updateWeight(num_k_count_e2, Fe_a, k, q, index_k_q, beta_k_q);
							carrier_scattering_phonon_index_e2b->updateWeight(num_k_count_e2, Fe_b, k, q, index_k_q, beta_k_q);
							num_k_count_e2 += 1;
						} 
					} else { // end acos_arg < 1
						
						cout << "problems with limits 2e..." << endl;
						exit(-1);
					}
				}
			}
		}
		
		
		//=========== h dynamics =========================
		// Expression 1: (-Q_plane*Q_plane + 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane) < 1
		phi_minus = asin((-K[k] + sqrt(K[k]*K[k] + phonon_energy_level_h))/K_max);
		phi_pluss = Pi/2.0;
		
		asin_arg_pluss = ( K[k] + sqrt(K[k]*K[k] + phonon_energy_level_h))/K_max;
		if (asin_arg_pluss < 1.0)
		{
			phi_pluss = asin(asin_arg_pluss);
		}
		
		// GRID 1
		phi_a = phi_minus; // Using half the angular interval. Multiply solution by 2.0
		phi_b = phi_pluss;
		
		for(int i = 1; i <= number_Phi_points; i++)
		{
			// Type 1 GC grid
			double zi = phi_zi_array[i-1];
			double wi = phi_wi_array[i-1];
			double omega_zi = phi_omega_zi_array[i-1];
			
			Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
			dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
		}
		
		for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
		{
			double q_minus = (-K[k] + sqrt(K[k]*K[k] + phonon_energy_level_h))/fabs(sin(Phi[phi_q]));
			//double q_pluss = ( K[k] + sqrt(K[k]*K[k] + phonon_energy_level_h))/fabs(sin(Phi[phi_q]));
			
			// Q: Type 1 Gauss-Chebychev Grid (a,b)
			double q_a = q_minus;
			double q_b = K_max;
			for(int i = 1; i <= number_Qph_points; i++)
			{
				// Type 1 GC grid
				double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
				double wi = Pi/((double)number_Qph_points);
				double omega_zi = 1.0/sqrt(1.0-zi*zi);
				
				
				Q_ph[number_Qph_points-1-(i-1)]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
				dQ_ph[number_Qph_points-1-(i-1)] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
			}
			
			for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
			{
				double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
				double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
				
				double acos_arg = (-Q_plane*Q_plane + 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
				double f_acos_arg = fabs(acos_arg);
				if (f_acos_arg < 1.0)
				{
					// Removed branch by multiplication of 2.0 below
					double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]
					
					double cos_k_q  = acos_arg;
					double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane + 2.0*K[k ] *Q_plane*cos_k_q);

					if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
					{
						// Calculate g: Tilmann Kuhn "Density matrix theory of coherent ultrafast dynamics" 1998 p.173
						// Calculate g: Inès Waldmüller "Intersubband Dynamics in Semiconductor Quantum Wells: Linear and Nonlinear Response of Quantum Confined Electrons" 2005 p.52
						// convert from gauss units by 1.0/(4*pi*eps0)
						double form_factor_h	= misc_interpolate_K_array_parabolic(Q_perp  , phonon_form_factor_Q, phonon_form_factor_hh, number_K_points);
						double Wh_q		= form_factor_h/(Q_ph[q]*Q_ph[q] + Q_perp*Q_perp);
						double g_h2		= phonon_matrix_element_constants*Wh_q;
						
						//==========================
						// h-P scattering part. 1
						//==========================
						// Energy difference
						double DF_hh   =  1.0/(K[k]*Q_plane*fabs(sin(angle_th_k_q)));

						int index_k_q;
						double beta_k_q;

						// Occupation number calculation
						//double fh_k_q	= interpolate_cubic_spline_nh->evaluate(k_q);
						double fh_k_q	= misc_interpolate_K_array_linear_index(k_q  , K, nh_k_old, number_K_points, &index_k_q, &beta_k_q);

						// Integral volume constants
						double integral_volume_h = delta_scale_h*Q_ph[q]*Q_ph[q]*dQ_ph[q]*sin(Phi[phi_q])*dPhi[phi_q]*common_constants_phonon; // Constants from q*q*dq*sin(th1)*dt1

						double Fh_a      = integral_volume_h*4.0*g_h2*DF_hh*carrier_scattering_phonon_n;
						double Fh_b      = integral_volume_h*4.0*g_h2*DF_hh*(1.0+carrier_scattering_phonon_n);
						
						carrier_scattering_phonon_index_h1a->updateWeight(num_k_count_h1, Fh_a, k, q, index_k_q, beta_k_q);
						carrier_scattering_phonon_index_h1b->updateWeight(num_k_count_h1, Fh_b, k, q, index_k_q, beta_k_q);
						num_k_count_h1 += 1;
					}
				} else { // end acos_arg < 1
				
					for(int i = 0; i < number_Qph_points; i++)
					{
						cout << "q = " << Q_ph[q] << endl;
					}
					cout << "problems with limits 1h..." << endl;
					exit(-1);
				} 
			}
		}
		
		// GRID 2
		phi_a = phi_pluss; // Using half the angular interval. Multiply solution by 2.0
		phi_b = Pi/2.0;
		
		for(int i = 1; i <= number_Phi_points; i++)
		{
			// Type 1 GC grid
			double zi = phi_zi_array[i-1];
			double wi = phi_wi_array[i-1];
			double omega_zi = phi_omega_zi_array[i-1];
			
			Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
			dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
		}
		
		for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
		{
			double q_minus = (-K[k] + sqrt(K[k]*K[k] + phonon_energy_level_h))/fabs(sin(Phi[phi_q]));
			double q_pluss = ( K[k] + sqrt(K[k]*K[k] + phonon_energy_level_h))/fabs(sin(Phi[phi_q]));
			
			// Q: Type 1 Gauss-Chebychev Grid (a,b)
			double q_a = q_minus;
			double q_b = q_pluss;
			for(int i = 1; i <= number_Qph_points; i++)
			{
				// Type 1 GC grid
				double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
				double wi = Pi/((double)number_Qph_points);
				double omega_zi = 1.0/sqrt(1.0-zi*zi);
				
				
				Q_ph[number_Qph_points-i]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
				dQ_ph[number_Qph_points-i] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
			}
			
			
			for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
			{
				double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
				double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
				
				
				double acos_arg = (-Q_plane*Q_plane + 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
				double f_acos_arg = fabs(acos_arg);
				if (f_acos_arg < 1.0)
				{
					// Removed branch by multiplication of 2.0 below
					double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]
					
					double cos_k_q  = acos_arg;
					double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane + 2.0*K[k ] *Q_plane*cos_k_q);

					if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
					{
						// Calculate g: Tilmann Kuhn "Density matrix theory of coherent ultrafast dynamics" 1998 p.173
						// Calculate g: Inès Waldmüller "Intersubband Dynamics in Semiconductor Quantum Wells: Linear and Nonlinear Response of Quantum Confined Electrons" 2005 p.52
						// convert from gauss units by 1.0/(4*pi*eps0)
						double form_factor_h	= misc_interpolate_K_array_parabolic(Q_perp  , phonon_form_factor_Q, phonon_form_factor_hh, number_K_points);
						double Wh_q		= form_factor_h/(Q_ph[q]*Q_ph[q] + Q_perp*Q_perp);
						double g_h2		= phonon_matrix_element_constants*Wh_q;
						
						//==========================
						// h-P scattering part. 1
						//==========================
						// Energy difference
						double DF_hh   =  1.0/(K[k]*Q_plane*fabs(sin(angle_th_k_q)));

						int index_k_q;
						double beta_k_q;

						// Occupation number calculation
						//double fh_k_q	= interpolate_cubic_spline_nh->evaluate(k_q);
						double fh_k_q	= misc_interpolate_K_array_linear_index(k_q  , K, nh_k_old, number_K_points, &index_k_q, &beta_k_q);

						// Integral volume constants
						double integral_volume_h = delta_scale_h*Q_ph[q]*Q_ph[q]*dQ_ph[q]*sin(Phi[phi_q])*dPhi[phi_q]*common_constants_phonon; // Constants from q*q*dq*sin(th1)*dt1

						double Fh_a      = integral_volume_h*4.0*g_h2*DF_hh*carrier_scattering_phonon_n;
						double Fh_b      = integral_volume_h*4.0*g_h2*DF_hh*(1.0+carrier_scattering_phonon_n);
						
						carrier_scattering_phonon_index_h1a->updateWeight(num_k_count_h1, Fh_a, k, q, index_k_q, beta_k_q);
						carrier_scattering_phonon_index_h1b->updateWeight(num_k_count_h1, Fh_b, k, q, index_k_q, beta_k_q);
						num_k_count_h1 += 1;
					}
				} else { // end acos_arg < 1
				
					for(int i = 0; i < number_Qph_points; i++)
					{
						cout << "q = " << Q_ph[q] << endl;
					}
					cout << "problems with limits 1h..." << endl;
					exit(-1);
				} 
			}
		}
		
		
		// Expression 2: (Q_plane*Q_plane + 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane) < 1
		if (K[k] > sqrt(phonon_energy_level_h))
		{
			double phi_minus = asin((K[k] - sqrt(K[k]*K[k] - phonon_energy_level_h))/K_max);
			double phi_pluss = Pi/2.0;
			
			double asin_arg_pluss = (K[k] + sqrt(K[k]*K[k] - phonon_energy_level_h))/K_max;
			if (asin_arg_pluss < 1.0)
			{
				phi_pluss = asin(asin_arg_pluss);
			}
			
			// GRID 1
			double phi_a = phi_minus; // Using half the angular interval. Multiply solution by 2.0
			double phi_b = phi_pluss;
			
			for(int i = 1; i <= number_Phi_points; i++)
			{
				// Type 1 GC grid
				double zi = phi_zi_array[i-1];
				double wi = phi_wi_array[i-1];
				double omega_zi = phi_omega_zi_array[i-1];
				
				Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
				dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
			}
			
			for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
			{
				double q_minus = (K[k] - sqrt(K[k]*K[k] - phonon_energy_level_h))/sin(Phi[phi_q]);
				//double q_pluss = (K[k] + sqrt(K[k]*K[k] - phonon_energy_level_h))/sin(Phi[phi_q]);
				
				// Q: Type 1 Gauss-Chebychev Grid (a,b)
				double q_a = q_minus;
				double q_b = K_max;
				for(int i = 1; i <= number_Qph_points; i++)
				{
					// Type 1 GC grid
					double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
					double wi = Pi/((double)number_Qph_points);
					double omega_zi = 1.0/sqrt(1.0-zi*zi);
					
					Q_ph[number_Qph_points-1-(i-1)]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
					dQ_ph[number_Qph_points-1-(i-1)] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
				}
				
				for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
				{
					double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
					double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
					
					// Second hh-phonon scattering
					double acos_arg = (Q_plane*Q_plane + 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
					double f_acos_arg = fabs(acos_arg);
					if (f_acos_arg < 1.0)
					{
						// Removed branch by multiplication of 2.0 below
						double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]
						
						double cos_k_q  = acos_arg;
						double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane - 2.0*K[k ] *Q_plane*cos_k_q);

						if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
						{
							// Calculate g: Tilmann Kuhn "Density matrix theory of coherent ultrafast dynamics" 1998 p.173
							// Calculate g: Inès Waldmüller "Intersubband Dynamics in Semiconductor Quantum Wells: Linear and Nonlinear Response of Quantum Confined Electrons" 2005 p.52
							// convert from gauss units by 1.0/(4*pi*eps0)
							double form_factor_h	= misc_interpolate_K_array_parabolic(Q_perp  , phonon_form_factor_Q, phonon_form_factor_hh, number_K_points);
							double Wh_q		= form_factor_h/(Q_ph[q]*Q_ph[q] + Q_perp*Q_perp);
							double g_h2		= phonon_matrix_element_constants*Wh_q;
							
							//==========================
							// h-P scattering part. 2
							//==========================
							// Energy difference
							double DF_hh   =  1.0/(K[k]*Q_plane*fabs(sin(angle_th_k_q)));

							int index_k_q;
							double beta_k_q;

							// Occupation number calculation
							//double fh_k_q	= interpolate_cubic_spline_nh->evaluate(k_q);
							double fh_k_q	= misc_interpolate_K_array_linear_index(k_q  , K, nh_k_old, number_K_points, &index_k_q, &beta_k_q);

							// Integral volume constants
							double integral_volume_h = delta_scale_h*Q_ph[q]*Q_ph[q]*dQ_ph[q]*sin(Phi[phi_q])*dPhi[phi_q]*common_constants_phonon; // Constants from q*q*dq*sin(th1)*dt1

							double Fh_a = integral_volume_h*4.0*g_h2*DF_hh*carrier_scattering_phonon_n;
							double Fh_b = integral_volume_h*4.0*g_h2*DF_hh*(1.0+carrier_scattering_phonon_n);
							
							carrier_scattering_phonon_index_h2a->updateWeight(num_k_count_h2, Fh_a, k, q, index_k_q, beta_k_q);
							carrier_scattering_phonon_index_h2b->updateWeight(num_k_count_h2, Fh_b, k, q, index_k_q, beta_k_q);
							num_k_count_h2 += 1;
						}
					} else { // end acos_arg < 1
					
						for(int i = 0; i < number_Qph_points; i++)
						{
							cout << "q = " << Q_ph[q] << endl;
						}
						cout << "problems with limits 1h..." << endl;
						exit(-1);
					} 
				}
			}
			
			// GRID 2
			phi_a = phi_pluss; // Using half the angular interval. Multiply solution by 2.0
			phi_b = Pi/2.0;
			
			for(int i = 1; i <= number_Phi_points; i++)
			{
				// Type 1 GC grid
				double zi = phi_zi_array[i-1];
				double wi = phi_wi_array[i-1];
				double omega_zi = phi_omega_zi_array[i-1];
				
				Phi[number_Phi_points-i]  = phi_a +  0.5*(phi_b-phi_a)*(zi+1.0);
				dPhi[number_Phi_points-i] = 0.5*(phi_b-phi_a)*wi*(1.0/omega_zi);
			}
			
			for(int phi_q = 0; phi_q < number_Phi_points; phi_q++)
			{
				double q_minus = (K[k] - sqrt(K[k]*K[k] - phonon_energy_level_h))/sin(Phi[phi_q]);
				double q_pluss = (K[k] + sqrt(K[k]*K[k] - phonon_energy_level_h))/sin(Phi[phi_q]);
				
				// Q: Type 1 Gauss-Chebychev Grid (a,b)
				double q_a = q_minus;
				double q_b = q_pluss;
				for(int i = 1; i <= number_Qph_points; i++)
				{
					// Type 1 GC grid
					double zi = cos((2.0*((double)i)-1.0)*Pi/(2.0*number_Qph_points));
					double wi = Pi/((double)number_Qph_points);
					double omega_zi = 1.0/sqrt(1.0-zi*zi);
					
					Q_ph[number_Qph_points-1-(i-1)]  = q_a +  0.5*(q_b-q_a)*(zi+1.0);
					dQ_ph[number_Qph_points-1-(i-1)] = 0.5*(q_b-q_a)*wi*(1.0/omega_zi);
				}
				
				for(int q = 0; q < number_Qph_points; q++) // Momentum magnitude q
				{
					double Q_perp  = fabs(Q_ph[q]*cos(Phi[phi_q])); // z- component of momentum
					double Q_plane = fabs(Q_ph[q]*sin(Phi[phi_q])); // In plane momentum
					
					// Second hh-phonon scattering
					double acos_arg = (Q_plane*Q_plane + 2.0*device_mh*hbar_wLO*a0*a0/(hbar*hbar))/(2.0*K[k]*Q_plane);
					double f_acos_arg = fabs(acos_arg);
					if (f_acos_arg < 1.0)
					{
						// Removed branch by multiplication of 2.0 below
						double angle_th_k_q     = acos(acos_arg); // Angle is between [0,pi] or [-pi,0]

						// For speed
						double cos_k_q  = acos_arg;
						double k_q      = sqrt(K[k ]*K[k ]   + Q_plane*Q_plane - 2.0*K[k ] *Q_plane*cos_k_q);

						if ((k_q <= K[number_K_points - 1])&&(k_q >= K[0]))
						{
							// Calculate g: Tilmann Kuhn "Density matrix theory of coherent ultrafast dynamics" 1998 p.173
							// Calculate g: Inès Waldmüller "Intersubband Dynamics in Semiconductor Quantum Wells: Linear and Nonlinear Response of Quantum Confined Electrons" 2005 p.52
							// convert from gauss units by 1.0/(4*pi*eps0)
							double form_factor_h	= misc_interpolate_K_array_parabolic(Q_perp  , phonon_form_factor_Q, phonon_form_factor_hh, number_K_points);
							double Wh_q		= form_factor_h/(Q_ph[q]*Q_ph[q] + Q_perp*Q_perp);
							double g_h2		= phonon_matrix_element_constants*Wh_q;
							
							//==========================
							// h-P scattering part. 2
							//==========================
							// Energy difference
							double DF_hh   =  1.0/(K[k]*Q_plane*fabs(sin(angle_th_k_q)));
							
							int index_k_q;
							double beta_k_q;

							// Occupation number calculation
							//double fh_k_q	= interpolate_cubic_spline_nh->evaluate(k_q);
							double fh_k_q	= misc_interpolate_K_array_linear_index(k_q  , K, nh_k_old, number_K_points, &index_k_q, &beta_k_q);

							// Integral volume constants
							double integral_volume_h = delta_scale_h*Q_ph[q]*Q_ph[q]*dQ_ph[q]*sin(Phi[phi_q])*dPhi[phi_q]*common_constants_phonon; // Constants from q*q*dq*sin(th1)*dt1

							double Fh_a = integral_volume_h*4.0*g_h2*DF_hh*carrier_scattering_phonon_n;
							double Fh_b = integral_volume_h*4.0*g_h2*DF_hh*(1.0+carrier_scattering_phonon_n);
							
							carrier_scattering_phonon_index_h2a->updateWeight(num_k_count_h2, Fh_a, k, q, index_k_q, beta_k_q);
							carrier_scattering_phonon_index_h2b->updateWeight(num_k_count_h2, Fh_b, k, q, index_k_q, beta_k_q);
							num_k_count_h2 += 1;
						}
					} else { // end acos_arg < 1
					
						for(int i = 0; i < number_Qph_points; i++)
						{
							cout << "q = " << Q_ph[q] << endl;
						}
						cout << "problems with limits 1h..." << endl;
						exit(-1);
					}
				}
			}
		}
	}
	// Sort arrays
	carrier_scattering_phonon_index_e1a->sortIndices();
	carrier_scattering_phonon_index_e1b->sortIndices();
	carrier_scattering_phonon_index_e2a->sortIndices();
	carrier_scattering_phonon_index_e2b->sortIndices();
	carrier_scattering_phonon_index_h1a->sortIndices();
	carrier_scattering_phonon_index_h1b->sortIndices();
	carrier_scattering_phonon_index_h2a->sortIndices();
	carrier_scattering_phonon_index_h2b->sortIndices();


	delete [] MAX_LENGTH_e1;
	delete [] MAX_LENGTH_e2;
	delete [] MAX_LENGTH_h1;
	delete [] MAX_LENGTH_h2;

	// Clean up temporary arrays
	for(int k = 0; k < number_K_points; k++)
	{
		delete [] MAX_LENGTH_ee_k_kp[k];
		delete [] MAX_LENGTH_eh_k_kp[k];
		delete [] MAX_LENGTH_he_k_kp[k];
	}
	delete [] MAX_LENGTH_ee_k_kp;
	delete [] MAX_LENGTH_eh_k_kp;
	delete [] MAX_LENGTH_he_k_kp;
}

/* Calculate the FAST carrier scattering in PARALLEL using the
 * pre-computed index sets: carrier_scattering_index_##
 * Calculations are done in parallel for different k-points
 * ensure that the initialization function is run once before this
 * 
 * fk Terms have been reorganized as compared to the other formulas
 * */
void Device::carrier_scattering_method2_parallel(double *ne_k_old, double *nh_k_old)
{	
	double mass_hole_electron = getHoleMass()/getElectronMass(); // Rato of the mass of hole / electron

	//=============================
	// CARRIER-CARRIER scattering
	//=============================
	#ifdef USE_OPENMP
	#pragma omp parallel for schedule(static,OMP_THREADS_LEVEL_2) num_threads(OMP_THREADS_LEVEL_2)
	#endif
	for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
	{
		double tmp_out_ee = 0;
		double tmp_in_ee = 0;
		double tmp_out_hh = 0;
		double tmp_in_hh = 0;
		double tmp_out_eh = 0;
		double tmp_in_eh = 0;
		double tmp_out_he = 0;
		double tmp_in_he = 0;
		
		for(int kp = 0; kp < number_Kp_points; kp++) // Momentum magnitude kp
		{
			//=====================
			// e-e / h-h scattering
			//=====================
			int index_ee_total = carrier_scattering_index_ee->num_indices(k,kp);

			for(int index_ee = 0; index_ee < index_ee_total; index_ee++)
			{
				int i_kp 	= carrier_scattering_index_ee->I_kp(k,kp,index_ee);
				int i_kq 	= carrier_scattering_index_ee->I_k_q(k,kp,index_ee);
				int i_kpq 	= carrier_scattering_index_ee->I_kp_q(k,kp,index_ee);
				
				double beta_kp  = carrier_scattering_index_ee->beta_kp(k,kp,index_ee);
				double beta_kq  = carrier_scattering_index_ee->beta_k_q(k,kp,index_ee);
				double beta_kpq = carrier_scattering_index_ee->beta_kp_q(k,kp,index_ee);

				double fe_kp	= ne_k_old[i_kp]  + (ne_k_old[i_kp+1]  - ne_k_old[i_kp] )*beta_kp;
				double fe_k_q 	= ne_k_old[i_kq]  + (ne_k_old[i_kq+1]  - ne_k_old[i_kq] )*beta_kq;
				double fe_kp_q  = ne_k_old[i_kpq] + (ne_k_old[i_kpq+1] - ne_k_old[i_kpq])*beta_kpq;
				
				double fh_kp 	= nh_k_old[i_kp]  + (nh_k_old[i_kp+1]  - nh_k_old[i_kp] )*beta_kp;
				double fh_k_q 	= nh_k_old[i_kq]  + (nh_k_old[i_kq+1]  - nh_k_old[i_kq] )*beta_kq;
				double fh_kp_q  = nh_k_old[i_kpq] + (nh_k_old[i_kpq+1] - nh_k_old[i_kpq])*beta_kpq;

				double weight_e = carrier_scattering_index_ee->F(k,kp,index_ee);
				double weight_h = carrier_scattering_index_hh->F(k,kp,index_ee);
				
				tmp_out_ee += weight_e*(fe_kp_q*(1.0 - fe_kp - fe_k_q) + fe_kp*fe_k_q);
				tmp_in_ee  += weight_e*fe_kp*fe_k_q*(1.0-fe_kp_q);
				
				tmp_out_hh += weight_h*(fh_kp_q*(1.0 - fh_kp - fh_k_q) + fh_kp*fh_k_q);
				tmp_in_hh  += weight_h*fh_kp*fh_k_q*(1.0-fh_kp_q);
			}

		
			//===============
			// e-h scattering
			//===============
			int index_eh_total = carrier_scattering_index_eh->num_indices(k,kp);

			for(int index_eh = 0; index_eh < index_eh_total; index_eh++)
			{
				int i_kp 	= carrier_scattering_index_eh->I_kp(k,kp,index_eh);
				int i_kq 	= carrier_scattering_index_eh->I_k_q(k,kp,index_eh);
				int i_kpq 	= carrier_scattering_index_eh->I_kp_q(k,kp,index_eh);

				double fh_kp 	= nh_k_old[i_kp]  + (nh_k_old[i_kp+1]  - nh_k_old[i_kp] )*carrier_scattering_index_eh->beta_kp(k,kp,index_eh);
				double fe_k_q 	= ne_k_old[i_kq]  + (ne_k_old[i_kq+1]  - ne_k_old[i_kq] )*carrier_scattering_index_eh->beta_k_q(k,kp,index_eh);
				double fh_kp_q  = nh_k_old[i_kpq] + (nh_k_old[i_kpq+1] - nh_k_old[i_kpq])*carrier_scattering_index_eh->beta_kp_q(k,kp,index_eh);

				double weight = carrier_scattering_index_eh->F(k,kp,index_eh);

				tmp_out_eh += weight*(fh_kp*(1.0 - fe_k_q - fh_kp_q) + fe_k_q*fh_kp_q);
				tmp_in_eh  += weight*fe_k_q*fh_kp_q*(1.0 - fh_kp);
			}

			
			//===============
			// h-e scattering
			//===============
			int index_he_total = carrier_scattering_index_he->num_indices(k,kp);

			for(int index_he = 0; index_he < index_he_total; index_he++)
			{
				int i_kp 	= carrier_scattering_index_he->I_kp(k,kp,index_he);
				int i_kq 	= carrier_scattering_index_he->I_k_q(k,kp,index_he);
				int i_kpq 	= carrier_scattering_index_he->I_kp_q(k,kp,index_he);

				double fe_kp 	= ne_k_old[i_kp]  + (ne_k_old[i_kp+1]  - ne_k_old[i_kp] )*carrier_scattering_index_he->beta_kp(k,kp,index_he);
				double fh_k_q 	= nh_k_old[i_kq]  + (nh_k_old[i_kq+1]  - nh_k_old[i_kq] )*carrier_scattering_index_he->beta_k_q(k,kp,index_he);
				double fe_kp_q  = ne_k_old[i_kpq] + (ne_k_old[i_kpq+1] - ne_k_old[i_kpq])*carrier_scattering_index_he->beta_kp_q(k,kp,index_he);

				double weight = carrier_scattering_index_he->F(k,kp,index_he);

				tmp_out_he += weight*(fe_kp*(1.0 - fh_k_q - fe_kp_q) + fh_k_q*fe_kp_q);
				tmp_in_he  += weight*fh_k_q*fe_kp_q*(1.0 - fe_kp);
			}
		}
		carrier_scattering_rates_e_out[k] = tmp_out_ee + tmp_out_eh;
		carrier_scattering_rates_e_in[k]  = tmp_in_ee  + tmp_in_eh;

		carrier_scattering_rates_h_out[k] = tmp_out_hh + tmp_out_he;
		carrier_scattering_rates_h_in[k]  = tmp_in_hh  + tmp_in_he;

		//========================================================================================================
		// PHONON-CARRIER scattering

		double tmp_e1a = 0;
		double tmp_e1b = 0;
		double tmp_e2a = 0;
		double tmp_e2b = 0;
		double tmp_h1a = 0;
		double tmp_h1b = 0;
		double tmp_h2a = 0;
		double tmp_h2b = 0;
		
		//===========
		// First e-e

		int index_ee_total = carrier_scattering_phonon_index_e1a->get_num_cols(k);

		for(int index_ee = 0; index_ee < index_ee_total; index_ee++)
		{
			int i_kq 	= carrier_scattering_phonon_index_e1a->I_kq(k,index_ee);
			double beta_kq  = carrier_scattering_phonon_index_e1a->beta_kq(k,index_ee);
			double fe_k_q 	= ne_k_old[i_kq]  + (ne_k_old[i_kq+1]  - ne_k_old[i_kq] )*beta_kq;

			double weight_ee1a = carrier_scattering_phonon_index_e1a->F(k,index_ee);
			double weight_ee1b = carrier_scattering_phonon_index_e1b->F(k,index_ee);
			
			tmp_e1a += weight_ee1a*(1.0-fe_k_q);
			tmp_e1b += weight_ee1b*(fe_k_q);
		}

		// Second e-e

		index_ee_total = carrier_scattering_phonon_index_e2a->get_num_cols(k);
		for(int index_ee = 0; index_ee < index_ee_total; index_ee++)
		{
			int i_kq 	= carrier_scattering_phonon_index_e2a->I_kq(k,index_ee);
			double beta_kq  = carrier_scattering_phonon_index_e2a->beta_kq(k,index_ee);
			double fe_k_q 	= ne_k_old[i_kq]  + (ne_k_old[i_kq+1]  - ne_k_old[i_kq] )*beta_kq;

			double weight_ee2a = carrier_scattering_phonon_index_e2a->F(k,index_ee);
			double weight_ee2b = carrier_scattering_phonon_index_e2b->F(k,index_ee);
			
			tmp_e2a += weight_ee2a*(fe_k_q);
			tmp_e2b += weight_ee2b*(1.0-fe_k_q);
		}
		//===========
		// First h-h

		int index_hh_total = carrier_scattering_phonon_index_h1a->get_num_cols(k);

		for(int index_hh = 0; index_hh < index_hh_total; index_hh++)
		{
			int i_kq 	= carrier_scattering_phonon_index_h1a->I_kq(k,index_hh);
			double beta_kq  = carrier_scattering_phonon_index_h1a->beta_kq(k,index_hh);
			double fh_k_q 	= nh_k_old[i_kq]  + (nh_k_old[i_kq+1]  - nh_k_old[i_kq] )*beta_kq;

			double weight_hh1a = carrier_scattering_phonon_index_h1a->F(k,index_hh);
			double weight_hh1b = carrier_scattering_phonon_index_h1b->F(k,index_hh);
			
			tmp_h1a += weight_hh1a*(1.0-fh_k_q);
			tmp_h1b += weight_hh1b*(fh_k_q);
		}

		// Second h-h

		index_hh_total = carrier_scattering_phonon_index_h2a->get_num_cols(k);
		for(int index_hh = 0; index_hh < index_hh_total; index_hh++)
		{
			int i_kq 	= carrier_scattering_phonon_index_h2a->I_kq(k,index_hh);
			double beta_kq  = carrier_scattering_phonon_index_h2a->beta_kq(k,index_hh);
			double fh_k_q 	= nh_k_old[i_kq]  + (nh_k_old[i_kq+1]  - nh_k_old[i_kq] )*beta_kq;

			double weight_hh2a = carrier_scattering_phonon_index_h2a->F(k,index_hh);
			double weight_hh2b = carrier_scattering_phonon_index_h2b->F(k,index_hh);
			
			tmp_h2a += weight_hh2a*(fh_k_q);
			tmp_h2b += weight_hh2b*(1.0-fh_k_q);
		}

		carrier_scattering_rates_phonon_e_in[k]  = tmp_e1b + tmp_e2a;
		carrier_scattering_rates_phonon_e_out[k] = tmp_e1a + tmp_e1b + tmp_e2a + tmp_e2b;
		carrier_scattering_rates_phonon_h_in[k]  = tmp_h1b + tmp_h2a;
		carrier_scattering_rates_phonon_h_out[k] = tmp_h1a + tmp_h1b + tmp_h2a + tmp_h2b;
	}
	

	for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
	{
		carrier_scattering_rates_e_total[k] = carrier_scattering_rates_e_in[k] - ne_k_old[k]*carrier_scattering_rates_e_out[k];
		carrier_scattering_rates_h_total[k] = carrier_scattering_rates_h_in[k] - nh_k_old[k]*carrier_scattering_rates_h_out[k];

		carrier_scattering_rates_phonon_e_total[k] = carrier_scattering_rates_phonon_e_in[k] - ne_k_old[k]*carrier_scattering_rates_phonon_e_out[k];
		carrier_scattering_rates_phonon_h_total[k] = carrier_scattering_rates_phonon_h_in[k] - nh_k_old[k]*carrier_scattering_rates_phonon_h_out[k];
	}


	//========================================
	// Conservation of particle numbers
	//========================================
	{
		double error_ne = 0.0;
		double error_ne_p = 0.0;
		double error_ne_m = 0.0;
		for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
		{
			double tmp = carrier_scattering_rates_e_total[k]*KdK[k];
			error_ne += tmp;
			if (carrier_scattering_rates_e_total[k] < 0)
			{
				error_ne_m += tmp;
			} else {
				error_ne_p += tmp;
			}
		}
		
		if (error_ne < 0.0)
		{
			if (error_ne_m < 0.0)
			{
				double tmp = 1.0-error_ne/error_ne_m;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_e_total[k] < 0.0)
					{
						carrier_scattering_rates_e_total[k] *= tmp;
					}
				}
			}
		} else {
			if (error_ne_p > 0.0)
			{
				double tmp = 1.0-error_ne/error_ne_p;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_e_total[k] > 0)
					{
						carrier_scattering_rates_e_total[k] *= tmp;
					}
				}
			}
		}
		
		double error_nh = 0.0;
		double error_nh_p = 0.0;
		double error_nh_m = 0.0;
		for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
		{
			double tmp = carrier_scattering_rates_h_total[k]*KdK[k];
			error_nh += tmp;
			if (carrier_scattering_rates_h_total[k] < 0)
			{
				error_nh_m += tmp;
			} else {
				error_nh_p += tmp;
			}
		}
		
		if (error_nh < 0.0)
		{
			if (error_nh_m < 0.0)
			{
				double tmp = 1.0-error_nh/error_nh_m;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_h_total[k] < 0.0)
					{
						carrier_scattering_rates_h_total[k] *= tmp;
					}
				}
			}
		} else {
			if (error_nh_p > 0.0)
			{
				double tmp = 1.0-error_nh/error_nh_p;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_h_total[k] > 0)
					{
						carrier_scattering_rates_h_total[k] *= tmp;
					}
				}
			}
		}
	}

	// Phonons
	{
		double error_ne = 0.0;
		double error_ne_p = 0.0;
		double error_ne_m = 0.0;
		for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
		{
			double tmp = carrier_scattering_rates_phonon_e_total[k]*KdK[k];
			error_ne += tmp;
			if (carrier_scattering_rates_phonon_e_total[k] < 0)
			{
				error_ne_m += tmp;
			} else {
				error_ne_p += tmp;
			}
		}
		
		if (error_ne < 0.0)
		{
			if (error_ne_m < 0.0)
			{
				double tmp = 1.0-error_ne/error_ne_m;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_phonon_e_total[k] < 0.0)
					{
						carrier_scattering_rates_phonon_e_total[k] *= tmp;
					}
				}
			}
		} else {
			if (error_ne_p > 0.0)
			{
				double tmp = 1.0-error_ne/error_ne_p;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_phonon_e_total[k] > 0)
					{
						carrier_scattering_rates_phonon_e_total[k] *= tmp;
					}
				}
			}
		}
		
		double error_nh = 0.0;
		double error_nh_p = 0.0;
		double error_nh_m = 0.0;
		for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
		{
			double tmp = carrier_scattering_rates_phonon_h_total[k]*KdK[k];
			error_nh += tmp;
			if (carrier_scattering_rates_phonon_h_total[k] < 0)
			{
				error_nh_m += tmp;
			} else {
				error_nh_p += tmp;
			}
		}
		
		if (error_nh < 0.0)
		{
			if (error_nh_m < 0.0)
			{
				double tmp = 1.0-error_nh/error_nh_m;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_phonon_h_total[k] < 0.0)
					{
						carrier_scattering_rates_phonon_h_total[k] *= tmp;
					}
				}
			}
		} else {
			if (error_nh_p > 0.0)
			{
				double tmp = 1.0-error_nh/error_nh_p;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_phonon_h_total[k] > 0)
					{
						carrier_scattering_rates_phonon_h_total[k] *= tmp;
					}
				}
			}
		}
	}

	// carrier-carrier scattering: Conservation of Energy
	{
		double sumE_e = 0.0;
		double sumE_h = 0.0;
		for(int k = 0; k < number_K_points; k++) 
		{
			sumE_e += Ek_e[k]*carrier_scattering_rates_e_total[k]*KdK[k];
			sumE_h += Ek_h[k]*carrier_scattering_rates_h_total[k]*KdK[k];
		}

		if ((fabs(sumE_e) >0.0)&&(fabs(sumE_h) > 0.0))
		{
			// Minimum energy
			if (fabs(sumE_e) < fabs(sumE_h))
			{
				double const1 = fabs(sumE_e/sumE_h);
				for(int k = 0; k < number_K_points; k++) 
				{
					carrier_scattering_rates_h_total[k] *= const1;
				}
			} else {
				double const1 = fabs(sumE_h/sumE_e);
				for(int k = 0; k < number_K_points; k++) 
				{
					carrier_scattering_rates_e_total[k] *= const1;
				}

			}
			
		} else if ((fabs(sumE_e) == 0.0)&&(fabs(sumE_h) > 0.0))
		{
			for(int k = 0; k < number_K_points; k++) 
			{
				carrier_scattering_rates_h_total[k] = 0.0; // set min value
			}
		} else if ((fabs(sumE_e) > 0.0)&&(fabs(sumE_h) == 0.0))
		{
			for(int k = 0; k < number_K_points; k++) 
			{
				carrier_scattering_rates_e_total[k] = 0.0; // set min value
			}
		}
	}
}

/* Assume rates are constant and only change ne_k, nh_k
 * */
void Device::carrier_scattering_method2_parallel_constant_rate(double *ne_k_old, double *nh_k_old)
{	
	double mass_hole_electron = getHoleMass()/getElectronMass(); // Rato of the mass of hole / electron
	//=============================
	// CARRIER-CARRIER scattering
	//=============================

	for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
	{
		// Carrier scattering
		carrier_scattering_rates_e_total[k] = carrier_scattering_rates_e_in[k] - ne_k_old[k]*carrier_scattering_rates_e_out[k];
		carrier_scattering_rates_h_total[k] = carrier_scattering_rates_h_in[k] - nh_k_old[k]*carrier_scattering_rates_h_out[k];

		// Phonon scattering
		carrier_scattering_rates_phonon_e_total[k] = carrier_scattering_rates_phonon_e_in[k] - ne_k_old[k]*carrier_scattering_rates_phonon_e_out[k];
		carrier_scattering_rates_phonon_h_total[k] = carrier_scattering_rates_phonon_h_in[k] - nh_k_old[k]*carrier_scattering_rates_phonon_h_out[k];
	}
	
	
	//========================================
	// Conservation of particle numbers
	
	// Jorg particle number balance
	{
		double error_ne = 0.0;
		double error_ne_p = 0.0;
		double error_ne_m = 0.0;
		for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
		{
			double tmp = carrier_scattering_rates_e_total[k]*KdK[k];
			error_ne += tmp;
			if (carrier_scattering_rates_e_total[k] < 0)
			{
				error_ne_m += tmp;
			} else {
				error_ne_p += tmp;
			}
		}
		
		if (error_ne < 0.0)
		{
			if (error_ne_m < 0.0)
			{
				double tmp = 1.0-error_ne/error_ne_m;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_e_total[k] < 0.0)
					{
						carrier_scattering_rates_e_total[k] *= tmp;
					}
				}
			}
		} else {
			if (error_ne_p > 0.0)
			{
				double tmp = 1.0-error_ne/error_ne_p;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_e_total[k] > 0)
					{
						carrier_scattering_rates_e_total[k] *= tmp;
					}
				}
			}
		}
		
		double error_nh = 0.0;
		double error_nh_p = 0.0;
		double error_nh_m = 0.0;
		for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
		{
			double tmp = carrier_scattering_rates_h_total[k]*KdK[k];
			error_nh += tmp;
			if (carrier_scattering_rates_h_total[k] < 0)
			{
				error_nh_m += tmp;
			} else {
				error_nh_p += tmp;
			}
		}
		
		if (error_nh < 0.0)
		{
			if (error_nh_m < 0.0)
			{
				double tmp = 1.0-error_nh/error_nh_m;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_h_total[k] < 0.0)
					{
						carrier_scattering_rates_h_total[k] *= tmp;
					}
				}
			}
		} else {
			if (error_nh_p > 0.0)
			{
				double tmp = 1.0-error_nh/error_nh_p;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_h_total[k] > 0)
					{
						carrier_scattering_rates_h_total[k] *= tmp;
					}
				}
			}
		}
	}

	//========================================
	// Conservation of particle numbers
	//========================================
	// Jorg particle number balance
	{
		double error_ne = 0.0;
		double error_ne_p = 0.0;
		double error_ne_m = 0.0;
		for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
		{
			double tmp = carrier_scattering_rates_phonon_e_total[k]*KdK[k];
			error_ne += tmp;
			if (carrier_scattering_rates_phonon_e_total[k] < 0)
			{
				error_ne_m += tmp;
			} else {
				error_ne_p += tmp;
			}
		}
		
		if (error_ne < 0.0)
		{
			if (error_ne_m < 0.0)
			{
				double tmp = 1.0-error_ne/error_ne_m;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_phonon_e_total[k] < 0.0)
					{
						carrier_scattering_rates_phonon_e_total[k] *= tmp;
					}
				}
			}
		} else {
			if (error_ne_p > 0.0)
			{
				double tmp = 1.0-error_ne/error_ne_p;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_phonon_e_total[k] > 0)
					{
						carrier_scattering_rates_phonon_e_total[k] *= tmp;
					}
				}
			}
		}
		
		double error_nh = 0.0;
		double error_nh_p = 0.0;
		double error_nh_m = 0.0;
		for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
		{
			double tmp = carrier_scattering_rates_phonon_h_total[k]*KdK[k];
			error_nh += tmp;
			if (carrier_scattering_rates_phonon_h_total[k] < 0)
			{
				error_nh_m += tmp;
			} else {
				error_nh_p += tmp;
			}
		}
		
		if (error_nh < 0.0)
		{
			if (error_nh_m < 0.0)
			{
				double tmp = 1.0-error_nh/error_nh_m;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_phonon_h_total[k] < 0.0)
					{
						carrier_scattering_rates_phonon_h_total[k] *= tmp;
					}
				}
			}
		} else {
			if (error_nh_p > 0.0)
			{
				double tmp = 1.0-error_nh/error_nh_p;
				for(int k = 0; k < number_K_points; k++) // Momentum magnitude k
				{
					if (carrier_scattering_rates_phonon_h_total[k] > 0)
					{
						carrier_scattering_rates_phonon_h_total[k] *= tmp;
					}
				}
			}
		}
	}
	// carrier-carrier scattering: Conservation of Energy
	{
		double sumE_e = 0.0;
		double sumE_h = 0.0;
		for(int k = 0; k < number_K_points; k++) 
		{
			sumE_e += Ek_e[k]*carrier_scattering_rates_e_total[k]*KdK[k];
			sumE_h += Ek_h[k]*carrier_scattering_rates_h_total[k]*KdK[k];
		}

		if ((fabs(sumE_e) >0.0)&&(fabs(sumE_h) > 0.0))
		{
			// Minimum energy
			if (fabs(sumE_e) < fabs(sumE_h))
			{
				double const1 = fabs(sumE_e/sumE_h);
				for(int k = 0; k < number_K_points; k++) 
				{
					carrier_scattering_rates_h_total[k] *= const1;
				}
			} else {
				double const1 = fabs(sumE_h/sumE_e);
				for(int k = 0; k < number_K_points; k++) 
				{
					carrier_scattering_rates_e_total[k] *= const1;
				}

			}
			
		} else if ((fabs(sumE_e) == 0.0)&&(fabs(sumE_h) > 0.0))
		{
			for(int k = 0; k < number_K_points; k++) 
			{
				carrier_scattering_rates_h_total[k] = 0.0; // set min value
			}
		} else if ((fabs(sumE_e) > 0.0)&&(fabs(sumE_h) == 0.0))
		{
			for(int k = 0; k < number_K_points; k++) 
			{
				carrier_scattering_rates_e_total[k] = 0.0; // set min value
			}
		}
	}
}


/* Set the background Fermi distribution from the instantanious variables
 * Helper function to update instantaneous fermi functions
 * */
void Device::carrier_scattering_rate_set_background_fermi()
{
	for(unsigned i = 0; i < number_K_points; i++)
	{
		fe_k_inst[i] = misc_get_fermi_distribution(device_inst_temp_e, device_inst_chem_pot_e, device_me, K[i]);
		fh_k_inst[i] = misc_get_fermi_distribution(device_inst_temp_h, device_inst_chem_pot_h, device_mh, K[i]);
	}
}

/* Set the instantanious temperature and chem. potential
 * device_inst_temp_e, device_inst_temp_h, device_inst_chem_pot_e, device_inst_chem_pot_h
 * 
 * Uses the current carriers ne_k nh_k
 * 
 * Find the best fit for a Fermi function to the current carriers.
 * The result is the temperature and the chem. pot.
 **/
void Device::carrier_scattering_rate_set_instant_temperature(double *ne_k_old, double *nh_k_old)
{
	//======================================
	// Find correct Temperature of carriers
	//======================================
	// Total number of carriers
	double Nsum_e = 0, Nsum_h = 0;
	for(unsigned i = 0; i < number_K_points; i++)
	{
		Nsum_e += ne_k_old[i]*KdK[i];
		Nsum_h += nh_k_old[i]*KdK[i];
	}
	Nsum_e /= (Pi*a0*a0);
	Nsum_h /= (Pi*a0*a0);
	//Nsum_h = Nsum_e;


//	double Nsum_e = misc_simpsons_quadrature_kdk(ne_k);
//	double Nsum_h = misc_simpsons_quadrature_kdk(nh_k);
//	double Nsum_h = Nsum_e;

	// Total energy in carriers
	double Esum_e = 0, Esum_h = 0;
	for(unsigned i = 0; i < number_K_points; i++)
	{
		Esum_e += Ek_e[i]*ne_k_old[i]*KdK[i];
		Esum_h += Ek_h[i]*nh_k_old[i]*KdK[i];
	}
	Esum_e /= (Pi*a0*a0);
	Esum_h /= (Pi*a0*a0);

//	double Esum_e = misc_simpsons_quadrature_Ek_kdk(Ek_e, ne_k);
//	double Esum_h = misc_simpsons_quadrature_Ek_kdk(Ek_h, nh_k);
	
	// Use as initial guess, the previous values
	double T0_e  = device_inst_temp_e;
	double T0_h  = device_inst_temp_h;
	double mu0_e = device_inst_chem_pot_e;
	double mu0_h = device_inst_chem_pot_h;

	//mu0_h = misc_get_chem_potential_analytic_2d(Nsum_h, device_mh, T0_h);
	//mu0_e = misc_get_chem_potential_analytic_2d(Nsum_e, device_me, T0_e);
	//carrier_scattering_rate_get_instant_temp_and_chemPot(&T0_h, &mu0_h, Nsum_h, Esum_h, device_mh, device_background_temp_h);
	//carrier_scattering_rate_get_instant_temp_and_chemPot(&T0_e, &mu0_e, Nsum_e, Esum_e, device_me, device_background_temp_e);
	carrier_scattering_rate_get_instant_temp_and_chemPot2(&T0_h, &mu0_h, Nsum_h, Esum_h, device_mh, device_background_temp_h);
	carrier_scattering_rate_get_instant_temp_and_chemPot2(&T0_e, &mu0_e, Nsum_e, Esum_e, device_me, device_background_temp_e);

	// Store for next iteration
	device_inst_temp_e = T0_e;
	device_inst_temp_h = T0_h;
	device_inst_chem_pot_h = mu0_h;
	device_inst_chem_pot_e = mu0_e;
}

/* Set the pair (T,mu) that fits with the total (N,E)
 * T -> Returned temperature
 * mu -> Returned chem. potential
 * N -> Total carrier density
 * E -> Total energy
 * m -> Effective mass of particle
 * t_bgr -> Background temperature of distribution
 * */
void Device::carrier_scattering_rate_get_instant_temp_and_chemPot(double *T, double *mu, double N, double E, double m, double T_bgr)
{
	// Assuming we have a 2D distribution we can use the analytic formaula for the Chem. Pot.

	double beta  = (kB/hbar)*sqrt(m/(E*Pi));
	double alpha = (hbar*N)*sqrt(Pi/(E*m));
	double y = beta*(*T);

	// Case that should not happen, but let it relax to backgrounid density
/*
	if ((E < 0) ||(N <0))
	{
		*T = device_background_temp;
		*mu = misc_get_chem_potential_analytic_2d(device_density, m, *T);
		return;
	}
*/

//	cout << "alpha = " << alpha << endl;
//	cout << "beta  = " << beta << endl;
//	cout << "NEWTON: T0 = " << y/beta << endl;

	// Solution of equations only exist if alpha<=sqrt(2), otherwise we get negative temperature
	// Keep old values in this case
	if (alpha < sqrt(2))
	{


		int MAXITER = 100;
		bool found_solution = false;
		double TOL = 0.5e-6;
		for(unsigned i = 0; i < MAXITER; i++)
		{
			double f,df;
			carrier_scattering_rate_get_instant_temp_and_chemPot_argMin(y, alpha, &f,&df);

			double newY = y - f/df;
			double y1;
			if (newY >= 0)
			{
				y1 = newY;
				
			} else {
				// Do one step of bisection between y=0 and current to try and find a zero	
				y1 = 0.5*(0.0 + y);
			}
			

			if (abs(y1-y)/abs(y1) < TOL)
			{
				found_solution = true;
				break;
			}

			y = y1;
		
		}

		// Find zero using newton
		if (found_solution)
		{
			// Force the estimated temperature to be Greater than the lattice temperature
			if (y/beta >= T_bgr)
			{
				*T  = y/beta;
				*mu = misc_get_chem_potential_analytic_2d(N, m, *T);
			}  else {
				//cout << "carrier_scattering_rate_get_instant_temp_and_chemPot: Newton found solution, but temp (" << y/beta << ") is less than lattice temperature -> Keeping t constant: " << *T << endl;
			}
			
		} else {
			//cout << "carrier_scattering_rate_get_instant_temp_and_chemPot: Newton did not converge.. keeping T constant" << *T << endl;
		}

	}
	else {
	//	cout << "NEWTON: T = " << *T << endl;
	}
}

/* Helper function that is used in the calculatino of the instantaneous temperature and density
 * Return minimizer and its derivative
 * */
void Device::carrier_scattering_rate_get_instant_temp_and_chemPot_argMin(double y, double alpha, double *f, double *df)
{
	// Evaluate the argument function -y^2 PolyLog(2,-1/A(alpha,y)) - beta

	double A = 1.0/(exp(alpha/y)-1.0);
	double MAX = 1.0/A;
	if (MAX < 0)
	{
		cout << "carrier_scattering_rate_get_instant_temp_and_chemPot_argMin(): Cannot use argument MAX < 0" << endl;
		cout << "alpha = " << alpha << endl;
		cout << "y     = " << y << endl;
		exit(-1);
	}

	// In the limit where y ->0 from above the argument function is
	if (y==0.0)
	{
		*f  = 0.5*(-2.0+alpha*alpha);
		*df = 0.0;
		return;
	}

	// Simpsons rule
	// f[x] = PolyLog[2, -MAX] = Int[-log(1+x)/x,{x,0,MAX}]
	double ASYMPTOTE_LIMIT = 1.0e9;

	if (MAX < ASYMPTOTE_LIMIT)
	{
		double sum = misc_polyLog2(MAX);

		*f  = -y*y*sum - 1.0;
		*df = alpha*alpha/(y*(exp(-alpha/y)-1.0)) - 2.0*y*sum;
		return;
	} else {

	
/*
		cout << "Cannot evaluate misc_polyLog2 here" << endl;
		cout << "alpha = " << alpha << endl;
		cout << "y     = " << y << endl;
		cout << "MAX   = " << MAX << endl;
		exit(-1);
*/		
		
		double sum = misc_polyLog2(ASYMPTOTE_LIMIT);
		double lm1 = log(MAX);
		double lm2 = log(ASYMPTOTE_LIMIT);
		double newSum = sum-0.5*(lm1*lm1-lm2*lm2);

		*f  = -y*y*newSum -1.0;
		*df = alpha*lm1/(exp(-alpha/y)-1.0) + y*lm1*lm1 - y*lm2*lm2 - 2.0*y*sum;
		return;
	} 
	
}

void Device::carrier_scattering_rate_get_instant_temp_and_chemPot2(double *T, double *mu, double N, double E, double m, double T_bgr)
{
	// Assuming we have a 2D distribution we can use the analytic formaula for the Chem. Pot.
	double beta = kB*m/(hbar*hbar*Pi*N);
	double alpha = E*beta/(kB*N);
	double y = beta*(*T);

/*
	// Case that should not happen, but let it relax to backgrounid density
	cout << "beta  = " << beta << endl;
	cout << "alpha = " << alpha << endl;
	cout << "y     = " << y << endl;
	cout << "NEWTON: T0 = " << y/beta << endl;
*/

	// Smallest possible energy that fits a given density and a POSITIVE temperature.
	// If the energy gets below this level, no Fermi function exists
 	// Assume constant temperature in this case! (Equal to last step)
	//double E_min = 0.5*kB*N/beta;
	double alpha_min = 0.5;

	// Solve system of equations
	if (alpha > alpha_min)
	{


		int MAXITER = 100;
		bool found_solution = false;
		double TOL = 0.5e-6;
		for(unsigned i = 0; i < MAXITER; i++)
		{
			double f,df;
			carrier_scattering_rate_get_instant_temp_and_chemPot2_argMin(y, alpha, beta, &f,&df);

			double newY = y - f/df;
			//cout << "y = " << y << ", f = " << f << ", df = " << df << ", newY = " << newY << endl;

			double y1;
			if (newY > 0)
			{
				y1 = newY;
				
			} else {
				// The alg. was equal to or missed y=0
				// Do one step of bisection between y=0 and current to try and find a zero	
				y1 = 0.5*(0.0 + y);
			}
			

			if (abs(y1-y)/abs(y1) < TOL)
			{
				found_solution = true;
				break;
			}

			y = y1;
		
		}

		// Did we find a solution?
		if (found_solution)
		{
			// Force the estimated temperature to be Greater than the lattice temperature
			if (y/beta >= T_bgr)
			{
				*T  = y/beta;
				*mu = misc_get_chem_potential_analytic_2d(N, m, *T);
			}  else {
				//cout << "carrier_scattering_rate_get_instant_temp_and_chemPot: Newton found solution, but temp (" << y/beta << ") is less than lattice temperature -> Keeping t constant: " << *T << endl;
			}
			
		} else {
			//cout << "carrier_scattering_rate_get_instant_temp_and_chemPot: Newton did not converge.. keeping T constant" << *T << endl;
		}

	}
	else {
	//	cout << "NEWTON: T = " << *T << endl;
	}
}

/* Helper function that is used in the calculation of the instantaneous temperature and density
 * Return minimizer and its derivative
 * */
void Device::carrier_scattering_rate_get_instant_temp_and_chemPot2_argMin(double y, double alpha, double beta, double *f, double *df)
{
	// In the limit where y ->0+ the function is
	if (y==0.0)
	{
		*f  = 0.5;
		*df = 0.0;
		return;
	}

	// Evaluate the argument function alpha -y^2 PolyLog(2,-exp(eps0))
	double y_inv = 1.0/y;
	double exp_inv_y = exp(y_inv);
	double A = exp_inv_y-1.0;

	// Simpsons rule
	// f[x] = PolyLog[2, -A] = Int[-log(1+x)/x,{x,0,A}]
	double sum = misc_polyLog2(A);
	*f  = alpha+y*y*sum;
	*df = 2.0*y*sum + (y_inv)*exp_inv_y/A;
	return;
	
}

//=======================================
// File IO functions and their helpers

/* Remove whitespaces from a string
 * Used when importing data from a file
 * */
void Device::file_removeWhitespace(std::string& str) 
{
    for (size_t i = 0; i < str.length(); i++) {
        if (str[i] == ' ' || str[i] == '\n' || str[i] == '\t') {
            str.erase(i, 1);
            i--;
        }
    }
}

/* Extract a double from a string
 * Removes anything after the pattern "//" inclusive
 * */
double Device::file_getDoubleFromLine(std::string line)
{
	string line2 = line;
	std::string mystr = line.substr(0, line.find("//", 0));
	//trim(mystr);		// Using boost
	file_removeWhitespace(mystr);
	
	if (mystr.empty())
	{
		cout << "Import failed, no value found in line" << endl;
		cout << "Line: |" << line2 << "|" << endl;
		exit(-1);
	}
	
	return atof(mystr.c_str());
}

/* Extract an int from a string
 * Removes anything after the pattern "//" inclusive
 * */
int Device::file_getIntFromLine(std::string line)
{
	string line2 = line;
	std::string mystr = line.substr(0, line.find("//", 0));
	//trim(mystr);		// Using boost
	file_removeWhitespace(mystr);
	
	if (mystr.empty())
	{
		cout << "Import failed, no value found in line" << endl;
		cout << "Line: |" << line2 << "|" << endl;
		exit(-1);
	}
	
	return atoi(mystr.c_str());
}

/* Read in QW config file from disk
 * device_name should be XX_T# where # is the transverse point
 * The device name must match the filename
 * material_XX.config where XX is the device name
 * */
void Device::file_readInConfigFromFile()
{
	std::ifstream device_config;
	std::stringstream fileName;
	fileName << "material_" << getName().substr(0, getName().find("_T",0)) << ".config";
	// Check if file exists
	if (!fileExists(fileName.str()))
	{
		cout << "file_readInConfigFromFile(): Cannot find file: " << fileName.str() << endl;
	}
	device_config.open((fileName.str()).c_str(), std::ofstream::in); 
	string line,line2;
	if (device_config)
	{
		getline(device_config,line); // Get first line
		number_K_points = file_getIntFromLine(line);
		
		getline(device_config,line);
		K_max = file_getDoubleFromLine(line);
		
		getline(device_config,line);
		device_me = file_getDoubleFromLine(line)*m0;
		
		getline(device_config,line);
		device_mh = file_getDoubleFromLine(line)*m0;
		
		getline(device_config,line);
		device_dipolemoment = file_getDoubleFromLine(line)*e;
		
		getline(device_config,line);
		device_background_temp_e = file_getDoubleFromLine(line);

		double backgroundTemp = 300.0;
		double maxTemp = device_background_temp_e;
		double deltaT = device_transverse_temp_scale*(maxTemp-backgroundTemp);

		// Modify according to pump scaling
		if (device_name.substr(0,2) == "QW")
		{
			device_background_temp_e = backgroundTemp + deltaT;
		}
		
		getline(device_config,line);
		device_background_temp_h = file_getDoubleFromLine(line);

		// Modify according to pump scaling
		if (device_name.substr(0,2) == "QW")
		{
			device_background_temp_h = backgroundTemp + deltaT;
		}
		
		getline(device_config,line);
		device_length_qw = file_getDoubleFromLine(line)*nm;
		
		getline(device_config,line);
		device_focus_E = file_getDoubleFromLine(line);
		
		getline(device_config,line);		
		device_effective_qw = file_getDoubleFromLine(line);

		getline(device_config,line);		
		device_deph_scale = file_getDoubleFromLine(line);

		getline(device_config,line);
		device_bandgap = file_getDoubleFromLine(line)*e;

		// Modify according to temperature
		if (device_name.substr(0,2) == "QW")
		{
			double lambda1 = 2.0*Pi*c0/(device_bandgap/hbar);
			double lambda2 = lambda1 + (4.0*nm/10.0)*deltaT; // 4nm/10K taken from Alex 2.Nov 2018
	
			device_bandgap = hbar*2.0*Pi*c0/lambda2;
		}
		
		getline(device_config,line);
		device_density = file_getDoubleFromLine(line);

		// Modify according to pump scaling
		if (device_name.substr(0,2) == "QW")
		{
			device_density = 5.0e14 + device_transverse_pump_scale*(device_density - 5.0e14);
		}
		
		getline(device_config,line);
		device_occ_relax_time = file_getDoubleFromLine(line)*ps;

		getline(device_config,line);
		device_occ_hole_relax_time = file_getDoubleFromLine(line)*ps;

		

		
		cout << "Load from file: " << fileName.str() << endl;
		
		cout << "number_K_points       = " << number_K_points << endl;
		cout << "K_max                 = " << K_max << " [1/a0] ( " << K_max/a0 << " [1/m])" << endl;
		cout << "me                    = " << device_me/m0 << " [m0]" << endl;
		cout << "mh                    = " << device_mh/m0 << " [m0]" << endl;
		cout << "device_dcv            = " << device_dipolemoment/e << " [eV]" << endl;
		cout << "T_background_e        = " << device_background_temp_e << " [K]" << endl;
		cout << "T_background_h        = " << device_background_temp_h << " [K]" << endl;
		cout << "device_length_qw      = " << device_length_qw/nm << " [nm]" << endl;
		cout << "device_focus_E        = " << device_focus_E << endl;
		cout << "device_eff_qw         = " << device_effective_qw << endl;
		cout << "device_deph_scale     = " << device_deph_scale << endl;
		cout << "device_bandgap        = " << device_bandgap/e << " [eV]" << endl;
		cout << "device_density        = " << device_density << " [1/m^2]" << endl;
		cout << "device_occ_relax      = " << device_occ_relax_time/ps << " [ps]" << endl;
		cout << "device_occ_hole_relax = " << device_occ_hole_relax_time/ps << " [ps]" << endl;
		
	} else {
	
		cout << "Device::file_readInConfigFromFile(): could not open file: " << fileName.str() << endl; 
		exit(-1);
	}
	
	device_config.close();
}

/* Write output variables to files
 * File name: OUTPUTKEY__OUTPUTNUMBER_VARIABLE.dat
 * OUTPUTKEY is set in the creation of the Device
 * OUTPUTNUMBER is set with "out_count"
 * VARIABLE is the spesific variable
 * 
 * The output is appended to each datafile
 * The output level is for control over what variables are written to file
 * The default is ALL variables
 * */
void Device::file_output_write(int output_level)
{
	double tmp = real(getElectricField());
	output_E_real.write(reinterpret_cast<const char*>(&tmp),sizeof(double));
	tmp = imag(getElectricField());
	output_E_imag.write(reinterpret_cast<const char*>(&tmp),sizeof(double));

/*
	tmp = ((double)device_effective_qw)*real(macroscopic_polarization);
	output_macpol_re.write(reinterpret_cast<const char*>(&tmp),sizeof(double));
	tmp = ((double)device_effective_qw)*imag(macroscopic_polarization);
	output_macpol_im.write(reinterpret_cast<const char*>(&tmp),sizeof(double));
*/
	if (output_level == 1)
	{
		double tmp_pol[number_K_points];
		for(int j=0; j<number_K_points; j++)
		{
			tmp_pol[j]=real(p_k[j]);
		}
		output_p_k_re.write(reinterpret_cast<const char*>(tmp_pol),number_K_points*sizeof(double));
		for(int j=0; j<number_K_points; j++)
		{
			tmp_pol[j]=imag(p_k[j]);
		}
		output_p_k_im.write(reinterpret_cast<const char*>(tmp_pol),number_K_points*sizeof(double));
		output_ne_k.write(reinterpret_cast<const char*>(ne_k),number_K_points*sizeof(double));
		output_nh_k.write(reinterpret_cast<const char*>(nh_k),number_K_points*sizeof(double));
		
	//	output_fe_k.write(reinterpret_cast<const char*>(fe_k_inst),number_K_points*sizeof(double));
	//	output_fh_k.write(reinterpret_cast<const char*>(fh_k_inst),number_K_points*sizeof(double));
		output_energy_renormalized_e_k.write(reinterpret_cast<const char*>(energy_renormalized_e),number_K_points*sizeof(double));
		output_energy_renormalized_h_k.write(reinterpret_cast<const char*>(energy_renormalized_h),number_K_points*sizeof(double));
	/*	
		output_scattering_out_e.write(reinterpret_cast<const char*>(carrier_scattering_rates_e_out),number_K_points*sizeof(double));
		output_scattering_in_e.write(reinterpret_cast<const char*>(carrier_scattering_rates_e_in),number_K_points*sizeof(double));
		output_scattering_out_h.write(reinterpret_cast<const char*>(carrier_scattering_rates_h_out),number_K_points*sizeof(double));
		output_scattering_in_h.write(reinterpret_cast<const char*>(carrier_scattering_rates_h_in),number_K_points*sizeof(double));
 	*/
		
		double tmp4_e[number_K_points];
		double tmp4_h[number_K_points];

	#if defined(USE_ISAK_HOLE_FILLING) || defined(USE_ISAK_HOLE_FILLING_TABLE)
		for(int i = 0; i < number_K_points; i++)
		{
			tmp4_e[i] = carrier_scattering_rates_e_total[i];
			tmp4_h[i] = carrier_scattering_rates_h_total[i];
		}
 		output_scattering_total_e.write(reinterpret_cast<const char*>(tmp4_e),number_K_points*sizeof(double));
		output_scattering_total_h.write(reinterpret_cast<const char*>(tmp4_h),number_K_points*sizeof(double));

		for(int i = 0; i < number_K_points; i++)
		{
			tmp4_e[i] = carrier_scattering_rates_phonon_e_total[i];
			tmp4_h[i] = carrier_scattering_rates_phonon_h_total[i];
		}
 		output_scattering_total_phonon_e.write(reinterpret_cast<const char*>(tmp4_e),number_K_points*sizeof(double));
		output_scattering_total_phonon_h.write(reinterpret_cast<const char*>(tmp4_h),number_K_points*sizeof(double));
	#else
 		output_scattering_total_e.write(reinterpret_cast<const char*>(carrier_scattering_rate_approximation_e),number_K_points*sizeof(double));
		output_scattering_total_h.write(reinterpret_cast<const char*>(carrier_scattering_rate_approximation_h),number_K_points*sizeof(double));
	#endif



		
		for(int i = 0; i < number_K_points; i++)
		{
			tmp4_e[i] = -2.0*imag((getElectricField() + electric_field_renormalized[i])*conj(p_k[i]));
		}
 		output_scattering_Rabi.write(reinterpret_cast<const char*>(tmp4_e),number_K_points*sizeof(double));

		for(int i = 0; i < number_K_points; i++)
		{
			tmp4_e[i] = -(ne_k[i]-fe_k[i])*SBE_OCC_PUMP;
			tmp4_h[i] = -(nh_k[i]-fh_k[i])*SBE_OCC_PUMP;
		}
 		output_scattering_Pump_e.write(reinterpret_cast<const char*>(tmp4_e),number_K_points*sizeof(double));
 		output_scattering_Pump_h.write(reinterpret_cast<const char*>(tmp4_h),number_K_points*sizeof(double));

		// SPONT
	#ifdef SBE_USE_SPONTAN_EMISSIONS
		for(int i = 0; i < number_K_points; i++)
		{
			tmp4_e[i] = device_spont_emission_wk[i]*ne_k[i]*nh_k[i];
		}
 		output_scattering_Spont_emission.write(reinterpret_cast<const char*>(tmp4_e),number_K_points*sizeof(double));
	#endif
	}

	double Nsum_e = 0, Nsum_h = 0, Esum_e = 0, Esum_h = 0,inversion_positive = 0;
	for(unsigned i = 0; i < number_K_points; i++)
	{
		Nsum_e += ne_k[i]*KdK[i];
		Nsum_h += nh_k[i]*KdK[i];
		Esum_e += ne_k[i]*KdK[i]*Ek_e[i];
		Esum_h += nh_k[i]*KdK[i]*Ek_h[i];
	}
	Nsum_e /= (Pi*a0*a0);
	Nsum_h /= (Pi*a0*a0);
	Esum_e /= (Pi*a0*a0);
	Esum_h /= (Pi*a0*a0);

	//inversion_positive /= (Pi*a0*a0);
//	double Nsum_e = misc_simpsons_quadrature_kdk(ne_k);
//	double Nsum_h = misc_simpsons_quadrature_kdk(nh_k);
//	double Esum_e = misc_simpsons_quadrature_Ek_kdk(Ek_e, ne_k);
//	double Esum_h = misc_simpsons_quadrature_Ek_kdk(Ek_h, nh_k);

	output_Nsum_e.write(reinterpret_cast<const char*>(&Nsum_e),sizeof(double));
	output_Nsum_h.write(reinterpret_cast<const char*>(&Nsum_h),sizeof(double));
/*
	output_Esum_e.write(reinterpret_cast<const char*>(&Esum_e),sizeof(double));
	output_Esum_h.write(reinterpret_cast<const char*>(&Esum_h),sizeof(double));
	//output_inversion_positive.write(reinterpret_cast<const char*>(&inversion_positive),sizeof(double));
	double Fsum_e = 0, Fsum_h = 0;
	for(unsigned i = 0; i < number_K_points; i++)
	{
		Fsum_e += fe_k_inst[i]*KdK[i];
		Fsum_h += fh_k_inst[i]*KdK[i];
	}
	Fsum_e /= (Pi*a0*a0);
	Fsum_h /= (Pi*a0*a0);

//	double Fsum_e = misc_simpsons_quadrature_kdk(fe_k_inst);
//	double Fsum_h = misc_simpsons_quadrature_kdk(fh_k_inst);

	output_Fsum_e.write(reinterpret_cast<const char*>(&Fsum_e),sizeof(double));
	output_Fsum_h.write(reinterpret_cast<const char*>(&Fsum_h),sizeof(double));
*/
	// Update inst_temp
	carrier_scattering_rate_set_instant_temperature(ne_k, nh_k);

	output_inst_temp_e.write(reinterpret_cast<const char*>(&device_inst_temp_e),sizeof(double));
	output_inst_temp_h.write(reinterpret_cast<const char*>(&device_inst_temp_h),sizeof(double));
//	output_inst_chem_pot_e.write(reinterpret_cast<const char*>(&device_inst_chem_pot_e),sizeof(double));
//	output_inst_chem_pot_h.write(reinterpret_cast<const char*>(&device_inst_chem_pot_h),sizeof(double));
}

/* Open output files for writing
 * File name: OUTPUTKEY__OUTPUTNUMBER_VARIABLE.dat
 * OUTPUTKEY is set in the creation of the Device
 * OUTPUTNUMBER is set with "out_count"
 * VARIABLE is the spesific variable
 * 
 * Creates and prepares the files for output as needed
 * Opens to append data if needed
 * */
void Device::file_output_open(int out_count, int output_level)
{


	// Setup output to files
	std::stringstream baseName;
	baseName << getToFileOutputKey() << out_count;
	
	std::stringstream fileName;
	fileName << baseName.str() << "_E_re_" << getName() << ".dat";
	openAppendBinary(&output_E_real, fileName.str());
	
	fileName.str(std::string());
	fileName << baseName.str() << "_E_im_" << getName() << ".dat";
	openAppendBinary(&output_E_imag, fileName.str());

/*
	fileName.str(std::string());
	fileName << baseName.str() << "_macpol_re_" << getName() << ".dat";
	openAppendBinary(&output_macpol_re, fileName.str());

	fileName.str(std::string());
	fileName << baseName.str() << "_macpol_im_" << getName() << ".dat";
	openAppendBinary(&output_macpol_im, fileName.str());
*/
	
	if (output_level == 1)
	{
		fileName.str(std::string());
		fileName << baseName.str() << "_p_re_" << getName() << ".dat";
		openAppendBinary(&output_p_k_re, fileName.str());
		
		fileName.str(std::string());
		fileName << baseName.str() << "_p_im_" << getName() << ".dat";
		openAppendBinary(&output_p_k_im, fileName.str());
	
		fileName.str(std::string());
		fileName << baseName.str() << "_ne_" << getName() << ".dat";
		openAppendBinary(&output_ne_k, fileName.str());
		
		fileName.str(std::string());
		fileName << baseName.str() << "_nh_" << getName() << ".dat";
		openAppendBinary(&output_nh_k, fileName.str());

		fileName.str(std::string());
		fileName << baseName.str() << "_fe_" << getName() << ".dat";
		openAppendBinary(&output_fe_k, fileName.str());
		
		fileName.str(std::string());
		fileName << baseName.str() << "_fh_" << getName() << ".dat";
		openAppendBinary(&output_fh_k, fileName.str());

		fileName.str(std::string());
		fileName << baseName.str() << "_energy_renormalized_e_" << getName() << ".dat";
		openAppendBinary(&output_energy_renormalized_e_k, fileName.str());
		
		fileName.str(std::string());
		fileName << baseName.str() << "_energy_renormalized_h_" << getName() << ".dat";
		openAppendBinary(&output_energy_renormalized_h_k, fileName.str());

		fileName.str(std::string());
		fileName << baseName.str() << "_scattering_rate_total_e_" << getName() << ".dat";
		openAppendBinary(&output_scattering_total_e, fileName.str());

		fileName.str(std::string());
		fileName << baseName.str() << "_scattering_rate_total_h_" << getName() << ".dat";
		openAppendBinary(&output_scattering_total_h, fileName.str());

		fileName.str(std::string());
		fileName << baseName.str() << "_scattering_rate_total_phonon_e_" << getName() << ".dat";
		openAppendBinary(&output_scattering_total_phonon_e, fileName.str());

		fileName.str(std::string());
		fileName << baseName.str() << "_scattering_rate_total_phonon_h_" << getName() << ".dat";
		openAppendBinary(&output_scattering_total_phonon_h, fileName.str());

		fileName.str(std::string());
		fileName << baseName.str() << "_Rabi_" << getName() << ".dat";
		openAppendBinary(&output_scattering_Rabi, fileName.str());

		fileName.str(std::string());
		fileName << baseName.str() << "_Pump_e_" << getName() << ".dat";
		openAppendBinary(&output_scattering_Pump_e, fileName.str());

		fileName.str(std::string());
		fileName << baseName.str() << "_Pump_h_" << getName() << ".dat";
		openAppendBinary(&output_scattering_Pump_h, fileName.str());

		fileName.str(std::string());
		fileName << baseName.str() << "_Spont_emission_" << getName() << ".dat";
		openAppendBinary(&output_scattering_Spont_emission, fileName.str());
	}
	fileName.str(std::string());
	fileName << baseName.str() << "_Nsum_e_" << getName() << ".dat";
	openAppendBinary(&output_Nsum_e, fileName.str());
	
	fileName.str(std::string());
	fileName << baseName.str() << "_Nsum_h_" << getName() << ".dat";
	openAppendBinary(&output_Nsum_h, fileName.str());

/*	
	fileName.str(std::string());
	fileName << baseName.str() << "_Esum_e_" << getName() << ".dat";
	openAppendBinary(&output_Esum_e, fileName.str());
	
	fileName.str(std::string());
	fileName << baseName.str() << "_Esum_h_" << getName() << ".dat";
	openAppendBinary(&output_Esum_h, fileName.str());
	
	fileName.str(std::string());
	fileName << baseName.str() << "_Fsum_e_" << getName() << ".dat";
	openAppendBinary(&output_Fsum_e, fileName.str());
	
	fileName.str(std::string());
	fileName << baseName.str() << "_Fsum_h_" << getName() << ".dat";
	openAppendBinary(&output_Fsum_h, fileName.str());
*/
	
	fileName.str(std::string());
	fileName << baseName.str() << "_inst_temp_e_" << getName() << ".dat";
	openAppendBinary(&output_inst_temp_e, fileName.str());
	
	fileName.str(std::string());
	fileName << baseName.str() << "_inst_temp_h_" << getName() << ".dat";
	openAppendBinary(&output_inst_temp_h, fileName.str());
/*		
	fileName.str(std::string());
	fileName << baseName.str() << "_inst_chem_pot_e_" << getName() << ".dat";
	openAppendBinary(&output_inst_chem_pot_e, fileName.str());
	
	fileName.str(std::string());
	fileName << baseName.str() << "_inst_chem_pot_h_" << getName() << ".dat";
	openAppendBinary(&output_inst_chem_pot_h, fileName.str());
*/

/*
	fileName.str(std::string());
	fileName << baseName.str() << "_inversion_positive" << getName() << ".dat";
	openAppendBinary(&output_inversion_positive, fileName.str());
*/
	
	/*
	fileName.str(std::string());
	fileName << baseName.str() << "_scattering_rate_out_e_" << getName() << ".dat";
	openAppendBinary(&output_scattering_out_e, fileName.str());

	fileName.str(std::string());
	fileName << baseName.str() << "_scattering_rate_in_e_" << getName() << ".dat";
	openAppendBinary(&output_scattering_in_e, fileName.str());

	fileName.str(std::string());
	fileName << baseName.str() << "_scattering_rate_out_h_" << getName() << ".dat";
	openAppendBinary(&output_scattering_out_h, fileName.str());

	fileName.str(std::string());
	fileName << baseName.str() << "_scattering_rate_in_h_" << getName() << ".dat";
	openAppendBinary(&output_scattering_in_h, fileName.str());
	*/

}

/* Close output file
 * File name: OUTPUTKEY__OUTPUTNUMBER_VARIABLE.dat
 * OUTPUTKEY is set in the creation of the Device
 * OUTPUTNUMBER is set with "out_count"
 * VARIABLE is the spesific variable
 * 
 * Close the output steam for all files
 * */
void Device::file_output_close(int output_level)
{
	output_E_real.close();
	output_E_imag.close();
/*
	output_macpol_re.close();
	output_macpol_im.close();
*/
	if (output_level == 1)
	{
		output_p_k_re.close();
		output_p_k_im.close();
		output_ne_k.close();
		output_nh_k.close();
		output_fe_k.close();
		output_fh_k.close();
		output_energy_renormalized_e_k.close();
		output_energy_renormalized_h_k.close();
		output_scattering_total_e.close();
		output_scattering_total_h.close();
		output_scattering_total_phonon_e.close();
		output_scattering_total_phonon_h.close();

		output_scattering_Rabi.close();
		output_scattering_Pump_e.close();
		output_scattering_Pump_h.close();
		output_scattering_Spont_emission.close();
	}
	output_Nsum_e.close();
	output_Nsum_h.close();
/*
	output_Esum_e.close();
	output_Esum_h.close();
	output_Fsum_e.close();
	output_Fsum_h.close();
*/
	output_inst_temp_e.close();
	output_inst_temp_h.close();
/*
	output_inst_chem_pot_e.close();
	output_inst_chem_pot_h.close();
*/
	//output_inversion_positive.close();
	
	/*
	output_scattering_out_e.close();
	output_scattering_in_e.close();
	output_scattering_out_h.close();
	output_scattering_in_h.close();
	*/
}

/* Save to file dynamic variables in order to restart from a savepoint 
 * save_count is the save identifier number
 * dev_num is the device identifier number
 * */
void Device::file_save_variables(int save_count, int dev_num)
{
	std::stringstream fileName;
	std::ofstream saveOut;

	double prev_tmp[18*number_K_points + 3*coulomb_matrix_red_size];
	for(int i = 0; i < number_K_points; i++)
	{
		prev_tmp[		     i] = carrier_scattering_prev_ne[i];
		prev_tmp[1*number_K_points + i] = carrier_scattering_prev_nh[i];
		prev_tmp[2*number_K_points + i] = renormalize_prev_ne[i];
		prev_tmp[3*number_K_points + i] = renormalize_prev_nh[i];
		prev_tmp[4*number_K_points + i] = coulomb_potential_epsilon_inv[i];
		prev_tmp[5*number_K_points + i] = fe_k_inst[i];
		prev_tmp[6*number_K_points + i] = fh_k_inst[i];

		prev_tmp[7*number_K_points  + i] = energy_renormalized_e[i];
		prev_tmp[8*number_K_points  + i] = energy_renormalized_h[i];
		prev_tmp[9*number_K_points  + i] = real(electric_field_renormalized[i]);
		prev_tmp[10*number_K_points + i] = imag(electric_field_renormalized[i]);

		prev_tmp[11*number_K_points + i] = real(p_k[i]);
		prev_tmp[12*number_K_points + i] = imag(p_k[i]);
		prev_tmp[13*number_K_points + i] = ne_k[i];
		prev_tmp[14*number_K_points + i] = nh_k[i];

		prev_tmp[15*number_K_points + i] = coulomb_potential_normalized_ee[i];
		prev_tmp[16*number_K_points + i] = coulomb_potential_normalized_hh[i];
		prev_tmp[17*number_K_points + i] = coulomb_potential_normalized_eh[i];
	}
	for(int i = 0; i < coulomb_matrix_red_size; i++)
	{
		prev_tmp[18*number_K_points + 				  i] = coulomb_matrix_ee_red[i];
		prev_tmp[18*number_K_points + 1*coulomb_matrix_red_size + i] = coulomb_matrix_hh_red[i];
		prev_tmp[18*number_K_points + 2*coulomb_matrix_red_size + i] = coulomb_matrix_eh_red[i];
	}
//	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_dynamic_vars.dat";
	saveBinary(fileName.str(), prev_tmp, 18*number_K_points + 3*coulomb_matrix_red_size);

/*
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_prev_scatt_e.dat";
	saveBinary(fileName.str(), carrier_scattering_prev_ne, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_prev_scatt_h.dat";
	saveBinary(fileName.str(), carrier_scattering_prev_nh, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_prev_renorm_e.dat";
	saveBinary(fileName.str(), renormalize_prev_ne, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_prev_renorm_h.dat";
	saveBinary(fileName.str(), renormalize_prev_nh, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_screen_eps_inv.dat";
	saveBinary(fileName.str(), coulomb_potential_epsilon_inv, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_fe.dat";
	saveBinary(fileName.str(), fe_k_inst, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_fh.dat";
	saveBinary(fileName.str(), fh_k_inst, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_energy_renormalized_e.dat";
	saveBinary(fileName.str(), energy_renormalized_e, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_energy_renormalized_h.dat";
	saveBinary(fileName.str(), energy_renormalized_h, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_electric_field_renormalized.dat";
	saveBinary(fileName.str(), electric_field_renormalized, number_K_points);

	// p_k real
//	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_p.dat";
	saveBinary(fileName.str(), p_k, number_K_points);

	// ne
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_ne.dat";
	saveBinary(fileName.str(), ne_k, number_K_points);
	
	// ne
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_nh.dat";
	saveBinary(fileName.str(), nh_k, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Vc_ee.dat";
	saveBinary(fileName.str(), coulomb_potential_normalized_ee, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Vc_hh.dat";
	saveBinary(fileName.str(), coulomb_potential_normalized_hh, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Vc_eh.dat";
	saveBinary(fileName.str(), coulomb_potential_normalized_eh, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Vc_ee_red.dat";
	saveBinary(fileName.str(), coulomb_matrix_ee_red, coulomb_matrix_red_size);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Vc_hh_red.dat";
	saveBinary(fileName.str(), coulomb_matrix_hh_red, coulomb_matrix_red_size);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Vc_eh_red.dat";
	saveBinary(fileName.str(), coulomb_matrix_eh_red, coulomb_matrix_red_size);
*/
	//-----


	#if defined(USE_ISAK_HOLE_FILLING)||defined(USE_ISAK_HOLE_FILLING_TABLE)

	double scatt_tmp[12*number_K_points];
	for(int i = 0; i < number_K_points; i++)
	{
		scatt_tmp[			i] = carrier_scattering_rates_e_in[i];
		scatt_tmp[1*number_K_points + 	i] = carrier_scattering_rates_e_out[i];
		scatt_tmp[2*number_K_points + 	i] = carrier_scattering_rates_h_in[i];
		scatt_tmp[3*number_K_points + 	i] = carrier_scattering_rates_h_out[i];
		scatt_tmp[4*number_K_points + 	i] = carrier_scattering_rates_phonon_e_in[i];
		scatt_tmp[5*number_K_points + 	i] = carrier_scattering_rates_phonon_e_out[i];
		scatt_tmp[6*number_K_points + 	i] = carrier_scattering_rates_phonon_h_in[i];
		scatt_tmp[7*number_K_points + 	i] = carrier_scattering_rates_phonon_h_out[i];
		scatt_tmp[8*number_K_points + 	i] = carrier_scattering_rates_e_total[i];
		scatt_tmp[9*number_K_points + 	i] = carrier_scattering_rates_h_total[i];
		scatt_tmp[10*number_K_points + 	i] = carrier_scattering_rates_phonon_e_total[i];
		scatt_tmp[11*number_K_points + 	i] = carrier_scattering_rates_phonon_h_total[i];
	}

	// Scattering
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_scattering_vars.dat";
	saveBinary(fileName.str(), scatt_tmp, 12*number_K_points);


/*
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cc_scatt_e_in.dat";
	saveBinary(fileName.str(), carrier_scattering_rates_e_in, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cc_scatt_e_out.dat";
	saveBinary(fileName.str(), carrier_scattering_rates_e_out, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cc_scatt_h_in.dat";
	saveBinary(fileName.str(), carrier_scattering_rates_h_in, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cc_scatt_h_out.dat";
	saveBinary(fileName.str(), carrier_scattering_rates_h_out, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cp_scatt_e_in.dat";
	saveBinary(fileName.str(), carrier_scattering_rates_phonon_e_in, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cp_scatt_e_out.dat";
	saveBinary(fileName.str(), carrier_scattering_rates_phonon_e_out, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cp_scatt_h_in.dat";
	saveBinary(fileName.str(), carrier_scattering_rates_phonon_h_in, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cp_scatt_h_out.dat";
	saveBinary(fileName.str(), carrier_scattering_rates_phonon_h_out, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cc_scatt_e_total.dat";
	saveBinary(fileName.str(), carrier_scattering_rates_e_total, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cc_scatt_h_total.dat";
	saveBinary(fileName.str(), carrier_scattering_rates_h_total, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cp_scatt_e_total.dat";
	saveBinary(fileName.str(), carrier_scattering_rates_phonon_e_total, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cp_scatt_h_total.dat";
	saveBinary(fileName.str(), carrier_scattering_rates_phonon_h_total, number_K_points);
*/

	#endif

	double misc_tmp[18];
	misc_tmp[0] = device_inst_temp_e_prev;
	misc_tmp[1] = device_inst_temp_h_prev;
	misc_tmp[2] = device_inst_temp_e;
	misc_tmp[3] = device_inst_temp_h;
	misc_tmp[4] = device_inst_chem_pot_e;
	misc_tmp[5] = device_inst_chem_pot_h;

	misc_tmp[6] = real(macroscopic_polarization);
	misc_tmp[7] = imag(macroscopic_polarization);
	misc_tmp[8] = real(macroscopic_polarization_m1);
	misc_tmp[9] = imag(macroscopic_polarization_m1);
	misc_tmp[10] = real(macroscopic_polarization_m2);
	misc_tmp[11] = imag(macroscopic_polarization_m2);

	misc_tmp[12] = real(electric_field);
	misc_tmp[13] = imag(electric_field);
	misc_tmp[14] = real(electric_field_tp05);
	misc_tmp[15] = imag(electric_field_tp05);
	misc_tmp[16] = real(electric_field_tp1);
	misc_tmp[17] = imag(electric_field_tp1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_misc_vars.dat";
	saveBinary(fileName.str(), misc_tmp, 18);

/*
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_temp_prev_e.dat";
	saveBinary(fileName.str(), &device_inst_temp_e_prev, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_temp_prev_h.dat";
	saveBinary(fileName.str(), &device_inst_temp_h_prev, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_temp_e.dat";
	saveBinary(fileName.str(), &device_inst_temp_e, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_temp_h.dat";
	saveBinary(fileName.str(), &device_inst_temp_h, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_chem_e.dat";
	saveBinary(fileName.str(), &device_inst_chem_pot_e, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_chem_h.dat";
	saveBinary(fileName.str(), &device_inst_chem_pot_h, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_macPol.dat";
	saveBinary(fileName.str(), &macroscopic_polarization, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_macPol_1.dat";
	saveBinary(fileName.str(), &macroscopic_polarization_m1, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_macPol_2.dat";
	saveBinary(fileName.str(), &macroscopic_polarization_m2, 1);

	// Electric field
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Et.dat";
	saveBinary(fileName.str(), &electric_field, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Et_05.dat";
	saveBinary(fileName.str(), &electric_field_tp05, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Et_1.dat";
	saveBinary(fileName.str(), &electric_field_tp1, 1);
*/

	#ifdef SBE_USE_SPONTAN_EMISSIONS
		#if defined(__ICC) || defined(__INTEL_COMPILER)
		//	VSLStreamStatePtr stream; // ICC or ICPC
			fileName.str(std::string());
			fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_rng_state.dat";
			vslSaveStreamF(stream,fileName.str().c_str());
		#elif defined(__GNUC__) || defined(__GNUG__)
		//	std::mt19937 device_generator; // G++
			fileName.str(std::string());
			fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_rng_state.dat";
			std::ofstream fout(fileName.str());
			fout << device_generator;
			fout.close();
		#endif
	#endif
	

	#ifdef USE_ISAK_HOLE_FILLING_TABLE
	
		// When saving variables, add points if at least 8 have been computed
		carrier_scattering_table_e_in->file_save();
		carrier_scattering_table_e_out->file_save();
		carrier_scattering_table_h_in->file_save();
		carrier_scattering_table_h_out->file_save();
	#endif
}

/* Load from files all variables from a savepoint
 * save_count is the save identifier number
 * dev_num is the device identifier number
 *  */
void Device::file_load_variables(int save_count, int dev_num)
{
	std::stringstream fileName;
	std::ifstream loadOut;

	double prev_tmp[18*number_K_points + 3*coulomb_matrix_red_size];
//	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_dynamic_vars.dat";
	loadBinary(fileName.str(), prev_tmp, 18*number_K_points + 3*coulomb_matrix_red_size);
	for(int i = 0; i < number_K_points; i++)
	{
		carrier_scattering_prev_ne[i] 	 	= prev_tmp[		       i];
		carrier_scattering_prev_nh[i] 	 	= prev_tmp[1*number_K_points + i];
		renormalize_prev_ne[i] 		 	= prev_tmp[2*number_K_points + i];
		renormalize_prev_nh[i]		 	= prev_tmp[3*number_K_points + i];
		coulomb_potential_epsilon_inv[i]	= prev_tmp[4*number_K_points + i];
		fe_k_inst[i]			 	= prev_tmp[5*number_K_points + i];
		fh_k_inst[i]			 	= prev_tmp[6*number_K_points + i];

		energy_renormalized_e[i]	 	= prev_tmp[7*number_K_points + i];
		energy_renormalized_h[i]	 	= prev_tmp[8*number_K_points + i];
		electric_field_renormalized[i]   	= prev_tmp[9*number_K_points + i]  + I*prev_tmp[10*number_K_points + i];

		
		p_k[i]  			 	= prev_tmp[11*number_K_points + i] + I*prev_tmp[12*number_K_points + i];
		ne_k[i] 			 	= prev_tmp[13*number_K_points + i];
		nh_k[i] 				= prev_tmp[14*number_K_points + i];

		coulomb_potential_normalized_ee[i] 	= prev_tmp[15*number_K_points + i];
		coulomb_potential_normalized_hh[i] 	= prev_tmp[16*number_K_points + i];
		coulomb_potential_normalized_eh[i] 	= prev_tmp[17*number_K_points + i];
	}
	for(int i = 0; i < coulomb_matrix_red_size; i++)
	{
		coulomb_matrix_ee_red[i] 		= prev_tmp[18*number_K_points + 			    i];
		coulomb_matrix_hh_red[i] 		= prev_tmp[18*number_K_points + 1*coulomb_matrix_red_size + i];
		coulomb_matrix_eh_red[i] 		= prev_tmp[18*number_K_points + 2*coulomb_matrix_red_size + i];
	}

/*
	// p_k
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_p.dat";
	loadBinary(fileName.str(), p_k, number_K_points);
	
	// ne
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_ne.dat";
	loadBinary(fileName.str(), ne_k, number_K_points);
	
	// nh
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_nh.dat";
	loadBinary(fileName.str(), nh_k, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_energy_renormalized_e.dat";
	loadBinary(fileName.str(), energy_renormalized_e, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_energy_renormalized_h.dat";
	loadBinary(fileName.str(), energy_renormalized_h, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_electric_field_renormalized.dat";
	loadBinary(fileName.str(), electric_field_renormalized, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Vc_ee.dat";
	loadBinary(fileName.str(), coulomb_potential_normalized_ee, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Vc_hh.dat";
	loadBinary(fileName.str(), coulomb_potential_normalized_hh, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Vc_eh.dat";
	loadBinary(fileName.str(), coulomb_potential_normalized_eh, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Vc_ee_red.dat";
	loadBinary(fileName.str(), coulomb_matrix_ee_red, coulomb_matrix_red_size);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Vc_hh_red.dat";
	loadBinary(fileName.str(), coulomb_matrix_hh_red, coulomb_matrix_red_size);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Vc_eh_red.dat";
	loadBinary(fileName.str(), coulomb_matrix_eh_red, coulomb_matrix_red_size);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_prev_scatt_e.dat";
	loadBinary(fileName.str(), carrier_scattering_prev_ne, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_prev_scatt_h.dat";
	loadBinary(fileName.str(), carrier_scattering_prev_nh, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_prev_renorm_e.dat";
	loadBinary(fileName.str(), renormalize_prev_ne, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_prev_renorm_h.dat";
	loadBinary(fileName.str(), renormalize_prev_nh, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_screen_eps_inv.dat";
	loadBinary(fileName.str(), coulomb_potential_epsilon_inv, number_K_points);
*/

	//-----

	#if defined(USE_ISAK_HOLE_FILLING)||defined(USE_ISAK_HOLE_FILLING_TABLE)
	double scatt_tmp[12*number_K_points];
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_scattering_vars.dat";
	loadBinary(fileName.str(), scatt_tmp, 12*number_K_points);
	for(int i = 0; i < number_K_points; i++)
	{
		carrier_scattering_rates_e_in[i]		= scatt_tmp[				i];
		carrier_scattering_rates_e_out[i] 		= scatt_tmp[1*number_K_points + 	i];
		carrier_scattering_rates_h_in[i] 		= scatt_tmp[2*number_K_points + 	i];
		carrier_scattering_rates_h_out[i] 		= scatt_tmp[3*number_K_points + 	i];
		carrier_scattering_rates_phonon_e_in[i] 	= scatt_tmp[4*number_K_points + 	i];
		carrier_scattering_rates_phonon_e_out[i] 	= scatt_tmp[5*number_K_points + 	i];
		carrier_scattering_rates_phonon_h_in[i] 	= scatt_tmp[6*number_K_points + 	i];
		carrier_scattering_rates_phonon_h_out[i] 	= scatt_tmp[7*number_K_points + 	i];
		carrier_scattering_rates_e_total[i] 		= scatt_tmp[8*number_K_points + 	i];
		carrier_scattering_rates_h_total[i] 		= scatt_tmp[9*number_K_points + 	i];
		carrier_scattering_rates_phonon_e_total[i] 	= scatt_tmp[10*number_K_points + 	i];
		carrier_scattering_rates_phonon_h_total[i] 	= scatt_tmp[11*number_K_points + 	i];
	}


/*
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cc_scatt_e_in.dat";
	loadBinary(fileName.str(), carrier_scattering_rates_e_in, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cc_scatt_e_out.dat";
	loadBinary(fileName.str(), carrier_scattering_rates_e_out, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cc_scatt_h_in.dat";
	loadBinary(fileName.str(), carrier_scattering_rates_h_in, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cc_scatt_h_out.dat";
	loadBinary(fileName.str(), carrier_scattering_rates_h_out, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cp_scatt_e_in.dat";
	loadBinary(fileName.str(), carrier_scattering_rates_phonon_e_in, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cp_scatt_e_out.dat";
	loadBinary(fileName.str(), carrier_scattering_rates_phonon_e_out, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cp_scatt_h_in.dat";
	loadBinary(fileName.str(), carrier_scattering_rates_phonon_h_in, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cp_scatt_h_out.dat";
	loadBinary(fileName.str(), carrier_scattering_rates_phonon_h_out, number_K_points);

	//--

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cc_scatt_e_total.dat";
	loadBinary(fileName.str(), carrier_scattering_rates_e_total, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cc_scatt_h_total.dat";
	loadBinary(fileName.str(), carrier_scattering_rates_h_total, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cp_scatt_e_total.dat";
	loadBinary(fileName.str(), carrier_scattering_rates_phonon_e_total, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_cp_scatt_h_total.dat";
	loadBinary(fileName.str(), carrier_scattering_rates_phonon_h_total, number_K_points);
*/
	#endif

	double misc_tmp[18];
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_misc_vars.dat";
	loadBinary(fileName.str(), misc_tmp, 18);
	device_inst_temp_e_prev 	= misc_tmp[0];
	device_inst_temp_h_prev 	= misc_tmp[1];
	device_inst_temp_e 		= misc_tmp[2];
	device_inst_temp_h 		= misc_tmp[3];
	device_inst_chem_pot_e 		= misc_tmp[4];
	device_inst_chem_pot_h 		= misc_tmp[5];

	macroscopic_polarization 	= misc_tmp[6]  + I*misc_tmp[7];
	macroscopic_polarization_m1 	= misc_tmp[8]  + I*misc_tmp[9];
	macroscopic_polarization_m2 	= misc_tmp[10] + I*misc_tmp[11];

	electric_field 			= misc_tmp[12] + I*misc_tmp[13];
	electric_field_tp05 		= misc_tmp[14] + I*misc_tmp[15];
	electric_field_tp1 		= misc_tmp[16] + I*misc_tmp[17];


	//----
/*
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_temp_prev_e.dat";
	loadBinary(fileName.str(), &device_inst_temp_e_prev, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_temp_prev_h.dat";
	loadBinary(fileName.str(), &device_inst_temp_h_prev, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_temp_e.dat";
	loadBinary(fileName.str(), &device_inst_temp_e, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_temp_h.dat";
	loadBinary(fileName.str(), &device_inst_temp_h, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_chem_e.dat";
	loadBinary(fileName.str(), &device_inst_chem_pot_e, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_chem_h.dat";
	loadBinary(fileName.str(), &device_inst_chem_pot_h, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_fe.dat";
	loadBinary(fileName.str(), fe_k_inst, number_K_points);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_inst_fh.dat";
	loadBinary(fileName.str(), fh_k_inst, number_K_points);

	// MacPol
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_macPol.dat";
	loadBinary(fileName.str(), &macroscopic_polarization, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_macPol_1.dat";
	loadBinary(fileName.str(), &macroscopic_polarization_m1, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_macPol_2.dat";
	loadBinary(fileName.str(), &macroscopic_polarization_m2, 1);

	// Electric field
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Et.dat";
	loadBinary(fileName.str(), &electric_field, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Et_05.dat";
	loadBinary(fileName.str(), &electric_field_tp05, 1);

	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_Et_1.dat";
	loadBinary(fileName.str(), &electric_field_tp1, 1);
*/

	#ifdef SBE_USE_SPONTAN_EMISSIONS
		#if defined(__ICC) || defined(__INTEL_COMPILER)
		//	VSLStreamStatePtr stream; // ICC or ICPC
			fileName.str(std::string());
			fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_rng_state.dat";
			vslDeleteStream(&stream);
			vslLoadStreamF(&stream,fileName.str().c_str());
			
		#elif defined(__GNUC__) || defined(__GNUG__)
		//	std::mt19937 device_generator; // G++
			fileName.str(std::string());
			fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_rng_state.dat";
			std::ifstream fin(fileName.str());
			fin >> device_generator;
			fin.close();
		#endif
	#endif

/*	
	// Macpol
	fileName.str(std::string());
	fileName << "save/save_" << save_count << "_DEV_" <<  dev_num << "_macPol.dat";
	loadBinary(fileName.str(), &macroscopic_polarization, 1);
*/
}

/* Load from files only carriers using the given save point 
 * out_count is the output file number
 * dev_num is the device identifier number
 * */
void Device::file_load_carriers(int out_count, int dev_num)
{
	std::stringstream baseName;
	baseName.str(std::string());
	baseName << getToFileOutputKey() << out_count;
	
	std::stringstream fileName;
	std::ifstream loadOut;
	/*
		Load dynamic variables from file
	*/
	double prev_tmp[18*number_K_points + 3*coulomb_matrix_red_size];
	fileName << "save/save_" << out_count << "_DEV_" <<  dev_num << "_dynamic_vars.dat";
	loadBinary(fileName.str(), prev_tmp, 18*number_K_points + 3*coulomb_matrix_red_size);

	for(int i = 0; i < number_K_points; i++)
	{
		p_k[i]  			 	= prev_tmp[11*number_K_points + i] + I*prev_tmp[12*number_K_points + i];
		ne_k[i] 			 	= prev_tmp[13*number_K_points + i];
		nh_k[i] 				= prev_tmp[14*number_K_points + i];
	}
	

/*
	// p_k
	fileName.str(std::string());
	fileName << "save/save_" << out_count << "_DEV_" <<  dev_num << "_p.dat";
	loadBinary(fileName.str(), p_k, number_K_points);
	
	// ne
	fileName.str(std::string());
	fileName << "save/save_" << out_count << "_DEV_" <<  dev_num << "_ne.dat";
	loadBinary(fileName.str(), ne_k, number_K_points);
	
	// nh
	fileName.str(std::string());
	fileName << "save/save_" << out_count << "_DEV_" <<  dev_num << "_nh.dat";
	loadBinary(fileName.str(), nh_k, number_K_points);
*/

	// Generate screening/scattering matrices
	double tmp = CARRIER_DEVIATION_TOL;
	CARRIER_DEVIATION_TOL = 0.0; // Force recalculation of screening
	renormalize_reducedCoulomb(p_k, ne_k, nh_k);
	CARRIER_DEVIATION_TOL = tmp;

	// Reset polarization
	for(int i =0 ; i < number_K_points; i++)
	{
		p_k[i] = 0.0;
	}
	//renormalize_reducedCoulomb(p_k, ne_k, nh_k); // Will ONLY recalcuate field_renormalize

/*	
	// Macpol
	fileName.str(std::string());
	fileName << "save/save_" << out_count << "_DEV_" <<  dev_num << "_macPol.dat";
	loadBinary(fileName.str(), &macroscopic_polarization, 1);

	// Electric field
	fileName.str(std::string());
	fileName << "save/save_" << out_count << "_DEV_" <<  dev_num << "_electric_field.dat";
	loadBinary(fileName.str(), &electric_field, 1);
*/

	//cout << "LOAD CARRIERS DONE! (" << out_count << ", " << dev_num << ")" << endl;
}

/* Load from files only carriers using an output file 
 * out_count is the output file number
 * line_num corresponds to the line number in the file
 * 
 * The target files correspond to the _ne_ and _nh_ files used for output
 * */
void Device::file_load_carriers_line(int out_count, int line_num)
{
	std::stringstream baseName;
	baseName.str(std::string());
	baseName << getToFileOutputKey() << out_count;
	
	std::stringstream fileName;
	std::ifstream loadOut;

	// p_k
	//fileName.str(std::string());
	//fileName << "save_" << save_count << "_DEV_" <<  dev_num << "_p.dat";
	//loadBinaryLine(fileName.str(), p_k, number_K_points, line_num);
	
	// ne
	fileName.str(std::string());
	fileName << baseName.str() << "_ne_" << getName() << ".dat";
	loadBinaryLine(fileName.str(), ne_k, number_K_points, line_num);
	
	// nh
	fileName.str(std::string());
	fileName << baseName.str() << "_nh_" << getName() << ".dat";
	loadBinaryLine(fileName.str(), nh_k, number_K_points, line_num);
	
	// Macpol
	//fileName.str(std::string());
	//fileName << "save_" << save_count << "_DEV_" <<  dev_num << "_macPol.dat";
	//loadBinaryLine(fileName.str(), &macroscopic_polarization, 1, line_num);
}

/* Load from files only carriers using an output file with a different name 
 * output_key is the base of the filenames to load from
 * out_count is the output file number
 * line_num corresponds to the line number in the file
 * 
 * The target files correspond to the _ne_ and _nh_ files used for output
 * */
void Device::file_load_carriers_line(const std::string &output_key, int out_count, int line_num)
{
	std::stringstream baseName;
	baseName.str(std::string());
	baseName << output_key << out_count;
	
	std::stringstream fileName;
	std::ifstream loadOut;

	// p_k
/*
	fileName.str(std::string());
	fileName << baseName.str() << "_p_" << getName() << ".dat";
	loadBinaryLine(fileName.str(), p_k, number_K_points, line_num);
*/
	// ne
	fileName.str(std::string());
	fileName << baseName.str() << "_ne_" << getName() << ".dat";
	loadBinaryLine(fileName.str(), ne_k, number_K_points, line_num);
	cout << "file_load_carriers_line: " << getName() << "  load from " << fileName.str() << endl;
	
	// nh
	fileName.str(std::string());
	fileName << baseName.str() << "_nh_" << getName() << ".dat";
	loadBinaryLine(fileName.str(), nh_k, number_K_points, line_num);
	cout << "file_load_carriers_line: " << getName() << "  load from " << fileName.str() << endl;
	
	// Macpol
	//fileName.str(std::string());
	//fileName << "save_" << save_count << "_DEV_" <<  dev_num << "_macPol.dat";
	//loadBinaryLine(fileName.str(), &macroscopic_polarization, 1, line_num);
}

/* Load from files only carriers using the ascii files stored in a directory 
 * device_num is the device identifier
 * line_num is the folder identifier
 * */
void Device::file_load_carriers_ascii(int device_num, int line_num)
{
	std::stringstream baseName;
	baseName.str(std::string());
	baseName << "neq_dist-line_" << line_num;
	
	if (!dirExists(baseName.str()))
	{
		cout << "file_load_carriers_ascii(): Could not find directory: " << baseName.str() << endl;
		exit(-1);
	}
	
	std::stringstream fileName;
	std::ifstream loadOut;
	double tmp;
	/*
		Load dynamic variables from file
	*/
	// p_k
	//fileName.str(std::string());
	//fileName << "save_" << save_count << "_DEV_" <<  dev_num << "_p.dat";
	//loadBinaryLine(fileName.str(), p_k, number_K_points, line_num);
	
	// ne
	fileName.str(std::string());
	fileName << baseName.str() << "/neq_dist_" << device_num << "_ne.txt";
	if (fileExists(fileName.str()))
	{
		loadOut.open(fileName.str().c_str(),std::ifstream::in);
		for(unsigned i = 0; i < number_K_points; i++)
		{
			loadOut >> tmp >> ne_k[i];
		}
		loadOut.close();
		cout << "file_load_carriers_ascii: " << getName() << "  load from " << fileName.str() << endl;
	} else {
		cout << "file_load_carriers_ascii: " << getName() << "  could NOT find file " << fileName.str() << endl;
	}
	
	// nh
	fileName.str(std::string());
	fileName << baseName.str() << "/neq_dist_" << device_num << "_nh.txt";
	if (fileExists(fileName.str()))
	{
		loadOut.open(fileName.str().c_str(),std::ifstream::in);
		for(unsigned i = 0; i < number_K_points; i++)
		{
			loadOut >> tmp >> nh_k[i];
		}
		loadOut.close();
		cout << "file_load_carriers_ascii: " << getName() << "  load from " << fileName.str() << endl;
	} else {
		cout << "file_load_carriers_ascii: " << getName() << "  could NOT find file " << fileName.str() << endl;;
	}
	
	// Macpol
	//fileName.str(std::string());
	//fileName << "save_" << save_count << "_DEV_" <<  dev_num << "_macPol.dat";
	//loadBinaryLine(fileName.str(), &macroscopic_polarization, 1, line_num);

}





//========================
// Misc helper functions

/* Integrate y[k]*K[k]*dK[0] using simpsons rule
 * */
double Device::misc_simpsons_quadrature_kdk(double *y)
{
	double total = 0.0;
	if (number_K_points % 2 == 0)
	{
		int SIMP_MAX = number_K_points;
		// Even
		total = y[0]*K[0] + y[SIMP_MAX-1]*K[SIMP_MAX-1];
		for(unsigned i = 1; i <= SIMP_MAX/2 -1; i++)
		{
			total += 2.0*y[2*i]*K[2*i];
		}
		for(unsigned i = 1; i <= SIMP_MAX/2; i++)
		{
			total += 4.0*y[2*i-1]*K[2*i-1];
		}
		total *= dK[0]/3.0; // Volume element

	} else {
		// Odd number in total
		int SIMP_MAX = number_K_points-1; // Use simpson on first N-1 points
		// Even
		total = y[0]*K[0] + y[SIMP_MAX-1]*K[SIMP_MAX-1];
		for(unsigned i = 1; i <= SIMP_MAX/2 -1; i++)
		{
			total += 2.0*y[2*i]*K[2*i];
		}
		for(unsigned i = 1; i <= SIMP_MAX/2; i++)
		{
			total += 4.0*y[2*i-1]*K[2*i-1];
		}
		total *= dK[0]/3.0; // Volume element

		// Integrate last point using Trapzoid
		total += y[number_K_points-1]*K[number_K_points-1];
	}
	total /= (Pi*a0*a0);
	return total;
}
/* Integrate Ek[k]*y[k]*K[k]*dK[0] using simpsons rule
 * */
double Device::misc_simpsons_quadrature_Ek_kdk(double *ek, double *y)
{
	double total = 0.0;
	if (number_K_points % 2 == 0)
	{
		int SIMP_MAX = number_K_points;
		// Even
		total = ek[0]*y[0]*K[0] + ek[SIMP_MAX-1]*y[SIMP_MAX-1]*K[SIMP_MAX-1];
		for(unsigned i = 1; i <= SIMP_MAX/2 -1; i++)
		{
			total += 2.0*ek[2*i]*y[2*i]*K[2*i];
		}
		for(unsigned i = 1; i <= SIMP_MAX/2; i++)
		{
			total += 4.0*ek[2*i-1]*y[2*i-1]*K[2*i-1];
		}
		total *= dK[0]/3.0; // Volume element

	} else {
		// Odd number in total
		int SIMP_MAX = number_K_points-1; // Use simpson on first N-1 points
		// Even
		total = ek[0]*y[0]*K[0] + ek[SIMP_MAX-1]*y[SIMP_MAX-1]*K[SIMP_MAX-1];
		for(unsigned i = 1; i <= SIMP_MAX/2 -1; i++)
		{
			total += 2.0*ek[2*i]*y[2*i]*K[2*i];
		}
		for(unsigned i = 1; i <= SIMP_MAX/2; i++)
		{
			total += 4.0*ek[2*i-1]*y[2*i-1]*K[2*i-1];
		}
		total *= dK[0]/3.0; // Volume element

		// Integrate last point using Trapzoid
		total += ek[number_K_points-1]*y[number_K_points-1]*K[number_K_points-1];
	}
	total /= (Pi*a0*a0);
	return total;
}

/* Return PolyLog[2,-x] to a given precision as found by adaptive Simpsons
 * */
double Device::misc_polyLog2(double x)
{
	double ERROR = 0.5e-4; 	// Requested error in evaluation integration
	return misc_adaptiveSimpsons(x,ERROR,50);
}

/* Argument for adaptiveSimpson routine for evaluation of misc_polyLog2()
 * */
double Device::misc_polyLog2_arg(double x)
{
	return -log(1.0+x)/x;
}

/* Auxiliary function for recursive misc_adaptiveSimpsons() below
 * */                                                                                                 
double Device::misc_adaptiveSimpsonsAux(double a, double b, double epsilon, double S, double fa, double fb, double fc, int bottom) 
{
	double c = (a + b)/2, h = b - a;                                                                  
	double d = (a + c)/2, e = (c + b)/2;                                                              
	double fd = misc_polyLog2_arg(d), fe = misc_polyLog2_arg(e);                                                                      
	double Sleft = (h/12)*(fa + 4*fd + fc);                                                           
	double Sright = (h/12)*(fc + 4*fe + fb);                                                          
	double S2 = Sleft + Sright;                                                                       
	if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)                                                    
	return S2 + (S2 - S)/15;                                                                        
	return misc_adaptiveSimpsonsAux(a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +                    
	misc_adaptiveSimpsonsAux(c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);                     
}         

/* Adaptive Simpson's main function
 * Integrate the function misc_polyLog2 from [0,b]
 * b -> Right endpoint
 * epsilon -> Error tolerance
 * maxRecursionDepth -> Maximum number of recursion calls
 * */
double Device::misc_adaptiveSimpsons(double b, double epsilon, int maxRecursionDepth)       
{
	double c = (0.0 + b)/2, h = b;                                                                  
	double fa = -1.0, fb = misc_polyLog2_arg(b), fc = misc_polyLog2_arg(c);                                                           
	double S = (h/6)*(fa + 4*fc + fb);                                                                
	return misc_adaptiveSimpsonsAux(0.0, b, epsilon, S, fa, fb, fc, maxRecursionDepth);                   
}    

/* Linear interpolation of function y-values stored in array
 * kp -> Where we want to evaluate the function
 * x0 -> The x-values of the function
 * dist -> The y-values of the function
 * numEl -> Number of elmeents in arrays x0, and dist
 * maxK_dist -> The value returned if kp > max(x0)
 * */
double Device::misc_interpolate_K_array(double kp, double *x0, double *dist, int numEl, double maxK_dist)
{

	// For when one tries to access at higher gridpoints
	if (kp >= x0[numEl-1])
	{
		cout << "WARNING: misc_interpolate_K_array():: k>x_end" << endl;
		exit(-1);
		return maxK_dist;
	}

	// For when one tries to access the lower gridpoints
	if (kp < x0[0])
	{
		cout << "WARNING: misc_interpolate_K_array():: k<x_0" << endl;
		exit(-1);
		return dist[0];
	}
	
	//double delta_x = x0[1]-x0[0]; // Assuming equidistant grids
	//double kp_indx = (kp-x0[0])/delta_x;
	//int i0 = floor(kp_indx);
	int i0 = misc_get_k_index(kp, x0, numEl); // In general grid
	int i1 = i0+1;

	return dist[i0] + (dist[i1]-dist[i0])*(kp-x0[i0])/(x0[i1]-x0[i0]);
}

// Same as misc_interpolate_K_array(), but also return the constants alpha, beta, and IND where: alpha*(f1 + beta*(f2-f1)) and f1 = f[IND]
double Device::misc_interpolate_K_array_index(double kp, double *x0, double *dist, int numEl, double maxK_dist, double *alpha, int *IND, double *beta)
{

	// For when one tries to access at higher gridpoints
	if (kp >= x0[numEl-1])
	{
		cout << "WARNING: misc_interpolate_K_array_index():: k>x_end" << endl;
		exit(-1);

	//	*IND   = 0;
	//	*alpha = 0.0;
	//	*beta  = 0.0;
	//	return maxK_dist;
	}

	// For when one tries to access the lower gridpoints
	if (kp < x0[0])
	{
		cout << "WARNING: misc_interpolate_K_array_index():: k<x_0" << endl;
		exit(-1);

	//	*IND = 0;
	//	*alpha = 1.0;
	//	*beta = 0.0;
	//	return dist[0];
	}
	
	//double delta_x = x0[1]-x0[0]; // Assuming equidistant grids
	//double kp_indx = (kp-x0[0])/delta_x;
	//int i0 = floor(kp_indx);
	int i0 = misc_get_k_index(kp, x0, numEl); // In general grid
	int i1 = i0+1;

	// Default settings
	*IND = i0;
	//*alpha = 1.0;
	*beta = (kp-x0[i0])/(x0[i1]-x0[i0]);

	return dist[i0] + (dist[i1]-dist[i0])*(kp-x0[i0])/(x0[i1]-x0[i0]);
}

// Same as misc_interpolate_K_array(), but also return the constants alpha, beta, and IND where: alpha*(f1 + beta*(f2-f1)) and f1 = f[IND]
double Device::misc_interpolate_K_array_linear_index(double kp, double *x0, double *dist, int numEl, int *IND, double *beta)
{

	// For when one tries to access at higher gridpoints
	if (kp > x0[numEl-1])
	{
		cout << "WARNING: misc_interpolate_K_array_index():: k>x_end" << endl;
		exit(-1);
	}

	// For when one tries to access the lower gridpoints
	if (kp < x0[0])
	{
		cout << "WARNING: misc_interpolate_K_array_index():: k<x_0" << endl;
		exit(-1);
	}
	
	int i0,i1;
	if (kp < x0[numEl-1])
	{
		//double delta_x = x0[1]-x0[0]; // Assuming equidistant grids
		//double kp_indx = (kp-x0[0])/delta_x;
		//i0 = floor(kp_indx);
		i0 = misc_get_k_index(kp, x0, numEl); // In general grid
		i1 = i0+1;
	} else {
		i0 = numEl-2;
		i1 = numEl-1;
	}

	// Default settings
	*IND = i0;
	//*alpha = 1.0;
	*beta = (kp-x0[i0])/(x0[i1]-x0[i0]);

	return dist[i0] + (dist[i1]-dist[i0])*(kp-x0[i0])/(x0[i1]-x0[i0]);
}

// Interpolate using parabolic approximation. Good for energy bands
double Device::misc_interpolate_K_array_parabolic(double kp, double *x0, double *dist, int numEl)
{

	// For when one tries to access at higher gridpoints
	if (kp > x0[numEl-1])
	{
		cout << "WARNING: misc_interpolate_K_array_parabolic():: k>x_end" << endl;
		exit(-1);
	}

	// For when one tries to access the lower gridpoints
	if (kp < x0[0])
	{
		cout << "WARNING: misc_interpolate_K_array_parabolic():: k<x_0" << endl;
		exit(-1);
	} 
	
	int i0,i1,i2;
	//double delta_x = x0[1]-x0[0]; // Assuming equidistant grids
	//double kp_indx = (kp-x0[0])/delta_x;
	//i1 = floor(kp_indx);
	i1 = misc_get_k_index(kp, x0, numEl); // In general grid
	i0 = i1-1;
	i2 = i1+1;
	
	if (i0 < 0)
	{
		// On left edge
		i0 += 1;
		i1 += 1;
		i2 += 1;
	} else if (i2 > numEl-1)
	{
		// On right edge
		i0 -= 1;
		i1 -= 1;
		i2 -= 1;
	}
	
	double y0 = dist[i0];
	double y1 = dist[i1];
	double y2 = dist[i2];
	double t0 = x0[i0];
	double t1 = x0[i1];
	double t2 = x0[i2];

	/*
	// Equidistant grid
	double delta_t = x0[1]-x0[0];
	double L0 =  0.5*(kp-t1)*(kp-t2);
	double L1 =     -(kp-t0)*(kp-t2);
	double L2 =  0.5*(kp-t0)*(kp-t1);
	
	return (y0*L0 + y1*L1 + y2*L2)/(delta_t*delta_t);
	*/
	
	// General grid
	double L0 = (kp-t1)*(kp-t2)/((t0-t1)*(t0-t2));
	double L1 = (kp-t0)*(kp-t2)/((t1-t0)*(t1-t2));
	double L2 = (kp-t0)*(kp-t1)/((t2-t0)*(t2-t1));
	return y0*L0 + y1*L1 + y2*L2;
}

// Interpolate using parabolic approximation. Good for energy bands
double Device::misc_interpolate_K_array_parabolic_index(double kp, double *x0, double *dist, int numEl, int *index, double *beta0, double *beta1, double *beta2)
{
	// For when one tries to access at higher gridpoints
	if (kp > x0[numEl-1])
	{
		cout << "WARNING: misc_interpolate_K_array_parabolic():: k>x_end" << endl;
		exit(-1);
	}

	// For when one tries to access the lower gridpoints
	if (kp < x0[0])
	{
		cout << "WARNING: misc_interpolate_K_array_parabolic():: k<x_0" << endl;
		exit(-1);
	} 
	
	int i0,i1,i2;
//	double delta_x = x0[1]-x0[0]; // Assuming equidistant grids
//	double kp_indx = (kp-x0[0])/delta_x;
//	i1 = floor(kp_indx);
	i1 = misc_get_k_index(kp, x0, numEl); // In general grid
	i0 = i1-1;
	i2 = i1+1;
	
	if (i0 < 0)
	{
		// On left edge
		i0 += 1;
		i1 += 1;
		i2 += 1;
	} else if (i2 > numEl-1)
	{
		// On right edge
		i0 -= 1;
		i1 -= 1;
		i2 -= 1;
	}
	
	double y0 = dist[i0];
	double y1 = dist[i1];
	double y2 = dist[i2];
	double t0 = x0[i0];
	double t1 = x0[i1];
	double t2 = x0[i2];
	
	// Equidistant grid
	/*
	double delta_t = x0[1]-x0[0];
	double inv_delta_t2 = 1.0/(delta_t*delta_t);
	double L0 =  0.5*(kp-t1)*(kp-t2);
	double L1 =     -(kp-t0)*(kp-t2);
	double L2 =  0.5*(kp-t0)*(kp-t1);

	// Interpolation works by giving the index i0 and constant beta1,beta2,beta3
	// From this one can find y0 = y[i0], y1 = y[i0+1], y2 = y[i0+2]
	// Then a given point is found by f(x) = beta0*y0 + beta1*y1 + beta2*y2
	*index = i0;
	*beta0 = L0*inv_delta_t2;
	*beta1 = L1*inv_delta_t2;
	*beta2 = L2*inv_delta_t2;
	
	return (y0*L0 + y1*L1 + y2*L2)*inv_delta_t2;
	*/

	// General grid
	double L0 = (kp-t1)*(kp-t2)/((t0-t1)*(t0-t2));
	double L1 = (kp-t0)*(kp-t2)/((t1-t0)*(t1-t2));
	double L2 = (kp-t0)*(kp-t1)/((t2-t0)*(t2-t1));


	// Interpolation works by giving the index i0 and constant beta1,beta2,beta3
	// From this one can find y0 = y[i0], y1 = y[i0+1], y2 = y[i0+2]
	// Then a given point is found by f(x) = beta0*y0 + beta1*y1 + beta2*y2
	*index = i0;
	*beta0 = L0;
	*beta1 = L1;
	*beta2 = L2;
	return y0*L0 + y1*L1 + y2*L2;
}
/* Return index 'i' of kp in array x such that x[i] <= kp < x[i+1]
 * Assumes array is ordered such that K[0] < K[1] < ... < K[numEl-1]
 * Does NOT check if kp is contained inside x[0] <= kp < x[numEl-1]
 */
int Device::misc_get_k_index(double kp, double *x, int numEl)
{
	int ia = 0;
	int ib = numEl-1;
	int ic = 0;
	int counter = 0;
	int counter_max = 10000;
	while (counter<counter_max)
	{
		ic = floor((ia+ib)/2.0);
		if ((x[ic] <= kp)&&(kp < x[ic+1]))
		{
			return ic;
		}
		
		if (x[ic] < kp)
		{
			ia = ic;
		} else {
			ib = ic;
		}

		counter++;
	}

	cout << "misc_get_k_index() Could not find element inside array!" << endl;
	cout << "kp = " << kp << endl;
	cout << "K[0] = " << K[0] << ", K[" << numEl-1 << "] = " << K[numEl-1] << endl;
	exit(-1);
}

/* Return the chemical potential for a 2D structure
 * The analytic formula can be found in:
 * Quamtum Theory of Optical and electronic properties of semiconductors
 * by: Haug and Koch
 * page 97
 * */
double Device::misc_get_chem_potential_analytic_2d(double dens, double mass, double temp)
{
	double beta = 1.0/(kB*temp);
	double cnt = hbar*hbar*beta*Pi*dens/mass;
	return log(exp(cnt) - 1.0)/beta;
}

/* Return the Fermi function at the given momentum value
 * for a species of given chemical potential, temperature, density and mass
 * */
double Device::misc_get_fermi_distribution(double temp, double chemPot, double mass, double waveNum)
{
	double kinetic_energy  = hbar*hbar*waveNum*waveNum/(2.0*mass*a0*a0);
	double exp_arg  = - ( chemPot - kinetic_energy)/(kB*temp);
	return 1.0/(exp(exp_arg) + 1.0);
	
}

// Set all dynamic variables to "Zero"
void Device::misc_zero_all_fields()
{
	// Zero electric field
	macroscopic_polarization = 0.0;
	macroscopic_polarization_m1 = 0.0;
	macroscopic_polarization_m2 = 0.0;
	electric_field = 0.0;
	electric_field_tp1 = 0.0;
	electric_field_tp05 = 0.0;

	device_inst_temp_e     = device_background_temp_e;
	device_inst_temp_h     = device_background_temp_h;
	device_inst_chem_pot_e = device_chemical_potential_e;
	device_inst_chem_pot_h = device_chemical_potential_h;
	
	device_inst_temp_e_prev = device_background_temp_e;
	device_inst_temp_h_prev = device_background_temp_h;
	for(unsigned i = 0; i < number_K_points; i++)
	{
		fe_k_inst[i] = fe_k[i];
		fh_k_inst[i] = fh_k[i];

		carrier_scattering_rate_approximation_e[i] = 0;
		carrier_scattering_rate_approximation_h[i] = 0;


		p_k[i] = 0.0;
		ne_k[i] = fe_k[i];
		nh_k[i] = fh_k[i];

		p_k_tmp[i] = 0;
		ne_k_tmp[i] = 0;
		nh_k_tmp[i] = 0;

		electric_field_renormalized[i] = 0.0;
		energy_renormalized_e[i] = 0.0;
		energy_renormalized_h[i] = 0.0;
	}
	for(unsigned i = 0; i < 4; i++)
	{
		for(unsigned j = 0; j < number_K_points; j++)
		{
			p_k_quad[i][j]  = 0;
			ne_k_quad[i][j] = 0;
			nh_k_quad[i][j] = 0;
		}
	}

	for(int i = 0; i < number_K_points; i++)
	{
		carrier_scattering_rates_e_in[i]  = 0.0;
		carrier_scattering_rates_e_out[i] = 0.0;
		carrier_scattering_rates_h_in[i]  = 0.0;
		carrier_scattering_rates_h_out[i] = 0.0;
		carrier_scattering_rates_e_total[i] = 0.0;
		carrier_scattering_rates_h_total[i] = 0.0;


		carrier_scattering_rates_phonon_e_in[i] = 0.0;
		carrier_scattering_rates_phonon_e_out[i] = 0.0;
		carrier_scattering_rates_phonon_e_total[i] = 0.0;
		carrier_scattering_rates_phonon_h_total[i] = 0.0;
		carrier_scattering_prev_ne[i] = 0.0;
		carrier_scattering_prev_nh[i] = 0.0;
		renormalize_prev_ne[i] = 0.0;
		renormalize_prev_nh[i] = 0.0;
	}

	#if defined(USE_ISAK_HOLE_FILLING_TABLE) || defined(USE_ISAK_HOLE_FILLING)


	// Calculation background scattering
	carrier_scattering_method2_parallel(ne_k, nh_k);
	#endif

	
	//=================================
	// Fill reduced Coulomb matrix
	#ifdef SBE_USE_FULL_SCREENING
	updateScreening_reducedCoulomb_full(ne_k, nh_k);
	#else
	updateScreening_reducedCoulomb_StaticPlasmonApproximation(ne_k, nh_k);
	#endif
	
	//==============================================
	// Calculate Energy and Field renormalization
	renormalize_reducedCoulomb(p_k, ne_k, nh_k);
	renormalize_reducedCoulomb_Energy(ne_k, nh_k); // Compute this once, in case it is not updated
}
