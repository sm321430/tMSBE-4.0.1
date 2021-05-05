
/*
	The DEVICE class is intended to be used on types of QW or ABSORBER.
	These objects contain the variable relevant for SBE solvers and such
	Isak Kilen @ 2016
*/

#ifndef __DEVICE_H_INCLUDED__
#define __DEVICE_H_INCLUDED__

#include <string>
#include <complex>
#include <fstream>
#include "classCyclic.h"
//#include "constantsAndMiscUnits.h"
#include "setup_simulation_variables.h"

#include "classIndexSet1d_coulomb.h" // For Coulomb matrix
#include "classIndexSet2d_screening.h" // For 2ndBorn code screening
#include "classIndexSet2d_phonon.h" // For 2ndBorn code screening
#include "classIndexSet3d_cc.h" // For 2ndBorn code indexing
#include "classDynamicInterpolation.cpp"

#ifdef MPI_BALANCE_WORKLOAD
	#include "myTimer.cpp"
#endif

#ifdef USE_DEVICE_TIMERS
	#include "statcenter.h"
#endif

#ifdef SBE_USE_SPONTAN_EMISSIONS
	#if defined(__ICC) || defined(__INTEL_COMPILER)
		#include "mkl_vsl.h"
	#elif defined(__GNUC__) || defined(__GNUG__)
		#include <random>
	#endif
#endif

class Device
{
	public:
		Device();
		~Device()
		{
			if (K != NULL)
			{
				for(int i =0; i < number_K_points; i++)
				{
					delete [] coulomb_matrix_ee[i];
					delete [] coulomb_matrix_hh[i];
					delete [] coulomb_matrix_eh[i];
				}
				delete [] coulomb_matrix_ee;
				delete [] coulomb_matrix_hh;
				delete [] coulomb_matrix_eh;
				delete [] K;

				delete [] dK;
				delete [] KdK;
				delete [] Ek_e;
				delete [] Ek_h;
				delete [] sumE_k;
				delete [] fe_k;
				delete [] fh_k;
				delete [] fe_k_inst;
				delete [] fh_k_inst;
				delete [] carrier_scattering_rate_approximation_e;
				delete [] carrier_scattering_rate_approximation_h;
				delete [] ne_k;
				delete [] nh_k;
				delete [] p_k;
				delete [] p_k_tmp;
				delete [] ne_k_tmp;
				delete [] nh_k_tmp;

				for(int i =0; i < 4; i++)
				{
					delete [] p_k_quad[i];
					delete [] ne_k_quad[i];
					delete [] nh_k_quad[i];
				}
				delete [] p_k_quad;
				delete [] ne_k_quad;
				delete [] nh_k_quad;

				delete p_k_ab4;
				delete ne_k_ab4;
				delete nh_k_ab4;
				delete [] electric_field_renormalized;
				delete [] energy_renormalized_e;
				delete [] energy_renormalized_h;
				delete [] coulomb_matrix_ee_red;
				delete [] coulomb_matrix_hh_red;
				delete [] coulomb_matrix_eh_red;
				delete [] coulomb_matrix_red_index;

				delete [] coulomb_potential_epsilon_inv;
				
				#ifdef SBE_USE_RESONANT_PUMP_MODEL
				delete [] device_pump_wk;
				#endif
				
				#if defined(USE_ISAK_HOLE_FILLING) || defined(USE_ISAK_HOLE_FILLING_TABLE)
				
				delete [] Q_sc;
				delete [] dQ_sc;
				delete [] Eta;
				delete [] dEta;
				
				delete [] theta_2B;
				delete [] dtheta_2B;
				delete [] Q;
				delete [] dQ;
				delete [] Kp;
				delete [] dKp;
				
				delete [] wavefunction_e;
				delete [] wavefunction_e_z;
				delete [] wavefunction_h;
				delete [] wavefunction_h_z;
				delete [] coulomb_potential_confinement_ee;
				delete [] coulomb_potential_confinement_hh;
				delete [] coulomb_potential_confinement_eh;
				delete [] coulomb_potential_normalized_ee;
				delete [] coulomb_potential_normalized_hh;
				delete [] coulomb_potential_normalized_eh;
				delete coulomb_potential_epsilon_index;
				delete carrier_scattering_index_ee;
				delete carrier_scattering_index_hh;
				delete carrier_scattering_index_eh;
				delete carrier_scattering_index_he;

				delete [] Q_ph;
				delete [] dQ_ph;
				delete [] Phi;
				delete [] dPhi;
				delete [] phonon_form_factor_ee;
				delete [] phonon_form_factor_hh;
				delete [] phonon_form_factor_Q;
				delete [] carrier_scattering_phonon_index_e1a;
				delete [] carrier_scattering_phonon_index_e1b;
				delete [] carrier_scattering_phonon_index_e2a;
				delete [] carrier_scattering_phonon_index_e2b;
				delete [] carrier_scattering_phonon_index_h1a;
				delete [] carrier_scattering_phonon_index_h1b;
				delete [] carrier_scattering_phonon_index_h2a;
				delete [] carrier_scattering_phonon_index_h2b;
				delete [] carrier_scattering_rates_phonon_e_in;
				delete [] carrier_scattering_rates_phonon_e_out;
				delete [] carrier_scattering_rates_phonon_h_in;
				delete [] carrier_scattering_rates_phonon_h_out;
				delete [] carrier_scattering_rates_phonon_e_total;
				delete [] carrier_scattering_rates_phonon_h_total;

				#endif
			
				#ifdef USE_ISAK_HOLE_FILLING_TABLE
				delete carrier_scattering_table_e_in;
				delete carrier_scattering_table_e_out;
				delete carrier_scattering_table_h_in;
				delete carrier_scattering_table_h_out;

				for(int i = 0; i < 8; i++)
				{
					delete [] carrier_scattering_table_gridpoints[i];
				}
				delete [] carrier_scattering_table_gridpoints;
				#endif

				#ifdef USE_HOLE_FILLING
				delete [] carrier_scattering_rates_e_in;
				delete [] carrier_scattering_rates_e_out;
				delete [] carrier_scattering_rates_h_in;
				delete [] carrier_scattering_rates_h_out;
				delete [] carrier_scattering_rates_e_total;
				delete [] carrier_scattering_rates_h_total;

				delete [] carrier_scattering_prev_ne;
				delete [] carrier_scattering_prev_nh;
				delete [] renormalize_prev_ne;
				delete [] renormalize_prev_nh;

				#endif

				#ifdef MPI_BALANCE_WORKLOAD
				delete MPI_load;
				#endif

				#ifdef SBE_USE_SPONTAN_EMISSIONS
				delete [] device_spont_emission_wk;
					#if defined(__ICC) || defined(__INTEL_COMPILER)
						vslDeleteStream(&stream);
					#else
						// for g++
					#endif
				#endif
			}			
		}

		void resetAllVariables(void);

		Device(const Device &obj)
		{
			resetAllVariables();

			device_name = obj.device_name;
			device_output_key = obj.device_output_key;
			device_position_x = obj.device_position_x;
			computationalDelay = obj.computationalDelay;
			analyticDelay = obj.analyticDelay;
			thresholdZeroFlag = obj.thresholdZeroFlag;
			
			macroscopic_polarization = obj.macroscopic_polarization;
			macroscopic_polarization_m1 = obj.macroscopic_polarization_m1;
			macroscopic_polarization_m2 = obj.macroscopic_polarization_m2;
			electric_field = obj.electric_field;
			electric_field_tp1  = obj.electric_field_tp1;
			electric_field_tp05 = obj.electric_field_tp05;
	
			device_focus_E = obj.device_focus_E;
			device_effective_qw = obj.device_effective_qw;
			device_deph_scale = obj.device_deph_scale;
			device_me = obj.device_me;
			device_mh = obj.device_mh;
			device_mr = obj.device_mr;
			device_Eb = obj.device_Eb;
			device_density = obj.device_density;
			device_bandgap = obj.device_bandgap;
			device_dipolemoment = obj.device_dipolemoment;
			device_dcv_hbar = obj.device_dcv_hbar;
			device_background_temp_e = obj.device_background_temp_e;
			device_background_temp_h = obj.device_background_temp_h;
			device_chemical_potential_e = obj.device_chemical_potential_e;
			device_chemical_potential_h = obj.device_chemical_potential_h;
			device_length_qw = obj.device_length_qw;
			device_deph_time = obj.device_deph_time;
			device_occ_relax_time = obj.device_occ_relax_time;
			device_occ_hole_relax_time = obj.device_occ_hole_relax_time;
			device_pol_frequency = obj.device_pol_frequency;
			K_max = obj.K_max;
			number_K_points = obj.number_K_points;

			coulomb_matrix_red_index_i = obj.coulomb_matrix_red_index_i;
			coulomb_matrix_red_index_j = obj.coulomb_matrix_red_index_j;

			device_inst_temp_e = obj.device_inst_temp_e;
			device_inst_temp_h = obj.device_inst_temp_h;
			device_inst_chem_pot_e = obj.device_inst_chem_pot_e;
			device_inst_chem_pot_h = obj.device_inst_chem_pot_h;
			device_inst_temp_e_prev = obj.device_inst_temp_e_prev;
			device_inst_temp_h_prev = obj.device_inst_temp_h_prev;

			SBE_TIME_SCALE = obj.SBE_TIME_SCALE;
			SBE_POL_DEPH = obj.SBE_POL_DEPH;
			SBE_OCC_PUMP = obj.SBE_OCC_PUMP;
			SBE_OCC_HOLE = obj.SBE_OCC_HOLE;

			ab4_precompute_counter = 0;
			ab4_prev_state = false;
			p_k_ab4 = NULL;
			ne_k_ab4 = NULL;
			nh_k_ab4 = NULL;
	
			hyst_tnext = obj.hyst_tnext;
			hyst_count = obj.hyst_count;
			hyst_dens_initial = obj.hyst_dens_initial;

			device_transverse_position = obj.device_transverse_position;
			device_transverse_pump_scale = obj.device_transverse_pump_scale;
			device_transverse_temp_scale = obj.device_transverse_temp_scale;

			#ifdef SBE_USE_RESONANT_PUMP_MODEL
			device_pump_W0 = obj.device_pump_W0;
			device_pump_E0 = obj.device_pump_E0;
			device_pump_ETA = obj.device_pump_ETA;
			device_pump_N0 = obj.device_pump_N0;

			if (obj.device_pump_wk != NULL)
			{
				device_pump_wk = new double[number_K_points];
			}
			#endif

			if (obj.K != NULL)
			{
				K  = new double[number_K_points]; 
				dK = new double[number_K_points];
				KdK = new double[number_K_points];
				fe_k = new double[number_K_points];
				fh_k = new double[number_K_points];
				Ek_e = new double[number_K_points];
				Ek_h = new double[number_K_points];
				sumE_k = new double[number_K_points];


				coulomb_matrix_ee = new double*[number_K_points];
				coulomb_matrix_hh = new double*[number_K_points];
				coulomb_matrix_eh = new double*[number_K_points];
				for(int i = 0; i < number_K_points; i++)
				{
					coulomb_matrix_ee[i] = new double[number_K_points];
					coulomb_matrix_hh[i] = new double[number_K_points];
					coulomb_matrix_eh[i] = new double[number_K_points];
				}


				coulomb_matrix_red_size = obj.coulomb_matrix_red_size;
				coulomb_matrix_ee_red =  new double[coulomb_matrix_red_size];
				coulomb_matrix_hh_red =  new double[coulomb_matrix_red_size];
				coulomb_matrix_eh_red =  new double[coulomb_matrix_red_size];
				coulomb_matrix_red_index = new indexSet1d_coulomb(*obj.coulomb_matrix_red_index);

				electric_field_renormalized = new std::complex<double>[number_K_points];
				energy_renormalized_e = new double[number_K_points];	
				energy_renormalized_h = new double[number_K_points];	

				fe_k_inst = new double[number_K_points];	
				fh_k_inst = new double[number_K_points];
				carrier_scattering_rate_approximation_e = new double[number_K_points];
				carrier_scattering_rate_approximation_h = new double[number_K_points];

				p_k = new std::complex<double>[number_K_points];
				ne_k = new double[number_K_points+2]; // 2 Extra point to speed up calculations in 2ndBorn code, should never be used elsewhere..
				nh_k = new double[number_K_points+2]; // 2 Extra point to speed up calculations in 2ndBorn code, should never be used elsewhere..
				
				p_k_tmp = new std::complex<double>[number_K_points];
				ne_k_tmp = new double[number_K_points+2]; // 2 Extra point to speed up calculations in 2ndBorn code, should never be used elsewhere..
				nh_k_tmp = new double[number_K_points+2]; // 2 Extra point to speed up calculations in 2ndBorn code, should never be used elsewhere..
	
				
				
				
				int NUM_TMP_STORAGE = 4;
				p_k_quad = new std::complex<double>*[NUM_TMP_STORAGE];
				for(unsigned i = 0; i < NUM_TMP_STORAGE; i++)
				{
					p_k_quad[i] = new std::complex<double>[number_K_points];
				}
				ne_k_quad = new double*[NUM_TMP_STORAGE];
				for(unsigned i = 0; i < NUM_TMP_STORAGE; i++)
				{
					ne_k_quad[i] = new double[number_K_points];
				}
				nh_k_quad = new double*[NUM_TMP_STORAGE];
				for(unsigned i = 0; i < NUM_TMP_STORAGE; i++)
				{
					nh_k_quad[i] = new double[number_K_points];
				}
			#ifdef SBE_USE_SPONTAN_EMISSIONS
				device_spont_emission_phasor = obj.device_spont_emission_phasor;
				device_spont_emission_wk = new double[number_K_points];
				device_spont_emission_on = obj.device_spont_emission_on;
			#endif
				CARRIER_DEVIATION_TOL = obj.CARRIER_DEVIATION_TOL;
				renormalize_prev_ne = new double[number_K_points];
				renormalize_prev_nh = new double[number_K_points];

				#ifdef USE_HOLE_FILLING
				carrier_scattering_rates_e_in = new double[number_K_points];
				carrier_scattering_rates_e_out = new double[number_K_points];
				carrier_scattering_rates_h_in = new double[number_K_points];
				carrier_scattering_rates_h_out = new double[number_K_points];
				carrier_scattering_rates_e_total = new double[number_K_points];
				carrier_scattering_rates_h_total = new double[number_K_points];

				CARRIER_DEVIATION_TOL = obj.CARRIER_DEVIATION_TOL;
				carrier_scattering_prev_ne = new double[number_K_points];
				carrier_scattering_prev_nh = new double[number_K_points];
				renormalize_prev_ne = new double[number_K_points];
				renormalize_prev_nh = new double[number_K_points];

				#endif

				coulomb_potential_epsilon_inv = new double[number_K_points];
			}


			
			number_Th_points = obj.number_Th_points;
			number_Z_points  = obj.number_Z_points;
			number_Q_points  = obj.number_Q_points;
			number_Kp_points = obj.number_Kp_points;
			number_Z_points  = obj.number_Z_points;
			
			number_Qsc_points = obj.number_Qsc_points;
			number_Eta_points = obj.number_Eta_points;
			if (obj.Q_sc != NULL)
			{
				Q_sc	= new double[number_Qsc_points];
				dQ_sc	= new double[number_Qsc_points];
				Eta		= new double[number_Eta_points];
				dEta	= new double[number_Eta_points];
			}


			if (obj.theta_2B != NULL)
			{
				theta_2B = new double[number_Th_points];
				dtheta_2B = new double[number_Th_points];
				Q = new double[number_Q_points];
				dQ = new double[number_Q_points];
				Kp = new double[number_Kp_points];
				dKp = new double[number_Kp_points];

				wavefunction_e = new double[number_Z_points];
				wavefunction_e_z = new double[number_Z_points];
				wavefunction_h = new double[number_Z_points];
				wavefunction_h_z = new double[number_Z_points];
				coulomb_potential_confinement_ee = new double[number_K_points];
				coulomb_potential_confinement_hh = new double[number_K_points];
				coulomb_potential_confinement_eh = new double[number_K_points];
				coulomb_potential_normalized_ee = new double[number_K_points];
				coulomb_potential_normalized_hh = new double[number_K_points];
				coulomb_potential_normalized_eh = new double[number_K_points];

				coulomb_potential_epsilon_index = new indexSet2d_screening(*obj.coulomb_potential_epsilon_index);
				carrier_scattering_index_ee = new indexSet3d_cc(*obj.carrier_scattering_index_ee);
				carrier_scattering_index_hh = new indexSet3d_cc(*obj.carrier_scattering_index_hh);
				carrier_scattering_index_eh = new indexSet3d_cc(*obj.carrier_scattering_index_eh);
				carrier_scattering_index_he = new indexSet3d_cc(*obj.carrier_scattering_index_he);

				number_Qph_points = obj.number_Qph_points;
				number_Phi_points = obj.number_Phi_points;
				
				Q_ph	= new double[number_Qph_points];
				dQ_ph	= new double[number_Qph_points];
				Phi	= new double[number_Phi_points];
				dPhi	= new double[number_Phi_points];
				
				phonon_form_factor_ee = new double[number_K_points];
				phonon_form_factor_hh = new double[number_K_points];
				phonon_form_factor_Q = new double[number_K_points];
				carrier_scattering_phonon_eps_zero = obj.carrier_scattering_phonon_eps_zero;
				carrier_scattering_phonon_eps_inf = obj.carrier_scattering_phonon_eps_inf;
				carrier_scattering_phonon_hbar_wLO = obj.carrier_scattering_phonon_hbar_wLO;
				carrier_scattering_phonon_temp = obj.carrier_scattering_phonon_temp;
				carrier_scattering_phonon_n    = obj.carrier_scattering_phonon_n;
				carrier_scattering_phonon_index_e1a = new indexSet2d_phonon(*obj.carrier_scattering_phonon_index_e1a);
				carrier_scattering_phonon_index_e1b = new indexSet2d_phonon(*obj.carrier_scattering_phonon_index_e1b);
				carrier_scattering_phonon_index_e2a = new indexSet2d_phonon(*obj.carrier_scattering_phonon_index_e2a);
				carrier_scattering_phonon_index_e2b = new indexSet2d_phonon(*obj.carrier_scattering_phonon_index_e2b);
				carrier_scattering_phonon_index_h1a = new indexSet2d_phonon(*obj.carrier_scattering_phonon_index_h1a);
				carrier_scattering_phonon_index_h1b = new indexSet2d_phonon(*obj.carrier_scattering_phonon_index_h1b);
				carrier_scattering_phonon_index_h2a = new indexSet2d_phonon(*obj.carrier_scattering_phonon_index_h2a);
				carrier_scattering_phonon_index_h2b = new indexSet2d_phonon(*obj.carrier_scattering_phonon_index_h2b);
				carrier_scattering_rates_phonon_e_in = new double[number_K_points];
				carrier_scattering_rates_phonon_e_out = new double[number_K_points];
				carrier_scattering_rates_phonon_h_in = new double[number_K_points];
				carrier_scattering_rates_phonon_h_out = new double[number_K_points];
				carrier_scattering_rates_phonon_e_total = new double[number_K_points];
				carrier_scattering_rates_phonon_h_total = new double[number_K_points];

				
			}

			if (obj.carrier_scattering_table_gridpoints != NULL)
			{
				carrier_scattering_table_e_in  = new DynamicInterpolation(*obj.carrier_scattering_table_e_in);
				carrier_scattering_table_e_out = new DynamicInterpolation(*obj.carrier_scattering_table_e_out);
				carrier_scattering_table_h_in  = new DynamicInterpolation(*obj.carrier_scattering_table_h_in);
				carrier_scattering_table_h_out = new DynamicInterpolation(*obj.carrier_scattering_table_h_out);
				
				carrier_scattering_table_gridpoints_numbers = obj.carrier_scattering_table_gridpoints_numbers;
				carrier_scattering_table_gridpoints = new double*[8];
				for(int i = 0; i < 8; i++)
				{
					carrier_scattering_table_gridpoints[i] = new double[6];
				}
			}
			#ifdef MPI_BALANCE_WORKLOAD
			if (MPI_load != NULL)
			{
				MPI_load = new myTimer(*obj.MPI_load);
			}
			#endif
		}
		
		void 		Print( void ) const;		// Print all device data to screen
		
		//======================================================
		// Functions to access various variables in the class
		
		std::complex<double> getMacroscopicPolarization(void) const;
		std::complex<double> getMacroscopicPolarization_m1(void) const;
		std::complex<double> getMacroscopicPolarization_m2(void) const;
		void setMacroscopicPolarization(std::complex<double>);
		void setMacroscopicPolarization_m1(std::complex<double>);
		void setMacroscopicPolarization_m2(std::complex<double>);
		std::complex<double> 	getElectricField(void) const;
		void 					setElectricField(std::complex<double>);
		std::complex<double> 	getElectricField_tp1(void) const;
		void 					setElectricField_tp1(std::complex<double>);
		std::complex<double> 	getElectricField_tp05(void) const;
		void 					setElectricField_tp05(std::complex<double>);
		
		
		void 		setName(const std::string &);		// Set name of device
		std::string getName(void) const;		// Return name of device
		void 		setToFileOutputKey(const std::string &);
		std::string getToFileOutputKey(void) const;
		double 	getPosition(void) const;		// Get 1D position
		void 		setPosition(double);		// Set 1D position
		void 		setDensity(double);		// Set density
		double 	getDensity(void) const;		// Get density
		void 		setBandgap(double);
		double 	getBandgap(void) const;

		double *getPointer_ne_k();
		double *getPointer_nh_k();
		std::complex<double> *getPointer_p_k();
		double *getPointer_s_ne_prev();
		double *getPointer_s_nh_prev();
		double *getPointer_cc_ne_prev();
		double *getPointer_cc_nh_prev();
		
		double getElectronMass(void) const;
		double getHoleMass(void) const;
		
		void setRelativeMass(void);
		double getRelativeMass(void) const;
		
		void setExitonBindEnergy(void);
		double getExitonBindEnergy(void) const;
		
		int getNumberKpoints(void) const;
		void setNumberKpoints(int);
		
		double getFocusE(void) const;
		void setFocusE(double);
		double getEffectiveQW(void) const;
		
		//===================
		// FILE IO functions

		// Read from .config files
		double file_getDoubleFromLine(std::string);
		int    	file_getIntFromLine(std::string);
		void 	file_removeWhitespace(std::string&);
		void 	file_readInConfigFromFile();
		
		// Output of simulation variables
		void file_output_write(int);
		void file_output_open(int,int);
		void file_output_close(int);
		
		// Save and load variables
		void file_save_variables(int, int);
		void file_load_variables(int, int);
		void file_load_carriers(int, int);
		void file_load_carriers_line(int, int);
		void file_load_carriers_line(const std::string &, int, int);
		void file_load_carriers_ascii(int, int);
		
		
		
		
		//=======================================
		// Various SBE functions
		
		void sbe_initialize(double);
		
		void sbe_delay_computations(double);
		void sbe_hysteresis_update_background_density(double);
		void sbe_set_spont_emissions(bool);
		void sbe_set_background_carrier_density(double);
		void sbe_set_real_pump_model(double, double, double, double);
		
		void sbe_RHS_pk(double, std::complex<double> *, std::complex<double> *, double *, double *,std::complex<double>);
		void sbe_RHS_ne(double, double *, std::complex<double> *, double *, double *,std::complex<double>);
		void sbe_RHS_nh(double, double *, std::complex<double> *, double *, double *,std::complex<double>);
		
		void sbe_RHS_pk_SmallAmp(double, std::complex<double> *, std::complex<double> *, double *, double *,std::complex<double>);
		void sbe_RHS_pk_debug(double, double, std::complex<double> *, std::complex<double> *, double *, double *,std::complex<double>);
		void sbe_RHS_ne_SmallAmp(double, double *, std::complex<double> *, double *, double *,std::complex<double>);
		void sbe_RHS_nh_SmallAmp(double, double *, std::complex<double> *, double *, double *,std::complex<double>);
		void sbe_RHS_pk_FieldFree(double, std::complex<double> *, std::complex<double> *, double *, double *);
		void sbe_RHS_ne_FieldFree(double, double *, std::complex<double> *, double *, double *);
		void sbe_RHS_nh_FieldFree(double, double *, std::complex<double> *, double *, double *);
		
		void sbe_iterate(double, double);
		void sbe_updateAnalytic(double);
		void sbe_iterate_rk4(double, double, bool, bool);
		void sbe_iterate_ab4(double, double, bool, bool);
		void sbe_iterate_abm4(double, double, bool, bool);
		
		void updateCoulombMatrix();
		void updateCoulombMatrixScreened();
		void updateScreening_reducedCoulomb_StaticPlasmonApproximation(double *, double *);
		void updateScreening_StaticPlasmonApproximation(double *, double *);
		void updateScreening_reducedCoulomb_full(double *, double *);
		void renormalize_Unscreened(std::complex<double> *, double*, double*);
		void zero_reducedCoulomb();
		void renormalize_reducedCoulomb(std::complex<double> *, double *, double *);
		void renormalize_reducedCoulomb_Energy(double *, double *);
		void renormalize_reducedCoulomb_Field(std::complex<double> *);
			
		//=======================================
		// Calculation of carrier scattering
		void carrier_scattering_init(double *, double *);
		void carrier_scattering_screening_full(double *, double *);
		void carrier_scattering_screening_full_fast(double *, double *);
		void carrier_scattering_method2_parallel(double *, double *);
		void carrier_scattering_method2_parallel_constant_rate(double *, double *);
		void carrier_scattering_method_old_full_serial(double *, double *);
		void carrier_scattering_method_old_reduced_serial(double *, double *);
		void carrier_scattering_calculate(long int, long int, double *, double *, bool);
		void carrier_scattering_method2_table(double *, double *);
		void carrier_scattering_rate_approx(double *, double *);
		
		void   carrier_scattering_rate_set_instant_temperature(double *, double *);
		void   carrier_scattering_rate_set_background_fermi();
		void   carrier_scattering_rate_get_instant_temp_and_chemPot(double *, double *, double, double, double, double);
		void   carrier_scattering_rate_get_instant_temp_and_chemPot2(double *, double *, double, double, double, double);
		void   carrier_scattering_rate_get_instant_temp_and_chemPot_argMin(double, double, double *, double *);
		void   carrier_scattering_rate_get_instant_temp_and_chemPot2_argMin(double, double, double, double *, double *);
		
		//=========================
		// Misc helper functions
		double misc_polyLog2(double);
		double misc_polyLog2_arg(double);
		double misc_simpsons_quadrature_kdk(double *);
		double misc_simpsons_quadrature_Ek_kdk(double *, double *);
		double misc_adaptiveSimpsonsAux(double a, double b, double epsilon,double S, double fa, double fb, double fc, int bottom);
		double misc_adaptiveSimpsons(double b,double epsilon,int maxRecursionDepth);
		double misc_interpolate_K_array(double, double *, double *, int,double);
		double misc_interpolate_K_array_index(double, double *, double *, int,double,double*, int *, double *);
		double misc_interpolate_K_array_linear_index(double kp, double *x0, double *dist, int numEl, int *IND, double *beta);
		double misc_interpolate_K_array_parabolic(double kp, double *x0, double *dist, int numEl);
		double misc_interpolate_K_array_parabolic_index(double kp, double *x0, double *dist, int numEl, int *, double *, double *, double *);
		int misc_get_k_index(double , double *, int);
		double misc_get_chem_potential_analytic_2d(double,double,double);
		double misc_get_fermi_distribution(double, double, double, double);

		void misc_zero_all_fields();

		void setTransversePosition(double);
		double getTransversePosition(void);
		void setTransverseBackgroundDensityScale(double);
		double getTransverseBackgroundDensityScale(void);
		void setTransverseBackgroundTemperatureScale(double);
		double get_MPI_load_mean_time(void);
		double get_MPI_load_mean_time_deviation(void);

	private:
		std::string device_name;		// Name of device
		std::string device_output_key;	// Start of output name of files
		double device_position_x;		// position of device, point position
		double computationalDelay;
		int analyticDelay;		//Counter for analytic delay thresholding
		int thresholdZeroFlag;	
	
		std::complex<double> macroscopic_polarization;	// Macroscopic polarization
		std::complex<double> macroscopic_polarization_m1;	// Macroscopic polarization previous timestep
		std::complex<double> macroscopic_polarization_m2;	// Macroscopic polarization previous timestep 2
		std::complex<double> electric_field;			// Electric field from Maxwell
		std::complex<double> electric_field_tp1;		// Electric field from Maxwell at t=t+dt
		std::complex<double> electric_field_tp05;		// Electric field from Maxwell at t=t+dt/2
		
		double device_focus_E;			// When focusing E on the device (Real number)
		double device_effective_qw;		// The feedback is multiplied to create effective QW's
		double device_deph_scale;		// The scaling of dephasing time =1 -> dep = 2.0Ebg

		double device_me;				// Electron mass
		double device_mh;				// Hole mass
		double device_mr;				// Relative mass
		double device_Eb;				// Excitonic binding energy
		double device_density;			// Density of device
		double device_bandgap;			// Bandgap energy
		double device_dipolemoment;		// Dipole moment in device
		double device_dcv_hbar;			// Dipole moment divided by hbar
		double device_background_temp_e;	// Background e-temperature in Kelvin
		double device_background_temp_h;	// Background h-temperature in Kelvin
		double device_chemical_potential_e;
		double device_chemical_potential_h;
		double device_length_qw;		// Length of QW (Approximated by single point else)
		
		double device_deph_time;		// Relaxation time for p_k -> 0
		double device_occ_relax_time;	// Relaxation time for n_k -> f_k
		double device_occ_hole_relax_time;	// Relaxation time for n_k -> f_k, instantanious
		double device_pol_frequency;	// Frequency of polarization vector p_k
		
		double 				K_max;		// Maximal wave number in units of a0^-1
		int 	  	number_K_points;		// Number of wave numbers
		double 				   *K;		// Wave number vector (midpoints in units of a0^-1) length is number_K_points
		double				  *dK;		// Wave number density vector (widths of intervals in units of a0^-1) length is number_K_points
		double				  *KdK;		// K[i]*dK[i]
		
		double 				*fe_k;		// Fermi-e distribution
		double 				*fh_k;		// Fermi-h distribution
		
		double 				*Ek_e;		// Kinetic energy electrons
		double 				*Ek_h;		// Kinetic energy holes
		double 			  *sumE_k;		// Sum of (E_g + Ek_e + Ek_h)/hbar
		
		double **coulomb_matrix_ee;
		double **coulomb_matrix_hh;
		double **coulomb_matrix_eh;
		int coulomb_matrix_red_size;
		double   *coulomb_matrix_ee_red;
		double   *coulomb_matrix_hh_red;
		double   *coulomb_matrix_eh_red;
		int 		coulomb_matrix_red_index_i; // Reduced k-grid length
		int 		coulomb_matrix_red_index_j; // Reduced k-grid length
		indexSet1d_coulomb *coulomb_matrix_red_index; // Index set for reduced Vc
		std::complex<double> *electric_field_renormalized;	// Renormalized electric field (scaled by 1/hbar)
		double  *energy_renormalized_e;	// Re-normalized energy (scaled by 1/hbar)		
		double  *energy_renormalized_h;	// Re-normalized energy (scaled by 1/hbar)		

		// Carrier scattering rate approximation variables
		double device_inst_temp_e;
		double device_inst_temp_h;
		double device_inst_chem_pot_e;
		double device_inst_chem_pot_h;
		double *fe_k_inst;
		double *fh_k_inst;
		double device_inst_temp_e_prev;
		double device_inst_temp_h_prev;
		double *carrier_scattering_rate_approximation_e;
		double *carrier_scattering_rate_approximation_h;

		
		
		// Dynamic variables
		std::complex<double> *p_k;		// Microscopic polarization
		double 				*ne_k;		// Electron occupation numbers
		double 				*nh_k;		// Hole occupation numbers
		
		
		// TMP STORAGE
		std::complex<double> *p_k_tmp;
		double 				*ne_k_tmp;
		double 				*nh_k_tmp;
		
		// Integrator storage space (RK4 gives 4 of each)
		std::complex<double> **p_k_quad;
		double 				**ne_k_quad;
		double				**nh_k_quad;
		
		// Integrator storage space (Adam-Bashford 4th order)
		int ab4_precompute_counter;
		bool ab4_prev_state;
		Cyclic< std::complex<double> > *p_k_ab4;
		Cyclic<double> *ne_k_ab4;
		Cyclic<double> *nh_k_ab4;
		
		// Scaling of equations
		double SBE_TIME_SCALE;
		double SBE_POL_DEPH;
		double SBE_OCC_PUMP;
		double SBE_OCC_HOLE;
		
		// Output variables to file
		std::ofstream output_E_real;
		std::ofstream output_E_imag;
		std::ofstream output_macpol_re;
		std::ofstream output_macpol_im;
		std::ofstream output_p_k_re;
		std::ofstream output_p_k_im;
		std::ofstream output_ne_k;
		std::ofstream output_nh_k;
		std::ofstream output_fe_k;
		std::ofstream output_fh_k;
		std::ofstream output_energy_renormalized_e_k;
		std::ofstream output_energy_renormalized_h_k;
		std::ofstream output_Nsum_e;
		std::ofstream output_Nsum_h;
		std::ofstream output_Esum_e;
		std::ofstream output_Esum_h;
		std::ofstream output_Fsum_e;
		std::ofstream output_Fsum_h;
		std::ofstream output_inst_temp_e;
		std::ofstream output_inst_temp_h;
		std::ofstream output_inst_chem_pot_e;
		std::ofstream output_inst_chem_pot_h;
		std::ofstream output_inversion_positive;
		std::ofstream output_scattering_out_e;
		std::ofstream output_scattering_in_e;
		std::ofstream output_scattering_out_h;
		std::ofstream output_scattering_in_h;
		std::ofstream output_scattering_total_e;
		std::ofstream output_scattering_total_h;
		std::ofstream output_scattering_total_phonon_e;
		std::ofstream output_scattering_total_phonon_h;
		std::ofstream output_scattering_Rabi;
		std::ofstream output_scattering_Pump_e;
		std::ofstream output_scattering_Pump_h;
 		std::ofstream output_scattering_Spont_emission;
		
		// Screening variables
		int number_Qsc_points;
		int number_Eta_points; // Screening angular grid
		double *Q_sc;
		double *dQ_sc;
		double *Eta;
		double *dEta;

		// Carrier-carrier scattering integrals
		double *theta_2B; // Angular integration points
		double *dtheta_2B; // Angular integration interval lengths
		int number_Th_points;
		int number_Z_points;
		int number_Q_points;
		int number_Kp_points;
		double *Q;
		double *dQ;
		double *Kp;
		double *dKp;
		
		
		double *wavefunction_e; // wavefunction for single particle
		double *wavefunction_e_z; // wavefunction for single particle
		double *wavefunction_h; // wavefunction for single particle
		double *wavefunction_h_z; // wavefunction for single particle
		double *coulomb_potential_confinement_ee; // Confinement functions for ee,eh and hh.
		double *coulomb_potential_confinement_hh; // Confinement functions for ee,eh and hh.
		double *coulomb_potential_confinement_eh; // Confinement functions for ee,eh and hh.
		double *coulomb_potential_epsilon_inv; // Screening of coulomb potential, V_screened = V/coulomb_epsilon_inv.
		double *coulomb_potential_normalized_ee; // Screened coulomb potential divided by a0*e*e/(2.0*eps0*eps).
		double *coulomb_potential_normalized_hh; // Screened coulomb potential divided by a0*e*e/(2.0*eps0*eps).
		double *coulomb_potential_normalized_eh; // Screened coulomb potential divided by a0*e*e/(2.0*eps0*eps).
		indexSet2d_screening *coulomb_potential_epsilon_index; // Index set for the coulomb screening.
		indexSet3d_cc *carrier_scattering_index_ee; // Index set of the carrier scattering indices e-e and h-h
		indexSet3d_cc *carrier_scattering_index_hh; // Index set of the carrier scattering indices e-e and h-h
		indexSet3d_cc *carrier_scattering_index_eh; // Index set of the carrier scattering indices e-h
		indexSet3d_cc *carrier_scattering_index_he; // Index set of the carrier scattering indices h-e
		double *carrier_scattering_rates_e_in;
		double *carrier_scattering_rates_e_out;
		double *carrier_scattering_rates_h_in;
		double *carrier_scattering_rates_h_out;
		double *carrier_scattering_rates_e_total;
		double *carrier_scattering_rates_h_total;

		DynamicInterpolation *carrier_scattering_table_e_in;
		DynamicInterpolation *carrier_scattering_table_e_out;
		DynamicInterpolation *carrier_scattering_table_h_in;
		DynamicInterpolation *carrier_scattering_table_h_out;
		double **carrier_scattering_table_gridpoints;
		int carrier_scattering_table_gridpoints_numbers;

		// Carrier-phonon scattering
		double carrier_scattering_phonon_eps_zero; // dielectric constant at low frequency
		double carrier_scattering_phonon_eps_inf;  // dielectric constant at high frequency
		double carrier_scattering_phonon_hbar_wLO; // [LO-Phonon-energy ] 
		double carrier_scattering_phonon_temp; // Phonon temperature
		double carrier_scattering_phonon_n; // Phonon distribution
		double *Q_ph;
		double *dQ_ph;
		double *Phi;
		double *dPhi;
		int number_Qph_points;
		int number_Phi_points;
		double *phonon_form_factor_ee; // Form factor for ee and hh.
		double *phonon_form_factor_hh; // Form factor for ee and hh.
		double *phonon_form_factor_Q; // Phonon Form factor Q grid
		indexSet2d_phonon *carrier_scattering_phonon_index_e1a; // Index set for the electron-phonon scattering
		indexSet2d_phonon *carrier_scattering_phonon_index_e1b; // Index set for the electron-phonon scattering
		indexSet2d_phonon *carrier_scattering_phonon_index_e2a; // Index set for the electron-phonon scattering
		indexSet2d_phonon *carrier_scattering_phonon_index_e2b; // Index set for the electron-phonon scattering
		indexSet2d_phonon *carrier_scattering_phonon_index_h1a; // Index set for the hole-phonon scattering
		indexSet2d_phonon *carrier_scattering_phonon_index_h1b; // Index set for the hole-phonon scattering
		indexSet2d_phonon *carrier_scattering_phonon_index_h2a; // Index set for the hole-phonon scattering
		indexSet2d_phonon *carrier_scattering_phonon_index_h2b; // Index set for the hole-phonon scattering
		double *carrier_scattering_rates_phonon_e_in;
		double *carrier_scattering_rates_phonon_e_out;
		double *carrier_scattering_rates_phonon_h_in;
		double *carrier_scattering_rates_phonon_h_out;
		double *carrier_scattering_rates_phonon_e_total;
		double *carrier_scattering_rates_phonon_h_total;

		// Storage for when skipping scattering calculation
		double CARRIER_DEVIATION_TOL;
		double *carrier_scattering_prev_ne;
		double *carrier_scattering_prev_nh;
		double *renormalize_prev_ne;
		double *renormalize_prev_nh;
		
		// Hysteresis variables
		double hyst_tnext;
		double hyst_count;
		double hyst_dens_initial;

		// Resonant pumping variables
		double device_pump_W0;	// Central frequency
		double device_pump_E0; 	// Field contributing to Carriers
		double device_pump_ETA; // Pump width
		double device_pump_N0;  // Initial density
		double *device_pump_wk; // Pump frequencies

		// Information about transverse domain
		double device_transverse_position;
		double device_transverse_pump_scale;
		double device_transverse_temp_scale;

		// Storage for MPI_load balancer
		#ifdef MPI_BALANCE_WORKLOAD
		myTimer *MPI_load;
		#endif

		#ifdef SBE_USE_SPONTAN_EMISSIONS
			std::complex<double> device_spont_emission_phasor;
			double *device_spont_emission_wk;
			bool device_spont_emission_on; // By default = true
		
			#if defined(__ICC) || defined(__INTEL_COMPILER)
				VSLStreamStatePtr stream; // ICC or ICPC
			#elif defined(__GNUC__) || defined(__GNUG__)
				std::mt19937 device_generator; // G++
			#endif
		#endif
};

/* Used to turn off spontanious emissions for devices
 * By default on if compiled with -D SBE_USE_SPONTAN_EMISSIONS
 * */
inline void Device::sbe_set_spont_emissions(bool newB)
{
	#ifdef SBE_USE_SPONTAN_EMISSIONS
	device_spont_emission_on = newB;
	if (device_spont_emission_on == false)
	{
		if (device_spont_emission_wk != NULL)
		{
			for (int i = 0; i < number_K_points; i++)
			{
				device_spont_emission_wk[i] = 0.0;
			}
		}
	}
	#endif
}

inline void Device::setName(const std::string & newName)
{
	device_name = newName;
}

inline std::string Device::getName(void) const
{
	return device_name;
}

inline void Device::setToFileOutputKey(const std::string & newName)
{
	device_output_key = newName;
}
inline std::string Device::getToFileOutputKey(void) const
{
	return device_output_key;
}

inline std::complex<double> Device::getMacroscopicPolarization(void) const
{
	return macroscopic_polarization;
}
inline std::complex<double> Device::getMacroscopicPolarization_m1(void) const
{
	return macroscopic_polarization_m1;
}
inline std::complex<double> Device::getMacroscopicPolarization_m2(void) const
{
	return macroscopic_polarization_m2;
}

inline void Device::setMacroscopicPolarization(std::complex<double> target)
{
	macroscopic_polarization = target;
}
inline void Device::setMacroscopicPolarization_m1(std::complex<double> target)
{
	macroscopic_polarization_m1 = target;
}
inline void Device::setMacroscopicPolarization_m2(std::complex<double> target)
{
	macroscopic_polarization_m2 = target;
}

inline std::complex<double> Device::getElectricField(void) const
{
	return electric_field;
}

inline void Device::setElectricField(std::complex<double> target)
{
	electric_field = target;
}

inline std::complex<double> Device::getElectricField_tp1(void) const
{
	return electric_field_tp1;
}

inline void Device::setElectricField_tp1(std::complex<double> target)
{
	electric_field_tp1 = target;
}

inline std::complex<double> Device::getElectricField_tp05(void) const
{
	return electric_field_tp05;
}

inline void Device::setElectricField_tp05(std::complex<double> target)
{
	electric_field_tp05 = target;
}

inline double Device::getPosition(void) const
{
	return device_position_x;
}

inline  void Device::setPosition(double newPosition_x)
{
	device_position_x = newPosition_x;
}

inline void Device::setDensity(double newDens)
{
	device_density = newDens;
}

inline double Device::getDensity() const
{
	return device_density;
}

inline double Device::getBandgap() const
{
	return device_bandgap;
}

inline void Device::setBandgap(double newEg)
{
	device_bandgap = newEg;
}

inline double Device::getElectronMass(void) const
{
	return device_me;
}

inline double Device::getHoleMass(void) const
{
	return device_mh;
}

inline void Device::setRelativeMass(void)
{
	device_mr = getElectronMass()*getHoleMass()/(getElectronMass() + getHoleMass());
}

inline double Device::getRelativeMass(void) const
{
	return device_mr;
}

inline void Device::setExitonBindEnergy(void)
{
	device_Eb = hbar*hbar/(2.0*getRelativeMass()*a0*a0);
}

inline double Device::getExitonBindEnergy(void) const
{
	return device_Eb;
}

inline int Device::getNumberKpoints(void) const
{
	return number_K_points;
}

inline void Device::setNumberKpoints(int target)
{
	number_K_points = target;
}

inline double Device::getFocusE(void) const
{
	return device_focus_E;
}

inline void Device::setFocusE(double newFocus)
{
	device_focus_E = newFocus;
}

inline double Device::getEffectiveQW(void) const
{
	return device_effective_qw;
}
inline double *Device::getPointer_ne_k(void)
{
	return ne_k;
}
inline double *Device::getPointer_nh_k(void)
{
	return nh_k;
}
inline std::complex<double> *Device::getPointer_p_k(void)
{
	return p_k;
}

inline double *Device::getPointer_s_ne_prev(void)
{
	return renormalize_prev_ne;
}
inline double *Device::getPointer_s_nh_prev(void)
{
	return renormalize_prev_nh;
}
inline double *Device::getPointer_cc_ne_prev(void)
{
	return carrier_scattering_prev_ne;
}
inline double *Device::getPointer_cc_nh_prev(void)
{
	return carrier_scattering_prev_nh;
}

inline void Device::setTransversePosition(double y_p)
{
	device_transverse_position = y_p;
}
inline double Device::getTransversePosition(void)
{
	return device_transverse_position;
}
inline void Device::setTransverseBackgroundDensityScale(double newP)
{
	device_transverse_pump_scale = newP;
}
inline double Device::getTransverseBackgroundDensityScale(void)
{
	return device_transverse_pump_scale;
}
inline void Device::setTransverseBackgroundTemperatureScale(double newP)
{
	device_transverse_temp_scale = newP;
}
inline double Device::get_MPI_load_mean_time()
{
	#ifdef MPI_BALANCE_WORKLOAD
	return MPI_load->get_real_time_mean();
	#else
	return 0;
	#endif
}

inline double Device::get_MPI_load_mean_time_deviation()
{
	#ifdef MPI_BALANCE_WORKLOAD
	return MPI_load->get_real_time_deviation();
	#else
	return 0;
	#endif

}


#endif
