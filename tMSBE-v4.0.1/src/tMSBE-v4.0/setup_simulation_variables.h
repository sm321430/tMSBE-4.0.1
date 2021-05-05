

#ifndef __SETUP_SIMULATION_VARIABLES_H_INCLUDED__
#define __SETUP_SIMULATION_VARIABLES_H_INCLUDED__

#include "constantsAndMiscUnits.h"


/*=================================================
 *  This file contains variables that define the
 *  simulation and will not change during runtime
 * */
 
// Should assume about 10% overhead while using timers
//#define USE_DEVICE_TIMERS
#define USE_MAIN_TIMERS
#ifdef USE_MAIN_TIMERS
	#include "myTimer.cpp"
	extern myTimerCentral *MainStat;
#endif

//#define MPI_BALANCE_WORKLOAD // experimental feature that attempts to estimate runtime of MPI nodes and use these times to optimize performance. ATM this is not better than a simpler estimate based on pump density.

#define USE_OPENMP		// USE OPENMP in general
#ifdef USE_OPENMP
	#include <omp.h>

	// Number of QWs to compute in parallel on a motherboard (=1 means each QW is computed in sequence)
	// OMP_THREADS_LEVEL_1 * OMP_THREADS_LEVEL_2 should NOT exceed the MAX number of threads for a node
	#define OMP_THREADS_LEVEL_1 2 

	// Number of threads to use INSIDE a single QW (=1 means that each QW is internally computed in serial)
	// OMP_THREADS_LEVEL_1 * OMP_THREADS_LEVEL_2 should NOT exceed the MAX number of threads for a node
	#define OMP_THREADS_LEVEL_2 1
#endif

//#define COMPILE_REFLECTION_CALCULATION // Use to disable QW carrier relaxation.

#define DUAL_PROP //Use to enable dual prop cavity which consists solely of twoArmStructures
//#define USE_EXPANDED_SBE  // Use to enable second order occupation numbers and third order polarization equations.
//#define CYCLIC_BOUNDARIES // Use to enable reflection from one boundary to another for cyclic operation i.e. ring cavities.
//#define TRANS_DELAY 	    // Use to enable transverse pulse delay for non-normal incidence on twoArmStructures. 

//#define SBE_CONFINEMENT_USE_FILE // Use confinment functions from files wavefunction_e_#.dat and wavefunction_h_#.dat
//#define SBE_USE_FULL_SCREENING	// Use Lindhard Formula, otherwise a static approximation is made

//#define SBE_UPDATE_SCREENING // Calculate Screening at every timestep, or only ONCE
#define SBE_UPDATE_BG_RENORM // Calculate Bg renormalization at every timestep, or only ONCE 
#define USE_HOLE_FILLING // WARNING: Dual rate hole filling WITH thresholding will underestimate the density. Can USE this to sim Spon. Emission.
//#define USE_ISAK_HOLE_FILLING
//#define USE_ISAK_HOLE_FILLING_TABLE

//#define SBE_USE_RESONANT_PUMP_MODEL // Pumping inverted carrier directly. Note: Remember to change the coefficients in the device for your simulation.
#define SBE_USE_SPONTAN_EMISSIONS // Causes the equilibrium density to be lower by up to 20%. Can be very time consuming, also check that n_vcsel_bgr is correct
#define ITERATE_QW // USE TO TURN ON QW ITERATION

//#define FFT_BPM_LENS // USE FFT_BPM when computing lens. This is required for proper lens functionality
//#define FFT_BPM_EXTRA_LAYERS // USE FFT_BPM in all normal (non-lens) material layers (SLOW)

 /*================================================
  * Global variables for both MAXWELL and SBE
  * */
  
const double t_max 		= 300.0*ns;
const double dt 		= 0.2*fs;				// Timestep
const double w0 		= 1.203740*e/hbar;			// Carrier frequency of light, 1030nm

const double out_DT_write 	= 100.0*ps;					// Print duration, when printing 
const double out_DT_freq  	= 0.5*fs;					// Write output every 1.0fs, when printing
const double out_DT_wait 	= 5.0*ns;					// Intervals between output strips
const double save_DT_global	= 10.0*ns;				// How often to save program state


//#define USE_CAVITY_SNAPSHOT // WARNING: LARGE SLOWDOWN.. (Up to a factor of 100, depending on # points and freq)
const double cav_snapshot_freq = 250.0*fs;			// How often to sample the cavity

const int SBE_COULOMB_THRESHOLDING 	= 1;					// Coulomb potential thresholding for SBE equations
const int SBE_FIELD_THRESHOLDING 	= 0;					// Field free thresholding for SBE equations
const int SBE_ANALYTIC_THRESHOLDING 	= 0;					// Analytic update thresholding for SBE equations

const double pulse_amplitude_threshold 		= 2.0e-2;			// Thresholding amplitude for SBE equations
const double pulse_amplitude_threshold_field 	= 1.0e-4;			// Thresholding for field free SBE equations
const double pulse_amplitude_threshold_analytic 	= 1.0e-8;			// Thresholding  amplitude for analytic SBE update


const int LOAD_STATE_NR 	= -1;							// If this is an int >=0 -> Will try to load that state
const int PL_MODE   		= 1;				// =1 to run Photoluminesence mode

const int OPTIMIZATION_MODE 	= 0;				// Run in optimization mode for parameters
const int ANALYSE_ENERGY_CURVE_ABS = 0;				// To find energy curve
const int ANALYSE_REFLECTION_SPEC  = 0;				// TO RUN AS REFLECTION WORKER

/*=================================================
 * Maxwell solver global variables
 * */

//Transverse VARIABLES
const double PUMP_RATIO      = 0.6; //w_0/w_p ratio
const double PUMP_SCALER     = PUMP_RATIO/0.6;
const double LENS_SPOT_GAIN  = 0.6*336*um;	
const double GAIN_SPOT_FWHM  = PUMP_SCALER*336*um; // number from Alex 2018
const double LENS_SPOT_SESAM = (1.0/5.0)*LENS_SPOT_GAIN; //Scale by a factor of  (should be 1/5 with current BPM)

const int NUM_TRANSVERSE_POINTS              = 1;//144; //256; // Should be even, faster with 2^N, fastest with machine multiple
const double TRANSVERSE_DOMAIN_SIZE_BG_RATIO = 0.8; // 75-79% gives supression of 1e-24 to 1e-10 for degree = 8
const double TRANSVERSE_DOMAIN_SIZE          = 5.0*GAIN_SPOT_FWHM;//5.0*GAIN_SPOT_FWHM; // Increase if noise is present in transverse profile.. // [-R/2,.. , 0,.. , (1-2/N)R]
const double TRANSVERSE_DOMAIN_OUTPUT_FWHM   = 0.8*TRANSVERSE_DOMAIN_SIZE; // FWHM of QWs where there should be output to file
const double TRANSVERSE_DOMAIN_OUTPUT_FREQ   = 11.668*um;//5.8334*um;//11.8*um; // Output at y=0, +/- f, +/- 2f,...


// Cavity setup
// Linear Cavity Setup cavity: Gain Chip -> Space1 -> Lens -> Space2 -> SESAM
const double real_cavity_length1    = 49053.8385107592694112*um;
const double real_cavity_length2    = 6622.1642764403850379*um;
const double real_cavity_lens_focus = 0.0061935615063968;

// Optical ABCD matrix: Gain -> SESAM
const double optical_abcd_A =  1.0-real_cavity_length1/real_cavity_lens_focus;
const double optical_abcd_B = real_cavity_length1+real_cavity_length2-real_cavity_length1*real_cavity_length2/real_cavity_lens_focus;
const double optical_abcd_C = -1.0/real_cavity_lens_focus;
const double optical_abcd_D =  1.0-real_cavity_length2/real_cavity_lens_focus;

//VCAV Setup cavity: OC->focus->length1->gain->length2->SESAM
//double VCAV_length2 = 10861.754972174*um; //Found in notes for F=10,
//const double VCAV_length2 = 12334.2089339*um;	//Computed for focus of 10 with w_0=0.6*336um
//const double VCAV_length2 = 45444.7898675878*um; //Smaller focus of 2.5 with w_0=0.6*336um
//const double VCAV_length2 = 17522.81432064076*um; //Smaller focus of 5 with w_0=0.6*336um (wrong)
//const double VCAV_length2 = 0.037189039199733; //FWHM intensity focus of 10 for w_sesam=(0.6*336/sqrt(10))um;
//const double VCAV_length2 = 0.024291778683655; //FWHM focus of 5 for w_sesam=40.32um;
//const double VCAV_length2 = 0.009917077119929; //FWHM focus of sqrt(5) for w_sesam=40.32um (paraxial propagator)
const double VCAV_length2 = 0.008768841103064; //FWHM focus of sqrt(5) for w_sesam=40.32um (modified for non-paraxial propagator)

const double VCAV_length1 = 100*um;
const double VCAV_focus   = (VCAV_length1+VCAV_length2)*1.04; 
const int VCAV_twoArmDeviceOutputLevel = 0;
const int VCAV_deviceOutputLevel = 0;
 
// DBR data
//const int num_dbr_layers = 61;				// Number of individual layers in DBR
const int num_dbr_layers = 61;				// Number of individual layers in DBR
//const double n_dbr_1 = 2.1;				// DBR Ref ind 1 [Ta2O5 @ 1040nm]
//const double n_dbr_2 = 3.63;				// DBR Ref ind 1 [AlGaAs @ 1040nm]
// Options...
//const double n_dbr_1 = 1.4499;				// DBR Ref ind 2 [SiO2 @ 1040nm]
//const double n_dbr_2 = 2.4813;				// DBR Ref ind 1 [TiO2 @ 1040nm]
//const double n_dbr_2 = 2.0972;				// DBR Ref ind 1 [Ta2O5 @ 1040nm]

const double n_dbr_2 = 3.3697;				// DBR Ref ind 1 [AlGaAs @ 1040nm]
const double n_dbr_1 = 2.9351;				// DBR Ref ind 2 [AlAs @ 1040nm]

// VCSEL data (RPG)
const int num_qw = 10;								// Number of QW's in VCSEL
const double length_qw_dbr = 0.5;			// Distance between QW 1 and DBR (or boundary) in units of lambda [0,1]
const double length_qw_qw = 0.5;			// Distance between each QW. In units of lambda [0,1]
const double length_qw_total = length_qw_dbr + (num_qw-1)*length_qw_qw + 0.25;	// Total length of cavity containing QW's in units of lambda [0,inf)
const double n_vcsel_bgr = 3.4453;			// VCSEL background refractive index [AlGaAs]
const double n_vcsel_cap_layer = 3.1778; 		// cap layer [InGaAs]
//const double n_vcsel_ar_coating = n_vcsel_cap_layer/sqrt(n_vcsel_bgr);		//IDEAL
//const double n_vcsel_ar_coating = 1.791989;		//IDEAL
const double n_vcsel_ar_coating = 1.45;		//IDEAL
const double n_vcsel_ar_coating2 = 2.0781;		//IDEAL

// Medium data
const double length_cav = 1600.25;//1600.25; //3200.25; 			// Length of medium in units of lambda (0,inf)
const double n_cav = 1.0;					// Ref ind of medium
//const double n_cav = 0.9397;					// Ref ind of medium

// Absorber data
const int num_abs = 1;					// Number of ABS's
const double length_abs_total = 0.0;			// Total length of cavity containing ABS's in units of lambda [0,inf)
const double length_abs_dbr = 0.75;			// Distance between ABS 1 and DBR (or boundary) in units of lambda [0,1]
//const double length_abs_dbr = 0.625;			// Distance between ABS 1 and DBR (or boundary) in units of lambda [0,1]
const double length_abs_abs = 0.5;			// Distance between each ABS. In units of lambda [0,1]
const double n_abs_bgr = 3.5006;				// ABS background refractive index
const double n_abs_cap_layer = 0; 			// NO absorber CAP layer
const double n_abs_ar_coating = 1.988;			// Optimal =sqrt(3.5145)=1.8747



// Initial pulse
const int    STARTL = 0;						// Start pulse at left side propagating right
const int    STARTR = 1;						// Start pulse at right side propagating left
const double pulse_initial_delay 		 	= 5.0*ps;	// Delay of pulse
//const double pulse_initial_amplitude 	 	= 1.0e7;		// amplitude of initial pulse

// DBR data
/* OLD DATA

const int num_dbr_layers = 46;				// Number of individual layers in DBR
//const double n_dbr_1 = 2.1;				// DBR Ref ind 1 [Ta2O5 @ 1040nm]
//const double n_dbr_2 = 3.63;				// DBR Ref ind 1 [AlGaAs @ 1040nm]
// Options...
//const double n_dbr_1 = 1.4499;				// DBR Ref ind 2 [SiO2 @ 1040nm]
//const double n_dbr_2 = 2.4813;				// DBR Ref ind 1 [TiO2 @ 1040nm]
//const double n_dbr_2 = 2.0972;				// DBR Ref ind 1 [Ta2O5 @ 1040nm]

const double n_dbr_1 = 2.9459;				// DBR Ref ind 2 [AlAs @ 1040nm]
const double n_dbr_2 = 3.4345;				// DBR Ref ind 1 [AlGaAs @ 1040nm]


// VCSEL data (RPG)
const int num_qw = 10;								// Number of QW's in VCSEL
const double length_qw_dbr = 0.5;			// Distance between QW 1 and DBR (or boundary) in units of lambda [0,1]
const double length_qw_qw = 0.5;			// Distance between each QW. In units of lambda [0,1]
const double length_qw_total = length_qw_dbr + (10-1.0)*length_qw_qw + 0.25 + 0.25;	// Total length of cavity containing QW's in units of lambda [0,inf)
const double n_vcsel_bgr = 3.4746;			// VCSEL background refractive index [AlGaAs]
const double n_vcsel_cap_layer = 3.1954; 		// cap layer [InGaAs]
const double n_vcsel_ar_coating = 1.988;		// HfO2

// Medium data
const double length_cav = 3200.25; 			// Length of medium in units of lambda (0,inf)
const double n_cav = 1.0;					// Ref ind of medium
//const double n_cav = 0.9397;					// Ref ind of medium

// Absorber data
const int num_abs = 1;					// Number of ABS's
const double length_abs_total = 0.0;			// Total length of cavity containing ABS's in units of lambda [0,inf)
const double length_abs_dbr = 0.75;			// Distance between ABS 1 and DBR (or boundary) in units of lambda [0,1]
const double length_abs_abs = 0.5;			// Distance between each ABS. In units of lambda [0,1]
const double n_abs_bgr = 3.5006;				// ABS background refractive index
const double n_abs_cap_layer = 0; 			// NO absorber CAP layer
const double n_abs_ar_coating = 1.88;			// Optimal =sqrt(3.5145)=1.8747



// Initial pulse
const int    STARTL = 0;						// Start pulse at left side propagating right
const int    STARTR = 1;						// Start pulse at right side propagating left
const double pulse_initial_delay 		 	= 5.0*ps;	// Delay of pulse
//const double pulse_initial_amplitude 	 	= 1.0e7;		// amplitude of initial pulse
const double pulse_amplitude_threshold 	= 2.0e2;			// Thresholdning amplitude for p_k equations
*/

#endif


