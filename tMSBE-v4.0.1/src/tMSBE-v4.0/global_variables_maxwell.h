
/*======================================================================
 * This file defines the global variables for Maxwell solver that are stored
 * in setup_simulation_variables.h
 * 
 * */
 
 
#ifndef __GLOBAL_VARIABLES_MAXWELL_H_INCLUDED__
#define __GLOBAL_VARIABLES_MAXWELL_H_INCLUDED__

#include "setup_simulation_variables.h"

extern const double t_max;					// Maximum simulation time
extern const double dt;						// Timestep of simulation

extern const double w0;						// Frequency of incoming light;

extern const int LOAD_STATE_NR;							// If this is an int >=0 -> Will try to load that state

extern const double LENS_SPOT_GAIN; // number from Alex 2018
extern const double PUMP_RATIO; //w_0/w_p ratio
extern const double PUMP_SCALER;
extern const double LENS_SPOT_GAIN;
extern const double GAIN_SPOT_FWHM;	
extern const double LENS_SPOT_SESAM; //Scale by a factor of 4

extern const int    NUM_TRANSVERSE_POINTS; //256; // Should be even, faster with 2^N, fastest with machine multiple
extern const double TRANSVERSE_DOMAIN_SIZE_BG_RATIO; // 75-79% gives supression of 1e-24 to 1e-10 for degree = 8
extern const double TRANSVERSE_DOMAIN_SIZE; // Increase if noise is present in transverse profile.. // [-R/2,.. , 0,.. , (1-2/N)R]
extern const double TRANSVERSE_DOMAIN_OUTPUT_FWHM; // FWHM of QWs where there should be output to file
extern const double TRANSVERSE_DOMAIN_OUTPUT_FREQ; // Output at y=0, +/- f, +/- 2f,...


// Cavity setup
// Linear Cavity Setup cavity: Gain Chip -> Space1 -> Lens -> Space2 -> SESAM
extern const double real_cavity_length1;
extern const double real_cavity_length2;
extern const double real_cavity_lens_focus;

// Optical ABCD matrix: Gain -> SESAM
extern const double optical_abcd_A;
extern const double optical_abcd_B;
extern const double optical_abcd_C;
extern const double optical_abcd_D;

//VCAV Setup cavity: OC->focus->length1->gain->length2->SESAM
extern const double VCAV_length2; 
extern const double VCAV_length1;
extern const double VCAV_focus;
extern const int VCAV_twoArmDeviceOutputLevel;
extern const int VCAV_deviceOutputLevel;

extern const int num_dbr_layers;
extern const double n_dbr_1;	// DBR refractive indices
extern const double n_dbr_2;	// DBR refractive indices

extern const int num_qw;
extern const double length_qw_total;			// No extra space in ABSORBER
extern const double length_qw_dbr;
extern const double length_qw_qw;
extern const double n_vcsel_bgr;	// VCSEL background refractive index
extern const double n_vcsel_cap_layer; 	// nInGaP cap layer
extern const double n_vcsel_ar_coating;	// ar coating cap
extern const double n_vcsel_ar_coating2;// ar coating cap

extern const double length_cav;
extern const double n_cav;

extern const int num_abs;						// Number of ABS's
extern const double length_abs_total;
extern const double length_abs_dbr;	// Distance between first ABS and DBR (or boundary) in units of lambda [0,1]
extern const double length_abs_abs;			// Distance between each ABS. In units of lambda [0,1]
extern const double n_abs_bgr;				// ABS background refractive index
extern const double n_abs_cap_layer; 			// NO absorber cap layer
extern const double n_abs_ar_coating;	// Perfect AR coating on ABS		

extern const int maxwell_num_saved_states; // Number of previous states to store in E+ and E-


#endif
