
/*
	The VECSEL class is used to design the device
	Can contain a series of devices from other classes
	Solvers for MAXWELL and for each device
*/

#include "classVECSEL.h"
#include "constantsAndMiscUnits.h"
#include "fileIO.cpp"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <unistd.h>

// For rounding
#include <string>
#include <math.h>
#include <iomanip>

#ifdef USE_THREADS
	#include <vector>
	#include <thread>
	#include <algorithm>
#endif

#include <mpi.h>

#ifdef USE_OPENMP
	#include <omp.h>
#endif

using namespace std;

/*
	Construction function
*/

VECSEL::VECSEL()
{
	cout << "Creating empty VECSEL" << endl;
	setName("New VECSEL");
	setToFileOutputKey("out_");
	setLambda(0.0);
	setNumberCavities(0);
	setNumberTwoArmInterfaces(0);
	setNumberTwoArmCavities(0);
	setNumberDevices(0);
	setNumberTwoArmDevices(0);
	setNumberBoundaries(0);

	VECSEL_transverse_points_number = -1;
	VECSEL_transverse_points_device_number = -1;
	VECSEL_transverse_points_R_max = 0.0;
	VECSEL_transverse_points_boundary_guard_ratio = 0.78;
	VECSEL_transverse_points_y = NULL;
	VECSEL_transverse_points_device_y = NULL;
	
	VECSEL_pulse_start_l = 0;
	VECSEL_pulse_start_r = 0;
	
	VECSEL_QW_FEEDBACK = 1.0;  
	
	VECSEL_initial_energy_shift = 0.0;
	VECSEL_cav_snapshot_num_points = 0;
	VECSEL_cav_snapshot_E = NULL;
	VECSEL_cav_snapshot_E_re = NULL;
	VECSEL_cav_snapshot_E_im = NULL;

	VECSEL_initial_transverse_FWHM = 500.0*um;
	VECSEL_initial_transverse_pulse_profile = NULL;
	VECSEL_initial_temp_profile_T = NULL;

	set_transverse_QW_pump_profile_SuperGaussian(13,2000.0*um); // Set initial SG pump profile
	//set_transverse_QW_temp_profile_SuperGaussian(13,2000.0*um); // Set initial SG temp profile

	cavity_trans_E_pl = NULL;
	cavity_trans_E_mi = NULL;
	cavity_trans_E_pl_k = NULL;
	cavity_trans_E_mi_k = NULL;
	cavity_trans_MacPol = NULL;
	
	// MPI related arrays
	MPI_MY_RANK = 0;
	MPI_WORK_DIST = NULL; // size(MPI_WORK_DIST) = number_of_workers x 2. Contains start and stop index of work
	MPI_WORK_DIST_E_OFFSET = NULL;
	MPI_WORK_DIST_P_OFFSET = NULL;
	MPI_WORK_DIST_E_SIZE = NULL;
	MPI_WORK_DIST_P_SIZE = NULL;
	MPI_WORK_DEVICE_SIZE = NULL;
	MPI_WORK_DEVICE_OFFSET = NULL;
	MPI_WORK_DIST_TOTAL = -1;
	MPI_FLAG;
	MPI_WORK_DIST_E_GLOBAL = NULL;
	MPI_WORK_DIST_P_GLOBAL = NULL;
	MPI_WORK_DIST_E_LOCAL = NULL;
	MPI_WORK_DIST_P_LOCAL = NULL;
	MPI_WORK_GANG = NULL;
	MPI_WORK_GANG_SIZE = -1;
	MPI_LoadBalancer = NULL;
	MPI_LoadBalancer_index_set = NULL;
	MPI_LoadBalancer_P_tmp = NULL;
	MPI_LoadBalancer_E_tmp = NULL;
	#ifdef MPI_BALANCE_WORKLOAD
	MPI_load = NULL; // timer
	#endif

	test_VECSEL_iteration = 0;
	init_VECSEL_iteration = 0;

	filter_diagnostics_prev_E = 0.0;
}

void VECSEL::Print() const 
{
	cout << "Print VECSEL:" << endl;
	cout << " name      = " << getName() << endl;
	cout << " filename  = " << getToFileOutputKey() << endl;
	cout << " lambda    = " << getLambda() << endl;
	cout << " # modules    = " << getNumberModules() << endl;
	cout << " # Boundaries = " << getNumberBoundaries() << endl;
	cout << " # Devices    = " << getNumberDevices() << endl;
	cout << " # Cavities   = " << getNumberCavities() << endl;
	cout << " # TwoArmDevices    = " << getNumberTwoArmDevices() << endl;
	cout << " # TwoArmCavities   = " << getNumberTwoArmCavities() << endl;
	cout << " # TwoArmInterfaces   = " << getNumberTwoArmInterfaces() << endl;
	if (VECSEL_transverse_points_y!=NULL)
	{
		cout << " Transverse dimension = ";
		for(int i = 0; i < VECSEL_transverse_points_number; i++)
		{
			cout << VECSEL_transverse_points_y[i]/um << " ";
		}

		cout << " [um]" << endl;
	}
	if (VECSEL_transverse_points_device_y != NULL)
	{
		cout << "Transverse device output points = ";
		for(int i = 0; i < VECSEL_transverse_points_device_number; i++)
		{
			cout << VECSEL_transverse_points_device_y[i]/um << " ";
		}

		cout << " [um]" << endl;

	}

	cout << " Gain chip SG pump profile degree = " << VECSEL_initial_pump_profile_SG_degree << endl;
	cout << " Gain chip SG pump profile FWHM = " << VECSEL_initial_pump_profile_SG_FWHM/um << " [um]" << endl;
	cout << " Gain chip SG temp profile degree = " << VECSEL_initial_temp_profile_SG_degree << endl;
	cout << " Gain chip SG temp profile FWHM = " << VECSEL_initial_temp_profile_SG_FWHM/um << " [um]" << endl;
	
	for(int i = 0; i < getNumberModules(); i++)
	{
		if (modules[i].isDevice()||modules[i].isTwoArmDevice())
		{
			cout << "<------- DEVICE_PRINT() ------------" << endl;
		} else {
			modules[i].Print();
		}
	}
	
	cout << "Cavity quick index:" << endl;
	for(unsigned i = 0; i < quick_index_cavity.size(); i++)
	{
		cout << quick_index_cavity[i] << " ";
	}
	cout << endl;

	cout << "Cavity quick index no BPM:" << endl;
	for(unsigned i = 0; i < quick_index_cavity_noBPM.size(); i++)
	{
		cout << quick_index_cavity_noBPM[i] << " ";
	}
	cout << endl;

	cout << "Cavity quick index free space:" << endl;
	for(unsigned i = 0; i < quick_index_cavity_freeSpace.size(); i++)
	{
		cout << quick_index_cavity_freeSpace[i] << " ";
	}
	cout << endl;

	cout << "Cavity quick index lens:" << endl;
	for(unsigned i = 0; i < quick_index_cavity_lens.size(); i++)
	{
		cout << quick_index_cavity_lens[i] << " ";
	}
	cout << endl;

	cout << "Cavity quick index halfCav lens:" << endl;
	for(unsigned i = 0; i < quick_index_cavity_lens_halfCav.size(); i++)
	{
		cout << quick_index_cavity_lens_halfCav[i] << " ";
	}
	cout << endl;

	cout << "Cavity quick index noQW:" << endl;
	for(unsigned i = 0; i < quick_index_cavity_noQW.size(); i++)
	{
		cout << quick_index_cavity_noQW[i] << " ";
	}
	cout << endl;
	cout << "Cavity quick index w/QW:" << endl;
	for(unsigned i = 0; i < quick_index_cavity_QW.size(); i++)
	{
		cout << quick_index_cavity_QW[i] << " ";
	}
	cout << endl;
	
	cout << "Boundary quick index:" << endl;
	for(unsigned i = 0; i < quick_index_boundary.size(); i++)
	{
		cout << quick_index_boundary[i] << " ";
	}
	cout << endl;
	
	cout << "Total Device quick index:" << endl;
	for(unsigned i = 0; i < quick_index_totalDevice.size(); i++)
	{
		cout << quick_index_totalDevice[i] << " ";
	}
	cout << endl;
	
	cout << "Device quick index:" << endl;
	for(unsigned i = 0; i < quick_index_device.size(); i++)
	{
		cout << quick_index_device[i] << " ";
	}
	cout << endl;
	
	cout << "TwoArmCavity quick index:" << endl;
	for(unsigned i = 0; i < quick_index_twoArmCavity.size(); i++)
	{
		cout << quick_index_twoArmCavity[i] << " ";
	}
	cout << endl;
	
	cout << "TwoArmCavity quick index w/QW:" << endl;
	for(unsigned i = 0; i < quick_index_twoArmCavity_QW.size(); i++)
	{
		cout << quick_index_twoArmCavity_QW[i] << " ";
	}
	cout << endl;
	
	cout << "TwoArmCavity quick index noQW:" << endl;
	for(unsigned i = 0; i < quick_index_twoArmCavity_noQW.size(); i++)
	{
		cout << quick_index_twoArmCavity_noQW[i] << " ";
	}
	cout << endl;
	
	cout << "TwoArmCavity quik index birefringent crystals:" << endl;
	for(unsigned i = 0; i < quick_index_birefringentCrystal.size(); i++)
	{
		cout << quick_index_birefringentCrystal[i] << " ";
	}
	cout << endl;
	
	cout << "TwoArmDevice quick index:" << endl;
	for(unsigned i = 0; i < quick_index_twoArmDevice.size(); i++)
	{
		cout << quick_index_twoArmDevice[i] << " ";
	}
	cout << endl;
	
	cout << "TwoArmInterface quick index:" << endl;
	for(unsigned i = 0; i < quick_index_twoArmInterface.size(); i++)
	{
		cout << quick_index_twoArmInterface[i] << " ";
	}
	cout << endl;
	
	cout << "TwoArmPostCav quick index:" << endl;
	for(unsigned i = 0; i < quick_index_twoArmPostCav.size(); i++)
	{
		cout << quick_index_twoArmPostCav[i] << " ";
	}
	cout << endl;
	
	cout << "Device prev cavity quick index:" << endl;
	for(unsigned i = 0; i < quick_index_totalDevice.size(); i++)
	{
		cout << quick_index_device_previous_cavity[i] << " ";
	}
	cout << endl;
	
}


/*
	n1 / n2 			-> Refractive indices in medium
	nNext				-> Refractive indices in the next layer AFTER the DBR (To ensure ordering)
	numLayers 			-> Number of layers in the DBR
	angle_of_incidence 	-> Angle in radians for TE
	external_index  	-> Index in incoming medium
*/
void VECSEL::addDBR_LEFT(double n1, double n2, int numLayers, double nNext, double angle_of_incidence, double external_index)
{
	double N_SUBSTRATE = n1;
	if (n2 > n1)
	{
		N_SUBSTRATE = n2;
	}
	
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 	= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	n1 		= n1*sqrt(1.0-a_sin2/(n1*n1));
	n2 		= n2*sqrt(1.0-a_sin2/(n2*n2));
	nNext 	= nNext*sqrt(1.0-a_sin2/(nNext*nNext));
	
	if (n1 == n2)
	{
		cout << "addDBR_LEFT: Refractive indices are equal, not funny..." << endl;
		exit(-1);
	} else if ((n1 <= 0.0)||(n2 <= 0.0)||(nNext <= 0.0))
	{
		cout << "addDBR_LEFT: Refractive indices have to be strictly positive" << endl;
		exit(-1);
	}
	
	if (numLayers < 0)
	{
		cout << "addDBR_LEFT: Have to use a positive number of layers" << endl;
		exit(-1);
	}

	// Get starting position
	int startPos = getNumberModules();
	if (startPos == 0)
	{
		// Add boundary to first element
		if (numLayers == 0)
		{
			cout << "addDBR_LEFT: Adding REFLECTING boundary on the left" << endl;
			addBoundary(1.0,N_SUBSTRATE);
		} else {
			cout << "addDBR_LEFT: Adding ABSORBING boundary on the left" << endl;
			addBoundary(0.0,N_SUBSTRATE);
		}
		
	} else if (startPos == 1)
	{
		// Starting on the left
		// Check if previous was Boundary
		if (!modules[0].isBoundary())
		{
			cout << "addDBR_LEFT: 1st module has to be a BOUNDARY..." << endl;
			exit(-1);
		}
	}
	
	// Ensure alternating reflectors
	double n_bragg[2] = {0,0};
	if (numLayers % 2 == 0)
	{
		// Pair number of layers
		if (n2 == nNext)
		{
			// Reverse ordering
			n_bragg[0] = n2;
			n_bragg[1] = n1;
		} else {
			n_bragg[0] = n1;
			n_bragg[1] = n2;
		}
	} else {
	
		// Odd number of layers
		if (n1 == nNext)
		{
			n_bragg[0] = n2;
			n_bragg[1] = n1;
		} else {
			n_bragg[0] = n1;
			n_bragg[1] = n2;
		}
	}
	
	Module *tmp;
	// Create numLayers cavities of length lambda/4/n_bragg[j]
	double deviceWidth = 0;
	for(int i = 0; i < numLayers; i++)
	{
		deviceWidth = getLambda()/(4.0*n_bragg[i%2]);
		tmp = addCavity(deviceWidth, n_bragg[i%2], angle_of_incidence, external_index);
	}
}

void VECSEL::addDBR_LEFT_nL(double n1, double L1, double n2, double L2, int numLayers, double angle_of_incidence, double external_index)
{
	double N_SUBSTRATE = n1;
	if (n2 > n1)
	{
		N_SUBSTRATE = n2;
	}
	
	if (n1 == n2)
	{
		cout << "addDBR_LEFT: Refractive indices are equal, not funny..." << endl;
		exit(-1);
	} else if ((n1 <= 0.0)||(n2 <= 0.0))
	{
		cout << "addDBR_LEFT: Refractive indices have to be strictly positive" << endl;
		exit(-1);
	}
	
	if (numLayers < 0)
	{
		cout << "addDBR_LEFT: Have to use a positive number of layers" << endl;
		exit(-1);
	}

	// Get starting position
	int startPos = getNumberModules();
	if (startPos == 0)
	{
		// Add boundary to first element
		if (numLayers == 0)
		{
			cout << "addDBR_LEFT: Adding REFLECTING boundary on the left" << endl;
			addBoundary(1.0,N_SUBSTRATE);
		} else {
			cout << "addDBR_LEFT: Adding ABSORBING boundary on the left" << endl;
			addBoundary(0.0,N_SUBSTRATE);
		}
		
	} else if (startPos == 1)
	{
		// Starting on the left
		// Check if previous was Boundary
		if (!modules[0].isBoundary())
		{
			cout << "addDBR_LEFT: 1st module has to be a BOUNDARY..." << endl;
			exit(-1);
		}
	}
	
	// Ensure alternating reflectors
	double n_bragg[2] = {n1,n2};
	double L_bragg[2] = {L1,L2};
	
	Module *tmp;
	// Create numLayers cavities of length lambda/4/n_bragg[j]
	double deviceWidth = 0;
	for(int i = 0; i < numLayers; i++)
	{
		tmp = addCavity(L_bragg[i%2], n_bragg[i%2], angle_of_incidence, external_index);
	}
}
// Based on the design idea by: Zhang. et al in Opt Quant Electron 2015 47:423-431
// Equations can be found in Calvez 2002 Photonics tech lett vol 14
// Make repetitions of ((n1n2)^D n1)^N where D and lambda0 are both solutions of an equation

void VECSEL::addDBM_LEFT(double n1, double n2, double nNext, int D, int N, double lambda0, double angle_of_incidence, double external_index)
{
	double N_SUBSTRATE = n1;
	if (n2 > n1)
	{
		N_SUBSTRATE = n2;
	}
	
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 	= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	n1 		= n1*sqrt(1.0-a_sin2/(n1*n1));
	n2 		= n2*sqrt(1.0-a_sin2/(n2*n2));
	
	if (n1 == n2)
	{
		cout << "addDBR_LEFT: Refractive indices are equal, not funny..." << endl;
		exit(-1);
	} else if ((n1 <= 0.0)||(n2 <= 0.0))
	{
		cout << "addDBR_LEFT: Refractive indices have to be strictly positive" << endl;
		exit(-1);
	}
	
	if ((N <= 0)||(D<=0))
	{
		cout << "addDBR_LEFT: Have to use a positive number of layers" << endl;
		exit(-1);
	}

	// Get starting position
	int startPos = getNumberModules();
	if (startPos == 0)
	{
		// Add boundary to first element		
		cout << "addDBR_LEFT: Adding ABSORBING boundary on the left" << endl;
		addBoundary(0.0,N_SUBSTRATE);
	
		
	} else if (startPos == 1)
	{
		// Starting on the left
		// Check if previous was Boundary
		if (!modules[0].isBoundary())
		{
			cout << "addDBR_LEFT: 1st module has to be a BOUNDARY..." << endl;
			exit(-1);
		}
	}
	
	// Ensure alternating reflectors
	double n_bragg[2] = {0,0};
	
	if (n1 == nNext)
	{
		n_bragg[0] = n2;
		n_bragg[1] = n1;
	} else {
		n_bragg[0] = n1;
		n_bragg[1] = n2;
	}
	
	
	Module *tmp;
	// Create numLayers cavities of length lambda/4/n_bragg[j]
	double deviceWidth = 0;
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < D; j++)
		{
			deviceWidth = lambda0/(4.0*n_bragg[0]);
			tmp = addCavity(deviceWidth, n_bragg[0]);
			
			deviceWidth = lambda0/(4.0*n_bragg[1]);
			tmp = addCavity(deviceWidth, n_bragg[1]);
		}
		
		deviceWidth = lambda0/(4.0*n_bragg[0]);
		tmp = addCavity(deviceWidth, n_bragg[0]);
	}
}


/*
	lambda  -> Target lambda for DBR
	n1 / n2 -> Refractive indices in medium
	nNext	-> Refractive indices in the next layer AFTER the DBR (To ensure ordering)
	numLayers -> Number of layers in the DBR
*/
void VECSEL::addDBR_LEFT_LAMBDA(double lambda, double n1, double n2, int numLayers, double nNext)
{
	double N_SUBSTRATE = n1;
	if (n2 > n1)
	{
		N_SUBSTRATE = n2;
	}
	
	if (n1 == n2)
	{
		cout << "addDBR_LEFT: Refractive indices are equal, not funny..." << endl;
		exit(-1);
	} else if ((n1 <= 0.0)||(n2 <= 0.0)||(nNext <= 0.0))
	{
		cout << "addDBR_LEFT: Refractive indices have to be strictly positive" << endl;
		exit(-1);
	}
	
	if (numLayers < 0)
	{
		cout << "addDBR_LEFT: Have to use a positive number of layers" << endl;
		exit(-1);
	}

	// Get starting position
	int startPos = getNumberModules();
	if (startPos == 0)
	{
		// Add boundary to first element
		if (numLayers == 0)
		{
			cout << "addDBR_LEFT: Adding REFLECTING boundary on the left" << endl;
			addBoundary(1.0,N_SUBSTRATE);
		} else {
			cout << "addDBR_LEFT: Adding ABSORBING boundary on the left" << endl;
			addBoundary(0.0,N_SUBSTRATE);
		}
		
	} else if (startPos == 1)
	{
		// Starting on the left
		// Check if previous was Boundary
		if (!modules[0].isBoundary())
		{
			cout << "addDBR_LEFT: 1st module has to be a BOUNDARY..." << endl;
			exit(-1);
		}
	}
	
	// Ensure alternating reflectors
	double n_bragg[2] = {0,0};
	if (numLayers % 2 == 0)
	{
		// Pair number of layers
		if (n2 == nNext)
		{
			// Reverse ordering
			n_bragg[0] = n2;
			n_bragg[1] = n1;
		} else {
			n_bragg[0] = n1;
			n_bragg[1] = n2;
		}
	} else {
	
		// Odd number of layers
		if (n1 == nNext)
		{
			n_bragg[0] = n2;
			n_bragg[1] = n1;
		} else {
			n_bragg[0] = n1;
			n_bragg[1] = n2;
		}
	}
	
	Module *tmp;
	// Create numLayers cavities of length lambda/4/n_bragg[j]
	double deviceWidth = 0;
	for(int i = 0; i < numLayers; i++)
	{
		deviceWidth = lambda/(4.0*n_bragg[i%2]);
		tmp = addCavity(deviceWidth, n_bragg[i%2]);
	}
}

/*
	Add in a DBR design for lambda_pump, that has a total length of lambda/4.
	
	The order of the layers is determined by nNext and n1,n2.
	 - default order is with: n1, n2, n1, ..., n2, n1, with numLayers
	   unless n1 == nNext or n2 == nNext, then the final layer is different
	tune - determines the amount of extra layer on the left / right.
	 - tune = 0 => all on the left
	 - tune = 1 => all on the right
*/
void VECSEL::addDBR_LEFT_PUMP(double lambda, double lambda_pump, double n1, double n2, int numLayers, double nTune, double nNext, double tune)
{
	double N_SUBSTRATE = n1;
	if (n2 > n1)
	{
		N_SUBSTRATE = n2;
	}
	
	if ((tune < 0.0)||(tune > 1.0))
	{
		cout << "addDBR_LEFT_PUMP: Tune has to be in the rage [0,1]" << endl;
		exit(-1);
	}
	if (n1 == n2)
	{
		cout << "addDBR_LEFT_PUMP: Refractive indices are equal, not funny..." << endl;
		exit(-1);
	} else if ((n1 <= 0.0)||(n2 <= 0.0)||(nNext <= 0.0))
	{
		cout << "addDBR_LEFT_PUMP: Refractive indices have to be strictly positive" << endl;
		exit(-1);
	}
	
	if (numLayers < 0)
	{
		cout << "addDBR_LEFT_PUMP: Have to use a positive number of layers" << endl;
		exit(-1);
	}

	// Get starting position
	int startPos = getNumberModules();
	if (startPos == 0)
	{
		// Add boundary to first element
		if (numLayers == 0)
		{
			cout << "addDBR_LEFT_PUMP: Adding REFLECTING boundary on the left" << endl;
			addBoundary(1.0,N_SUBSTRATE);
		} else {
			cout << "addDBR_LEFT_PUMP: Adding ABSORBING boundary on the left" << endl;
			addBoundary(0.0,N_SUBSTRATE);
		}
		
	} else if (startPos == 1)
	{
		// Starting on the left
		// Check if previous was Boundary
		if (!modules[0].isBoundary())
		{
			cout << "addDBR_LEFT_PUMP: 1st module has to be a BOUNDARY..." << endl;
			exit(-1);
		}
	}
	
	// Ensure alternating reflectors
	double n_bragg[2] = {0,0};
	if (numLayers % 2 == 0)
	{
		// Pair number of layers
		if (n2 == nNext)
		{
			// Reverse ordering
			n_bragg[0] = n2;
			n_bragg[1] = n1;
		} else {
			n_bragg[0] = n1;
			n_bragg[1] = n2;
		}
	} else {
	
		// Odd number of layers
		if (n1 == nNext)
		{
			n_bragg[0] = n2;
			n_bragg[1] = n1;
		} else {
			n_bragg[0] = n1;
			n_bragg[1] = n2;
		}
	}
	
	//======================================
	// Adding in extra layer  for stability
	//======================================
	double total_pump_length = lambda/4.0;
	double MIN_LAYER_SIZE = 20*nm;
	double diff_length = 0;
	for(unsigned i = 0; i < numLayers; i++)
	{
		diff_length += lambda_pump/(4.0*n_bragg[i%2]);
	}
	
	int mult = floor(diff_length/total_pump_length);
	
	if (mult == 0)
	{
		cout << "addDBR_LEFT_PUMP: Number of layers is too few for lambda/4 " << endl;
	}
	
	double res  = diff_length - ((double)mult)*total_pump_length; // Residual
	
	
	Module *tmp;
	
	// Add in residual layer
	if (tune < 1)
	{
		if (MIN_LAYER_SIZE > res*(1.0-tune))
		{
			cout << "addDBR_LEFT_PUMP: Extra cavity (first) is smaller than LIMIT of " << MIN_LAYER_SIZE/nm << " [nm]" << endl;
			exit(-1);
		}
		tmp = addCavity(res*(1.0-tune), nTune);
		cout << "addDBR_LEFT_PUMP: Extra cavity (first) added of length = " << res*(1-tune)/nm << " [nm]" << endl;
	}
	
	// Create numLayers cavities of length lambda/4/n_bragg[j]
	double deviceWidth = 0;
	for(int i = 0; i < numLayers; i++)
	{
		deviceWidth = lambda_pump/(4.0*n_bragg[i%2]);
		tmp = addCavity(deviceWidth, n_bragg[i%2]);
	}
	
	// Add last residual layer
	if (tune > 0)
	{
		if (MIN_LAYER_SIZE > res*tune)
		{
			cout << "addDBR_LEFT_PUMP: Extra cavity (last ) is smaller than LIMIT of " << MIN_LAYER_SIZE/nm << " [nm]" << endl;
			exit(-1);
		}
		tmp = addCavity(res*tune, nTune);
		cout << "addDBR_LEFT_PUMP: Extra cavity (last ) added of length = " << res*tune/nm << " [nm]" << endl;
	}
}


/*
	n1 / n2 -> Refractive indices in medium
	nPrev	-> Refractive indices in the previous layer BEFORE the DBR (To ensure ordering)
	numLayers -> Number of layers in the DBR
*/

void VECSEL::addDBR_RIGHT(double n1, double n2, int numLayers, double nPrev)
{
	if (n1 == n2)
	{
		cout << "addDBR_RIGHT: Refractive indices are equal, not funny..." << endl;
		exit(-1);
	} else if ((n1 <= 0.0)||(n2 <= 0.0)||(nPrev <= 0.0))
	{
		cout << "addDBR_RIGHT: Refractive indices have to be strictly positive" << endl;
		exit(-1);
	}
	
	if (numLayers < 0)
	{
		cout << "addDBR_RIGHT: Have to use a positive number of layers" << endl;
		exit(-1);
	}
	
	int startInd = getNumberModules();
	if ((startInd == 0)||(startInd == 1))
	{
		// Should use addDBR_LEFT
		cout << "addDBR_RIGHT: Use addDBR_RIGHT when starting on the left?" << endl;
		exit(-1);
	}
	
	// Ensure alternating reflectors
	double n_bragg[2] = {0,0};
	if (n1 == nPrev)
	{
	//	cout << "1: " << n1 << " == " << nPrev << endl;
	//	cout << "First is = " << n2 << endl;
		// Reverse ordering
		n_bragg[0] = n2;
		n_bragg[1] = n1;
	} else {
	//	cout << "2: " << n1 << " != " << nPrev << endl;
	//	cout << "First is = " << n1 << endl;
		n_bragg[0] = n1;
		n_bragg[1] = n2;
	}
	
	Module *tmp;
	// Create numLayers cavities of length lambda/4/n_bragg[j]
	double deviceWidth = 0;
	for(int i = 0; i < numLayers; i++)
	{
		deviceWidth = getLambda()/(4.0*n_bragg[i%2]);
		tmp = addCavity(deviceWidth, n_bragg[i%2]);
	}
	
}


/*
	n1 / n2 			-> Refractive indices in medium
	nNext				-> Refractive indices in the next layer AFTER the DBR (To ensure ordering)
	numLayers 			-> Number of layers in the DBR
	angle_of_incidence 	-> Angle in radians for TE
	external_index  	-> Index in incoming medium
*/
void VECSEL::addTwoArmDBR_FRONT(double n1, double n2, int numLayers, double nBack, double angle_of_incidence, double external_index)
{
	double N_SUBSTRATE = n1;
	if (n2 > n1)
	{
		N_SUBSTRATE = n2;
	}
	
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 		= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	double intIndex 	= nBack;
	n1 			= n1*sqrt(1.0-a_sin2/(n1*n1));
	n2 			= n2*sqrt(1.0-a_sin2/(n2*n2));
	nBack 			= nBack*sqrt(1.0-a_sin2/(nBack*nBack));
	
	if (n1 == n2)
	{
		cout << "addTwoArmDBR_FRONT: Refractive indices are equal, not funny...but letting it slide " << endl;
		//exit(-1);
	} else if ((n1 <= 0.0)||(n2 <= 0.0)||(nBack <= 0.0))
	{
		cout << "addTwoArmDBR_FRONT: Refractive indices have to be strictly positive" << endl;
		exit(-1);
	}
	
	if (numLayers < 0)
	{
		cout << "addTwoArmDBR_FRONT: Have to use a positive number of layers" << endl;
		exit(-1);
	}

	// Get starting position
	int startPos = getNumberModules();
	Module *tmp;
		
	if (startPos == 0)
	{
		// Add boundary to first element
		if (numLayers == 0)
		{
			cout << "addTwoArmDBR_FRONT: Adding REFLECTING boundary on the left" << endl;
			addBoundary(1.0,N_SUBSTRATE);
		} else {
			cout << "addTwoArmDBR_FRONT: Adding ABSORBING boundary on the left" << endl;
			addBoundary(0.0,N_SUBSTRATE);
		}
	}
	
	if (modules[startPos-1].isCavity())
	{
		tmp = addTwoArmInterface(1.0, angle_of_incidence, 1.0, intIndex);	
	}

	// Ensure alternating reflectors
	double n_bragg[2] = {0,0};
	if (numLayers % 2 == 0)
	{
		// Even number of layers
		if (n2 == nBack)
		{
			n_bragg[0] = n1;
			n_bragg[1] = n2;
		} else {
			n_bragg[0] = n2;
			n_bragg[1] = n1;
		}
	} else {
	
		// Odd number of layers
		if (n1 == nBack)
		{
			n_bragg[0] = n2;
			n_bragg[1] = n1;
		} else {
			n_bragg[0] = n1;
			n_bragg[1] = n2;
		}
	}
	
	// Create numLayers cavities of length lambda/4/n_bragg[j]
	double deviceWidth = 0;
	for(int i = 0; i < numLayers; i++)
	{
		deviceWidth = getLambda()/(4.0*n_bragg[i%2]);
		tmp = addTwoArmCavity(deviceWidth, n_bragg[i%2], angle_of_incidence, external_index);
	}
}


/*
 *  Create a gain medium on the left with QW. 
 * 
 * 	|---o-o-o-o-o-------------|--CAP--|--AR--|
 * 
 *   . If there is nothing on the left of this, a boundary is added
 * 
 *  Where:
 *  numQW			-> Number of QW's in the medium [0,inf)
 * 						- If numQW = 0		=> Empty medium of length cavityLength
 * 	dx0_qw			-> Distance from LEFTMOST QW's to the edge in units of WHOLE PEAKS
 *  dx1_qw			-> Distance from RIGHTMOST QW's to right edge in units of WHOLE PEAKS
 * 	cavityIndex		-> Refractive background index in medium
 * 	capIndex		-> Refractive index of CAP layer, if capIndex = 0 then NO CAP layer
 *  arIndex			-> Refractive index of AR coating, if arIndex = 0 then NO ar coating
 *  fillerLength	-> Length of cavity to fill in AFTER AR and CAP layer in units of m
 *  angle_of_incidence -> Angle in radians for TE
 *  external_index  -> Index in incoming medium
 * */

void VECSEL::addCUSTOM_CLUSTER_QW_LEFT(int numClusters, int *clusterQW, int dx0_qw, double dx1_qw, double cavityIndex, double capIndex, double arIndex, double angle_of_incidence, double external_index)
{
	// Distance between QW that are clustered
	//double BARRIER_LENGTH = 16.35*nm;
	//double BARRIER_LENGTH = 13.0*nm;
	double BARRIER_LENGTH = 10.0*nm;
	
	
	
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 	= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	cavityIndex 	= cavityIndex*sqrt(1.0-a_sin2/(cavityIndex*cavityIndex));
	capIndex 		= capIndex*sqrt(1.0-a_sin2/(capIndex*capIndex));
	arIndex 		= arIndex*sqrt(1.0-a_sin2/(arIndex*arIndex));
	
	
	int numQW = 0;
	for(unsigned i = 0; i < numClusters; i++)
	{
		numQW += clusterQW[i];
		if (clusterQW[i] < 0)
		{
			cout << "addCUSTOM_CLUSTER_QW_LEFT(): Cannot have NEGATIVE number of QW's" << endl;
			exit(-1);
		}
		if ((clusterQW[i]-1)*BARRIER_LENGTH + BARRIER_LENGTH > getLambda()*0.5/cavityIndex)
		{
			cout << "addCUSTOM_CLUSTER_QW_LEFT(): Too many QW's in a node" << endl;
			cout << "Num QW's = " << clusterQW[i] << endl;
			exit(-1);
		}
	}
	
	if (numClusters > 0)
	{
		if ((clusterQW[0] == 0)||(clusterQW[numClusters-1] == 0))
		{
			cout << "addCUSTOM_CLUSTER_QW_LEFT(): First and/or last element of clusterQW[-] cannot be zero" << endl;
			exit(-1);
		}
	}
	
	// Check the input for consistencty
	 if (dx0_qw < 0) {
		cout << "addCUSTOM_CLUSTER_QW_LEFT(): dx0_qw < 0, quitting" << endl;
		exit(-1);
	} else if (dx1_qw < 0) {
		cout << "addCUSTOM_CLUSTER_QW_LEFT(): dx1_qw <= 0, quitting" << endl;
		exit(-1);
	} else if (numQW < 0) {
		cout << "addCUSTOM_CLUSTER_QW_LEFT(): numQW < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the QW, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addCUSTOM_CLUSTER_QW_LEFT: No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}
	
	if (((numQW == 0)&&(dx0_qw== 0))&&(dx1_qw == 0))
	{
		cout << "addCUSTOM_CLUSTER_QW_LEFT(): No QW's and NO empty nodes, refusing to make empty gain region" << endl;
		exit(-1);
	}
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;
	
	double QW_START = 0;
	
	double EXTRA_PADDING = 0.5;
	// Count number of empty peaks before first QW
	for(unsigned i = 0; i < numClusters; i++)
	{
		if (clusterQW[i] == 0)
		{
			EXTRA_PADDING += 0.5;
		} else {
			// Find distance to first peak
			QW_START 	= BARRIER_LENGTH*(clusterQW[i]-1)/2;
			break;
		}
	}
	
	
	
	
	
	tmpWidth = getLambda()*(dx0_qw*0.5 + EXTRA_PADDING)/cavityIndex - QW_START;
	tmp = addCavity(tmpWidth, cavityIndex);
	
	double QW_MID1 = QW_START;	// Distance from last QW to node
	double QW_MID2 = 0;			// Distance from node to first QW
	
	if (numQW == 0)
	{
		// CASE 1: EMPTY Cavity
		if (dx1_qw > 0)
		{
			tmpWidth = getLambda()*(dx1_qw)/cavityIndex;
			tmp = addCavity(tmpWidth, cavityIndex);
		}
		
	} else {
		
		// CASE 2: numQW > 0
		int QW_NUM = 1;
		// Plan structure
		for(int i = 0; i < numClusters; i++)
		{
			if (clusterQW[i] > 0)
			{
				EXTRA_PADDING = 0.0;
				for(unsigned j = 0; j < clusterQW[i]; j++)
				{
					tmp = addDevice();
					tmpName.str("");
					tmpName << "QW" << QW_NUM;
					tmp->getDevice()->setName(tmpName.str());
					tmp->setOutputToFile(1);
					tmp->getDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
	
					// Add cavity of length
					if (j < clusterQW[i]-1)
					{
						tmp = addCavity(BARRIER_LENGTH, cavityIndex);
					}
	
					QW_NUM += 1;
				}
		
		
				QW_MID1 = getLambda()*0.25/cavityIndex - BARRIER_LENGTH*(clusterQW[i]-1)/2;

				if (i < numClusters-1)
				{
					int kk = i+1;
					while (clusterQW[kk] == 0)
					{
						EXTRA_PADDING += 0.5;
						kk++;
					}
					QW_MID2 = getLambda()*0.25/cavityIndex - BARRIER_LENGTH*(clusterQW[kk]-1)/2;
					tmp = addCavity(getLambda()*EXTRA_PADDING/cavityIndex + QW_MID1 + QW_MID2, cavityIndex);
				} 
			} else {
				//
			}
		}


		// Add last layer
		tmpWidth = getLambda()*(dx1_qw*0.5)/cavityIndex + QW_MID1;
		tmp = addCavity(tmpWidth, cavityIndex);
	}
	

	// Add CAP layer
/*		
	if (capIndex > 0)
	{
		tmpWidth = (getLambda()/2.0)/capIndex;
		tmp = addCavity(tmpWidth, capIndex);
		tmp->getCavity()->setName("Cap");
	} 
	// Add AR coating			
	if (arIndex > 0)
	{
		tmpWidth = (getLambda()/4.0)/arIndex;
		tmp = addCavity(tmpWidth, arIndex);
		tmp->getCavity()->setName("Ar");
	}
*/
	tmpWidth = 0.503490*getLambda()/3.1778;
        tmp = addCavity(tmpWidth, 3.1778);
        tmp->getCavity()->setName("Cap");

        tmpWidth = 0.190713*getLambda()/2.0781;
        tmp = addCavity(tmpWidth, 2.0781);
        tmp->getCavity()->setName("Ta2O5");

        tmpWidth = 0.111396*getLambda()/1.45;
        tmp = addCavity(tmpWidth, 1.45);
        tmp->getCavity()->setName("SiO2");


/*
	tmpWidth = 0.117717*(getLambda())/2.963510;
	tmp = addCavity(tmpWidth, 2.963510);
	tmp->getCavity()->setName("Ar AlAs");

	tmpWidth = 0.351447*(getLambda())/3.396376;
	tmp = addCavity(tmpWidth, 3.396376);
	tmp->getCavity()->setName("AR AlGaAs");

	tmpWidth = 0.183256*(getLambda())/2.963510;
	tmp = addCavity(tmpWidth, 2.963510);
	tmp->getCavity()->setName("Ar AlAs");

	tmpWidth = 0.099180*(getLambda())/3.396376;
	tmp = addCavity(tmpWidth, 3.396376);
	tmp->getCavity()->setName("AR AlGaAs");

	tmpWidth = 0.322405*(getLambda())/2.963510;
	tmp = addCavity(tmpWidth, 2.963510);
	tmp->getCavity()->setName("Ar AlAs");

	tmpWidth = 0.217130*(getLambda())/3.396376;
	tmp = addCavity(tmpWidth, 3.396376);
	tmp->getCavity()->setName("AR AlGaAs");

	tmpWidth = 0.219079*(getLambda())/2.963510;
	tmp = addCavity(tmpWidth, 2.963510);
	tmp->getCavity()->setName("Ar AlAs");

	tmpWidth = 0.080894*(getLambda())/3.195400;
	tmp = addCavity(tmpWidth, 3.195400);
	tmp->getCavity()->setName("CAP");

	tmpWidth = 0.253982*(getLambda())/1.45;
	tmp = addCavity(tmpWidth, 1.45);
	tmp->getCavity()->setName("AR - SiN");
	*/

}

// Special clustering of QWs where two are centered near the peak of the antinode and two are on the wings
void VECSEL::addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(int numClusters, int *clusterQW, bool USE_BOOSTER, int dx0_qw, double dx1_qw, double cavityIndex, double capIndex, double arIndex, double angle_of_incidence, double external_index)
{
	// Distance between QW that are clustered
	//double BARRIER_LENGTH = 16.35*nm;
	//double BARRIER_LENGTH = 13.0*nm;
	double BARRIER_LENGTH_CENTER = 16.35*nm;
	double BARRIER_LENGTH_WINGS = 24.0*nm;
	double booster_layer_space = 16.35*nm;
	

	// Length of region that has to be centered on the antinode
	double qw_span = BARRIER_LENGTH_CENTER + 2.0*BARRIER_LENGTH_WINGS;
	
	
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 	= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	cavityIndex 	= cavityIndex*sqrt(1.0-a_sin2/(cavityIndex*cavityIndex));
	capIndex 		= capIndex*sqrt(1.0-a_sin2/(capIndex*capIndex));
	arIndex 		= arIndex*sqrt(1.0-a_sin2/(arIndex*arIndex));
	
	
	int numQW = 0;
	for(unsigned i = 0; i < numClusters; i++)
	{
		numQW += clusterQW[i];
		if (clusterQW[i] < 0)
		{
			cout << "addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(): Cannot have NEGATIVE number of QW's" << endl;
			exit(-1);
		}
		if (BARRIER_LENGTH_CENTER + 2.0*BARRIER_LENGTH_WINGS > getLambda()*0.5/cavityIndex)
		{
			cout << "addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(): Too much space for antinode" << endl;
			cout << "BARRIER_CENTER = " << BARRIER_LENGTH_CENTER/nm << " [nm]" << endl;
			cout << "BARRIER_WINGS  = " << BARRIER_LENGTH_WINGS/nm  << " [nm]" << endl;
			exit(-1);
		}
		

		if (clusterQW[i]>4)
		{
			cout << "addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(): Need 4 QWs per antinode. No more. No less." << endl;
			cout << "clusterQW[i] = " << clusterQW[i] << endl;
			exit(-1);
		} else if (clusterQW[i]==1) {
			cout << "addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(): Need 4 QWs per antinode. No more. No less." << endl;
			cout << "clusterQW[i] = " << clusterQW[i] << endl;
			exit(-1);
		} else if (clusterQW[i]==2) {
			cout << "addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(): Need 4 QWs per antinode. No more. No less." << endl;
			cout << "clusterQW[i] = " << clusterQW[i] << endl;
			exit(-1);
		} else if (clusterQW[i]==3) {
			cout << "addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(): Need 4 QWs per antinode. No more. No less." << endl;
			cout << "clusterQW[i] = " << clusterQW[i] << endl;
			exit(-1);
		}
	}
	
	if (numClusters > 0)
	{
		if ((clusterQW[0] == 0)||(clusterQW[numClusters-1] == 0))
		{
			cout << "addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(): First and/or last element of clusterQW[-] cannot be zero" << endl;
			exit(-1);
		}
	}
	
	// Check the input for consistencty
	 if (dx0_qw < 0) {
		cout << "addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(): dx0_qw < 0, quitting" << endl;
		exit(-1);
	} else if (dx1_qw < 0) {
		cout << "addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(): dx1_qw <= 0, quitting" << endl;
		exit(-1);
	} else if (numQW < 0) {
		cout << "addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(): numQW < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the QW, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT: No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}
	
	if (((numQW == 0)&&(dx0_qw== 0))&&(dx1_qw == 0))
	{
		cout << "addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(): No QW's and NO empty nodes, refusing to make empty gain region" << endl;
		exit(-1);
	}
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;
	
	double QW_START = 0;
	
	double EXTRA_PADDING = 0.5;
	// Count number of empty peaks before first QW
	for(unsigned i = 0; i < numClusters; i++)
	{
		if (clusterQW[i] == 0)
		{
			EXTRA_PADDING += 0.5;
		} else {
			// Find distance to first peak
			QW_START 	= qw_span/2.0;
			break;
		}
	}
	
	if (USE_BOOSTER) // == true
	{
		tmp = addCavity(booster_layer_space, cavityIndex);

		// Include BOOSTER QW near DBR
		tmp = addDevice();
		tmpName.str("");
		tmpName << "QW1";
		tmp->getDevice()->setName(tmpName.str());
		tmp->setOutputToFile(1);

		tmpWidth = getLambda()*(dx0_qw*0.5 + EXTRA_PADDING)/cavityIndex - QW_START - booster_layer_space;
		tmp = addCavity(tmpWidth, cavityIndex);
		
	} else {
	
		tmpWidth = getLambda()*(dx0_qw*0.5 + EXTRA_PADDING)/cavityIndex - QW_START;
		tmp = addCavity(tmpWidth, cavityIndex);
	}
	
	double QW_MID1 = QW_START;	// Distance from last QW to node
	double QW_MID2 = 0;			// Distance from node to first QW
	
	if (numQW == 0)
	{
		// CASE 1: EMPTY Cavity
		if (dx1_qw > 0)
		{
			tmpWidth = getLambda()*(dx1_qw)/cavityIndex;
			tmp = addCavity(tmpWidth, cavityIndex);
		}
		
	} else {
		
		// CASE 2: numQW > 0
		int QW_NUM = 1;
		if (USE_BOOSTER) // == true
		{
			QW_NUM = 2;
		}
		// Plan structure
		for(int i = 0; i < numClusters; i++)
		{
			if (clusterQW[i] > 0)
			{
				EXTRA_PADDING = 0.0;

				if (clusterQW[i]==4)
				{
					tmp = addDevice();
					tmpName.str("");
					tmpName << "QW" << QW_NUM;
					tmp->getDevice()->setName(tmpName.str());
					tmp->setOutputToFile(1);
					QW_NUM += 1;


					tmp = addCavity(BARRIER_LENGTH_WINGS, cavityIndex);

					tmp = addDevice();
					tmpName.str("");
					tmpName << "QW" << QW_NUM;
					tmp->getDevice()->setName(tmpName.str());
					tmp->setOutputToFile(1);
					QW_NUM += 1;

					tmp = addCavity(BARRIER_LENGTH_CENTER, cavityIndex);
					
					tmp = addDevice();
					tmpName.str("");
					tmpName << "QW" << QW_NUM;
					tmp->getDevice()->setName(tmpName.str());
					tmp->setOutputToFile(1);
					QW_NUM += 1;

					tmp = addCavity(BARRIER_LENGTH_WINGS, cavityIndex);

					tmp = addDevice();
					tmpName.str("");
					tmpName << "QW" << QW_NUM;
					tmp->getDevice()->setName(tmpName.str());
					tmp->setOutputToFile(1);
					QW_NUM += 1;
				}

		
		
				QW_MID1 = getLambda()*0.25/cavityIndex - qw_span/2.0;

				if (i < numClusters-1)
				{
					int kk = i+1;
					while (clusterQW[kk] == 0)
					{
						EXTRA_PADDING += 0.5;
						kk++;
					}
					QW_MID2 = getLambda()*0.25/cavityIndex - qw_span/2.0;
					tmp = addCavity(getLambda()*EXTRA_PADDING/cavityIndex + QW_MID1 + QW_MID2, cavityIndex);
				} 
			} else {
				//
			}
		}


		// Add last layer
		tmpWidth = getLambda()*(dx1_qw*0.5)/cavityIndex + QW_MID1;
		tmp = addCavity(tmpWidth, cavityIndex);
	}
	

/*		
	// Add CAP layer
	if (capIndex > 0)
	{
		tmpWidth = (getLambda()/2.0)/capIndex;
		tmp = addCavity(tmpWidth, capIndex);
		tmp->getCavity()->setName("Cap");
	} 
	// Add AR coating			
	if (arIndex > 0)
	{
		tmpWidth = (getLambda()/4.0)/arIndex;
		tmp = addCavity(tmpWidth, arIndex);
		tmp->getCavity()->setName("Ar");
	}
*/

	tmpWidth = 0.486678*getLambda()/3.1778;
	tmp = addCavity(tmpWidth, 3.1778);
	tmp->getCavity()->setName("Cap");

	tmpWidth = 0.175184*getLambda()/2.0781;
	tmp = addCavity(tmpWidth, 2.0781);
	tmp->getCavity()->setName("Ar Perfect");

	tmpWidth = 0.143084*getLambda()/1.45;
	tmp = addCavity(tmpWidth, 1.45);
	tmp->getCavity()->setName("Ar Perfect");


/*
	tmpWidth = 0.117717*(getLambda())/2.963510;
	tmp = addCavity(tmpWidth, 2.963510);
	tmp->getCavity()->setName("Ar AlAs");

	tmpWidth = 0.351447*(getLambda())/3.396376;
	tmp = addCavity(tmpWidth, 3.396376);
	tmp->getCavity()->setName("AR AlGaAs");

	tmpWidth = 0.183256*(getLambda())/2.963510;
	tmp = addCavity(tmpWidth, 2.963510);
	tmp->getCavity()->setName("Ar AlAs");

	tmpWidth = 0.099180*(getLambda())/3.396376;
	tmp = addCavity(tmpWidth, 3.396376);
	tmp->getCavity()->setName("AR AlGaAs");

	tmpWidth = 0.322405*(getLambda())/2.963510;
	tmp = addCavity(tmpWidth, 2.963510);
	tmp->getCavity()->setName("Ar AlAs");

	tmpWidth = 0.217130*(getLambda())/3.396376;
	tmp = addCavity(tmpWidth, 3.396376);
	tmp->getCavity()->setName("AR AlGaAs");

	tmpWidth = 0.219079*(getLambda())/2.963510;
	tmp = addCavity(tmpWidth, 2.963510);
	tmp->getCavity()->setName("Ar AlAs");

	tmpWidth = 0.080894*(getLambda())/3.195400;
	tmp = addCavity(tmpWidth, 3.195400);
	tmp->getCavity()->setName("CAP");

	tmpWidth = 0.253982*(getLambda())/1.45;
	tmp = addCavity(tmpWidth, 1.45);
	tmp->getCavity()->setName("AR - SiN");
	*/

}


/*
 *  Create a gain medium on the left with QW. 
 * 
 * 	|---o-o-o-o-o-------------|--CAP--|--AR--|
 * 
 *   . If there is nothing on the left of this, a boundary is added
 * 
 *  Where:
 *  numQW			-> Number of QW's in the medium [0,inf)
 * 						- If numQW = 0		=> Empty medium of length cavityLength
 * 	dx0_qw			-> Distance from LEFTMOST QW to the edge in units of m
 *  dx_qw			-> Distance between each QW in units of m
 *  cavityLength	-> Total length of cavity (not including AR and CAP)
 * 						- If cavityLength = 0    => Fit cavity tight on QW's (If numQW = 0 => ERROR)
 * 						- If cavityLength < space for QW with spacing   => ERROR
 * 	cavityIndex		-> Refractive background index in medium
 * 	capIndex		-> Refractive index of CAP layer, if capIndex = 0 then NO CAP layer
 *  arIndex			-> Refractive index of AR coating, if arIndex = 0 then NO ar coating
 *  fillerLength	-> Length of cavity to fill in AFTER AR and CAP layer in units of m
 * */


void VECSEL::addCUSTOM_QW_LEFT_OPTIMAL_ANGLE(int numQW, double *width, double cavityLength, double cavityIndex, double capIndex, double arIndex, double arLength,bool PASSIVE_STRUCTURE,double BAD_GROWTH_FACTOR,double BAD_GROWTH_FACTOR_CAP, double angle_of_incidence, double external_index, double ar_temperature_index_diff)
{

	double BARRIER = 5*nm;
	// Check the input for consistencty
	if (cavityLength < 0)
	{
		cout << "addCUSTOM_QW_LEFT(): cavityLength < 0, quitting" << endl;
		exit(-1);
	} else if (numQW < 0) {
		cout << "addCUSTOM_QW_LEFT(): numQW < 0, quitting" << endl;
		exit(-1);
	}

	for(unsigned i = 0; i < numQW; i++)
	{
		if (width[i] <=0.0)
		{
			cout << "addCUSTOM_QW_LEFT(): width[i] <= 0" << endl;
			exit(-1);
		} else if (width[i] < BARRIER)
		{
			cout << "addCUSTOM_QW_LEFT(): width[i] < min BARRIER length" << endl;
			exit(-1);
		}
	}
	
	// Check a few logical inconsistencies
	if ((numQW == 0) && (cavityLength == 0))
	{
		cout << "addCUSTOM_QW_LEFT(): numQW =0 && cavityLength = 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the QW, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addCUSTOM_QW_LEFT: No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}

	
	Module *tmp;
	double tmpWidth;
	
	
	
	std::stringstream tmpName;
	if (numQW == 0)
	{
		// CASE 1: EMPTY Cavity
		tmpWidth = cavityLength;
		tmp = addCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		
	} else {
		
		// CASE 2: numQW > 0
		// First cavity
		tmpWidth = BAD_GROWTH_FACTOR*width[0];
		tmp = addCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		
		
		// Add first device
		tmp = addDevice();
		if (PASSIVE_STRUCTURE)
		{
			tmp->getDevice()->setName("QW1_PASSIVE");
		} else {
			tmp->getDevice()->setName("QW1");
		}
		tmp->setOutputToFile(1);

		double totalWidth = BAD_GROWTH_FACTOR*width[0];
		// Add rest in dx_qw distance
		for(int i = 1; i < numQW; i++)
		{
			// Add Cavity
			totalWidth += BAD_GROWTH_FACTOR*width[i];
			tmp = addCavity(BAD_GROWTH_FACTOR*width[i], cavityIndex, angle_of_incidence, external_index);

			
			// Add device
			tmp = addDevice();
			if (PASSIVE_STRUCTURE)
			{
				tmp->getDevice()->setName("QW1_PASSIVE");
				
			} else {
				tmpName.str("");
				tmpName << "QW" << i+1;
				tmp->getDevice()->setName(tmpName.str());
			}
			tmp->setOutputToFile(1);
		}
		
		if (cavityLength > 0)
		{
			// Add rest of Cavity
			tmpWidth = BAD_GROWTH_FACTOR*(cavityLength - totalWidth);
			tmp = addCavity( tmpWidth , cavityIndex, angle_of_incidence, external_index);

		} else {
			// Add final barrier
			//tmpWidth = BARRIER;
			tmp = addCavity( BAD_GROWTH_FACTOR*width[numQW] , cavityIndex, angle_of_incidence, external_index);
		
		
		}
	}
/*
	// Add CAP layer
	if (capIndex > 0)
	{
		tmpWidth = BAD_GROWTH_FACTOR_CAP*(getLambda()/2.0)/capIndex;
		//tmpWidth = BAD_GROWTH_FACTOR*(980*nm/2.0)/capIndex;
		tmp = addCavity(tmpWidth, capIndex, angle_of_incidence, external_index);
		tmp->getCavity()->setName("Cap");
	} 
		
	// Add AR coating	
	if (arIndex > 0)
	{
		//tmpWidth = BAD_GROWTH_FACTOR_AR*(getLambda()/4.0)/arIndex;
		//tmpWidth = BAD_GROWTH_FACTOR_AR*(140.0*nm); // Target length = 140*nm, Alex estimates width = 125*nm
		tmp = addCavity(arLength, arIndex, angle_of_incidence, external_index);
		tmp->getCavity()->setName("Ar");
	} 
*/
	double n_AlGaAs = 3.4583 +ar_temperature_index_diff;
	double n_AlAs   = 2.9601 +ar_temperature_index_diff;
	double n_cap    = 3.1778 +ar_temperature_index_diff;
	double n_si02   = 1.4708; // Lemarchand 2013

	// CAP
	tmpWidth = 128.96*nm; // 127.38
	tmp = addCavity(tmpWidth, 3.1978, angle_of_incidence, external_index);
	tmp->getCavity()->setName("CAP");

	// AR 2
	tmpWidth = 87.88*nm; // 87.77nm
	tmp = addCavity(tmpWidth, 1.9826, angle_of_incidence, external_index);
	tmp->getCavity()->setName("Ar");

	// AR 2
	tmpWidth = 118.98*nm; // 119.23nm
	tmp = addCavity(tmpWidth, 1.477, angle_of_incidence, external_index);
	tmp->getCavity()->setName("Ar");
/*
	// CAP
	tmpWidth = 0.273099*getLambda()/3.1978;
	tmp = addCavity(tmpWidth, 3.1978, angle_of_incidence, external_index);
	tmp->getCavity()->setName("CAP");

	// AR 2
	tmpWidth = 0.117960*getLambda()/1.9826;
	tmp = addCavity(tmpWidth, 1.9826, angle_of_incidence, external_index);
	tmp->getCavity()->setName("Ar");

	// AR 2
	tmpWidth = 0.206474*getLambda()/1.477;
	tmp = addCavity(tmpWidth, 1.477, angle_of_incidence, external_index);
	tmp->getCavity()->setName("Ar");
*/

	//====================================================================
	// Add in angle dependence for ALL previous layers (DBR/gold/...)
	// First cavity
	std::complex<double> n_curr = modules[quick_index_cavity[quick_index_cavity.size()-1]].getRefInd() + I*modules[quick_index_cavity[quick_index_cavity.size()-1]].getRefInd_im();
	std::complex<double> n_prev = external_index;

	// Convert to standard form
	double n2_tilde = real(n_prev)*real(n_curr) + imag(n_prev)*imag(n_curr);
	double k2_tilde = real(n_prev)*imag(n_curr) - real(n_curr)*imag(n_prev);
	double n1_tilde = abs(n_prev)*abs(n_prev);

	// Compute temporary terms
	double n1_sin_th_2 = n1_tilde*sin(angle_of_incidence); // na*sin(theta) 
	n1_sin_th_2 *= n1_sin_th_2; // ^2
	double norm_n_2 = n2_tilde*n2_tilde + k2_tilde*k2_tilde;
	double term_a = 1.0+n1_sin_th_2/norm_n_2;

	// Kovalenko 2001: Descartes-Snell law of refraction with absorption
	// compute sin(theta)^2
	double sin_th_2 = 0.5*(term_a - sqrt(term_a*term_a - 4.0*n2_tilde*n2_tilde*n1_sin_th_2/(norm_n_2*norm_n_2)));
	double cos_th = sqrt(1.0-sin_th_2);
	modules[quick_index_cavity[quick_index_cavity.size()-1]].setCosTh(cos_th,cos_th);

	// Iterate over all cavities
	for(int i=quick_index_cavity.size()-2; i>=0; i--)
	{
		n_curr = modules[quick_index_cavity[i]].getRefInd() + I*modules[quick_index_cavity[i]].getRefInd_im();
		n_prev = modules[quick_index_cavity[i+1]].getRefInd() + I*modules[quick_index_cavity[i+1]].getRefInd_im();

		// Convert to standard form
		n2_tilde = real(n_prev)*real(n_curr) + imag(n_prev)*imag(n_curr);
		k2_tilde = real(n_prev)*imag(n_curr) - real(n_curr)*imag(n_prev);
		n1_tilde = abs(n_prev)*abs(n_prev);

		// Compute temporary terms
		n1_sin_th_2 = n1_tilde*n1_tilde*sin_th_2; // (n1 sin(th))^2
		norm_n_2 = n2_tilde*n2_tilde + k2_tilde*k2_tilde;
		term_a = 1.0+n1_sin_th_2/norm_n_2;

		// Kovalenko 2001: Descartes-Snell law of refraction with absorption
		// compute sin(theta)^2
		sin_th_2 = 0.5*(term_a - sqrt(term_a*term_a - 4.0*n2_tilde*n2_tilde*n1_sin_th_2/(norm_n_2*norm_n_2)));
		cos_th = sqrt(1.0-sin_th_2);
		modules[quick_index_cavity[i]].setCosTh(cos_th,cos_th);
	}
}

/*
 *  Create a gain medium on the left with QW. 
 * 
 * 	|---o-o-o-o-o-------------|--CAP--|--AR--|
 * 
 *   . If there is nothing on the left of this, a boundary is added
 * 
 *  Where:
 *  numQW			-> Number of QW's in the medium [0,inf)
 * 						- If numQW = 0		=> Empty medium of length cavityLength
 * 	dx0_qw			-> Distance from LEFTMOST QW to the edge in units of lambda [0,1]
 *  dx_qw			-> Distance between each QW in units of lambda (0,1]
 *  cavityLength	-> Total length of cavity (not including AR and CAP)
 * 						- If cavityLength = 0    => Fit cavity tight on QW's (If numQW = 0 => ERROR)
 * 						- If cavityLength < space for QW with spacing   => ERROR
 * 	cavityIndex		-> Refractive background index in medium
 * 	capIndex		-> Refractive index of CAP layer, if capIndex = 0 then NO CAP layer
 *  arIndex			-> Refractive index of AR coating, if arIndex = 0 then NO ar coating
 *  angle_of_incidence -> Angle in radians for TE
 *  external_index  -> Index in incoming medium
 * */

void VECSEL::addRPG_QW_LEFT(int numQW, double dx0_qw, double dx_qw, double cavityLength, double cavityIndex, double capIndex, double arIndex, double angle_of_incidence, double external_index)
{
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 	= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	cavityIndex 	= cavityIndex*sqrt(1.0-a_sin2/(cavityIndex*cavityIndex));
	capIndex 		= capIndex*sqrt(1.0-a_sin2/(capIndex*capIndex));
	arIndex 		= arIndex*sqrt(1.0-a_sin2/(arIndex*arIndex));
	
	
	// Check the input for consistencty
	if (cavityLength < 0)
	{
		cout << "addRPG_QW_LEFT(): cavityLength < 0, quitting" << endl;
		exit(-1);
	} else if (dx0_qw < 0) {
		cout << "addRPG_QW_LEFT(): dx0_qw < 0, quitting" << endl;
		exit(-1);
	} else if (dx_qw < 0) {
		cout << "addRPG_QW_LEFT(): dx_qw <= 0, quitting" << endl;
		exit(-1);
	} else if (numQW < 0) {
		cout << "addRPG_QW_LEFT(): numQW < 0, quitting" << endl;
		exit(-1);
	}
	
	// Check a few logical inconsistencies
	if ((numQW == 0) && (cavityLength == 0))
	{
		cout << "addRPG_QW_LEFT(): numQW =0 && cavityLength = 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the QW, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addRPG_QW_LEFT: No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}
	
	if ((numQW > 0)&&(cavityLength>0))
	{
		double tmpWidth = dx0_qw + ((double)numQW -1.0)*dx_qw;
		if (tmpWidth > cavityLength)
		{
			cout << "addRPG_QW_LEFT(): cavityLength too short, quitting" << endl;
			exit(-1);
		}
		
		if ((dx_qw == 0)&&(numQW>1))
		{
			cout << "addRPG_QW_LEFT(): cannot have numQW>0 && dx_qw = 0, quitting" << endl;
			exit(-1);
		}
	}
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;
	if (numQW == 0)
	{
		// CASE 1: EMPTY Cavity
		tmpWidth = (getLambda()*cavityLength)/cavityIndex;
		tmp = addCavity(tmpWidth, cavityIndex);
		
	} else {
		
		// CASE 2: numQW > 0
		// First cavity
		if (dx0_qw > 0)
		{
			tmpWidth = (getLambda()*dx0_qw)/cavityIndex;
			tmp = addCavity(tmpWidth, cavityIndex);
		}
		
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add first device
			tmp = addDevice();
			tmpName.str("");
			tmpName << "QW1_T" << j+1;
			tmp->getDevice()->setName(tmpName.str());
			tmp->setOutputToFile(1);
			tmp->getDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
		}

		// Add rest in dx_qw distance
		double deviceDistance = (getLambda()*dx_qw)/(cavityIndex);
		for(int i = 1; i < numQW; i++)
		{
			// Add Cavity
			tmp = addCavity(deviceDistance, cavityIndex);
			
			for(int j = 0; j < VECSEL_transverse_points_number; j++)
			{
				// Add device
				tmp = addDevice();
				tmpName.str("");
				tmpName << "QW" << i+1 << "_T" << j+1;
				tmp->getDevice()->setName(tmpName.str());
				tmp->setOutputToFile(1);
				tmp->getDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
			}
		}
		
		if (cavityLength > 0)
		{
			// Add rest of Cavity
			tmpWidth = getLambda()*(cavityLength - dx0_qw - ((double)numQW -1.0)*dx_qw)/cavityIndex;
			tmp = addCavity( tmpWidth , cavityIndex);
		}
	}
	
/*
	// Add CAP layer
	if (capIndex > 0)
	{
		tmpWidth = (getLambda()/2.0)/capIndex;
		tmp = addCavity(tmpWidth, capIndex);
		tmp->getCavity()->setName("Cap");
	} 

	// Add AR coating	
	if (arIndex > 0)
	{
		tmpWidth = (0.25*getLambda())/arIndex;
		tmp = addCavity(tmpWidth, arIndex);
		tmp->getCavity()->setName("Ar");
	} 
*/
	// CAP
	tmpWidth = 0.432421*(getLambda())/3.1778;
	tmp = addCavity(tmpWidth, 3.1778);
	tmp->getCavity()->setName("cap");
	// AR 1
	tmpWidth = 0.181938*(getLambda())/2.0781;
	tmp = addCavity(tmpWidth, 2.0781);
	tmp->getCavity()->setName("Ta2O5");

	// AR 2
	tmpWidth = 0.146560*(getLambda())/1.45;
	tmp = addCavity(tmpWidth, 1.45);
	tmp->getCavity()->setName("SiO2");

/*
	// AR 4
	tmpWidth = (getLambda()/2.0)/3.756268;
	tmp = addCavity(tmpWidth, 3.756268);
	tmp->getCavity()->setName("Ar4");
	
	// AR 3
	tmpWidth = (getLambda()/4.0)/2.870939;
	tmp = addCavity(tmpWidth, 2.870939);
	tmp->getCavity()->setName("Ar3");
	
	// AR 2
	tmpWidth = (getLambda()/4.0)/1.694157;
	tmp = addCavity(tmpWidth, 1.694157);
	tmp->getCavity()->setName("Ar2");
	
	// AR 1
	tmpWidth = (getLambda()/4.0)/1.126036;
	tmp = addCavity(tmpWidth, 1.126036);
	tmp->getCavity()->setName("Ar1");
*/
	
}

/* Similar to addRPG_QW_LEFT() however, this function will only add ONE qw at a given antinode and fill the rest with empty barrier material
*/

void VECSEL::addRPG_QW_LEFT_EFF(int anti_node, int numQW, double dx0_qw, double dx_qw, double cavityLength, double cavityIndex, double capIndex, double arIndex, double angle_of_incidence, double external_index)
{
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 	= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	cavityIndex 	= cavityIndex*sqrt(1.0-a_sin2/(cavityIndex*cavityIndex));
	capIndex 		= capIndex*sqrt(1.0-a_sin2/(capIndex*capIndex));
	arIndex 		= arIndex*sqrt(1.0-a_sin2/(arIndex*arIndex));
	
	
	// Check the input for consistencty
	if (cavityLength < 0)
	{
		cout << "addRPG_QW_LEFT_EFF(): cavityLength < 0, quitting" << endl;
		exit(-1);
	} else if (dx0_qw < 0) {
		cout << "addRPG_QW_LEFT_EFF(): dx0_qw < 0, quitting" << endl;
		exit(-1);
	} else if (dx_qw < 0) {
		cout << "addRPG_QW_LEFT_EFF(): dx_qw <= 0, quitting" << endl;
		exit(-1);
	} else if (numQW < 0) {
		cout << "addRPG_QW_LEFT_EFF(): numQW < 0, quitting" << endl;
		exit(-1);
	} else if ((anti_node < 1) || (anti_node > numQW)) {
		cout << "addRPG_QW_LEFT_EFF(): anti_node must be in range [1, numQW]" << endl;
		exit(-1);
	}
	
	// Check a few logical inconsistencies
	if ((numQW == 0) && (cavityLength == 0))
	{
		cout << "addRPG_QW_LEFT(): numQW =0 && cavityLength = 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the QW, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addRPG_QW_LEFT: No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}
	
	if ((numQW > 0)&&(cavityLength>0))
	{
		double tmpWidth = dx0_qw + ((double)numQW -1.0)*dx_qw;
		if (tmpWidth > cavityLength)
		{
			cout << "addRPG_QW_LEFT(): cavityLength too short, quitting" << endl;
			exit(-1);
		}
		
		if ((dx_qw == 0)&&(numQW>1))
		{
			cout << "addRPG_QW_LEFT(): cannot have numQW>0 && dx_qw = 0, quitting" << endl;
			exit(-1);
		}
	}
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;
	if (numQW == 0)
	{
		// CASE 1: EMPTY Cavity
		tmpWidth = (getLambda()*cavityLength)/cavityIndex;
		tmp = addCavity(tmpWidth, cavityIndex);
		
	} else {
		
		// CASE 2: numQW > 0
		// First cavity
		if (dx0_qw + (anti_node-1)*dx_qw > 0)
		{
			tmpWidth = getLambda()*(dx0_qw + (anti_node-1)*dx_qw)/cavityIndex;
			tmp = addCavity(tmpWidth, cavityIndex);
			tmp->getCavity()->setName("FrontMatQW");
			//tmp->setOutputToFile(1);
		}
		
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add first device
			tmp = addDevice();
			tmpName.str("");
			tmpName << "QW" << anti_node << "_T" << j+1;
			tmp->getDevice()->setName(tmpName.str());
			//tmp->setOutputToFile(1);
			tmp->getDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
			
		}
		
		double delta_L = 0.0;
		if (cavityLength > 0)
		{
			delta_L = cavityLength - dx0_qw - ((double)numQW -1.0)*dx_qw;
		} else {
			delta_L = 0.0;
		}

		if (anti_node < numQW)
		{
			tmpWidth = getLambda()*((numQW-anti_node)*dx_qw + delta_L)/cavityIndex;
			tmp = addCavity(tmpWidth, cavityIndex);
		} else {
			if (delta_L > 0)
			{
				// Add rest of Cavity
				tmpWidth = getLambda()*(delta_L)/cavityIndex;
				tmp = addCavity( tmpWidth , cavityIndex);
			}
		}
	}
	
/*
	// Add CAP layer
	if (capIndex > 0)
	{
		tmpWidth = (getLambda()/2.0)/capIndex;
		tmp = addCavity(tmpWidth, capIndex);
		tmp->getCavity()->setName("Cap");
	} 

	// Add AR coating	
	if (arIndex > 0)
	{
		tmpWidth = (0.25*getLambda())/arIndex;
		tmp = addCavity(tmpWidth, arIndex);
		tmp->getCavity()->setName("Ar");
	} 
*/
	// CAP
	tmpWidth = 0.432421*(getLambda())/3.1778;
	tmp = addCavity(tmpWidth, 3.1778);
	tmp->getCavity()->setName("cap");
	//tmp->setOutputToFile(1);
	
	// AR 1
	tmpWidth = 0.181938*(getLambda())/2.0781;
	tmp = addCavity(tmpWidth, 2.0781);
	tmp->getCavity()->setName("Ta2O5");
	//tmp->setOutputToFile(1);

	// AR 2
	tmpWidth = 0.146560*(getLambda())/1.45;
	tmp = addCavity(tmpWidth, 1.45);
	tmp->getCavity()->setName("SiO2");
//	tmp->setOutputToFile(1);

/*
	// AR 4
	tmpWidth = (getLambda()/2.0)/3.756268;
	tmp = addCavity(tmpWidth, 3.756268);
	tmp->getCavity()->setName("Ar4");
	
	// AR 3
	tmpWidth = (getLambda()/4.0)/2.870939;
	tmp = addCavity(tmpWidth, 2.870939);
	tmp->getCavity()->setName("Ar3");
	
	// AR 2
	tmpWidth = (getLambda()/4.0)/1.694157;
	tmp = addCavity(tmpWidth, 1.694157);
	tmp->getCavity()->setName("Ar2");
	
	// AR 1
	tmpWidth = (getLambda()/4.0)/1.126036;
	tmp = addCavity(tmpWidth, 1.126036);
	tmp->getCavity()->setName("Ar1");
*/
	
}

/* Similar to addRPG_QW_LEFT() however, this function will only add ONE qw at a given antinode and fill the rest with empty barrier material
*/

void VECSEL::addRPG_QW_simple(int numQW, double cavityLength1, double cavityLength2, double cavityIndex, double angle_of_incidence, double external_index)
{
		
	// Check the input for consistencty
	if (cavityLength1<0 || cavityLength2<0)
	{
		cout<< "addRPG_QW_simple(): cavity lengths cannot be negative. Quitting." <<endl;
		exit(-1);	
	}
	// Check a few logical inconsistencies
	if ((numQW == 0) && (cavityLength1+cavityLength2 == 0))
	{
		cout << "addRPG_QW_simple(): numQW =0 && cavityLength = 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the simple of the QW, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addRPG_QW_simple: No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}
	
	
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;
/*	tmpWidth = getLambda();
	tmp = addCavity(tmpWidth, 100.0);
	tmp->setOutputToFile(1); 
	//tmp->getCavity()->setName("BPMFS2");
	tmp->getCavity()->setName("TESTCAV");*/
		
	if (numQW == 0)
	{
		// CASE 1: EMPTY Cavity with probe
		tmpWidth = getLambda()*cavityLength1/cavityIndex;
		tmp = addCavity(tmpWidth, cavityIndex);
		tmp->getCavity()->setName("CAVBACK");
		tmp->setOutputToFile(2); 
		
		tmpWidth = getLambda()*cavityLength2/cavityIndex;
		tmp = addCavity(tmpWidth, cavityIndex);
		tmp->getCavity()->setName("CAVFRONT");
		tmp->setOutputToFile(2); 
	} else {
		
		// CASE 2: numQW > 0
		// First cavity
			
		tmpWidth = getLambda()*cavityLength1/cavityIndex;
		tmp = addCavity(tmpWidth, cavityIndex);
		tmp->getCavity()->setName("CAVBACK");
		tmp->setOutputToFile(2); 
		
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add first device
			tmp = addDevice();
			tmpName.str("");
			tmpName << "QW6" << "_T" << j+1;
			tmp->getDevice()->setName(tmpName.str());
			tmp->setOutputToFile(2); 
			tmp->getDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
		}
			tmpWidth = getLambda()*cavityLength2/cavityIndex;
			tmp = addCavity(tmpWidth, cavityIndex);
			tmp->getCavity()->setName("CAVFRONT");	
			tmp->setOutputToFile(2); 
		

	}	
}

/* Adds a basic twoArm cavity primarily for debugging purpoposes. Also adds a two arm interface in front if needed*/
/*void VECSEL::addTwoArmCavity(double cavityLength, double cavityIndex, double angle_of_incidence, double external_index)
{
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 	= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	cavityIndex 	= cavityIndex*sqrt(1.0-a_sin2/(cavityIndex*cavityIndex));

	if (cavityLength < 0)
	{
		cout << "addTwoArmCavity(): cavityLength < 0, quitting" << endl;
		exit(-1);
	} else if (!(modules.back().isCavity()))
	{
		cout << "addTwoArmCavity(): Previous module must be a cavity, quitting" << endl;
	}

	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;

	double indx=getNumberModules();
	
	tmpWidth = getLambda();
	tmp = addTwoArmInterface(tmpWidth, cavityIndex, angle_of_incidence, 1.0);	
		
	tmpWidth = cavityLength*getLambda()/cavityIndex;
	tmp = addTwoArmCavity(tmpWidth, cavityIndex);
	tmp->getTwoArmCavity()->setName("TACAV");	
}*/

/* Similar to addRPG_QW_LEFT() however, this function will only add ONE qw at a given antinode and fill the rest with empty barrier material. This is the two arm equivalent.
*/

void VECSEL::addTwoArmRPG_QW_FRONT_EFF(int material, int anti_node, int numQW, double dx0_qw, double dx_qw, double cavityLength, double cavityIndex, double capIndex, double arIndex, double arIndex2, double angle_of_incidence, double external_index)
{
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 		= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	double intIndex 	= cavityIndex;
	cavityIndex 		= cavityIndex*sqrt(1.0-a_sin2/(cavityIndex*cavityIndex));
	capIndex 		= capIndex*sqrt(1.0-a_sin2/(capIndex*capIndex));
	arIndex 		= arIndex*sqrt(1.0-a_sin2/(arIndex*arIndex));
	arIndex2                 = arIndex2*sqrt(1.0-a_sin2/(arIndex2*arIndex2));	

	// Check the input for consistencty
	if (cavityLength < 0)
	{
		cout << "addTwoArmRPG_QW_front_EFF(): cavityLength < 0, quitting" << endl;
		exit(-1);
	} else if (dx0_qw < 0) {
		cout << "addTwoArmRPG_QW_front_EFF(): dx0_qw < 0, quitting" << endl;
		exit(-1);
	} else if (dx_qw < 0) {
		cout << "addTwoArmRPG_QW_front_EFF(): dx_qw <= 0, quitting" << endl;
		exit(-1);
	} else if (numQW < 0) {
		cout << "addTwoArmRPG_QW_front_EFF(): numQW < 0, quitting" << endl;
		exit(-1);
	} else if (anti_node > numQW) {
		cout << "addTwoArmRPG_QW_front_EFF(): anti_node must be <numQW" << endl;
		exit(-1);
	} else if ((anti_node < 1) && (numQW>0))
	{
		cout << "addTwoArmRPG_QW_front_EFF(): anti_node must be positive" << endl;
		exit(-1);
	}
	
	// Check a few logical inconsistencies
	if ((numQW == 0) && (cavityLength == 0))
	{
		cout << "addTwoArmRPG_QW_front(): numQW =0 && cavityLength = 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the QW, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addTwoArmRPG_QW_front_EFF: No LEFT boundary detected, exiting for now. Not tested." << endl;
		exit(-1);
	} else if (modules[startInd-1].isBoundary())
	{
		cout<< "addTwoArmRPG_QW_front_EFF: TwoArm cavity proceeded by a boundary is untested. Exiting." << endl;
		exit(-1);
	}
	

	if ((numQW > 0)&&(cavityLength>0))
	{
		double tmpWidth = dx0_qw + ((double)numQW -1.0)*dx_qw;
		if (tmpWidth > cavityLength)
		{
			cout << "addTwoArmRPG_QW_front_EFF(): cavityLength too short, quitting" << endl;
			exit(-1);
		}
		
		if ((dx_qw == 0)&&(numQW>1))
		{
			cout << "addTwoArmRPG_QW_front_EFF(): cannot have numQW>0 && dx_qw = 0, quitting" << endl;
			exit(-1);
		}
	}
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;
	//WTF-Interface to start. Requires external index to be 1.0 and does not yet include forced delay. Need to make nonuniform cyclic structures.

	if (modules[startInd-1].isCavity())
	{
		tmpWidth = 1.0*getLambda();
		tmp = addTwoArmInterface(1.0, angle_of_incidence, 1.0, intIndex);	
		//tmp->getTwoArmInterface()->setReflect(-1.0); //Absorber has perfectly reflecting back boundary
	}
	
	// AR 1
	tmpWidth = 0.146560*(getLambda())/arIndex;
	tmp = addTwoArmCavity(tmpWidth, arIndex, angle_of_incidence, external_index);
	tmp->getTwoArmCavity()->setName("SiO2");

	// AR 2
	tmpWidth = 0.181938*(getLambda())/arIndex2;
	tmp = addTwoArmCavity(tmpWidth, arIndex2, angle_of_incidence, external_index);
	tmp->getTwoArmCavity()->setName("Ta2O5");
	
	// CAP
	tmpWidth = 0.432421*(getLambda())/capIndex;
	tmp = addTwoArmCavity(tmpWidth, capIndex, angle_of_incidence, external_index);
	tmp->getTwoArmCavity()->setName("cap");

	if (numQW == 0)
	{
		// CASE 1: EMPTY Cavity
		tmpWidth = (getLambda()*cavityLength)/cavityIndex;
		tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		
	} else {
		
		// CASE 2: numQW > 0
		//Front cavity
		double delta_L = 0.0;
		if (cavityLength > 0)
		{
			delta_L = cavityLength - dx0_qw - ((double)numQW -1.0)*dx_qw;
		} else {
			delta_L = 0.0;
		}

		if (anti_node < numQW)
		{
			tmpWidth = getLambda()*((numQW-anti_node)*dx_qw + delta_L)/cavityIndex;
			tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		} else {
			if (delta_L > 0)
			{
				// Add rest of Cavity
				tmpWidth = getLambda()*(delta_L)/cavityIndex;
				tmp = addTwoArmCavity( tmpWidth , cavityIndex, angle_of_incidence, external_index);
			}
		}
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add first device
			tmp = addTwoArmDevice();
			tmpName.str("");
			tmpName << "QW" << material << "_T" << j+1;
			tmp->getTwoArmDevice()->setName(tmpName.str());
			//WTFF-tmp->setOutputToFile(2);
			tmp->getTwoArmDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
		}
		// Back cavity
		if (dx0_qw + (anti_node-1)*dx_qw > 0)
		{
			tmpWidth = getLambda()*(dx0_qw + (anti_node-1)*dx_qw)/cavityIndex;
			tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		}
		
	}	
}

/* Generates a QW structure with two cavitiies separated by a QW.
*/

void VECSEL::addQW_STRUCT(int numQW, double cavityLength1, double cavityLength2, double cavityIndex, double angle_of_incidence, double external_index)
{
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 		= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	double intIndex 	= cavityIndex;
	cavityIndex 		= cavityIndex*sqrt(1.0-a_sin2/(cavityIndex*cavityIndex));
	
	
	// Check the input for consistencty
	if (cavityLength1 < 0)
	{
		cout << "addQW_STRUCT(): cavityLength1 < 0, quitting" << endl;
		exit(-1);
	} else if (cavityLength2 < 0)
	{
		cout << "addQW_STRUCT(): cavityLength2 < 0, quitting" << endl;
		exit(-1);
	} 	
	// Check a few logical inconsistencies
	if ((numQW == 0) && (cavityLength1 == 0 || cavityLength2 ==0))
	{
		cout << "addQW_STRUCT(): numQW =0 && cavityLength = 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the QW, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addQW_STRUCT: No LEFT boundary detected. TwoArmStructure as first cavity is untested. Exiting" << endl;
		exit(-1);
	}
	
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;
	tmpWidth = getLambda();
	
	if (numQW == 0)
	{
		// CASE 1: EMPTY Cavity
		tmpWidth = (getLambda()*(cavityLength1+cavityLength2))/cavityIndex;
		tmp = addCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		
	} else if (numQW < 0)
	{
		// CASE 2: numQW<0
		//ABSORBING QW
		//First cavity
				
		tmpWidth = getLambda()*cavityLength1/cavityIndex;
		tmp = addCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		//tmp->setOutputToFile(2);
		tmp->getCavity()->setName("CFRONT");
		
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add first device
			tmp = addDevice();
			tmpName.str("");
			tmpName << "ABS1" << "_T" << j+1;
			tmp->getDevice()->setName(tmpName.str());
			tmp->setOutputToFile(2);
			tmp->getDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
		}
			tmpWidth = getLambda()*cavityLength2/cavityIndex;
			tmp = addCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
			tmp->setOutputToFile(1);
			tmp->getCavity()->setName("CBACK");
	} else 
	{	
		// CASE 3: numQW > 0
		// Gain QW
		//First cavity
				
		tmpWidth = getLambda()*cavityLength1/cavityIndex;
		tmp = addCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		//tmp->setOutputToFile(2);
		tmp->getCavity()->setName("CFRONT");
		
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add first device
			tmp = addDevice();
			tmpName.str("");
			tmpName << "QW6" << "_T" << j+1;
			tmp->getDevice()->setName(tmpName.str());
			tmp->setOutputToFile(2);
			tmp->getDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
		}
			tmpWidth = getLambda()*cavityLength2/cavityIndex;
			tmp = addCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
			tmp->setOutputToFile(1);
			tmp->getCavity()->setName("TACBACK");
	}	
}


/* Generates a two arm structure with an interface and two twoArmCavities separated by a QW. Will eventually support full angular implementation 
*/
void VECSEL::addTwoArmQW_STRUCT(int numQW, double cavityLength1, double cavityLength2, double cavityIndex, double angle_of_incidence, double external_index, double reflect)
{
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 		= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	double intIndex 	= cavityIndex;
	cavityIndex 		= cavityIndex*sqrt(1.0-a_sin2/(cavityIndex*cavityIndex));
	
	
	// Check the input for consistencty
	if (cavityLength1 < 0)
	{
		cout << "addTwoArmQW_STRUCT(): cavityLength1 < 0, quitting" << endl;
		exit(-1);
	} else if (cavityLength2 < 0)
	{
		cout << "addTwoArmQW_STRUCT(): cavityLength2 < 0, quitting" << endl;
		exit(-1);
	} 	
	// Check a few logical inconsistencies
	if ((numQW == 0) && (cavityLength1 == 0 || cavityLength2 ==0))
	{
		cout << "addTwoArmQW_STRUCT(): numQW =0 && cavityLength = 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the QW, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addTwoArmQW_STRUCT: No LEFT boundary detected. TwoArmStructure as first cavity is untested. Exiting" << endl;
		exit(-1);
	}
	
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;
	tmpWidth = getLambda();
	tmp = addTwoArmInterface(1.0, angle_of_incidence, 1.0, intIndex);	
	tmp->getTwoArmInterface()->setReflect(reflect); //Absorber has perfectly reflecting back boundary
	//tmp->getTwoArmInterface()->setReflect(0.0); //Absorber has perfectly reflecting back boundary
	
	if (numQW == 0)
	{
		// CASE 1: EMPTY Cavity
		tmpWidth = (getLambda()*(cavityLength1+cavityLength2))/cavityIndex;
		tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		
	} else if (numQW < 0)
	{
		// CASE 2: numQW<0
		//ABSORBING QW
		//First cavity
				
		tmpWidth = getLambda()*cavityLength1/cavityIndex;
		tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		//tmp->setOutputToFile(2);
		tmp->getTwoArmCavity()->setName("TACFRONT");
		
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add first device
			tmp = addTwoArmDevice();
			tmpName.str("");
			tmpName << "ABS1" << "_T" << j+1;
			tmp->getTwoArmDevice()->setName(tmpName.str());
			tmp->setOutputToFile(2);
			tmp->getTwoArmDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
		}
			tmpWidth = getLambda()*cavityLength2/cavityIndex;
			tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
			tmp->setOutputToFile(1);
			tmp->getTwoArmCavity()->setName("TACBACK");
	} else 
	{	
		// CASE 3: numQW > 0
		// Gain QW
		//First cavity
				
		tmpWidth = getLambda()*cavityLength1/cavityIndex;
		tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		//tmp->setOutputToFile(2);
		tmp->getTwoArmCavity()->setName("TACFRONT");
		
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add first device
			tmp = addTwoArmDevice();
			tmpName.str("");
			tmpName << "QW6" << "_T" << j+1;
			tmp->getTwoArmDevice()->setName(tmpName.str());
			tmp->setOutputToFile(2);
			tmp->getTwoArmDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
		}
			tmpWidth = getLambda()*cavityLength2/cavityIndex;
			tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
			tmp->setOutputToFile(1);
			tmp->getTwoArmCavity()->setName("TACBACK");
	}	
}

/*
 *  Create an absorber on the RIGHT
 * 
 *  |--AR--|--CAP--|---------------o-o-o-o---|
 * 
 *  Where:
 *  numABS			-> Number of ABS's in the medium [0,inf)
 * 						- If numABS = 0		=> Empty medium of length cavityLength
 * 	dx0_abs			-> Distance from RIGHTMOST ABS to the edge in units of lambda [0,1]
 *  dx_abs			-> Distance between each ABS in units of lambda (0,1]
 *  cavityLength	-> Total length of cavity (not including AR and CAP)
 * 						- If cavityLength = 0    => Fit cavity tight on ABS's (If numABS = 0 => ERROR)
 * 						- If cavityLength < space for ABS with spacing   => ERROR
 * 	cavityIndex		-> Refractive background index in medium
 * 	capIndex		-> Refractive index of CAP layer, if capIndex = 0 then NO CAP layer
 *  arIndex			-> Refractive index of AR coating, if arIndex = 0 then NO ar coating
 * */

void VECSEL::addRPG_ABS_RIGHT(int numABS, double dx0_abs, double dx_abs, double cavityLength, double cavityIndex, double capIndex, double arIndex)
{
	// Check the input for consistencty
	if (cavityLength < 0)
	{
		cout << "addRPG_ABS_RIGHT(): cavityLength < 0, quitting" << endl;
		exit(-1);
	} else if (dx0_abs < 0) {
		cout << "addRPG_ABS_RIGHT(): dx0_abs < 0, quitting" << endl;
		exit(-1);
	} else if (dx_abs < 0) {
		cout << "addRPG_ABS_RIGHT(): dx_abs <= 0, quitting" << endl;
		exit(-1);
	} else if (numABS < 0) {
		cout << "addRPG_ABS_RIGHT(): numABS < 0, quitting" << endl;
		exit(-1);
	}
	
	// Check a few logical inconsistencies
	if ((numABS == 0) && (cavityLength == 0))
	{
		cout << "addRPG_ABS_RIGHT(): numABS =0 && cavityLength = 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the ABS, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addRPG_ABS_RIGHT: No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}
	
	if ((numABS > 0)&&(cavityLength>0))
	{
		double tmpWidth = dx0_abs + ((double)numABS -1.0)*dx_abs;
		if (tmpWidth > cavityLength)
		{
			cout << "addRPG_ABS_RIGHT(): cavityLength too short, quitting" << endl;
			cout << "dx0_abs = " << dx0_abs << endl;
			cout << "cavityL = " << cavityLength << endl;
			exit(-1);
		}
		
		if ((dx_abs == 0)&&(numABS > 1))
		{
			cout << "addRPG_ABS_RIGHT(): cannot have numABS>0 && dx_qw = 0, quitting" << endl;
			exit(-1);
		}
	}
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;
/*	
	// AR 1
	tmpWidth = 0.262484*getLambda()/1.45;
	tmp = addCavity(tmpWidth, 1.45);
	tmp->getCavity()->setName("Ar");
	
	// AR 2
	tmpWidth = 0.132304*getLambda()/3.195400;
	tmp = addCavity(tmpWidth, 3.195400);
	tmp->getCavity()->setName("cap");
	
	// AR 1
	tmpWidth = 0.180177*getLambda()/2.963510;
	tmp = addCavity(tmpWidth, 2.963510);
	tmp->getCavity()->setName("Ar1");
	
	// AR 2
	tmpWidth = 0.149296*getLambda()/3.396376;
	tmp = addCavity(tmpWidth, 3.396376);
	tmp->getCavity()->setName("Ar2");

	// AR 1
	tmpWidth = 0.193327*getLambda()/2.963510;
	tmp = addCavity(tmpWidth, 2.963510);
	tmp->getCavity()->setName("Ar1");
	
	// AR 2
	tmpWidth = 0.291200*getLambda()/3.396376;
	tmp = addCavity(tmpWidth, 3.396376);
	tmp->getCavity()->setName("Ar2");

	// AR 1
	tmpWidth = 0.358399*getLambda()/2.963510;
	tmp = addCavity(tmpWidth, 2.963510);
	tmp->getCavity()->setName("Ar1");
	
	// AR 2
	tmpWidth = 0.222876*getLambda()/3.396376;
	tmp = addCavity(tmpWidth, 3.396376);
	tmp->getCavity()->setName("Ar2");

	// AR 1
	tmpWidth = 0.431467*getLambda()/2.963510;
	tmp = addCavity(tmpWidth, 2.963510);
	tmp->getCavity()->setName("Ar1");
*/


	// Add AR coating	
	if (arIndex > 0)
	{
		tmpWidth = (getLambda()/4.0)/arIndex;
		tmp = addCavity(tmpWidth, arIndex);
		tmp->getCavity()->setName("Ar");
	}


	// Add CAP layer
	if (capIndex > 0)
	{
		tmpWidth = (getLambda()/2.0)/capIndex;
		tmp = addCavity(tmpWidth, capIndex);
		tmp->getCavity()->setName("Cap");
	}
	
	if (numABS == 0)
	{
		// CASE 1: EMPTY Cavity
		tmpWidth = (getLambda()*cavityLength)/cavityIndex;
		tmp = addCavity(tmpWidth, cavityIndex);
		
	} else {
		
		// CASE 2: numABS > 0
		if (cavityLength > 0)
		{
			// Add first part of cavity
			tmpWidth = getLambda()*(cavityLength - dx0_abs - ((double)numABS -1.0)*dx_abs)/cavityIndex;
			tmp = addCavity( tmpWidth , cavityIndex);
		}

		// Add first device
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add device
			tmp = addDevice();
			tmpName.str("");
			tmpName << "ABS" << numABS << "_T" << j+1;
			tmp->getDevice()->setName(tmpName.str());
			//tmp->setOutputToFile(2);
			tmp->getDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
		}
		
		// Add rest
		double deviceDistance = (getLambda()*dx_abs)/cavityIndex;
		for(int i = numABS-1; i>=1; i--)
		{
			// Add Cavity
			tmp = addCavity(deviceDistance, cavityIndex);
			
			// Add device
			for(int j = 0; j < VECSEL_transverse_points_number; j++)
			{
				// Add device
				tmp = addDevice();
				tmpName.str("");
				tmpName << "ABS" << i << "_T" << j+1;
				tmp->getDevice()->setName(tmpName.str());
				//tmp->setOutputToFile(1);
				tmp->getDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
			}
		}
		
		// Add final cavity
		if (dx0_abs > 0)
		{
			tmpWidth = getLambda()*dx0_abs/cavityIndex;
			tmp = addCavity(tmpWidth, cavityIndex);
			tmp->getCavity()->setName("CBACK");
			//tmp->setOutputToFile(1);
		}
	}
}

/* Create a twoArm simple quantum well structure without an interface
 *
 * |---|--QW--|---|
 * */

void VECSEL::addTwoArmQW_noInterface(int numQW, double cavityLength1, double cavityLength2, double cavityIndex, double angle_of_incidence, double external_index)
{
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 		= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	double intIndex 	= cavityIndex;
	cavityIndex 		= cavityIndex*sqrt(1.0-a_sin2/(cavityIndex*cavityIndex));
	
	
	// Check the input for consistencty
	if (cavityLength1 < 0)
	{
		cout << "addTwoArmQW_STRUCT(): cavityLength1 < 0, quitting" << endl;
		exit(-1);
	} else if (cavityLength2 < 0)
	{
		cout << "addTwoArmQW_STRUCT(): cavityLength2 < 0, quitting" << endl;
		exit(-1);
	} 	
	// Check a few logical inconsistencies
	if ((numQW == 0) && (cavityLength1 == 0 || cavityLength2 ==0))
	{
		cout << "addTwoArmQW_STRUCT(): numQW =0 && cavityLength = 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the QW, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addTwoArmQW_STRUCT: No LEFT boundary detected. TwoArmStructure as first cavity is untested. Exiting" << endl;
		exit(-1);
	}
	
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;
	tmpWidth = getLambda();
	
	if (numQW == 0)
	{
		// CASE 1: EMPTY Cavity
		tmpWidth = (getLambda()*(cavityLength1+cavityLength2))/cavityIndex;
		tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		
	} else if (numQW < 0)
	{
		// CASE 2: numQW<0
		//ABSORBING QW
		//First cavity
				
		tmpWidth = getLambda()*cavityLength1/cavityIndex;
		tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		//tmp->setOutputToFile(2);
		tmp->getTwoArmCavity()->setName("TACFRONT");
		
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add first device
			tmp = addTwoArmDevice();
			tmpName.str("");
			tmpName << "ABS1" << "_T" << j+1;
			tmp->getTwoArmDevice()->setName(tmpName.str());
			tmp->setOutputToFile(2);
			tmp->getTwoArmDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
		}
			tmpWidth = getLambda()*cavityLength2/cavityIndex;
			tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
			tmp->setOutputToFile(1);
			tmp->getTwoArmCavity()->setName("TACBACK");
	} else 
	{	
		// CASE 3: numQW > 0
		// Gain QW. Use QW6 material file. Hard-coded
		//First cavity
				
		tmpWidth = getLambda()*cavityLength1/cavityIndex;
		tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		//tmp->setOutputToFile(2);
		tmp->getTwoArmCavity()->setName("TACFRONT");
		
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add first device
			tmp = addTwoArmDevice();
			tmpName.str("");
			tmpName << "QW6" << "_T" << j+1;
			tmp->getTwoArmDevice()->setName(tmpName.str());
			tmp->setOutputToFile(2);
			tmp->getTwoArmDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
		}
			tmpWidth = getLambda()*cavityLength2/cavityIndex;
			tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
			tmp->setOutputToFile(1);
			tmp->getTwoArmCavity()->setName("TACBACK");
	}	
}

/* Create a twoArm simple structure without an interface
 * */

void VECSEL::addTwoArm_Space(double cavityLength, double cavityIndex, double angle_of_incidence, double external_index, const std::string &name)
{	
	
	// Check the input for consistencty
	if (cavityLength < 0)
	{
		cout << "addTwoArm_Space(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the QW, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addTwoArm_Space: No LEFT boundary detected. TwoArmStructure as first cavity is untested. Exiting" << endl;
		exit(-1);
	}
	
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;
	tmpWidth = (getLambda()*cavityLength)/cavityIndex;
	tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
	tmp->setOutputToFile(1);
	tmp->getTwoArmCavity()->setName(name);
		
}

/* Create a birefringent crystal structure without an interface
 * */

void VECSEL::addBirefringentCrystal(double width, double cavityIndex, double extraordinary_n, double angle_of_incidence)
{	
	
	// Check the input for consistencty
	if (width < 0)
	{
		cout << "addTwoArmQW_STRUCT(): cavityLength1 < 0, quitting" << endl;
		exit(-1);
	}	
	// If there is nothing on the LEFT of the QW, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addTwoArmQW_STRUCT: No LEFT boundary detected. TwoArmStructure as first cavity is untested. Exiting" << endl;
		exit(-1);
	}
	
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;
	tmpWidth = (getLambda()*width)/cavityIndex;
	tmp = addBirefringentCrystal(tmpWidth, cavityIndex, extraordinary_n);
	tmp->setOutputToFile(1);
	tmp->getBirefringentCrystal()->setName("BRC");
		
}

/*
 *  Create a twoArm non-normal incidence absorber
 * 
 *  |--AR--|--CAP--|---------------o-o-o-o---|
 * 
 *  Where:
 *  numABS			-> Number of ABS's in the medium [0,inf)
 * 						- If numABS = 0		=> Empty medium of length cavityLength
 * 	dx0_abs			-> Distance from RIGHTMOST ABS to the edge in units of lambda [0,1]
 *  dx_abs			-> Distance between each ABS in units of lambda (0,1]
 *  cavityLength	-> Total length of cavity (not including AR and CAP)
 * 						- If cavityLength = 0    => Fit cavity tight on ABS's (If numABS = 0 => ERROR)
 * 						- If cavityLength < space for ABS with spacing   => ERROR
 * 	cavityIndex		-> Refractive background index in medium
 * 	capIndex		-> Refractive index of CAP layer, if capIndex = 0 then NO CAP layer
 *  arIndex			-> Refractive index of AR coating, if arIndex = 0 then NO ar coating
 * */

void VECSEL::addTwoArmRPG_ABS_LEFT(int numABS, double dx0_abs, double dx_abs, double cavityLength, double cavityIndex, double capIndex, double arIndex, double external_index, double angle_of_incidence)
{
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 		= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	double intIndex 	= cavityIndex;
	cavityIndex 		= cavityIndex*sqrt(1.0-a_sin2/(cavityIndex*cavityIndex));
	capIndex 		= capIndex*sqrt(1.0-a_sin2/(capIndex*capIndex));
	arIndex 		= arIndex*sqrt(1.0-a_sin2/(arIndex*arIndex));

	// Check the input for consistencty
	if (cavityLength < 0)
	{
		cout << "addTwoArmRPG_ABS(): cavityLength < 0, quitting" << endl;
		exit(-1);
	} else if (dx0_abs < 0) {
		cout << "addTwoArmRPG_ABS(): dx0_abs < 0, quitting" << endl;
		exit(-1);
	} else if (dx_abs < 0) {
		cout << "addTwoArmRPG_ABS(): dx_abs <= 0, quitting" << endl;
		exit(-1);
	} else if (numABS < 0) {
		cout << "addTwoArmRPG_ABS(): numABS < 0, quitting" << endl;
		exit(-1);
	}
	
	// Check a few logical inconsistencies
	if ((numABS == 0) && (cavityLength == 0))
	{
		cout << "addTwoArmRPG_ABS(): numABS =0 && cavityLength = 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the ABS, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addTwoArmRPG_ABS: No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}
	
	if ((numABS > 0)&&(cavityLength>0))
	{
		double tmpWidth = dx0_abs + ((double)numABS -1.0)*dx_abs;
		if (tmpWidth > cavityLength)
		{
			cout << "addTwoArmRPG_ABS(): cavityLength too short, quitting" << endl;
			cout << "dx0_abs = " << dx0_abs << endl;
			cout << "cavityL = " << cavityLength << endl;
			exit(-1);
		}
		
		if ((dx_abs == 0)&&(numABS > 1))
		{
			cout << "addTwoArmRPG_ABS(): cannot have numABS>0 && dx_qw = 0, quitting" << endl;
			exit(-1);
		}
	}
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;

	if (modules[startInd-1].isCavity())
	{
		tmp = addTwoArmInterface(1.0, angle_of_incidence, 1.0, intIndex);	
		tmp->getTwoArmInterface()->setReflect(-1.0); //Absorber has perfectly reflecting back boundary
		//tmp->getTwoArmInterface()->setReflect(0.0); //Absorber has perfectly reflecting back boundary
	}


	if (numABS == 0)
	{
		// CASE 1: EMPTY Cavity
		tmpWidth = (getLambda()*cavityLength)/cavityIndex;
		tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		
	} else {
		// CASE 2: numABS > 0
		
		// Add final cavity
		if (dx0_abs > 0)
		{
			tmpWidth = getLambda()*dx0_abs/cavityIndex;
			tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
			//tmp->getTwoArmCavity()->setName("TACBACK");
			//tmp->setOutputToFile(1);
		}
		
		// Add rest
		double deviceDistance = (getLambda()*dx_abs)/cavityIndex;
		for(int i = numABS-1; i>=1; i--)
		{
			// Add device
			for(int j = 0; j < VECSEL_transverse_points_number; j++)
			{
				// Add device
				tmp = addTwoArmDevice();
				tmpName.str("");
				tmpName << "ABS" << i << "_T" << j+1;
				tmp->getTwoArmDevice()->setName(tmpName.str());
				//WTFFF-tmp->setOutputToFile(1);
				tmp->getTwoArmDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
			}
			// Add Cavity
			tmp = addTwoArmCavity(deviceDistance, cavityIndex, angle_of_incidence, external_index);
		}

		// Add first device
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add device
			tmp = addTwoArmDevice();
			tmpName.str("");
			tmpName << "ABS" << numABS << "_T" << j+1;
			tmp->getTwoArmDevice()->setName(tmpName.str());
			//tmp->setOutputToFile(0);
			tmp->getTwoArmDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
		}
		
		if (cavityLength > 0)
		{
			// Add first part of cavity
			tmpWidth = getLambda()*(cavityLength - dx0_abs - ((double)numABS -1.0)*dx_abs)/cavityIndex;
			tmp = addTwoArmCavity( tmpWidth , cavityIndex, angle_of_incidence, external_index);
		}
		
	}
	

	// Add CAP layer
	if (capIndex > 0)
	{
		tmpWidth = (getLambda()/2.0)/capIndex;
		tmp = addTwoArmCavity(tmpWidth, capIndex, angle_of_incidence, external_index);
		tmp->getTwoArmCavity()->setName("Cap");
	}
	
	// Add AR coating	
	if (arIndex > 0)
	{
		tmpWidth = (getLambda()/4.0)/arIndex;
		tmp = addTwoArmCavity(tmpWidth, arIndex, angle_of_incidence, external_index);
		tmp->getTwoArmCavity()->setName("Ar");
	}

}

/*
 *  Create a twoArm non-normal incidence absorber
 * 
 *  |--AR--|--CAP--|---------------o-o-o-o---|
 * 
 *  Where:
 *  numABS			-> Number of ABS's in the medium [0,inf)
 * 						- If numABS = 0		=> Empty medium of length cavityLength
 * 	dx0_abs			-> Distance from RIGHTMOST ABS to the edge in units of lambda [0,1]
 *  dx_abs			-> Distance between each ABS in units of lambda (0,1]
 *  cavityLength	-> Total length of cavity (not including AR and CAP)
 * 						- If cavityLength = 0    => Fit cavity tight on ABS's (If numABS = 0 => ERROR)
 * 						- If cavityLength < space for ABS with spacing   => ERROR
 * 	cavityIndex		-> Refractive background index in medium
 * 	capIndex		-> Refractive index of CAP layer, if capIndex = 0 then NO CAP layer
 *  arIndex			-> Refractive index of AR coating, if arIndex = 0 then NO ar coating
 * */

void VECSEL::addTwoArmRPG_ABS(int numABS, double dx0_abs, double dx_abs, double cavityLength, double cavityIndex, double capIndex, double arIndex, double external_index, double angle_of_incidence)
{
	//==========================
	// With angle of incidence
	//==========================
	double a_sin2 		= external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence);
	double intIndex 	= cavityIndex;
	cavityIndex 		= cavityIndex*sqrt(1.0-a_sin2/(cavityIndex*cavityIndex));
	capIndex 		= capIndex*sqrt(1.0-a_sin2/(capIndex*capIndex));
	arIndex 		= arIndex*sqrt(1.0-a_sin2/(arIndex*arIndex));

	// Check the input for consistencty
	if (cavityLength < 0)
	{
		cout << "addTwoArmRPG_ABS(): cavityLength < 0, quitting" << endl;
		exit(-1);
	} else if (dx0_abs < 0) {
		cout << "addTwoArmRPG_ABS(): dx0_abs < 0, quitting" << endl;
		exit(-1);
	} else if (dx_abs < 0) {
		cout << "addTwoArmRPG_ABS(): dx_abs <= 0, quitting" << endl;
		exit(-1);
	} else if (numABS < 0) {
		cout << "addTwoArmRPG_ABS(): numABS < 0, quitting" << endl;
		exit(-1);
	}
	
	// Check a few logical inconsistencies
	if ((numABS == 0) && (cavityLength == 0))
	{
		cout << "addTwoArmRPG_ABS(): numABS =0 && cavityLength = 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the ABS, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addTwoArmRPG_ABS: No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}
	
	if ((numABS > 0)&&(cavityLength>0))
	{
		double tmpWidth = dx0_abs + ((double)numABS -1.0)*dx_abs;
		if (tmpWidth > cavityLength)
		{
			cout << "addTwoArmRPG_ABS(): cavityLength too short, quitting" << endl;
			cout << "dx0_abs = " << dx0_abs << endl;
			cout << "cavityL = " << cavityLength << endl;
			exit(-1);
		}
		
		if ((dx_abs == 0)&&(numABS > 1))
		{
			cout << "addTwoArmRPG_ABS(): cannot have numABS>0 && dx_qw = 0, quitting" << endl;
			exit(-1);
		}
	}
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;

	if (modules[startInd-1].isCavity())
	{
		tmp = addTwoArmInterface(1.0, angle_of_incidence, 1.0, intIndex);	
		tmp->getTwoArmInterface()->setReflect(-1.0); //Absorber has perfectly reflecting back boundary
		//tmp->getTwoArmInterface()->setReflect(0.0); //Absorber has perfectly reflecting back boundary
	}

/*	// AR 1
	tmpWidth = 0.146560*(getLambda())/arIndex;
	tmp = addTwoArmCavity(tmpWidth, arIndex, angle_of_incidence, external_index);
	tmp->getTwoArmCavity()->setName("SiO2");
*/

	// Add AR coating	
	if (arIndex > 0)
	{
		tmpWidth = (getLambda()/4.0)/arIndex;
		tmp = addTwoArmCavity(tmpWidth, arIndex, angle_of_incidence, external_index);
		tmp->getTwoArmCavity()->setName("Ar");
	}


	// Add CAP layer
	if (capIndex > 0)
	{
		tmpWidth = (getLambda()/2.0)/capIndex;
		tmp = addTwoArmCavity(tmpWidth, capIndex, angle_of_incidence, external_index);
		tmp->getTwoArmCavity()->setName("Cap");
	}
	
	if (numABS == 0)
	{
		// CASE 1: EMPTY Cavity
		tmpWidth = (getLambda()*cavityLength)/cavityIndex;
		tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
		
	} else {
		
		// CASE 2: numABS > 0
		if (cavityLength > 0)
		{
			// Add first part of cavity
			tmpWidth = getLambda()*(cavityLength - dx0_abs - ((double)numABS -1.0)*dx_abs)/cavityIndex;
			tmp = addTwoArmCavity( tmpWidth , cavityIndex, angle_of_incidence, external_index);
		}

		// Add first device
		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			// Add device
			tmp = addTwoArmDevice();
			tmpName.str("");
			tmpName << "ABS" << numABS << "_T" << j+1;
			tmp->getTwoArmDevice()->setName(tmpName.str());
			//tmp->setOutputToFile(0);
			tmp->getTwoArmDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
		}
		
		// Add rest
		double deviceDistance = (getLambda()*dx_abs)/cavityIndex;
		for(int i = numABS-1; i>=1; i--)
		{
			// Add Cavity
			tmp = addTwoArmCavity(deviceDistance, cavityIndex, angle_of_incidence, external_index);
			
			// Add device
			for(int j = 0; j < VECSEL_transverse_points_number; j++)
			{
				// Add device
				tmp = addTwoArmDevice();
				tmpName.str("");
				tmpName << "ABS" << i << "_T" << j+1;
				tmp->getTwoArmDevice()->setName(tmpName.str());
				//WTFFF-tmp->setOutputToFile(1);
				tmp->getTwoArmDevice()->setTransversePosition(VECSEL_transverse_points_y[j]);
			}
		}
		
		// Add final cavity
		if (dx0_abs > 0)
		{
			tmpWidth = getLambda()*dx0_abs/cavityIndex;
			tmp = addTwoArmCavity(tmpWidth, cavityIndex, angle_of_incidence, external_index);
			//tmp->getTwoArmCavity()->setName("TACBACK");
			//tmp->setOutputToFile(1);
		}
	}
}


void VECSEL::addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(double arIndex, double arLength, double dx0, double dx1, double cavityIndex, double DBR_n1, double DBR_n1_length, double DBR_n2, double DBR_n2_length, int num_dbr_layers, double PHASE_n, double PHASE_length, double angle_of_incidence, double external_index)
{
	// Check the input for consistencty
	if (arIndex <= 1.0)
	{
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): arIndex <= 1.0, quitting" << endl;
		exit(-1);

	} else if (arLength < 0.0) {
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): arLength < 0.0, quitting" << endl;
		exit(-1);

	} else if (dx0 < 0.0) {
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): dx0 < 0.0, quitting" << endl;
		exit(-1);

	} else if (dx1 < 0.0) {
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): dx1 < 0.0, quitting" << endl;
		exit(-1);

	} else if (cavityIndex <= 1.0) {
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): cavityIndex <= 1.0, quitting" << endl;
		exit(-1);
	} else if  (DBR_n1 <= 1.0) {
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): DBR_n1 <= 1.0, quitting" << endl;
		exit(-1);

	} else if (DBR_n1_length <= 0.0) {
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): DBR_n1_length <= 0.0, quitting" << endl;
		exit(-1);

	} else if (DBR_n2 <= 1.0) {
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): DBR_n2 <= 1.0, quitting" << endl;
		exit(-1);

	} else if (DBR_n2_length <= 0.0) {
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): DBR_n2_length <= 0.0, quitting" << endl;
		exit(-1);

	} else if (num_dbr_layers <= 0) {
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): num_dbr_layers <= 0.0, quitting" << endl;
		exit(-1);

	} else if (PHASE_n <= 1.0) {
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): PHASE_n <= 1.0, quitting" << endl;
		exit(-1);
	
	} else if (PHASE_length <= 0.0) {
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): PHASE_length <= 0.0, quitting" << endl;
		exit(-1);

	}

	// Illogical statements
	if (DBR_n1 == DBR_n2)
	{
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): DBR layers cannot be equal, quitting" << endl;
		exit(-1);
	}

	
	// If there is nothing on the LEFT of the ABS, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(): No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}
	
	Module *tmp;
	double tmpWidth;
	std::stringstream tmpName;



	double ar1_l = 112.60*nm; // 111.1nm
	double ar1_n = 1.4770;
	// AR 1
	tmp = addCavity(ar1_l, ar1_n, angle_of_incidence, external_index);
	tmp->getCavity()->setName("Ar1");

	int index_initial_cavity = quick_index_cavity.size()-1; // Store initial cavity (AR1)

	// AR 2
	double ar2_l = 113.54*nm; //114.07nm
	double ar2_n = 1.9826;
	tmp = addCavity(ar2_l, ar2_n, angle_of_incidence, external_index);
	tmp->getCavity()->setName("Ar2");

	// Surface layer
	tmp = addCavity(dx0 , cavityIndex, angle_of_incidence, external_index);
	
	// Add device
	tmp = addDevice();	
	tmpName.str("");
	tmpName << "ABS1";
	tmp->getDevice()->setName(tmpName.str());
	tmp->setOutputToFile(1);
	
	// Add final cavity
	tmp = addCavity(dx1 , cavityIndex, angle_of_incidence, external_index);

	// ADD DBR
	// Ensure alternating reflectors
	double n_bragg[2] = {DBR_n1,DBR_n2};
	double L_bragg[2] = {DBR_n1_length,DBR_n2_length};
	
	// Create numLayers cavities of given length
	for(int i = 0; i < num_dbr_layers; i++)
	{
		tmp = addCavity(L_bragg[i%2], n_bragg[i%2], angle_of_incidence, external_index);
	}


	// ADD PHASE LAYER
	tmp = addCavity(PHASE_length , PHASE_n, angle_of_incidence, external_index);


	
	//====================================================================
	// Add in angle dependence for ALL layers in this function
	// First cavity
	std::complex<double> n_curr = modules[quick_index_cavity[index_initial_cavity]].getRefInd() + I*modules[quick_index_cavity[index_initial_cavity]].getRefInd_im();
	std::complex<double> n_prev = external_index;

	// Convert to standard form
	double n2_tilde = real(n_prev)*real(n_curr) + imag(n_prev)*imag(n_curr);
	double k2_tilde = real(n_prev)*imag(n_curr) - real(n_curr)*imag(n_prev);
	double n1_tilde = abs(n_prev)*abs(n_prev);

	// Compute temporary terms
	double n1_sin_th_2 = n1_tilde*sin(angle_of_incidence); // na*sin(theta) 
	n1_sin_th_2 *= n1_sin_th_2; // ^2
	double norm_n_2 = n2_tilde*n2_tilde + k2_tilde*k2_tilde;
	double term_a = 1.0+n1_sin_th_2/norm_n_2;

	// Kovalenko 2001: Descartes-Snell law of refraction with absorption
	// compute sin(theta)^2
	double sin_th_2 = 0.5*(term_a - sqrt(term_a*term_a - 4.0*n2_tilde*n2_tilde*n1_sin_th_2/(norm_n_2*norm_n_2)));
	double cos_th = sqrt(1.0-sin_th_2);
	modules[quick_index_cavity[index_initial_cavity]].setCosTh(cos_th,cos_th);

	// Iterate over all cavities
	for(int i=index_initial_cavity+1; i<quick_index_cavity.size(); i++)
	{
		n_curr = modules[quick_index_cavity[i]].getRefInd() + I*modules[quick_index_cavity[i]].getRefInd_im();
		n_prev = modules[quick_index_cavity[i-1]].getRefInd() + I*modules[quick_index_cavity[i-1]].getRefInd_im();

		// Convert to standard form
		n2_tilde = real(n_prev)*real(n_curr) + imag(n_prev)*imag(n_curr);
		k2_tilde = real(n_prev)*imag(n_curr) - real(n_curr)*imag(n_prev);
		n1_tilde = abs(n_prev)*abs(n_prev);

		// Compute temporary terms
		n1_sin_th_2 = n1_tilde*n1_tilde*sin_th_2; // (n1 sin(th))^2
		norm_n_2 = n2_tilde*n2_tilde + k2_tilde*k2_tilde;
		term_a = 1.0+n1_sin_th_2/norm_n_2;

		// Kovalenko 2001: Descartes-Snell law of refraction with absorption
		// compute sin(theta)^2
		sin_th_2 = 0.5*(term_a - sqrt(term_a*term_a - 4.0*n2_tilde*n2_tilde*n1_sin_th_2/(norm_n_2*norm_n_2)));
		cos_th = sqrt(1.0-sin_th_2);
		modules[quick_index_cavity[i]].setCosTh(cos_th,cos_th);
	}


}


/*
 * Add an empty space with a seperator in the middle for output
 *  -> probe_tune is the tuning of space. Range(probe_tine) = [-1,1] in units of lambda
 * */

void VECSEL::addCUSTOM_SPACE(double cavityLength, double leftDeviceCavityLength, double cavityIndex, double probe_tune, const std::string &name)
{
	if (cavityLength < 0)
	{
		cout << "addRPG_SPACE(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the SPACE, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addRPG_SPACE: No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}
	
	if ((probe_tune > 1.0)||(probe_tune < -1.0))
	{
			cout << "addRPG_SPACE: Probe is out of bounds, should be in [-1,1] but is = " << probe_tune << endl;
			exit(-1);
	}
	
	double length1 = getLambda()*(cavityLength/2.0 + probe_tune)/cavityIndex - leftDeviceCavityLength;
	double length2 = getLambda()*(cavityLength)/cavityIndex - length1;
	
	Module *tmp;
	
	tmp = addCavity(length1, cavityIndex);
	tmp->setOutputToFile(1);
	tmp->getCavity()->setName(name);
	
	tmp = addCavity(length2, cavityIndex);
}

/*
 * Add an empty space with a seperator in the middle for output
 *  -> probe_tune is the tuning of space. Range(probe_tine) = [-1,1] in units of lambda
 * */

void VECSEL::addRPG_SPACE(double cavityLength, double cavityIndex, double probe_tune, const std::string &name)
{
	if (cavityLength < 0)
	{
		cout << "addRPG_SPACE(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the SPACE, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		//cout << "addRPG_SPACE: No LEFT boundary detected, adding perfect reflection" << endl;
		//addBoundary(1.0,1.0);
	}
	
	if ((probe_tune > 1.0)||(probe_tune < -1.0))
	{
			cout << "addRPG_SPACE: Probe is out of bounds, should be in [-1,1] but is = " << probe_tune << endl;
			exit(-1);
	}
	
	double length1 = cavityLength/2.0 + probe_tune;
	double length2 = cavityLength - length1;
	
	Module *tmp;
	
	tmp = addCavity((getLambda()*length1)/cavityIndex, cavityIndex);
	tmp->setOutputToFile(1);
	tmp->getCavity()->setName(name);
	
	tmp = addCavity((getLambda()*length2)/cavityIndex, cavityIndex);
}


/*
 * Add an ETALON + AIR AROUND with a seperator in the middle for output
 *  -> probe_tune is the tuning of space. Range(probe_tine) = [-1,1] in units of lambda
 * */

void VECSEL::addRPG_SPACE_ETALON(double cavityLength, double cavityIndex, double probe_tune, const std::string &name, double ETALON_INDEX, double ETALON_WIDTH)
{
	if (cavityLength < 0)
	{
		cout << "addRPG_SPACE_ETALON(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the SPACE, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		//cout << "addRPG_SPACE: No LEFT boundary detected, adding perfect reflection" << endl;
		//addBoundary(1.0,1.0);
	}
	
	if ((probe_tune > 1.0)||(probe_tune < -1.0))
	{
			cout << "addRPG_SPACE_ETALON: Probe is out of bounds, should be in [-1,1] but is = " << probe_tune << endl;
			exit(-1);
	}

	double length1 = (cavityLength-0.5)/2.0 + probe_tune;
	
	Module *tmp;
	
	tmp = addCavity((getLambda()*length1)/cavityIndex, cavityIndex);
	tmp->setOutputToFile(1);
	tmp->getCavity()->setName(name);

	tmp = addCavity((getLambda()*length1)/cavityIndex, cavityIndex);

	tmp = addCavity(ETALON_WIDTH, ETALON_INDEX);
	tmp->getCavity()->setName("ETALON1");
	
	tmp = addCavity((getLambda()*0.5)/cavityIndex, cavityIndex);
}
void VECSEL::addRPG_SPACE_ETALON2(double cavityLength, double cavityIndex, double probe_tune, const std::string &name, double ETALON_INDEX1, double ETALON_WIDTH1, double ETALON_INDEX2, double ETALON_WIDTH2)
{
	if (cavityLength < 0)
	{
		cout << "addRPG_SPACE_ETALON(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the SPACE, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		//cout << "addRPG_SPACE: No LEFT boundary detected, adding perfect reflection" << endl;
		//addBoundary(1.0,1.0);
	}
	
	if ((probe_tune > 1.0)||(probe_tune < -1.0))
	{
			cout << "addRPG_SPACE_ETALON: Probe is out of bounds, should be in [-1,1] but is = " << probe_tune << endl;
			exit(-1);
	}

	Module *tmp;
	
	tmp = addCavity((getLambda()*cavityLength)/cavityIndex, cavityIndex);
	tmp->setOutputToFile(1);
	tmp->getCavity()->setName(name);

	if (ETALON_WIDTH1 < ETALON_WIDTH2)
	{
		cout<< "VECSEL::addRPG_SPACE_ETALON2(): Are you sure you're designing the correct cavity etalon ORDER???" << endl;
		exit(-1);
	}

	tmp = addCavity(ETALON_WIDTH1, ETALON_INDEX1);
	tmp->getCavity()->setName("ETALON1");
	tmp = addCavity((0.25*getLambda())/cavityIndex, cavityIndex);

	tmp = addCavity(ETALON_WIDTH2, ETALON_INDEX2);
	tmp->getCavity()->setName("ETALON2");
	tmp = addCavity((0.25*getLambda())/cavityIndex, cavityIndex);
}

// wa_s -> is wa*dt
// wb_s -> wb*dt
// w0_s -> w0*dt
// width_s -> spectral_width*dt
// filter_length is number of timesteps in filter

void VECSEL::addRPG_SPACE_FILTER(double cavityLength, double cavityIndex, double probe_tune, const std::string &name, double wa_s, double wb_s, double w0_s, double width_s, int filter_length)
{
	if (cavityLength < 0)
	{
		cout << "addRPG_SPACE_FILTER(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the SPACE, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		//cout << "addRPG_SPACE: No LEFT boundary detected, adding perfect reflection" << endl;
		//addBoundary(1.0,1.0);
	}
	
	if ((probe_tune > 1.0)||(probe_tune < -1.0))
	{
			cout << "addRPG_SPACE_FILTER: Probe is out of bounds, should be in [-1,1] but is = " << probe_tune << endl;
			exit(-1);
	}
	
	double length1 = cavityLength/2.0 + probe_tune;
	double length2 = cavityLength - length1 - 10.0;
	double length3 = 10.0;
	
	Module *tmp;
	
	tmp = addCavity((getLambda()*(0.1*length1))/cavityIndex, cavityIndex);
	tmp = addFilter(filter_length);
	//tmp->getFilter()->setdoubleExp_flatTopWindow(wa_s, wb_s);
	tmp->getFilter()->setFilter_pluss_doubleGauss_flatTopWindow(wa_s, wb_s, width_s, w0_s); // Gauss filter
	//tmp->getFilter()->setFilter_pluss_diagnostic_square(0.5);
	//tmp->getFilter()->setFilter_pluss_PassThroughDelay(); // Delay "filter"
	tmp->getFilter()->setFilter_minus_PassThroughDelay(); // Delay "filter"

	tmp = addCavity((getLambda()*(0.9*length1))/cavityIndex, cavityIndex);
	tmp->setOutputToFile(1);
	tmp->getCavity()->setName(name);
	
	tmp = addCavity((getLambda()*length2)/cavityIndex, cavityIndex);

	tmp = addCavity((getLambda()*length3)/cavityIndex, cavityIndex); // Added layer to make the filter work better
}

void VECSEL::addRPG_SPACE_ANGLE(double cavityLength, double cavityIndex, double probe_tune, const std::string &name, double angle_of_incidence, double angle_of_incidence_transfer)
{
	if (cavityLength < 0)
	{
		cout << "addRPG_SPACE(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the SPACE, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		//cout << "addRPG_SPACE: No LEFT boundary detected, adding perfect reflection" << endl;
		//addBoundary(1.0,1.0);
	}
	
	if ((probe_tune > 1.0)||(probe_tune < -1.0))
	{
			cout << "addRPG_SPACE: Probe is out of bounds, should be in [-1,1] but is = " << probe_tune << endl;
			exit(-1);
	}
	
	double length1 = cavityLength/2.0 + probe_tune; 
	double length2 = cavityLength - length1;
	
	Module *tmp;

	if (angle_of_incidence != angle_of_incidence_transfer)
	{
	
		tmp = addCavity(0.01*(getLambda()*length1)/cavityIndex, cavityIndex, angle_of_incidence, cavityIndex);
		double cos_left, cos_right; tmp->getCosTh(&cos_left,&cos_right);
		tmp->setCosTh(cos_left,1.0);

		// Change left angle
		tmp = addCavity(0.99*(getLambda()*length1)/cavityIndex, cavityIndex, angle_of_incidence_transfer, cavityIndex);
	} else {
		tmp = addCavity((getLambda()*length1)/cavityIndex, cavityIndex, angle_of_incidence, cavityIndex);
	}
	tmp->setOutputToFile(2);
	tmp->getCavity()->setName(name);
	
	tmp = addCavity((getLambda()*length2)/cavityIndex, cavityIndex, angle_of_incidence_transfer, cavityIndex);
}


void VECSEL::addRPG_SPACE_ANGLE_LEFT(double cavityLength, double cavityIndex, double probe_tune, const std::string &name, double angle_of_incidence, double angle_of_incidence_transfer)
{
	if (cavityLength < 0)
	{
		cout << "addRPG_SPACE(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the SPACE, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		//cout << "addRPG_SPACE: No LEFT boundary detected, adding perfect reflection" << endl;
		//addBoundary(1.0,1.0);
	}
	
	if ((probe_tune > 1.0)||(probe_tune < -1.0))
	{
			cout << "addRPG_SPACE: Probe is out of bounds, should be in [-1,1] but is = " << probe_tune << endl;
			exit(-1);
	}
	
	double length2 = cavityLength/2.0 + probe_tune; 
	double length1 = cavityLength - length2;
	
	Module *tmp;

	if (angle_of_incidence != angle_of_incidence_transfer)
	{
	
		tmp = addCavity(0.01*(getLambda()*length1)/cavityIndex, cavityIndex, angle_of_incidence, cavityIndex);
		double cos_left, cos_right; tmp->getCosTh(&cos_left,&cos_right);
		tmp->setCosTh(cos_left,1.0);

		// Change left angle
		tmp = addCavity(0.99*(getLambda()*length1)/cavityIndex, cavityIndex, angle_of_incidence_transfer, cavityIndex);
	} else {
		tmp = addCavity((getLambda()*length1)/cavityIndex, cavityIndex, angle_of_incidence, cavityIndex);
	}
	
	tmp = addCavity((getLambda()*length2)/cavityIndex, cavityIndex, angle_of_incidence_transfer, cavityIndex);
	tmp->setOutputToFile(1);
	tmp->getCavity()->setName(name);
}

void VECSEL::addRPG_SPACE_LENS(double cavityLength, double cavityIndex, double lens_1, double armLength, double lens_2)
{
	if (cavityLength < 0)
	{
		cout << "addRPG_SPACE_LENS(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	Module *tmp;
	tmp = addCavity((getLambda()*cavityLength)/cavityIndex, cavityIndex, 0.0, cavityIndex);
	tmp->getCavity()->setName("BPMHC");
	tmp->getCavity()->set_equivalent_optical_cavity(lens_1, armLength, lens_2);
}

void VECSEL::addRPG_SPACE_LENS_ANGLE(double cavityLength, double cavityIndex, double lens_1, double armLength, double lens_2, double angle_of_incidence)
{
	if (cavityLength < 0)
	{
		cout << "addRPG_SPACE_LENS(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	Module *tmp;
	tmp = addCavity((getLambda()*cavityLength)/cavityIndex, cavityIndex, angle_of_incidence, cavityIndex);
	tmp->getCavity()->setName("BPMHC");
	tmp->getCavity()->set_equivalent_optical_cavity(lens_1, armLength, lens_2);
}

void VECSEL::addRPG_SPACE_BPM(double cavityLength, double opticalLength, double cavityIndex, double probe_tune, const std::string &name)
{
	if (cavityLength < 0)
	{
		cout << "addRPG_SPACE(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the SPACE, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addRPG_SPACE: No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}
	
	if ((probe_tune > 1.0)||(probe_tune < -1.0))
	{
			cout << "addRPG_SPACE: Probe is out of bounds, should be in [-1,1] but is = " << probe_tune << endl;
			exit(-1);
	}
	
	double length1 = cavityLength/2.0 + probe_tune;
	double length2 = cavityLength - length1;
	
	Module *tmp;
	
	tmp = addCavity((getLambda()*length1)/cavityIndex, cavityIndex);
	//tmp->setOutputToFile(2); 
	//tmp->getCavity()->setName("BPMFS2");
	tmp->getCavity()->setName(name);
	
	tmp = addCavity((getLambda()*length2)/cavityIndex, cavityIndex);
	tmp->getCavity()->setName("BPMFS");
	tmp->getCavity()->set_equivalent_optical_cavity(0.0, opticalLength, 0.0);
}

void VECSEL::addRPG_SPACE_BPM_ANGLE(double cavityLength, double opticalLength, double cavityIndex, double probe_tune, const std::string &name, double angle_of_incidence)
{
	if (cavityLength < 0)
	{
		cout << "addRPG_SPACE(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the SPACE, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addRPG_SPACE: No LEFT boundary detected, adding perfect reflection" << endl;
		addBoundary(1.0,1.0);
	}
	
	if ((probe_tune > 1.0)||(probe_tune < -1.0))
	{
			cout << "addRPG_SPACE: Probe is out of bounds, should be in [-1,1] but is = " << probe_tune << endl;
			exit(-1);
	}
	
	double length1 = cavityLength/2.0 + probe_tune;
	double length2 = cavityLength - length1;
	
	Module *tmp;
	
	tmp = addCavity((getLambda()*length1)/cavityIndex, cavityIndex, angle_of_incidence, cavityIndex);
	//tmp->setOutputToFile(2); 
	//tmp->getCavity()->setName("BPMFS2");
	tmp->getCavity()->setName(name);
	
	tmp = addCavity((getLambda()*length2)/cavityIndex, cavityIndex, angle_of_incidence, cavityIndex);
	tmp->getCavity()->setName("BPMFS");
	tmp->getCavity()->set_equivalent_optical_cavity(0.0, opticalLength, 0.0);	
}


void VECSEL::addRPG_SPACE_LENS_transform(double cavityLength, double cavityIndex, double A, double B, double C, double D)
{
	if (cavityLength < 0)
	{
		cout << "addRPG_SPACE(): cavityLength < 0, quitting" << endl;
		exit(-1);
	}
	
	// If there is nothing on the LEFT of the SPACE, add a boundary for consistency
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		//cout << "addRPG_SPACE: No LEFT boundary detected, adding perfect reflection" << endl;
		//addBoundary(1.0,1.0);
	}

	// Construct Lengths/Focus from ABCD matrix formulation
	double LENS_1 = B/(1.0-D);
	double LENS_2 = B/(1.0-A);
	double LENGTH = B;
	
	Module *tmp;

	tmp = addCavity((getLambda()*cavityLength)/cavityIndex, cavityIndex, 0.0, cavityIndex);
	tmp->getCavity()->setName("BPMTRANS");
	tmp->getCavity()->set_equivalent_optical_cavity(LENS_1, LENGTH, LENS_2);
}


void VECSEL::addQW_SPACE_ABS(int numQW, int numABS, double dx0_qw, double dx0_abs, double cavityLength, double cavityIndex)
{
	// Check the input
	if ((numQW < 0)||(numABS < 0))
	{
		cout << "addQW_SPACE_ABS: Need numQW, numABS >= 0" << endl;
		exit(-1);
	}
	if ((dx0_qw < 0)||(dx0_abs < 0))
	{
		cout << "addQW_SPACE_ABS: Need dx0_qw, dx0_abs >= 0" << endl;
		exit(-1);
	}
	if (cavityLength <= 0)
	{
		cout << "addQW_SPACE_ABS: Need cavityLength > 0" << endl;
		exit(-1);
	}
	
	// Find position of previous element
	int startInd = getNumberModules();
	if (startInd == 0)
	{
		cout << "addQW_SPACE_ABS: Adding perfect absorbtion on the left" << endl;
		addBoundary(0.0,1.0);
		startInd += 1;
	}
	
	// Store starting position
	double startPos = modules.back().getPosition1();
	
	// Start adding cavity / devices
	double deviceDistance = getLambda()/(2.0*cavityIndex);
	double tmpWidth = 0;
	
	// Start adding cavity
	Module *tmp;
	std::stringstream tmpName;
	if ((numQW == 0)&&( numABS == 0))
	{
		// CASE 1: EMPTY Cavity
		tmpWidth = getLambda()*cavityLength;
		tmp = addCavity(tmpWidth, cavityIndex);
		
	} else if ((numQW > 0)&&( numABS == 0)) {
	
		// Check if too many QW's
		if (getLambda()*cavityLength - ((double)numQW-1.0)*deviceDistance - dx0_qw < 0.0)
		{
			cout << "addQW_SPACE_ABS: Cannot fit all QW's into cavity. Either add length or reduce numQW" << endl;
			exit(-1);
		}
		
		// CASE 2: No absorbers
		
		// First cavity
		tmpWidth = getLambda()*dx0_qw;
		tmp = addCavity(tmpWidth, cavityIndex);
		
		// Add first device
		tmp = addDevice();	
		tmp->getDevice()->setName("QW 1");

		// Add rest in lambda/2 distance
		for(int i = 1; i < numQW; i++)
		{
			// Add Cavity
			tmp = addCavity(deviceDistance, cavityIndex);
			
			// Add device
			tmp = addDevice();
			tmpName.str("");
			tmpName << "QW " << i+1;
			tmp->getDevice()->setName(tmpName.str());
		}
		
		// Add rest of Cavity
		tmpWidth = getLambda()*(cavityLength - modules.back().getPosition1());
		tmp = addCavity( tmpWidth , cavityIndex);

	} else if ((numQW == 0)&&( numABS > 0))
	{
		// Check if too many ABS's
		if (getLambda()*cavityLength - ((double)numABS-1.0)*deviceDistance - dx0_abs < 0.0)
		{
			cout << "addQW_SPACE_ABS: Cannot fit all ABS's into cavity. Either add length or reduce numABS" << endl;
			exit(-1);
		}
		// CASE 3: No QW
		
		// First cavity
		tmpWidth = getLambda()*cavityLength - ( getLambda()*dx0_abs + deviceDistance*((double)numABS-1.0) );
		tmp = addCavity(tmpWidth, cavityIndex);
		
		// Add first device
		tmp = addDevice();	
		tmp->getDevice()->setName("ABS 1");
		
		
		// Add rest in lambda/2 distance
		for(int i = 1; i < numABS; i++)
		{
			// Add Cavity
			tmp = addCavity(deviceDistance, cavityIndex);
			
			// Add device
			tmp = addDevice();
			tmpName.str("");
			tmpName << "ABS " << i+1;
			tmp->getDevice()->setName(tmpName.str());
		}
		
		// Add final cavity
		tmpWidth = getLambda()*dx0_abs;
		tmp = addCavity(tmpWidth, cavityIndex);
		
	} else {
		
		// Check if too many ABS's
		if (getLambda()*cavityLength - ((double)numQW-1.0)*deviceDistance - dx0_qw - ((double)numABS-1.0)*deviceDistance - dx0_abs < 0.0)
		{
			cout << "addQW_SPACE_ABS: Cannot fit all QW's + ABS's into cavity. Either add length or reduce numQW and numABS" << endl;
			exit(-1);
		}
	
		// CASE 4: Both QW and ABS
		
		// First cavity
		tmpWidth = getLambda()*dx0_qw;
		tmp = addCavity(tmpWidth, cavityIndex);
		
		// Add 1st QW
		tmp = addDevice();	
		tmp->getDevice()->setName("QW 1");
		
		// Add rest in lambda/2 distance
		for(int i = 1; i < numQW; i++)
		{
			// Add Cavity
			tmp = addCavity(deviceDistance, cavityIndex);

			// Add device
			tmp = addDevice();
			tmpName.str("");
			tmpName << "QW " << i+1;
			tmp->getDevice()->setName(tmpName.str());
		}
		
		// Add middle cavity
		double middle_space = getLambda()*cavityLength - getLambda()*dx0_qw - deviceDistance*((double)numQW-1.0) - getLambda()*dx0_abs - deviceDistance*((double)numABS-1.0);
		
		tmp = addCavity(middle_space,cavityIndex);
		
		// Add first device
		tmp = addDevice();	
		tmp->getDevice()->setName("ABS 1");
		
		// Add ABSORBERS
		for(int i = 1; i < numABS; i++)
		{
			// Add Cavity
			tmp = addCavity(deviceDistance, cavityIndex);

			// Add device
			tmp = addDevice();
			tmpName.str("");
			tmpName << "ABS " << i+1;
			tmp->getDevice()->setName(tmpName.str());
		}
		
		// Add final cavity
		tmpWidth = getLambda()*dx0_abs;
		tmp = addCavity(tmpWidth, cavityIndex);
	}
	
	cout << "addQW_SPACE_ABS: Adding REFLECTING boundary on the right" << endl;
	addBoundary(1.0,1.0);
}

// Move QWs 1 to num to new positions given in newPos[0] -> newPos[num-1]
// Should only be used in a linear cavity
// Not set up for twoArmDevices
void VECSEL::moveQWpositions(int num, double *newPos, double DT)
{
	
	// Find target cavities
	int target_cav[num+1];
	for(int i = 0; i < num; i++)
	{
		target_cav[i] = quick_index_device_previous_cavity[i];
	}
	
	// Find next cavity after the last QW
	for(int i = target_cav[num-1]+1; i < modules.size(); i++)
	{
		if (modules[i].isCavity())
		{
			target_cav[num] = i;
			break;
		}
	}
	
	
/*	
	// TODO: REMOVE fake placement
	for(int i = 0; i < num; i++)
	{
		double oldPos = modules[target_cav[i]+1].getPosition1();
		newPos[i] = oldPos - 10.0*nm;
	}

	cout << "Cav positions.. before" << endl;
	for(int i = 0; i < num+1; i++)
	{
		double x0 = modules[target_cav[i]].getPosition0();
		double x1 = modules[target_cav[i]].getPosition1();
		cout << "[" << target_cav[i] << "]: x = [ " << x0/1.0e-6 << " / " << x1/1.0e-6 << "] [um] width = " << (x1-x0)/1.0e-9 << " [nm]" << endl;
	}
*/
	
	// check that new points are consistent
	double x0 = modules[target_cav[0  ]].getPosition0();
	double x1 = modules[target_cav[num]].getPosition1();
	
	for(int i = 0; i < num; i++)
	{
		if ((newPos[i] <= x0)||(newPos[i] >= x1))
		{
			cout << "moveQWpositions():: Positions out of bound" << endl;
			cout << "Boundaries = [ " << x0 << " / " << x1 << " ]" << endl;
			cout << "Position:" << endl;
			for(int j = 0; j < num; j++)
			{
				cout << "[" << j << "]: " << newPos[j] << endl;
			}
			exit(-1);
		}
	}
	
	
	// Set new limits for cavities, keep boundaries fixed
	modules[target_cav[0]].setPosition(modules[target_cav[0]].getPosition0(),newPos[0]);	// Cav
	modules[target_cav[0]+1].setPosition(newPos[0],newPos[0]); // Device
	
	for(int i = 1; i < num; i++)
	{
		int indx = target_cav[i];
		modules[target_cav[i]].setPosition(modules[target_cav[i]-1].getPosition1(),newPos[i]); // Cav
		modules[target_cav[i]+1].setPosition(newPos[i],newPos[i]); // Device
	}
	modules[target_cav[num]].setPosition(newPos[num-1],modules[target_cav[num]].getPosition1());	// Cav
	
/*	
	cout << "Cav positions.. after" << endl;
	for(int i = 0; i < num+1; i++)
	{
		double x0 = modules[target_cav[i]].getPosition0();
		double x1 = modules[target_cav[i]].getPosition1();
		cout << "[" << target_cav[i] << "]: x = [ " << x0/1.0e-6 << " / " << x1/1.0e-6 << "] [um] width = " << (x1-x0)/1.0e-9 << " [nm]" << endl;
	}
*/
	
	// Reset new storage with the initialize function
	for(int i =0; i < num+1; i++)
	{
		modules[target_cav[i]].getCavity()->initializeZero(DT,VECSEL_transverse_points_number, VECSEL_transverse_points_y, VECSEL_transverse_points_R_max, VECSEL_transverse_points_boundary_guard_ratio);
	}
	//cout << "sucsess??" << endl;
	//exit(-1);
}

void VECSEL::getQWpositions(int num, double **newPos)
{	
	// Export QW positnions
	for(int i = 0; i < num; i++)
	{
		(*newPos)[i] =  modules[quick_index_device[i]].getPosition0();
	}
}

/*
	Add EMPTY CAVITY 
*/

Module * VECSEL::addCavity()
{	
	// Update number
	setNumberCavities(getNumberCavities() + 1);
	
	// Update Quick refrence
	quick_index_cavity.push_back(getNumberModules());


	// Check if previous cavity is end of two arm structure	
	if(modules.back().isTwoArmCavity())
	{
		quick_index_twoArmPostCav.push_back(getNumberModules());
		for(unsigned i=0; i<getNumberModules()-1; i++)
		{
			if(modules[getNumberModules()-1-i].isTwoArmInterface())
			{
				modules[getNumberModules()-1-i].getTwoArmInterface()->setPostCav(getNumberModules());
				break;
			}
		}
	}
	// Add new module
	Module *tmp = addModule();
	
	// Add Cavity
	modules.back().addCavity();
	modules.back().setToFileOutputKey(getToFileOutputKey());
	return tmp;
}

/*
	Copy structure into our currenct structure
*/
void VECSEL::appendStructure(VECSEL *baseStruct)
{
	double offset;
	if (getNumberModules()>0)
	{
		offset = modules.back().getPosition1();
	} else {
		offset = 0.0;
	}

	// Set transverse data to old cavity
	set_transverse_dimensions(baseStruct->getTransversePoints(), baseStruct->getTransverseRmax(), baseStruct->getTransverseBoundaryGuardRatio());

	
	for(int i =0; i < baseStruct->getNumberModules(); i++)
	{
		// Create a copy
		//Module *tmp = addModule();
		Module *base = baseStruct->getModule(i);
		
		if (base->isDevice())
		{
			std::stringstream oldName;
			oldName << base->getDevice()->getName();
			
			Module *tmp = addDevice();
			tmp->getDevice()->setName(oldName.str());
			tmp->setOutputToFile(1);

			double pos_y = base->getDevice()->getTransversePosition();
			tmp->getDevice()->setTransversePosition(pos_y);
			
		} else if (base->isCavity())
		{
			double width = base->getWidth();
			double n1    = base->getCavity()->getRefInd();
			double n2    = base->getCavity()->getRefInd_im();
			std::complex<double> nc = n1 + I*n2;
			double cos_th_left, cos_th_right;
			base->getCavity()->getCosTh(&cos_th_left, &cos_th_right);

			addCavity_cos(width, nc, cos_th_left, cos_th_right);
			
		} else if (base->isBoundary())
		{
			// Boundary info
			double ref = base->getBoundary()->getRefCoeff();
			double nextN = base->getBoundary()->getNextCavityIndex();
			addBoundary(ref,nextN);
			
		} else if (base->isLossElement())
		{
			// Boundary info
			double loss_plus  = 0;
			double loss_minus = 0;
			base->getLossElement()->getLossCoeff(&loss_minus, &loss_plus);
			addBoundary(loss_minus,loss_plus);
			
		} else {
				cout << "VECSEL::appendStructure(): Unknown module. Not configured for VCAV" << endl;
				exit(-1);
		}
		
		double x0 = modules.back().getPosition0();
		double x1 = modules.back().getPosition1();
		modules.back().setPosition(x0+offset,x1+offset);
	}

	// Set all QWs with the given transverse gain profile
	set_transverse_QW_pump_profile_SuperGaussian(baseStruct->getSGPumpDegree(), baseStruct->getSGPumpFWHM());
	//set_transverse_QW_temp_profile_SuperGaussian(baseStruct->getSGTempDegree(), baseStruct->getSGTempFWHM());
	set_transverse_QW_temp_profile(baseStruct->VECSEL_initial_temp_profile_T);
}

/*
	Add TwoArmCAVITY with width and refractive index n1
	Position is automatricly set from previous device
*/

Module * VECSEL::addTwoArmCavity(double width, std::complex<double> n1)
{
	
	if (width<=0)
	{
		cout<<"addTwoArmCavity: Cannot have width<=0. Quitting"<<endl;
		exit(-1);
	}
	if (abs(n1)<=0)
	{
		cout<<"addTwoArmCavity: Cannot have abs(n1)<=0. Quitting"<<endl;
		exit(-1);
	}
	// Update number
	setNumberTwoArmCavities(getNumberTwoArmCavities() + 1);
	
	// Update Quick refrence
	quick_index_twoArmCavity.push_back(getNumberModules());

	// Add new module
	Module *tmp = addModule();

	// Add Cavity
	modules.back().addTwoArmCavity();
	
	double pos = 0;
	if (getNumberModules() > 1) // It self + others
	{
		// Set position to previous device
		pos = modules[getNumberModules()-2].getPosition1();
	}
	
	modules.back().setRefInd(real(n1));
	modules.back().setRefInd_im(imag(n1));
	modules.back().setWidth(width);
	modules.back().setToFileOutputKey(getToFileOutputKey());
	
	modules.back().setPosition(pos, pos + width);
	return tmp;
}

Module * VECSEL::addTwoArmInterface(std::complex<double> n1, double angle_of_incidence, double external_index, double interference_index)
{
	if (abs(n1)<=0)
	{
		cout<<"addTwoArmInterface: Cannot have abs(n1)<=0. Quitting"<<endl;
		exit(-1);
	}
	if (external_index<=0)
	{
		cout<<"addTwoArmInterface: Cannot have external_index<=0. Quitting"<<endl;
		exit(-1);
	}
	if (interference_index<=0)
	{
		cout<<"addTwoArmInterface: Cannot have interferenceIndex<=0. Quitting"<<endl;
		exit(-1);
	}

	if (angle_of_incidence<0)
	{
		cout<<"addTwoArmInterface: Negative angle of incidence untested. Quitting."<<endl;
		exit(-1);
	}
	

	int indx=getNumberModules();
	if (!modules[indx-1].isCavity())
	{
		cout<<"addTwoArmInterface: Previous modules must be a cavity. Coincident TwoArmStructs untested. Quitting."<<endl;
		exit(-1);
	}
	// Snell's Law. Assumes real refractive indices
	double cos_th = sqrt(1.0-external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence)/(real(n1)*real(n1)));

	// Update number
	setNumberTwoArmInterfaces(getNumberTwoArmInterfaces() + 1);

	// Update Quick refrence
	quick_index_twoArmInterface.push_back(getNumberModules());
	
	// Add new module
	Module *tmp = addModule();

	// Add Two Arm Interface
	modules.back().addTwoArmInterface();
	
	
	double pos = 0;
	if (getNumberModules() > 1) // It self + others
	{
		// Set position to previous device
		pos = modules[getNumberModules()-2].getPosition1();
	}

	double width = getLambda();
	#ifdef TRANS_DELAY
		if(angle_of_incidence==0)
		{
			cout<<"addTwoArmInterface: Cannot use transverse delay with zero angle of incidence. Quitting."<<endl;
			exit(-1);
		} else	
		{
			width = tan(angle_of_incidence)*(VECSEL_transverse_points_y[VECSEL_transverse_points_number-1]-VECSEL_transverse_points_y[0]);
		}
	#endif

	double angle_of_interference=asin(sin(angle_of_incidence)*external_index/interference_index); //Snell's law	

	modules.back().setRefInd(real(n1));
	modules.back().setRefInd_im(imag(n1));
	modules.back().setCosTh(cos_th,cos_th);
	modules.back().setWidth(width);
	modules.back().setAngle(angle_of_incidence); //Angle at interface
	modules.back().setIntAngle(angle_of_interference); //Angle for interference fringes (scaled by QW refractive index)
	modules.back().setPrevCav(getNumberModules()-1);
	modules.back().setToFileOutputKey(getToFileOutputKey());	
	modules.back().setPosition(pos, pos + width);
	return tmp;
}

Module * VECSEL::addBirefringentCrystal(double width, std::complex<double> n1, double extraordinary_n)
{
	if (width<=0)
	{
		cout<<"addBirefringentCrystal: Cannot have width<=0. Quitting"<<endl;
		exit(-1);
	}
	if (abs(n1)<=0)
	{
		cout<<"addBirefringentCrystal: Cannot have abs(n1)<=0. Quitting"<<endl;
		exit(-1);
	}
	if (getNumberTwoArmCavities()==0)
	{
		cout<<"addBirefringentCrystal: BRC as first element is untested. Quitting."<<endl;
		exit(-1);
	}

	// Update number
	setNumberBirefringentCrystals(getNumberBirefringentCrystals() + 1);
	
	// Update Quick refrence
	quick_index_birefringentCrystal.push_back(getNumberModules());
	
	// Add new module
	Module *tmp = addModule();
	
	// Add Cavity
	modules.back().addBirefringentCrystal();
	
	
	double pos = 0;
	if (getNumberModules() > 1) // It self + others
	{
		// Set position to previous device
		pos = modules[getNumberModules()-2].getPosition1();
	}
	
	modules.back().setRefInd(real(n1));
	modules.back().setRefInd_extraAxis(extraordinary_n);
	modules.back().setRefInd_im(imag(n1));
	modules.back().setWidth(width);
	modules.back().setToFileOutputKey(getToFileOutputKey());
	
	modules.back().setPosition(pos, pos + width);
	return tmp;
}
		

Module * VECSEL::addTwoArmCavity(double width, std::complex<double> n1, double angle_of_incidence, double external_index)
{
	if (width<=0)
	{
		cout<<"addTwoArmCavity: Cannot have width<=0. Quitting"<<endl;
		exit(-1);
	}
	if (abs(n1)<=0)
	{
		cout<<"addTwoArmCavity: Cannot have abs(n1)<=0. Quitting"<<endl;
		exit(-1);
	}
	if (external_index<=0)
	{
		cout<<"addTwoArmCavity: Cannot have external_index<=0. Quitting"<<endl;
		exit(-1);
	}
	if (angle_of_incidence<0)
	{
		cout<<"addTwoArmCavity: Negative angle of incidence untested. Quitting."<<endl;
		exit(-1);
	}
	// Snell's Law. Assumes real refractive indices
	double cos_th = sqrt(1.0-external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence)/(real(n1)*real(n1)));
	//width = width/cos_th;

	// Update number
	setNumberTwoArmCavities(getNumberTwoArmCavities() + 1);
	
	// Update Quick refrence
	quick_index_twoArmCavity.push_back(getNumberModules());
	
	// Add new module
	Module *tmp = addModule();
	
	// Add Cavity
	modules.back().addTwoArmCavity();
	
	
	double pos = 0;
	if (getNumberModules() > 1) // It self + others
	{
		// Set position to previous device
		pos = modules[getNumberModules()-2].getPosition1();
	}
	
	
	modules.back().setRefInd(real(n1));
	modules.back().setRefInd_im(imag(n1));
	modules.back().setWidth(width);
	modules.back().setCosTh(cos_th,cos_th);
	modules.back().setToFileOutputKey(getToFileOutputKey());	
	modules.back().setPosition(pos, pos + width);

	return tmp;
}

Module * VECSEL::addTwoArmCavity_cos(double width, std::complex<double> n1, double cos_th_left, double cos_th_right)
{
	if (width<=0)
	{
		cout<<"addTwoArmCavity: Cannot have width<=0. Quitting"<<endl;
		exit(-1);
	}
	if (abs(n1)<=0)
	{
		cout<<"addTwoArmCavity: Cannot have abs(n1)<=0. Quitting"<<endl;
		exit(-1);
	}
	
	// Update number
	setNumberTwoArmCavities(getNumberTwoArmCavities() + 1);
	
	// Update Quick refrence
	quick_index_twoArmCavity.push_back(getNumberModules());
	
	// Add new module
	Module *tmp = addModule();
	
	// Add Cavity
	modules.back().addTwoArmCavity();
	
	
	double pos = 0;
	if (getNumberModules() > 1) // It self + others
	{
		// Set position to previous device
		pos = modules[getNumberModules()-2].getPosition1();
	}
	
	
	modules.back().setRefInd(real(n1));
	modules.back().setRefInd_im(imag(n1));
	modules.back().setWidth(width);
	modules.back().setCosTh(cos_th_left, cos_th_right);
	modules.back().setToFileOutputKey(getToFileOutputKey());
	
	modules.back().setPosition(pos, pos + width);
	return tmp;
}

/*
	Add CAVITY with width and refractive index n1
	Position is automatricly set from previous device
*/

Module * VECSEL::addCavity(double width, std::complex<double> n1)
{
	
	// Update number
	setNumberCavities(getNumberCavities() + 1);
	
	// Update Quick refrence
	quick_index_cavity.push_back(getNumberModules());
	
	// Check if previous cavity is end of two arm structure	
	if(modules.back().isTwoArmCavity())
	{
		quick_index_twoArmPostCav.push_back(getNumberModules());
		for(unsigned i=0; i<getNumberModules()-1; i++)
		{
			if(modules[getNumberModules()-1-i].isTwoArmInterface())
			{
				modules[getNumberModules()-1-i].getTwoArmInterface()->setPostCav(getNumberModules());
				break;
			}
		}
	}

	// Add new module
	Module *tmp = addModule();

	// Add Cavity
	modules.back().addCavity();
	
	double pos = 0;
	if (getNumberModules() > 1) // It self + others
	{
		// Set position to previous device
		pos = modules[getNumberModules()-2].getPosition1();
	}
	
	modules.back().setRefInd(real(n1));
	modules.back().setRefInd_im(imag(n1));
	modules.back().setWidth(width);
	modules.back().setToFileOutputKey(getToFileOutputKey());
	
	modules.back().setPosition(pos, pos + width);
	return tmp;
}

Module * VECSEL::addCavity(double width, std::complex<double> n1, double angle_of_incidence, double external_index)
{
	// Snell's Law. Assumes real refractive indices
	double cos_th = sqrt(1.0-external_index*external_index*sin(angle_of_incidence)*sin(angle_of_incidence)/(real(n1)*real(n1)));
	//width = width/cos_th;

	// Update number
	setNumberCavities(getNumberCavities() + 1);
	
	// Update Quick refrence
	quick_index_cavity.push_back(getNumberModules());
	
	// Check if previous cavity is end of two arm structure	
	if(modules.back().isTwoArmCavity())
	{
		quick_index_twoArmPostCav.push_back(getNumberModules());
		for(unsigned i=0; i<getNumberModules()-1; i++)
		{
			if(modules[getNumberModules()-1-i].isTwoArmInterface())
			{
				modules[getNumberModules()-1-i].getTwoArmInterface()->setPostCav(getNumberModules());
				break;
			}
		}
	}
	
	// Add new module
	Module *tmp = addModule();
	
	// Add Cavity
	modules.back().addCavity();
	
	
	double pos = 0;
	if (getNumberModules() > 1) // It self + others
	{
		// Set position to previous device
		pos = modules[getNumberModules()-2].getPosition1();
	}
	
	
	modules.back().setRefInd(real(n1));
	modules.back().setRefInd_im(imag(n1));
	modules.back().setWidth(width);
	modules.back().setCosTh(cos_th,cos_th);
	modules.back().setToFileOutputKey(getToFileOutputKey());
	
	modules.back().setPosition(pos, pos + width);
	return tmp;
}

Module * VECSEL::addCavity_cos(double width, std::complex<double> n1, double cos_th_left, double cos_th_right)
{
	// Update number
	setNumberCavities(getNumberCavities() + 1);
	
	// Update Quick refrence
	quick_index_cavity.push_back(getNumberModules());
	
	// Check if previous cavity is end of two arm structure	
	if(modules.back().isTwoArmCavity())
	{
		quick_index_twoArmPostCav.push_back(getNumberModules());
		for(unsigned i=0; i<getNumberModules()-1; i++)
		{
			if(modules[getNumberModules()-1-i].isTwoArmInterface())
			{
				modules[getNumberModules()-1-i].getTwoArmInterface()->setPostCav(getNumberModules());
				break;
			}
		}
	}
	
	// Add new module
	Module *tmp = addModule();
	
	// Add Cavity
	modules.back().addCavity();
	
	
	double pos = 0;
	if (getNumberModules() > 1) // It self + others
	{
		// Set position to previous device
		pos = modules[getNumberModules()-2].getPosition1();
	}
	
	
	modules.back().setRefInd(real(n1));
	modules.back().setRefInd_im(imag(n1));
	modules.back().setWidth(width);
	modules.back().setCosTh(cos_th_left, cos_th_right);
	modules.back().setToFileOutputKey(getToFileOutputKey());
	
	modules.back().setPosition(pos, pos + width);
	return tmp;
}

Module * VECSEL::addDevice()
{
	// Update Number of devices
	setNumberDevices(getNumberDevices() + 1);
	
	// Update Quick refrence
	quick_index_device.push_back(getNumberModules());
	
	// Update Quick refrence
	quick_index_totalDevice.push_back(getNumberModules());
	
	// Update Quick refrence 
	int prev_cav_index;
	for(int i = getNumberModules()-1; i>=0; i--)
	{
		if (modules[i].isCavity())
		{
			prev_cav_index = i;
			break;
		}
	}
	quick_index_device_previous_cavity.push_back(prev_cav_index);

	// Add new module
	Module *tmp = addModule();

	// Add device
	modules.back().addDevice();
	
	if (getNumberModules() <= 1)
	{
		cout << "VECSEL::addDevice(): Need boundary & cavity on the far left of VECSEL" << endl;
		exit(-1);
		
	} else if (getNumberModules() > 1)
	{
		// Set position to previous device
		double pos = modules[getNumberModules()-2].getPosition1();
		
		modules.back().setPosition(pos, pos);
		modules.back().setToFileOutputKey(getToFileOutputKey());
	}
	return tmp;
}

Module * VECSEL::addTwoArmDevice()
{
	// Update Number of devices
	setNumberTwoArmDevices(getNumberTwoArmDevices() + 1);
	
	// Update Quick refrence
	quick_index_twoArmDevice.push_back(getNumberModules());
	
	// Update Quick refrence
	quick_index_totalDevice.push_back(getNumberModules());
	
	// Update Quick refrence 
	int prev_cav_index;
	for(int i = getNumberModules()-1; i>=0; i--)
	{
		if (modules[i].isTwoArmCavity())
		{
			prev_cav_index = i;
			break;
		}
	}
	quick_index_device_previous_cavity.push_back(prev_cav_index);
	
	// Add new module
	Module *tmp = addModule();

	// Add device
	modules.back().addTwoArmDevice();
	
	if (getNumberModules() <= 1)
	{
		cout << "VECSEL::addTwoArmDevice(): Need boundary & cavity on the far left of VECSEL" << endl;
		exit(-1);
		
	} else if (getNumberModules() > 1)
	{
		// Set position to previous device
		double pos = modules[getNumberModules()-2].getPosition1();
		
		modules.back().setPosition(pos, pos);
		modules.back().setToFileOutputKey(getToFileOutputKey());
	}
	return tmp;
}

Module * VECSEL::addBoundary()
{
	// Update number
	setNumberBoundaries(getNumberBoundaries() + 1);
	
	// Update Quick refrence
	quick_index_boundary.push_back(getNumberModules());
	// Add new module
	Module *tmp = addModule();
	
	// Add Cavity
	modules.back().addBoundary();
	
	if (getNumberModules() > 1)
	{
		// Set position to previous device
		double pos = modules[getNumberModules()-2].getPosition1();
		modules.back().setPosition(pos, pos);
	}
	return tmp;
}

Module * VECSEL::addBoundary(double refCoeff, double next_cavity_index)
{
	// Update number
	setNumberBoundaries(getNumberBoundaries() + 1);
	
	// Update Quick refrence
	quick_index_boundary.push_back(getNumberModules());
	
	// Add new module
	Module *tmp = addModule();
	
	// Add Cavity
	modules.back().addBoundary(refCoeff);
	modules.back().getBoundary()->setNextCavityIndex(next_cavity_index);
	
	if (getNumberModules() > 1)
	{
		// Set position to previous device
		double pos = modules[getNumberModules()-2].getPosition1();
		modules.back().setPosition(pos, pos);
	}
	return tmp;
}

Module * VECSEL::addLossElement()
{
	// Update number
	//setNumberBoundaries(getNumberBoundaries() + 1);
	
	// Update Quick refrence
	//quick_index_boundary.push_back(getNumberModules());
	
	// Add new module
	Module *tmp = addModule();
	
	// Add Cavity
	modules.back().addLossElement();
	
	if (getNumberModules() > 1)
	{
		// Set position to previous device
		double pos = modules[getNumberModules()-2].getPosition1();
		modules.back().setPosition(pos, pos);
	}
	return tmp;
}

Module * VECSEL::addLossElement(double loss_minus, double loss_plus)
{
	// Update number
	//setNumberBoundaries(getNumberBoundaries() + 1);
	
	// Update Quick refrence
	//quick_index_boundary.push_back(getNumberModules());
	
	// Add new module
	Module *tmp = addModule();
	
	// Add Cavity
	modules.back().addLossElement(loss_minus, loss_plus);
	
	if (getNumberModules() > 1)
	{
		// Set position to previous device
		double pos = modules[getNumberModules()-2].getPosition1();
		modules.back().setPosition(pos, pos);
	}
	return tmp;
}

Module * VECSEL::addFilter(int numEl)
{
	// Update number
	//setNumberBoundaries(getNumberBoundaries() + 1);
	
	// Update Quick refrence
	//quick_index_boundary.push_back(getNumberModules());
	
	// Add new module
	Module *tmp = addModule();
	
	// Add Cavity
	modules.back().addFilter(numEl);
	
	if (getNumberModules() > 1)
	{
		// Set position to previous device
		double pos = modules[getNumberModules()-2].getPosition1();
		modules.back().setPosition(pos, pos);
	}
	return tmp;
}

Module * VECSEL::addModule()
{
	//Module *tmp = new Module(); // Initialize to NULL
	//modules.push_back(*tmp);
	modules.push_back(Module()); // Temporary object is copied in
	return &modules.back();
}



//=======================================
// File IO functions and their helpers

/* Remove all output files from devices */
void VECSEL::file_output_device_disable()
{
    #ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		modules[indx].setOutputToFile(0);
	}
	#endif
    
}

/* Reduce output files from devices by turning off output outside a given region */
void VECSEL::file_output_device_reduce(double x_max)
{
    #ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		double trans_x;
		if (modules[indx].isDevice())
		{
			trans_x = modules[indx].getDevice()->getTransversePosition();
		} else
		{
			trans_x = modules[indx].getTwoArmDevice()->getTransversePosition();
		}

		if (abs(trans_x) < x_max/2.0)
		{
			modules[indx].setOutputToFileLevel(1);
			modules[indx].setOutputToFile(1);
		} else {
			modules[indx].setOutputToFile(0); // No output outside of range
		}
	}
	#endif   
}

/* Set output frequency for devices y = 0, y = +/- dx, y +/- 2*dx, ... inside the rangle [-x_max/2, x_max/2] */
void VECSEL::file_output_device_set_freq(double output_frequency, double x_max, int twoArmDevice_OutputLevel, int device_OutputLevel)
{
	if (output_frequency <= 0.0)
	{
		cout << "VECSEL::file_output_device_set_freq() Error. Cannot have 'output_frequency' <= 0" << endl;
		cout << "output_frequency = " << output_frequency << endl;
		exit(-1);
	
	} else if (x_max <= 0.0)
	{
		cout << "VECSEL::file_output_device_set_freq() Error. Cannot have 'x_max' <= 0" << endl;
		cout << "x_max = " << x_max << endl;
		exit(-1);

	} else if (x_max >= VECSEL_transverse_points_R_max)
	{
		cout << "VECSEL::file_output_device_set_freq() Error. Cannot have 'x_max' > domain size" << endl;
		cout << "x_max = " << x_max << endl;
		cout << "domain size = " << VECSEL_transverse_points_R_max << endl;
		exit(-1);

	} else if ((output_frequency < VECSEL_transverse_points_y[1]-VECSEL_transverse_points_y[0]) && VECSEL_transverse_points_number>1)
	{
		cout << "VECSEL::file_output_device_set_freq() Error. Cannot have 'output_frequency' faster than resolution" << endl;
		cout << "output_frequency = " << output_frequency << endl;
		cout << "domain resolution = " << VECSEL_transverse_points_y[1]-VECSEL_transverse_points_y[0] << endl;
		exit(-1);

	} else if  (0.5*x_max <= output_frequency)
	{
		cout << "VECSEL::file_output_device_set_freq() Error. Cannot have 'output_frequency' longer than '0.5*x_max'" << endl;
		cout << "0.5*x_max = " << 0.5*x_max << endl;
		cout << "output_frequency = " << output_frequency << endl;
		exit(-1);
	}

    #ifdef ITERATE_QW
	// Turn off all output
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		if (modules[indx].isDevice())
		{
			double trans_x = modules[indx].getDevice()->getTransversePosition();
		} else
		{
			double trans_x = modules[indx].getTwoArmDevice()->getTransversePosition();
		}
		modules[indx].setOutputToFile(0);
	}

	// Turn on output only for given devices
	int MAX_OUTPUT = floor(0.5*x_max/output_frequency);

	VECSEL_transverse_points_device_number = 1+2*(MAX_OUTPUT-1);
	if (VECSEL_transverse_points_device_y == NULL)
	{
		VECSEL_transverse_points_device_y = new double[VECSEL_transverse_points_device_number];
	}
	for (int i  = 0; i < MAX_OUTPUT; i++)
	{
		double output_x_p = +i*output_frequency;
		int min_indx_p, min_indx_m;
		double min_v;
		if (abs(output_x_p) <= x_max/2.0)
		{
			min_indx_p = -1;
			min_v = 1e99;

			for(int j = 0; j < VECSEL_transverse_points_number; j++)
			{
				if (abs(VECSEL_transverse_points_y[j]-output_x_p) < min_v)
				{
					min_indx_p = j;
					min_v = abs(VECSEL_transverse_points_y[j]-output_x_p);
				}
			}

			if (min_indx_p >= 0)
			{
				for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
				{
					int indx = quick_index_totalDevice[j]; // Index of device
					double trans_x;
					if (modules[indx].isDevice())
					{
						trans_x = modules[indx].getDevice()->getTransversePosition();
						if (trans_x-VECSEL_transverse_points_y[min_indx_p] == 0)
						{
							// Set output to file
							modules[indx].setOutputToFile(device_OutputLevel);
							modules[indx].setOutputToFileLevel(device_OutputLevel);
						}
					} else
					{
						trans_x = modules[indx].getTwoArmDevice()->getTransversePosition();
						if (trans_x-VECSEL_transverse_points_y[min_indx_p] == 0)
						{
							// Set output to file
							modules[indx].setOutputToFile(twoArmDevice_OutputLevel);
							modules[indx].setOutputToFileLevel(twoArmDevice_OutputLevel);
						}
					}
				}
			} else {
				cout << "VECSEL::file_output_device_set_freq() Cannot find a QW in range [0, x_max/2] closer than max" << endl;
				exit(-1);
			}
		}

		double output_x_m = -i*output_frequency;
		if (abs(output_x_m) <= x_max/2.0)
		{
			min_indx_m = -1;
			min_v = 1e99;

			for(int j = 0; j < VECSEL_transverse_points_number; j++)
			{
				if (abs(VECSEL_transverse_points_y[j]-output_x_m) < min_v)
				{
					min_indx_m = j;
					min_v = abs(VECSEL_transverse_points_y[j]-output_x_m);
				}
			}


			if (min_indx_m >= 0)
			{
				for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
				{
					int indx = quick_index_totalDevice[j]; // Index of device
					double trans_x;	
					if (modules[indx].isDevice())
					{
						trans_x = modules[indx].getDevice()->getTransversePosition();
						if (trans_x-VECSEL_transverse_points_y[min_indx_m] == 0)
						{
							// Set output to file
							modules[indx].setOutputToFile(device_OutputLevel);
							modules[indx].setOutputToFileLevel(device_OutputLevel);
						}
					} else
					{
						trans_x = modules[indx].getTwoArmDevice()->getTransversePosition();
						if (trans_x-VECSEL_transverse_points_y[min_indx_m] == 0)
						{
							// Set output to file
							modules[indx].setOutputToFile(twoArmDevice_OutputLevel);
							modules[indx].setOutputToFileLevel(twoArmDevice_OutputLevel);
						}
					}
				}
			} else {
				cout << "VECSEL::file_output_device_set_freq() Cannot find a QW in range [-x_max/2, 0] closer than max" << endl;
				exit(-1);
			}
		}

		if (i==0)
		{
			VECSEL_transverse_points_device_y[MAX_OUTPUT-1] = VECSEL_transverse_points_y[min_indx_p];
		} else{
			if (abs(output_x_p) <= x_max/2.0)
			{
				VECSEL_transverse_points_device_y[MAX_OUTPUT-1+i] = VECSEL_transverse_points_y[min_indx_p];
			}
			if (abs(output_x_m) <= x_max/2.0)
			{
				VECSEL_transverse_points_device_y[MAX_OUTPUT-1-i] = VECSEL_transverse_points_y[min_indx_m];
			}
		}
	}

	std::stringstream fileName;
	fileName.str("");
	fileName << getToFileOutputKey() << "transverse_grid_device_y.dat";
	saveBinary(fileName.str(), VECSEL_transverse_points_device_y, VECSEL_transverse_points_device_number);
	

	#endif
}

/* Write the given cavity structure to file
 * Each row represents a module
 * 
 * Output format:
 * 1st column: What type are we talking about
 *  -> 0 is boundary
 *  -> 1 is cavity
 *  -> 2 is device
 *  -> 3 is everything else
 * 2nd column: Position of first vertex
 * 3rd column: Position of second vertex
 * 4th column: Refractive index if a cavity. else -1
 * */
void VECSEL::file_output_structure(void)
{
	if (MPI_MY_RANK==0)
	{
		std::ofstream printStruc;
		
		std::cout.precision(5);
		
		std::stringstream fileName;
		fileName << getToFileOutputKey() << "system_structure.dat";
		printStruc.open((fileName.str()).c_str(),std::ofstream::out);
			
		for(unsigned i = 0; i < getNumberModules(); i++)
		{		
			if (modules[i].isBoundary())
			{
				printStruc << std::fixed << std::setprecision(16) << 0 << " " << modules[i].getPosition0() << " " << modules[i].getPosition1() << " " << "-1" << endl;
			}
			if (modules[i].isCavity())
			{
				printStruc << std::fixed << std::setprecision(16) << 1 << " " << modules[i].getPosition0() << " " << modules[i].getPosition1() << " " << modules[i].getRefInd() << endl;
			}
			if (modules[i].isTwoArmInterface())
			{
				printStruc << std::fixed << std::setprecision(16) << 1 << " " << modules[i].getPosition0() << " " << modules[i].getPosition1() << " " << modules[i].getRefInd() << endl;
			}
			if (modules[i].isTwoArmCavity())
			{
				printStruc << std::fixed << std::setprecision(16) << 1 << " " << modules[i].getPosition0() << " " << modules[i].getPosition1() << " " << modules[i].getRefInd() << endl;
			}
			if (modules[i].isDevice())
			{
				int indx0 = i;
				int indx1 = i;
				while (modules[indx0].isCavity())
				{
					indx0--;
					if (indx0 < 0)
					{
						cout << "VECSEL::file_output_structure() cannot find previous cavity" << endl;
					}
				}
				while (modules[indx1].isCavity())
				{
					indx1++;
					if (indx1 < 0)
					{
						cout << "VECSEL::file_output_structure() cannot find next cavity" << endl;
					}
				}
				double position_x0 = modules[indx0].getPosition1();
				double position_x1 = modules[indx1].getPosition0();
				printStruc << std::fixed << std::setprecision(16) << 2 << " " << position_x0 << " " << position_x1 << " " << "-1" << endl; 
				
			}
			if (modules[i].isTwoArmDevice())
			{
				int indx0 = i;
				int indx1 = i;
				while (modules[indx0].isTwoArmCavity())
				{
					indx0--;
					if (indx0 < 0)
					{
						cout << "VECSEL::file_output_structure() cannot find previous cavity" << endl;
					}
				}
				while (modules[indx1].isTwoArmCavity())
				{
					indx1++;
					if (indx1 < 0)
					{
						cout << "VECSEL::file_output_structure() cannot find next cavity" << endl;
					}
				}
				double position_x0 = modules[indx0].getPosition1();
				double position_x1 = modules[indx1].getPosition0();
				printStruc << std::fixed << std::setprecision(16) << 2 << " " << position_x0 << " " << position_x1 << " " << "-1" << endl; 
				
			}
		//	if (modules[i].isDevice())
		//	{
		//		printStruc << std::fixed << std::setprecision(16) << 2 << " " << modules[i].getPosition0() << " " << modules[i].getPosition1() << " " << "-1" << endl; 
		//	}
			if (modules[i].isLossElement())
			{
				printStruc << std::fixed << std::setprecision(16) << 3 << " " << modules[i].getPosition0() << " " << modules[i].getPosition1() << " " << "-1" << endl; 
			}
		}
		printStruc.close();
	}
}

/* Open all output files everywhere */
void VECSEL::file_output_open_all(int out_count)
{
	if (MPI_MY_RANK==0)
	{
		// Open global stuff
		file_output_open(out_count);
		
		for(unsigned i=0; i<getNumberModules(); i++)
		{
			if (!(modules[i].isDevice()||modules[i].isTwoArmDevice()))
			{
				modules[i].file_output_open(out_count);
			}
		}
	}
	
	#ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		modules[indx].file_output_open(out_count);
	}
	#endif
    
}

/* Close all output files everywhere */
void VECSEL::file_output_close_all()
{
	if (MPI_MY_RANK==0)
	{
		// Close global stuff
		file_output_close();

		for(unsigned i=0; i<getNumberModules(); i++)
		{
			if (!(modules[i].isDevice()||modules[i].isTwoArmDevice()))
			{
				modules[i].file_output_close();
			}
		}
	}
	
	#ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		modules[indx].file_output_close();
	}
	#endif
	
}

/* Write to all output files everywhere */
void VECSEL::file_output_write_all(double t_sim)
{
	if (MPI_MY_RANK==0)
	{
		// Output global stuff
		file_output_write(t_sim);
		
		// For output of Cavities
		for(unsigned i=0; i<getNumberModules(); i++)
		{
			if (!(modules[i].isDevice()||modules[i].isTwoArmDevice()))
			{
				modules[i].file_output_write();
			}
		}
	}
	#ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device	
		modules[indx].file_output_write();
	}
	#endif
	
}

/* Save dynamic variables to files everywhere */
void VECSEL::file_save_variables(int save_count, int offset)
{	
	// Iterate cavities
	int cav_num = 1;
	
	int filt_num = 1;
	
	if (MPI_MY_RANK==0)
	{
		for(unsigned i=0; i<getNumberModules(); i++)
		{
			if (modules[i].isCavity())
			{
				modules[i].getCavity()->file_save_variables(save_count, offset+cav_num); // Simulation variables
				cav_num++;
				
			} else if (modules[i].isTwoArmInterface())
			{
				modules[i].getTwoArmInterface()->file_save_variables(save_count, offset+cav_num); // Simulation variables
				cav_num++;
				
			} else if (modules[i].isTwoArmCavity())
			{
				modules[i].getTwoArmCavity()->file_save_variables(save_count, offset+cav_num); // Simulation variables
				cav_num++;
				
			} else if (modules[i].isFilter())
			{
				modules[i].getFilter()->file_save_variables(save_count,offset+filt_num); // Simulation variables
				filt_num++;
			} 
		}	
	}
	
	#ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		int dev_num = j+1;
		if (modules[indx].isDevice())
		{
			modules[indx].getDevice()->file_save_variables(save_count,offset+dev_num); // Simulation variables
		} else
		{
			modules[indx].getTwoArmDevice()->file_save_variables(save_count,offset+dev_num); // Simulation variables
		}
		//dev_num++;
	}
	#endif
}

/* Load dyanmic variables from files everywhere */
void VECSEL::file_load_variables(int save_count, int offset)
{	
	// Iterate cavities
	int cav_num = 1;
	int filt_num = 1;

	if (MPI_MY_RANK==0)
	{
		for(unsigned i=0; i<getNumberModules(); i++)
		{
			if (modules[i].isCavity())
			{
				modules[i].getCavity()->file_load_variables(save_count, offset+cav_num); // Simulation variables
				cav_num++;
				
			} else if (modules[i].isTwoArmCavity())
			{
				modules[i].getTwoArmCavity()->file_load_variables(save_count, offset+cav_num); // Simulation variables
				cav_num++;
				
			} else if (modules[i].isTwoArmInterface())
			{
				modules[i].getTwoArmInterface()->file_load_variables(save_count, offset+cav_num); // Simulation variables
				cav_num++;
				
			} else if (modules[i].isFilter())
			{
				modules[i].getFilter()->file_load_variables(save_count, offset+filt_num); // Simulation variables
				filt_num++;
			} 
		}
	}
	
	#ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		int dev_num = j+1;
		if (modules[indx].isDevice())
		{
			modules[indx].getDevice()->file_load_variables(save_count, offset+dev_num); // Simulation variables
		} else
		{
			modules[indx].getTwoArmDevice()->file_load_variables(save_count, offset+dev_num); // Simulation variables
		}
		//dev_num++;
	}
	#endif
}

/* Load carriers from files at a spesific time in all devices
 * output_name gives the target output name
 * t_T is the requested time to load carriers at
 *  */
void VECSEL::file_load_device_carriers_time(const std::string &output_name, double t_T)
{
	// Find line corresponding to target time from files with output_name
	// Look through the files: output_name__*_t.dat for the correct file
	
	int fileNum = -1; // Start looking at file = 0
	int lineNum = -1;
	
	std::stringstream fileName;
	fileName << output_name << "0_t.dat";
	
	double tmp;
	double *t_int = new double[2];
	
	while (fileExists(fileName.str()))
	{
		// Load file and check if it containes t_T
		int tmpLines = detectNumberOfLines(fileName.str(), &tmp, 1);
		
		loadBinaryLine(fileName.str(), &t_int[0], 1, 1);
		loadBinaryLine(fileName.str(), &t_int[1], 1, tmpLines);
		
		if ((t_int[0] <= t_T)&&(t_T <= t_int[1]))
		{
			double *allT = new double[tmpLines];
			loadBinary(fileName.str(), allT, tmpLines);

			int targetLine = 0;
			double minT = 1.0;
			for(int j = 0; j < tmpLines; j++)
			{
				if (abs(t_T - allT[j]) < minT)
				{
					targetLine = j;
					minT = abs(t_T - allT[j]);
				}
			}
			
			targetLine++; // Compensate for array index to line #

			cout << "Found in file  = " << fileName.str() << endl;
			cout << "Target time    = " << t_T/ps << " [ps]" << endl;
			cout << "Found index    = " << targetLine << " / " << tmpLines << endl;
			cout << "Deviation time = " << abs(allT[targetLine-1] - t_T)/fs << " [fs]" << endl;


			delete allT;
			
			lineNum = targetLine;
		}
		
		// Check next file
		fileNum++;
		fileName.str("");
		fileName << output_name << fileNum << "_t.dat";
		
		
		if (lineNum != -1)
		{
			break;
		}
	}
	
	if (fileNum == -1)
	{
		cout << "file_load_device_carriers_time(): Could not find requested file " << endl;
		cout << "fileName = " << fileName.str() << endl;
		exit(-1);
	}
	
	if (lineNum==-1)
	{
		cout << "file_load_device_carriers_time(): Could not find requested time " << endl;
		cout << "Time = " << t_T << endl;
		exit(-1);
	}

	file_load_device_carriers_line(output_name, fileNum, lineNum);
}

/* Load device carriers from a given line of output file  */
void VECSEL::file_load_device_carriers_line(int out_count, int line_num)
{	
    #ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		if (modules[indx].isDevice())
		{
			modules[indx].getDevice()->file_load_carriers_line(out_count,line_num); // Simulation variables
		} else
		{
			modules[indx].getTwoArmDevice()->file_load_carriers_line(out_count,line_num); // Simulation variables
		}
	}
	#endif    
}

/* Load device carriers from save point  */
void VECSEL::file_load_device_carriers(int out_count)
{	
    #ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		int dev_num = j+1;
		if (modules[indx].isDevice())
		{
			modules[indx].getDevice()->file_load_carriers(out_count,dev_num); // Simulation variables
		} else
		{
			modules[indx].getTwoArmDevice()->file_load_carriers(out_count,dev_num); // Simulation variables
		}
		dev_num++;
	}
	#endif 
}

/* Load device carriers from a given line of output file  */
void VECSEL::file_load_device_carriers_line(const std::string &output_name, int out_count, int line_num)
{	
    #ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		if (modules[indx].isDevice())
		{
			modules[indx].getDevice()->file_load_carriers_line(output_name,out_count,line_num); // Simulation variables
		} else
		{
			modules[indx].getTwoArmDevice()->file_load_carriers_line(output_name,out_count,line_num); // Simulation variables
		}
	}
	#endif
    
}

/* Load device carriers from an ascii file in a directory 
 * line_num is the folder identifier
 * */
void VECSEL::file_load_device_carriers_ascii(int line_num)
{	
	// Iterate cavities
	#ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		int dev_num = j+1;
		if (modules[indx].isDevice())
		{
			modules[indx].getDevice()->file_load_carriers_ascii(dev_num, line_num); // Simulation variables
		} else
		{
			modules[indx].getTwoArmDevice()->file_load_carriers_ascii(dev_num, line_num); // Simulation variables
		}
		dev_num++;
	}
	#endif


	cout << "file_load_device_carriers_ascii: Load line " << line_num << " Complete" << endl;
}

/* Write output to output files 
 * t_sim is the current simulation time
 * */
void VECSEL::file_output_write(double t_sim)
{
	if (MPI_MY_RANK==0)
	{
		output_simulation_time.write(reinterpret_cast<const char*>(&t_sim),sizeof(double));
		
		
		#ifdef USE_CAVITY_SNAPSHOT
			if (VECSEL_cav_snapshot_output_count >= VECSEL_cav_snapshot_output_wait)
			{
				file_snapshot_write(t_sim);
				VECSEL_cav_snapshot_output_count = 0;
			} else {
				VECSEL_cav_snapshot_output_count++;
			}
		#endif
		
		#ifdef DUAL_PROP
			if(true) //(modules[0].isBoundary() && VECSEL_pulse_start_l == 1)
			{
			/*==============================
				 Output of device on LEFT
			==============================*/
				int indx = quick_index_twoArmCavity[0]; // Final cavity
				double rR = modules[0].getBoundary()->getRefCoeff(); // Right boundary
				double ni = modules[indx].getRefInd();

				std::complex<double> E_bp[VECSEL_transverse_points_number];
				std::complex<double> E_bm[VECSEL_transverse_points_number];
				std::complex<double> E_prop;
				modules[indx].getTwoArmCavity()->getEbp_front_wall(E_bp);
				modules[indx].getTwoArmCavity()->getEbm_front_wall(E_bm);
				
				for(int i = 0; i < VECSEL_transverse_points_number; i++)
				{
					// output
					E_prop=E_bp[i];
					double tmp = real(sqrt(ni*(1.0-rR))*E_prop);
					output_E_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
				
					tmp = imag(sqrt(ni*(1.0-rR))*E_prop);
					output_E_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
					
					E_prop=E_bm[i];
					tmp = real(sqrt(ni*(1.0-rR))*E_prop);
					output_back_E_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
				
					tmp = imag(sqrt(ni*(1.0-rR))*E_prop);
					output_back_E_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
				}
			}

		#else
			if(true) //(modules[0].isBoundary() && VECSEL_pulse_start_l == 1)
			{
			/*==============================
				 Output of device on LEFT
			==============================*/
				int indx = quick_index_cavity[0]; // Final cavity
				double rR = modules[0].getBoundary()->getRefCoeff(); // Right boundary
				double ni = modules[indx].getRefInd();

				std::complex<double> E_prop[VECSEL_transverse_points_number];
				//std::complex<double> E_plus[VECSEL_transverse_points_number];
				//std::complex<double> E_minus[VECSEL_transverse_points_number];
				modules[indx].getCavity()->getEminus_left_wall(E_prop);
				//modules[indx].getCavity()->getEminus_right_wall(E_minus);
				for(int i = 0; i < VECSEL_transverse_points_number; i++)
				{
					// output
					double tmp = real(sqrt(ni*(1.0-rR))*E_prop[i]);
					output_E_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
				
					tmp = imag(sqrt(ni*(1.0-rR))*E_prop[i]);
					output_E_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
				}
			}

			if(true) //(modules.back().isBoundary() && VECSEL_pulse_start_r ==1) 
			{
			/*==============================
				 Output of device on RIGHT
			==============================*/
				int indx = quick_index_cavity.back(); // Final cavity
				double rR = modules.back().getBoundary()->getRefCoeff(); // Right boundary
				double ni = modules[indx].getRefInd();

				std::complex<double> E_prop[VECSEL_transverse_points_number];
				//std::complex<double> E_plus[VECSEL_transverse_points_number];
				//std::complex<double> E_minus[VECSEL_transverse_points_number];
				modules[indx].getCavity()->getEpluss_right_wall(E_prop);
				//modules[indx].getCavity()->getEminus_right_wall(E_minus);
				for(int i = 0; i < VECSEL_transverse_points_number; i++)
				{
					// output
					double tmp = real(sqrt(ni*(1.0-rR))*E_prop[i]);
					output_back_E_real[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
				
					tmp = imag(sqrt(ni*(1.0-rR))*E_prop[i]);
					output_back_E_imag[i].write(reinterpret_cast<const char*>(&tmp),sizeof(double));
				}
			}
		#endif
	}
}

/* Initialize output files, open with correct numbering 
 * out_count is the output identification number
 * */
void VECSEL::file_output_open(int out_count)
{
	if (MPI_MY_RANK==0)
	{
		// Setup output to files
		std::stringstream baseName;
		baseName << getToFileOutputKey() << out_count;
		
		std::stringstream fileName;
		fileName << baseName.str() << "_t.dat";
		//output_simulation_time.open((fileName.str()).c_str(), std::ofstream::out|std::ofstream::binary|std::ofstream::app);
		openAppendBinary(&output_simulation_time, fileName.str());

		for(int i = 0; i < VECSEL_transverse_points_number; i++)
		{
			if(true)//( VECSEL_pulse_start_l == 1 )
			{
				fileName.str("");
				fileName << baseName.str() << "_E_re_OUTPUT_T" << i << ".dat";
				openAppendBinary(&output_E_real[i], fileName.str());
		
				fileName.str("");
				fileName << baseName.str() << "_E_im_OUTPUT_T" << i << ".dat";
				openAppendBinary(&output_E_imag[i], fileName.str());
			}
	
			if(true)//( VECSEL_pulse_start_r == 1)
			{
				fileName.str("");
				fileName << baseName.str() << "_E_re_OUTPUTBACK_T" << i << ".dat";
				openAppendBinary(&output_back_E_real[i], fileName.str());
		
				fileName.str("");
				fileName << baseName.str() << "_E_im_OUTPUTBACK_T" << i << ".dat";
				openAppendBinary(&output_back_E_imag[i], fileName.str());
			}
		}
		
	
		#ifdef USE_CAVITY_SNAPSHOT
		VECSEL_cav_snapshot_output_count = VECSEL_cav_snapshot_output_wait;
		fileName.str("");
		fileName << baseName.str() << "_cav_snapshot_t.dat";
		openAppendBinary(&output_Cav_Snapshot_t, fileName.str());
		
		for(int i = 0; i < VECSEL_transverse_points_number; i++)
		{
			fileName.str("");
			fileName << baseName.str() << "_cav_snapshot_E_re_T"<< i <<".dat";
			openAppendBinary(&output_Cav_Snapshot_E_real[i], fileName.str());
			
			fileName.str("");
			fileName << baseName.str() << "_cav_snapshot_E_im_T"<< i <<".dat";
			openAppendBinary(&output_Cav_Snapshot_E_imag[i], fileName.str());
		}
		#endif
	}
}

/* Close all output files */
void VECSEL::file_output_close()
{
	if (MPI_MY_RANK==0)
	{
		output_simulation_time.close();
		for(int i = 0; i < VECSEL_transverse_points_number; i++)
		{
		
			if(true) //( VECSEL_pulse_start_l == 1)
			{
				output_E_real[i].close();
				output_E_imag[i].close();
			}

			if(true) //( VECSEL_pulse_start_r == 1)
			{
				output_back_E_real[i].close();
				output_back_E_imag[i].close();
			}
		}
		
		#ifdef USE_CAVITY_SNAPSHOT
		output_Cav_Snapshot_t.close();
		for(int i = 0; i < VECSEL_transverse_points_number; i++)
		{
			output_Cav_Snapshot_E_real[i].close();
			output_Cav_Snapshot_E_imag[i].close();
		}
		#endif
	}
}

/* Print out the field in the z-direction in the cavity
 * Requires alot of calculations, so takes a long time
 * VECSEL_cav_snapshot_num_points: Sets number of z-points to print at
 * */
void VECSEL::file_snapshot_write(double t)
{
	cout << "VECSEL::file_snapshot_write() Snapshot not updated for transverse dimensions.." << endl;
	exit(-1);
/*
	for(int i = 0; i < VECSEL_cav_snapshot_num_points; i++)
	{
		int indx = VECSEL_cav_snapshot_index[i]; // Find correct cavity
		
		VECSEL_cav_snapshot_E[i] = modules[indx].getCavity()->evaluateEprop(VECSEL_cav_snapshot_x[i]);
		
		VECSEL_cav_snapshot_E_re[i] = real(VECSEL_cav_snapshot_E[i]);
		VECSEL_cav_snapshot_E_im[i] = imag(VECSEL_cav_snapshot_E[i]);
		
	}
	
	output_Cav_Snapshot_E_real.write(reinterpret_cast<const char*>(VECSEL_cav_snapshot_E_re),VECSEL_cav_snapshot_num_points*sizeof(double));
	output_Cav_Snapshot_E_imag.write(reinterpret_cast<const char*>(VECSEL_cav_snapshot_E_im),VECSEL_cav_snapshot_num_points*sizeof(double));
	output_Cav_Snapshot_t.write(reinterpret_cast<const char*>(&t),sizeof(double));
*/
}

/* Print out ONLY ONCE the field in the z-direction in the cavity 
 * Requires alot of calculations, so takes a long time
 * VECSEL_cav_snapshot_num_points: Sets number of z-points to print at 
 * */
void VECSEL::file_snapshot_write_single(int ID, int NUM_POINTS)
{
	cout << "VECSEL::file_snapshot_write_single() Snapshot not updated for transverse dimensions.." << endl;
	exit(-1);
/*
	std::stringstream fileName;
	std::stringstream baseName;
	baseName << getToFileOutputKey() << ID;
	
	fileName << baseName.str() << "_cav_snapshot_E_re.dat";
	openAppendBinary(&output_Cav_Snapshot_E_real, fileName.str());
	
	fileName.str("");
	fileName << baseName.str() << "_cav_snapshot_E_im.dat";
	openAppendBinary(&output_Cav_Snapshot_E_imag, fileName.str());
	
	//=================
	// Compute x, data
	//=================
	baseName.str("");
	baseName << getToFileOutputKey() << ID;

	fileName.str("");
	fileName << baseName.str() << "_cav_snapshot_x.dat";
	openAppendBinary(&output_Cav_Snapshot_x, fileName.str());

	double DX = (modules.back().getPosition1() - modules[0].getPosition0())/((double)NUM_POINTS);
	#ifndef USE_CAVITY_SNAPSHOT
	//VECSEL_cav_snapshot_x = new double[NUM_POINTS];
	#endif
	for(int i =0; i < NUM_POINTS; i++)
	{
		VECSEL_cav_snapshot_x[i] = modules[0].getPosition0() + ((double)i+0.5)*DX;
	}
	output_Cav_Snapshot_x.write(reinterpret_cast<const char*>(&VECSEL_cav_snapshot_x[0]),NUM_POINTS*sizeof(double));
	output_Cav_Snapshot_x.close();
	
	
	
	// Find correct indices
	#ifndef USE_CAVITY_SNAPSHOT
	//VECSEL_cav_snapshot_index = new int[NUM_POINTS];
	#endif
	for(int i =0; i < NUM_POINTS; i++)
	{
		// Find correct cavity
		for(int j = 0; j < modules.size(); j++)
		{
			if (((VECSEL_cav_snapshot_x[i] >= modules[j].getPosition0())&&(VECSEL_cav_snapshot_x[i] <= modules[j].getPosition1()))&&(modules[j].isCavity()))
			{
				VECSEL_cav_snapshot_index[i] = j;
				break;
			} 
		}
	}

	// Initialize E
	#ifndef USE_CAVITY_SNAPSHOT
	VECSEL_cav_snapshot_E 	= new std::complex<double>[NUM_POINTS];
	VECSEL_cav_snapshot_E_re = new double[NUM_POINTS];
	VECSEL_cav_snapshot_E_im = new double[NUM_POINTS];
	#endif
	
	
	for(int i = 0; i < NUM_POINTS; i++)
	{
		int indx = VECSEL_cav_snapshot_index[i]; // Find correct cavity
		
		VECSEL_cav_snapshot_E[i] = modules[indx].getCavity()->evaluateEprop(VECSEL_cav_snapshot_x[i]);
		VECSEL_cav_snapshot_E_re[i] = real(VECSEL_cav_snapshot_E[i]);
		VECSEL_cav_snapshot_E_im[i] = imag(VECSEL_cav_snapshot_E[i]);
	}
	
	output_Cav_Snapshot_E_real.write(reinterpret_cast<const char*>(VECSEL_cav_snapshot_E_re),NUM_POINTS*sizeof(double));
	output_Cav_Snapshot_E_imag.write(reinterpret_cast<const char*>(VECSEL_cav_snapshot_E_im),NUM_POINTS*sizeof(double));
	//output_Cav_Snapshot_t.write(reinterpret_cast<const char*>(&t),sizeof(double));
	
	
	//output_Cav_Snapshot_t.close();
	output_Cav_Snapshot_E_real.close();
	output_Cav_Snapshot_E_imag.close();
*/
}



//===================
// Maxwells function

/**
	Update an edge 'i'  on to the RIGHT of a cavity, inside the domain using the iteration scheme
	The contribution from any devices on the edge is precomputed in 'device_MacPol'
	
	device_MacPol -> Any contribution from the device itself = MacPol*(eff_qw/focusE)
	
	Should not be used on boundary edges such as Vcav, boundary or others
*/

void VECSEL::iterateModules_updateSingleSurface_transfer_matrix(int i)
{
	std::complex<double> a11, a12, a21, a22, *E_pl, *E_mi, *E_pl_k, *E_mi_k, *MacPol;
	
	
	int indx   = quick_index_cavity[i];			// Current cavity
	int indx_k = quick_index_cavity[i+1];     	// Next cavity


	modules[indx].getCavity()->get_transfer_matrix(&a11,&a12,&a21,&a22,&MacPol);
	
	E_pl   = modules[indx  ].getCavity()->interpolateEpluss_x1();
	E_mi_k = modules[indx_k].getCavity()->interpolateEminus_x0();

	E_pl_k = modules[indx_k].getCavity()->setEpluss();
	E_mi   = modules[indx].getCavity()->setEminus();

	cblas_zcopy(VECSEL_transverse_points_number, MacPol, 1  , E_pl_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11  , E_pl  , 1  , E_pl_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a12  , E_mi_k, 1  , E_pl_k, 1);
	
	cblas_zcopy(VECSEL_transverse_points_number, MacPol, 1  , E_mi, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21  , E_pl  , 1, E_mi, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a22  , E_mi_k, 1, E_mi, 1);
}

void VECSEL::iterateModules_updateSingleSurface_transfer_matrix(Cavity *cav0, Cavity *cav1)
{
	std::complex<double> a11, a12, a21, a22, *E_pl, *E_mi, *E_pl_k, *E_mi_k, *MacPol;
	
	cav0->get_transfer_matrix(&a11,&a12,&a21,&a22,&MacPol);
	
	E_pl   = cav0->interpolateEpluss_x1();
	E_mi_k = cav1->interpolateEminus_x0();

	E_pl_k = cav1->setEpluss();
	E_mi   = cav0->setEminus();
	
	cblas_zcopy(VECSEL_transverse_points_number, MacPol, 1  , E_pl_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11  , E_pl  , 1  , E_pl_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a12  , E_mi_k, 1  , E_pl_k, 1);
	
	cblas_zcopy(VECSEL_transverse_points_number, MacPol, 1  , E_mi, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21  , E_pl  , 1, E_mi, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a22  , E_mi_k, 1, E_mi, 1);
}

void VECSEL::iterateModules_updateSingleSurface_transfer_matrix_noQW(Cavity *cav0, Cavity *cav1)
{
	std::complex<double> a11, a12, a21, a22, *E_pl, *E_mi, *E_mi_k, *E_pl_k, *MacPol;

	cav0->get_transfer_matrix(&a11,&a12,&a21,&a22,&MacPol);
	
	
	E_pl   = cav0->interpolateEpluss_x1();
	E_mi_k = cav1->interpolateEminus_x0();

	E_pl_k = cav1->setEpluss();
	E_mi   = cav0->setEminus();
	
	memset(E_pl_k, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a11  , E_pl  , 1  , E_pl_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a12  , E_mi_k, 1  , E_pl_k, 1);
	
	memset(E_mi, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a21  , E_pl  , 1, E_mi, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a22  , E_mi_k, 1, E_mi, 1);
}

void VECSEL::iterateModules_updateSingleSurface_transfer_matrix_noQW_debug(Cavity *cav0, Cavity *cav1)
{
	std::complex<double> a11, a12, a21, a22, *E_pl, *E_mi, *E_pl_k, *MacPol;
	
	std::complex<double> *E_mi_k;

	cav0->get_transfer_matrix(&a11,&a12,&a21,&a22,&MacPol);
	
	
	E_pl   = cav0->interpolateEpluss_x1();
	E_mi_k = cav1->interpolateEminus_x0();

	E_pl_k = cav1->setEpluss();
	E_mi   = cav0->setEminus();

	
	memset(E_pl_k, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a11  , E_pl  , 1  , E_pl_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a12  , E_mi_k, 1  , E_pl_k, 1);
	
	memset(E_mi, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a21  , E_pl  , 1, E_mi, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a22  , E_mi_k, 1, E_mi, 1);
}

void VECSEL::iterateModules_updateSingleSurface_TwoArm_transfer_matrix(int i)
{
	std::complex<double> a11, a12, a21, a22, *E_fp, *E_fp_k, *E_fm, *E_fm_k, *E_bp, *E_bp_k, *E_bm, *E_bm_k, *MacPol_fp, *MacPol_fm, *MacPol_bp, *MacPol_bm;
	
	int indx   = quick_index_twoArmCavity[i];	// Current cavity
	int indx_k = quick_index_twoArmCavity[i+1];     	// Next cavity
	
	std::complex<double> dummy=1.0;

	modules[indx ].getTwoArmCavity()->get_transfer_matrix(&a11,&a12,&a21,&a22,&MacPol_fp, &MacPol_fm, &MacPol_bp, &MacPol_bm);
	
	E_fp   = modules[indx ].getTwoArmCavity()->interpolateEfp_x1();
	E_fm   = modules[indx ].getTwoArmCavity()->interpolateEfm_x1();
	E_bp_k = modules[indx_k].getTwoArmCavity()->interpolateEbp_x0();
	E_bm_k = modules[indx_k].getTwoArmCavity()->interpolateEbm_x0();

	E_fp_k = modules[indx_k].getTwoArmCavity()->setEfp();
	E_fm_k = modules[indx_k].getTwoArmCavity()->setEfm();
	E_bp   = modules[indx ].getTwoArmCavity()->setEbp();
	E_bm   = modules[indx ].getTwoArmCavity()->setEbm();


	cblas_zcopy(VECSEL_transverse_points_number, MacPol_fp, 1  , E_fp_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &dummy,  MacPol_bp, 1  , E_fp_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a12  , E_bp_k  , 1  , E_fp_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11  , E_fp, 1  , E_fp_k, 1);
	
	cblas_zcopy(VECSEL_transverse_points_number, MacPol_fm, 1  , E_fm_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &dummy, MacPol_bm, 1  , E_fm_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a12  , E_bm_k, 1, E_fm_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11  , E_fm  , 1, E_fm_k, 1);

	cblas_zcopy(VECSEL_transverse_points_number, MacPol_fp, 1  , E_bp, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &dummy, MacPol_bp, 1  , E_bp, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21  , E_fp  , 1  , E_bp, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a22  , E_bp_k, 1  , E_bp, 1);
	
	cblas_zcopy(VECSEL_transverse_points_number, MacPol_fm, 1  , E_bm, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &dummy, MacPol_bm, 1  , E_bm, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21  , E_fm  , 1, E_bm, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a22  , E_bm_k, 1, E_bm, 1);
}

void VECSEL::iterateModules_updateSingleSurface_TwoArm_transfer_matrix(TwoArmCavity *cav0, TwoArmCavity *cav1)
{
	std::complex<double> a11, a12, a21, a22, *E_fp, *E_fp_k, *E_fm, *E_fm_k, *E_bp, *E_bp_k, *E_bm, *E_bm_k, *MacPol_fp, *MacPol_fm, *MacPol_bp, *MacPol_bm;
	std::complex<double> dummy=1.0;
	
	cav0->get_transfer_matrix(&a11,&a12,&a21,&a22,&MacPol_fp, &MacPol_fm, &MacPol_bp, &MacPol_bm);

	E_fp   = cav0->interpolateEfp_x1();
	E_bp_k = cav1->interpolateEbp_x0();
	E_bm_k = cav1->interpolateEbm_x0();
	E_fm   = cav0->interpolateEfm_x1();

	E_fp_k = cav1->setEfp();
	E_fm_k = cav1->setEfm();
	E_bp   = cav0->setEbp();
	E_bm   = cav0->setEbm();
	
	cblas_zcopy(VECSEL_transverse_points_number, MacPol_fp, 1  , E_fp_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &dummy, MacPol_bp, 1  , E_fp_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a12  , E_bp_k, 1  , E_fp_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11  , E_fp  , 1  , E_fp_k, 1);

	cblas_zcopy(VECSEL_transverse_points_number, MacPol_fm, 1  , E_fm_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &dummy, MacPol_bm, 1  , E_fm_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a12  , E_bm_k, 1, E_fm_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11  , E_fm  , 1, E_fm_k, 1);
	
	cblas_zcopy(VECSEL_transverse_points_number, MacPol_bp, 1  , E_bp, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &dummy, MacPol_fp, 1  , E_bp, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a22  , E_bp_k  , 1  , E_bp, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21  , E_fp   , 1  , E_bp, 1);
	
	cblas_zcopy(VECSEL_transverse_points_number, MacPol_bm, 1  , E_bm, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &dummy, MacPol_fm, 1  , E_bm, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a22  , E_bm_k, 1, E_bm, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21  , E_fm  , 1, E_bm, 1);
}

void VECSEL::iterateModules_updateSingleSurface_TwoArm_transfer_matrix_interface_delay(Cavity *cav0, TwoArmInterface *cavInt, Cavity *cav1, TwoArmCavity *cavFront)
{
	std::complex<double> a11_left, a12_left, a21_left, a22_left, a11_right, a12_right, a21_right, a22_right, *E_p_in, *E_m_in, *E_bm, *E_bp, *E_p_out, *E_m_out, *E_fp, *E_fm, *phase_in, *phase_out;

	std::complex<double>  E_fp_int[VECSEL_transverse_points_number], E_fm_int[VECSEL_transverse_points_number], E_bp_int[VECSEL_transverse_points_number], E_bm_int[VECSEL_transverse_points_number];
	
	std::complex<double> dummy=1.0;
	std::complex<double> E_tmp;

	cavInt->get_transfer_matrix(&a11_left,&a12_left,&a21_left,&a22_left, &a11_right, &a12_right, &a21_right, &a22_right);
	
	E_p_in  = cav0->interpolateEpluss_x1();
	E_m_in  = cav1->interpolateEminus_x0();
	E_bm    = cavFront->interpolateEbm_x0();
	E_bp    = cavFront->interpolateEbp_x0();

	E_m_out = cav1->setEpluss();
	E_p_out = cav0->setEminus();
	E_fp    = cavFront->setEfp();
	E_fm    = cavFront->setEfm();

	phase_in=cavInt->getPhaseIn();
	phase_out=cavInt->getPhaseOut();
	
	cavInt->getBoundary_Efp(E_fp_int);	
	cavInt->getBoundary_Efm(E_fm_int);	
	cavInt->getBoundary_Ebp(E_bp_int);	
	cavInt->getBoundary_Ebm(E_bm_int);	

	for( int i =0; i < VECSEL_transverse_points_number; i++)
	{
		//Forward field into interface
		cavInt->setEfp(i, &E_p_in[i]);
		cavInt->setEfm(i, &E_m_in[i]);
		
		//Backward pluss field	
		E_tmp=a21_right*E_fm_int[i]+a22_right*phase_in[i]*E_bm[i];
		//E_tmp=a21_right*E_fm_int[i]+a22_right*E_bm[i];
		cavInt->setEbp(i, &E_tmp);
	
		//Backward minus field
		E_tmp=a21_left*E_fp_int[i]+a22_left*phase_out[i]*E_bp[i];
		//E_tmp=a21_left*E_fp_int[i]+a22_left*E_bp[i];
		cavInt->setEbm(i, &E_tmp);	
		
		E_fp_int[i]*=phase_in[i];
		E_fm_int[i]*=phase_out[i];
	}
	
	
	memset(E_fp, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a12_left  , E_bp      , 1 , E_fp, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11_left  , E_fp_int  , 1 , E_fp, 1);
	
	memset(E_fm, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a12_right  , E_bm     , 1 , E_fm, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11_right  , E_fm_int , 1 , E_fm, 1);
	
	memset(E_m_out, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &dummy  , E_bm_int , 1 , E_m_out, 1);

	memset(E_p_out, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));	
	cblas_zaxpy(VECSEL_transverse_points_number, &dummy  , E_bp_int , 1 , E_p_out, 1);

}

void VECSEL::iterateModules_updateSingleSurface_TwoArm_transfer_matrix_interface(Cavity *cav0, TwoArmInterface *cavInt, Cavity *cav1, TwoArmCavity *cavFront)
{
	std::complex<double> a11_left, a12_left, a21_left, a22_left, a11_right, a12_right, a21_right, a22_right, *E_m_in, *E_p_in, *E_bm, *E_bp, *E_p_out, *E_m_out, *E_fp, *E_fm;
	cavInt->get_transfer_matrix(&a11_left,&a12_left,&a21_left,&a22_left, &a11_right, &a12_right, &a21_right, &a22_right);

	E_p_in  = cav0->interpolateEpluss_x1();
	E_m_in  = cav1->interpolateEminus_x0();
	E_bm    = cavFront->interpolateEbm_x0();
	E_bp    = cavFront->interpolateEbp_x0();
	
	E_m_out = cav1->setEpluss();
	E_p_out = cav0->setEminus();
	E_fp    = cavFront->setEfp();
	E_fm    = cavFront->setEfm();

	
	memset(E_fp, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a12_left  , E_bp  , 1  , E_fp, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11_left  , E_p_in, 1  , E_fp, 1);
	
	memset(E_fm, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a12_right  , E_bm  , 1, E_fm, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11_right  , E_m_in, 1, E_fm, 1);
	
	memset(E_m_out, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	memset(E_p_out, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	
	//VCAV implementation
	cblas_zaxpy(VECSEL_transverse_points_number, &a22_right  , E_bp  , 1  , E_m_out, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21_left   , E_p_in, 1  , E_m_out, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a22_left , E_bm  , 1, E_p_out, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21_right, E_m_in , 1, E_p_out, 1);
	
	//MockLinear implementation
	//cblas_zaxpy(VECSEL_transverse_points_number, &a22_right  , E_bm  , 1  , E_m_out, 1);
	//cblas_zaxpy(VECSEL_transverse_points_number, &a21_right  , E_m_in, 1  , E_m_out, 1);	
	//cblas_zaxpy(VECSEL_transverse_points_number, &a22_left  , E_bp  , 1  , E_p_out, 1);
	//cblas_zaxpy(VECSEL_transverse_points_number, &a21_left  , E_p_in, 1  , E_p_out, 1);

}

void VECSEL::iterateModules_updateSingleSurface_TwoArm_transfer_matrix_back(TwoArmCavity *cav0, std::complex<double> Reflection)
{	
	std::complex<double> *E_bp, *E_bm;
	std::complex<double>  E_fp[VECSEL_transverse_points_number], E_fm[VECSEL_transverse_points_number];

	cav0->getEfp_back_wall(E_fp);
	cav0->getEfm_back_wall(E_fm);
	
	E_bp   = cav0->setEbp();
	E_bm   = cav0->setEbm();

	memset(E_bp, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &Reflection, E_fp, 1  , E_bp, 1);
	memset(E_bm, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &Reflection, E_fm, 1  , E_bm, 1);
}

void VECSEL::iterateModules_updateSingleSurface_TwoArm_transfer_matrix_noQW(TwoArmCavity *cav0, TwoArmCavity *cav1)
{
	std::complex<double> a11, a12, a21, a22, *E_fp, *E_fp_k, *E_fm, *E_fm_k, *E_bp, *E_bp_k, *E_bm, *E_bm_k, *MacPol_fp, *MacPol_fm, *MacPol_bp, *MacPol_bm;
	
	cav0->get_transfer_matrix(&a11,&a12,&a21,&a22,&MacPol_fp, &MacPol_fm, &MacPol_bp, &MacPol_bm);

	E_fp   = cav0->interpolateEfp_x1();
	E_bp_k = cav1->interpolateEbp_x0();
	E_fm   = cav0->interpolateEfm_x1();
	E_bm_k = cav1->interpolateEbm_x0();

	E_fp_k = cav1->setEfp();
	E_bp   = cav0->setEbp();
	E_fm_k = cav1->setEfm();
	E_bm   = cav0->setEbm();
	
	memset(E_fp_k, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a12  , E_bp_k , 1  , E_fp_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11  , E_fp   , 1  , E_fp_k, 1);
	
	memset(E_fm_k, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a12  , E_bm_k , 1, E_fm_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11  , E_fm   , 1, E_fm_k, 1);

	memset(E_bp, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a22  , E_bp_k , 1  , E_bp, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21  , E_fp   , 1  , E_bp, 1);
	
	memset(E_bm, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a22  , E_bm_k  , 1, E_bm, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21  , E_fm , 1, E_bm, 1);
}

void VECSEL::iterateModules_updateSingleSurface_TwoArm_transfer_matrix_birefringentCrystal_front(TwoArmCavity *cav0, TwoArmCavity *brc)
{
	std::complex<double> a11, a12, a21, a22, a11_ex, a12_ex, a21_ex, a22_ex, *E_fp, *E_fp_k, *E_fm, *E_fm_k, *E_bp, *E_bp_k, *E_bm, *E_bm_k, *MacPol_fp, *MacPol_fm, *MacPol_bp, *MacPol_bm;
	
	cav0->get_transfer_matrix(&a11,&a12,&a21,&a22,&MacPol_fp, &MacPol_fm, &MacPol_bp, &MacPol_bm);
	brc->get_transfer_matrix_extraAxis(&a11_ex,&a12_ex,&a21_ex,&a22_ex);

	E_fp     = cav0->interpolateEfp_x1();
	E_bp_k   = brc->interpolateEbp_x0();
	E_fm     = cav0->interpolateEfm_x1();
	E_bm_k   = brc->interpolateEbm_x0();
	
	E_fp_k   = brc->setEfp();
	E_bp     = cav0->setEbp();
	E_fm_k   = brc->setEfm();
	E_bm     = cav0->setEbm();


	memset(E_fp_k, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a12     , E_bp_k , 1  , E_fp_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11     , E_fp   , 1  , E_fp_k, 1);
	
	memset(E_fm_k, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a12_ex  , E_bm_k , 1, E_fm_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11     , E_fm   , 1, E_fm_k, 1);

	memset(E_bp, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a22     , E_bp_k , 1  , E_bp, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21     , E_fp   , 1  , E_bp, 1);
	
	memset(E_bm, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a22_ex  , E_bm_k  , 1, E_bm, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21     , E_fm , 1, E_bm, 1);
	
}

void VECSEL::iterateModules_updateSingleSurface_TwoArm_transfer_matrix_birefringentCrystal_back(TwoArmCavity *brc, TwoArmCavity *cav1)
{
	std::complex<double> a11, a12, a21, a22, a11_ex, a12_ex, a21_ex, a22_ex, *E_fp, *E_fp_k, *E_fm, *E_fm_k, *E_bp, *E_bp_k, *E_bm, *E_bm_k, *MacPol_fp, *MacPol_fm, *MacPol_bp, *MacPol_bm;
	
	brc->get_transfer_matrix(&a11,&a12,&a21,&a22,&MacPol_fp, &MacPol_fm, &MacPol_bp, &MacPol_bm);
	brc->get_transfer_matrix_extraAxis(&a11_ex,&a12_ex,&a21_ex,&a22_ex);

	E_fp     = brc->interpolateEfp_x1();
	E_bp_k   = cav1->interpolateEbp_x0();
	E_fm     = brc->interpolateEfm_x1();
	E_bm_k   = cav1->interpolateEbm_x0();
	
	E_fp_k   = cav1->setEfp();
	E_bp     = brc->setEbp();
	E_fm_k   = cav1->setEfm();
	E_bm     = brc->setEbm();


	memset(E_fp_k, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a12     , E_bp_k , 1  , E_fp_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11     , E_fp   , 1  , E_fp_k, 1);
	
	memset(E_fm_k, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a12     , E_bm_k , 1, E_fm_k, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a11_ex  , E_fm   , 1, E_fm_k, 1);

	memset(E_bp, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a22     , E_bp_k , 1  , E_bp, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21     , E_fp   , 1  , E_bp, 1);
	
	memset(E_bm, 0.0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	cblas_zaxpy(VECSEL_transverse_points_number, &a22     , E_bm_k  , 1, E_bm, 1);
	cblas_zaxpy(VECSEL_transverse_points_number, &a21_ex  , E_fm , 1, E_bm, 1);
	
}


void VECSEL::iterateModules_updateAllSurface(VECSEL *model, double t_sim)
{
	//================================================================
	// 1. Apply boundary conditions: Reflecting/Absorbing periodic boundaries: 1st boundary
	//================================================================
	if (model->modules[0].isBoundary())
	{
	
		#ifdef DUAL_PROP
			int indx   = model->quick_index_twoArmCavity[0];			// Current cavity
			int indx_k = model->quick_index_twoArmCavity.back();     		// Last cavity
			
			double r1 = model->modules[0].getBoundary()->getRefCoeff();
			double r2 = model->modules.back().getBoundary()->getRefCoeff();
			double x0 = model->modules[indx].getPosition0();
			double ni = model->modules[indx].getRefInd();
			double n0 = model->modules[0].getBoundary()->getNextCavityIndex();

			// Reflected from cavity
			model->modules[indx].getTwoArmCavity()->getEbp_front_wall(cavity_trans_E_mi);
			model->modules[0].getBoundary()->getEpluss(cavity_trans_E_pl);

			// Transmission through boundary for output
			std::complex<double> *transmitted = model->modules[0].getBoundary()->setEminus();
			std::complex<double> scale = sqrt((ni/n0)*(1.0-r1));
			memcpy(transmitted, cavity_trans_E_mi, VECSEL_transverse_points_number*sizeof(std::complex<double>));
			cblas_zscal(VECSEL_transverse_points_number, &scale, transmitted, 1);
			model->modules[indx].getTwoArmCavity()->getEbm_front_wall(cavity_trans_E_mi);
			cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_mi , 1, transmitted, 1);
			
			// Feedback into cavity
			std::complex<double> *feedback = model->modules[indx].getTwoArmCavity()->setEfp();
			
			#ifdef CYCLIC_BOUNDARIES
				memcpy(feedback, cavity_trans_E_pl, VECSEL_transverse_points_number*sizeof(std::complex<double>));
				model->modules[indx_k].getTwoArmCavity()->getEfp_right_wall(cavity_trans_E_pl);
				scale = -sqrt(r2);
				cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_pl , 1, feedback, 1);
			#else
				scale = -sqrt(r1);
				memcpy(feedback, cavity_trans_E_pl, VECSEL_transverse_points_number*sizeof(std::complex<double>));
				model->modules[indx].getTwoArmCavity()->getEbp_front_wall(cavity_trans_E_mi);
				cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_mi , 1, feedback, 1);
			#endif
				
			feedback = model->modules[indx].getTwoArmCavity()->setEfm();
			#ifdef CYCLIC_BOUNDARIES
				memcpy(feedback, cavity_trans_E_pl, VECSEL_transverse_points_number*sizeof(std::complex<double>));
				model->modules[indx_k].getTwoArmCavity()->getEfm_right_wall(cavity_trans_E_pl);
				scale = -sqrt(r2);
				cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_pl , 1, feedback, 1);
			#else
				scale = -sqrt(r1);
				model->modules[indx].getTwoArmCavity()->getEbm_front_wall(cavity_trans_E_mi);
				memcpy(feedback, cavity_trans_E_pl, VECSEL_transverse_points_number*sizeof(std::complex<double>));
				cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_mi , 1, feedback, 1);
			#endif
		
		#else
			int indx   = model->quick_index_cavity[0];			// Current cavity
			int indx_k = model->quick_index_cavity.back();     		// Last cavity
		
			double r1 = model->modules[0].getBoundary()->getRefCoeff();
			double r2 = model->modules.back().getBoundary()->getRefCoeff();
			double x0 = model->modules[indx].getPosition0();
			double ni = model->modules[indx].getRefInd();
			double n0 = model->modules[0].getBoundary()->getNextCavityIndex();
			
			// Reflected from cavity
			model->modules[indx].getCavity()->getEminus_left_wall(cavity_trans_E_mi);
			model->modules[0].getBoundary()->getEpluss(cavity_trans_E_pl);

			// Transmission through boundary for output
			std::complex<double> *transmitted = model->modules[0].getBoundary()->setEminus();
			std::complex<double> scale = sqrt((ni/n0)*(1.0-r1));
			memcpy(transmitted, cavity_trans_E_mi, VECSEL_transverse_points_number*sizeof(std::complex<double>));
			cblas_zscal(VECSEL_transverse_points_number, &scale, transmitted, 1);
			
			// Feedback into cavity
			std::complex<double> *feedback = model->modules[indx].getCavity()->setEpluss();
			#ifdef CYCLIC_BOUNDARIES
				memcpy(feedback, cavity_trans_E_pl, VECSEL_transverse_points_number*sizeof(std::complex<double>));
				model->modules[indx_k].getCavity()->getEpluss_right_wall(cavity_trans_E_pl);
				scale = -sqrt(r2);
				cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_pl , 1, feedback, 1);
			#else
				scale = -sqrt(r1);
				memcpy(feedback, cavity_trans_E_pl, VECSEL_transverse_points_number*sizeof(std::complex<double>));
				cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_mi , 1, feedback, 1);
			#endif
		#endif
		
	}
	
	//================================================================
	// 2. Apply boundary conditions: Reflecting/Absorbing periodic boundaries: Last boundary
	//================================================================
	if (model->modules.back().isBoundary())
	{
	
		#ifdef DUAL_PROP
			int indx   = model->quick_index_twoArmCavity[0];			// Current cavity
			int indx_k = model->quick_index_twoArmCavity.back();     		// Last cavity
			
			// Final boundary
			double r1 = model->modules[0].getBoundary()->getRefCoeff();
			double r2 = model->modules.back().getBoundary()->getRefCoeff();
			double x1 = model->modules[indx_k].getPosition1();
			double ni = model->modules[indx_k].getRefInd();
			double n0 = model->modules.back().getBoundary()->getNextCavityIndex();

			// Reflected from cavity
			model->modules[indx_k].getTwoArmCavity()->getEfp_back_wall(cavity_trans_E_pl);
			model->modules.back().getBoundary()->getEminus(cavity_trans_E_mi);

			// Transmission through boundary for output
			std::complex<double> *transmitted = model->modules.back().getBoundary()->setEpluss();
			std::complex<double> scale = sqrt((ni/n0)*(1.0-r2));
			memcpy(transmitted, cavity_trans_E_pl, VECSEL_transverse_points_number*sizeof(std::complex<double>));
			cblas_zscal(VECSEL_transverse_points_number, &scale, transmitted, 1);
			model->modules[indx_k].getTwoArmCavity()->getEfm_back_wall(cavity_trans_E_pl);
			cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_pl , 1, transmitted, 1);
		
			// Feedback into cavity
			std::complex<double> *feedback = model->modules[indx_k].getTwoArmCavity()->setEbm();
			#ifdef CYCLIC_BOUNDARIES
				memcpy(feedback, cavity_trans_E_mi, VECSEL_transverse_points_number*sizeof(std::complex<double>));
				model->modules[indx].getTwoArmCavity()->getEbm_front_wall(cavity_trans_E_pl);
				scale = -sqrt(r1);
				cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_pl , 1, feedback, 1);
			#else
				scale = -sqrt(r2);
				memcpy(feedback, cavity_trans_E_mi, VECSEL_transverse_points_number*sizeof(std::complex<double>));
				cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_pl , 1, feedback, 1);
			#endif
			
			feedback = model->modules[indx_k].getTwoArmCavity()->setEbp();
			#ifdef CYCLIC_BOUNDARIES
				memcpy(feedback, cavity_trans_E_mi, VECSEL_transverse_points_number*sizeof(std::complex<double>));
				model->modules[indx].getTwoArmCavity()->getEbp_front_wall(cavity_trans_E_mi);
				scale = -sqrt(r1);
				cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_mi , 1, feedback, 1);
			#else
				scale = -sqrt(r2);
				model->modules[indx_k].getTwoArmCavity()->getEfp_back_wall(cavity_trans_E_pl);
				memcpy(feedback, cavity_trans_E_mi, VECSEL_transverse_points_number*sizeof(std::complex<double>));
				cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_pl , 1, feedback, 1);
			#endif
		
		#else
			int indx   = model->quick_index_cavity[0];			// Current cavity
			int indx_k = model->quick_index_cavity.back();     		// Last cavity
			
			// Final boundary
			double r1 = model->modules[0].getBoundary()->getRefCoeff();
			double r2 = model->modules.back().getBoundary()->getRefCoeff();
			double x1 = model->modules[indx_k].getPosition1();
			double ni = model->modules[indx_k].getRefInd();
			double n0 = model->modules.back().getBoundary()->getNextCavityIndex();

			// Reflected from cavity
			model->modules[indx_k].getCavity()->getEpluss_right_wall(cavity_trans_E_pl);
			model->modules.back().getBoundary()->getEminus(cavity_trans_E_mi);

			// Transmission through boundary for output
			std::complex<double> *transmitted = model->modules.back().getBoundary()->setEpluss();
			std::complex<double> scale = sqrt((ni/n0)*(1.0-r2));
			memcpy(transmitted, cavity_trans_E_pl, VECSEL_transverse_points_number*sizeof(std::complex<double>));
			cblas_zscal(VECSEL_transverse_points_number, &scale, transmitted, 1);
		
			// Feedback into cavity
			std::complex<double> *feedback = model->modules[indx_k].getCavity()->setEminus();
			#ifdef CYCLIC_BOUNDARIES
				memcpy(feedback, cavity_trans_E_mi, VECSEL_transverse_points_number*sizeof(std::complex<double>));
				model->modules[indx].getCavity()->getEminus_left_wall(cavity_trans_E_mi);
				scale = -sqrt(r1);
				cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_mi , 1, feedback, 1);
			#else
				scale = -sqrt(r2);
				memcpy(feedback, cavity_trans_E_mi, VECSEL_transverse_points_number*sizeof(std::complex<double>));
				cblas_zaxpy(VECSEL_transverse_points_number, &scale, cavity_trans_E_pl , 1, feedback, 1);
			#endif
		#endif	
	}

	if( model->quick_index_cavity_noQW.size() > 0 )
	{
		// Iterate over all (noQW) cavities
		#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
		for(int i=0; i<model->quick_index_cavity_noQW.size()-1; i++)
		{
			int indx   = quick_index_cavity_noQW[i];	// Current cavity
			int indx_k = indx+1;     			// Next cavity
		
			if(modules[indx_k].isTwoArmInterface())
			{
				int indx_end = modules[indx_k].getTwoArmInterface()->getPostCav();
				#ifdef TRANS_DELAY
					model->iterateModules_updateSingleSurface_TwoArm_transfer_matrix_interface_delay(modules[indx].getCavity(),modules[indx_k].getTwoArmInterface(),modules[indx_end].getCavity(), modules[indx_k+1].getTwoArmCavity()); //Update interface with transversal delay
				#else
					model->iterateModules_updateSingleSurface_TwoArm_transfer_matrix_interface(modules[indx].getCavity(),modules[indx_k].getTwoArmInterface(),modules[indx_end].getCavity(), modules[indx_k+1].getTwoArmCavity()); //Update interface without any pulse delay
				#endif
				std::complex<double> tmp_reflect=modules[indx_k].getTwoArmInterface()->getReflect(); 	
				model->iterateModules_updateSingleSurface_TwoArm_transfer_matrix_back(modules[indx_end-1].getTwoArmCavity(),tmp_reflect);
			}else
			{
				model->iterateModules_updateSingleSurface_transfer_matrix_noQW(modules[indx].getCavity(), modules[indx_k].getCavity());
			}
		}
	}
	if( model->quick_index_twoArmCavity_noQW.size() > 0 )
	{
		// Iterate over all (noQW) two arm cavities
		#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
		for(int i=0; i<model->quick_index_twoArmCavity_noQW.size()-1; i++)
		{
			int indx   = quick_index_twoArmCavity_noQW[i];	// Current cavity
			int indx_k = indx+1;			     	// Next cavity
			if(modules[indx_k].isTwoArmCavity())
			{
				model->iterateModules_updateSingleSurface_TwoArm_transfer_matrix_noQW(modules[indx].getTwoArmCavity(), modules[indx_k].getTwoArmCavity());
			}	
		}
	}
	
	if( model->quick_index_birefringentCrystal.size() > 0 )
	{
		// Iterate over all (noQW) two arm cavities
		#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
		for(int i=0; i<model->quick_index_birefringentCrystal.size(); i++)
		{
			int indx   = quick_index_birefringentCrystal[i];	// Current cavity
			int indx_mk = indx-1;			     		// Next cavity
			int indx_pk = indx+1;
			model->iterateModules_updateSingleSurface_TwoArm_transfer_matrix_birefringentCrystal_front(modules[indx_mk].getTwoArmCavity(), modules[indx].getBirefringentCrystal());
			model->iterateModules_updateSingleSurface_TwoArm_transfer_matrix_birefringentCrystal_back(modules[indx].getBirefringentCrystal(), modules[indx_pk].getTwoArmCavity());
		}
	}
}

void VECSEL::iterateModules(double t_sim, double DT)
{
	/* ===========================================
	 * SOLVE MAXWELL's EQUATIONS
	 * */
	if (MPI_MY_RANK == 0)
	{
		#ifdef MPI_BALANCE_WORKLOAD
		MPI_load->start();
		#endif

		#ifdef USE_MAIN_TIMERS
		MainStat->start("VECSEL::transfer matrix");
		#endif
//					// Normal start conditions on boundaries
			if (VECSEL_pulse_start_l == 1)
			{
				std::complex<double> *tmp = modules[0].getBoundary()->setEpluss();
				maxwell_initial_E(t_sim,tmp);
			}

			if (VECSEL_pulse_start_r == 1)
			{
				std::complex<double> *tmp = modules.back().getBoundary()->setEminus();
				maxwell_initial_E(t_sim,tmp);
			}



		#ifdef USE_MAIN_TIMERS
		MainStat->start("TM Update Storage");
		#endif
		#pragma omp parallel num_threads(OMP_THREADS_LEVEL_1)
		{
			for(unsigned i=0; i<quick_index_cavity_lens.size(); i++)	
			{	
				#pragma omp single nowait
				{
					int indx = quick_index_cavity_lens[i];
					modules[indx].getCavity()->updateStorage_lens_pluss();
				}
				#pragma omp single nowait
				{
					int indx = quick_index_cavity_lens[i];
					modules[indx].getCavity()->updateStorage_lens_minus();
				}
			}
	
			for(unsigned i=0; i<quick_index_cavity_lens_halfCav.size(); i++)
			{
				#pragma omp single nowait
				{
					int indx = quick_index_cavity_lens_halfCav[i];
					modules[indx].getCavity()->updateStorage_lens_pluss();
				}
				#pragma omp single nowait
				{
					int indx = quick_index_cavity_lens_halfCav[i];
					modules[indx].getCavity()->updateStorage_lens_minus();
				}
			}
	
			#pragma omp for
			for(unsigned i = 0; i < quick_index_cavity_freeSpace.size(); i++)
			{
				int indx = quick_index_cavity_freeSpace[i];
				modules[indx].getCavity()->updateStorage_freeSpace_forcedBPM(); // Update storage of fields
			}
			#pragma omp for
			for(unsigned i = 0; i < quick_index_cavity_noBPM.size(); i++)
			{
				int indx = quick_index_cavity_noBPM[i];
				modules[indx].getCavity()->updateStorage_freeSpace(); // Update storage of fields
			}
			#pragma omp for
			for(unsigned i = 0; i < quick_index_twoArmCavity.size(); i++)
			{
				int indx = quick_index_twoArmCavity[i];
				modules[indx].getTwoArmCavity()->updateStorage_freeSpace(); // Update storage of fields
			}
			#pragma omp for
			for(unsigned i = 0; i < quick_index_birefringentCrystal.size(); i++)
			{
				int indx = quick_index_birefringentCrystal[i];
				modules[indx].getBirefringentCrystal()->updateStorage_freeSpace(); // Update storage of fields
			}
			#pragma omp for
			for(unsigned i = 0; i < quick_index_twoArmInterface.size(); i++)
			{
				int indx = quick_index_twoArmInterface[i];
				modules[indx].getTwoArmInterface()->updateStorage_freeSpace(); // Update storage of fields
			}
		}

		#ifdef USE_MAIN_TIMERS
		MainStat->stop("TM Update Storage");
		#endif
		

		#ifdef USE_MAIN_TIMERS
		MainStat->start("TM B.C.");
		#endif
		//============
		// Update Cav
		//============
		iterateModules_updateAllSurface(this,t_sim);

		#ifdef USE_MAIN_TIMERS
		MainStat->stop("TM B.C.");
		#endif
		

		#ifdef USE_MAIN_TIMERS
		MainStat->stop("VECSEL::transfer matrix");
		#endif

		#ifdef MPI_BALANCE_WORKLOAD
		MPI_load->stop();
		#endif


		//=============================================
		// Code moved from iterateModules_updateAllSurface
		// Collect macroscopic polarizations
		#ifdef USE_MAIN_TIMERS
		MainStat->start("VECSEL::MPI_comm_gath");
		#endif

		#ifdef ITERATE_QW	
		MPI_Gatherv( MPI_WORK_DIST_P_LOCAL, MPI_WORK_DIST_P_SIZE[MPI_MY_RANK],MPI_DOUBLE_COMPLEX, // Where I store my stuff
				 MPI_WORK_DIST_P_GLOBAL[0], MPI_WORK_DIST_P_SIZE, MPI_WORK_DIST_P_OFFSET, MPI_DOUBLE_COMPLEX, // Distribution of work
				0, *MPI_WORK_GANG); // Who is recieving the data
		#endif

		#ifdef USE_MAIN_TIMERS
		MainStat->stop("VECSEL::MPI_comm_gath");
		#endif


		#ifdef USE_MAIN_TIMERS
		MainStat->start("VECSEL::prepare MPI");
		#endif


		// unscramble
		#ifdef ITERATE_QW
		#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
		for(int i = 0; i < quick_index_totalDevice.size(); i++)
		{
			MPI_LoadBalancer_P_tmp[4*MPI_LoadBalancer_index_set[i]] = MPI_WORK_DIST_P_GLOBAL[0][4*i];
			MPI_LoadBalancer_P_tmp[1+4*MPI_LoadBalancer_index_set[i]] = MPI_WORK_DIST_P_GLOBAL[0][1+4*i];
			MPI_LoadBalancer_P_tmp[2+4*MPI_LoadBalancer_index_set[i]] = MPI_WORK_DIST_P_GLOBAL[0][2+4*i];
			MPI_LoadBalancer_P_tmp[3+4*MPI_LoadBalancer_index_set[i]] = MPI_WORK_DIST_P_GLOBAL[0][3+4*i];
		}
		#endif		

		
		#pragma omp parallel num_threads(OMP_THREADS_LEVEL_1)
		{
			// Set all MacPol into correct cavities
			#ifdef ITERATE_QW	
			#pragma omp for 
			for(int j = 0; j < quick_index_totalDevice.size(); j = j + VECSEL_transverse_points_number)
			{
				int indx_k = quick_index_device_previous_cavity[j];
				if (modules[indx_k].isCavity())
				{
					modules[indx_k].getCavity()->set_transfer_matrix_macPol(&(MPI_LoadBalancer_P_tmp[4*j]));	
				} else
				{
					modules[indx_k].getTwoArmCavity()->set_transfer_matrix_macPol_fp(&(MPI_LoadBalancer_P_tmp[4*j]));
					modules[indx_k].getTwoArmCavity()->set_transfer_matrix_macPol_fm(&(MPI_LoadBalancer_P_tmp[4*j+1]));
					modules[indx_k].getTwoArmCavity()->set_transfer_matrix_macPol_bp(&(MPI_LoadBalancer_P_tmp[4*j+2]));
					modules[indx_k].getTwoArmCavity()->set_transfer_matrix_macPol_bm(&(MPI_LoadBalancer_P_tmp[4*j+3]));
				}
			}
			
			#endif
			
			// Iterate final parts of transfer matrix for cavities with QWs
			if( quick_index_cavity_QW.size() > 0)
			{
				#pragma omp for 
				for(int i=0; i<quick_index_cavity_QW.size()-1; i = i + 2)
				{
					int indx   = quick_index_cavity_QW[i];		// Current cavity
					int indx_k = quick_index_cavity_QW[i+1];     	// Next cavity
					iterateModules_updateSingleSurface_transfer_matrix(modules[indx].getCavity(), modules[indx_k].getCavity());
				}
			}
	
			// Iterate final parts of transfer matrix for two arm cavities with QWs
			if( quick_index_twoArmCavity_QW.size() > 0)
			{	
				#pragma omp for
				for(int i=0; i<quick_index_twoArmCavity_QW.size()-1; i = i + 2)
				{
					int indx   = quick_index_twoArmCavity_QW[i];		// Current cavity
					int indx_k = quick_index_twoArmCavity_QW[i+1];     	// Next cavity
					iterateModules_updateSingleSurface_TwoArm_transfer_matrix(modules[indx].getTwoArmCavity(), modules[indx_k].getTwoArmCavity());
				}
			}
		}
		
		#ifdef ITERATE_QW	
		
		#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
		for(int i = 0; i < quick_index_device_previous_cavity.size(); i = i + VECSEL_transverse_points_number)
		{
			
			int indx_k = quick_index_device_previous_cavity[i];
			if(modules[indx_k].isCavity())
			{
				modules[indx_k].getCavity()->evaluateEprop_x1_fast(    &MPI_LoadBalancer_E_tmp[8*i  ], 8);
				modules[indx_k].getCavity()->evaluateEprop_x1_tp1_fast(&MPI_LoadBalancer_E_tmp[8*i+1], 8);
			}
			else
			{
				//WTFF-modules[indx_k].getTwoArmCavity()->evaluateEprop_back_wall(    &MPI_LoadBalancer_E_tmp[8*i]);
				modules[indx_k+VECSEL_transverse_points_number+1].getTwoArmCavity()->evaluateEprop_front_wall(    &MPI_LoadBalancer_E_tmp[8*i]);
			}
		}

		// scramble
		#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
		for(int i = 0; i < quick_index_totalDevice.size(); i++)
		{
			for(int j=0; j < 8; j++)
			{
				MPI_WORK_DIST_E_GLOBAL[8*i+j] = MPI_LoadBalancer_E_tmp[8*MPI_LoadBalancer_index_set[i]+j];
			}
		}
		

		// Update the storage of polarization vectors
		std::complex<double> *dummy = MPI_WORK_DIST_P_GLOBAL[2];
		MPI_WORK_DIST_P_GLOBAL[2] = MPI_WORK_DIST_P_GLOBAL[1];
		MPI_WORK_DIST_P_GLOBAL[1] = MPI_WORK_DIST_P_GLOBAL[0];
		MPI_WORK_DIST_P_GLOBAL[0] = dummy; // Will be overwritten by QWs, does not need to be zeroed
		#endif


		#ifdef USE_MAIN_TIMERS
		MainStat->stop("VECSEL::prepare MPI");
		#endif
	} else {

		// WORKERS WAIT HERE
		#ifdef USE_MAIN_TIMERS
		MainStat->start("VECSEL::MPI_comm_gath");
		#endif

		#ifdef ITERATE_QW	
		MPI_Gatherv( MPI_WORK_DIST_P_LOCAL, MPI_WORK_DIST_P_SIZE[MPI_MY_RANK],MPI_DOUBLE_COMPLEX, // Where I store my stuff
				 MPI_WORK_DIST_P_GLOBAL[0], MPI_WORK_DIST_P_SIZE, MPI_WORK_DIST_P_OFFSET, MPI_DOUBLE_COMPLEX, // Distribution of work
				0, *MPI_WORK_GANG); // Who is recieving the data
		#endif

		#ifdef USE_MAIN_TIMERS
		MainStat->stop("VECSEL::MPI_comm_gath");
		#endif

	}
	

    /* ===========================================
     * SOLVE SBE
     * */
	
	#ifdef ITERATE_QW

//	MPI_Barrier(*MPI_WORK_GANG);

/*
	if (MPI_MY_RANK==0)
	{
		for(int i = 0; i < MPI_WORK_DIST_TOTAL; i++)
		{
			MPI_WORK_DIST_E_GLOBAL[i] = std::complex<double>(i+1,i+1);
		}
	}

	//=====
	MPI_Comm new_world;
        MPI_Comm_dup(*MPI_WORK_GANG,&new_world);
        for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
        {
                MPI_Barrier(new_world);
                if (i == MPI_MY_RANK)
                {
                        cout << "Process[" << i << "]: a)" << endl;
                }
        }
        MPI_Barrier(new_world);
*/
	//=====
	#ifdef USE_MAIN_TIMERS
	MainStat->start("VECSEL::MPI_comm_scatt");
	#endif

	// Send E(t) at the given QW to workers
	MPI_Scatterv(MPI_WORK_DIST_E_GLOBAL, MPI_WORK_DIST_E_SIZE, MPI_WORK_DIST_E_OFFSET, MPI_DOUBLE_COMPLEX, // Distribution of work
			MPI_WORK_DIST_E_LOCAL, MPI_WORK_DIST_E_SIZE[MPI_MY_RANK],MPI_DOUBLE_COMPLEX, // Where I store my stuff
			0, *MPI_WORK_GANG); // Who is sending the data
	#ifdef USE_MAIN_TIMERS
	MainStat->stop("VECSEL::MPI_comm_scatt");
	#endif
/*
	for(int i = 0; i < 2*MPI_WORK_DIST_TOTAL; i++)
	{
		MPI_WORK_DIST_E_LOCAL[i] = MPI_WORK_DIST_E_GLOBAL[i];
	}
*/
	//=====
/*
        for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
        {
                MPI_Barrier(new_world);
                if (i == MPI_MY_RANK)
                {
                        cout << "Process[" << i << "]: b) E_local = ";
			for(int j = 0; j < MPI_WORK_DIST_E_SIZE[MPI_MY_RANK]; j++)
			{
				cout << MPI_WORK_DIST_E_LOCAL[j] << " ";
			}
			cout << endl;
                }
        }
        MPI_Barrier(new_world);
*/
	//=====
	#ifdef USE_MAIN_TIMERS
	MainStat->start("VECSEL::compute all QWs");
	#endif

	// Compute P(t) for my own QWs using MPI_WORK_DIST_E_LOCAL
	std::complex<double> *openmp_MPI_WORK_DIST_E_LOCAL = MPI_WORK_DIST_E_LOCAL;
	int **openmp_MPI_WORK_DIST = MPI_WORK_DIST;
	std::complex<double> *openmp_MPI_WORK_DIST_P_LOCAL = MPI_WORK_DIST_P_LOCAL;
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1) shared(openmp_MPI_WORK_DIST_E_LOCAL, openmp_MPI_WORK_DIST) firstprivate(openmp_MPI_WORK_DIST_P_LOCAL)
	for(int j = openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < openmp_MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		// do work inside QW 
		//#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_2)
		//for(int j = 0; j < number_K_points; j++)
		


		int indx = quick_index_totalDevice[j]; // Index of device
		if(modules[indx].isDevice())
		{	
			std::complex<double> E_prop     = openmp_MPI_WORK_DIST_E_LOCAL[  8*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))]; // Set propagating field at device
			std::complex<double> E_prop_tp1 = openmp_MPI_WORK_DIST_E_LOCAL[1+8*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))];

			Device *dev0 = modules[indx].getDevice();

			// Electric field in device
			double focusElectricField = dev0->getFocusE();
			dev0->setElectricField(E_prop*focusElectricField);

			std::complex<double> E_prop_tp05 = 0.5*(E_prop + E_prop_tp1);
			dev0->setElectricField_tp1(E_prop_tp1*focusElectricField);
			dev0->setElectricField_tp05(E_prop_tp05*focusElectricField);
		
			// Iterate device
			dev0->sbe_iterate(t_sim, DT);

		

			double focusE 	= dev0->getFocusE();
			double eff_qw 	= dev0->getEffectiveQW();
			openmp_MPI_WORK_DIST_P_LOCAL[4*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))] = VECSEL_QW_FEEDBACK*dev0->getMacroscopicPolarization()*(eff_qw/focusE);
		} else
		{
			std::complex<double> E_prop_fp     = openmp_MPI_WORK_DIST_E_LOCAL[  8*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))]; // Set propagating field at device
			std::complex<double> E_prop_fp_tp1 = openmp_MPI_WORK_DIST_E_LOCAL[1+8*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))];
			std::complex<double> E_prop_fm     = openmp_MPI_WORK_DIST_E_LOCAL[2+8*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))]; // Set propagating field at device
			std::complex<double> E_prop_fm_tp1 = openmp_MPI_WORK_DIST_E_LOCAL[3+8*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))];
			std::complex<double> E_prop_bp     = openmp_MPI_WORK_DIST_E_LOCAL[4+8*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))]; // Set propagating field at device
			std::complex<double> E_prop_bp_tp1 = openmp_MPI_WORK_DIST_E_LOCAL[5+8*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))];
			std::complex<double> E_prop_bm     = openmp_MPI_WORK_DIST_E_LOCAL[6+8*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))]; // Set propagating field at device
			std::complex<double> E_prop_bm_tp1 = openmp_MPI_WORK_DIST_E_LOCAL[7+8*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))];
	
			std::complex<double> E_prop_fp_tp05 = 0.5*(E_prop_fp + E_prop_fp_tp1);
			std::complex<double> E_prop_fm_tp05 = 0.5*(E_prop_fm + E_prop_fm_tp1);
			std::complex<double> E_prop_bp_tp05 = 0.5*(E_prop_bp + E_prop_bp_tp1);
			std::complex<double> E_prop_bm_tp05 = 0.5*(E_prop_bm + E_prop_bm_tp1);
	
			TwoArmDevice *dev0 = modules[indx].getTwoArmDevice();
			// Electric field in device
			double focusElectricField = dev0->getFocusE();
			dev0->setElectricField_fp(E_prop_fp*focusElectricField);
			dev0->setElectricField_fm(E_prop_fm*focusElectricField);
			dev0->setElectricField_bp(E_prop_bp*focusElectricField);
			dev0->setElectricField_bm(E_prop_bm*focusElectricField);
			
			dev0->setElectricField_fp_tp1(E_prop_fp_tp1*focusElectricField);
			dev0->setElectricField_fm_tp1(E_prop_fm_tp1*focusElectricField);
			dev0->setElectricField_bp_tp1(E_prop_bp_tp1*focusElectricField);
			dev0->setElectricField_bm_tp1(E_prop_bm_tp1*focusElectricField);
				
			dev0->setElectricField_fp_tp05(E_prop_fp_tp05*focusElectricField);
			dev0->setElectricField_fm_tp05(E_prop_fm_tp05*focusElectricField);
			dev0->setElectricField_bp_tp05(E_prop_bp_tp05*focusElectricField);
			dev0->setElectricField_bm_tp05(E_prop_bm_tp05*focusElectricField);	
				
			// Iterate device
			dev0->sbe_iterate(t_sim, DT);

			double eff_qw 	= dev0->getEffectiveQW();
			openmp_MPI_WORK_DIST_P_LOCAL[4*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))] = VECSEL_QW_FEEDBACK*dev0->getMacroscopicPolarization_fp1()*(eff_qw/focusElectricField);
			openmp_MPI_WORK_DIST_P_LOCAL[1+4*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))] = VECSEL_QW_FEEDBACK*dev0->getMacroscopicPolarization_fm1()*(eff_qw/focusElectricField);
			openmp_MPI_WORK_DIST_P_LOCAL[2+4*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))] = VECSEL_QW_FEEDBACK*dev0->getMacroscopicPolarization_bp1()*(eff_qw/focusElectricField);
			openmp_MPI_WORK_DIST_P_LOCAL[3+4*(j-(openmp_MPI_WORK_DIST[MPI_MY_RANK][0]-1))] = VECSEL_QW_FEEDBACK*dev0->getMacroscopicPolarization_bm1()*(eff_qw/focusElectricField);	
		}
	}

	#ifdef USE_MAIN_TIMERS
	MainStat->stop("VECSEL::compute all QWs");
	#endif

	//=====
/*
        for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
        {
                MPI_Barrier(new_world);
                if (i == MPI_MY_RANK)
                {
                        cout << "Process[" << i << "]: d) P_local";
			for(int j = 0; j < MPI_WORK_DIST_P_SIZE[MPI_MY_RANK]; j++)
			{
				cout << MPI_WORK_DIST_P_LOCAL[j] << " ";
			}
			cout << endl;
                }
        }
        MPI_Barrier(new_world);
*/
	//=====

	// Send P(t) back to master
/*
	if (MPI_MY_RANK>0)
	{
		#ifdef USE_MAIN_TIMERS
		MainStat->start("VECSEL::MPI_comm_gath");
		#endif
		MPI_Gatherv( MPI_WORK_DIST_P_LOCAL, MPI_WORK_DIST_P_SIZE[MPI_MY_RANK],MPI_DOUBLE_COMPLEX, // Where I store my stuff
				 MPI_WORK_DIST_P_GLOBAL[0], MPI_WORK_DIST_P_SIZE, MPI_WORK_DIST_P_OFFSET, MPI_DOUBLE_COMPLEX, // Distribution of work
				0, *MPI_WORK_GANG); // Who is recieving the data
		#ifdef USE_MAIN_TIMERS
		MainStat->stop("VECSEL::MPI_comm_gath");
		#endif
	}
*/
	
/*
	for(int i = 0; i < MPI_WORK_DIST_TOTAL; i++)
	{
		MPI_WORK_DIST_P_GLOBAL[0][i] = MPI_WORK_DIST_P_LOCAL[i];
	}
*/	
	//=====
/*
	if (MPI_MY_RANK==0)
	{
		for(int j = 0; j < 3; j++)
		{
			cout << "P_global["<< j <<"][...] = ";
			for(int i = 0; i < MPI_WORK_DIST_TOTAL; i++)
			{
				cout << MPI_WORK_DIST_P_GLOBAL[j][i] << " ";
			}
			cout << endl;
		}
		
	}
        for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
        {
                MPI_Barrier(new_world);
                if (i == MPI_MY_RANK)
                {
                        cout << "Process[" << i << "]: e)" << endl;
                }
        }

        MPI_Barrier(new_world);
	exit(-1);
*/
	
	
	
	#endif
}


//===============================
// Maxwell functions

/* Initialize all simulation variables
 * DT is the requested time step for the simulation
 * Calles initialize on all objects
 * */
void VECSEL::maxwell_initialize(double DT, int mpi_rank)
{
	#ifdef MPI_BALANCE_WORKLOAD
	MPI_load = new myTimer("Maxwell timer");
	#endif
	
	mpi_initialize_rank(mpi_rank);
	
	if (MPI_MY_RANK==0)
	{
		// Initialize Transverse Dimension and output
		output_E_real = new std::ofstream[VECSEL_transverse_points_number];
		output_E_imag = new std::ofstream[VECSEL_transverse_points_number];
		output_back_E_real = new std::ofstream[VECSEL_transverse_points_number];
		output_back_E_imag = new std::ofstream[VECSEL_transverse_points_number];
		output_Cav_Snapshot_E_real = new std::ofstream[VECSEL_transverse_points_number];
		output_Cav_Snapshot_E_imag = new std::ofstream[VECSEL_transverse_points_number];
		
		// Set initial pulse transverse profile
		if (VECSEL_initial_transverse_pulse_profile ==  NULL)
		{
			VECSEL_initial_transverse_pulse_profile = new std::complex<double>[VECSEL_transverse_points_number];
		}
		double waist0 = (1.0/sqrt(2.0))*VECSEL_initial_transverse_FWHM/sqrt(2.0*log(2.0)); // For Gaussian shape with beam waist, fwhm of |E(t)|
		//double waist0 = VECSEL_initial_transverse_FWHM/sqrt(2.0*log(2.0)); // For Gaussian shape with beam waist, fwhm of |E(t)|^2
		for(int i = 0; i < VECSEL_transverse_points_number; i++)
		{
			VECSEL_initial_transverse_pulse_profile[i] = exp(-VECSEL_transverse_points_y[i]*VECSEL_transverse_points_y[i]/(waist0*waist0));
		}

		// Update delays in Cavities
		double maximal_index = 0;
		double minTime = 1.0;
		double totalTime = 0;
		bool toLargeDT = false;
		
		VECSEL_DT = DT;
			int ext_rtt_steps = 0; // Number of timesteps in a round-trip  outside of the filters averaging interface length
			int ext_rtt_steps2 = 0; // Number of timesteps in a round-trip  outside of the filters for Initial condition

			// Boundaries
			modules[0].getBoundary()->initializeZero(VECSEL_transverse_points_number);
			modules.back().getBoundary()->initializeZero(VECSEL_transverse_points_number);

			// Cavities
			if(quick_index_cavity.size()>0)
			{
			for(unsigned i=0; i<quick_index_cavity.size(); i++)
			{
				int indx = quick_index_cavity[i];
				modules[indx].getCavity()->initializeZero(DT,VECSEL_transverse_points_number, VECSEL_transverse_points_y, VECSEL_transverse_points_R_max, VECSEL_transverse_points_boundary_guard_ratio); // Simulation variables
				//Find next cavity's refractive index
				double leftRefInd = modules[indx].getCavity()->getRefInd();
				double leftWidth  = modules[indx].getCavity()->getWidth();
				
				if (leftRefInd > maximal_index)
				{
					maximal_index = leftRefInd;
				}
			
				// Get total length
				totalTime += leftWidth/(c0/leftRefInd);
			
				// Get minimal width
				if (leftWidth/(c0/leftRefInd) < minTime)
				{
					minTime = leftWidth/(c0/leftRefInd);
				}
				ext_rtt_steps += 2*(modules[indx].getCavity()->getNumberOfTimesteps()-1); // CORRECT
				ext_rtt_steps2 += 2*(modules[indx].getCavity()->getNumberOfTimesteps()-2); // CORRECT
				
			}
			}
		
		if(quick_index_twoArmInterface.size()>0)
			{	
			for(unsigned i=0; i<quick_index_twoArmInterface.size(); i++)
			{
				int indx = quick_index_twoArmInterface[i];
				modules[indx].getTwoArmInterface()->initializeZero(DT,VECSEL_transverse_points_number, VECSEL_transverse_points_y, VECSEL_transverse_points_R_max, VECSEL_transverse_points_boundary_guard_ratio); // Simulation variables
			
				//Find next cavity's refractive index
				double leftRefInd = modules[indx].getTwoArmInterface()->getRefInd();
				double leftWidth  = modules[indx].getTwoArmInterface()->getWidth();
				
				if (leftRefInd > maximal_index)
				{
					maximal_index = leftRefInd;
				}
			
				// Get total length
				totalTime += leftWidth/(c0/leftRefInd);
			
				// Get minimal width
				if (leftWidth/(c0/leftRefInd) < minTime)
				{
					minTime = leftWidth/(c0/leftRefInd);
				}
				ext_rtt_steps += 2*(modules[indx].getTwoArmInterface()->getNumberOfTimesteps()-1); // CORRECT
				ext_rtt_steps2 += 2*(modules[indx].getTwoArmInterface()->getNumberOfTimesteps()-2); // CORRECT
				
			}
			}		

			if(quick_index_twoArmCavity.size()>0)
			{
			for(unsigned i=0; i<quick_index_twoArmCavity.size(); i++)
			{
				int indx = quick_index_twoArmCavity[i];
		// Simulation variables
				double leftRefInd=0.0;
				double leftWidth=0.0;
				//Find next cavity's refractive index
				modules[indx].getTwoArmCavity()->initializeZero(DT,VECSEL_transverse_points_number, VECSEL_transverse_points_y, VECSEL_transverse_points_R_max, VECSEL_transverse_points_boundary_guard_ratio); 
				leftRefInd = modules[indx].getTwoArmCavity()->getRefInd();
				leftWidth  = modules[indx].getTwoArmCavity()->getWidth();			
				ext_rtt_steps += 2*(modules[indx].getTwoArmCavity()->getNumberOfTimesteps()-1); // CORRECT
				ext_rtt_steps2 += 2*(modules[indx].getTwoArmCavity()->getNumberOfTimesteps()-2); // CORRECT
				
				if (leftRefInd > maximal_index)
				{
					maximal_index = leftRefInd;
				}
			
				// Get total length. Doubled accounting for correct number of passes
				totalTime += 2.0*leftWidth/(c0/leftRefInd);
			
				// Get minimal width
				if (leftWidth/(c0/leftRefInd) < minTime)
				{
					minTime = leftWidth/(c0/leftRefInd);
				}	
			}
			}
			if(quick_index_birefringentCrystal.size()>0)
			{
			for(unsigned i=0; i<quick_index_birefringentCrystal.size(); i++)
			{
				int indx = quick_index_birefringentCrystal[i];
		// Simulation variables
				double leftRefInd=0.0;
				double leftWidth=0.0;
				//Find next cavity's refractive index
				modules[indx].getBirefringentCrystal()->initializeZero(DT,VECSEL_transverse_points_number, VECSEL_transverse_points_y, VECSEL_transverse_points_R_max, VECSEL_transverse_points_boundary_guard_ratio); 
				leftRefInd = modules[indx].getBirefringentCrystal()->getRefInd();
				leftWidth  = modules[indx].getBirefringentCrystal()->getWidth();			
				ext_rtt_steps += 2*(modules[indx].getBirefringentCrystal()->getNumberOfTimesteps()-1); // CORRECT
				ext_rtt_steps2 += 2*(modules[indx].getBirefringentCrystal()->getNumberOfTimesteps()-2); // CORRECT
				
				if (leftRefInd > maximal_index)
				{
					maximal_index = leftRefInd;
				}
			
				// Get total length. Doubled accounting for correct number of passes
				totalTime += 2.0*leftWidth/(c0/leftRefInd);
			
				// Get minimal width
				if (leftWidth/(c0/leftRefInd) < minTime)
				{
					minTime = leftWidth/(c0/leftRefInd);
				}	
			}
			}
			#ifndef CYCLIC_BOUNDARIES	
				VECSEL_ROUND_TRIP_TIME = 2.0*totalTime;
			#endif
		
			// Set FILTER DELAY
			int filter_rtt_steps = 0;
			for(unsigned i=0; i<modules.size(); i++)
			{
				if (modules[i].isFilter())
				{
					modules[i].getFilter()->setFilter_pluss_active_steps(ext_rtt_steps);
					VECSEL_ROUND_TRIP_TIME += (modules[i].getFilter()->getFilterLength()-1)*VECSEL_DT;
					filter_rtt_steps += (modules[i].getFilter()->getFilterLength()-1); // Subtract 2 to avoid double counting
				}
			}

			VECSEL_ROUND_TRIP_ITERATIONS = ext_rtt_steps2 + filter_rtt_steps;
		

			//cout << "STRUCTURE: LINEAR CAVITY" << endl;
			cout << " Round trip time        = " << VECSEL_ROUND_TRIP_TIME/ps << " [ps]" << endl;
			cout << " Round trip #timesteps  = " << ext_rtt_steps + filter_rtt_steps << " ("<< ext_rtt_steps <<", "<< filter_rtt_steps <<")" << endl;
			cout << " Largest timestep       = " << (minTime)/fs << " [fs]"<< endl;

			// Output round trip time to file
			std::stringstream fileName;
			fileName << getToFileOutputKey() << "round_trip_time" << ".dat";
			saveBinary(fileName.str(), &VECSEL_ROUND_TRIP_TIME, 1);



			double cos_th_k_left, cos_th_k_right;
			if(quick_index_cavity.size()>0)
			{
			// Set up transfer matrix
			for(unsigned i=0; i<quick_index_cavity.size()-1; i++)
			{
				int indx   = quick_index_cavity[i];
				std::complex<double> nk    = 1.0;
				std::complex<double> T_Em  = 1.0;
				if(modules[indx+1].isTwoArmInterface())
				{
				int indx_k=indx+1;
				nk = modules[indx_k].getTwoArmInterface()->getRefInd() + I*modules[indx_k].getTwoArmInterface()->getRefInd_im();
				modules[indx_k].getCosTh(&cos_th_k_left, &cos_th_k_right);
				nk = nk*cos_th_k_left;

				T_Em = modules[indx_k].getTwoArmInterface()->get_transport_Eb_x0();
				} else
				{
				int indx_k = quick_index_cavity[i+1];
				nk = modules[indx_k].getCavity()->getRefInd() + I*modules[indx_k].getCavity()->getRefInd_im();
				modules[indx_k].getCosTh(&cos_th_k_left, &cos_th_k_right);
				nk = nk*cos_th_k_left;

				T_Em = modules[indx_k].getCavity()->get_transport_Em_x0();
				}

				double loss_p = 0.0;
				double loss_m = 0.0;
				if (modules[indx+1].isLossElement())
				{
					modules[indx+1].getLossElement()->getLossCoeff(&loss_p, &loss_m);
				}

				modules[indx].getCavity()->set_transfer_matrix(nk, T_Em, loss_p, loss_m);
				std::complex<double> a11,a12,a21,a22,*M;
				modules[indx].getCavity()->get_transfer_matrix(&a11,&a12,&a21,&a22,&M);
			}
			}
			if(quick_index_twoArmCavity.size()>0)
			{
			for(unsigned i=0; i<quick_index_twoArmCavity.size()-1; i++)
			{
				int indx   = quick_index_twoArmCavity[i];
				int indx_k = quick_index_twoArmCavity[i+1];

				std::complex<double> nk    = 1.0;
				std::complex<double> T_Em  = 1.0;
				if(!modules[indx+1].isCavity())
				{
				nk = modules[indx_k].getTwoArmCavity()->getRefInd() + I*modules[indx_k].getTwoArmCavity()->getRefInd_im();
				modules[indx_k].getCosTh(&cos_th_k_left, &cos_th_k_right);
				nk = nk*cos_th_k_left;
				
				T_Em = modules[indx_k].getTwoArmCavity()->get_transport_Eb_x0();
				}
					

				double loss_p = 0.0;
				double loss_m = 0.0;
				if (modules[indx+1].isLossElement())
				{
					modules[indx+1].getLossElement()->getLossCoeff(&loss_p, &loss_m);
				}

				
				modules[indx].getTwoArmCavity()->set_transfer_matrix(nk, T_Em, loss_p, loss_m);
				std::complex<double> a11,a12,a21,a22,*M;
				modules[indx].getTwoArmCavity()->get_transfer_matrix(&a11,&a12,&a21,&a22,&M, &M, &M, &M);
			}
			}
			if(quick_index_birefringentCrystal.size()>0)
			{
			for(unsigned i=0; i<quick_index_birefringentCrystal.size(); i++)
			{
				int indx   = quick_index_birefringentCrystal[i];
				int indx_pre = indx-1;
				int indx_post = indx+1;
			
				std::complex<double> n_post    = 1.0;
				std::complex<double> n_pre    = 1.0;
				std::complex<double> T_Em     = 1.0;
				double loss_p = 0.0;
				double loss_m = 0.0;
				
				T_Em = modules[indx_post].getTwoArmCavity()->get_transport_Eb_x0();
				n_pre = modules[indx_pre].getTwoArmCavity()->getRefInd() + I*modules[indx_pre].getTwoArmCavity()->getRefInd_im();
				n_post = modules[indx_post].getTwoArmCavity()->getRefInd() + I*modules[indx_pre].getTwoArmCavity()->getRefInd_im();
	
				modules[indx].getBirefringentCrystal()->set_transfer_matrix(n_post, T_Em, loss_p, loss_m);
				modules[indx].getBirefringentCrystal()->set_transfer_matrix_extraAxis(n_pre, n_post);
				std::complex<double> a11_ex,a12_ex,a21_ex,a22_ex,*M;
				modules[indx].getBirefringentCrystal()->get_transfer_matrix_extraAxis(&a11_ex,&a12_ex,&a21_ex,&a22_ex);
			}
			}
			if(quick_index_twoArmInterface.size()>0)
			{
			for(unsigned i=0; i<quick_index_twoArmInterface.size(); i++)
			{
				int indx   = quick_index_twoArmInterface[i];
				int indx_k = indx+1;
				int indx_end = quick_index_twoArmPostCav[i];			

				std::complex<double> nk = modules[indx_k].getTwoArmCavity()->getRefInd() + I*modules[indx_k].getTwoArmCavity()->getRefInd_im();
				modules[indx_k].getCosTh(&cos_th_k_left, &cos_th_k_right);
				nk = nk*cos_th_k_left;
				
				std::complex<double> T_Em = modules[indx_k].getTwoArmCavity()->get_transport_Eb_x0();
				std::complex<double> T_Ep_left = modules[indx-1].getCavity()->get_transport_Ep_x1();
				std::complex<double> T_Ep_right = modules[indx_end].getCavity()->get_transport_Em_x0();
				
				double loss_p = 0.0;
				double loss_m = 0.0;
				if (modules[indx+1].isLossElement())
				{
					modules[indx+1].getLossElement()->getLossCoeff(&loss_p, &loss_m);
				}
		
				modules[indx].getTwoArmInterface()->set_transfer_matrix(nk, T_Ep_left, T_Ep_right, T_Em, loss_p, loss_m);
				
				std::complex<double> a11_left,a12_left,a21_left,a22_left, a11_right, a12_right, a21_right, a22_right;
				modules[indx].getTwoArmInterface()->get_transfer_matrix(&a11_left,&a12_left,&a21_left,&a22_left, &a11_right, &a12_right, &a21_right, &a22_right);
			}
			}
		// Output central frequency to file.
		fileName.str("");
		fileName << getToFileOutputKey() << "w0.dat";
		double tmp = (2.0*Pi*c0)/getLambda();
		saveBinary(fileName.str(), &tmp, 1);

		// Output transverse grid to file
		fileName.str("");
		fileName << getToFileOutputKey() << "transverse_grid_y.dat";
		saveBinary(fileName.str(), VECSEL_transverse_points_y, VECSEL_transverse_points_number);

		double rL = modules.front().getBoundary()->getRefCoeff(); // Left boundary
		
		fileName.str("");
		fileName << getToFileOutputKey() << "reflection_left.dat";
		saveBinary(fileName.str(), &rL, 1);
			double rR = modules.back().getBoundary()->getRefCoeff(); // Right boundary
			fileName.str("");
			fileName << getToFileOutputKey() << "reflection_right" << ".dat";
			saveBinary(fileName.str(), &rR, 1);
			
			cout << "Intensity reflection  coeff [L/R] = " << rL << " / " << rR << endl;
			cout << "Intensity tansmission coeff [L/R] = " << 1.0-rL << " / " << 1.0-rR << endl;
			
//		} // end: VCAV
	
	} // if (MPI_MY_RANK==0)

	
	
	#ifdef ITERATE_QW	
	
	// Find total #QWs
	int total_work = quick_index_totalDevice.size();
	// 1. Init. schedule for all nodes
	mpi_initialize_work_arrays(total_work);
	
	// 2. Initialize QWs located in MPI_WORK_DIST[MPI_MY_RANK][0] = 1, MPI_WORK_DIST[MPI_MY_RANK][1] = #QWs
	// Serial	
		#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
		for(int i = MPI_WORK_DIST[MPI_MY_RANK][0]-1; i < MPI_WORK_DIST[MPI_MY_RANK][1]; i++)
		{
			int index = quick_index_totalDevice[i];
			if (modules[index].isDevice())
			{	
				// Hole filling
				modules[index].getDevice()->sbe_initialize(DT); // Simulation variables
			} else
			{
				// Hole filling
				modules[index].getTwoArmDevice()->sbe_initialize(DT); // Simulation variables
			}
		}	
	#endif	
	
	if (MPI_MY_RANK==0)
	{
		//=================================
		// Cavity snapshot - Also in LOAD
		//=================================
		// Setup output to files
		#ifdef USE_CAVITY_SNAPSHOT
			//Not setup for TwoArmCavity
			std::stringstream baseName;
			baseName << getToFileOutputKey();
	
			fileName.str("");
			fileName << baseName.str() << "cav_snapshot_x.dat";
			openAppendBinary(&output_Cav_Snapshot_x, fileName.str());
			
			int inda = quick_index_cavity[0];
			int indb = quick_index_cavity[quick_index_cavity.size()-1];
			
		
			// Find correct indices
			//VECSEL_cav_snapshot_index = new int[VECSEL_cav_snapshot_num_points];
			double xi = modules[inda].getPosition0();
			int count = 0;
			while(xi < modules[indb].getPosition1())
			{
				VECSEL_cav_snapshot_x.push_back(xi);
				
				// Find correct cavity
				for(int k = 0; k < modules.size(); k++)
				{
					if (((VECSEL_cav_snapshot_x[count] >= modules[k].getPosition0())&&(VECSEL_cav_snapshot_x[count] < modules[k].getPosition1()))&&(modules[k].isCavity()))
					{
						VECSEL_cav_snapshot_index.push_back(k);
						break;
					} 
				}
				int indx = VECSEL_cav_snapshot_index[count];
				
				double DX;
				if (modules[indx].getCavity()->getRefInd() > 1.0)
				{
					DX = 0.2*0.25*getLambda()/modules[indx].getCavity()->getRefInd(); // 5 points per L/4
					
				} else {
					// In AIR we will sample at one rate
					DX = 0.5*0.25*getLambda()/modules[indx].getCavity()->getRefInd(); // dx <=0.25 L/n
				}
				
				
				if (xi+DX < modules[indx].getPosition1())
				{
					xi += DX;
				} else {
					xi = modules[indx].getPosition1();
				}
				count += 1;
			}
			
			VECSEL_cav_snapshot_num_points = count;
			cout << "CAVITY SNAPSHOT: num points = " << VECSEL_cav_snapshot_num_points << endl;
			
			output_Cav_Snapshot_x.write(reinterpret_cast<const char*>(&VECSEL_cav_snapshot_x[0]),VECSEL_cav_snapshot_num_points*sizeof(double));
			output_Cav_Snapshot_x.close();
		
			// Initialize E
			VECSEL_cav_snapshot_E 	= new std::complex<double>[VECSEL_cav_snapshot_num_points];
			VECSEL_cav_snapshot_E_re = new double[VECSEL_cav_snapshot_num_points];
			VECSEL_cav_snapshot_E_im = new double[VECSEL_cav_snapshot_num_points];
		
			// Initialize counters
			VECSEL_cav_snapshot_output_wait  = ceil(cav_snapshot_freq/DT);
			VECSEL_cav_snapshot_output_count = VECSEL_cav_snapshot_output_wait;	
		#endif

	}

	// Find LENS cavity and move to correct quick index list
	// Above the round-trip time includes the lens cavity
	// thus: this should only change the Maxwell update iteration
	quick_index_cavity_noBPM = quick_index_cavity;
	int BPM_ctr=0;
	for(unsigned i=0; i<quick_index_cavity.size(); i++)
	{
		int indx   = quick_index_cavity[i];
		if (modules[indx].getCavity()->getName().substr(0,3) == "BPM")
		{
			if(modules[indx].getCavity()->getName().substr(0,5) == "BPMHC")
			{
				// Add to quick_index_cavity_halfCav
				quick_index_cavity_lens_halfCav.push_back(indx);
			} else if (modules[indx].getCavity()->getName().substr(0,5) == "BPMFS")
			{
			  	quick_index_cavity_freeSpace.push_back(indx);
			} else if (modules[indx].getCavity() -> getName().substr(0,8) == "BPMTRANS")
			{
				// Add to quick_index_cavity_lens
				quick_index_cavity_lens.push_back(indx);
			} else
			{
				cout<<"BPM cavity not recognized"<<endl;
				exit(-1);
			}
			// Remove from quick_index_cavity_noBPM
			quick_index_cavity_noBPM.erase(quick_index_cavity_noBPM.begin()+i-BPM_ctr);
			BPM_ctr++;
		}
	}

	if (quick_index_birefringentCrystal.size()>0 && quick_index_birefringentCrystal.back()==getNumberModules())
	{
		cout<<"Birefringent crystal at last element. Untested. Quitting."<<endl;
		exit(-1);
	}

	if (quick_index_twoArmCavity.size()>0)
	{
	// Find two arm cavities with NO QWs between them
	for(unsigned i=0; i<quick_index_twoArmCavity.size()-1; i++)
	{
		int i1   = quick_index_twoArmCavity[i];
		int i2   = quick_index_twoArmCavity[i+1];
		
		if (i2-i1 == 1)
		{
			// No distance => No QWs
			quick_index_twoArmCavity_noQW.push_back(i1);
		} else {	
			// Can be QWs between
			bool has_qw = false;
			for(int j = i1; j <= i2; j++)
			{
				if (modules[j].isTwoArmDevice())
				{
					has_qw = true;
					break;
				}
			}
			if (has_qw == true)
			{
				quick_index_twoArmCavity_QW.push_back(i1);
				quick_index_twoArmCavity_QW.push_back(i2);
	
			} else {
				quick_index_twoArmCavity_noQW.push_back(i1);
			}
		}
	}
	quick_index_twoArmCavity_noQW.push_back(quick_index_twoArmCavity.back());
	}
	if(quick_index_cavity.size()>0)
	{
	// Find cavities with NO QWs between them
	for(unsigned i=0; i<quick_index_cavity.size()-1; i++)
	{
		int i1   = quick_index_cavity[i];
		int i2   = quick_index_cavity[i+1];
	
		
		if (i2-i1 == 1)
		{
			// No distance => No QWs
			quick_index_cavity_noQW.push_back(i1);
		} else {
			// Can be QWs between
			bool has_qw = false;
			for(int j = i1; j <= i2; j++)
			{
				if (modules[j].isDevice())
				{
					has_qw = true;
					break;
				}
			}
			if (has_qw == true)
			{
				quick_index_cavity_QW.push_back(i1);
				quick_index_cavity_QW.push_back(i2);
	
			} else {
				quick_index_cavity_noQW.push_back(i1);
			}
		}
	}
	quick_index_cavity_noQW.push_back(quick_index_cavity.back());
	}
}


/* The initaial pulse in system 
 * Will create a seed pulse for a given amount of time.
 * There are multiple options for the type of pulse
 *  -> Bump function
 *  -> Sinc function
 *  -> (default) Sech()^2
 * Or one can use
 *  -> dual frequency input
 *  -> Single frequency input
 *  -> Triangular pulse input 
 * Or one can turn it off for spont. emission. as input.
 * */
void VECSEL::maxwell_initial_E(double x_t , std::complex<double> *tmp)
{
//	if (init_VECSEL_iteration < VECSEL_ROUND_TRIP_ITERATIONS)
	if (abs(x_t-VECSEL_initial_delay) < 1000*ps)
	{
		init_VECSEL_iteration++;

		// Energy: int(|E|^2) = 4*(amp^2)/(3*const1)
		double const1 = 1.76275/VECSEL_initial_fwhm;// Correct constant , fwhm of |E(t)|
		//double const1 = 1.21169/VECSEL_initial_fwhm;// Correct constant , fwhm of |E(t)|^2

		double mc = cosh((x_t - VECSEL_initial_delay)*const1);
		std::complex<double> shift = exp(-I*(VECSEL_initial_energy_shift)*(x_t - VECSEL_initial_delay));
		std::complex<double> z_profile = VECSEL_initial_amplitude*shift/(mc*mc);

		// Scale transverse profile with the z_profile
		memset(tmp,0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
		cblas_zaxpy(VECSEL_transverse_points_number, &z_profile, VECSEL_initial_transverse_pulse_profile, 1, tmp, 1);


	} else {
		memset(tmp,0,VECSEL_transverse_points_number*sizeof(std::complex<double>));
	}
}



//=================
// Misc functions

/* Return the cavity field evaluatated at the given device nr. Not setup for VCAV */
void VECSEL::misc_getFieldAtDevice(int device_nr, std::complex<double> *tmp)
{
	if(quick_index_device.size() > 0)
	{
		if (device_nr >= quick_index_totalDevice.size())
		{
			cout << "misc_getFieldAtDevice(): Requesting device not in list" << endl;
			cout << "Requesting   = " << device_nr << endl;
			cout << "# of devices = " << quick_index_device.size() << endl;
		}
		
		// Index of previous cavity
		int indx_k 					= quick_index_device_previous_cavity[device_nr];
		modules[indx_k].getCavity()->evaluateEprop_x1(tmp);
		
	} else {
		// Use boundary
		
		if (modules.back().isBoundary())
		{
			// Index of previous cavity
			int indx_k 					= quick_index_cavity.back();

			modules[indx_k].getCavity()->evaluateEprop_x1(tmp);
			
		} else {
			cout << "misc_getFieldAtDevice():: Cannot detect any cavities.." << endl;
			exit(-1);
		}
	}
}


//=====================
// Filter functions

void VECSEL::setFilter_pluss_minus_pass()
{
	for(unsigned i = 0; i < modules.size(); i++)
	{
		if (modules[i].isFilter()){
			modules[i].getFilter()->setFilter_pluss_PassThroughDelay(); // Delay "filter"
			modules[i].getFilter()->setFilter_minus_PassThroughDelay(); // Delay "filter"
		}
	}
}

void VECSEL::setFilter_pluss_gauss_minus_pass(double wa_s, double wb_s, double width_s, double w0_s)
{
	for(unsigned i = 0; i < modules.size(); i++)
	{
		if (modules[i].isFilter()){
			modules[i].getFilter()->setFilter_minus_PassThroughDelay(); // Delay "filter"
			modules[i].getFilter()->setFilter_pluss_doubleGauss_flatTopWindow(wa_s, wb_s, width_s, w0_s); // Gauss filter
		}
	}
}


//=====================================================
// Functions for access to devices and their variables

/* Disable the spontaneous emissions from all devices
 * If one wants to do this all the time, it is faster to turn it off from the flag.
 * */
void VECSEL::device_disable_spont_emissions_all()
{
	#ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		if (modules[indx].isDevice())
		{
			modules[indx].getDevice()->sbe_set_spont_emissions(false);
		} else
		{
			modules[indx].getTwoArmDevice()->sbe_set_spont_emissions(false);
		}
	}
	#endif
}

/* When using a realistic pump and you want to change the parameters
 * Set the pump parameters for each device
 * */
void VECSEL::device_set_qw_real_pump_parameters(double W0, double E0, double ETA, double nCavity)
{
	// TODO: Modify Pump strength based on QW placement and pump absorbtion...
	int ind_qw_cav = 0;
	int ind_air_cav = 0;
	
	// Possible expansion: Set strength different for each QW based on pump frequency
	#ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		if (modules[indx].isDevice())
		{
			modules[indx].getDevice()->sbe_set_real_pump_model(W0, E0, ETA,VECSEL_DT);
		} else
		{
			modules[indx].getTwoArmDevice()->sbe_set_real_pump_model(W0, E0, ETA,VECSEL_DT);
		}
	}
	#endif	
}

/* Change background carrier density in all devices
 * */
void VECSEL::device_set_background_carrier_density(double new_density,int NUM_QW)
{
	#ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		if (modules[indx].isDevice())
		{
			modules[indx].getDevice()->sbe_set_background_carrier_density(new_density);
		} else
		{
			modules[indx].getTwoArmDevice()->sbe_set_background_carrier_density(new_density);
		}
	}
	#endif
}

/* Set a delay for when computations are supposed to start
 * Used on in special situations
 * */
void VECSEL::device_set_computational_delay(double delay)
{
	// For a given delay in computation
	// computations will start after the delay
	#ifdef ITERATE_QW
	#pragma omp parallel for num_threads(OMP_THREADS_LEVEL_1)
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		if(modules[indx].isDevice())
		{
			modules[indx].getDevice()->sbe_delay_computations(delay);
		} else
		{	
			modules[indx].getTwoArmDevice()->sbe_delay_computations(delay);
		}
	}
	#endif

	cout << "device_set_computational_delay: " << delay/ps << " [ps] Complete" << endl;
}



//=================================
// Diagnostics function

/* Remove all field content in the cavity
 * */
void VECSEL::diagnostics_zero_all_cavity_fields()
{
	for(unsigned i = 0; i < modules.size(); i++)
	{
		if (modules[i].isCavity())
		{
			modules[i].getCavity()->clearFields();
		} else if (modules[i].isFilter()){
			modules[i].getFilter()->clearFields();
		}
	}
}

/* Set this in order to REMOVE the polarization feedback into the pulse
 * */
void VECSEL::diagnostics_set_qw_feedback(double newFeedback)
{
	VECSEL_QW_FEEDBACK = newFeedback;
}

/* Return the maximally allowed timestep based on cavity lengths
 * */
double VECSEL::diagnostics_findMaxTimestep()
{
	double dt = 1000;
	
	if(quick_index_cavity.size()>0)
	{
	for(unsigned i=0; i<quick_index_cavity.size(); i++)
	{
		int indx = quick_index_cavity[i];

		//Find next cavity's refractive index
		double leftRefInd = modules[indx].getCavity()->getRefInd();
		double leftWidth  = modules[indx].getCavity()->getWidth();
	
		// Get minimal width
		if (leftWidth/(c0/leftRefInd) < dt)
		{
			dt = leftWidth/(c0/leftRefInd);
		}
		
	}
	}

	if(quick_index_twoArmCavity.size()>0)
	{
	for(unsigned i=0; i<quick_index_twoArmCavity.size(); i++)
	{
		int indx = quick_index_twoArmCavity[i];

		//Find next cavity's refractive index
		double leftRefInd = modules[indx].getTwoArmCavity()->getRefInd();
		double leftWidth  = modules[indx].getTwoArmCavity()->getWidth();
	
		// Get minimal width
		if (leftWidth/(c0/leftRefInd) < dt)
		{
			dt = leftWidth/(c0/leftRefInd);
		}
		
	}
	}

	return dt;

}

/* To remove everything from this device
 * */
void VECSEL::diagnostics_clear_VECSEL()
{
	cout << "VECSEL: Remove all modules" << endl;
	for(unsigned i=0; i<getNumberModules(); i++)
    {
    	modules[i].Remove();
    }
    modules.clear();
    
    cout << "VECSEL: Clearing variables" << endl;
    setName("");
    setLambda(0.0);
    
}

/* Initialize MPI worker arrays
 * distribute the total work as equal as possible over all nodes
 * Note: Will quit if the number of nodes > number of devices
 * */
void VECSEL::mpi_initialize_work_arrays(int total_devices)
{
	// Set up QW distribution depending on #QWs and #NODES
	if (MPI_WORK_GANG_SIZE < 1)
	{
		cout << "mpi_initialize_work_arrays()::ERROR cannot have number of nodes < 1.." << endl;
		cout << "total devices = " << total_devices << endl;
		cout << "# of nodes  = " << MPI_WORK_GANG_SIZE << endl;
		cout << " Please run program with number of nodes >= 1" << endl << endl;
		MPI_Finalize();
		exit(-1);	
	}

	// Set up work distribution
	MPI_WORK_DIST = new int*[MPI_WORK_GANG_SIZE];
	for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
	{
		MPI_WORK_DIST[i] = new int[2]; // Store QW # from-to
	}

	if (total_devices < 0)
	{
		cout << "mpi_initialize_work_arrays()::ERROR cannot process < 0 devices.." << endl;
		cout << "total devices = " << total_devices << endl;
		cout << "# of nodes  = " << MPI_WORK_GANG_SIZE << endl;
		cout << " Please run program with total number of devices >= 0" << endl << endl;
		MPI_Finalize();
		exit(-1);	
	}

	if ((MPI_WORK_GANG_SIZE>1)&&(total_devices < MPI_WORK_GANG_SIZE))
	{
		cout << "mpi_initialize_work_arrays()::ERROR Too MANY nodes used!" << endl;
		cout << "total devices = " << total_devices << endl;
		cout << "# of nodes  = " << MPI_WORK_GANG_SIZE << endl;
		cout << " Please run program with number of nodes <= number of devices" << endl << endl;
		MPI_Finalize();
		exit(-1);	
	}

	MPI_WORK_DIST_TOTAL = total_devices;

	MPI_Barrier(*MPI_WORK_GANG);
	//===============================
	// Work scheduling
	MPI_WORK_DEVICE_SIZE = new int[MPI_WORK_GANG_SIZE];
	if (total_devices > 0)
	{
		MPI_LoadBalancer_index_set = new int[total_devices];
		for(int i = 0; i < total_devices; i++)
		{
			MPI_LoadBalancer_index_set[i] = i;
		}
		#ifdef MPI_BALANCE_WORKLOAD

			// Use uniform distribution to allow for easy timing
			//====================================================================
			double work_per_thread = ((double)MPI_WORK_DIST_TOTAL/(double)MPI_WORK_GANG_SIZE)/((double)OMP_THREADS_LEVEL_1);
			if (work_per_thread >= 10)
			{
				// Safe region
			} else {
				cout << "Hybrid MPI-OpenMP issue.." << endl;
				cout << "The # work/thread is too low and speedup will saturate above work/thread < 10" << endl;
				cout << "Parallel speedup ratio is saturating, but total time will still improve" << endl;
				cout << "Total work   = " << MPI_WORK_DIST_TOTAL << endl;
				cout << "MPI_nodes    = " << MPI_WORK_GANG_SIZE << endl;
				cout << "Thread/node  = " << OMP_THREADS_LEVEL_1 << endl;
				cout << "work/thread  = " << work_per_thread << endl;
			//	exit(-1);
			}

			// MPI_WORK_DIST[n][0] with elements in [1,MPI_WORK_DIST_TOTAL]
			// MPI_WORK_DIST[n][1] with elements in [1,MPI_WORK_DIST_TOTAL]
			// such that MPI_WORK_DIST[n][0] <= MPI_WORK_DIST[n][1]
			// for each MPI_node: n
			int base_val = floor((double)MPI_WORK_DIST_TOTAL/(double)MPI_WORK_GANG_SIZE);
			int remain = MPI_WORK_DIST_TOTAL - base_val*MPI_WORK_GANG_SIZE;
			int remove = MPI_WORK_GANG_SIZE-1-remain;

			if (remain > 0)
			{
				// Work is NOT multiple of workforce
				MPI_WORK_DEVICE_SIZE[0] = base_val - remove;
				for(int i = 1; i < MPI_WORK_GANG_SIZE; i++)
				{
					MPI_WORK_DEVICE_SIZE[i] = base_val + 1;
				}
			} else {
				// Work is multiple of workforce
				for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
				{
					MPI_WORK_DEVICE_SIZE[i] = base_val;
				}
			}


			//========================================================
			// Finetune hybrid MPI node usage for better performance.
			// Only tuned for FAST QWs i.e. no 2nd Born scattering
			//#if !defined(USE_ISAK_HOLE_FILLING) && !defined(USE_ISAK_HOLE_FILLING_TABLE)
			#if defined(__ICC) || defined(__INTEL_COMPILER)
			if (1.0/(double)MPI_WORK_GANG_SIZE < 0.1)
			#else
			if (false)
			#endif
			{
				// Reduce MASTER workload
				// This reduces memory load times for the master node AND reduces wait times
				// but only usefull for large number of MPI nodes.
				if (MPI_WORK_DEVICE_SIZE[0] >= MPI_WORK_GANG_SIZE-1)
				{
			
					// Option 1: NO QWs for master
					int cnt = MPI_WORK_DEVICE_SIZE[0];
					MPI_WORK_DEVICE_SIZE[0] = 0;
					int counter = 1;
					while (cnt>0)
					{
						if( counter >= MPI_WORK_GANG_SIZE)
						{
							counter = 1;
						}
		
						MPI_WORK_DEVICE_SIZE[counter] += 1;
						cnt--;
						counter++;
					}

			
	/*				
					// Option 2: All workers have SAME number of QWs, master takes remaining or 0
					int mult = floor((double)MPI_WORK_DEVICE_SIZE[0]/((double)MPI_WORK_GANG_SIZE-1));
					int scale = +mult; // Give master less work
					MPI_WORK_DEVICE_SIZE[0] = MPI_WORK_DEVICE_SIZE[0] - scale*(MPI_WORK_GANG_SIZE-1);
					for(int i = 1; i < MPI_WORK_GANG_SIZE; i++)
					{
						MPI_WORK_DEVICE_SIZE[i] = MPI_WORK_DEVICE_SIZE[i] + scale;
					}

	*/
				}
			} else if ((MPI_WORK_GANG_SIZE > 1)&&((work_per_thread < 200))) {

				// For cases where there is a medium amount of work, but not enough MPI workers

				// Reduce work on MASTER by 10%, if possible
				// This will compensate for memory load times for very high number of QWs..
				// This number (10%) can be changed if you think this is a problem..
				if (MPI_WORK_DEVICE_SIZE[0] >= MPI_WORK_GANG_SIZE-1)
				{
					int mult = floor(0.1*MPI_WORK_DEVICE_SIZE[0]/((double)MPI_WORK_GANG_SIZE-1));
					int scale = +mult; // Give master less work
					MPI_WORK_DEVICE_SIZE[0] = MPI_WORK_DEVICE_SIZE[0] - scale*(MPI_WORK_GANG_SIZE-1);
					for(int i = 1; i < MPI_WORK_GANG_SIZE; i++)
					{
						MPI_WORK_DEVICE_SIZE[i] = MPI_WORK_DEVICE_SIZE[i] + scale;
					}
				}
			}
		#else

			// Use preset from file  OR attempt to estimate best setting
			if (fileExists("MPI_CONFIG_WEIGHTS.dat"))
			{
				// IF MPI preset exists, load from file
				// Format:
				// Line 1: Number of subgroups, Sub-group-size
				// Line 2: Total number of devices, Index set of jobs

				FILE *mpi_setup = fopen("MPI_CONFIG_WEIGHTS.dat","r+");
				if (mpi_setup != NULL)
				{
					double maxwell_time, maxwell_std;
					double device_times[total_devices];

					// Load MPI work sizes
					int row_size = -1;
					fscanf(mpi_setup, "%d", &row_size);
					if (row_size == total_devices)
					{
						
						fscanf(mpi_setup, "%lf %lf", &maxwell_time, &maxwell_std);
						for(int i = 0; i < row_size; i++)
						{
							double time, std;
							fscanf(mpi_setup, "%lf %lf", &time, &std);

							device_times[i] = time + 3.0*std;
						}
						fclose(mpi_setup);
						
					} else {
						cout << "mpi_initialize_work_arrays()::ERROR file MPI_CONFIG_WEIGHTS.dat is set for a different total number of devices" << endl;
						cout << "total number of devices  = " << total_devices << endl;
						cout << "value in file = " << row_size << endl;
						exit(-1);
					}

					MPI_LoadBalancer = new parSchedule();
					if (1.0/(double)MPI_WORK_GANG_SIZE < 0.05)
					{
						MPI_LoadBalancer->optimize_schedule(total_devices, device_times, MPI_WORK_GANG_SIZE-1, 0.0, 1.0);
						MPI_WORK_DEVICE_SIZE[0] = 0;
						MPI_LoadBalancer->get_sub_group_info(&(MPI_WORK_DEVICE_SIZE[1]));

					} else {
						MPI_LoadBalancer->optimize_schedule(total_devices, device_times, MPI_WORK_GANG_SIZE, 2.0*(maxwell_time+3.0*maxwell_std), 1.25);
						MPI_LoadBalancer->get_sub_group_info(&(MPI_WORK_DEVICE_SIZE[0]));
					}
					MPI_LoadBalancer->get_sub_group_index(MPI_LoadBalancer_index_set);

					// Apply balancing
					int current_device_index[total_devices];
					int current_device_index_prev_cav[total_devices];
					for(int i = 0; i < total_devices; i++)
					{
						current_device_index[i] = quick_index_totalDevice[i];
						current_device_index_prev_cav[i] = quick_index_device_previous_cavity[i];
					}
					for(int i = 0; i < total_devices; i++)
					{
						quick_index_totalDevice[i] 		= current_device_index[MPI_LoadBalancer_index_set[i]];
						quick_index_device_previous_cavity[i] 	= current_device_index_prev_cav[MPI_LoadBalancer_index_set[i]];
					}

					if (MPI_MY_RANK == 0)
					{
						MPI_LoadBalancer->file_write_structure("MPI_CONFIG_BALANCED.dat");
					}
					delete MPI_LoadBalancer;
				}
				
			} else {
				// No MPI preset found, using density data to esitmate weights and attempt balancing
				// Weights for boole: 1.86-2.02
				// Weights for hamilton: 2.49-2.72
				double est_weights[total_devices];
				double max_weight = 1.0;
				for(int i = 0; i < total_devices; i++)
				{
					int index = quick_index_totalDevice[i];
					if (modules[index].isDevice())
					{
						double qw_density_scale = modules[index].getDevice()->getTransverseBackgroundDensityScale();
						if (qw_density_scale > 0.75)
						{
							est_weights[i] = 1.0;
						} else 
						{
						est_weights[i] = 1.94; // boole
						//est_weights[i] = 2.6; // hamilton, but gives only a minor correction
						}
					} else
					{
						double qw_density_scale = modules[index].getTwoArmDevice()->getTransverseBackgroundDensityScale();
						if (qw_density_scale > 0.75)
						{
							est_weights[i] = 1.0;
						} else 
						{
						est_weights[i] = 1.94; // boole
						//est_weights[i] = 2.6; // hamilton, but gives only a minor correction
						}
						est_weights[i]*=3.0;
					}	

					if (est_weights[i] > max_weight)
					{
						max_weight = est_weights[i];
					}
				}

				MPI_LoadBalancer = new parSchedule();
				if (1.0/(double)MPI_WORK_GANG_SIZE < 0.05)
				{
					MPI_LoadBalancer->optimize_schedule(total_devices, est_weights, MPI_WORK_GANG_SIZE-1, 0.0, 1.0);
					MPI_WORK_DEVICE_SIZE[0] = 0;
					MPI_LoadBalancer->get_sub_group_info(&(MPI_WORK_DEVICE_SIZE[1]));
				} else {
					MPI_LoadBalancer->optimize_schedule(total_devices, est_weights, MPI_WORK_GANG_SIZE, 2*max_weight, 1.25);
					MPI_LoadBalancer->get_sub_group_info(&(MPI_WORK_DEVICE_SIZE[0]));
				}
				MPI_LoadBalancer->get_sub_group_index(MPI_LoadBalancer_index_set);

				if (MPI_MY_RANK == 0)
				{
					MPI_LoadBalancer->file_write_structure("MPI_CONFIG_BALANCED_est.dat");
				}
				delete MPI_LoadBalancer;
				
				// Apply balancing
				int current_device_index[total_devices];
				int current_device_index_prev_cav[total_devices];
				for(int i = 0; i < total_devices; i++)
				{
					current_device_index[i] = quick_index_totalDevice[i];
					current_device_index_prev_cav[i] = quick_index_device_previous_cavity[i];
				}
				for(int i = 0; i < total_devices; i++)
				{
					quick_index_totalDevice[i] 			= current_device_index[MPI_LoadBalancer_index_set[i]];
					//WTFF-quick_index_device_previous_cavity[i] 	= current_device_index_prev_cav[MPI_LoadBalancer_index_set[i]];
				}

				
			}

		
		#endif
	} else {
		for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
		{
			MPI_WORK_DEVICE_SIZE[i] = 0;
		}

	}

	MPI_Barrier(*MPI_WORK_GANG);

/*

*/
	//====================================================================
	MPI_WORK_DEVICE_TIMING_SIZE = new int[MPI_WORK_GANG_SIZE];
	for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
	{
		MPI_WORK_DEVICE_TIMING_SIZE[i] = 2*MPI_WORK_DEVICE_SIZE[i];
	}

	MPI_WORK_DEVICE_TIMING_OFFSET = new int[MPI_WORK_GANG_SIZE];
	MPI_WORK_DEVICE_TIMING_OFFSET[0] = 0;
	for(int i = 1; i < MPI_WORK_GANG_SIZE; i++)
	{
		MPI_WORK_DEVICE_TIMING_OFFSET[i] = MPI_WORK_DEVICE_TIMING_OFFSET[i-1] + MPI_WORK_DEVICE_TIMING_SIZE[i-1];
	}
	
	// Transfer size data to workforce
	if (MPI_WORK_DEVICE_SIZE[0]>0)
	{
		MPI_WORK_DIST[0][0] = 1;
		MPI_WORK_DIST[0][1] = MPI_WORK_DEVICE_SIZE[0];
	} else {
		MPI_WORK_DIST[0][0] = 1;
		MPI_WORK_DIST[0][1] = 0;
	}
	for(int i = 1; i < MPI_WORK_GANG_SIZE; i++)
	{
		MPI_WORK_DIST[i][0] = MPI_WORK_DIST[i-1][1] + 1;
		MPI_WORK_DIST[i][1] = MPI_WORK_DIST[i-1][1] + MPI_WORK_DEVICE_SIZE[i];
	}

	
/*
	// Simple work schedule
	MPI_WORK_DIST[0][0] = 1;
	MPI_WORK_DIST[0][1] = floor((double)MPI_WORK_DIST_TOTAL/(double)MPI_WORK_GANG_SIZE); // Work of master node, floor
	//MPI_WORK_DIST[0][1] = ceil((double)MPI_WORK_DIST_TOTAL/(double)MPI_WORK_GANG_SIZE); // Work of master node, ceil
	for(int i = 1; i < MPI_WORK_GANG_SIZE; i++)
	{
		int tmp1 = floor((double)(MPI_WORK_DIST_TOTAL-MPI_WORK_DIST[i-1][1])/(double)(MPI_WORK_GANG_SIZE-i));
		MPI_WORK_DIST[i][0] = MPI_WORK_DIST[i-1][1] + 1;
		MPI_WORK_DIST[i][1] = MPI_WORK_DIST[i-1][1] + tmp1;
	}
*/
	

	// MPI scatterv and gatherv organization arrays
	
	MPI_WORK_DIST_E_OFFSET = new int[MPI_WORK_GANG_SIZE];
	MPI_WORK_DIST_E_OFFSET[0] = 0;
	MPI_WORK_DIST_E_SIZE   = new int[MPI_WORK_GANG_SIZE];
	MPI_WORK_DIST_E_SIZE[0] = 8*(MPI_WORK_DIST[0][1] - MPI_WORK_DIST[0][0]+1);
	for(int i = 1; i < MPI_WORK_GANG_SIZE; i++)
	{
		MPI_WORK_DIST_E_OFFSET[i] = MPI_WORK_DIST_E_OFFSET[i-1] + MPI_WORK_DIST_E_SIZE[i-1];
		MPI_WORK_DIST_E_SIZE[i]   = 8*(MPI_WORK_DIST[i][1] - MPI_WORK_DIST[i][0]+1);
	}
	MPI_WORK_DIST_P_OFFSET = new int[MPI_WORK_GANG_SIZE];
	MPI_WORK_DIST_P_OFFSET[0] = 0;
	MPI_WORK_DIST_P_SIZE   = new int[MPI_WORK_GANG_SIZE];
	MPI_WORK_DIST_P_SIZE[0] = 4*(MPI_WORK_DIST[0][1] - MPI_WORK_DIST[0][0]+1);
	for(int i = 1; i < MPI_WORK_GANG_SIZE; i++)
	{
		MPI_WORK_DIST_P_OFFSET[i] = MPI_WORK_DIST_P_OFFSET[i-1] + MPI_WORK_DIST_P_SIZE[i-1];
		MPI_WORK_DIST_P_SIZE[i]   = 4*(MPI_WORK_DIST[i][1] - MPI_WORK_DIST[i][0]+1);
	}

	// Find max number of WORK for any worker
	int MPI_WORK_DIST_E_SIZE_MAX = 0;
	int MPI_WORK_DIST_P_SIZE_MAX = 0;
	for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
	{
		if (MPI_WORK_DIST_E_SIZE[i] > MPI_WORK_DIST_E_SIZE_MAX)
		{
			MPI_WORK_DIST_E_SIZE_MAX = MPI_WORK_DIST_E_SIZE[i];
		}
		if (MPI_WORK_DIST_P_SIZE[i] > MPI_WORK_DIST_P_SIZE_MAX)
		{
			MPI_WORK_DIST_P_SIZE_MAX = MPI_WORK_DIST_P_SIZE[i];
		}
	}
	
	// Local memory of work
	MPI_WORK_DIST_E_LOCAL = new std::complex<double>[MPI_WORK_DIST_E_SIZE_MAX];
	MPI_WORK_DIST_P_LOCAL = new std::complex<double>[MPI_WORK_DIST_P_SIZE_MAX];

	for(int i =0; i < MPI_WORK_DIST_E_SIZE_MAX; i++)
	{
		MPI_WORK_DIST_E_LOCAL[i] = 0;
	}
	for(int i =0; i < MPI_WORK_DIST_P_SIZE_MAX; i++)
	{
		MPI_WORK_DIST_P_LOCAL[i] = 0;
	}
	MPI_WORK_DIST_P_GLOBAL = new std::complex<double>*[3];
	for(int i =0; i < 3; i++)
	{
		MPI_WORK_DIST_P_GLOBAL[i] = NULL;
	}

	if (MPI_MY_RANK == 0)
	{
		// In the master thread, the polarizations for each QW are all in one array
		MPI_WORK_DIST_E_GLOBAL = new std::complex<double>[8*MPI_WORK_DIST_TOTAL];
		for(int j = 0; j < 8*MPI_WORK_DIST_TOTAL; j++)
		{
			MPI_WORK_DIST_E_GLOBAL[j] = 0;
		}
		
		for(int i = 0; i < 3; i++)
		{
			MPI_WORK_DIST_P_GLOBAL[i] = new std::complex<double>[4*MPI_WORK_DIST_TOTAL];
			for(int j = 0; j < 4*MPI_WORK_DIST_TOTAL; j++)
			{
				MPI_WORK_DIST_P_GLOBAL[i][j] = 0;
			}
		}
		MPI_LoadBalancer_P_tmp = new std::complex<double>[4*MPI_WORK_DIST_TOTAL];
		MPI_LoadBalancer_E_tmp = new std::complex<double>[8*MPI_WORK_DIST_TOTAL];
		
		cout << "MPI WORK DISTRIBUTED OVER " << MPI_WORK_GANG_SIZE << " NODES" << endl;
		cout << "MPI_WORK_DIST[][]: " << endl;
		for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
		{
			cout << "MPI_WORK_DIST[" << i << "] = [" << MPI_WORK_DIST[i][0] << ", " << MPI_WORK_DIST[i][1] << "]" << endl;
		}
		cout << "MPI_WORK_DIST_E_OFFSET[] = ";
		for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
		{
			cout << MPI_WORK_DIST_E_OFFSET[i] << " ";
		}
		cout << endl;
		cout << "MPI_WORK_DIST_E_SIZE[]: ";
		for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
		{
			cout << MPI_WORK_DIST_E_SIZE[i] << " ";
		}
		cout << endl;
		cout << "MPI_WORK_DIST_P_OFFSET[] = ";
		for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
		{
			cout << MPI_WORK_DIST_P_OFFSET[i] << " ";
		}
		cout << endl;
		cout << "MPI_WORK_DIST_P_SIZE[]: ";
		for(int i = 0; i < MPI_WORK_GANG_SIZE; i++)
		{
			cout << MPI_WORK_DIST_P_SIZE[i] << " ";
		}
		cout << endl;
	}

	// SET UP MPI TIMERS
	#ifdef USE_MAIN_TIMERS
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		std::stringstream oldName;
		int indx = quick_index_totalDevice[j]; // Index of device
		if(modules[indx].isDevice())
		{
		oldName << modules[indx].getDevice()->getName();
		} else
		{
		oldName << modules[indx].getTwoArmDevice()->getName();
		}		

		#ifdef USE_DEVICE_TIMERS
		MainStat->newSubTimer("VECSEL::compute all QWs",oldName.str());
		#endif
	}
	#endif
	
}

/* Initialize MPI rank
 */
void VECSEL::mpi_initialize_rank(int my_rank)
{
	MPI_MY_RANK = my_rank;	
}

/* Set the communication group for MPI workers
 */
void VECSEL::mpi_set_work_gang(MPI_Comm *work_gang)
{
	MPI_WORK_GANG = work_gang;	
	MPI_Comm_size(*MPI_WORK_GANG,&MPI_WORK_GANG_SIZE);
}


