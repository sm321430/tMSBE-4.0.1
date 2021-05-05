#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h> // For sleep(1)
#include <stdlib.h> // For system()

#include <mpi.h>

#include "constantsAndMiscUnits.h"
#include "classVECSEL.h"
#include "fileIO.cpp"
#include "classCommunication.h"


//! A class to keep track of file-output simulation meta-data
/*! Data such as number of output files, when to write to output, when to write to save, etc. is stored in counter. This class collects these counters and basic printing of information */
class out_slice_param
{
	public:
		//! A constructor
		out_slice_param()
		{
			C_wait 		= 0;
			C_write		= 0;
			C_freq		= 0;
			COUNT_WAIT  = 0;
			COUNT_WRITE = 0;
			COUNT_FREQ 	= 0;
			out_counter = 0;
			NEW_OUTPUT_FILE = 1;

			save_counter = 0;
			C_save     = 0;
			COUNT_SAVE = 0;

		}

		//! A constructor
		out_slice_param(double DT, double OUT_DT_WAIT, double OUT_DT_WRITE, double OUT_DT_FREQ, double SAVE_DT_WAIT)
		{
			//C_wait 		= 0;
			C_write		= 0;
			C_freq		= 0;
			COUNT_WAIT  	= ceil(OUT_DT_WAIT/DT);
			COUNT_WRITE 	= ceil(OUT_DT_WRITE/DT);
			COUNT_FREQ  	= ceil(OUT_DT_FREQ/DT);
			C_wait 		= COUNT_WAIT;
			out_counter 	= 0;
			NEW_OUTPUT_FILE = 1;
			
			COUNT_SAVE 	= ceil(SAVE_DT_WAIT/DT);
			C_save 		= COUNT_SAVE;
			save_counter 	= 0;
		}

		//! Return true if writing can begin
		bool isDoneWaiting()
		{
			return (C_wait >= COUNT_WAIT);
		}

		//! Return true if currently writing to file
		bool isWriting()
		{
			return (C_write < COUNT_WRITE);
		}

		//! Return true if writing to file should start
		bool writeNow()
		{
			return (C_freq % COUNT_FREQ == 0);
		}

		//! Return true if currently saving to file
		bool isSaving()
		{
			return (C_save >= COUNT_SAVE);
		}

		//! Return true if a new output file should be generated
		int get_new_output_file()
		{
			return NEW_OUTPUT_FILE;
		}
		//! Set output file name
		void set_new_output_file(int newFile)
		{
			NEW_OUTPUT_FILE = newFile;
		}
		//! Return counter for number of output-files
		int get_out_counter()
		{
			return out_counter;
		}
		//! Increment counter for number of output-files
		void update_out_counter()
		{
			out_counter++;
		}

		//! Return counter for number of save-files
		int get_save_counter()
		{
			return save_counter;
		}
		//! Increment counter for number of save files
		void update_save_counter()
		{
			save_counter++;
		}
		
		//! Increment wait time counter
		void update_c_wait()
		{
			C_wait++;
		}
		//! Increment write counter
		void update_c_write()
		{
			C_write++;
		}
		//! Increment write counter
		void update_c_freq()
		{
			C_freq++;
		}
		//! Increment save counter
		void update_c_save()
		{
			C_save++;
		}
		//! Set wait counter
		void set_c_wait(int newCount)
		{
			C_wait = newCount;
		}
		//! Set write counter
		void set_c_write(int newCount)
		{
			C_write = newCount;
		}
		//! Set write counter
		void set_c_freq(int newCount)
		{
			C_freq = newCount;
		}
		//! Set save counter
		void set_c_save(int newCount)
		{
			C_save = newCount;
		}

		//! Reset write counter, wait counter, output file counter, and save_counter. Create new output-file and save-file
		void reset()
		{
			//C_wait 		= 0;
			C_write		= 0;
			C_freq		= 0;
			//COUNT_WAIT  = ceil(OUT_DT_WAIT/DT);
			//COUNT_WRITE = ceil(OUT_DT_WRITE/DT);
			//COUNT_FREQ  = ceil(OUT_DT_FREQ/DT);
			C_wait 		= COUNT_WAIT;
			out_counter 	= 0;
			NEW_OUTPUT_FILE = 1;
			
			//COUNT_SAVE = 10000000;
			C_save = COUNT_SAVE;
			save_counter = 0;
	
		}

		//! Print output data to screen
		void Print()
		{
			cout << "C_wait      = " << C_wait << endl;
			cout << "C_write     = " << C_write << endl;
			cout << "C_freq      = " << C_freq << endl;
			cout << "COUNT_WAIT  = " << COUNT_WAIT << endl;
			cout << "COUNT_WRITE = " << COUNT_WRITE << endl;
			cout << "COUNT_FREQ  = " << COUNT_FREQ << endl;
			cout << "out_counter = " << out_counter << endl;
			cout << "NEW_OUTPUT_FILE = " << NEW_OUTPUT_FILE << endl;
			cout << "save_counter = " << save_counter << endl;
			cout << "C_save      = " << C_save << endl;
			cout << "COUNT_SAVE  = " << COUNT_SAVE << endl;
		}
	private:
		//! Counter for waiting/ no output time
		int C_wait;

		//! Counter for writing /output time
		int C_write;

		//! Counter for writing (might be every n'th timestep)
		int C_freq;
		
		//! Number of # timesteps to wait
		int COUNT_WAIT;
		//! Number of # timesteps to write
		int COUNT_WRITE;
		//! Number of # timesteps to write
		int COUNT_FREQ;
		//! Output file counter
		int out_counter;
		//! Flag to write a new output file or append to old file
		int NEW_OUTPUT_FILE;

		//! Save-file counter
		int save_counter;

		//! Counter for when to save
		int C_save;

		//! Counter for how long between each save
		int COUNT_SAVE;
	
};


void out_slice(double t, VECSEL *, int);
void out_slice(double t, VECSEL *, out_slice_param *);
void runSimulation(VECSEL *, double, double,double,double,double,double,int);
void readInValuesFromFile(double *, double *, double *, double *, double *);

void saveMiscVariables(int, double  );
void loadMiscVariables(int, double *);
void save_slice(double, VECSEL *,int);
void save_slice(double, VECSEL *, out_slice_param *);

using namespace std;


#ifdef USE_MAIN_TIMERS
	#include "myTimer.cpp"
	myTimerCentral *MainStat = new myTimerCentral();
#endif


// Output counters
int 	C_wait 		= 0;
int 	C_write		= 0;
int 	C_freq		= 0;
int 	COUNT_WAIT  = 0;
int 	COUNT_WRITE = 0;
int		COUNT_FREQ 	= 0;
int 	out_counter = 0;
int NEW_OUTPUT_FILE = 1;

int     save_counter = 0;
int     C_save     = 0;
int		COUNT_SAVE = 0;

int LOAD_FROM_FILE_NR = -1;


int main(int argc, char *argv[])
{

#ifdef USE_MAIN_TIMERS
	MainStat->newTimer("runSimulation");
	MainStat->newSubTimer("runSimulation","main file i/o");
	MainStat->newSubTimer("runSimulation","VECSEL::iterateModules");
	MainStat->newSubTimer("VECSEL::iterateModules","VECSEL::transfer matrix");
	MainStat->newSubTimer("VECSEL::transfer matrix","TM I.C.");
	MainStat->newSubTimer("VECSEL::transfer matrix","TM Update Storage");
	MainStat->newSubTimer("VECSEL::transfer matrix","TM B.C.");
	MainStat->newSubTimer("VECSEL::iterateModules","VECSEL::prepare MPI");
	MainStat->newSubTimer("VECSEL::iterateModules","VECSEL::MPI_comm_scatt");
	MainStat->newSubTimer("VECSEL::iterateModules","VECSEL::MPI_comm_gath");
	MainStat->newSubTimer("VECSEL::iterateModules","VECSEL::compute all QWs");
#endif
	// MPI info
	int   MPI_MY_RANK = 0; // Default, will also take care of case where MPI is not used
	char  MPI_MY_HOSTNAME[256];
	
	// Start MPI
	int mpi_required_threading = MPI_THREAD_FUNNELED;
	int mpi_provided_threading;
	MPI_Init_thread(&argc, &argv, mpi_required_threading, &mpi_provided_threading);
	//MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&MPI_MY_RANK); // Get the rank of this current process
	gethostname(MPI_MY_HOSTNAME,255);


	if (mpi_provided_threading < mpi_required_threading)
	{
		if (MPI_MY_RANK==0)
		{
			cout << "MPI: Warning: This MPI implementation provides insufficient threading support!" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}
	cout << "MPI: Hello world! I am process number: " << MPI_MY_RANK << " on host: " << MPI_MY_HOSTNAME << endl;
	
	#ifdef USE_OPENMP
	omp_set_nested(1); // Enable nested parallelism
	
	if (MPI_MY_RANK==0)
	{
		cout << endl;
		// Running program with OpenMP
		#pragma omp parallel
		{
			#pragma omp master
			{               
				//cout << "OMP USING T  = " << NUM_THREADS << " OpenMP threads" << endl;
				// Get environment information 
				int procs = omp_get_num_procs();
				int nthreads = omp_get_num_threads();
				int maxt = omp_get_max_threads();
				int inpar = omp_in_parallel();
				int dynamic = omp_get_dynamic();
				int nested = omp_get_nested();

				// Print environment information 
				printf("OpenMP info\n");
				printf("Number of processors = %d\n", procs);
				printf("Number of threads = %d\n", nthreads);
				printf("Max threads = %d\n", maxt);
				printf("In parallel? = %d\n", inpar);
				printf("Dynamic threads enabled? = %d\n", dynamic);
				printf("Nested parallelism supported? = %d\n", nested);  


				if (OMP_THREADS_LEVEL_1*OMP_THREADS_LEVEL_2 > nthreads)
				{
					cout << "OpenMP: Requesting that we use: " << endl;
					cout << "OMP_THREADS_LEVEL_1 = " << OMP_THREADS_LEVEL_1 << endl;
					cout << "OMP_THREADS_LEVEL_2 = " << OMP_THREADS_LEVEL_2 << endl;
					cout << "but the max number of threads is = " << nthreads << endl;
					exit(-1);
				}
			}
		}
		
		#pragma omp parallel
		{
			#pragma omp critical
			{
				cout << "Thread[" << omp_get_thread_num() << "]: Reporting for duty.." << endl;
			}
		}
		cout << endl;
	}
	#endif
	
	if (out_DT_write >= t_max)
	{
		cout << "Write time has to be less than the simulation time" << endl;
		exit(-1);
	}

	
	LOAD_FROM_FILE_NR = LOAD_STATE_NR; // Load from file

	/*
		Construction manual: (!!)
		
		1. A VECSEL is the primary class of the program. Create one of these and add more modules inside
		2. A MODULE is anything inside the VECSEL, QW, ABS, BOUNDARY,... whatever, add as many as you like
		3. Configure each module with the correct parameters: density ... whatever is in the module
		4. There HAS to be TWO boundaries, one on each end.
	*/
	
	Cavity *cav;
	Device *dev;
	Boundary *bnd;
	Module *mod;
	
	
	//double lambda = 2.0*Pi*hbar*c0/(1.298221484*e);
	//double lambda = 2.0*Pi*hbar*c0/(1.230*e);
	double lambda = 2.0*Pi*c0/w0;
	cout << "lambda = " << lambda/nm << " [nm]" << endl;
	cout << "w0     = " << hbar*w0/e << " [eV]" << endl;
	cout << "SVEA limit = " << (lambda/c0)/fs << " [fs]" << endl;
	
	// BASE ABSORBTION
	//double OUTPUT_COUPLING = 0.99*0.99;
	//double ANGLE_OF_INCIDENCE = Pi*0.0/180.0;
	double OUTPUT_COUPLING = 0.0;
	double ANGLE_OF_INCIDENCE = 0.0;
	double pulse_initial_fwhm = 0.0;		 // Full width half max
	double pulse_initial_amplitude = 0.0;		 // amplitude of initial pulse
	double pulse_initial_energy_shift = 0.00*e/hbar; // Shift in w

	readInValuesFromFile(&OUTPUT_COUPLING, &ANGLE_OF_INCIDENCE,&pulse_initial_fwhm, &pulse_initial_amplitude,&pulse_initial_energy_shift);

	// BASE GAIN
	VECSEL *baseGain = new VECSEL;
	baseGain->setName("dummyGain");
	baseGain->setToFileOutputKey("dummyGain__");
	baseGain->setLambda(lambda);

	baseGain->set_transverse_dimensions(NUM_TRANSVERSE_POINTS,TRANSVERSE_DOMAIN_SIZE, TRANSVERSE_DOMAIN_SIZE_BG_RATIO); // Transverse domain [-R/2, (1 - 2/N)R/2] with super gaussian boundary guard FWHM

	baseGain->set_transverse_QW_pump_profile_SuperGaussian(13,GAIN_SPOT_FWHM); // order, FWHM. // Corresponds to a waist w0 = 175um
//	baseGain->set_transverse_QW_temp_profile_SuperGaussian(2 ,287.3*um); // order, FWHM. NOTE: sets QW temperature profile relative to a lattice (300K) (material__QW#.config temperature at center of SG)
	baseGain->set_transverse_QW_temp_profile_file(11249, "transverse_temperature_profile_11249.txt",PUMP_SCALER); // NOTE: sets QW temperature profile relative to a lattice (300K) (material__QW#.config temperature at peak)

//	baseGain->set_transverse_QW_pump_profile_SuperGaussian(13,8000.0*um); // order, FWHM. // Corresponds to a waist w0 = 175um
//baseGain->set_transverse_QW_temp_profile_SuperGaussian(13,8000.0*um); // order, FWHM. NOTE: sets QW temperature profile relative to a lattice (300K) (material__QW#.config temperature at center of SG)
	
	// ===== MQW 3x 121 + 1
/*
	std::complex<double> n_AU = 0.2426102 + 6.7059991*I; // FIT TO GROWN STRUCTURE
	baseGain->addBoundary(0.0,1.0);
	baseGain->addCavity(200.0*nm, n_AU); //Au
	baseGain->addCavity(62.0*nm, 3.1978, ANGLE_OF_INCIDENCE, 1.0);
	baseGain->addDBR_LEFT_nL(3.3998,76.93*nm,2.9551,88.8*nm,2*16, ANGLE_OF_INCIDENCE, 1.0);

	double qwWidth[13] = { 8.215*nm,
			      98.84*nm,
			      32.51*nm,
			      16.43*nm,
			      99.95*nm,
			      32.51*nm,
			      16.43*nm,
			      32.51*nm,
			      68.79*nm,
			      32.51*nm,
			      16.43*nm,
			      32.51*nm,
			      32.635*nm};


	baseGain->addCUSTOM_QW_LEFT_OPTIMAL_ANGLE(12, qwWidth, 0.0, 3.4653, 3.1978, n_vcsel_ar_coating,GAIN_AR_LENGTH,false, 1.0, 1.0, ANGLE_OF_INCIDENCE,1.0, GAIN_TEMPERATURE_FACTOR*BAD_GAIN_TEMP);
	// ...
*/
	// ==== MQW Nx121+1
/*
	baseGain->addBoundary(0.0,1.0);
	baseGain->addDBR_LEFT(n_dbr_1,n_dbr_2,num_dbr_layers, n_vcsel_bgr, ANGLE_OF_INCIDENCE, n_cav);
	int numClusters[6] = {4,4,4,4,4,4};
	baseGain->addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(6, numClusters, true, 0, 0, n_vcsel_bgr, n_vcsel_cap_layer, n_vcsel_ar_coating, ANGLE_OF_INCIDENCE, n_cav);
	// ...
*/

	// ===== RPG
//	baseGain->addBoundary(0.0,1.0);
//	baseGain->addDBR_LEFT(n_dbr_1,n_dbr_2,num_dbr_layers, n_vcsel_bgr, ANGLE_OF_INCIDENCE, n_cav);
//	baseGain->addRPG_QW_LEFT(num_qw, length_qw_dbr, length_qw_qw,length_qw_total,n_vcsel_bgr,n_vcsel_cap_layer,n_vcsel_ar_coating,ANGLE_OF_INCIDENCE, n_cav);
//	baseGain->addRPG_QW_LEFT_EFF(6, num_qw, length_qw_dbr, length_qw_qw,length_qw_total,n_vcsel_bgr,n_vcsel_cap_layer,n_vcsel_ar_coating,ANGLE_OF_INCIDENCE, n_cav);
	// ....

	// ===== MQW10
//	baseGain->addBoundary(0.0,1.0);
//	baseGain->addDBR_LEFT(n_dbr_1,n_dbr_2,num_dbr_layers, n_vcsel_bgr, ANGLE_OF_INCIDENCE, n_cav);
//	int numClusters[1] = {10};
//	baseGain->addCUSTOM_CLUSTER_QW_LEFT(1, numClusters, 0, 0, n_vcsel_bgr, n_vcsel_cap_layer, n_vcsel_ar_coating, ANGLE_OF_INCIDENCE, n_cav);
	// ...

	// ===== VCAV
	// ----- OC -> Length1 -> RPG
/*	baseGain->addBoundary(1.0,1.0);
	baseGain->addRPG_SPACE_LENS(length_cav, n_cav,VCAV_focus, VCAV_length1,  0.0);
	baseGain->addTwoArmRPG_QW_FRONT_EFF(5, num_qw, length_qw_dbr, length_qw_qw,length_qw_total,n_vcsel_bgr,n_vcsel_cap_layer,n_vcsel_ar_coating,ANGLE_OF_INCIDENCE, n_cav);
	baseGain->addTwoArmDBR_FRONT(n_dbr_1,n_dbr_2,num_dbr_layers, n_vcsel_bgr, ANGLE_OF_INCIDENCE, n_cav);
*/


	int NUM_QW = 8; // For optimization alg..
	



#ifdef COMPILE_REFLECTION_CALCULATION
	if (true)
#else
	if(false)
#endif
	{
		//======================================================
		//  Equilibrium gain / Absorbtion
		//======================================================
		// REFLECTION SPECTRUM
		cout << "Calculating reflection spectrum" << endl;
		VECSEL *refSpecQW = new VECSEL;
		refSpecQW->setName("vcselSpectrum");
		refSpecQW->setToFileOutputKey("refSpecQW__");
		refSpecQW->setLambda(lambda);
		refSpecQW->setInitialPulseAmplitude(2.0e4);	// Weak probe pulse
		refSpecQW->setInitialPulseFwhm(100.0*fs);		// Short probe
		refSpecQW->setInitialPulseDelay(0.5*ps);		// Short enough for 10*fs pulse

		double probe_pulse_transverse_FWHM = LENS_SPOT_GAIN; // 

		//Rescaled cavity length to include half of pulse delay in original length
		refSpecQW->set_transverse_dimensions(NUM_TRANSVERSE_POINTS,TRANSVERSE_DOMAIN_SIZE, TRANSVERSE_DOMAIN_SIZE_BG_RATIO); // Transverse domain [-R/2, (1 - 2/N)R/2] with super gaussian boundary guard FWHM
		// RPG STRUCTURE		
		if ((STARTL == 0)&&(STARTR==1))
		{
			double length_cav_refSpec=500.25;
			refSpecQW->setStartingPoint(0,1); 
			//refSpecQW->appendStructure(baseGain);
			refSpecQW->addBoundary(0.0,1.0);
			//refSpecQW->addDBR_LEFT(n_dbr_1,n_dbr_2,num_dbr_layers, n_vcsel_bgr, ANGLE_OF_INCIDENCE, n_cav);
			//refSpecQW->addRPG_QW_LEFT_EFF(6, num_qw, length_qw_dbr, length_qw_qw,length_qw_total,n_vcsel_bgr,n_vcsel_cap_layer,n_vcsel_ar_coating,ANGLE_OF_INCIDENCE, n_cav);	
			//refSpecQW->addRPG_SPACE_ANGLE_LEFT(10.0, n_cav,-0.125,"CAVL",ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);
			
			//refSpecQW->addRPG_SPACE_LENS_ANGLE(length_cav_refSpec, n_cav,VCAV_focus, VCAV_length1,  0.0, ANGLE_OF_INCIDENCE);
			refSpecQW->addRPG_SPACE_ANGLE(500.0, n_cav,-0.125,"CAVL", ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);
			refSpecQW->addRPG_SPACE_ANGLE(1.0, n_cav,-0.125,"CAVM", ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);
			refSpecQW->addTwoArmRPG_QW_FRONT_EFF(6, 6, num_qw, length_qw_dbr, length_qw_qw,length_qw_total,n_vcsel_bgr,n_vcsel_cap_layer,n_vcsel_ar_coating, n_vcsel_ar_coating2, ANGLE_OF_INCIDENCE, n_cav);
			refSpecQW->addTwoArmDBR_FRONT(n_dbr_1,n_dbr_2,num_dbr_layers, n_vcsel_bgr, ANGLE_OF_INCIDENCE, n_cav);
			//refSpecQW->addTwoArmQW_STRUCT(6, 800.0, 1.25, 1.0, ANGLE_OF_INCIDENCE, 1.0);
			refSpecQW->addRPG_SPACE_BPM_ANGLE(length_cav_refSpec, VCAV_length2,  n_cav,0.375,"CAVC", ANGLE_OF_INCIDENCE);
			//refSpecQW->addRPG_SPACE_ANGLE(500.0, n_cav,-0.125,"CAVR", ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);
			refSpecQW->addBoundary(0.0,1.0); // Reflecting on the right

			//refSpecQW->addRPG_QW_simple(0, 1.0, 500.0, 1.0, 0.0, 1.0);
			//refSpecQW->addRPG_SPACE_LENS_transform(length_cav, n_cav, optical_abcd_A, optical_abcd_B, optical_abcd_C, optical_abcd_D);
			//refSpecQW->addTwoArmRPG_QW_FRONT_EFF(5, 6, num_qw, length_qw_dbr, length_qw_qw,length_qw_total,n_vcsel_bgr,n_vcsel_cap_layer,n_vcsel_ar_coating, n_vcsel_ar_coating2, ANGLE_OF_INCIDENCE, n_cav);
			//refSpecQW->addTwoArmDBR_FRONT(n_dbr_1,n_dbr_2,num_dbr_layers, n_vcsel_bgr, ANGLE_OF_INCIDENCE, n_cav);
			//refSpecQW->addRPG_SPACE_BPM(length_cav, VCAV_length2,  n_cav,0.375,"CAVR");
			//refSpecQW->addRPG_SPACE_ANGLE(10.0, n_cav,-0.125,"CAVP", ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);
			//refSpecQW->addRPG_ABS_RIGHT(num_abs, length_abs_dbr,length_abs_abs,length_abs_total,n_abs_bgr,n_abs_cap_layer,n_abs_ar_coating);
		
			refSpecQW->set_transverse_QW_pump_profile_SuperGaussian(13,GAIN_SPOT_FWHM); // order, FWHM. // Corresponds to a waist w0 = 175um
			refSpecQW->set_transverse_QW_temp_profile_file(11249, "transverse_temperature_profile_11249.txt",PUMP_SCALER); // NOTE: sets QW temperature profile relative to a lattice (300K) (material__QW#.config temperature at peak)
			refSpecQW->setInitialPulseTransverseFWHM(LENS_SPOT_GAIN); // Initial pulse transverse FWHM
		
		} else if ((STARTL == 1)&&(STARTR==0))
		{
			cout << " Needs fixing" << endl;
			exit(-1);
		} else {
			cout << "Reflection spectrum, cannot have a pulse from both sides" << endl;
			exit(-1);
		}
		//refSpecQW->setInitialPulseDualPulse(-0.0780*e/hbar,0.0800*e/hbar); // FOR DBM
		double DOMAIN_MAX_DT = refSpecQW->diagnostics_findMaxTimestep();
		double SIM_DT = 0.1*fs;
		if (DOMAIN_MAX_DT < SIM_DT)
		{
			SIM_DT = DOMAIN_MAX_DT*(0.99);
			cout << "Timestep reset to dt = " << SIM_DT/fs << " [fs]" << endl;
		}

		int max_work_qw = (refSpecQW->getNumberDevices() + refSpecQW->getNumberTwoArmDevices()); // Get number of QWs
		if (max_work_qw==0)
		{
			max_work_qw = 1; // Always let 0 work be done by rank=0 worker
		}
		if (MPI_MY_RANK < max_work_qw) // These are the workers that go into a problem
		{
			// Create new work gang
			MPI_Comm work_gang;
			MPI_Comm_split(MPI_COMM_WORLD, 0, MPI_MY_RANK, &work_gang);

			// Do work with new gang
			refSpecQW->mpi_set_work_gang(&work_gang);
			refSpecQW->maxwell_initialize(SIM_DT,MPI_MY_RANK); // Set all simulation variables to zero
			//refSpecQW->device_set_qw_real_pump_parameters(2.0*Pi*c0/(808*nm), 1.0e7, 1.0e12,n_cav); // Must be after init.
			//refSpecQW->file_load_device_carriers_line(0, 1);
			//refSpecQW->file_load_device_carriers_ascii(37832);

			refSpecQW->device_disable_spont_emissions_all();
	//		refSpecQW->file_output_device_disable();
	//		refSpecQW->file_output_device_reduce(TRANSVERSE_DOMAIN_OUTPUT_FWHM);
			refSpecQW->file_output_device_set_freq(TRANSVERSE_DOMAIN_OUTPUT_FREQ, TRANSVERSE_DOMAIN_OUTPUT_FWHM, 2, 0);

			MPI_Barrier(work_gang); // To ensure syncronization of timers

			if (MPI_MY_RANK==0)
			{
				refSpecQW->Print();
				refSpecQW->file_output_structure(); // Print structure to file
			}

			MPI_Barrier(work_gang); // To ensure syncronization of timers
			runSimulation(refSpecQW, 16*ps, SIM_DT, 15*ps, SIM_DT, 10*ps, 100*ps,MPI_MY_RANK);
			#ifdef USE_MAIN_TIMERS
			for(int i = 0; i < max_work_qw; i++) {
				MPI_Barrier(work_gang);
				if (i == MPI_MY_RANK) {
					MainStat->Print_short();
				}
			}
			#endif

			// Free work gang
			MPI_Comm_free(&work_gang);
		} else {
			MPI_Comm slacking_gang;
			MPI_Comm_split(MPI_COMM_WORLD, 1, MPI_MY_RANK, &slacking_gang);
			MPI_Comm_free(&slacking_gang);
		}
		MPI_Barrier(MPI_COMM_WORLD); // Wait for all workers to finish jobs
		//MPI_Finalize();
		//exit(-1); 
	
		// ABSORBTION SPECTRUM
		cout << "Calculating absorbtion spectrum" << endl;
		VECSEL *refSpecABS = new VECSEL;
		refSpecABS->setName("AbsorberSpectrum");
		refSpecABS->setToFileOutputKey("refSpecABS__");
		refSpecABS->setLambda(lambda);
		refSpecABS->setInitialPulseAmplitude(2.0e4);	// Weak probe pulse
		refSpecABS->setInitialPulseFwhm(100*fs);		// Short probe
		refSpecABS->setInitialPulseDelay(0.5*ps);		// Short enough for 10*fs pulse

		refSpecABS->set_transverse_dimensions(NUM_TRANSVERSE_POINTS,TRANSVERSE_DOMAIN_SIZE, TRANSVERSE_DOMAIN_SIZE_BG_RATIO); // Transverse domain
		
		if ((STARTL == 0)&&(STARTR==1))
		{
			refSpecABS->setStartingPoint(1, 0);
			refSpecABS->addBoundary(0.0,1.0);
			//refSpecABS->addRPG_SPACE_LENS_transform(600.25, n_cav, optical_abcd_A, optical_abcd_B, optical_abcd_C, optical_abcd_D);
			//refSpecABS->addRPG_SPACE_ANGLE(1.0, n_cav,-0.125,"CAVOC",0, 0);
			//refSpecABS->addRPG_SPACE_BPM(750.0, VCAV_length2,  n_cav,0.375,"CAVL");	
			refSpecABS->addRPG_SPACE_ANGLE(1.0, n_cav,-0.125,"CAVM",0, 0);
			//refSpecABS->addRPG_SPACE_BPM(500.0, VCAV_length2,  n_cav,0.375,"CAVL");	
			//refSpecABS->addTwoArmQW_STRUCT(6, 1.0, 1.25, 1.0, ANGLE_OF_INCIDENCE, 1.0);
			//refSpecABS->addTwoArmQW_STRUCT(-1, 1.0, 1.25, 1.0, ANGLE_OF_INCIDENCE, 1.0, -1.0);
			//refSpecABS->addQW_STRUCT(-1, 1.0, 1.25, 1.0, ANGLE_OF_INCIDENCE, 1.0);
			//refSpecABS->addTwoArmRPG_ABS(num_abs, length_abs_dbr,length_abs_abs,length_abs_total,n_abs_bgr,n_abs_cap_layer,n_abs_ar_coating, n_cav, ANGLE_OF_INCIDENCE);
			refSpecABS->addRPG_SPACE_BPM(1000.0, VCAV_length2,  n_cav,0.375,"CAVR");
			refSpecABS->addRPG_ABS_RIGHT(num_abs, length_abs_dbr,length_abs_abs,length_abs_total,n_abs_bgr,n_abs_cap_layer,n_abs_ar_coating);
			//refSpecABS->addBoundary(OUTPUT_COUPLING,1.0);
			refSpecABS->addBoundary(1.0,1.0);

/*
			// SESAM WITH DBR ON RIGHT
			refSpecABS->addRPG_SPACE_ANGLE(600.25, n_cav,-0.125,"CAVP",ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);
			refSpecABS->addLossElement(1.0-OUTPUT_COUPLING,0.0); // Reduce left moving wave, give beter GAIN/ABS picture, NO EXTRA GDD
			refSpecABS->addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(1.99037079,130*nm,  // Ar index/length
								     9.175*nm, 19.175*nm,3.521966, // SESAM length dx0/dx1/index
								     2.969682, 87.17243*nm,3.526966, 73.5332039*nm, 2*30, // DBR n1/L1/n2/L2/#
								     3.52196644, 100.0*nm, // Phase layer n/length
								     ANGLE_OF_INCIDENCE, 1.0); // theta/external index
			refSpecABS->addBoundary(0.0,1.0);
*/
		refSpecABS->setInitialPulseTransverseFWHM(LENS_SPOT_SESAM); // Initial pulse transverse FWHM
			
		} else if ((STARTL == 1)&&(STARTR==0))
		{
			cout << "Need to check" << endl;
			exit(-1);
			
		} else {
			cout << "Reflection spectrum, cannot have a pulse from both sides" << endl;
			exit(-1);
		}

		DOMAIN_MAX_DT = refSpecABS->diagnostics_findMaxTimestep();
		SIM_DT = 0.1*fs;
		if (DOMAIN_MAX_DT < SIM_DT)
		{
			SIM_DT = DOMAIN_MAX_DT*(0.99);
			cout << "Timestep reset to dt = " << SIM_DT/fs << " [fs]" << endl;
		}

		int max_work_abs = refSpecABS->getNumberDevices(); // Get number of QWs
		if (max_work_abs==0)
		{
			max_work_abs = 1; // Always let 0 work be done by rank=0 worker
		}
		if (MPI_MY_RANK < max_work_abs) // These are the workers that go into a problem
		{
			// Create new work gang
			MPI_Comm work_gang;
			MPI_Comm_split(MPI_COMM_WORLD, 0, MPI_MY_RANK, &work_gang);

			// Do work with new gang
			refSpecABS->mpi_set_work_gang(&work_gang);
		
			refSpecABS->maxwell_initialize(SIM_DT,MPI_MY_RANK); // Set all simulation variables to zero
		//	refSpecABS->file_output_device_disable();
		//	refSpecABS->file_output_device_reduce(TRANSVERSE_DOMAIN_OUTPUT_FWHM);
			refSpecABS->file_output_device_set_freq(TRANSVERSE_DOMAIN_OUTPUT_FREQ, TRANSVERSE_DOMAIN_OUTPUT_FWHM, 0, 0);
		//	refSpecABS->set_transverse_ABS_fake_focus(3.1623,LENS_SPOT_SESAM);

			if (MPI_MY_RANK==0)
			{
				refSpecABS->Print();
				refSpecABS->file_output_structure(); // Print structure to file
			}

			MPI_Barrier(work_gang); // To ensure syncronization of timers
			
			runSimulation(refSpecABS, 21*ps, SIM_DT, 20*ps, SIM_DT,30.0*ps,30.0*ps,MPI_MY_RANK);

			#ifdef USE_MAIN_TIMERS
			for(int i = 0; i < max_work_qw; i++) {
				MPI_Barrier(work_gang);
				if (i == MPI_MY_RANK) {
					MainStat->Print_short();
				}
			}
			#endif

			// Free work gang
			MPI_Comm_free(&work_gang);
		} else {
			MPI_Comm slacking_gang;
			MPI_Comm_split(MPI_COMM_WORLD, 1, MPI_MY_RANK, &slacking_gang);
			MPI_Comm_free(&slacking_gang);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		exit(-1);

	
	}
	
	/*========================
		Run main simulation
	========================*/
	{
		// MAIN DEVICE
		VECSEL *primary = new VECSEL;
		primary->setName("primary test");
		primary->setToFileOutputKey("out__");
		primary->setLambda(lambda);
		primary->setInitialPulseAmplitude(pulse_initial_amplitude);
		primary->setInitialPulseFwhm(pulse_initial_fwhm);
		primary->setInitialPulseDelay(pulse_initial_delay);
		primary->setInitialPulseEnergyShift(pulse_initial_energy_shift);

		primary->set_transverse_dimensions(NUM_TRANSVERSE_POINTS,TRANSVERSE_DOMAIN_SIZE, TRANSVERSE_DOMAIN_SIZE_BG_RATIO); // Transverse domain [-R/2, (1 - 2/N)R/2] with super gaussian boundary guard FWHM
		//primary->setStartingPoint(STARTL, STARTR);
		primary->setStartingPoint(1, 0);
		if ((STARTL == 0)&&(STARTR==1))
		{
			primary->addBoundary(OUTPUT_COUPLING,1.0);
			//primary->addBoundary(0.0,1.0);

			//LINEAR CAVITY STRUCTURE
			//primary->addDBR_LEFT(n_dbr_1,n_dbr_2,num_dbr_layers, n_vcsel_bgr, ANGLE_OF_INCIDENCE, n_cav);
			//primary->addRPG_QW_LEFT_EFF(6, num_qw, length_qw_dbr, length_qw_qw,length_qw_total,n_vcsel_bgr,n_vcsel_cap_layer,n_vcsel_ar_coating,ANGLE_OF_INCIDENCE, n_cav);
			//primary->addRPG_SPACE_LENS_ANGLE(1200.25, n_cav,VCAV_focus, VCAV_length1,  0.0, ANGLE_OF_INCIDENCE);
			//primary->addRPG_SPACE_BPM_ANGLE(length_cav, VCAV_length2,  n_cav,0.375,"CAVC", ANGLE_OF_INCIDENCE);
			//primary->addRPG_SPACE_ANGLE(1.0, n_cav,-0.125,"CAVP", ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);


			//THIS IS THE PRIMARY STRUCTURE
			//primary->addRPG_SPACE_ANGLE(1.0, n_cav,-0.125,"CAVOC", ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);
			//primary->addRPG_SPACE_LENS_ANGLE(1600.25, n_cav,VCAV_focus, VCAV_length1,  0.0, ANGLE_OF_INCIDENCE);
			//primary->addRPG_SPACE_BPM_ANGLE(1600.25, VCAV_length2,  n_cav,0.375,"CAVC", ANGLE_OF_INCIDENCE);
			primary->addTwoArmRPG_ABS_LEFT(num_abs, length_abs_dbr,length_abs_abs,length_abs_total,n_abs_bgr,n_abs_cap_layer,n_abs_ar_coating, n_cav, ANGLE_OF_INCIDENCE);
			primary->addTwoArm_Space(1500.0,1.0,ANGLE_OF_INCIDENCE,1.0, "CAV1");
			//primary->addBirefringentCrystal(3190.48543,1.6431,1.4799, ANGLE_OF_INCIDENCE); //Calcite crystal
			primary->addBirefringentCrystal(200.0,1.0,1.165, ANGLE_OF_INCIDENCE); //Calcite crystal
			primary->addTwoArm_Space(1500.0,1.0,ANGLE_OF_INCIDENCE,1.0, "CAV2");
			primary->addTwoArmRPG_QW_FRONT_EFF(6, 6, num_qw, length_qw_dbr, length_qw_qw,length_qw_total,n_vcsel_bgr,n_vcsel_cap_layer,n_vcsel_ar_coating, n_vcsel_ar_coating2, ANGLE_OF_INCIDENCE, n_cav);
			primary->addTwoArmDBR_FRONT(n_dbr_1,n_dbr_2,num_dbr_layers, n_vcsel_bgr, ANGLE_OF_INCIDENCE, n_cav);
			//primary->addTwoArmQW_noInterface(6, 50.0, 300.125, 1.0, ANGLE_OF_INCIDENCE, 1.0);
			//primary->addRPG_SPACE_LENS_ANGLE(1600.25, n_cav,VCAV_focus, VCAV_length1,  0.0, ANGLE_OF_INCIDENCE);
			//primary->addQW_STRUCT(6, 1.0, 1.25, 1.0, ANGLE_OF_INCIDENCE, 1.0);
			//primary->addRPG_SPACE_ANGLE(1.0, n_cav,-0.125,"CAVGAIN2", ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);
			//primary->addRPG_SPACE_ANGLE(1.0, n_cav,-0.125,"CAVGAIN1", ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);
			//primary->addRPG_SPACE_BPM_ANGLE(1600.0, VCAV_length2,  n_cav,0.375,"CAVC", ANGLE_OF_INCIDENCE);
			//primary->addRPG_SPACE_ANGLE(1.0, n_cav,-0.125,"CAVP", ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);
			//primary->addRPG_SPACE_LENS_ANGLE(1700.0, n_cav,VCAV_focus, VCAV_length1,  0.0, ANGLE_OF_INCIDENCE);
			//primary->addTwoArmQW_STRUCT(6, 5.0, 1.25, 1.0, ANGLE_OF_INCIDENCE, 1.0, 1.0);
			//primary->addTwoArmRPG_QW_FRONT_EFF(6, 6, num_qw, length_qw_dbr, length_qw_qw,length_qw_total,n_vcsel_bgr,n_vcsel_cap_layer,n_vcsel_ar_coating, n_vcsel_ar_coating2, ANGLE_OF_INCIDENCE, n_cav);
			//primary->addTwoArmDBR_FRONT(n_dbr_1,n_dbr_2,num_dbr_layers, n_vcsel_bgr, ANGLE_OF_INCIDENCE, n_cav);
			//primary->addRPG_SPACE_BPM_ANGLE(1500.0, VCAV_length2,  n_cav,0.375,"CAVM", ANGLE_OF_INCIDENCE);
			//primary->addTwoArmQW_STRUCT(6, 5.0, 1.25, 1.0, ANGLE_OF_INCIDENCE, 1.0, -1.0);
			//primary->addTwoArmRPG_QW_FRONT_EFF(6, 6, num_qw, length_qw_dbr, length_qw_qw,length_qw_total,n_vcsel_bgr,n_vcsel_cap_layer,n_vcsel_ar_coating, n_vcsel_ar_coating2, ANGLE_OF_INCIDENCE, n_cav);
			//primary->addTwoArmDBR_FRONT(n_dbr_1,n_dbr_2,num_dbr_layers, n_vcsel_bgr, ANGLE_OF_INCIDENCE, n_cav);
			
			//primary->addRPG_SPACE_LENS_ANGLE(3200.0, n_cav, 0.0, VCAV_length2+VCAV_length1,  VCAV_focus, ANGLE_OF_INCIDENCE);
			//primary->addTwoArmRPG_ABS(num_abs, length_abs_dbr,length_abs_abs,length_abs_total,n_abs_bgr,n_abs_cap_layer,n_abs_ar_coating, n_cav, ANGLE_OF_INCIDENCE);
			//primary->addRPG_SPACE_BPM_ANGLE(3200.0, VCAV_length2,  n_cav,0.375,"CAVR", ANGLE_OF_INCIDENCE);
			//primary->addRPG_ABS_RIGHT(num_abs, length_abs_dbr,length_abs_abs,length_abs_total,n_abs_bgr,n_abs_cap_layer,n_abs_ar_coating);
			
			//SIMPLE AIR V-CAV

			//EXTRA STRUCTURES FOR REFERENCE
			//primary->addBoundary(0.0,1.0);
			//primary->addRPG_SPACE_ANGLE(1.0, n_cav,-0.125,"CAVOC", ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);
			//primary->addTwoArmQW_STRUCT(6, 800.0, 1.25, 1.0, ANGLE_OF_INCIDENCE, 1.0);
			//primary->addRPG_SPACE_ANGLE(1.0, n_cav,-0.125,"CAVP", ANGLE_OF_INCIDENCE, ANGLE_OF_INCIDENCE);
			//primary->addBoundary(0.0,1.0);
		 

			//REFERENCE FUNCTIONS	
			//primary->addTwoArmQW_STRUCT(6, 5.0, 1.25, 1.0, ANGLE_OF_INCIDENCE, 1.0, 0.0);
			//primary->appendStructure(baseGain);	
			//primary->addQW_STRUCT(6, 800.0, 1.25, 1.0, ANGLE_OF_INCIDENCE, 1.0);
			//primary->addRPG_SPACE(length_cav, n_cav,0.375,"CAVP");
			//primary->addRPG_SPACE_LENS_transform(length_cav, n_cav, optical_abcd_A, optical_abcd_B, optical_abcd_C, optical_abcd_D);
			//primary->addBoundary(OUTPUT_COUPLING,1.0);
			//primary->addRPG_SPACE_ANGLE(length_cav, n_cav,0.375,"CAVP",ANGLE_OF_INCIDENCE, 0.0);
			//primary->addRPG_SPACE_ANGLE_LEFT(10.0, n_cav,-0.125,"CAVL",0.0, 0.0);
			//primary->addRPG_SPACE_LENS(length_cav/2.0, n_cav,VCAV_focus, VCAV_length1,  0.0);
			//primary->addTwoArmDBR_FRONT(n_dbr_1,n_dbr_2,num_dbr_layers, n_vcsel_bgr, ANGLE_OF_INCIDENCE, n_cav);	
			//primary->addRPG_SPACE_BPM(length_cav, VCAV_length2,  n_cav,0.375,"CAVR");
			//primary->addRPG_SPACE_ANGLE(10.0, n_cav,-0.125,"CAVP",0.0, 0.0);
			//primary->addTwoArmRPG_ABS(num_abs, length_abs_dbr,length_abs_abs,length_abs_total,n_abs_bgr,n_abs_cap_layer,n_abs_ar_coating, ANGLE_OF_INCIDENCE, n_cav);
			//primary->addRPG_SPACE_LENS(length_cav/2.0, n_cav,VCAV_focus, VCAV_length1,  0.0);
			//primary->addBoundary(OUTPUT_COUPLING,1.0);
			primary->addBoundary(1.0,1.0);
			//primary->addBoundary(0.0,1.0);
			primary->set_transverse_QW_pump_profile_SuperGaussian(13,GAIN_SPOT_FWHM); // order, FWHM. // Corresponds to a waist w0 = 175um
			primary->set_transverse_QW_temp_profile_file(11249, "transverse_temperature_profile_11249.txt",PUMP_SCALER); // NOTE: sets QW temperature profile relative to a lattice (300K) (material__QW#.config temperature at peak)
			primary->setInitialPulseTransverseFWHM(LENS_SPOT_GAIN); // Initial pulse transverse FWHM

		} else if ((STARTL == 1)&&(STARTR==0))
		{
			cout << "Not tested " << endl;
			exit(-1);
		}
	


		//================================
		// Initialize simulation variables
		//
		int max_work_qw = (primary->getNumberDevices() + primary->getNumberTwoArmDevices()); // Get number of QWs
		if (max_work_qw==0)
		{
			max_work_qw = 1; // Always let 0 work be done by rank=0 worker
		}
		if (MPI_MY_RANK < max_work_qw) // These are the workers that go into a problem
		{
			// Create new work gang
			MPI_Comm work_gang;
			MPI_Comm_split(MPI_COMM_WORLD, 0, MPI_MY_RANK, &work_gang);

			// Do work with new gang
			primary->mpi_set_work_gang(&work_gang);
			
			primary->maxwell_initialize(dt,MPI_MY_RANK); // Set all simulation variables to zero

	//		primary->file_output_device_disable();
	//		primary->file_output_device_reduce(TRANSVERSE_DOMAIN_OUTPUT_FWHM);
			primary->file_output_device_set_freq(TRANSVERSE_DOMAIN_OUTPUT_FREQ, TRANSVERSE_DOMAIN_OUTPUT_FWHM, VCAV_twoArmDeviceOutputLevel, VCAV_deviceOutputLevel);
	//		primary->set_transverse_ABS_fake_focus(3.1623,LENS_SPOT_SESAM);

			if (MPI_MY_RANK==0)
			{
				primary->Print();
				primary->file_output_structure(); // Print structure to file
			}
		
			MPI_Barrier(work_gang); // To ensure syncronization of timers

			//=====================================
			#ifdef MPI_BALANCE_WORKLOAD
			// ! At this moment this is not optimized enough. Default estimates are just as good.
			// Test computational demand of VECSEL nodes before running simulation
			// Will output file 'MPI_CONFIG_WEIGHTS.dat' with information to be used in second run of main program
			if (MPI_MY_RANK==0)
			{
				cout << endl << "Start MPI performance analysis.." << endl << endl;
			}
			primary->setInitialPulseAmplitude(2.0e7);	// A stronger pulse
			primary->setInitialPulseFwhm(100*fs);		// about mode-locked length
			runSimulation(primary, 120*ps, dt, 1.0*ps, out_DT_freq, 130.0*ps, 130.0*ps, MPI_MY_RANK);	

			// Gather and process performance data
			if (MPI_MY_RANK==0)
			{
				cout << "test simulation done, printing timing data to file: MPI_CONFIG_WEIGHTS.dat" << endl;
				double counter_device_data[2*primary->getNumberDevices()];
				primary->get_device_MPI_balance_counters(counter_device_data); // get all device data

				// get Maxwell data
				double maxwell_time, maxwell_time_sdev;
				primary->get_maxwell_MPI_balance_counters(&maxwell_time, &maxwell_time_sdev);
				
				FILE *mpi_balance = fopen("MPI_CONFIG_WEIGHTS.dat","w+");
				if (mpi_balance != NULL)
				{
					fprintf(mpi_balance, "%d\n", primary->getNumberDevices());
					fprintf(mpi_balance,"%.6f %.6f\n",maxwell_time/1e-6, maxwell_time_sdev/1e-6);
					for(int i = 0; i < primary->getNumberDevices(); i++)
					{
						fprintf(mpi_balance,"%.6f %.6f\n",counter_device_data[2*i]/1e-6, counter_device_data[2*i+1]/1e-6);
					}
					fclose(mpi_balance);
					cout << "Performance analysis done, restart without the flag MPI_BALANCE_WORKLOAD" << endl;
				} else {
					cout << "main.cpp:: ERROR writing to MPI_CONFIG_WEIGHTS.dat.." << endl;
					exit(-1);
				}
			} else {
				primary->get_device_MPI_balance_counters(NULL);
			}
			MPI_Barrier(work_gang); // To ensure syncronization of timers
			exit(-1);
			#endif
			//=====================================
			// Run regular program
			runSimulation(primary, t_max, dt, out_DT_write, out_DT_freq,out_DT_wait, save_DT_global,MPI_MY_RANK);	

			#ifdef USE_MAIN_TIMERS
			for(int i = 0; i < max_work_qw; i++) {
				MPI_Barrier(work_gang);
				if (i == MPI_MY_RANK) {
					MainStat->Print_short();
				}
			}
			#endif

			// Free work gang
			MPI_Comm_free(&work_gang);
		} else {
			MPI_Comm slacking_gang;
			MPI_Comm_split(MPI_COMM_WORLD, 1, MPI_MY_RANK, &slacking_gang);
			MPI_Comm_free(&slacking_gang);
		}
		MPI_Barrier(MPI_COMM_WORLD); // Wait for all workers to finish jobs
	}
}

void runSimulation(VECSEL *system, double T_MAX, double DT, double OUT_DT_WRITE, double OUT_DT_FREQ, double OUT_DT_WAIT, double OUT_DT_SAVE,int my_rank)
{
	// Set output data
	C_wait 		= 0;
	C_write		= 0;
	C_freq		= 0;
	COUNT_WAIT  = ceil(OUT_DT_WAIT/DT);
	COUNT_WRITE = ceil(OUT_DT_WRITE/DT);
	COUNT_FREQ  = ceil(OUT_DT_FREQ/DT);
	C_wait 		= COUNT_WAIT;
	out_counter = 0;
	NEW_OUTPUT_FILE = 1;
	
	// Set save variables
	
	COUNT_SAVE = ceil(OUT_DT_SAVE/DT);
	C_save = COUNT_SAVE;
	save_counter = 0;

	// Output time variables
	ofstream calcTime("output_calculation_time.dat",ios::app);
	
	double t = 0;
	int tmp_o = 0;

	#ifdef USE_MAIN_TIMERS
	MainStat->start("runSimulation");
	#endif
	
	struct timeval t1, t2;
	gettimeofday(&t1,NULL);
	
	if (LOAD_FROM_FILE_NR >=0)
	{
		// Load controll variables
		loadMiscVariables(LOAD_FROM_FILE_NR, &t);
		
		// Load Maxwell and Bloch
		system->file_load_variables(LOAD_FROM_FILE_NR,0);

		// Write to output at once
		C_wait 		= COUNT_WAIT;
		C_write		= 0;
		C_freq		= 0;
		//out_counter = 0;
		NEW_OUTPUT_FILE = 1;
		COUNT_WAIT  = ceil(OUT_DT_WAIT/DT);
		COUNT_WRITE = ceil(OUT_DT_WRITE/DT);
		COUNT_FREQ  = ceil(OUT_DT_FREQ/DT);
		COUNT_SAVE  = ceil(OUT_DT_SAVE/DT);
		
		
		if (my_rank==0)
		{
			cout << "LOAD DONE!" << endl;
		}
		//exit(-1);
		
		//calcTime.open("out_calculation_time.dat",ios::app);
	} else {
		// Create file for keeping track of time in iterations only
		// Simulation time recorded | Calculation time used in total
		//calcTime.open("out_calculation_time.dat", ios::trunc | ios::app);
	}

	double num_perc = 100.0;
	double t_curr = 0;
	while (t <= T_MAX)
	{
		#ifdef USE_MAIN_TIMERS
		MainStat->start("VECSEL::iterateModules");
		#endif
		
		system->iterateModules(t,DT);
		#ifdef USE_MAIN_TIMERS
		MainStat->stop("VECSEL::iterateModules");
		#endif
		#ifdef USE_MAIN_TIMERS
		MainStat->start("main file i/o");
		#endif
		
		// Output data
		out_slice(t,system,my_rank);

		#ifdef USE_MAIN_TIMERS
		MainStat->stop("main file i/o");
		#endif


		// Single step done
//		system->file_output_close_all();
//		exit(-1);

		if (my_rank==0)
		{
			if ((t/T_MAX) > tmp_o/num_perc)
			{
				gettimeofday(&t2,NULL);
				double diff 	 = ((double)t2.tv_sec + (double)t2.tv_usec * .000001 - (double)t1.tv_sec - (double)t1.tv_usec * .000001);
				//printf("[%3.1f%%] t_sim = %.2f [ps] t_real = %.2f [s]\n",100.0*tmp_o/10.0,t/ps,diff);
				
				double avgTime   = 100.0*diff/tmp_o;
				//double est_t   = diff/(t/T_MAX);        // Estimate of total time
				double est_t     = avgTime; // Estimate of total time
				double est_t_day = est_t/(60.0*60.0*24.0);
				double est_t_hr  = 24.0*(est_t_day - floor(est_t_day));
				double est_t_min = 60.0*(est_t_hr - floor(est_t_hr));
				double est_t_sec = 60.0*(est_t_min - floor(est_t_min));

				double rem_t     = est_t - diff;
				double rem_t_day = rem_t/(60.0*60.0*24.0);
				double rem_t_hr  = 24.0*(rem_t_day - floor(rem_t_day));
				double rem_t_min = 60.0*(rem_t_hr - floor(rem_t_hr));
				double rem_t_sec = 60.0*(rem_t_min - floor(rem_t_min));


				printf("[%3.1f%%] t_sim = %.2f [ps] t_real = %.2f [s] EST: REMAIN = %.0f:%.0f:%.0f:%.0f [d:hr:min:s] (%.0f [s]) TOTAL = %.0f:%.0f:%.0f:%.0f [d:hr:min:s] (%.0f [s])\n",100.0*tmp_o/num_perc,t/ps,diff,floor(rem_t_day), floor(rem_t_hr),floor(rem_t_min),floor(rem_t_sec), rem_t, floor(est_t_day), floor(est_t_hr),floor(est_t_min),floor(est_t_sec),est_t);
				
				tmp_o 		+= 1;
				calcTime << t/ps << "\t" << diff << endl;
			}
		}
		
		t += DT;

		#ifdef USE_MAIN_TIMERS
		MainStat->start("main file i/o");
		#endif
		
		// Save program to file
		save_slice(t, system, my_rank);
		#ifdef USE_MAIN_TIMERS
		MainStat->stop("main file i/o");
		#endif

	}

	if (my_rank==0)	
	{
		gettimeofday(&t2,NULL);
		double diff = ((double)t2.tv_sec + (double)t2.tv_usec * .000001 - (double)t1.tv_sec - (double)t1.tv_usec * .000001);
		printf("Simulation finished in %.2f [s]\n",diff);
		calcTime << t/ps << "\t" << diff << endl;
		calcTime.close();
	}
	
	// Final snapshot of Cavity
	//system->toFileWriteCavitySnapshot_single(800, 10000);

	#ifdef USE_MAIN_TIMERS
	MainStat->stop("runSimulation");
	#endif
	
	
	if (my_rank==0)
	{
		cout << "Program finished" << endl;
		std::stringstream tmp;
		tmp << system->getToFileOutputKey() << "run_completed.dat";
		std::string tmp2 = tmp.str();
		fclose(fopen(tmp2.c_str(),"w+"));// Write a small file that will signal that the output is done
	}

	system->file_output_close_all();
}


void out_slice(double t, VECSEL *primary,int my_rank)
{
	// If waiting check if done
	if (C_wait >= COUNT_WAIT)
	{
			if (NEW_OUTPUT_FILE == 1)
			{
				primary->file_output_open_all(out_counter);
				NEW_OUTPUT_FILE = 0;
				if (my_rank==0)
				{
					cout << " -> Output " << out_counter << " start" << endl;
				}
			}
			
			// If wringing check how long it has been writing
			if (C_write < COUNT_WRITE)
			{
				// Output every out_t_freq if we are writing
				if (C_freq % COUNT_FREQ == 0)
				{
				primary->file_output_write_all(t);
				}
				C_freq += 1;
				
			} else {
				if (my_rank==0)
				{
					cout << " -> Output " << out_counter << " done" << endl;
				}
				
				// We are done writing, reset timers
				out_counter += 1;				 // Next output number
				
				primary->file_output_close_all();
				NEW_OUTPUT_FILE = 1;
				//primary->file_output_open_all(out_counter);
				
				C_write = 0;
				C_wait = 0;
				
				if (C_freq % COUNT_FREQ == 0)
				{
					C_freq = 0;
				}
			
				if (my_rank==0)
				{
					std::stringstream tmp;
					tmp << primary->getToFileOutputKey() << out_counter-1 << "_output_completed.dat";
					std::string tmp2 = tmp.str();
					fclose(fopen(tmp2.c_str(),"w+"));// Write a small file that will signal that the output is done
				}
		
			}
			
			C_write += 1;
	} 
	C_wait += 1;
}



void out_slice(double t, VECSEL *primary, out_slice_param *primary_output_param)
{
	// If waiting check if done
	//if (C_wait >= COUNT_WAIT)
	if (primary_output_param->isDoneWaiting())
	{
			//if (NEW_OUTPUT_FILE == 1)
			if (primary_output_param->get_new_output_file() ==1)
			{
				primary->file_output_open_all(primary_output_param->get_out_counter());
				primary_output_param->set_new_output_file(0);
				//NEW_OUTPUT_FILE = 0;
				cout << " -> Output " << primary_output_param->get_out_counter() << " start" << endl;
			}
			
			// If writing check how long it has been writing
			//if (C_write < COUNT_WRITE)
			if (primary_output_param->isWriting())
			{
				// Output every out_t_freq if we are writing
				//if (C_freq % COUNT_FREQ == 0)
				if (primary_output_param->writeNow())
				{
					primary->file_output_write_all(t);
				}
				primary_output_param->update_c_freq();
				//C_freq += 1;
				
			} else {
				cout << " -> Ouput " << primary_output_param->get_out_counter() << " done" << endl;
				
				// We are done writing, reset timers
				primary_output_param->update_out_counter();
				//out_counter += 1;				 // Next output number
				
				primary->file_output_close_all();
				primary_output_param->set_new_output_file(1);
				//NEW_OUTPUT_FILE = 1;
				//primary->toFileOpenOutputAll(out_counter);
				
				primary_output_param->set_c_write(0);
				primary_output_param->set_c_wait(0);
				//C_write = 0;
				//C_wait = 0;
				
				//if (C_freq % COUNT_FREQ == 0)
				if (primary_output_param->writeNow())
				{
					//C_freq = 0;
					primary_output_param->set_c_freq(0);
				}

				std::stringstream tmp;
				tmp << primary->getToFileOutputKey() << primary_output_param->get_out_counter()-1 << "_output_completed.dat";
				std::string tmp2 = tmp.str();
				fclose(fopen(tmp2.c_str(),"w+"));// Write a small file that will signal that the output is done
			}
			
			//C_write += 1;
			primary_output_param->update_c_write();
	} 
	//C_wait += 1;
	primary_output_param->update_c_wait();
}


void saveMiscVariables(int save_count, double t_sim)
{
	std::stringstream fileName;
	/*
		Save Misc variables
	*/
	
	// Time
	fileName.str("");
	fileName << "save/save_" << save_count << "_t_sim.dat";
	saveBinary(fileName.str(), &t_sim, 1);
	cout << "SAVE: T_sim           = " << t_sim << endl;
	
	// Counters
	fileName.str("");
	fileName << "save/save_" << save_count << "_C_wait.dat";
	saveBinary(fileName.str(), &C_wait, 1);
	cout << "SAVE: C_wait          = " << C_wait << endl;
	
	fileName.str("");
	fileName << "save/save_" << save_count << "_C_write.dat";
	saveBinary(fileName.str(), &C_write, 1);
	cout << "SAVE: C_write         = " << C_write << endl;
	
	fileName.str("");
	fileName << "save/save_" << save_count << "_C_freq.dat";
	saveBinary(fileName.str(), &C_freq, 1);
	cout << "SAVE: C_freq          = " << C_freq << endl;
	
	fileName.str("");
	fileName << "save/save_" << save_count << "_COUNT_WAIT.dat";
	saveBinary(fileName.str(), &COUNT_WAIT, 1);
	cout << "SAVE: COUNT_WAIT      = " << COUNT_WAIT << endl;
	
	fileName.str("");
	fileName << "save/save_" << save_count << "_COUNT_WRITE.dat";
	saveBinary(fileName.str(), &COUNT_WRITE, 1);
	cout << "SAVE: COUNT_WRITE     = " << COUNT_WRITE << endl;
	
	fileName.str("");
	fileName << "save/save_" << save_count << "_COUNT_FREQ.dat";
	saveBinary(fileName.str(), &COUNT_FREQ, 1);
	cout << "SAVE: COUNT_FREQ      = " << COUNT_FREQ << endl;
	
	fileName.str("");
	fileName << "save/save_" << save_count << "_out_counter.dat";
	saveBinary(fileName.str(), &out_counter, 1);
	cout << "SAVE: out_counter     = " << out_counter << endl;
	
	fileName.str("");
	fileName << "save/save_" << save_count << "_NEW_OUTPUT_FILE.dat";
	saveBinary(fileName.str(), &NEW_OUTPUT_FILE, 1);
	cout << "SAVE: NEW_OUTPUT_FILE = " << NEW_OUTPUT_FILE << endl;

	
	fileName.str("");
	fileName << "save/save_" << save_count << "_save_counter.dat";
	saveBinary(fileName.str(), &save_counter, 1);
	cout << "SAVE: save_counter    = " << save_counter << endl;

	
	fileName.str("");
	fileName << "save/save_" << save_count << "_C_save.dat";
	saveBinary(fileName.str(), &C_save, 1);
	cout << "SAVE: C_save          = " << C_save << endl;

	
	fileName.str("");
	fileName << "save/save_" << save_count << "_COUNT_SAVE.dat";
	saveBinary(fileName.str(), &COUNT_SAVE, 1);
	cout << "SAVE: save_count      = " << save_count << endl;

}

void save_slice(double t, VECSEL *primary,int my_rank)
{
	// If waiting check if done
	if (C_save >= COUNT_SAVE)
	{
		// If save folder does not exist
		if (!dirExists("save"))
		{
			dirMake("save");
		}
		
		if (my_rank==0)
		{	
			// Save controll program
			saveMiscVariables(save_counter, t);
		}
			
		// Save Maxwell and Bloch
		primary->file_save_variables(save_counter,0);
			
		if (my_rank==0)
		{
			cout << " -> Save " << save_counter << endl;
		}
		save_counter++;
		C_save = 0;
	} 
	C_save += 1;
}

void save_slice(double t, VECSEL *primary, out_slice_param *primary_save_param)
{
	
	// If waiting check if done
	//if (C_save >= COUNT_SAVE)
	if (primary_save_param->isSaving())
	{
			// Save controll program
			saveMiscVariables(primary_save_param->get_save_counter(), t);
			
			// Save Maxwell and Bloch
			primary->file_save_variables(primary_save_param->get_save_counter(),0);
			
			cout << " -> Save " << primary_save_param->get_save_counter() << endl;
			
			//save_counter++;
			//C_save = 0;
			primary_save_param->update_save_counter();
			primary_save_param->set_c_save(0);
	} 
	//C_save += 1;
	primary_save_param->update_c_save();
}

void loadMiscVariables(int save_count, double *t_sim)
{
	if (!dirExists("save"))
	{
		cout << "Cannot find save directory: save" << endl;
		exit(-1);
	}
	std::stringstream fileName;
	double tmpd = 0;
	int tmpi = 0;
	//
	//	Load Misc variables
	//
	// Time
	cout << "LOAD MAIN: Load file = " << save_count << endl;
	tmpd = *t_sim;
	fileName << "save/save_" << save_count << "_t_sim.dat";
	loadBinary(fileName.str(), t_sim, 1);
	cout << "LOAD: T_sim           = " << *t_sim << "(" << tmpd << ")" << endl;
	
	// Counters
	tmpi = C_wait;
	fileName.str("");
	fileName << "save/save_" << save_count << "_C_wait.dat";
	loadBinary(fileName.str(), &C_wait, 1);
	cout << "LOAD: C_wait          = " << C_wait << "(" << tmpi << ")" << endl;
	
	tmpi = C_write;
	fileName.str("");
	fileName << "save/save_" << save_count << "_C_write.dat";
	loadBinary(fileName.str(), &C_write, 1);
	cout << "LOAD: C_write         = " << C_write << "(" << tmpi << ")" << endl;
	
	tmpi = C_freq;
	fileName.str("");
	fileName << "save/save_" << save_count << "_C_freq.dat";
	loadBinary(fileName.str(), &C_freq, 1);
	cout << "LOAD: C_freq          = " << C_freq << "(" << tmpi << ")" << endl;
	
	tmpi = COUNT_WAIT;
	fileName.str("");
	fileName << "save/save_" << save_count << "_COUNT_WAIT.dat";
	loadBinary(fileName.str(), &COUNT_WAIT, 1);
	cout << "LOAD: COUNT_WAIT      = " << COUNT_WAIT << "(" << tmpi << ")" << endl;
	
	tmpi = COUNT_WRITE;
	fileName.str("");
	fileName << "save/save_" << save_count << "_COUNT_WRITE.dat";
	loadBinary(fileName.str(), &COUNT_WRITE, 1);
	cout << "LOAD: COUNT_WRITE     = " << COUNT_WRITE << "(" << tmpi << ")" << endl;
	
	tmpi = COUNT_FREQ;
	fileName.str("");
	fileName << "save/save_" << save_count << "_COUNT_FREQ.dat";
	loadBinary(fileName.str(), &COUNT_FREQ, 1);
	cout << "LOAD: COUNT_FREQ      = " << COUNT_FREQ << "(" << tmpi << ")" << endl;
	
	tmpi = out_counter;
	fileName.str("");
	fileName << "save/save_" << save_count << "_out_counter.dat";
	loadBinary(fileName.str(), &out_counter, 1);
	cout << "LOAD: out_counter     = " << out_counter << "(" << tmpi << ")" << endl;
	
	tmpi = NEW_OUTPUT_FILE;
	fileName.str("");
	fileName << "save/save_" << save_count << "_NEW_OUTPUT_FILE.dat";
	loadBinary(fileName.str(), &NEW_OUTPUT_FILE, 1);
	cout << "LOAD: NEW_OUTPUT_FILE = " << NEW_OUTPUT_FILE << "(" << tmpi << ")" << endl;
	NEW_OUTPUT_FILE = 1;
	
	tmpi = save_counter;
	fileName.str("");
	fileName << "save/save_" << save_count << "_save_counter.dat";
	loadBinary(fileName.str(), &save_counter, 1);
	cout << "LOAD: save_counter    = " << save_counter << "(" << tmpi << ")" << endl;
	
	tmpi = C_save;
	fileName.str("");
	fileName << "save/save_" << save_count << "_C_save.dat";
	loadBinary(fileName.str(), &C_save, 1);
	cout << "LOAD: C_save          = " << C_save << "(" << tmpi << ")" << endl;
	
	tmpi = COUNT_SAVE;
	fileName.str("");
	fileName << "save/save_" << save_count << "_COUNT_SAVE.dat";
	loadBinary(fileName.str(), &COUNT_SAVE, 1);
	cout << "LOAD: save_count      = " << save_count << "(" << tmpi << ")" << endl;

}




/*
	Remove whitespaces
	From SelectricSimian @ stackOverflow
	http://stackoverflow.com/questions/14233065/remove-whitespace-in-stdstring
*/

void removeWhitespace(std::string& str) 
{
    for (size_t i = 0; i < str.length(); i++) {
        if (str[i] == ' ' || str[i] == '\n' || str[i] == '\t') {
            str.erase(i, 1);
            i--;
        }
    }
}

double getDoubleFromLine(std::string line)
{
	string line2 = line;
	std::string mystr = line.substr(0, line.find("//", 0));
	//trim(mystr);		// Using boost
	removeWhitespace(mystr);
	
	if (mystr.empty())
	{
		cout << "Import failed, no value found in line" << endl;
		cout << "Line: |" << line2 << "|" << endl;
		exit(-1);
	}
	
	return atof(mystr.c_str());
}

int getIntFromLine(std::string line)
{
	string line2 = line;
	std::string mystr = line.substr(0, line.find("//", 0));
	//trim(mystr);		// Using boost
	removeWhitespace(mystr);
	
	if (mystr.empty())
	{
		cout << "Import failed, no value found in line" << endl;
		cout << "Line: |" << line2 << "|" << endl;
		exit(-1);
	}
	
	return atoi(mystr.c_str());
}

void readInValuesFromFile(double *REFLECTION, double *ANGLE, double *INIT_FWHM, double *INIT_AMP, double *INIT_ENR_SHIFT)
{
	std::ifstream cavity_config;
	std::stringstream fileName;
	fileName << "cavity.config";
	// Check if file exists
	if (!fileExists(fileName.str()))
	{
		cout << "readInValuesFromFile(): Cannot find file: " << fileName.str() << endl;
	}
	cavity_config.open((fileName.str()).c_str(), std::ofstream::in); 
	string line,line2;
	if (cavity_config)
	{
		getline(cavity_config,line); // Get first line
		*REFLECTION = getDoubleFromLine(line); 		// Read in amplitude reflection
		*REFLECTION = (*REFLECTION)*(*REFLECTION);	// Make into intensity reflection
		
		getline(cavity_config,line);
		*ANGLE = getDoubleFromLine(line)*Pi/180.0;	// Read in angle of incidence
		
		getline(cavity_config,line);
		*INIT_FWHM = getDoubleFromLine(line)*fs;	// Read in fwhm of initial pulse
		
		getline(cavity_config,line);
		*INIT_AMP = getDoubleFromLine(line);		// Read in amplitude of pulse
		
		getline(cavity_config,line);
		*INIT_ENR_SHIFT = getDoubleFromLine(line);		// Read in energy shift of pulse [meV]
		*INIT_ENR_SHIFT = *INIT_ENR_SHIFT*e/(1000.0*hbar);
		
		cout << "Load from file: " << fileName.str() << endl;
		
		cout << "reflection        		= " << sqrt(*REFLECTION) << " (amplitude)" << endl;
		cout << "angle of incidence		= " << (*ANGLE)*180.0/Pi << " [deg]" << endl;
		cout << "initial fwhm			= " << (*INIT_FWHM)/fs << " [fs]" << endl;
		cout << "initial amplitude		= " << (*INIT_AMP)*um << " [V/um]" << endl;
		double energy = (0.5*eps0*c0)*(4.0*(*INIT_AMP)*(*INIT_AMP)/(3.0*1.76275/(*INIT_FWHM)));
		cout << "=> pulse energy (air)	= " << energy*cm*cm/um << " [uJ/cm^2]" << endl;
		cout << "initial energy shift	= " << 1000.0*hbar*(*INIT_ENR_SHIFT)/e << " [meV]"  << endl;
		
	} else {
	
		cout << "readInValuesFromFile(): could not open file: " << fileName.str() << endl; 
		exit(-1);
	}
	
	cavity_config.close();
}


