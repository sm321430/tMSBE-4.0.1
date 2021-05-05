
#ifndef __VECSEL_H_INCLUDED__
#define __VECSEL_H_INCLUDED__

#include <fstream>
#include <string>
#include <vector>
#include "classModule.h"
#include <mpi.h>
#include "setup_simulation_variables.h"
#include "classParallelScheduler.h"
#ifdef USE_MAIN_TIMERS
	#include "myTimer.cpp"
#endif


#if defined(__ICC) || defined(__INTEL_COMPILER)
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#elif defined(__GNUC__) || defined(__GNUG__)
extern "C"{
        // Y= a*X + Y
        void cblas_zaxpy(const int N, std::complex<double>* alpha, std::complex<double>* X, const int LDX, std::complex<double>* Y, const int INCY);

        // Y = X
        void cblas_zcopy(const int N, std::complex<double>* X, const int LDX, std::complex<double>* Y, const int INCY);

        // X = a*X
        void cblas_zscal(const int N, std::complex<double>* alpha, std::complex<double>* X, const int LDX);
}
#endif

//! The primary class for a laser cavity
/*! This class hold all structural information about the simulation domain and iteration of the transfer matrix method. \n
    A VECSEL cavity is buildt by stacking modules in sequence. Each module can be a cavity, a QW, a loss element, a boundary, etc.. \n
    Examples of how to run a simple simulation can be found in main.cpp, however the VECSEL class contains information about various gain structure, DBRs, etc. that can be included into a cavity. \n


   \sa Module
   \sa Device
   \sa TwoArmDevice
   \sa TwoArmInterface
   \sa TwoArmCavity
   \sa Cavity
   \sa Boundary
   \sa LossElement \n

   2018: Added more comments \n
   Isak Kilen

   2020: Added two-arm non-normal implementation \n
   Sam McLaren

*/
class VECSEL
{
	public:
		//! A constructor
		VECSEL();
		//! A destructor
		~VECSEL()
		{
			if (MPI_WORK_DIST != NULL)
			{
				for(int i =0; i < MPI_WORK_GANG_SIZE; i++)
				{
					delete [] MPI_WORK_DIST[i];
				}
				delete [] MPI_WORK_DIST; // size(MPI_WORK_DIST) = number_of_workers x 2. Contains start and stop index of work
				delete [] MPI_WORK_DIST_E_OFFSET;
				delete [] MPI_WORK_DIST_E_SIZE ;
				delete [] MPI_WORK_DIST_P_OFFSET;
				delete [] MPI_WORK_DIST_P_SIZE ;
				delete [] MPI_WORK_DIST_E_LOCAL;
				delete [] MPI_WORK_DIST_P_LOCAL;
				delete [] MPI_WORK_DEVICE_SIZE;
				delete [] MPI_LoadBalancer_index_set;
				delete [] MPI_LoadBalancer_P_tmp;
				delete [] MPI_LoadBalancer_E_tmp;
				delete [] MPI_WORK_DEVICE_OFFSET;
				delete [] MPI_WORK_DEVICE_TIMING_SIZE;
				delete [] MPI_WORK_DEVICE_TIMING_OFFSET;

				if (MPI_MY_RANK==0)
				{
					for(int i =0; i < 3; i++)
					{
						delete [] MPI_WORK_DIST_P_GLOBAL[i];
					}
					delete [] MPI_WORK_DIST_E_GLOBAL;
				}
				delete [] MPI_WORK_DIST_P_GLOBAL;

				delete MPI_LoadBalancer;

				delete [] VECSEL_transverse_points_y;
				delete [] VECSEL_transverse_points_device_y;
				delete [] VECSEL_initial_transverse_pulse_profile;
				delete [] VECSEL_initial_temp_profile_T;

				delete [] cavity_trans_E_pl;
				delete [] cavity_trans_E_mi;
				delete [] cavity_trans_E_pl_k;
				delete [] cavity_trans_E_mi_k;
				delete [] cavity_trans_MacPol;
			}

			#ifdef USE_CAVITY_SNAPSHOT
			if (VECSEL_cav_snapshot_E != NULL)
			{
				delete [] VECSEL_cav_snapshot_E;
				delete [] VECSEL_cav_snapshot_E_re;
				delete [] VECSEL_cav_snapshot_E_im;
			}
			#endif

			#ifdef MPI_BALANCE_WORKLOAD
			delete MPI_load;
			#endif
			
			
			
		}
		//! A copy-constructor
		VECSEL(const VECSEL &obj)
		{
			VECSEL_name = obj.VECSEL_name;
			VECSEL_output_key = obj.VECSEL_output_key;
			modules = obj.modules;
			
			
			VECSEL_cav_snapshot_E = NULL;
			VECSEL_cav_snapshot_E_re = NULL;
			VECSEL_cav_snapshot_E_im = NULL;
			VECSEL_cav_snapshot_index = obj.VECSEL_cav_snapshot_index;
			VECSEL_cav_snapshot_x = obj.VECSEL_cav_snapshot_x;
			VECSEL_cav_snapshot_num_points = obj.VECSEL_cav_snapshot_num_points;
			VECSEL_cav_snapshot_output_count = obj.VECSEL_cav_snapshot_output_count;
			VECSEL_cav_snapshot_output_wait = obj.VECSEL_cav_snapshot_output_wait;
			
			number_cavities = obj.number_cavities;
			number_twoArmCavities = obj.number_twoArmCavities;
			number_birefringentCrystals = obj.number_twoArmCavities;
			number_twoArmInterfaces = obj.number_twoArmInterfaces;
			number_devices = obj.number_devices;
			number_twoArmDevices = obj.number_twoArmDevices;
			number_boundaries = obj.number_boundaries;

			VECSEL_initial_fwhm = obj.VECSEL_initial_fwhm;
			VECSEL_initial_delay = obj.VECSEL_initial_delay;
			VECSEL_initial_amplitude = obj.VECSEL_initial_amplitude;
			VECSEL_initial_energy_shift = obj.VECSEL_initial_energy_shift;
			VECSEL_initial_transverse_FWHM = obj.VECSEL_initial_transverse_FWHM;
			VECSEL_DT = obj.VECSEL_DT;
			VECSEL_ROUND_TRIP_TIME = obj.VECSEL_ROUND_TRIP_TIME;
			VECSEL_ROUND_TRIP_ITERATIONS = obj.VECSEL_ROUND_TRIP_ITERATIONS;
			VECSEL_QW_FEEDBACK = obj.VECSEL_QW_FEEDBACK;
			VECSEL_initial_pump_profile_SG_degree = obj.VECSEL_initial_pump_profile_SG_degree;
			VECSEL_initial_pump_profile_SG_sigma = obj.VECSEL_initial_pump_profile_SG_sigma;
			VECSEL_initial_pump_profile_SG_FWHM = obj.VECSEL_initial_pump_profile_SG_FWHM;
			VECSEL_initial_temp_profile_SG_degree = obj.VECSEL_initial_temp_profile_SG_degree;
			VECSEL_initial_temp_profile_SG_sigma = obj.VECSEL_initial_temp_profile_SG_sigma;
			VECSEL_initial_temp_profile_SG_FWHM = obj.VECSEL_initial_temp_profile_SG_FWHM;

			test_VECSEL_iteration = obj.test_VECSEL_iteration;
			init_VECSEL_iteration = obj.init_VECSEL_iteration;
			filter_diagnostics_prev_E = obj.filter_diagnostics_prev_E;

			quick_index_cavity = obj.quick_index_cavity;
			quick_index_twoArmCavity = obj.quick_index_twoArmCavity;
			quick_index_twoArmInterface = obj.quick_index_twoArmInterface;
			quick_index_twoArmPostCav = obj.quick_index_twoArmPostCav;
			quick_index_cavity_noBPM = obj.quick_index_cavity_noBPM;
			quick_index_cavity_lens = obj.quick_index_cavity_lens;
			quick_index_cavity_lens_halfCav = obj.quick_index_cavity_lens_halfCav;
			quick_index_boundary = obj.quick_index_boundary;
			quick_index_totalDevice = obj.quick_index_totalDevice;
			quick_index_device = obj.quick_index_device;
			quick_index_twoArmDevice = obj.quick_index_twoArmDevice;
			quick_index_device_previous_cavity = obj.quick_index_device_previous_cavity;
			quick_index_twoArmDevice_previous_cavity = obj.quick_index_twoArmDevice_previous_cavity;
			quick_index_birefringentCrystal = obj.quick_index_birefringentCrystal;

			VECSEL_lambda = obj.VECSEL_lambda;
			VECSEL_pulse_start_l = obj.VECSEL_pulse_start_l;
			VECSEL_pulse_start_r = obj.VECSEL_pulse_start_r;
			
			

			cavity_trans_E_pl = NULL;
			cavity_trans_E_mi = NULL;
			cavity_trans_E_pl_k = NULL;
			cavity_trans_E_mi_k = NULL;
			cavity_trans_MacPol = NULL;


			if (obj.MPI_WORK_DIST != NULL)
			{
				mpi_initialize_work_arrays(obj.MPI_WORK_DIST_TOTAL);
				MPI_MY_RANK = obj.MPI_MY_RANK;
				mpi_set_work_gang(obj.MPI_WORK_GANG);
			}

			if (obj.VECSEL_transverse_points_y != NULL)
			{
				VECSEL_transverse_points_y = new double[obj.VECSEL_transverse_points_number];
				for(int i = 0; i < obj.VECSEL_transverse_points_number; i++)
				{
					VECSEL_transverse_points_y[i] = obj.VECSEL_transverse_points_y[i];
				}
				
				cavity_trans_E_pl = new std::complex<double>[obj.VECSEL_transverse_points_number];
				cavity_trans_E_mi = new std::complex<double>[obj.VECSEL_transverse_points_number];
				cavity_trans_E_pl_k = new std::complex<double>[obj.VECSEL_transverse_points_number];
				cavity_trans_E_mi_k = new std::complex<double>[obj.VECSEL_transverse_points_number];
				cavity_trans_MacPol = new std::complex<double>[obj.VECSEL_transverse_points_number];
			}
		
			if (obj.VECSEL_transverse_points_device_y != NULL)
			{
				VECSEL_transverse_points_device_y = new double[obj.VECSEL_transverse_points_device_number];
				for(int i = 0; i < obj.VECSEL_transverse_points_device_number; i++)
				{
					VECSEL_transverse_points_device_y[i] = obj.VECSEL_transverse_points_device_y[i];
				}
			}

			if (obj.VECSEL_initial_transverse_pulse_profile != NULL)
			{
				VECSEL_initial_transverse_pulse_profile = new std::complex<double>[obj.VECSEL_transverse_points_number];
				for(int i = 0; i < obj.VECSEL_transverse_points_number; i++)
				{
					VECSEL_initial_transverse_pulse_profile[i] = obj.VECSEL_initial_transverse_pulse_profile[i];
				}

				
			}
			if (obj.VECSEL_initial_temp_profile_T != NULL)
			{
				VECSEL_initial_temp_profile_T = new double[obj.VECSEL_transverse_points_number];
				for(int i = 0; i < obj.VECSEL_transverse_points_number; i++)
				{
					VECSEL_initial_temp_profile_T[i] = obj.VECSEL_initial_temp_profile_T[i];
				}
			}

			#ifdef USE_CAVITY_SNAPSHOT
			if (VECSEL_cav_snapshot_E != NULL)
			{
				VECSEL_cav_snapshot_E 	= new std::complex<double>[VECSEL_cav_snapshot_num_points];
				VECSEL_cav_snapshot_E_re = new double[VECSEL_cav_snapshot_num_points];
				VECSEL_cav_snapshot_E_im = new double[VECSEL_cav_snapshot_num_points];

				for(int i = 0; i < VECSEL_cav_snapshot_num_points; i++)
				{
					VECSEL_cav_snapshot_E[i] = obj.VECSEL_cav_snapshot_E[i];
					VECSEL_cav_snapshot_E_re[i] = obj.VECSEL_cav_snapshot_E_re[i];
					VECSEL_cav_snapshot_E_im[i] = obj.VECSEL_cav_snapshot_E_im[i];
				}
			}
			#endif

			#ifdef MPI_BALANCE_WORKLOAD
			if (MPI_load != NULL)
			{
				MPI_load = new myTimer(*obj.MPI_load);
			}
			#endif

		}
		//! Print VECSEL information to screen
		void Print( void ) const;
		
		//! Set name of VECSEL cavity
		void setName(const std::string &);
	
		//! Return name of VECSEL cavity
		std::string getName(void) const;

		//! Set output file name key
		/*! Appears at the beginning of the output file name of ALL VECSEL output.*/
		void setToFileOutputKey(const std::string &);

		//! Return output file name key
		std::string getToFileOutputKey(void) const;

		//! Return number of modules in VECSEL cavity
		int getNumberModules(void) const;

		//! Return a pointer to a given module #
		Module * getModule(int);
		
		//!  Return the central wavelength for the simulation
		double getLambda(void) const;
		//! Set the central wavelength for the simulation
		void setLambda(double);
		
		//! Set counter for number of cavities in VECSEL
		void setNumberCavities(int);
		//! Return counter for number of cavities in VECSEL
		int  getNumberCavities(void) const;
		//! Set counter for number of devices in VECSEL
		void setNumberDevices(int);
		//! Return counter for number of devices in VECSEL
		int  getNumberDevices(void) const;
		//! Set counter for number of two arm interfaces in VECSEL
		void setNumberTwoArmInterfaces(int);
		//! Return counter for number of two arm interfaces in VECSEL
		int  getNumberTwoArmInterfaces(void) const;
		//! Set counter for number of two arm cavities in VECSEL
		void setNumberTwoArmCavities(int);
		//! Return counter for number of two arm cavities in VECSEL
		int  getNumberBirefringentCrystals(void) const;
		//! Set counter for number of  birefringent crystals in VECSEL
		void setNumberBirefringentCrystals(int);
		//! Return counter for number of  birefringent crystals in VECSEL
		int  getNumberTwoArmCavities(void) const;
		//! Set counter for number of two arm devices in VECSEL
		void setNumberTwoArmDevices(int);
		//! Return counter for number of two arm devices in VECSEL
		int  getNumberTwoArmDevices(void) const;
		//! Set counter for number of boundaries in VECSEL
		void setNumberBoundaries(int);
		//! Return counter for number of boundaries in VECSEL
		int  getNumberBoundaries(void) const;
		//! Set reflection coefficient for right boundary
		void setRightBoundaryRef(double);

		//! Return nuber of transverse points
		int getTransversePoints();

		//! Return size of transverse dimsion
		double getTransverseRmax();

		//! Return ratio of domain used for boundary guard (super gaussian)
		double getTransverseBoundaryGuardRatio();
		
		//===========================
		// Cavity functions
		
		//! Create a DBR on the left side with lambda/4 thicknesses
		void addDBR_LEFT(double, double, int, double,double,double);
		//! Create a DBM on the left side with lambda/4 thicknesses
		void addDBM_LEFT(double, double, double, int, int,double, double, double);
		//! Create a DBR on the left side with given thicknesses
		void addDBR_LEFT_nL(double, double, double, double, int, double, double);
		//! Create a DBR on the right side
		void addDBR_RIGHT(double, double, int, double);
		//! Create a signal DBR on the left side with a pump DBR on the back
		/*! Add in a DBR design for lambda_pump, that has a total length of lambda/4.*/
		void addDBR_LEFT_PUMP(double, double, double, double, int, double, double, double);
		//! Create a DBR on the left side with custom_lambda/4 thickneses
		void addDBR_LEFT_LAMBDA(double, double, double, int, double);
		//! Create a non-normal TwoArmDBR from the front (left to right)
		void addTwoArmDBR_FRONT(double, double, int, double,double,double);
		// Custom gain structures
		//! Create a clustered QW structure on the left side
		/*! This structure clusters the QWs with equal spacing on top of antinodes e.x. MQW10, MQW44 \n 
		    No DBR is added on the left! \n
		 \param numQW Number of QW's in the medium [0,inf). If numQW = 0 => Empty medium of length cavityLength
		 \param dx0_qw Distance from LEFTMOST QW's to the edge in units of WHOLE PEAKS
		 \param dx1_qw Distance from RIGHTMOST QW's to right edge in units of WHOLE PEAKS
		 \param cavityIndex Refractive background index in medium
		 \param capIndex Refractive index of CAP layer, if capIndex = 0 then NO CAP layer
		 \param arIndex Refractive index of AR coating, if arIndex = 0 then NO ar coating
		 \param fillerLength Length of cavity to fill in AFTER AR and CAP layer in units of m
		 \param angle_of_incidence Angle in radians for TE
		 \param external_index Index in incoming medium
		*/
		void addCUSTOM_CLUSTER_QW_LEFT(int, int *, int, double, double, double, double, double, double);
		//! Create a clustered QW structure on the left side with pattern 121+booster
		void addCUSTOM_CLUSTER_121_BOOSTER_QW_LEFT(int, int *, bool, int, double, double, double, double, double, double);
		//! Add custom QW placement on the left with dimensions given in the array AND include ANGLE OF INCIDENCE
		void addCUSTOM_QW_LEFT_OPTIMAL_ANGLE(int, double *, double, double, double, double,double,bool,double,double,double,double,double);
		
		// Resonant Periodic Gain structures
		//! Create RPG gain structure on the left side (debug function)
		void addRPG_QW_LEFT_debug(int,double,double,double,double,double,double,double,double);
		//! Create RPG gain structure on the left side
		void addRPG_QW_LEFT(int,double,double,double,double,double,double,double,double);
		//! Create RPG gain structure on the left side with only 1 qw at a given location
		void addRPG_QW_LEFT_EFF(int, int,double,double,double,double,double,double,double,double);
		//! Create RPG gain structure on the left side with only 1 qw at a given location and no cap layers. Intended for debugging purposes.
		void addRPG_QW_simple(int, double, double, double, double, double);
		
		//! Creates simple two arm structure with interface and single cavity. Can be added to.
		//void addTwoArmCavity(double, double, double, double);
		//! Create basic Structure with two cavities separated by a single QW.
		void addQW_STRUCT(int, double, double, double, double, double); 
		//! Create basic twoArmStructure with interface, and two twoArmCavities separate by a single QW.
		void addTwoArmQW_STRUCT(int, double, double, double, double, double, double); 
		//! Create basic twoArmStructure without interface, and two twoArmCavities separate by a single QW.
		void addTwoArmQW_noInterface(int, double, double, double, double, double); 
		//! Create RPG gain structure with single effective QW for non-normal two arm cavity. Placement and material file are distinct
		void addTwoArmRPG_QW_FRONT_EFF(int, int, int, double, double, double, double, double, double, double, double, double);		
		//! Create simple single ABS on the right, no DBR
		void addRPG_ABS_RIGHT(int,double,double,double,double,double,double);
		//! Create simple two arm structure with angle of incidence
		void addBirefringentCrystal(double, double,double,double);
		//! Create simple two arm structure with angle of incidence
		void addTwoArm_Space(double,double,double,double, const std::string &);
		//! Create simple single ABS as two arm structure with angle of incidence
		void addTwoArmRPG_ABS(int,double,double,double,double,double,double,double,double);
		//! Create simple single ABS as two arm structure with angle of incidence starting from left side
		void addTwoArmRPG_ABS_LEFT(int,double,double,double,double,double,double,double,double);
		//! Create a single ABS on the right with a DBR and ANGLE OF INCIDENCE
		void addRPG_ABS_DBR_PHASE_RIGHT_ANGLE(double, double, double, double, double, double, double,double,double, int, double, double, double, double);
		//! Create a solid cavity with all QW on the left, all ABS on the right
		void addQW_SPACE_ABS(int,int,double,double,double,double);
		
		// Some space functions
		//! Add space with a probe in the middle of the region to print output to file.
		void addRPG_SPACE(double, double,double,const std::string &);
		//! Add space with a probe in the middle of the region to print output to file and forced BPM.
		void addRPG_SPACE_BPM(double,double, double,double,const std::string &);
		//! Add space with a probe in the middle of the region to print output to file (disabled currently) and forced BPM at an angle.
		void addRPG_SPACE_BPM_ANGLE(double,double, double,double,const std::string &, double);
		//! Add space+ETALON with a probe to print output to file.
		void addRPG_SPACE_ETALON(double, double,double,const std::string &, double, double);
		//! Add space+2xETALON with a probe to print output to file.
		void addRPG_SPACE_ETALON2(double, double,double,const std::string &, double, double, double, double);
		//! Add space at an ANGLE with a probe to print output to file.
		void addRPG_SPACE_ANGLE(double, double,double,const std::string &,double, double);
		//! Add space at an ANGLE with a probe to print output to file. Name goes to second cavity.
		void addRPG_SPACE_ANGLE_LEFT(double, double,double,const std::string &,double, double);
		//! Add cavity with a single lens
		/*! Construct a cavity with: Gain -> Lens1 -> Space -> Lens2 -> SESAM\n
		    \param cavety_length Length of air cavity, should be equal to 2*f1 where f1 is lens focus
		    \param n_bagr Refractive index of cavity
		    \param focus1, armLength, focus2 are cavity elements 	
		*/
		void addRPG_SPACE_LENS(double, double, double, double, double);
		//! Add cavity with a single lens
		/*! Construct a cavity with: Gain -> Lens1 -> Space -> Lens2 -> SESAM\n
		    \param cavety_length Length of air cavity, should be equal to 2*f1 where f1 is lens focus
		    \param n_bgr Refractive index of cavity
		    \param focus1, armLength, focus2 are cavity elements 	
	            \param angle_of_incidence is angle
		*/
		void addRPG_SPACE_LENS_ANGLE(double, double, double, double, double, double);
		//! Add space with a THIN focusing lens
		/*! Construct a cavity with: Gain -> Lens1 -> Space -> Lens2 -> SESAM\n
		    This is an equvalent formulation
		    \param cavety_length Length of air cavity, should be equal to where f1 is lens focus
		    \param n_bagr Refractive index of cavity
		    \param ABCD components from matrix of optical system
		*/
		void addRPG_SPACE_LENS_transform(double, double, double, double, double, double);
		//! Add space+FILTER with a probe to print output to file.
		void addRPG_SPACE_FILTER(double, double,double,const std::string &,double,double,double,double,int);
		//! Add custom space with a probe to print output to file. Can subtract part of the first cavity to adjust total length.
		void addCUSTOM_SPACE(double,double, double,double,const std::string &);
		
		//! Set transverse pump profile to super gaussian shape with amplitude 1 at the peak
		/*! Changes the QW background carrier density of all CURRENT QWs by scaling according to shape function. WARNING: Will be retroactivly applied with appendStructure() command.\n
			\param degree Degree of super gaussian profile. Degree=1 gives regular gaussian.
			\param fwhm FWHM of super guassian [m]
			\sa getSGPumpDegree
			\sa getSGPumpFWHM
		*/
		void set_transverse_QW_pump_profile_SuperGaussian(int, double);


		//! Set transverse temperature profile to super gaussian shape with amplitude 1 at the peak
		/*! Changes the QW background carrier density of all CURRENT QWs by scaling according to shape function. WARNING: Will be retroactivly applied with appendStructure() command.\n
			\param degree Degree of super gaussian profile. Degree=1 gives regular gaussian.
			\param fwhm FWHM of super guassian [m]
			\sa getSGTempDegree
			\sa getSGTempFWHM
		*/
		void set_transverse_QW_temp_profile_SuperGaussian(int, double);

		//! Set transverse temperature profile form file with data in the range [0,1]
		/*! Changes the QW background carrier density of all CURRENT QWs by scaling according to shape function. File is formated with two columns, first column contains transverse coordinate [units of m] and second column contains temperature scale [in range 0->1] WARNING: Will be retroactivly applied with appendStructure() command.\n
			\param size is the length of the file containing information
			\param name of file to load, should be located in run dir
		*/
		void set_transverse_QW_temp_profile_file(int, const std::string &,double);

		//! Set the focus on the ABS with a given FWHM
		/*! This is an alternaltive to testing the focusing on the SESAM by the lens
			\param fake_focus the new fake focus that is set on the ABS
			\param the fake lasing spot FWHM
		*/
		void set_transverse_ABS_fake_focus(double fake_focus, double fake_spot_fwhm);

		//! Return Super Gaussian pump degree
		/*!
			\sa set_transverse_QW_pump_profile_SuperGaussian
		*/
		int getSGPumpDegree();

		//! Return Super Gaussian pump FWHM
		/*!
			\sa set_transverse_QW_pump_profile_SuperGaussian
		*/
		double getSGPumpFWHM();

		//! Return Super Gaussian temperature degree
		/*!
			\sa set_transverse_QW_temp_profile_SuperGaussian
		*/
		int getSGTempDegree();

		//! Return Super Gaussian temperature FWHM
		/*!
			\sa set_transverse_QW_temp_profile_SuperGaussian
		*/
		double getSGTempFWHM();

		//! Add empty cavity to the right side of the VECSEL domain
		Module * addCavity();
		//! Add Cavity of givien width and refractive index to the right side of the VECSEL domain
		Module * addCavity(double, std::complex<double>);
		//! Add Cavity of given width, refractive index, and angle of incidence on the rigth side of the VECSEL domain
		Module * addCavity(double, std::complex<double>, double, double);
		//! Add Cavity of given width, refractive index, and cos(th)
		Module * addCavity_cos(double, std::complex<double>, double, double);
		//! Add empty Device on the right side of the VECSEL domain
		Module * addDevice();
		//! Add empty two arm cavity to the right side of the VECSEL domain
		Module * addTwoArmCavity();
		//! Add TwoArmCavity of given width and refractive index to the right side of the VECSEL domain
		Module * addTwoArmCavity(double, std::complex<double>);
		//! Add birefringent crystal of given width, refractive indices, and angle of incidence on the rigth side of the VECSEL domain
		Module * addBirefringentCrystal(double, std::complex<double>, double);
		//! Add TwoArmCavity of given width, refractive index, and angle of incidence on the rigth side of the VECSEL domain
		Module * addTwoArmCavity(double, std::complex<double>, double, double);
		//! Add TwoArmCavity of given width, refractive index, and cos(th)
		Module * addTwoArmCavity_cos(double, std::complex<double>, double, double);
		//! Add empty TwoArmDevice on the right side of the VECSEL domain
		Module * addTwoArmDevice();

		//! Add empty two arm interface to the right side of the VECSEL domain
		Module * addTwoArmInterface();
		//! Add TwoArmInterface of given angle and refractive index to the right side of the VECSEL domain. Utilizes interference phases for given refractive index and angles
		Module * addTwoArmInterface(std::complex<double>, double, double, double);
		
		//! Add empty Boundary on the right side of the VECSEL domain
		Module * addBoundary();
		//! Add a Boundary on the right side of the VECSEL domain with given reflection and refractive index
		Module * addBoundary(double,double);
		//! Add empty LossElement on the right side of the VECSEL domain
		Module * addLossElement();
		//! Add a LossElement on the right side of the VECSEL domain with given loss in left/right directions
		Module * addLossElement(double,double);
		//! Add a Filter with given number of elements
		Module * addFilter(int);

		//! Add a given GAIN structrue to the current setup. This is usefull when re-using same design multiple times.
		/*! Is designed to be used with GAIN structures and will call set_QW_transverse_pump_profile_SuperGaussian() with the current degree/FWHM parameters on all QWs.\n */
		void appendStructure(VECSEL *);
		
		
		//====================
		// Maxwell functions
		
		//! Initialize Maxwell solver and all cavity elements
		void maxwell_initialize(double,int);
		
		//! Iterate all modules (fields, QWs, etc..) one timestep forward
		void iterateModules(double,double);
		//! Iterate a single surface one timestep forward
		/*! Applies the Maxwell boundary conditions to a single surface */
		void iterateModules_updateSingleSurface_transfer_matrix(int);
		void iterateModules_updateSingleSurface_transfer_matrix(Cavity *a, Cavity *b);

		//! Iterate a single surface one timestep forward WITHOUT QWs between
		/*! Applies the Maxwell boundary conditions to a single surface WITHOUT QWs between */
		void iterateModules_updateSingleSurface_transfer_matrix_noQW(Cavity *a, Cavity *b);

		void iterateModules_updateSingleSurface_transfer_matrix_noQW_debug(Cavity *a, Cavity *b);
		/*! Applies the Maxwell boundary conditions to a single two arm surface */
		void iterateModules_updateSingleSurface_TwoArm_transfer_matrix(int);
		void iterateModules_updateSingleSurface_TwoArm_transfer_matrix(TwoArmCavity *a, TwoArmCavity *b);

		//! Iterate a single two arm surface one timestep forward WITHOUT QWs between
		/*! Applies the Maxwell boundary conditions to a single surface WITHOUT QWs between */
		void iterateModules_updateSingleSurface_TwoArm_transfer_matrix_noQW(TwoArmCavity *a, TwoArmCavity *b);

		//! Iterate a front two arm surface one timestep forward WITHOUT QWs between and WITH a birefringent crystal
		/*! Applies the Maxwell boundary conditions to a dual surface WITHOUT QWs between */
		void iterateModules_updateSingleSurface_TwoArm_transfer_matrix_birefringentCrystal_front(TwoArmCavity *a, TwoArmCavity *brc);
	
		//! Iterate a front two arm surface one timestep forward WITHOUT QWs between and WITH a birefringent crystal
		/*! Applies the Maxwell boundary conditions to a dual surface WITHOUT QWs between */
		void iterateModules_updateSingleSurface_TwoArm_transfer_matrix_birefringentCrystal_back(TwoArmCavity *brc, TwoArmCavity *a);

		//! Iterate the back on a two arm cavity structure 
		/*! Applies given reflection to field for two arm cavity*/
		void iterateModules_updateSingleSurface_TwoArm_transfer_matrix_back(TwoArmCavity *a, std::complex<double>);

		//! Iterate a single two arm interface one timestep forward excluding the interface width
		/*! Applies the Maxwell boundary conditions to a single surface WITHOUT QWs between */
		void iterateModules_updateSingleSurface_TwoArm_transfer_matrix_interface(Cavity *a, TwoArmInterface *b, Cavity *c, TwoArmCavity *d);

		//! Iterate a single two arm interface one timestep forward including time delay
		/*! Applies the Maxwell boundary conditions to a single surface WITHOUT QWs between */
		void iterateModules_updateSingleSurface_TwoArm_transfer_matrix_interface_delay(Cavity *a, TwoArmInterface *b, Cavity *c, TwoArmCavity *d);
		
		//! Iterate all fields one timestep forward.
		/*! Applies the Maxwell boundary conditions to each material surface */
		void iterateModules_updateAllSurface(VECSEL *,double);

		//! Set starting point for input field.
		/*! Input pulse can come from left or right in cavity. More configation is possible by editing iterateModules(). Parameter: (1,0) gives input from left, (0,1) gives input from right.*/
		void setStartingPoint(int, int);
		
		//! Set initial field amplitude in units [V/m]
		void setInitialPulseAmplitude(double);
		//! Set initial pulse width in units [s]
		/*! Width is measured as the FWHM of the intenstiy profile. Edit maxwell_initial_E() to change.*/
		void setInitialPulseFwhm(double);
		//! Set initial pulse delay in units of [s]
		void setInitialPulseDelay(double);
		//! Optional: Give inital pulse a energy shift [eV]. Default = 0
		void setInitialPulseEnergyShift(double);
		//! Set initla pulse transverse FWHM [m]
		/*! Gaussian profile assumed */
		void setInitialPulseTransverseFWHM(double);
		
		//! Initial field in cavity.
		/*! The input field is by default a pulse that can be configured with width and amplitude. It is also possible to switch the input field to a dual-wavelength input. \n
			\sa setStartingPoint()
			\sa setInitialPulseAmplitude()
			\sa setInitialPulseFwhm()
			\sa setInitialPulseDelay()
			\sa setInitialPulseEnergyShift()
			\sa setInitialPulseDualFreq()
			\sa setInitialPulseOff()
		*/
		void 	maxwell_initial_E(double, std::complex<double>*);
		
		
		//====================
		// File IO functions
		
		//! File operation: Save all variables to file in ./save/
		/*! \param save_count ID for output file.
		    \param offset summed with save ID. default = 0 
		*/
		void file_save_variables(int,int);
		//! File opeation: Load all saved variables (QW,field,...) from file in ./save/
		/*! \param save_count ID for output file.
		    \param offset summed with save ID. default = 0 
		*/
		void file_load_variables(int,int);
		//! File operation: load QW carriers from a file. With the time code closest to the input
		/*! \param output_name File name that matches a file of the form "output_name__*_t.dat"
		    \param t_T the simulation time [s] that should be loaded. If the given time cannot be found, the closest match is used.
		*/
		void file_load_device_carriers_time(const std::string &, double);
		//! File operation: load QW carriers from the given line from an output file
		/*! \param out_count the output # that should be used
		    \param line_num the line number in the given output file.
		*/
		void file_load_device_carriers_line(int, int);
		//! File operation: load QW carriers from the given line of a given file
		/*! \param output_name base file name to load from
		    \param out_count output file # to load from.
		    \param line_num line # in given file to load from
		*/
		void file_load_device_carriers_line(const std::string &, int, int);
		//! File operation: load QW carriers from an ASCII in folder "neq_dist-line_#" where # is the input line_number. Each QW carrier file should have the name "neq_dist_#_ne.txt" and "neq_dist_#_nh.txt", where # is the QW#=1,2,3,..
		/*! \param line_num the spesific file ID# that should be loaded. */
		void file_load_device_carriers_ascii(int);
		//! File operation: load QW carriers from a file in ./save/
		/*! \param out_count the save file # to load from. */
		void file_load_device_carriers(int);

		//! Disable all output to file from all devices (QWs,filters,..)
		void file_output_device_disable(void);
		//! Reduce the number of variables output to file from all devices (QWs, filters,..)
		/*!
			\param x_max Turn off output from all devices outside this range. Turn on all devices inside region.
		*/
		void file_output_device_reduce(double);
		//! Reduce the number of variables output to file from all devices (QWs, filters,..) by only allowing output at for a given spatial frequency dy. Will output at transverse coordinates y = 0, +/- dy, +/- 2*dy, ...  inside the range [-x_max/2, x_max/2]
		/*!
			\param freq Output frequency in transverse domain
			\param x_max turn off all output outside range [-x_max/2, x_max/2]
			\param twoArmDevice_OutputLevel set at positive for full output
			\param device_OutputLevel set at positive for full output
		*/
		void file_output_device_set_freq(double, double, int, int);
		//! File operation: Write output to file 
		/*! \param t_sim current simulation time*/
		void file_output_write(double);
		//! File operation: Open all output files for writing
		/*! \param out_count file # to open for writing */
		void file_output_open(int);
		//! File operation: Close all output files
		void file_output_close(void);
		//! Ask all components to output to file
		/*! \param t_sim current simulation time */
		void file_output_write_all(double);
		//! Ask all components to open files to prepare for writing
		/*! \param out_count file # to open for writing */
		void file_output_open_all(int);
		//! Ask all components to close all output files
		void file_output_close_all(void);
		//! File operation: Print current cavity to file: 
		/*! Write the given cavity structure to file. Each row represents a module\n
		   Output format:\n
		   1st column: What type are we talking about\n
		    -> 0 is boundary\n
		    -> 1 is cavity\n
		    -> 2 is device\n
		    -> 3 is everything else\n
		   2nd column: Position of first vertex\n
		   3rd column: Position of second vertex\n
		   4th column: Refractive index if a cavity. else -1\n
		 */
		void file_output_structure(void);
		//! Experimental: Write a snapshot of the EM field at ALL locations inside the cavity over time. Is VERY computationaly intensive.
		void file_snapshot_write(double);
		//! Experimental: Write a snapshot of the EM field at ALL locations inside the cavity.
		void file_snapshot_write_single(int, int);

		//================
		// Misc functions
		
		//! Return the cavity field evaluatated at the given device nr
		void 	misc_getFieldAtDevice(int,std::complex<double>*);
		//! Return the cavity field evaluatated at the given TwoArm device nr
		void 	misc_getFieldAtTwoArmDevice(int,std::complex<double>*);
		//! Initialize arrays for MPI workers. Should first run mpi_initialize_rank().
		void 					mpi_initialize_work_arrays(int);
		//! Initialize the MPI rank of the current worker
		void mpi_initialize_rank(int);
		//! Initialize MPI work group. Add current worker to the communication group.
		void mpi_set_work_gang(MPI_Comm *);

		//! Set transverse dimension.
		/*! Will create a position vector in range [-R_max/2, (1/2 - 1/N)R_max] with a given # elements. Each transverse point will be simulated.\n
		   \param numP # of transverse points
		   \param R_max Size of transverse domain
		   \param FWHM of super gaussian Boundary guard in FFT BPM
		*/ 
		void set_transverse_dimensions(int, double, double);
		
		//=====================================================
		// Functions for access to devices and their variables
		
		//! Create a delay in QW computations. Will keep SBE constant until a given time has passed.
		void device_set_computational_delay(double);
		//! Diasble all QWs spontaneous emissions.
		void device_disable_spont_emissions_all();
		//! Set all QWs background carrier density
		void device_set_background_carrier_density(double,int);
		//! Set all QW pump parameters for resonant pump
		/*! \param W0 central frequency
		    \param E0 pump amplitye
		    \param ETA width parameter
  		*/
		void device_set_qw_real_pump_parameters(double, double, double,double);
		
		//! Move a QW to a new position. Should only be used with a pulse AND when no field is present. Will zero out all fields in the surrounding material cavities
		/*! \param num # of QWs to change
		    \param newPos array containingnew positions of QWs
		    \param DT timestep in simulation
		*/
		void moveQWpositions(int num, double *,double);
		//! Return the position [m] of QWs relative to zero.
		/*! \param num # of QWs to retun for 0,1,..,num-1
		    \param newPos array to store the QW positions
		*/
		void getQWpositions(int num, double **);
		
		//=====================
		// Filter functions
		
		//! Set Filter to pass through in both direction. (delay only)
		void setFilter_pluss_minus_pass();
		//! Set Filter for left moving field to pass through (delay) and the right moving field to a double gaussian with given parameters.
		/*! \param wa_s Gaussian 1 central frequency. Relative to central frequency of field.
		    \param wb_s Gaussian 2 central frequency. Relative to central frequency of field.
		    \param width_s FWHM of each Gaussian.
		    \param w0_s 
		*/
		void setFilter_pluss_gauss_minus_pass(double, double, double, double);

		//=====================
		// Diagnostics functions

		//! Copy device performance data to array
		void get_device_MPI_balance_counters(double *);
		void get_maxwell_MPI_balance_counters(double *, double *);

		//! Return the longest possible timestep in the given VECSEL cavity.
		double diagnostics_findMaxTimestep(void);
		//! Set multiplier for Mac.Pol. from QWs. default = 1 
		void diagnostics_set_qw_feedback(double);
		//! Set all cavity fields to zero.
		void diagnostics_zero_all_cavity_fields();
		//! Remove all components in the VECSEL cavity.
		void diagnostics_clear_VECSEL();
	
	private:
		Module * addModule();				// Add device to VECSEL
		std::string VECSEL_name;			// Name of VECSEL
		std::string VECSEL_output_key;	// Start of output name of files
		std::vector<Module> modules;	// All modules

		// TRANSVERSE POINTS
		int VECSEL_transverse_points_number;
		int VECSEL_transverse_points_device_number;
		double VECSEL_transverse_points_R_max;
		double VECSEL_transverse_points_boundary_guard_ratio;
		double *VECSEL_transverse_points_y;
		double *VECSEL_transverse_points_device_y;
		std::complex<double> *VECSEL_initial_transverse_pulse_profile;
		std::complex<double> *cavity_trans_E_pl;
		std::complex<double> *cavity_trans_E_mi;
		std::complex<double> *cavity_trans_E_pl_k;
		std::complex<double> *cavity_trans_E_mi_k;
		std::complex<double> *cavity_trans_MacPol;
		
		// Output time data to file
		std::ofstream output_simulation_time;
		std::ofstream *output_E_real;
		std::ofstream *output_E_imag;
		std::ofstream *output_back_E_real;
		std::ofstream *output_back_E_imag;
		
		std::complex<double> *VECSEL_cav_snapshot_E;
		double *VECSEL_cav_snapshot_E_re;
		double *VECSEL_cav_snapshot_E_im;
		std::vector<int> VECSEL_cav_snapshot_index;
		std::vector<double> VECSEL_cav_snapshot_x;
		int 	VECSEL_cav_snapshot_num_points;
		int 	VECSEL_cav_snapshot_output_count;
		int 	VECSEL_cav_snapshot_output_wait;
		std::ofstream output_Cav_Snapshot_x;
		std::ofstream output_Cav_Snapshot_t;
		std::ofstream *output_Cav_Snapshot_E_real;
		std::ofstream *output_Cav_Snapshot_E_imag;
		
		int number_cavities;
		int number_devices;
		int number_twoArmCavities;
		int number_birefringentCrystals;
		int number_twoArmInterfaces;
		int number_twoArmDevices;
		int number_boundaries;
		
		// Settings for initial pulse
		double VECSEL_initial_fwhm;
		double VECSEL_initial_delay;
		double VECSEL_initial_amplitude;
		double VECSEL_initial_energy_shift;
		double VECSEL_initial_transverse_FWHM;
		double VECSEL_DT;
		double VECSEL_ROUND_TRIP_TIME;
		double VECSEL_ROUND_TRIP_ITERATIONS;
		double VECSEL_QW_FEEDBACK;
		int VECSEL_initial_pump_profile_SG_degree;
		double VECSEL_initial_pump_profile_SG_sigma;
		double VECSEL_initial_pump_profile_SG_FWHM;
		int VECSEL_initial_temp_profile_SG_degree;
		double VECSEL_initial_temp_profile_SG_sigma;
		double VECSEL_initial_temp_profile_SG_FWHM;
		double *VECSEL_initial_temp_profile_T;

		double test_VECSEL_iteration;
		double init_VECSEL_iteration;
		std::complex<double> filter_diagnostics_prev_E;
		
		std::vector<int> quick_index_cavity;
		std::vector<int> quick_index_cavity_noQW;
		std::vector<int> quick_index_cavity_QW;
		std::vector<int> quick_index_cavity_noBPM;
		std::vector<int> quick_index_cavity_freeSpace;
		std::vector<int> quick_index_cavity_lens;
		std::vector<int> quick_index_cavity_lens_halfCav;
		std::vector<int> quick_index_twoArmCavity;
		std::vector<int> quick_index_twoArmInterface;
		std::vector<int> quick_index_twoArmPostCav;
		std::vector<int> quick_index_twoArmCavity_noQW;
		std::vector<int> quick_index_twoArmCavity_QW;
		std::vector<int> quick_index_birefringentCrystal; 
		std::vector<int> quick_index_boundary;
		std::vector<int> quick_index_totalDevice;
		std::vector<int> quick_index_device;
		std::vector<int> quick_index_device_previous_cavity;
		std::vector<int> quick_index_twoArmDevice;
		std::vector<int> quick_index_twoArmDevice_previous_cavity;
		
		int **MPI_WORK_DIST; // size(MPI_WORK_DIST) = number_of_workers x 2. Contains start and stop index of work
		int  *MPI_WORK_DIST_E_OFFSET;
		int  *MPI_WORK_DIST_E_SIZE;
		int  *MPI_WORK_DIST_P_OFFSET;
		int  *MPI_WORK_DIST_P_SIZE;
		int  *MPI_WORK_DEVICE_SIZE;
		int  *MPI_WORK_DEVICE_OFFSET;
		int  *MPI_WORK_DEVICE_TIMING_SIZE;
		int  *MPI_WORK_DEVICE_TIMING_OFFSET;
		int   MPI_WORK_DIST_TOTAL;
		int MPI_MY_RANK;
		int MPI_FLAG;
		std::complex<double> *MPI_WORK_DIST_E_GLOBAL;
		std::complex<double> **MPI_WORK_DIST_P_GLOBAL;
		std::complex<double> *MPI_WORK_DIST_E_LOCAL;
		std::complex<double> *MPI_WORK_DIST_P_LOCAL;
		MPI_Comm *MPI_WORK_GANG;
		int MPI_WORK_GANG_SIZE;
		parSchedule *MPI_LoadBalancer;
		int *MPI_LoadBalancer_index_set;
		std::complex<double> *MPI_LoadBalancer_P_tmp;
		std::complex<double> *MPI_LoadBalancer_E_tmp;

		// Properties of the VECSEL
		double VECSEL_lambda;			// Wavelength of VECSEL
		
		int VECSEL_pulse_start_l;
		int VECSEL_pulse_start_r;

		//! Set QW temperature scale to all current QWs
		/*! Transfer transverse temperature profile to QWs
		*/
		void set_transverse_QW_temp_profile(double *);

		// Storage for MPI_load balancer
		#ifdef MPI_BALANCE_WORKLOAD
		myTimer *MPI_load;
		#endif
};

inline void VECSEL::setName( const std::string & newName)
{
	VECSEL_name = newName;
}

inline std::string VECSEL::getName(void) const
{
	return VECSEL_name;
}

inline void VECSEL::setToFileOutputKey(const std::string & newName)
{
	VECSEL_output_key = newName;
}
inline std::string VECSEL::getToFileOutputKey(void) const
{
	return VECSEL_output_key;
}


inline int VECSEL::getNumberModules() const
{
	return modules.size();
}
inline Module * VECSEL::getModule(int i)
{
	return &(modules[i]);
}

inline double VECSEL::getLambda(void) const
{
	return VECSEL_lambda;
}

inline void VECSEL::setLambda(double newLambda)
{
	VECSEL_lambda = newLambda;
}

inline void VECSEL::setStartingPoint(int start_l, int start_r)
{
	VECSEL_pulse_start_l = start_l;
	VECSEL_pulse_start_r = start_r;
}

inline void VECSEL::setInitialPulseAmplitude(double nA)
{
	VECSEL_initial_amplitude = nA;
}

inline void VECSEL::setInitialPulseFwhm(double nF)
{
	VECSEL_initial_fwhm = nF;
}

inline void VECSEL::setInitialPulseDelay(double nD)
{
	VECSEL_initial_delay = nD;
}

inline void VECSEL::setInitialPulseEnergyShift(double nE)
{
	VECSEL_initial_energy_shift = nE;
}


inline void VECSEL::setNumberCavities(int nc)
{
	number_cavities = nc;
}
inline int  VECSEL::getNumberCavities(void) const
{
	return number_cavities;
}
inline void VECSEL::setNumberDevices(int nd)
{
	number_devices = nd;
}
inline int  VECSEL::getNumberDevices(void) const
{
	return number_devices;
}

inline void VECSEL::setNumberTwoArmCavities(int nc)
{
	number_twoArmCavities = nc;
}
inline int  VECSEL::getNumberTwoArmCavities(void) const
{
	return number_twoArmCavities;
}

inline void VECSEL::setNumberBirefringentCrystals(int nc)
{
	number_birefringentCrystals = nc;
}
inline int  VECSEL::getNumberBirefringentCrystals(void) const
{
	return number_birefringentCrystals;
}

inline void VECSEL::setNumberTwoArmInterfaces(int nc)
{
	number_twoArmInterfaces = nc;
}
inline int  VECSEL::getNumberTwoArmInterfaces(void) const
{
	return number_twoArmInterfaces;
}
inline void VECSEL::setNumberTwoArmDevices(int nd)
{
	number_twoArmDevices = nd;
}
inline int  VECSEL::getNumberTwoArmDevices(void) const
{
	return number_twoArmDevices;
}
inline void VECSEL::setNumberBoundaries(int nb)
{
	number_boundaries = nb;
}
inline int  VECSEL::getNumberBoundaries(void) const
{
	return number_boundaries;
}

inline void VECSEL::setRightBoundaryRef(double newRef)
{
	if (modules.back().isBoundary())
	{
		modules.back().getBoundary()->setRefCoeff(newRef);
	}
}

inline void VECSEL::setInitialPulseTransverseFWHM(double newFWHM)
{
	VECSEL_initial_transverse_FWHM = newFWHM;
}

inline void VECSEL::get_device_MPI_balance_counters(double *data)
{
	double data_local[MPI_WORK_DEVICE_TIMING_SIZE[MPI_MY_RANK]];
	int cnt = 0;
	#ifdef ITERATE_QW	
	for(int j = MPI_WORK_DIST[MPI_MY_RANK][0]-1; j < MPI_WORK_DIST[MPI_MY_RANK][1]; j++)
	{
		int indx = quick_index_totalDevice[j]; // Index of device
		if (modules[indx].isDevice())
		{	
			data_local[2*cnt  ] = modules[indx].getDevice()->get_MPI_load_mean_time();
			data_local[2*cnt+1] = modules[indx].getDevice()->get_MPI_load_mean_time_deviation();
			cnt++;
		} else
		{
			data_local[2*cnt  ] = modules[indx].getTwoArmDevice()->get_MPI_load_mean_time();
			data_local[2*cnt+1] = modules[indx].getTwoArmDevice()->get_MPI_load_mean_time_deviation();
			cnt++;
		}
	}
	#endif



	MPI_Gatherv( data_local, MPI_WORK_DEVICE_TIMING_SIZE[MPI_MY_RANK], MPI_DOUBLE, // Where I store my stuff
		     data      , MPI_WORK_DEVICE_TIMING_SIZE             , MPI_WORK_DEVICE_TIMING_OFFSET, MPI_DOUBLE, // Distribution of work
			0, *MPI_WORK_GANG); // Who is recieving the data
}

inline void VECSEL::get_maxwell_MPI_balance_counters(double *timer, double *sdev)
{
	#ifdef MPI_BALANCE_WORKLOAD
	*timer = MPI_load->get_real_time_mean();
	*sdev = MPI_load->get_real_time_deviation();
	#else
	*timer = 0.0;
	*sdev = 0.0;
	#endif
}

inline void VECSEL::set_transverse_QW_pump_profile_SuperGaussian(int degree, double newFWHM)
{
	VECSEL_initial_pump_profile_SG_degree = degree;
	VECSEL_initial_pump_profile_SG_FWHM = newFWHM;
	VECSEL_initial_pump_profile_SG_sigma = newFWHM/(2.0*sqrt(2.0)*std::pow(log(2.0),1.0/(2.0*degree)));

	// QWs are already in the system, apply gaussian shape function to ALL QWs.
	if (quick_index_device.size()>0)
	{
		for(int i = 0; i < quick_index_device.size(); i++)
		{
			int indx = quick_index_device[i]; // Index of device

			double qw_y_pos = modules[indx].getDevice()->getTransversePosition();

			double b = std::pow(sqrt(2.0)*VECSEL_initial_pump_profile_SG_sigma,2.0*VECSEL_initial_pump_profile_SG_degree);
			double dx = std::pow(qw_y_pos,2.0*VECSEL_initial_pump_profile_SG_degree);
			double pump_profile_scale = exp(-dx/b);

			modules[indx].getDevice()->setTransverseBackgroundDensityScale(pump_profile_scale);
		}
	}
	// QWs are already in the system, apply gaussian shape function to ALL QWs.
	if (quick_index_twoArmDevice.size()>0)
	{
		for(int i = 0; i < quick_index_twoArmDevice.size(); i++)
		{
			int indx = quick_index_twoArmDevice[i]; // Index of device

			double qw_y_pos = modules[indx].getTwoArmDevice()->getTransversePosition();

			double b = std::pow(sqrt(2.0)*VECSEL_initial_pump_profile_SG_sigma,2.0*VECSEL_initial_pump_profile_SG_degree);
			double dx = std::pow(qw_y_pos,2.0*VECSEL_initial_pump_profile_SG_degree);
			double pump_profile_scale = exp(-dx/b);

			modules[indx].getTwoArmDevice()->setTransverseBackgroundDensityScale(pump_profile_scale);
		}
	}
}

inline void VECSEL::set_transverse_QW_temp_profile_SuperGaussian(int degree, double newFWHM)
{
	VECSEL_initial_temp_profile_SG_degree = degree;
	VECSEL_initial_temp_profile_SG_FWHM = newFWHM;
	VECSEL_initial_temp_profile_SG_sigma = newFWHM/(2.0*sqrt(2.0)*std::pow(log(2.0),1.0/(2.0*degree)));

	// Transfer to stored array
	if (VECSEL_initial_temp_profile_T == NULL)
	{
		VECSEL_initial_temp_profile_T = new double[VECSEL_transverse_points_number];
	}


	for(int i = 0; i < VECSEL_transverse_points_number; i++)
	{

		double b = std::pow(sqrt(2.0)*VECSEL_initial_temp_profile_SG_sigma,2.0*VECSEL_initial_temp_profile_SG_degree);
		double dx = std::pow(VECSEL_transverse_points_y[i],2.0*VECSEL_initial_temp_profile_SG_degree);
		VECSEL_initial_temp_profile_T[i] = exp(-dx/b);
	}
	

	// QWs are already in the system, apply gaussian shape function to ALL QWs.
	if (quick_index_device.size()>0||quick_index_twoArmDevice.size()>0)
	{
		// Send QW temperature scale to QWs
		set_transverse_QW_temp_profile(VECSEL_initial_temp_profile_T);
	}
}

inline void VECSEL::set_transverse_QW_temp_profile_file(int size, const std::string & fileName,double scale)
{
	double temp_array_z[size];
	double temp_array_T[size];

	// Load from file
	FILE *fid = fopen(fileName.c_str(),"r+");
	if (fid != NULL)
	{
		for(int i =0 ; i < size; i++)
		{
			int id = fscanf(fid,"%le %le",&(temp_array_z[i]),&(temp_array_T[i]));
			 temp_array_z[i]/=scale; //scale the transverse points to align with pump scaling
		}
		fclose(fid);
	} else {
		cout << "VECSEL::set_transverse_QW_temp_profile_file() Error: failed to read file: " << fileName.c_str() << endl;
		exit(-1);
	}

	if (VECSEL_transverse_points_y[0] < temp_array_z[0])
	{
		cout << "VECSEL::set_transverse_QW_temp_profile_file() Lower transverse temperature outside of range in file.." << endl;
		cout << "min_y grid = " << VECSEL_transverse_points_y[0]/um << " [um]" << endl;
		cout << "min_y file = " << temp_array_z[0]/um << " [um]" << endl;
		exit(-1);
	} else if (VECSEL_transverse_points_y[VECSEL_transverse_points_number-1] > temp_array_z[size-1])
	{
		cout << "VECSEL::set_transverse_QW_temp_profile_file() Upper transverse temperature outside of range in file.." << endl;
		cout << "max_y grid = " << VECSEL_transverse_points_y[VECSEL_transverse_points_number-1] << " [um]" << endl;
		cout << "max_y file = " << temp_array_z[VECSEL_transverse_points_number-1] << " [um]" << endl;
		exit(-1);
	}

	if (VECSEL_initial_temp_profile_T == NULL)
	{
		VECSEL_initial_temp_profile_T = new double[VECSEL_transverse_points_number];
	}

	// Interpolate to transverse grid
	for(int i = 0; i < VECSEL_transverse_points_number; i++)
	{
		double target = VECSEL_transverse_points_y[i];
		for(int j = 1; j < size; j++)
		{
			if ((temp_array_z[j-1] < target)&&(target <= temp_array_z[j]))
			{
				VECSEL_initial_temp_profile_T[i] = temp_array_T[j-1] + (temp_array_T[j]-temp_array_T[j-1])*(target-temp_array_z[j-1])/(temp_array_z[j]-temp_array_z[j-1]);
				break;
			}
		}
	}
	
	if (quick_index_device.size()>0 || quick_index_twoArmDevice.size()>0)
	{
		// Send QW temperature scale to QWs
		set_transverse_QW_temp_profile(VECSEL_initial_temp_profile_T);
		
	}
}


inline void VECSEL::set_transverse_QW_temp_profile(double *newScale)
{
	// Set Temperature scale to QWs
	for(int i = 0; i < quick_index_device.size(); i++)
	{
		int indx = quick_index_device[i]; // Index of device
		double qw_y_pos = modules[indx].getDevice()->getTransversePosition();

		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			if (fabs(qw_y_pos - VECSEL_transverse_points_y[j]) < 1e-12)
			{
				modules[indx].getDevice()->setTransverseBackgroundTemperatureScale(newScale[j]);
				break;
			}
		}
	}
	for(int i = 0; i < quick_index_twoArmDevice.size(); i++)
	{
		int indx = quick_index_twoArmDevice[i]; // Index of device
		double qw_y_pos = modules[indx].getTwoArmDevice()->getTransversePosition();

		for(int j = 0; j < VECSEL_transverse_points_number; j++)
		{
			if (fabs(qw_y_pos - VECSEL_transverse_points_y[j]) < 1e-12)
			{
				modules[indx].getTwoArmDevice()->setTransverseBackgroundTemperatureScale(newScale[j]);
				break;
			}
		}
	}

}


inline void VECSEL::set_transverse_dimensions(int numP, double R_max, double bg_ratio)
{
	if (R_max <= 0.0)
	{
		cout << "VECSEL():: Maximal transverse point: R_max <= 0" << endl;
		cout << "trying to use R_max = " << R_max << endl;
		exit(-1);
	}
	if (numP <= 0)
	{
		cout << "VECSEL():: Number of transverse points should be >= 1" << endl;
		cout << "trying to use VECSEL_transverse_points_number = " << VECSEL_transverse_points_number << endl;
		exit(-1);
	}

	VECSEL_transverse_points_number = numP;
	VECSEL_transverse_points_R_max = R_max;
	VECSEL_transverse_points_boundary_guard_ratio = bg_ratio;
	if (VECSEL_transverse_points_y==NULL)
	{
		VECSEL_transverse_points_y = new double[VECSEL_transverse_points_number];
		cavity_trans_E_pl = new std::complex<double>[VECSEL_transverse_points_number];
		cavity_trans_E_mi = new std::complex<double>[VECSEL_transverse_points_number];
		cavity_trans_E_pl_k = new std::complex<double>[VECSEL_transverse_points_number];
		cavity_trans_E_mi_k = new std::complex<double>[VECSEL_transverse_points_number];
		cavity_trans_MacPol = new std::complex<double>[VECSEL_transverse_points_number];

	} else {
		delete [] VECSEL_transverse_points_y;
		VECSEL_transverse_points_y = new double[VECSEL_transverse_points_number];

		delete [] cavity_trans_E_pl;
		delete [] cavity_trans_E_mi;
		delete [] cavity_trans_E_pl_k;
		delete [] cavity_trans_E_mi_k;
		delete [] cavity_trans_MacPol;

		cavity_trans_E_pl = new std::complex<double>[VECSEL_transverse_points_number];
		cavity_trans_E_mi = new std::complex<double>[VECSEL_transverse_points_number];
		cavity_trans_E_pl_k = new std::complex<double>[VECSEL_transverse_points_number];
		cavity_trans_E_mi_k = new std::complex<double>[VECSEL_transverse_points_number];
		cavity_trans_MacPol = new std::complex<double>[VECSEL_transverse_points_number];
	}
	
	if (VECSEL_transverse_points_number == 1)
	{
		// With a single point we only simulate r=0
		VECSEL_transverse_points_y[0] = 0.0;
	} else if (VECSEL_transverse_points_number > 1) {
		// only allow even number of points
		// i = N/2 will give r = 0 (center point)
		if (VECSEL_transverse_points_number%2 != 0)
		{
			cout << "VECSEL():: Only use even number of gridpoints.." << endl;
			cout << "trying to use VECSEL_transverse_points_number = " << VECSEL_transverse_points_number << endl;
			exit(-1);
		}

		double dx = R_max/((double)VECSEL_transverse_points_number);
		for(int i = 0; i < VECSEL_transverse_points_number; i++)
		{
			VECSEL_transverse_points_y[i] = dx*( (double)i - (double)VECSEL_transverse_points_number/2.0);
		}
	} else {
		cout  << "VECSEL::set_transverse_dimensions() cannot have # of transverse points <= 0 .." << endl;
		exit(-1);
	}
}
inline int VECSEL::getTransversePoints()
{
	return VECSEL_transverse_points_number;
}

inline double VECSEL::getTransverseRmax()
{
	return VECSEL_transverse_points_R_max;
}

inline double VECSEL::getTransverseBoundaryGuardRatio()
{
	return VECSEL_transverse_points_boundary_guard_ratio;
}


inline int VECSEL::getSGPumpDegree()
{
	return VECSEL_initial_pump_profile_SG_degree;
}
inline double VECSEL::getSGPumpFWHM()
{
	return VECSEL_initial_pump_profile_SG_FWHM;
}
inline int VECSEL::getSGTempDegree()
{
	return VECSEL_initial_temp_profile_SG_degree;
}
inline double VECSEL::getSGTempFWHM()
{
	return VECSEL_initial_temp_profile_SG_FWHM;
}

#endif






