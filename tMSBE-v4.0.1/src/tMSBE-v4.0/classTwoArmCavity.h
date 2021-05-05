
	

#ifndef __TWOARMCAVITY_H_INCLUDED__
#define __TWOARMCAVITY_H_INCLUDED__

#include <string>
#include <complex>
#include <fstream>
#include "classCyclic.h"
#include "classBPM.h"
#include "setup_simulation_variables.h"
#include "constantsAndMiscUnits.h"
#include <stdlib.h>
#include <iostream>
using namespace std;

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



//! The CAVITY class is placed inside the 'Module()' container between other objects.
/*!
	The material TwoArmCavity class stores the electromagnetic (EM) fields inside each material layer
	for a non-normal incience cavity with just a single angle.
the number of points that are depend on the given timestep.

	In order to iterate the cavity fields one has to call: updateStorate() \n
	There are various methods here that can be used to recover information about the EM field
	or set the EM field to a desired value

	Standalone useage: \n
	TwoArmCavity *a = new TwoArmCavity(refractive_index, cavity_width, position_1, position_2, dt); // Create object \n
	a->initializeZero(dt, NUM_TRANS, *transv_y); \n
	// .. \n
	// For each timestep \n
	a->setEfp(array1); \n
	a->setEfm(array1); \n
	a->setEbp(array1); \n
	a->setEbm(array1); \n
	a->updateStorage(); // Prepare for next timestep, clear up single space \n
	\n

	\sa Module \n
	2020: Module Created	\n
	Sam McLaren
*/

class TwoArmCavity
{
	public:
		//! A constructor
		TwoArmCavity()
		{
			setName("TACAV");
			setToFileOutputKey("out_");
			setRefInd(0.0);
			setRefInd_extraAxis(0.0);
			setRefInd_im(0.0);
			setCosTh(1.0, 1.0);
			setWidth(0.0);
			setPosition0(0.0);
			setPosition1(0.0);
			setNumberOfTimesteps(0);
			setNumberOfTimesteps_extraAxis(0);
			time_delay = -1;
			time_delay_ex = -1;
			cavity_dt = -1;
			electric_field_fp = NULL;
			electric_field_fm = NULL;
			electric_field_bp = NULL;
			electric_field_bm = NULL;
			transfer_matrix_a11 = 0.0;
			transfer_matrix_a12 = 0.0;
			transfer_matrix_a21 = 0.0;
			transfer_matrix_a22 = 0.0;
			transfer_matrix_a11_ex = 0.0;
			transfer_matrix_a12_ex = 0.0;
			transfer_matrix_a21_ex = 0.0;
			transfer_matrix_a22_ex = 0.0;

			transfer_matrix_b = 0.0;
			transfer_matrix_MacPol_fp = NULL;
			transfer_matrix_MacPol_fm = NULL;
			transfer_matrix_MacPol_bp = NULL;
			transfer_matrix_MacPol_bm = NULL;

			cavity_transverse_points_num = 0;
			cavity_transverse_points_y = NULL;
			temp_transverse_array_fp = NULL;
			temp_transverse_array_fm = NULL;
			temp_transverse_array_bp = NULL;
			temp_transverse_array_bm = NULL;
			output_E_real = NULL;
			output_E_imag = NULL;
			output_E_fp_real = NULL;
			output_E_fp_imag = NULL;
			output_E_fm_real = NULL;
			output_E_fm_imag = NULL;
			output_E_bp_real = NULL;
			output_E_bp_imag = NULL;
			output_E_bm_real = NULL;
			output_E_bm_imag = NULL;

			freeSpace_propagator = NULL;

		}
		//! A destructor
		~TwoArmCavity()
		{
			if (electric_field_fp != NULL)
			{
				delete  electric_field_fp;
				delete  electric_field_fm;
				delete  electric_field_bp;
				delete  electric_field_bm;
				delete [] cavity_transverse_points_y;
				delete [] transfer_matrix_MacPol_fp;
				delete [] transfer_matrix_MacPol_fm;
				delete [] transfer_matrix_MacPol_bp;
				delete [] transfer_matrix_MacPol_bm;

				delete [] output_E_real;
				delete [] output_E_imag;
				delete [] output_E_fp_real;
				delete [] output_E_fp_imag;
				delete [] output_E_fm_real;
				delete [] output_E_fm_imag;
				delete [] output_E_bp_real;
				delete [] output_E_bp_imag;
				delete [] output_E_bm_real;
				delete [] output_E_bm_imag;

				delete [] temp_transverse_array_fp;
				delete [] temp_transverse_array_fm;
				delete [] temp_transverse_array_bp;
				delete [] temp_transverse_array_bm;
			}

			if (freeSpace_propagator != NULL)
			{
				delete freeSpace_propagator;
			}
		}
		
		//! A copy constructor
		TwoArmCavity(const TwoArmCavity &obj)
		{
			cavity_name = obj.cavity_name;	
			cavity_output_key = obj.cavity_output_key;
			cavity_refractive_index = obj.cavity_refractive_index;
			cavity_refractive_index_ex = obj.cavity_refractive_index_ex;
			cavity_refractive_index_im = obj.cavity_refractive_index_im;
			cavity_left_cos_th = obj.cavity_left_cos_th;
			cavity_right_cos_th = obj.cavity_right_cos_th;
			cavity_width = obj.cavity_width;
			cavity_position_x0 = obj.cavity_position_x0;
			cavity_position_x1 = obj.cavity_position_x1;
			numberOfTimesteps = obj.numberOfTimesteps;
			numberOfTimesteps_extraAxis = obj.numberOfTimesteps_extraAxis;
			time_delay = obj.time_delay;
			time_delay_ex = obj.time_delay_ex;
			cavity_dt = obj.cavity_dt;

			
			cavity_transverse_points_num = obj.cavity_transverse_points_num;
			cavity_transverse_points_y = NULL;

			temp_transverse_array_fp = NULL;
			temp_transverse_array_fm = NULL;
			temp_transverse_array_bp = NULL;
			temp_transverse_array_bm = NULL;
			
			electric_field_fp = NULL;
			electric_field_fm = NULL;
			electric_field_bp = NULL;
			electric_field_bm = NULL;

			freeSpace_propagator = NULL;

			if (obj.electric_field_fp != NULL)
			{
				cavity_transverse_points_y = new double[cavity_transverse_points_num];
				
				electric_field_fp = new Cyclic<std::complex<double> >(*obj.electric_field_fp);
				electric_field_fm = new Cyclic<std::complex<double> >(*obj.electric_field_fm);
				electric_field_bp = new Cyclic<std::complex<double> >(*obj.electric_field_fp);
				electric_field_bm = new Cyclic<std::complex<double> >(*obj.electric_field_fm);
				
				for(int i = 0; i < cavity_transverse_points_num; i++)
				{
					cavity_transverse_points_y[i] = obj.cavity_transverse_points_y[i];
				}
				transfer_matrix_MacPol_fp = new std::complex<double>[cavity_transverse_points_num];
				transfer_matrix_MacPol_fm = new std::complex<double>[cavity_transverse_points_num];
				transfer_matrix_MacPol_bp = new std::complex<double>[cavity_transverse_points_num];
				transfer_matrix_MacPol_bm = new std::complex<double>[cavity_transverse_points_num];
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
				
			}
			if (obj.freeSpace_propagator != NULL)
			{
				freeSpace_propagator = new BPM(*obj.freeSpace_propagator);
			}
		}

		//! A constructor
		TwoArmCavity(double refInd, double width, double x0, double x1, double dt)
		{
			setName("TACAV");
			setToFileOutputKey("out__");
			setRefInd(refInd);
			setRefInd_extraAxis(refInd);
			setRefInd_im(0.0);
			setCosTh(1.0,1.0);
			setWidth(width);
			setPosition0(x0);
			setPosition1(x1);
			setNumberOfTimesteps(0);
			setNumberOfTimesteps_extraAxis(0);
			time_delay = -1;
			time_delay_ex = -1;
			cavity_dt = dt;

			transfer_matrix_a11 = 0.0;
			transfer_matrix_a12 = 0.0;
			transfer_matrix_a21 = 0.0;
			transfer_matrix_a22 = 0.0;
			transfer_matrix_a11_ex = 0.0;
			transfer_matrix_a12_ex = 0.0;
			transfer_matrix_a21_ex = 0.0;
			transfer_matrix_a22_ex = 0.0;

			transfer_matrix_b = 0.0;
			transfer_matrix_MacPol_fp = NULL;
			transfer_matrix_MacPol_fm = NULL;
			transfer_matrix_MacPol_bm = NULL;
			transfer_matrix_MacPol_bm = NULL;

			electric_field_fp = NULL;
			electric_field_fm = NULL;
			electric_field_bp = NULL;
			electric_field_bm = NULL;

			cavity_transverse_points_num = -1;
			cavity_transverse_points_y = NULL;
	
			temp_transverse_array_fp = NULL;
			temp_transverse_array_fm = NULL;
			temp_transverse_array_bp = NULL;
			temp_transverse_array_bm = NULL;
		}

		//! Print an overview of the cavity data to screen
		void Print(void)
		{
			cout << "Print two arm cavity:" << endl;
			cout << " -> name    = " << getName() << endl;
			cout << " -> n       = " << getRefInd() << " +i " << getRefInd_im() << endl;
			cout << " -> n_ex       = " << getRefInd_extraAxis() << endl;
			cout << " -> width   = " << getWidth()/um << " [um]" << endl;
			cout << " -> x0      = " << getPosition0() << endl;
			cout << " -> x1      = " << getPosition1() << endl;
			int dL = getDelay();
			cout << " -> delay   = " << dL << endl;
			dL = getDelay_extraAxis();
			cout << " -> delay extra axis  = " << dL << endl;
			double left, right; getCosTh(&left, &right);
			cout << " ->cos(th)  = " << left << ", " << right << endl;
			cout << " -> # trans = " << cavity_transverse_points_num << endl;
			cout << " -> r [um]  = ";
			for(int i = 0; i < cavity_transverse_points_num; i++)
			{
				cout << cavity_transverse_points_y[i]/um << ", ";
			}
			cout << endl;
		}


		//! Set the name of the cavity
		void setName(const std::string &);

		//! Return the cavity name
		std::string getName(void) const;

		//! Set the file output name
		void setToFileOutputKey(const std::string &);

		//! Return the file output name
		std::string getToFileOutputKey(void) const;

		//! Return the REAL refractive index
		double getRefInd(void) const;

		//! Set the REAL refractive index (plus arm in dual prop)
		void   setRefInd(double);

		//! Set the IMAG refractive index (plus arm in dual prop)
		void   setRefInd_im(double);

		//! Return the REAL extraordinary refractive index (minus arm in dual prop)
		double getRefInd_extraAxis(void) const;

		//! Return the IMAG refractive index
		double getRefInd_im(void) const;

		//! Set the extraordinary axis refractive index
		void   setRefInd_extraAxis(double);

		//! Return the cavity width
		double getWidth(void) const;

		//! Set the cavity width
		void   setWidth(double);

		//! Return the z-position of the left boundary
		double getPosition0(void) const;

		//! Set the z-position of the left boundary
		void   setPosition0(double);

		//! Return the z-position of the right boundary
		double getPosition1(void) const;
	
		//! Set the z-position of the right boundary
		void   setPosition1(double);

		//! Return the number cos(theta) for this cavity
		void   getCosTh(double*, double*) const;

		//! Set the number cos(theta) for this cavity
		void   setCosTh(double, double);

		//! Return number of timesteps used in cavity
		int    getNumberOfTimesteps(void) const;

		//! Set number of timesteps used in cavity
		void   setNumberOfTimesteps(int);

		//! Return number of timesteps used in cavity along extraordinary axis
		int    getNumberOfTimesteps_extraAxis(void) const;

		//! Set number of timesteps used in cavity along extraordinary axis
		void   setNumberOfTimesteps_extraAxis(int);

		//! Return cavity resolution in terms of timesteps along extraordinary axis
		/*! There are #timesteps +2 stored in each container 
		    with the extra two helping with boundary effects
		*/
		int    getDelay_extraAxis(void) const;

		//! Return cavity resolution in terms of timesteps
		/*! There are #timesteps +2 stored in each container 
		    with the extra two helping with boundary effects
		*/
		int    getDelay(void) const;

		//! Set cavity resolution in terms of timesteps along extraordinary axis
		/*! There are #timesteps +2 stored in each container 
		    with the extra two helping with boundary effects
		*/
		void   setDelay_extraAxis(int);

		//! Set cavity resolution in terms of timesteps
		/*! There are #timesteps +2 stored in each container 
		    with the extra two helping with boundary effects
		*/
		void   setDelay(int);

		//! Return timestep 'dt' used in simulation
		int    getCavityDT(void) const;

		//! Set timestep 'dt' used in simulation
		void   setCavityDT(double);
		
		//! Set first value of E forward pluss
		void setEfp(std::complex<double>*);

		//! Set first value of E forward minus
		void setEfm(std::complex<double>*);

		//! Set first value of E backward pluss
		void setEbp(std::complex<double>*);

		//! Set first value of E backward minus
		void setEbm(std::complex<double>*);
		
		//! Return pointer to memory of first value of E forward pluss
		/*! Overloaded function for setting memory*/
		std::complex<double> * setEfp(void);

		//! Return pointer to memory of first value of E forward minus
		/*! Overloaded function for setting memory*/
		std::complex<double> * setEfm(void);

		//! Return pointer to memory of first value of E backward pluss
		/*! Overloaded function for setting memory*/
		std::complex<double> * setEbp(void);

		//! Return pointer to memory of first value of E backward minus
		/*! Overloaded function for setting memory*/
		std::complex<double> * setEbm(void);

		//! Return first value of E forward pluss
		void getEfp(int, std::complex<double>*) const;

		//! Return first value of E forward minus
		void getEfm(int, std::complex<double>*) const;

		//! Return first value of E backward pluss
		void getEbp(int, std::complex<double>*) const;

		//! Return first value of E backward minus
		void getEbm(int, std::complex<double>*) const;
		
		//! Initialize all storage in this class
		/*! This function initializes the class using the simulation setup+timestep
		    It is assumed that the cavity has already been setup with required dimensions
		    \param DT timestep of simulation
		    \param Nx Number of transverse points
		    \param x Transverse points x in [-R_max/2, (1/2 - 1/N)R_max]
		    \param R_max Width of transverse grid
		    \param bg_ratio Ratio of domain used as boundary guard (SuperGaussian) for FFT-BPM propagator
		*/
		void initializeZero(double DT, int Nx, double*x, double R_max, double boundary_guard_rato);

		//! Save four E fields to binary files in folder 'save/'
		/*! Saves four E fields  to binary files named with 'cavity_name' and transverse dimension number
			If folder 'save/' does not exist in directory nothing can be saved
			If files already exist, they will be overwritten.
		 */

		void file_save_variables(int, int);
		//! Load four E fields from binary files in folder 'save/'
		/*! Load four E field files from files named with 'cavity_name' and transverse dimension number
			This function assumes that the files in 'save/' have same number of timesteps as cavity
			If files cannot be found, the program will stop 
 		 */
		void file_load_variables(int, int);

		//! Update storage of four E fields
		/*! This function iterates the storage structures such that a new timestep of E+/E- is ready
		    to be stored. The new storage space is NOT set to zero..\n
		    Note: Computes FFT BPM if needed
		 */
		void updateStorage_freeSpace();

		//! Set all components of four E fields to zero
		void clearFields();
		
		//! Write all simulation variables to binary file(s) at left wall
		/*! The specific output data can be changed in this function.
		    Complex numbers are saved in two files named with 're' and 'im'
		    If you change the output data, also change the functions: file_output_open and file_output_close
		 */
		void file_output_write_backWall(int);
		
		//! Write all simulation variables to binary file(s)
		/*! The specific output data can be changed in this function.
		    Complex numbers are saved in two files named with 're' and 'im'
		    If you change the output data, also change the functions: file_output_open and file_output_close
		 */
		void file_output_write(int);

		//! Open new binary output files with new names if required
		/*! Open new binary output files that are ready to record output
			The files are named according to the 'cavity_name' and given output names
			The new files will be created or overwrite old files with the same name
		        If you change the output data, also change the functions: file_output_write and file_output_close
		 */
		void file_output_open(int,int);

		//! Close output files
		/*! If you change the output data, also change the functions: file_output_write and file_output_close
		 */
		void file_output_close(int);
		
		//! Return E forward pluss evaluated on the front boundary
		void getEfp_front_wall(std::complex<double> *);
		
		//! Return E forward minus evaluated on the front boundary
		void getEfm_front_wall(std::complex<double> *);

		//! Return E backward pluss evaluated on the front boundary
		/*! Note: calls 'interpolateEbp_x0()' to evaluate on the boundary */
		void getEbp_front_wall(std::complex<double> *);

		//! Return E backward minus evaluated on the front boundary
		/*! Note: calls 'interpolateEbp_x0()' to evaluate on the boundary */
		void getEbm_front_wall(std::complex<double> *);

		//! Return E backward minus evaluated on the front boundary for extraodinary axis propgation
		/*! Note: calls 'interpolateEbp_x0()' to evaluate on the boundary */
		void getEbm_ex_front_wall(std::complex<double> *);

		//! Return E backward pluss evaluated on the back boundary
		void getEbp_back_wall(std::complex<double> *);
		
		//! Return E backward minus evaluated on the back boundary
		void getEbm_back_wall(std::complex<double> *);

		//! Return E forward pluss evaluated on the back boundary
		/*! Note: calls 'interpolateEfp_x1()' to evaluate on the boundary */
		void getEfp_back_wall(std::complex<double> *);

		//! Return E forward minus evaluated on the back boundary
		/*! Note: calls 'interpolateEfm_x0()' to evaluate on the boundary */
		void getEfm_back_wall(std::complex<double> *);

		//! Return E forward minus evaluated on the back boundary for extraodinary axis propagation
		/*! Note: calls 'interpolateEfm_x0()' to evaluate on the boundary */
		void getEfm_ex_back_wall(std::complex<double> *);
		
		//! Return E forward pluss evaluated on the front boundary at the previous timestep
		void getEfp_front_wall_tp1(std::complex<double> *);
		
		//! Return E forward minus evaluated on the front boundary at the previous timestep
		void getEfm_front_wall_tp1(std::complex<double> *);

		//! Return E backward pluss evaluated on the front boundary at the previous timestep
		/*! Note: calls 'interpolateEbp_x0()' to evaluate on the boundary at the previous timestep */
		void getEbp_front_wall_tp1(std::complex<double> *);

		//! Return E backward minus evaluated on the front boundary at the previous timestep
		/*! Note: calls 'interpolateEbm_x0()' to evaluate on the boundary at the previous timestep */
		void getEbm_front_wall_tp1(std::complex<double> *);

		//! Return E backward minus evaluated on the front boundary at the previous timestepi for extraodinary axis propagation
		/*! Note: calls 'interpolateEbm_x0()' to evaluate on the boundary at the previous timestep */
		void getEbm_ex_front_wall_tp1(std::complex<double> *);

		//! Return E backward pluss evaluated on the back boundary at the previous timestep
		void getEbp_back_wall_tp1(std::complex<double> *);
		
		//! Return E backward minus evaluated on the back boundary at the previous timestep
		void getEbm_back_wall_tp1(std::complex<double> *);

		//! Return E forward pluss evaluated on the back boundary at the previous timestep
		/*! Note: calls 'interpolateEfp_x0()' to evaluate on the boundary at the previous timestep */
		void getEfp_back_wall_tp1(std::complex<double> *);

		//! Return E forward minus evaluated on the back boundary at the previous timestep
		/*! Note: calls 'interpolateEfm_x0()' to evaluate on the boundary at the previous timestep */
		void getEfm_back_wall_tp1(std::complex<double> *);

		//! Return E forward minus evaluated on the back boundary at the previous timestep for extraordinary axis propagation
		/*! Note: calls 'interpolateEfm_x0()' to evaluate on the boundary at the previous timestep */
		void getEfm_ex_back_wall_tp1(std::complex<double> *);

		//! Return E all eight fields on the back boundary at the previous timestep
		/*! Note: fixed at 8. For future updates this may be changed */
		void evaluateEprop_back_wall(std::complex<double> *tmp);

		//! Return E all eight fields on the front boundary at the previous timestep
		/*! Note: fixed at 8. For future updates this may be changed */
		void evaluateEprop_front_wall(std::complex<double> *tmp);

		//! Return E all eight fields on the back boundary at the previous timestep for extraodinary axis included
		/*! Note: fixed at 8. For future updates this may be changed */
		void evaluateEprop_back_wall_ex(std::complex<double> *tmp);

		//! Return E all eight fields on the front boundary at the previous timestep for extraodinary axis included
		/*! Note: fixed at 8. For future updates this may be changed */
		void evaluateEprop_front_wall_ex(std::complex<double> *tmp);

		//! Return E backward pluss on the front boundary (possibly) interpolated between two timesteps
		void interpolateEbp_x0(std::complex<double> *);

		std::complex<double> * interpolateEbp_x0(void);
		
		//! Return E backward minus on the front boundary (possibly) interpolated between two timesteps
		void interpolateEbm_x0(std::complex<double> *);

		//! Return E backward minus on the front boundary (possibly) interpolated between two timesteps extraodinary axis
		void interpolateEbm_ex_x0(std::complex<double> *);

		std::complex<double> * interpolateEbm_x0(void);

		//! Return E forward pluss on the back boundary linearly interpolated between two timesteps
		void interpolateEfp_x1(std::complex<double> *);

		std::complex<double> * interpolateEfp_x1(void);

		//! Return E forward minus on the back boundary linearly interpolated between two timesteps
		void interpolateEfm_x1(std::complex<double> *);

		//! Return E forward minus on the back boundary linearly interpolated between two timesteps extraordinary axis
		void interpolateEfm_ex_x1(std::complex<double> *);

		std::complex<double> * interpolateEfm_x1(void);

		//! Return E backward pluss from the past timestep evaluated on the back boundary linearly interpolated between two timesteps
		void interpolateEbp_x0_tp1(std::complex<double> *);

		std::complex<double> * interpolateEbp_x0_tp1(void);

		//! Return E forward pluss from the past timestep evaluated on the right boundary linearly interpolated between two timesteps
		void interpolateEfp_x1_tp1(std::complex<double> *);

		std::complex<double> * interpolateEfp_x1_tp1(void);

		//! Return E backward minus from the past timestep evaluated on the back boundary linearly interpolated between two timesteps
		void interpolateEbm_x0_tp1(std::complex<double> *);

		//! Return E backward minus from the past timestep evaluated on the back boundary linearly interpolated between two timesteps extraordinary axis
		void interpolateEbm_ex_x0_tp1(std::complex<double> *);

		std::complex<double> * interpolateEbm_x0_tp1(void);

		//! Return E forward minus from the past timestep evaluated on the right boundary linearly interpolated between two timesteps
		void interpolateEfm_x1_tp1(std::complex<double> *);

		//! Return E forward minus from the past timestep evaluated on the right boundary linearly interpolated between two timesteps extraordinary axis
		void interpolateEfm_ex_x1_tp1(std::complex<double> *);

		std::complex<double> * interpolateEfm_x1_tp1(void);

		//! Return E forward pluss at a given 'z' position inside the materal layer. Uses linear interpolation.
		void interpolateEfp(double, std::complex<double> *);

		//! Return E forward minus at a given 'z' position inside the materal layer. Uses linear interpolation.
		void interpolateEfm(double, std::complex<double> *);

		//! Return E backward pluss at a given 'z' position inside the materal layer. Uses linear interpolation.
		void interpolateEbp(double, std::complex<double> *);

		//! Return E backward minus at a given 'z' position inside the materal layer. Uses linear interpolation.
		void interpolateEbm(double, std::complex<double> *);
		
		//! Return the element needed to compute transfer matrix method.
		/*! This method requires that the internal variables have been initialized using 'set_transfer_matrix()' and 'set_transfer_matrix_macPol()' for each macPol. The matrix elements a11, a12, a21, a22 and the macroscopic polarization is returned. */
		void get_transfer_matrix(std::complex<double> *a11, std::complex<double> *a12, std::complex<double> *a21, std::complex<double> *a22, std::complex<double> **b_fp, std::complex<double> **b_fm, std::complex<double> **b_bp, std::complex<double> **b_bm);
		
		//! Return the element needed to compute transfer matrix method. 
		/*! This method requires that the internal variables have been initialized using 'set_transfer_matrix_extraAxis()'. The matrix elements a11_ex, a12_ex, a21_ex, a22_ex */
		void get_transfer_matrix_extraAxis(std::complex<double> *, std::complex<double> *, std::complex<double> *, std::complex<double> *);

		//! Initialize the transfer matrix method variables for the given cavity
		/*! Computes the matrix elements a11, a12, a21, a22 and the macroscopic polarization constants. Includes any angular dependence. */
		void set_transfer_matrix(std::complex<double>, std::complex<double>, double, double);

		//! Initialize the transfer matrix method variables for the given cavity for birefringent crystal with extraordinary axis.
		/*! Computes the matrix elements a11_ex, a12_ex, a21_ex, a22_ex. */
		void set_transfer_matrix_extraAxis(std::complex<double>, std::complex<double>);

		//! Set the macroscopic polarization for forward pluss direction that will later be used in the transfer matrix method
		void set_transfer_matrix_macPol_fp(std::complex<double> *);

		//! Set the macroscopic polarization for forward minus direction that will later be used in the transfer matrix method
		void set_transfer_matrix_macPol_fm(std::complex<double> *);

		//! Set the macroscopic polarization for backward pluss direction that will later be used in the transfer matrix method
		void set_transfer_matrix_macPol_bp(std::complex<double> *);

		//! Set the macroscopic polarization for backward minus direction that will later be used in the transfer matrix method
		void set_transfer_matrix_macPol_bm(std::complex<double> *);

		//! Returns the transport phase factor for E forward in this cavity
                std::complex<double> get_transport_Ef_x1(void);

                //! Returns the transport phase factor for E backward in this cavity
                std::complex<double> get_transport_Eb_x0(void);

		//! Returns the transport phase factor for E forward in this cavity along extraordinary axis
                std::complex<double> get_transport_Ef_ex_x1(void);

                //! Returns the transport phase factor for E backward in this cavity along extraordinary axis
                std::complex<double> get_transport_Eb_ex_x0(void);
	
	private:
		//! Name of device
		std::string cavity_name;

		//! Start of output name of files
		std::string cavity_output_key;

		//! Refractive index of cavity (REAL) (plus arm in dual prop)
		double cavity_refractive_index;

		//! Refractive index of cavity (REAL) extraoindary axis (minus arm in dual prop)
		double cavity_refractive_index_ex;

		//! Refractive index of cavity (IMAG)
		double cavity_refractive_index_im;

		//! Length of Cavity [m]
		double cavity_width;

		//! Position of Cavity left boundary [m]
		double cavity_position_x0;

		//! Position of cavity right boundary [m]
		double cavity_position_x1;

		//! # timestesp the cavity delays the electric field when propagating
		int time_delay;

		//! # timesteps the cavity delays the electric field when propagating along extraordinary axis
		int time_delay_ex;

		//! Timestep from the external solver [s]
		double cavity_dt;

		//! The factor for the angular dependence given by the angle of incidence cos(th) on the left
		double cavity_left_cos_th;

		//! The factor for the angular dependence given by the angle of incidence cos(th) on the right
		double cavity_right_cos_th;

		//! # of transverse points in cavity
		int cavity_transverse_points_num;

		//! Containes the position of each transverse points starting [-R_max, R_max]
		double *cavity_transverse_points_y;
		
		//! Interpolation weight for left boundary
		double cavity_weight_x0;

		//! Interpolation weight for right boundary
		double cavity_weight_x1;
		
		//! Interpolation weight for left boundary extraordinary axis
		double cavity_weight_ex_x0;

		//! Interpolation weight for right boundary extraordinary axis
		double cavity_weight_ex_x1;

		//! Temporary array to avoid alot of heap-creation/destruction
		std::complex<double> *temp_transverse_array_fp;

		//! Temporary array to avoid alot of heap-creation/destruction
		std::complex<double> *temp_transverse_array_fm;

		//! Temporary array to avoid alot of heap-creation/destruction
		std::complex<double> *temp_transverse_array_bp;

		//! Temporary array to avoid alot of heap-creation/destruction
		std::complex<double> *temp_transverse_array_bm;

		//! Transport phase factor for E forward to transfer from front to back boundary (1)
		std::complex<double> Ef_transp_x0;

		//! Transport phase factor for E forward to transfer from back to front boundary
		std::complex<double> Ef_transp_x1;

		//! Transport phase factor for E backward to transfer from front to back boundary (1)
		std::complex<double> Eb_transp_x0;

		//! Transport phase factor for E backward to transfer from back to front boundary
		std::complex<double> Eb_transp_x1;

		//! Transport phase factor for E forward to transfer from back to front boundary extraordinary axis
		std::complex<double> Ef_ex_transp_x1;

		//! Transport phase factor for E backward to transfer from front to back boundary (1) extraordinary axis
		std::complex<double> Eb_ex_transp_x0;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a11;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a12;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a21;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a22;

		//! Transfer matrix element precomputed with refractive indices for extraordinary axis
		std::complex<double> transfer_matrix_a11_ex;

		//! Transfer matrix element precomputed with refractive indices for extraordinary axis
		std::complex<double> transfer_matrix_a12_ex;

		//! Transfer matrix element precomputed with refractive indices for extraordinary axis
		std::complex<double> transfer_matrix_a21_ex;

		//! Transfer matrix element precomputed with refractive indices for extraordinary axis
		std::complex<double> transfer_matrix_a22_ex;

		//! Transfer matrix  macroscopic polarization factors
		std::complex<double> transfer_matrix_b;

		//! Transfer matrix storage current for macroscopic polarization for forward pluss direction
		std::complex<double> *transfer_matrix_MacPol_fp;

		//! Transfer matrix storage current for macroscopic polarization for forward minus direction
		std::complex<double> *transfer_matrix_MacPol_fm;

		//! Transfer matrix storage current for macroscopic polarization for backward pluss direction
		std::complex<double> *transfer_matrix_MacPol_bp;

		//! Transfer matrix storage current for macroscopic polarization for backward minus direction
		std::complex<double> *transfer_matrix_MacPol_bm;
		
		//! Electric field E forward pluss propagating from front boundary to back boundary
		Cyclic<std::complex<double> > *electric_field_fp;

		//! Electric field E forward minus propagating from forward boundary to back boundary
		Cyclic<std::complex<double> > *electric_field_fm;
		
		//! Electric field E backward pluss propagating from back boundary to front boundary
		Cyclic<std::complex<double> > *electric_field_bp;

		//! Electric field E backward minus propagating from back boundary to front boundary
		Cyclic<std::complex<double> > *electric_field_bm;

		//! FFT BPM free space propagator
		BPM *freeSpace_propagator;
	
		//! Number of stored timesteps in E+ and E- (plus arm also ordinary axis)
		int numberOfTimesteps;
	
		//! Number of stored timesteps in E+ and E- (minus arm also extraordinary axis)
		int numberOfTimesteps_extraAxis;
		
		//! Output to file of real(E)
		std::ofstream *output_E_real;

		//! Output to file of imag(E)
		std::ofstream *output_E_imag;

		//! Output to file of real(E_fp)
		std::ofstream *output_E_fp_real;

		//! Output to file of imag(E_fp)
		std::ofstream *output_E_fp_imag;

		//! Output to file of real(E_fm)
		std::ofstream *output_E_fm_real;

		//! Output to file of imag(E_fm)
		std::ofstream *output_E_fm_imag;

		//! Output to file of real(E_bp)
		std::ofstream *output_E_bp_real;

		//! Output to file of imag(E_bp)
		std::ofstream *output_E_bp_imag;

		//! Output to file of real(E_bm)
		std::ofstream *output_E_bm_real;

		//! Output to file of imag(E_bm)
		std::ofstream *output_E_bm_imag;
};

inline void TwoArmCavity::setName(const std::string & newName)
{
	cavity_name = newName;
}

inline std::string TwoArmCavity::getName(void) const
{
	return cavity_name;
}

inline void TwoArmCavity::setCosTh(double new_cos_th_left, double new_cos_th_right)
{
	cavity_left_cos_th  = new_cos_th_left;
	cavity_right_cos_th = new_cos_th_right;
}

inline void TwoArmCavity::getCosTh(double *left, double *right) const
{
	*left  = cavity_left_cos_th;
	*right = cavity_right_cos_th;
}

inline void TwoArmCavity::setToFileOutputKey(const std::string & newName)
{
	cavity_output_key = newName;
}
inline std::string TwoArmCavity::getToFileOutputKey(void) const
{
	return cavity_output_key;
}

inline void TwoArmCavity::setRefInd(double newInd)
{
	cavity_refractive_index = newInd;
}

inline void TwoArmCavity::setRefInd_extraAxis(double newInd)
{
	cavity_refractive_index_ex = newInd;
}

inline void TwoArmCavity::setRefInd_im(double newInd)
{
	cavity_refractive_index_im = newInd;
}

inline double TwoArmCavity::getRefInd(void) const
{
	return cavity_refractive_index;
}

inline double TwoArmCavity::getRefInd_extraAxis(void) const
{
	return cavity_refractive_index_ex;
}

inline double TwoArmCavity::getRefInd_im(void) const
{
	return cavity_refractive_index_im;
}

inline void TwoArmCavity::setWidth(double newWidth)
{
	cavity_width = newWidth;
}

inline double TwoArmCavity::getWidth(void) const
{
	return cavity_width;
}

inline void TwoArmCavity::setPosition0(double newX0)
{
	cavity_position_x0 = newX0;
}

inline double TwoArmCavity::getPosition0(void) const
{
	return cavity_position_x0;
}

inline void TwoArmCavity::setPosition1(double newX1)
{
	cavity_position_x1 = newX1;
}

inline double TwoArmCavity::getPosition1(void) const
{
	return cavity_position_x1;
}

inline void TwoArmCavity::setNumberOfTimesteps(int nt)
{
	numberOfTimesteps = nt;
}

inline int TwoArmCavity::getNumberOfTimesteps_extraAxis(void) const
{
	return numberOfTimesteps_extraAxis;
}

inline void TwoArmCavity::setNumberOfTimesteps_extraAxis(int nt)
{
	numberOfTimesteps_extraAxis = nt;
}

inline int TwoArmCavity::getNumberOfTimesteps(void) const
{
	return numberOfTimesteps;
}

inline int TwoArmCavity::getDelay_extraAxis(void) const
{
	return time_delay_ex;
}

inline int TwoArmCavity::getDelay(void) const
{
	return time_delay;
}

inline void TwoArmCavity::setDelay(int delay)
{
	time_delay = delay;
}

inline void TwoArmCavity::setDelay_extraAxis(int delay)
{
	time_delay_ex = delay;
}

inline int TwoArmCavity::getCavityDT(void) const
{
	return cavity_dt;
}

inline void TwoArmCavity::setCavityDT(double newDT)
{
	cavity_dt = newDT;
}

inline void TwoArmCavity::setEfp(std::complex<double> *target)
{
	electric_field_fp->setFirstContainer(target);
}

inline void TwoArmCavity::setEfm(std::complex<double> *target)
{
	electric_field_fm->setFirstContainer(target);
}

inline void TwoArmCavity::setEbp(std::complex<double> *target)
{
	electric_field_bp->setFirstContainer(target);
}

inline void TwoArmCavity::setEbm(std::complex<double> *target)
{
	electric_field_bm->setFirstContainer(target);
}

inline std::complex<double> * TwoArmCavity::setEfp()
{
	return electric_field_fp->getFirstContainer();
}

inline std::complex<double> * TwoArmCavity::setEfm()
{
	return electric_field_fm->getFirstContainer();
}

inline std::complex<double> * TwoArmCavity::setEbp()
{
	return electric_field_bp->getFirstContainer();
}

inline std::complex<double> * TwoArmCavity::setEbm()
{
	return electric_field_bm->getFirstContainer();
}

inline void TwoArmCavity::getEfp(int k, std::complex<double> *tmp) const
{
	memcpy(tmp, electric_field_fp->getContainerNr(k), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void TwoArmCavity::getEfm(int k, std::complex<double> *tmp) const
{
	memcpy(tmp, electric_field_fm->getContainerNr(k), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void TwoArmCavity::getEbp(int k, std::complex<double> *tmp) const
{
	memcpy(tmp, electric_field_bp->getContainerNr(k), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void TwoArmCavity::getEbm(int k, std::complex<double> *tmp) const
{
	memcpy(tmp, electric_field_bm->getContainerNr(k), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void TwoArmCavity::getEfp_front_wall(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_fp->getFirstContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void TwoArmCavity::getEbp_back_wall(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_bp->getFirstContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void TwoArmCavity::getEfm_front_wall(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_fm->getFirstContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void TwoArmCavity::getEbm_back_wall(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_bm->getFirstContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void TwoArmCavity::getEfp_back_wall(std::complex<double> *tmp)
{
	interpolateEfp_x1(tmp);
	cblas_zscal(cavity_transverse_points_num, &Ef_transp_x1, tmp , 1);
}

inline void TwoArmCavity::getEbp_front_wall(std::complex<double> *tmp)
{
	interpolateEbp_x0(tmp);
	cblas_zscal(cavity_transverse_points_num, &Eb_transp_x0, tmp , 1);
}

inline void TwoArmCavity::getEfm_back_wall(std::complex<double> *tmp)
{
	interpolateEfm_x1(tmp);
	cblas_zscal(cavity_transverse_points_num, &Ef_transp_x1, tmp , 1);
}

inline void TwoArmCavity::getEbm_front_wall(std::complex<double> *tmp)
{
	interpolateEbm_x0(tmp);
	cblas_zscal(cavity_transverse_points_num, &Eb_transp_x0, tmp , 1);
}

inline void TwoArmCavity::getEfm_ex_back_wall(std::complex<double> *tmp)
{
	interpolateEfm_x1(tmp);
	cblas_zscal(cavity_transverse_points_num, &Ef_ex_transp_x1, tmp , 1);
}

inline void TwoArmCavity::getEbm_ex_front_wall(std::complex<double> *tmp)
{
	interpolateEbm_x0(tmp);
	cblas_zscal(cavity_transverse_points_num, &Eb_ex_transp_x0, tmp , 1);
}


inline void TwoArmCavity::getEfp_front_wall_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_fp->get2ndContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void TwoArmCavity::getEbp_back_wall_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_bp->get2ndContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void TwoArmCavity::getEfm_front_wall_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_fm->get2ndContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void TwoArmCavity::getEbm_back_wall_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_bm->get2ndContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void TwoArmCavity::getEfp_back_wall_tp1(std::complex<double> *tmp)
{
	interpolateEfp_x1_tp1(tmp);
	cblas_zscal(cavity_transverse_points_num, &Ef_transp_x1, tmp , 1);
}

inline void TwoArmCavity::getEbp_front_wall_tp1(std::complex<double> *tmp)
{
	interpolateEbp_x0_tp1(tmp);
	cblas_zscal(cavity_transverse_points_num, &Eb_transp_x0, tmp , 1);
}

inline void TwoArmCavity::getEfm_back_wall_tp1(std::complex<double> *tmp)
{
	interpolateEfm_x1_tp1(tmp);
	cblas_zscal(cavity_transverse_points_num, &Ef_transp_x1, tmp , 1);
}

inline void TwoArmCavity::getEbm_front_wall_tp1(std::complex<double> *tmp)
{
	interpolateEbm_x0_tp1(tmp);
	cblas_zscal(cavity_transverse_points_num, &Eb_transp_x0, tmp , 1);
}

inline void TwoArmCavity::getEfm_ex_back_wall_tp1(std::complex<double> *tmp)
{
	interpolateEfm_x1_tp1(tmp);
	cblas_zscal(cavity_transverse_points_num, &Ef_ex_transp_x1, tmp , 1);
}

inline void TwoArmCavity::getEbm_ex_front_wall_tp1(std::complex<double> *tmp)
{
	interpolateEbm_x0_tp1(tmp);
	cblas_zscal(cavity_transverse_points_num, &Eb_ex_transp_x0, tmp , 1);
}

inline void TwoArmCavity::evaluateEprop_back_wall(std::complex<double> *tmp)
{
	std::complex<double> *tmp_Eprop_forward;
        
	tmp_Eprop_forward=interpolateEfp_x1();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[0], 8);
	cblas_zscal(cavity_transverse_points_num, &Ef_transp_x1, &tmp[0], 8);
        tmp_Eprop_forward=interpolateEfp_x1_tp1();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[1], 8);
	cblas_zscal(cavity_transverse_points_num, &Ef_transp_x1, &tmp[1], 8);
        tmp_Eprop_forward=interpolateEfm_x1();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[2], 8);
	cblas_zscal(cavity_transverse_points_num, &Ef_transp_x1, &tmp[2], 8);
        tmp_Eprop_forward=interpolateEfm_x1_tp1();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[3], 8);
	cblas_zscal(cavity_transverse_points_num, &Ef_transp_x1, &tmp[3], 8);
	
	std::complex<double> tmp_Eprop_backward[cavity_transverse_points_num];	
	getEbp_back_wall(tmp_Eprop_backward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[4], 8);
	getEbp_back_wall_tp1(tmp_Eprop_backward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[5], 8);
	getEbm_back_wall(tmp_Eprop_backward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[6], 8);
	getEbm_back_wall_tp1(tmp_Eprop_backward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[7], 8);
}

inline void TwoArmCavity::evaluateEprop_front_wall(std::complex<double> *tmp)
{
	std::complex<double> dummy=0.0;
	std::complex<double> tmp_Eprop_forward[cavity_transverse_points_num];	
	getEfp_front_wall(tmp_Eprop_forward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[0], 8);
	getEfp_front_wall_tp1(tmp_Eprop_forward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[1], 8);
	getEfm_front_wall(tmp_Eprop_forward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[2], 8);
	getEfm_front_wall_tp1(tmp_Eprop_forward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[3], 8);
	
	std::complex<double> *tmp_Eprop_backward; 

	tmp_Eprop_backward=interpolateEbp_x0();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[4], 8);
	cblas_zscal(cavity_transverse_points_num, &Eb_transp_x0, &tmp[4], 8);
        tmp_Eprop_backward=interpolateEbp_x0_tp1();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[5], 8);
	cblas_zscal(cavity_transverse_points_num, &Eb_transp_x0, &tmp[5], 8);
        tmp_Eprop_backward=interpolateEbm_x0();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[6], 8);
	cblas_zscal(cavity_transverse_points_num, &Eb_transp_x0, &tmp[6], 8);
        tmp_Eprop_backward=interpolateEbm_x0_tp1();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[7], 8);
	cblas_zscal(cavity_transverse_points_num, &Eb_transp_x0, &tmp[7], 8);
}

inline void TwoArmCavity::evaluateEprop_back_wall_ex(std::complex<double> *tmp)
{
	std::complex<double> *tmp_Eprop_forward;
        
	tmp_Eprop_forward=interpolateEfp_x1();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[0], 8);
	cblas_zscal(cavity_transverse_points_num, &Ef_transp_x1, &tmp[0], 8);
        tmp_Eprop_forward=interpolateEfp_x1_tp1();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[1], 8);
	cblas_zscal(cavity_transverse_points_num, &Ef_transp_x1, &tmp[1], 8);
        tmp_Eprop_forward=interpolateEfm_x1();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[2], 8);
	cblas_zscal(cavity_transverse_points_num, &Ef_ex_transp_x1, &tmp[2], 8);
        tmp_Eprop_forward=interpolateEfm_x1_tp1();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[3], 8);
	cblas_zscal(cavity_transverse_points_num, &Ef_ex_transp_x1, &tmp[3], 8);
	
	std::complex<double> tmp_Eprop_backward[cavity_transverse_points_num];	
	getEbp_back_wall(tmp_Eprop_backward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[4], 8);
	getEbp_back_wall_tp1(tmp_Eprop_backward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[5], 8);
	getEbm_back_wall(tmp_Eprop_backward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[6], 8);
	getEbm_back_wall_tp1(tmp_Eprop_backward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[7], 8);
}

inline void TwoArmCavity::evaluateEprop_front_wall_ex(std::complex<double> *tmp)
{
	std::complex<double> dummy=0.0;
	std::complex<double> tmp_Eprop_forward[cavity_transverse_points_num];	
	getEfp_front_wall(tmp_Eprop_forward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[0], 8);
	getEfp_front_wall_tp1(tmp_Eprop_forward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[1], 8);
	getEfm_front_wall(tmp_Eprop_forward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[2], 8);
	getEfm_front_wall_tp1(tmp_Eprop_forward);
        cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_forward, 1, &tmp[3], 8);
	
	std::complex<double> *tmp_Eprop_backward; 

	tmp_Eprop_backward=interpolateEbp_x0();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[4], 8);
	cblas_zscal(cavity_transverse_points_num, &Eb_transp_x0, &tmp[4], 8);
        tmp_Eprop_backward=interpolateEbp_x0_tp1();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[5], 8);
	cblas_zscal(cavity_transverse_points_num, &Eb_transp_x0, &tmp[5], 8);
        tmp_Eprop_backward=interpolateEbm_x0();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[6], 8);
	cblas_zscal(cavity_transverse_points_num, &Eb_ex_transp_x0, &tmp[6], 8);
        tmp_Eprop_backward=interpolateEbm_x0_tp1();
	cblas_zcopy(cavity_transverse_points_num, tmp_Eprop_backward, 1, &tmp[7], 8);
	cblas_zscal(cavity_transverse_points_num, &Eb_ex_transp_x0, &tmp[7], 8);
}

inline void TwoArmCavity::get_transfer_matrix_extraAxis(std::complex<double> *a11_ex, std::complex<double> *a12_ex, std::complex<double> *a21_ex, std::complex<double> *a22_ex)
{
		*a11_ex = transfer_matrix_a11_ex;
		*a12_ex = transfer_matrix_a12_ex;
		*a21_ex = transfer_matrix_a21_ex;
		*a22_ex = transfer_matrix_a22_ex;
}

inline void TwoArmCavity::get_transfer_matrix(std::complex<double> *a11, std::complex<double> *a12, std::complex<double> *a21, std::complex<double> *a22, std::complex<double> **MPol_fp_t, std::complex<double> **MPol_fm_t, std::complex<double> **MPol_bp_t, std::complex<double> **MPol_bm_t)
{
		*a11 = transfer_matrix_a11;
		*a12 = transfer_matrix_a12;
		*a21 = transfer_matrix_a21;
		*a22 = transfer_matrix_a22;

		*MPol_fp_t = transfer_matrix_MacPol_fp;
		*MPol_fm_t = transfer_matrix_MacPol_fm;
		*MPol_bp_t = transfer_matrix_MacPol_bp;
		*MPol_bm_t = transfer_matrix_MacPol_bm;
}

// nk -> Refractive index of next material, if there is an angle, then it should already be included!
// Tk -> Transfer function for next material
inline void TwoArmCavity::set_transfer_matrix(std::complex<double> nk, std::complex<double> T_Em, double loss_pluss, double loss_minus)
{
	std::complex<double> ni = getRefInd() + I*getRefInd_im();
	double cos_th_i_left, cos_th_i_right; getCosTh(&cos_th_i_left, &cos_th_i_right);
	ni = ni*cos_th_i_right;

	std::complex<double> T_Ep = get_transport_Ef_x1();

        transfer_matrix_a11 = 2.0*ni*T_Ep/(ni+nk);
        transfer_matrix_a12 = (nk-ni)*T_Em/(ni+nk);
        transfer_matrix_a21 = (ni-nk)*T_Ep/(ni+nk);
        transfer_matrix_a22 = 2.0*nk*T_Em/(ni+nk);
	
	// Loss is never during a QW, thus does not need to be inside MacPol
	double lp = sqrt(1.0-loss_pluss);
	double lm = sqrt(1.0-loss_minus);

	transfer_matrix_a11 *= lp;
	transfer_matrix_a12 *= lp;
	transfer_matrix_a21 *= lm;
	transfer_matrix_a22 *= lm;

	// Macroscpic polarization
	transfer_matrix_b = I*mu0*c0*w0/(ni+nk);

}

// nk -> Refractive index of next material, if there is an angle, then it should already be included!
// Tk -> Transfer function for next material
inline void TwoArmCavity::set_transfer_matrix_extraAxis(std::complex<double> n_pre, std::complex<double> n_post)
{
	std::complex<double> ni = getRefInd() + I*getRefInd_im();
	std::complex<double> T_Ep_ex = get_transport_Ef_ex_x1();
	std::complex<double> T_Em_ex = get_transport_Eb_ex_x0();

        transfer_matrix_a11_ex = 2.0*n_pre*T_Ep_ex/(ni+n_pre);
        transfer_matrix_a12_ex = (n_post-ni)*T_Em_ex/(ni+n_post);
        transfer_matrix_a21_ex = (ni-n_pre)*T_Ep_ex/(ni+n_pre);
        transfer_matrix_a22_ex = 2.0*n_post*T_Em_ex/(ni+n_post);
}

inline void TwoArmCavity::set_transfer_matrix_macPol_fp(std::complex<double> *mpol)
{
        cblas_zcopy(cavity_transverse_points_num, mpol, 4, transfer_matrix_MacPol_fp, 1);
	// Multiply by constant
	cblas_zscal(cavity_transverse_points_num, &transfer_matrix_b, transfer_matrix_MacPol_fp , 1);
}

inline void TwoArmCavity::set_transfer_matrix_macPol_fm(std::complex<double> *mpol)
{
        cblas_zcopy(cavity_transverse_points_num, mpol, 4, transfer_matrix_MacPol_fm, 1);
	// Multiply by constant
	cblas_zscal(cavity_transverse_points_num, &transfer_matrix_b, transfer_matrix_MacPol_fm , 1);
}

inline void TwoArmCavity::set_transfer_matrix_macPol_bp(std::complex<double> *mpol)
{
        cblas_zcopy(cavity_transverse_points_num, mpol, 4, transfer_matrix_MacPol_bp, 1);
	// Multiply by constant
	cblas_zscal(cavity_transverse_points_num, &transfer_matrix_b, transfer_matrix_MacPol_bp , 1);
}

inline void TwoArmCavity::set_transfer_matrix_macPol_bm(std::complex<double> *mpol)
{
        cblas_zcopy(cavity_transverse_points_num, mpol, 4, transfer_matrix_MacPol_bm, 1);
	// Multiply by constant
	cblas_zscal(cavity_transverse_points_num, &transfer_matrix_b, transfer_matrix_MacPol_bm , 1);
}
inline std::complex<double> TwoArmCavity::get_transport_Ef_x1(void)
{
        return Ef_transp_x1;
}
inline std::complex<double> TwoArmCavity::get_transport_Eb_x0(void)
{
        return Eb_transp_x0;
}
inline std::complex<double> TwoArmCavity::get_transport_Ef_ex_x1(void)
{
        return Ef_ex_transp_x1;
}
inline std::complex<double> TwoArmCavity::get_transport_Eb_ex_x0(void)
{
        return Eb_ex_transp_x0;
}
inline std::complex<double> * TwoArmCavity::interpolateEbp_x0()
{
	return electric_field_bp->get2ndLastContainer();
}

inline std::complex<double> * TwoArmCavity::interpolateEfp_x1()
{
	return electric_field_fp->get2ndLastContainer();
}

inline std::complex<double> * TwoArmCavity::interpolateEbm_x0()
{
	return electric_field_bm->get2ndLastContainer();
}

inline std::complex<double> * TwoArmCavity::interpolateEfm_x1()
{
	return electric_field_fm->get2ndLastContainer();
}

inline std::complex<double> * TwoArmCavity::interpolateEbp_x0_tp1()
{
	return electric_field_bp->get3rdLastContainer();
}

inline std::complex<double> * TwoArmCavity::interpolateEfp_x1_tp1()
{
	return electric_field_fp->get3rdLastContainer();
}

inline std::complex<double> * TwoArmCavity::interpolateEbm_x0_tp1()
{
	return electric_field_bm->get3rdLastContainer();
}

inline std::complex<double> * TwoArmCavity::interpolateEfm_x1_tp1()
{
	return electric_field_fm->get3rdLastContainer();
}

#endif









