#ifndef __TWOARMINT_H_INCLUDED__
#define __TWOARMINT_H_INCLUDED__

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



//! The TwoArmInterface class is placed inside the 'Module()' container between other objects.
/*!
	The material TwoArmInterface class stores the electromagnetic (EM) fields inside each material layer
	for a non-normal incience cavity with just a single angle.
the number of points that are depend on the given timestep.

	This class acts as an intermediary between traditional single arm Cavity and Two-Arm Cavities.

	In order to iterate the cavity fields one has to call: updateStorate() \n
	There are various methods here that can be used to recover information about the EM field
	or set the EM field to a desired value

	
	Standalone useage: \n
	Cavity *a = new TwoArmInterface(refractive_index, cavity_width, position_1, position_2, dt); // Create object \n
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

class TwoArmInterface
{
	public:
		//! A constructor
		TwoArmInterface()
		{
			setName("INT");
			setToFileOutputKey("out_");
			setRefInd(0.0);
			setRefInd_im(0.0);
			setReflect(0.0);
			setCosTh(1.0, 1.0);
			setWidth(0.0);
			setPosition0(0.0);
			setPosition1(0.0);
			setNumberOfTimesteps(0);
			time_delay = -1;
			cavity_dt = -1;
			electric_field_fp = NULL;
			electric_field_fm = NULL;
			electric_field_bp = NULL;
			electric_field_bm = NULL;
			transfer_matrix_a11_left = 0.0;
			transfer_matrix_a12_left = 0.0;
			transfer_matrix_a21_left = 0.0;
			transfer_matrix_a22_left = 0.0;
			transfer_matrix_a11_right = 0.0;
			transfer_matrix_a12_right = 0.0;
			transfer_matrix_a21_right = 0.0;
			transfer_matrix_a22_right = 0.0;

			cavity_weights_Efp = NULL;
			cavity_weights_Efm = NULL;
			cavity_weights_Ebp = NULL;
			cavity_weights_Ebm = NULL;
			Efp_trans_delay	   = NULL;
			Efm_trans_delay	   = NULL;
			Ebp_trans_delay	   = NULL;
			Ebm_trans_delay	   = NULL;
			delay_index_Efp	   = NULL;
			delay_index_Efm	   = NULL;
			delay_index_Ebp	   = NULL;
			delay_index_Ebm	   = NULL;

			transfer_matrix_b = 0.0;

			cavity_transverse_points_num = 0;
			cavity_transverse_points_y = NULL;
			temp_transverse_array = NULL;
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

			angle = NULL;		
			intAngle = NULL;		
			prev_cavity = NULL;
			post_cavity = NULL;
			cavity_transverse_phase_in = NULL;
			cavity_transverse_phase_out = NULL;

			freeSpace_propagator = NULL;

		}
		//! A destructor
		~TwoArmInterface()
		{
			if (electric_field_fp != NULL)
			{
				delete  electric_field_fp;
				delete  electric_field_fm;
				delete  electric_field_bp;
				delete  electric_field_bm;
				
				delete  cavity_weights_Efp;
				delete 	cavity_weights_Efm;
				delete 	cavity_weights_Ebp;
				delete 	cavity_weights_Ebm;
				delete 	Efp_trans_delay;
				delete 	Efm_trans_delay;
				delete 	Ebp_trans_delay;
				delete 	Ebm_trans_delay;
				delete 	delay_index_Efp;
				delete 	delay_index_Efm;
				delete 	delay_index_Ebp;
				delete 	delay_index_Ebm;
				
				delete [] cavity_transverse_points_y;
				delete [] cavity_transverse_phase_in;
				delete [] cavity_transverse_phase_out;

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

				delete [] temp_transverse_array;
			}

			if (freeSpace_propagator != NULL)
			{
				delete freeSpace_propagator;
			}
		}
		
		//! A copy constructor
		TwoArmInterface(const TwoArmInterface &obj)
		{
			cavity_name = obj.cavity_name;	
			cavity_output_key = obj.cavity_output_key;
			cavity_refractive_index = obj.cavity_refractive_index;
			cavity_refractive_index_im = obj.cavity_refractive_index_im;
			cavity_left_cos_th = obj.cavity_left_cos_th;
			cavity_right_cos_th = obj.cavity_right_cos_th;
			cavity_width = obj.cavity_width;
			cavity_position_x0 = obj.cavity_position_x0;
			cavity_position_x1 = obj.cavity_position_x1;
			numberOfTimesteps = obj.numberOfTimesteps;
			time_delay = obj.time_delay;
			cavity_dt = obj.cavity_dt;

			reflect=obj.reflect;
			angle=obj.angle;
			intAngle=obj.intAngle;
			prev_cavity=obj.prev_cavity;
			post_cavity=obj.post_cavity;
			cavity_transverse_phase_in=obj.cavity_transverse_phase_in;
			cavity_transverse_phase_out=obj.cavity_transverse_phase_out;
			
			cavity_transverse_points_num = obj.cavity_transverse_points_num;
			cavity_transverse_points_y = NULL;

			temp_transverse_array = NULL;
			
			electric_field_fp = NULL;
			electric_field_fm = NULL;
			electric_field_bp = NULL;
			electric_field_bm = NULL;

			cavity_weights_Efp = NULL;
			cavity_weights_Efm = NULL;
			cavity_weights_Ebp = NULL;
			cavity_weights_Ebm = NULL;
			Efp_trans_delay	   = NULL;
			Efm_trans_delay	   = NULL;
			Ebp_trans_delay	   = NULL;
			Ebm_trans_delay	   = NULL;
			delay_index_Efp	   = NULL;
			delay_index_Efm	   = NULL;
			delay_index_Ebp	   = NULL;
			delay_index_Ebm	   = NULL;

			freeSpace_propagator = NULL;

			if (obj.electric_field_fp != NULL)
			{
				cavity_transverse_points_y = new double[cavity_transverse_points_num];
			
				for( int i = 0; i < cavity_transverse_points_num; i++)
				{	
					electric_field_fp[i] = new Cyclic<std::complex<double> >(*obj.electric_field_fp[i]);
					electric_field_fm[i] = new Cyclic<std::complex<double> >(*obj.electric_field_fm[i]);
					electric_field_bp[i] = new Cyclic<std::complex<double> >(*obj.electric_field_bp[i]);
					electric_field_bm[i] = new Cyclic<std::complex<double> >(*obj.electric_field_bm[i]);
				}		
		
				for(int i = 0; i < cavity_transverse_points_num; i++)
				{
					cavity_transverse_points_y[i] = obj.cavity_transverse_points_y[i];
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
				
			}
			if (obj.freeSpace_propagator != NULL)
			{
				freeSpace_propagator = new BPM(*obj.freeSpace_propagator);
			}
		}

		//! A constructor
		TwoArmInterface(double refInd, double width, double x0, double x1, double dt)
		{
			setName("INT");
			setToFileOutputKey("out_");
			setReflect(0.0);
			setRefInd(refInd);
			setRefInd_im(0.0);
			setCosTh(1.0,1.0);
			setWidth(width);
			setPosition0(x0);
			setPosition1(x1);
			setNumberOfTimesteps(0);
			time_delay = -1;
			cavity_dt = dt;

			transfer_matrix_a11_left = 0.0;
			transfer_matrix_a12_left = 0.0;
			transfer_matrix_a21_left = 0.0;
			transfer_matrix_a22_left = 0.0;
			transfer_matrix_a11_right = 0.0;
			transfer_matrix_a12_right = 0.0;
			transfer_matrix_a21_right = 0.0;
			transfer_matrix_a22_right = 0.0;

			transfer_matrix_b = 0.0;

			electric_field_fp = NULL;
			electric_field_fm = NULL;
			electric_field_bp = NULL;
			electric_field_bm = NULL;

			cavity_weights_Efp = NULL;
			cavity_weights_Efm = NULL;
			cavity_weights_Ebp = NULL;
			cavity_weights_Ebm = NULL;
			Efp_trans_delay	   = NULL;
			Efm_trans_delay	   = NULL;
			Ebp_trans_delay	   = NULL;
			Ebm_trans_delay	   = NULL;
			delay_index_Efp	   = NULL;
			delay_index_Efm	   = NULL;
			delay_index_Ebp	   = NULL;
			delay_index_Ebm	   = NULL;

			cavity_transverse_points_num = -1;
			cavity_transverse_points_y = NULL;
	
			temp_transverse_array = NULL;
			
			angle = NULL;		
			intAngle = NULL;		
			prev_cavity = NULL;
			post_cavity = NULL;
			cavity_transverse_phase_in = NULL;
			cavity_transverse_phase_out = NULL;
		}

		//! Print an overview of the cavity data to screen
		void Print(void)
		{
			cout << "Print two arm cavity:" << endl;
			cout << " -> name    = " << getName() << endl;
			cout << " -> n       = " << getRefInd() << " +i " << getRefInd_im() << endl;
			cout << " -> width   = " << getWidth()/um << " [um]" << endl;
			cout << " -> angle   = " << getAngle() << " rad" <<endl;
			cout << " -> intAngle= " << getIntAngle() << " rad" <<endl;
			cout << " -> cavity  = " << getPrevCav() <<", "<< getPostCav() << endl;
			cout << " -> x0      = " << getPosition0() << endl;
			cout << " -> x1      = " << getPosition1() << endl;
			cout << " -> refl.   = " << getReflect() << endl;
			int dL = getDelay();
			cout << " -> delay   = " << dL << endl;
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

		//! Set the REAL refractive index
		void   setRefInd(double);

		//! Return the reflection coefficient
		std::complex<double> getReflect(void) const;

		//! Set the reflection coefficient
		void   setReflect(std::complex<double>);

		//! Return the IMAG refractive index
		double getRefInd_im(void) const;

		//! Set the IMAG refractive index
		void   setRefInd_im(double);

		//! Return the cavity width
		double getWidth(void) const;

		//! Set the cavity width
		void   setWidth(double);

		//! Return the z-position of the left boundary
		double getPosition0(void) const;

		//! Set the angle of incidence
		void setAngle(double);

		//! Set the angle for introducing transverse phase oscillations
		void setIntAngle(double);

		//! Set the previous cavity
		void setPrevCav(int);
		
		//! Set post cavity
		void setPostCav(int);

		//! Return the angle of incidence
		double getAngle(void) const;

		//! Return the angle of interference fringes
		double getIntAngle(void) const;

		//! Return the previous cavity
		int getPrevCav(void) const;
		
		//! Return post cavity
		int getPostCav(void) const;

		//! Return entering transverse phase increment
		std::complex<double>* getPhaseIn(void) const;

		//! Return exit transverse phase increment
		std::complex<double>* getPhaseOut(void) const;

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

		//! Return cavity resolution in terms of timesteps
		/*! There are #timesteps +2 stored in each container 
		    with the extra two helping with boundary effects
		*/
		int    getDelay(void) const;

		//! Set cavity resolution in terms of timesteps
		/*! There are #timesteps +2 stored in each container 
		    with the extra two helping with boundary effects
		*/
		void   setDelay(int);

		//! Return timestep 'dt' used in simulation
		int    getCavityDT(void) const;

		//! Set timestep 'dt' used in simulation
		void   setCavityDT(double);
		
		//! Return pointer to memory of first value of E forward pluss
		void setEfp(int, std::complex<double>*);

		//! Return pointer to memory of first value of E forward minus
		void setEfm(int, std::complex<double>*);

		//! Return pointer to memory of first value of E backward pluss
		void setEbp(int, std::complex<double>*);

		//! Return pointer to memory of first value of E backward minus
		void setEbm(int, std::complex<double>*);
		
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
	
		//! Return front of forward pluss field
		void getBoundary_Efp(std::complex<double>*);
	
		//! Return front of forward minus field
		void getBoundary_Efm(std::complex<double>*);
	
		//! Return front of backward pluss field
		void getBoundary_Ebp(std::complex<double>*);
	
		//! Return front of backward minus field
		void getBoundary_Ebm(std::complex<double>*);
		
		//! Return forward pluss field at given point
		std::complex<double> getEfp_singlePoint(int, int, std::complex<double>, double);
		
		//! Return forward minus field at given point
		std::complex<double> getEfm_singlePoint(int, int, std::complex<double>, double);
		
		//! Return backward pluss field at given point
		std::complex<double> getEbp_singlePoint(int, int, std::complex<double>, double);
		
		//! Return backward minus field at given point
		std::complex<double> getEbm_singlePoint(int, int, std::complex<double>, double);
			
		//! Return the element needed to compute transfer matrix method.
		/*! This method requires that the internal variables have been initialized using 'set_transfer_matrix()' and 'set_transfer_matrix_macPol()' for each macPol. The matrix elements a11, a12, a21, a22 and the macroscopic polarization is returned. */
		void get_transfer_matrix(std::complex<double> *, std::complex<double> *, std::complex<double> *, std::complex<double> *, std::complex<double> *, std::complex<double> *, std::complex<double> *, std::complex<double>*);

		//! Initialize the transfer matrix method variables for the given cavity
		/*! Computes the matrix elements a11, a12, a21, a22 and the macroscopic polarization constants. Includes any angular dependence. */
		void set_transfer_matrix(std::complex<double>, std::complex<double>, std::complex<double>, std::complex<double>, double, double);

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
	
	private:
		//! Name of device
		std::string cavity_name;

		//! Start of output name of files
		std::string cavity_output_key;

		//! Refractive index of cavity (REAL)
		double cavity_refractive_index;

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

		//! Timestep from the external solver [s]
		double cavity_dt;

		//! Angle of incidence in degrees
		double angle;

		//! Angle for interference frings in degrees
		double intAngle;

		//! Phase induced on entering two arm cavity by translating reference frame
		std::complex<double> *cavity_transverse_phase_in;

		//! Phase induced on exiting two arm cavity by translating reference frame
		std::complex<double> *cavity_transverse_phase_out;

		//! Previous cavity before two arm semiconductor structure
		int prev_cavity;

		//! Post cavity after two arm semiconductor structure
		int post_cavity;

		//! The factor for the angular dependence given by the angle of incidence cos(th) on the left
		double cavity_left_cos_th;

		//! The factor for the angular dependence given by the angle of incidence cos(th) on the right
		double cavity_right_cos_th;

		//! # of transverse points in cavity
		int cavity_transverse_points_num;

		//! Containes the position of each transverse points starting [-R_max, R_max]
		double *cavity_transverse_points_y;
	
		//! Interpolation weight for outgoing fields
		double cavity_weight_x0;

		//! Interpolation weight for right boundary
		double cavity_weight_x1;

		//! Interpolation weights for forward pluss field
		double *cavity_weights_Efp;

		//! Interpolation weights for forward minus field
		double *cavity_weights_Efm;

		//! Interpolation weights for backward pluss field
		double *cavity_weights_Ebp;

		//! Interpolation weights for backward minus field
		double *cavity_weights_Ebm;
		
		//! Phases for transverse delay interpolation in forward pluss direction
		std::complex<double> *Efp_trans_delay;

		//! Phases for transverse delay interpolation in forward minus direction
		std::complex<double> *Efm_trans_delay;

		//! Phases for transverse delay interpolation in backward pluss direction
		std::complex<double> *Ebp_trans_delay;

		//! Phases for transverse delay interpolation in backward minus direction
		std::complex<double> *Ebm_trans_delay;


		//! Indices for transverse delay interpolation in forward pluss direction
		int *delay_index_Efp;

		//! Indices for transverse delay interpolation in forward minus direction
		int *delay_index_Efm;

		//! Indices for transverse delay interpolation in backward pluss direction
		int *delay_index_Ebp;

		//! Indices for transverse delay interpolation in backward minus direction
		int *delay_index_Ebm;

		//! Reflection coefficient
		std::complex<double> reflect;

		//! Temporary array to avoid alot of heap-creation/destruction
		std::complex<double> *temp_transverse_array;

		//! Transport phase factor for E forward to transfer from front to back boundary (1)
		std::complex<double> Ef_transp_x0;

		//! Transport phase factor for E forward to transfer from back to front boundary
		std::complex<double> Ef_transp_x1;

		//! Transport phase factor for E backward to transfer from front to back boundary (1)
		std::complex<double> Eb_transp_x0;

		//! Transport phase factor for E backward to transfer from back to front boundary
		std::complex<double> Eb_transp_x1;

		//! Transport phase factors for E forward plus to transfer a transversally dependent distance
		std::complex<double> *Efp_transp_delay;

		//! Transport phase factors for E forward minus to transfer a transversally dependent distance
		std::complex<double> *Efm_transp_delay;

		//! Transport phase factors for E backward pluss to transfer a transversally dependent distance
		std::complex<double> *Ebp_transp_delay;

		//! Transport phase factors for E backward minus to transfer a transversally dependent distance
		std::complex<double> *Ebm_transp_delay;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a11_left;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a12_left;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a21_left;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a22_left;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a11_right;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a12_right;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a21_right;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a22_right;

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
		Cyclic<std::complex<double> > ** electric_field_fp;

		//! Electric field E forward minus propagating from forward boundary to back boundary
		Cyclic<std::complex<double> > ** electric_field_fm;
		
		//! Electric field E backward pluss propagating from back boundary to front boundary
		Cyclic<std::complex<double> > ** electric_field_bp;

		//! Electric field E backward minus propagating from back boundary to front boundary
		Cyclic<std::complex<double> > ** electric_field_bm;

		//! FFT BPM free space propagator
		BPM *freeSpace_propagator;
	
		//! Number of stored timesteps in E+ and E-
		int numberOfTimesteps;
		
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

inline void TwoArmInterface::setName(const std::string & newName)
{
	cavity_name = newName;
}

inline std::string TwoArmInterface::getName(void) const
{
	return cavity_name;
}

inline void TwoArmInterface::setCosTh(double new_cos_th_left, double new_cos_th_right)
{
	cavity_left_cos_th  = new_cos_th_left;
	cavity_right_cos_th = new_cos_th_right;
}

inline void TwoArmInterface::getCosTh(double *left, double *right) const
{
	*left  = cavity_left_cos_th;
	*right = cavity_right_cos_th;
}

inline void TwoArmInterface::setToFileOutputKey(const std::string & newName)
{
	cavity_output_key = newName;
}
inline std::string TwoArmInterface::getToFileOutputKey(void) const
{
	return cavity_output_key;
}

inline void TwoArmInterface::setRefInd(double newInd)
{
	cavity_refractive_index = newInd;
}

inline void TwoArmInterface::setReflect(std::complex<double> newReflection)
{
	reflect = newReflection;
}

inline std::complex<double> TwoArmInterface::getReflect(void) const
{
	return reflect;
}

inline void TwoArmInterface::setRefInd_im(double newInd)
{
	cavity_refractive_index_im = newInd;
}

inline double TwoArmInterface::getRefInd(void) const
{
	return cavity_refractive_index;
}


inline std::complex<double>* TwoArmInterface::getPhaseIn(void) const
{
	return cavity_transverse_phase_in;
}

inline std::complex<double>* TwoArmInterface::getPhaseOut(void) const
{
	return cavity_transverse_phase_out;
}

inline void TwoArmInterface::setPostCav(int cav_ind)
{
	post_cavity=cav_ind;
}

inline void TwoArmInterface::setPrevCav(int cav_ind)
{
	prev_cavity=cav_ind;
}

inline int TwoArmInterface::getPrevCav(void) const
{
	return prev_cavity;
}

inline int TwoArmInterface::getPostCav(void) const
{
	return post_cavity;
}

inline double TwoArmInterface::getAngle(void) const
{
	return angle;
}

inline double TwoArmInterface::getIntAngle(void) const
{
	return intAngle;
}

inline void TwoArmInterface::setAngle(double ang)
{
	angle=ang;
}

inline void TwoArmInterface::setIntAngle(double ang)
{
	intAngle=ang;
}

inline double TwoArmInterface::getRefInd_im(void) const
{
	return cavity_refractive_index_im;
}

inline void TwoArmInterface::setWidth(double newWidth)
{
	cavity_width = newWidth;
}

inline double TwoArmInterface::getWidth(void) const
{
	return cavity_width;
}

inline void TwoArmInterface::setPosition0(double newX0)
{
	cavity_position_x0 = newX0;
}

inline double TwoArmInterface::getPosition0(void) const
{
	return cavity_position_x0;
}

inline void TwoArmInterface::setPosition1(double newX1)
{
	cavity_position_x1 = newX1;
}

inline double TwoArmInterface::getPosition1(void) const
{
	return cavity_position_x1;
}

inline void TwoArmInterface::setNumberOfTimesteps(int nt)
{
	numberOfTimesteps = nt;
}

inline int TwoArmInterface::getNumberOfTimesteps(void) const
{
	return numberOfTimesteps;
}

inline int TwoArmInterface::getDelay(void) const
{
	return time_delay;
}

inline void TwoArmInterface::setDelay(int delay)
{
	time_delay = delay;
}

inline int TwoArmInterface::getCavityDT(void) const
{
	return cavity_dt;
}

inline void TwoArmInterface::setCavityDT(double newDT)
{
	cavity_dt = newDT;
}

inline void TwoArmInterface::setEfp(int ind, std::complex<double> * value)
{
	electric_field_fp[ind]->setFirstContainer(value);
}

inline void TwoArmInterface::setEfm(int ind, std::complex<double> * value)
{
	electric_field_fm[ind]->setFirstContainer(value);
}

inline void TwoArmInterface::setEbp(int ind, std::complex<double> * value)
{
	electric_field_bp[ind]->setFirstContainer(value);
}

inline void TwoArmInterface::setEbm(int ind, std::complex<double> * value)
{
	electric_field_bm[ind]->setFirstContainer(value);
}


inline void TwoArmInterface::get_transfer_matrix(std::complex<double> *a11_left, std::complex<double> *a12_left, std::complex<double> *a21_left, std::complex<double> *a22_left, std::complex<double> *a11_right, std::complex<double> *a12_right, std::complex<double> *a21_right, std::complex<double> *a22_right)
{
		*a11_left = transfer_matrix_a11_left;
		*a12_left = transfer_matrix_a12_left;
		*a21_left = transfer_matrix_a21_left;
		*a22_left = transfer_matrix_a22_left;
		*a11_right = transfer_matrix_a11_right;
		*a12_right = transfer_matrix_a12_right;
		*a21_right = transfer_matrix_a21_right;
		*a22_right = transfer_matrix_a22_right;
}

// nk -> Refractive index of next material, if there is an angle, then it should already be included!
// Tk -> Transfer function for next material
inline void TwoArmInterface::set_transfer_matrix(std::complex<double> nk, std::complex<double> T_Ep_left, std::complex<double> T_Ep_right, std::complex<double> T_Em, double loss_pluss, double loss_minus)
{
	std::complex<double> ni = getRefInd() + I*getRefInd_im();
	double cos_th_i_left, cos_th_i_right;
	getCosTh(&cos_th_i_left, &cos_th_i_right);
	ni = ni*cos_th_i_right;

	transfer_matrix_a11_left = 2.0*ni*T_Ep_left/(ni+nk);
        transfer_matrix_a12_left = (nk-ni)*T_Em/(ni+nk);
        transfer_matrix_a21_left = (ni-nk)*T_Ep_left/(ni+nk);
        transfer_matrix_a22_left = 2.0*nk*T_Em/(ni+nk);
	transfer_matrix_a11_right = 2.0*ni*T_Ep_right/(ni+nk);
        transfer_matrix_a12_right = (nk-ni)*T_Em/(ni+nk);
        transfer_matrix_a21_right = (ni-nk)*T_Ep_right/(ni+nk);
        transfer_matrix_a22_right = 2.0*nk*T_Em/(ni+nk);
	
	// Loss is never during a QW, thus does not need to be inside MacPol
	double lp = sqrt(1.0-loss_pluss);
	double lm = sqrt(1.0-loss_minus);

	transfer_matrix_a11_left *= lp;
	transfer_matrix_a12_left *= lp;
	transfer_matrix_a21_left *= lm;
	transfer_matrix_a22_left *= lm;
	transfer_matrix_a11_right *= lp;
	transfer_matrix_a12_right *= lp;
	transfer_matrix_a21_right *= lm;
	transfer_matrix_a22_right *= lm;

	// Macroscpic polarization
	transfer_matrix_b = I*mu0*c0*w0/(ni+nk);


}
inline void TwoArmInterface::set_transfer_matrix_macPol_fp(std::complex<double> *mpol)
{
	memcpy(transfer_matrix_MacPol_fp, mpol, cavity_transverse_points_num*sizeof(std::complex<double>));

	// Multiply by constant
	cblas_zscal(cavity_transverse_points_num, &transfer_matrix_b, transfer_matrix_MacPol_fp , 1);
}

inline void TwoArmInterface::set_transfer_matrix_macPol_fm(std::complex<double> *mpol)
{
	memcpy(transfer_matrix_MacPol_fm, mpol, cavity_transverse_points_num*sizeof(std::complex<double>));

	// Multiply by constant
	cblas_zscal(cavity_transverse_points_num, &transfer_matrix_b, transfer_matrix_MacPol_fm , 1);
}

inline void TwoArmInterface::set_transfer_matrix_macPol_bp(std::complex<double> *mpol)
{
	memcpy(transfer_matrix_MacPol_bp, mpol, cavity_transverse_points_num*sizeof(std::complex<double>));

	// Multiply by constant
	cblas_zscal(cavity_transverse_points_num, &transfer_matrix_b, transfer_matrix_MacPol_bp , 1);
}

inline void TwoArmInterface::set_transfer_matrix_macPol_bm(std::complex<double> *mpol)
{
	memcpy(transfer_matrix_MacPol_bm, mpol, cavity_transverse_points_num*sizeof(std::complex<double>));

	// Multiply by constant
	cblas_zscal(cavity_transverse_points_num, &transfer_matrix_b, transfer_matrix_MacPol_bm , 1);
}
inline std::complex<double> TwoArmInterface::get_transport_Ef_x1(void)
{
        return Ef_transp_x1;
}
inline std::complex<double> TwoArmInterface::get_transport_Eb_x0(void)
{
        return Eb_transp_x0;
}

#endif









