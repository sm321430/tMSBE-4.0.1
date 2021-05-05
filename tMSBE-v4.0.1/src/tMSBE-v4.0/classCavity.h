
	

#ifndef __CAVITY_H_INCLUDED__
#define __CAVITY_H_INCLUDED__

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
	The material CAVITY class stores the electromagnetic (EM) fields inside each material layer
	the number of points that are depend on the given timestep.

	In order to iterate the cavity fileds one has to call: updateStorate() \n
	There are various methods here that can be used to recover information about the EM field
	or set the EM field to a desired value

	Standalone useage: \n
	Cavity *a = new Cavity(refractive_index, cavity_width, position_1, position_2, dt); // Create object \n
	a->initializeZero(dt, NUM_TRANS, *transv_y); \n
	// .. \n
	// For each timestep \n
	a->setEpluss(array1); \n
	a->setEminus(array1); \n
	a->updateStorage(); // Prepare for next timestep, clear up single space \n
	\n

	\sa Module \n
	2018: Transverse fields added	\n
	2016: Initial version \n
	Isak Kilen
*/

class Cavity
{
	public:
		//! A constructor
		Cavity()
		{
			setName("CAV");
			setToFileOutputKey("out_");
			setRefInd(0.0);
			setRefInd_im(0.0);
			setCosTh(1.0, 1.0);
			setWidth(0.0);
			setPosition0(0.0);
			setPosition1(0.0);
			setNumberOfTimesteps(0);
			time_delay = -1;
			cavity_dt = -1;
			electric_field_pluss = NULL;
			electric_field_minus = NULL;
			transfer_matrix_a11 = 0.0;
			transfer_matrix_a12 = 0.0;
			transfer_matrix_a21 = 0.0;
			transfer_matrix_a22 = 0.0;

			transfer_matrix_b = 0.0;
			transfer_matrix_MacPol = NULL;

			cavity_transverse_points_num = 0;
			cavity_transverse_points_y = NULL;
			output_E_real = NULL;
			output_E_imag = NULL;
			output_E_pluss_real = NULL;
			output_E_pluss_imag = NULL;
			output_E_minus_real = NULL;
			output_E_minus_imag = NULL;

			freeSpace_propagator = NULL;
			lens_propagator_pluss = NULL;
			lens_propagator_minus = NULL;

			temp_transverse_array_pluss = NULL;
			temp_transverse_array_minus = NULL;

			cavity_eq_opt_LENS1_focus = 0.0;
			cavity_eq_opt_LENS2_focus = 0.0;
			cavity_eq_opt_length = 0.0;
		}
		//! A destructor
		~Cavity()
		{
			if (electric_field_pluss != NULL)
			{
				delete  electric_field_pluss;
				delete  electric_field_minus;
				delete [] cavity_transverse_points_y;
				delete [] transfer_matrix_MacPol;

				delete [] output_E_real;
				delete [] output_E_imag;
				delete [] output_E_pluss_real;
				delete [] output_E_pluss_imag;
				delete [] output_E_minus_real;
				delete [] output_E_minus_imag;

				delete [] temp_transverse_array_pluss;
				delete [] temp_transverse_array_minus;

			}

			if (freeSpace_propagator != NULL)
			{
				delete freeSpace_propagator;
			}

			if (lens_propagator_pluss != NULL)
			{
				delete lens_propagator_pluss;
				delete lens_propagator_minus;
			}
		}
		//! A copy constructor
		Cavity(const Cavity &obj)
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

			
			cavity_transverse_points_num = obj.cavity_transverse_points_num;
			cavity_transverse_points_y = NULL;

			electric_field_pluss = NULL;
			electric_field_minus = NULL;

			freeSpace_propagator = NULL;
			lens_propagator_pluss = NULL;
			lens_propagator_minus = NULL;
			temp_transverse_array_pluss = NULL;
			temp_transverse_array_minus = NULL;

			if (obj.electric_field_pluss != NULL)
			{
				cavity_transverse_points_y = new double[cavity_transverse_points_num];
				
				electric_field_pluss = new Cyclic<std::complex<double> >(*obj.electric_field_pluss);
				electric_field_minus = new Cyclic<std::complex<double> >(*obj.electric_field_minus);
				
				for(int i = 0; i < cavity_transverse_points_num; i++)
				{
					cavity_transverse_points_y[i] = obj.cavity_transverse_points_y[i];
				}
				transfer_matrix_MacPol = new std::complex<double>[cavity_transverse_points_num];
				output_E_real = new std::ofstream[cavity_transverse_points_num];
				output_E_imag = new std::ofstream[cavity_transverse_points_num];
				output_E_pluss_real = new std::ofstream[cavity_transverse_points_num];
				output_E_pluss_imag = new std::ofstream[cavity_transverse_points_num];
				output_E_minus_real = new std::ofstream[cavity_transverse_points_num];
				output_E_minus_imag = new std::ofstream[cavity_transverse_points_num];

				temp_transverse_array_pluss = new std::complex<double>[cavity_transverse_points_num];
				temp_transverse_array_minus = new std::complex<double>[cavity_transverse_points_num];
				
			}
			if (obj.freeSpace_propagator != NULL)
			{
				freeSpace_propagator = new BPM(*obj.freeSpace_propagator);
			}

			if (obj.lens_propagator_pluss != NULL)
			{
				lens_propagator_pluss = new BPM(*obj.lens_propagator_pluss);
				lens_propagator_minus = new BPM(*obj.lens_propagator_minus);
			}

			cavity_eq_opt_LENS1_focus = obj.cavity_eq_opt_LENS1_focus;
			cavity_eq_opt_LENS2_focus = obj.cavity_eq_opt_LENS2_focus;
			cavity_eq_opt_length = obj.cavity_eq_opt_length;
		}

		//! A constructor
		Cavity(double refInd, double width, double x0, double x1, double dt)
		{
			setName("CAV");
			setToFileOutputKey("out_");
			setRefInd(refInd);
			setRefInd_im(0.0);
			setCosTh(1.0,1.0);
			setWidth(width);
			setPosition0(x0);
			setPosition1(x1);
			setNumberOfTimesteps(0);
			time_delay = -1;
			cavity_dt = dt;

			transfer_matrix_a11 = 0.0;
			transfer_matrix_a12 = 0.0;
			transfer_matrix_a21 = 0.0;
			transfer_matrix_a22 = 0.0;

			transfer_matrix_b = 0.0;
			transfer_matrix_MacPol = NULL;

			electric_field_pluss = NULL;
			electric_field_minus = NULL;

			cavity_transverse_points_num = -1;
			cavity_transverse_points_y = NULL;
			temp_transverse_array_pluss = NULL;
			temp_transverse_array_minus = NULL;
		}

		//! Print an overview of the cavity data to screen
		void Print(void)
		{
			cout << "Print cavity:" << endl;
			cout << " -> name    = " << getName() << endl;
			cout << " -> n       = " << getRefInd() << " +i " << getRefInd_im() << endl;
			cout << " -> width   = " << getWidth()/um << " [um]" << endl;
			cout << " -> x0      = " << getPosition0() << endl;
			cout << " -> x1      = " << getPosition1() << endl;
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
/*
			cout << "E+ (r,z):" << endl;
			// Print field (diagnostics only)
			std::complex<double> ep[cavity_transverse_points_num];
			for(int i = 0; i < cavity_transverse_points_num; i++)
			{
				cout << "r = " << cavity_transverse_points_y[i]/um << "\t: ";

				for(int j = 0; j < getDelay()+2; j++)
				{
					getEpluss(j,ep);
					cout << ep[i] << " ";
				}
				cout << endl;
			}

			cout << "E- (r,z):" << endl;
			// Print field (diagnostics only)
			for(int i = 0; i < cavity_transverse_points_num; i++)
			{
				cout << "r = " << cavity_transverse_points_y[i]/um << "\t: ";

				for(int j = 0; j < getDelay()+2; j++)
				{
					getEminus(j,ep);
					cout << ep[i] << " ";
				}
				cout << endl;
			}
*/
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
		
		//! Set first value of E+
		void setEpluss(std::complex<double>*);

		//! Set first value of E-
		void setEminus(std::complex<double>*);

		//! Return pointer to memory of first value of E+
		/*! Overloaded function for setting memory*/
		std::complex<double> * setEpluss(void);

		//! Return pointer to memory of first value of E-
		/*! Overloaded function for setting memory*/
		std::complex<double> * setEminus(void);

		//! Return first value of E+
		void getEpluss(int, std::complex<double>*) const;

		//! Return first value of E-
		void getEminus(int, std::complex<double>*) const;
		
		//! Initialize all storage in this class
		/*! This function initializes the class using the simulation setup+timestep
		    It is assumed that the cavity has already been setup width required dimensions
		    \param DT timestep of simulation
		    \param Nx Number of transverse points
		    \param x Transverse points x in [-R_max/2, (1/2 - 1/N)R_max]
		    \param R_max Width of transverse grid
		    \param bg_ratio Ratio of domain used as boundary guard (SuperGaussian) for FFT-BPM propagator
		*/
		void initializeZero(double DT, int Nx, double*x, double R_max, double boundary_guard_rato);

		//! Save E+/E- to binary files in folder 'save/'
		/*! Saves E+/E- to binary files named with 'cavity_name' and transverse dimension number
			If folder 'save/' does not exist in directory nothing can be saved
			If files already exist, they will be overwritten.
		 */
		void file_save_variables(int, int);
		//! Load E+/E-  from binary files in folder 'save/'
		/*! Load E+/E- files from files named with 'cavity_name' and transverse dimension number
			This function assumes that the files in 'save/' have same number of timesteps as cavity
			If files cannot be found, the program will stop 
 		 */
		void file_load_variables(int, int);

		//! Update storage of E+/E-
		/*! This function iterates the storage structures such that a new timestep of E+/E- is ready
		    to be stored. The new storage space is NOT set to zero..\n
		    Note: Computes FFT BPM if needed
		 */
		void updateStorage_freeSpace();

		//! Update storage of E+/E- with forced BPM
		/*! This function iterates the storage structures such that a new timestep of E+/E- is ready
		    to be stored. The new storage space is NOT set to zero..\n
		    Note: Computes FFT BPM if needed
		 */
		void updateStorage_freeSpace_forcedBPM();

		//! Update storage of E+
		/*! This function iterates the storage structures such that a new timestep of E+ is ready
		    to be stored. The new storage space is NOT set to zero..\n
		    Note: Computes FFT BPM if needed
		 */
		void updateStorage_lens_pluss();

		//! Update storage of E-
		/*! This function iterates the storage structures such that a new timestep of E+ is ready
		    to be stored. The new storage space is NOT set to zero..\n
		    Note: Computes FFT BPM if needed
		 */
		void updateStorage_lens_minus();

		//! Set all components of E+/E- to zero
		void clearFields();
		
		//! Write all simulation variables to binary file(s)
		/*! The specific output data can be changed in this function.
		    Complex numbers are saved in two files named with 're' and 'im'
		    If you change the output data, also change the functions: file_output_open and file_output_close
		 */
		void file_output_write(int);
		
		//! Write all simulation variables to binary file(s) at left wall
		/*! The specific output data can be changed in this function.
		    Complex numbers are saved in two files named with 're' and 'im'
		    If you change the output data, also change the functions: file_output_open and file_output_close
		 */
		void file_output_write_leftWall(int);

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
		
		//! Return E+ evaluated on the left boundary
		void getEpluss_left_wall(std::complex<double> *);

		//! Return E- evaluated on the left boundary
		/*! Note: calls 'interpolateEminus_x0()' to evaluate on the boundary */
		void getEminus_left_wall(std::complex<double> *);

		//! Return E+ evaluated on the right boundary
		/*! Note: calls 'interpolateEpluss_x1()' to evaluate on the boundary */
		void getEpluss_right_wall(std::complex<double> *);

		//! Return E- evaluated on the right boundary
		void getEminus_right_wall(std::complex<double> *);

		//! Return E+ evaluated on the left boundary the previous timestep
		void getEpluss_left_wall_tp1(std::complex<double> *);

		//! Return E- evaluated on the left boundary the previous timestep
		/*! Note: calls 'interpolateEminus_x0()' to evaluate on the boundary */
		void getEminus_left_wall_tp1(std::complex<double> *);

		//! Return E+ evaluated on the right boundary the previous timestep
		/*! Note: calls 'interpolateEpluss_x1()' to evaluate on the boundary */
		void getEpluss_right_wall_tp1(std::complex<double> *);

		//! Return E- evaluated on the right boundary the previous timestep
		void getEminus_right_wall_tp1(std::complex<double> *);
		
		//! Return E- on the left boundary (possibly) interpolated between two timesteps
		void interpolateEminus_x0(std::complex<double> *);

		std::complex<double> * interpolateEminus_x0(void);

		//! Return E+ on the right boundary linearly interpolated between two timesteps
		void interpolateEpluss_x1(std::complex<double> *);

		std::complex<double> * interpolateEpluss_x1(void);

		//! Return E- from the past timestep evaluated on the left boundary linearly interpolated between two timesteps
		void interpolateEminus_x0_tp1(std::complex<double> *);

		std::complex<double> * interpolateEminus_x0_tp1(void);

		//! Return E+ from the past timestep evaluated on the right boundary linearly interpolated between two timesteps
		void interpolateEpluss_x1_tp1(std::complex<double> *);

		std::complex<double> * interpolateEpluss_x1_tp1(void);

		//! Return E- at a given 'z' position inside the materal layer. Uses linear interpolation.
		void interpolateEminus(double, std::complex<double> *);

		//! Return E+ at a given 'z' position inside the material layer. Uses linear interpolation.
		void interpolateEpluss(double, std::complex<double> *);

		//! Return E(z,t) at a given 'z' position inside the material layer. Uses linear interpolation.
		void evaluateEprop(double, std::complex<double> *);

		//! Return E(z,t) at the right boundary.
		void evaluateEprop_x1(std::complex<double> *);

		//! Return E(z,t) at the right boundary.
		void evaluateEprop_x1_fast(std::complex<double> *, int stride);

		//! Return E(z,t) at the left boundary.
		void evaluateEprop_x0(std::complex<double> *);

		//! Return E(z,t) at the right boundary from the previous timestep
		void evaluateEprop_x1_tp1(std::complex<double> *);

		//! Return E(z,t) at the right boundary from the previous timestep
		void evaluateEprop_x1_tp1_fast(std::complex<double> *, int stride);

		//! Return E(z,t) at the left boundary from the previous timestep
		void evaluateEprop_x0_tp1(std::complex<double> *);

		//! Return the element needed to compute transfer matrix method.
		/*! This method requires that the internal variables have been initialized using 'set_transfer_matrix()' and 'set_transfer_matrix_macPol()'. The matrix elements a11, a12, a21, a22 and the macroscopic polarization is returned. */
		void get_transfer_matrix(std::complex<double> *a11, std::complex<double> *a12, std::complex<double> *a21, std::complex<double> *a22, std::complex<double> **b);

		//! Initialize the transfer matrix method variables for the given cavity
		/*! Computes the matrix elements a11, a12, a21, a22 and the macroscopic polarization constants. Includes any angular dependence. */
		void set_transfer_matrix(std::complex<double>, std::complex<double>, double, double);

		//! Set the macroscopic polarization that will later be used in the transfer matrix method
		void set_transfer_matrix_macPol(std::complex<double> *);
		
		//! Set equivalent optical cavity paramters
		/*! Only used when a LENS is in the cavity\n
		    \param focus1 focus on gain chip side
		    \param length free space propagation length
		    \param focus focus on SESAM side*/
		void set_equivalent_optical_cavity(double FOCUS_1, double LENGTH, double FOCUS_2);

		//! Returns the transport phase factor for E+ in this cavity
                std::complex<double> get_transport_Ep_x1(void);

                //! Returns the transport phase factor for E- in this cavity
                std::complex<double> get_transport_Em_x0(void);
	
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

		//! Transport phase factor for E+ to transfer from left to left boundary (1)
		std::complex<double> Epluss_transp_x0;

		//! Transport phase factor for E- to transfer from right to left boundary
		std::complex<double> Eminus_transp_x0;

		//! Transport phase factor for E+ to transfer from left to right boundary
		std::complex<double> Epluss_transp_x1;

		//! Transport phase factor for E- to transfer from right to right boundary (1)
		std::complex<double> Eminus_transp_x1;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a11;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a12;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a21;

		//! Transfer matrix element precomputed with refractive indices
		std::complex<double> transfer_matrix_a22;

		//! Transfer matrix  macroscopic polarization factors
		std::complex<double> transfer_matrix_b;

		//! Transfer matrix storage current for macroscopic polarization
		std::complex<double> *transfer_matrix_MacPol;
		
		//! Electric field E+ propagating from left boundary to right boundary
		Cyclic<std::complex<double> > *electric_field_pluss;

		//! Electric field E- propagating from right boundary to left boundary
		Cyclic<std::complex<double> > *electric_field_minus;

		//! FFT BPM free space propagator
		BPM *freeSpace_propagator;
		//! FFT BPM lens space propagator for E+
		BPM *lens_propagator_pluss;
		//! FFT BPM lens space propagator for E-
		BPM *lens_propagator_minus;

		//! Equivalent cavity Lens1 focus
		double cavity_eq_opt_LENS1_focus;

		//! Equivalent cavity Lens2 focus
		double cavity_eq_opt_LENS2_focus;

		//! Equivalent cavity air gap length
		double cavity_eq_opt_length;
		
		//! Number of stored timesteps in E+ and E-
		int numberOfTimesteps;

		//! Temporary array to avoid alot of heap-creation/destruction
		std::complex<double> *temp_transverse_array_pluss;

		//! Temporary array to avoid alot of heap-creation/destruction
		std::complex<double> *temp_transverse_array_minus;
		
		//! Output to file of real(E)
		std::ofstream *output_E_real;

		//! Output to file of imag(E)
		std::ofstream *output_E_imag;

		//! Output to file of real(E+)
		std::ofstream *output_E_pluss_real;

		//! Output to file of imag(E+)
		std::ofstream *output_E_pluss_imag;

		//! Output to file of real(E-)
		std::ofstream *output_E_minus_real;

		//! Output to file of imag(E-)
		std::ofstream *output_E_minus_imag;
};

inline void Cavity::setName(const std::string & newName)
{
	cavity_name = newName;
}

inline std::string Cavity::getName(void) const
{
	return cavity_name;
}

inline void Cavity::setCosTh(double new_cos_th_left, double new_cos_th_right)
{
	cavity_left_cos_th  = new_cos_th_left;
	cavity_right_cos_th = new_cos_th_right;
}

inline void Cavity::getCosTh(double *left, double *right) const
{
	*left  = cavity_left_cos_th;
	*right = cavity_right_cos_th;
}

inline void Cavity::setToFileOutputKey(const std::string & newName)
{
	cavity_output_key = newName;
}
inline std::string Cavity::getToFileOutputKey(void) const
{
	return cavity_output_key;
}

inline void Cavity::setRefInd(double newInd)
{
	cavity_refractive_index = newInd;
}

inline void Cavity::setRefInd_im(double newInd)
{
	cavity_refractive_index_im = newInd;
}

inline double Cavity::getRefInd(void) const
{
	return cavity_refractive_index;
}

inline double Cavity::getRefInd_im(void) const
{
	return cavity_refractive_index_im;
}

inline void Cavity::setWidth(double newWidth)
{
	cavity_width = newWidth;
}

inline double Cavity::getWidth(void) const
{
	return cavity_width;
}

inline void Cavity::setPosition0(double newX0)
{
	cavity_position_x0 = newX0;
}

inline double Cavity::getPosition0(void) const
{
	return cavity_position_x0;
}

inline void Cavity::setPosition1(double newX1)
{
	cavity_position_x1 = newX1;
}

inline double Cavity::getPosition1(void) const
{
	return cavity_position_x1;
}

inline void Cavity::setNumberOfTimesteps(int nt)
{
	numberOfTimesteps = nt;
}

inline int Cavity::getNumberOfTimesteps(void) const
{
	return numberOfTimesteps;
}

inline int Cavity::getDelay(void) const
{
	return time_delay;
}

inline void Cavity::setDelay(int delay)
{
	time_delay = delay;
}

inline int Cavity::getCavityDT(void) const
{
	return cavity_dt;
}

inline void Cavity::setCavityDT(double newDT)
{
	cavity_dt = newDT;
}

inline void Cavity::setEpluss(std::complex<double> *target)
{
	electric_field_pluss->setFirstContainer(target);
}

inline void Cavity::setEminus(std::complex<double> *target)
{
	electric_field_minus->setFirstContainer(target);
}

inline std::complex<double> * Cavity::setEpluss()
{
	return electric_field_pluss->getFirstContainer();
}

inline std::complex<double> * Cavity::setEminus()
{
	return electric_field_minus->getFirstContainer();
}

inline void Cavity::getEpluss(int k, std::complex<double> *tmp) const
{
	memcpy(tmp, electric_field_pluss->getContainerNr(k), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void Cavity::getEminus(int k, std::complex<double> *tmp) const
{	
	memcpy(tmp, electric_field_minus->getContainerNr(k), cavity_transverse_points_num*sizeof(std::complex<double>));
}

inline void Cavity::getEpluss_left_wall(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_pluss->getFirstContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
/*
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i]*Epluss_transp_x0; // *1
	}
*/
}

inline void Cavity::getEminus_left_wall(std::complex<double> *tmp)
{
	interpolateEminus_x0(tmp);
	cblas_zscal(cavity_transverse_points_num, &Eminus_transp_x0, tmp , 1);
/*
	for(int i = 0; i < cavity_transverse_points_num; i++)
        {
                tmp[i] = tmp[i]*Eminus_transp_x0;
        }
*/
}

inline void Cavity::getEpluss_right_wall(std::complex<double> *tmp)
{
	interpolateEpluss_x1(tmp);
	cblas_zscal(cavity_transverse_points_num, &Epluss_transp_x1, tmp , 1);
/*
	for(int i = 0; i < cavity_transverse_points_num; i++)
        {
                tmp[i] = tmp[i]*Epluss_transp_x1;
        }
*/
}

inline void Cavity::getEminus_right_wall(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_minus->getFirstContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
/*
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i]*Eminus_transp_x1; // *1
	}
*/
}


inline void Cavity::getEpluss_left_wall_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_pluss->get2ndContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
/*
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i]*Epluss_transp_x0; // *1
	}
*/
}

inline void Cavity::getEminus_left_wall_tp1(std::complex<double> *tmp)
{
	interpolateEminus_x0_tp1(tmp);
	cblas_zscal(cavity_transverse_points_num, &Eminus_transp_x0, tmp , 1);
/*
	for(int i = 0; i < cavity_transverse_points_num; i++)
        {
                tmp[i] = tmp[i]*Eminus_transp_x0;
        }
*/
}

inline void Cavity::getEpluss_right_wall_tp1(std::complex<double> *tmp)
{
	interpolateEpluss_x1_tp1(tmp);
	cblas_zscal(cavity_transverse_points_num, &Epluss_transp_x1, tmp , 1);
/*
	for(int i = 0; i < cavity_transverse_points_num; i++)
        {
                tmp[i] = tmp[i]*Epluss_transp_x1;
        }
*/
}

inline void Cavity::getEminus_right_wall_tp1(std::complex<double> *tmp)
{
	memcpy(tmp, electric_field_minus->get2ndContainer(), cavity_transverse_points_num*sizeof(std::complex<double>));
/*
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i]*Eminus_transp_x1; // *1
	}
*/
}


inline void Cavity::evaluateEprop_x1(std::complex<double> *tmp)
{
	getEpluss_right_wall(temp_transverse_array_pluss);
	getEminus_right_wall(tmp);

	std::complex<double> dummy = 1.0;
	cblas_zaxpy(cavity_transverse_points_num, &dummy, temp_transverse_array_pluss  , 1, tmp, 1);
/*
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + temp_transverse_array_pluss[i];
	}
*/
}

inline void Cavity::evaluateEprop_x1_fast(std::complex<double> *tmp, int stride)
{
	cblas_zcopy(cavity_transverse_points_num, electric_field_minus->getFirstContainer(), 1, tmp, stride);

	std::complex<double> *tmp_Ep_x1 = interpolateEpluss_x1();
	cblas_zaxpy(cavity_transverse_points_num, &Epluss_transp_x1, tmp_Ep_x1, 1, tmp, stride);
}


inline void Cavity::evaluateEprop_x0(std::complex<double> *tmp)
{
	getEpluss_left_wall(tmp);
	getEminus_left_wall(temp_transverse_array_minus);

	std::complex<double> dummy = 1.0;
	cblas_zaxpy(cavity_transverse_points_num, &dummy, temp_transverse_array_minus  , 1, tmp, 1);
/*
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + temp_transverse_array_minus[i];
	}
*/
}

inline void Cavity::evaluateEprop_x1_tp1(std::complex<double> *tmp)
{
	getEpluss_right_wall_tp1(temp_transverse_array_pluss);
	getEminus_right_wall_tp1(tmp);

	std::complex<double> dummy = 1.0;
	cblas_zaxpy(cavity_transverse_points_num, &dummy, temp_transverse_array_pluss  , 1, tmp, 1);
/*
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + temp_transverse_array_pluss[i];
	}
*/
}

inline void Cavity::evaluateEprop_x1_tp1_fast(std::complex<double> *tmp, int stride)
{
	cblas_zcopy(cavity_transverse_points_num, electric_field_minus->get2ndContainer(), 1, tmp, stride);

	std::complex<double> *tmp_Ep_x1 = interpolateEpluss_x1_tp1();
	cblas_zaxpy(cavity_transverse_points_num, &Epluss_transp_x1, tmp_Ep_x1, 1, tmp, stride);
}

inline void Cavity::evaluateEprop_x0_tp1(std::complex<double> *tmp)
{
	getEpluss_left_wall_tp1(tmp);
	getEminus_left_wall_tp1(temp_transverse_array_minus);

	std::complex<double> dummy = 1.0;
	cblas_zaxpy(cavity_transverse_points_num, &dummy, temp_transverse_array_minus  , 1, tmp, 1);
/*
	for(int i = 0; i < cavity_transverse_points_num; i++)
	{
		tmp[i] = tmp[i] + temp_transverse_array_minus[i];
	}
*/
}

inline void Cavity::get_transfer_matrix(std::complex<double> *a11, std::complex<double> *a12, std::complex<double> *a21, std::complex<double> *a22, std::complex<double> **MPol_t)
{
		*a11 = transfer_matrix_a11;
		*a12 = transfer_matrix_a12;
		*a21 = transfer_matrix_a21;
		*a22 = transfer_matrix_a22;

		*MPol_t = transfer_matrix_MacPol;
/*
		for(int i = 0; i < cavity_transverse_points_num; i++)
		{
			MPol_t[i] = transfer_matrix_b*transfer_matrix_MacPol[i];
		}
*/
}

// nk -> Refractive index of next material, if there is an angle, then it should already be included!
// Tk -> Transfer function for next material
inline void Cavity::set_transfer_matrix(std::complex<double> nk, std::complex<double> T_Em, double loss_pluss, double loss_minus)
{
	std::complex<double> ni = getRefInd() + I*getRefInd_im();
	double cos_th_i_left, cos_th_i_right; getCosTh(&cos_th_i_left, &cos_th_i_right);
	ni = ni*cos_th_i_right;

	std::complex<double> T_Ep = get_transport_Ep_x1();

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
inline void Cavity::set_transfer_matrix_macPol(std::complex<double> *mpol)
{
        cblas_zcopy(cavity_transverse_points_num, mpol, 4, transfer_matrix_MacPol, 1);
	//memcpy(transfer_matrix_MacPol, mpol, cavity_transverse_points_num*sizeof(std::complex<double>));

	// Multiply by constant
	cblas_zscal(cavity_transverse_points_num, &transfer_matrix_b, transfer_matrix_MacPol , 1);
}

inline void Cavity::set_equivalent_optical_cavity(double FOCUS_1, double LENGTH, double FOCUS_2)
{
	cavity_eq_opt_LENS1_focus  = FOCUS_1;
	cavity_eq_opt_LENS2_focus  = FOCUS_2;
	cavity_eq_opt_length = LENGTH;
}

inline std::complex<double> Cavity::get_transport_Ep_x1(void)
{
        return Epluss_transp_x1;
}
inline std::complex<double> Cavity::get_transport_Em_x0(void)
{
        return Eminus_transp_x0;
}

inline std::complex<double> * Cavity::interpolateEminus_x0()
{
	return electric_field_minus->get2ndLastContainer();
}

inline std::complex<double> * Cavity::interpolateEpluss_x1()
{
	return electric_field_pluss->get2ndLastContainer();
}

inline std::complex<double> * Cavity::interpolateEminus_x0_tp1()
{
	return electric_field_minus->get3rdLastContainer();
}

inline std::complex<double> * Cavity::interpolateEpluss_x1_tp1()
{
	return electric_field_pluss->get3rdLastContainer();
}

#endif









