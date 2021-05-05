#ifndef CLASSBPM_H
#define CLASSBPM_H

#include <stdio.h>

#include "constantsAndMiscUnits.h"

#include <iostream>
using namespace std;

#if defined(__ICC) || defined(__INTEL_COMPILER)
//#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#include "mkl_dfti.h"

#elif defined(__GNUC__) || defined(__GNUG__)
#include <fftw3.h>
extern "C"{
	// X = a*X
	void cblas_zscal(const int N, std::complex<double>* alpha, std::complex<double>* X, const int LDX);
}
#endif



//! Beam propagation method
/*! FFT propagator for BPM\n
    Free space propagator: Propagate a function u(x)exp(I k z) for a distance Z. Will ZERO pad domain with a hard coded zero for added resolution.\n
    Lens propagator: Use equivalent ABCD cavity with components: (gain chip) Lens1 -> Free space -> Lens2 (SESAM). The relevant parameters are set outside of this code.\n
    Note: When comparing with GAUSSIAN beam solution we have to use PARAXIAL propagator, turn off BOUNDARY GUARD, and might have to set z=-z in Gaussian solution\n
	\n
	\n
    // Use this x grid to get optimal error\n
    double dx = x_max/((double)num_x);\n
    double x[num_x];\n
    for(int i = 0; i < num_x; i++)\n
    {\n
	double index = ((double)i)*((double)num_x-1.0)/((double)num_x-1.0);\n
	x[i] = dx*( (double)i - (double)num_x/2.0);\n
    }\n
    Note: Optimized for using FFTW3 with g++ and MLK-FFT with icc compilers. GSL-FFT is sub-optimal.\n
*/
class BPM
{
	public:
		//! A constructor
		BPM()
		{
			prop = NULL;
			lens = NULL;
			lens2 = NULL;
			boundary_guard_2 = NULL;
			k_r = NULL;
			
			ZERO_PADDING = 0;
			ZERO_PADDING_DOMAIN_SIZE = 0.0;

			N = 0;
			boundary_guard_ratio = 0.8;
			

		#if defined(__ICC) || defined(__INTEL_COMPILER)
			A_fft = NULL;
			
		#elif defined(__GNUC__) || defined(__GNUG__)
			dfft_x = NULL;
			dfft_w = NULL;
		#endif
		}

		//! A constructor
		BPM(int numEl, double R_max)
		{
			prop = NULL;
			lens = NULL;
			lens2 = NULL;
			boundary_guard_2 = NULL;
			k_r = NULL;

			if (numEl <= 0)
			{
				printf("BPM::initialize() num_x must be > 0\n");
				exit(-1);
			} else if ((numEl > 1)&&(numEl % 2 != 0))
			{
				printf("BPM::initialize() Use num_x even, not odd\n");
				exit(-1);
			}

			N = numEl;
			ZERO_PADDING = 0*N; // Added to both sides i.e. # of exta elements are 2*ZERO_PADDING


			// Minimal zero padding to ensure a well sized domain
			// Also ensures that the TOTAL number of elements N + 2*ZERO_PADDING = 2^a
			ZERO_PADDING_DOMAIN_SIZE = 3000.0*um;
			if  (ZERO_PADDING > 0)
			{
				if (R_max < ZERO_PADDING_DOMAIN_SIZE)
				{
					int a = ceil(log2((double)N*ZERO_PADDING_DOMAIN_SIZE/R_max));
					ZERO_PADDING = floor(0.5*(pow(2,a) - N));
				} else {
					// Domain size is ok, fix N such that the TOTAL number of elements N + 2*ZERO_PADDING = 2^a
					int a = ceil(log2((double)N));
					ZERO_PADDING = floor(0.5*(pow(2,a) - N));
				}
			}


			boundary_guard_ratio = 0.8;

			prop = new std::complex<double>[N+2*ZERO_PADDING];
			lens = new std::complex<double>[N+2*ZERO_PADDING];
			lens2 = new std::complex<double>[N+2*ZERO_PADDING];
			boundary_guard_2 = new double[N];

			k_r = new double[N+2*ZERO_PADDING];

			for(int i = 0; i < N; i++)
			{
				boundary_guard_2[i] = 1.0;
			}
			
			for(int i = 0; i < N+2*ZERO_PADDING; i++)
			{
				prop[i] = 0.0;
				lens[i] = 0.0;
				lens2[i] = 0.0;
			}
			
		#if defined(__ICC) || defined(__INTEL_COMPILER)
			A_fft = NULL;
			A_fft = new double[2*N+ 2*2*ZERO_PADDING];
			for(int i = 0; i < 2*N + 2*2*ZERO_PADDING; i++)
			{
				A_fft[i] = 0;
			}

			status = DftiCreateDescriptor( &mkl_fft_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, N+2*ZERO_PADDING);
			if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
			{
				cout << "classBPM::BPM(int,double) DftiCreateDescriptor: Message = " << DftiErrorMessage(status) << endl;
				exit(-1);
			}
			status = DftiSetValue( mkl_fft_handle, DFTI_THREAD_LIMIT, 0);
			if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
			{
				cout << "classBPM::BPM(int,double) DftiSetValue: Message = " << DftiErrorMessage(status) << endl;
				exit(-1);
			}
			status = DftiSetValue( mkl_fft_handle, DFTI_NUMBER_OF_USER_THREADS, 0);
			if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
			{
				cout << "classBPM::BPM(copy) DftiSetValue: Message = " << DftiErrorMessage(status) << endl;
				exit(-1);
			}
			status = DftiCommitDescriptor( mkl_fft_handle);
			if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
			{
				cout << "classBPM::BPM(int,double) DftiCommitDescriptor: Message = " << DftiErrorMessage(status) << endl;
				exit(-1);
			}
			
		#elif defined(__GNUC__) || defined(__GNUG__)
			dfft_x = NULL;
			dfft_w = NULL;
			dfft_x = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N+2*ZERO_PADDING));
			dfft_w = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N+2*ZERO_PADDING));
			plan_fft  = fftw_plan_dft_1d(N+2*ZERO_PADDING, dfft_x, dfft_w, FFTW_FORWARD, FFTW_MEASURE);
			plan_ifft = fftw_plan_dft_1d(N+2*ZERO_PADDING, dfft_w, dfft_x, FFTW_BACKWARD, FFTW_MEASURE);
		#endif
		}

		//! A destructor
		~BPM()
		{
			if (prop!=NULL)
			{
				delete [] prop;
				delete [] lens;
				delete [] lens2;
				delete [] boundary_guard_2;
				delete [] k_r;
				
			}
		#if defined(__ICC) || defined(__INTEL_COMPILER)
			if (A_fft != NULL)
			{
				delete [] A_fft;
				status = DftiFreeDescriptor(&mkl_fft_handle);
			}
			
		#elif defined(__GNUC__) || defined(__GNUG__)
			if (dfft_x != NULL)
			{
				fftw_destroy_plan(plan_fft);
				fftw_destroy_plan(plan_ifft);
				fftw_free(dfft_x);
				fftw_free(dfft_w);
			}
		#endif
			
		}

		//! A copy-constructor
		BPM(const BPM & obj)
		{
			prop = NULL;
			lens = NULL;
			lens2 = NULL;
			boundary_guard_2 = NULL;
			k_r = NULL;

			ZERO_PADDING_DOMAIN_SIZE = obj.ZERO_PADDING_DOMAIN_SIZE;

			if (obj.prop!=NULL)
			{
				N = obj.N;
				ZERO_PADDING = obj.ZERO_PADDING; // Added to both sides i.e. # of exta elements are 2*ZERO_PADDING
	//			carrier_freq = obj.carrier_freq;
				boundary_guard_ratio = obj.boundary_guard_ratio;
				
				
				prop = new std::complex<double>[N + 2*ZERO_PADDING];
				lens = new std::complex<double>[N + 2*ZERO_PADDING];
				lens2 = new std::complex<double>[N + 2*ZERO_PADDING];
				boundary_guard_2 = new double[N];

				k_r = new double[N+2*ZERO_PADDING];
				
				
				for(int i = 0; i < N; i++)
				{
					boundary_guard_2[i] = 0.0;
				}

				for(int i = 0; i < N+2*ZERO_PADDING; i++)
				{
					prop[i] = obj.prop[i];
					lens[i] = obj.lens[i];
					lens2[i] = obj.lens2[i];
				}
			}

		#if defined(__ICC) || defined(__INTEL_COMPILER)
			A_fft = NULL;

			if (obj.A_fft!=NULL)
			{
				A_fft = new double[2*N + 2*2*ZERO_PADDING];
				for(int i = 0; i < 2*N + 2*2*ZERO_PADDING; i++)
				{
					A_fft[i] = obj.A_fft[i];
				}

				status = DftiCreateDescriptor( &mkl_fft_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, N+2*ZERO_PADDING);
				if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
				{
					cout << "classBPM::BPM(copy) DftiCreateDescriptor: Message = " << DftiErrorMessage(status) << endl;
					exit(-1);
				}
				status = DftiSetValue( mkl_fft_handle, DFTI_THREAD_LIMIT, 0);
				if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
				{
					cout << "classBPM::BPM(int,double) DftiSetValue: Message = " << DftiErrorMessage(status) << endl;
					exit(-1);
				}
				status = DftiSetValue( mkl_fft_handle, DFTI_NUMBER_OF_USER_THREADS, 0);
				if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
				{
					cout << "classBPM::BPM(copy) DftiSetValue: Message = " << DftiErrorMessage(status) << endl;
					exit(-1);
				}
				status = DftiCommitDescriptor( mkl_fft_handle);
				if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
				{
					cout << "classBPM::BPM(copy) DftiCommitDescriptor: Message = " << DftiErrorMessage(status) << endl;
					exit(-1);
				}
			}
			
		#elif defined(__GNUC__) || defined(__GNUG__)
			dfft_x = NULL;
			dfft_w = NULL;

			if (obj.dfft_x != NULL)
			{
				dfft_x = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N+2*ZERO_PADDING));
				dfft_w = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N+2*ZERO_PADDING));
				plan_fft  = fftw_plan_dft_1d(N+2*ZERO_PADDING, dfft_x, dfft_w, FFTW_FORWARD, FFTW_MEASURE);
				plan_ifft = fftw_plan_dft_1d(N+2*ZERO_PADDING, dfft_w, dfft_x, FFTW_BACKWARD, FFTW_MEASURE);
			}
		#endif
		}

		//! Construct FreeSpace propagator
		void initialize_freeSpace_FFT_BPM(double target_Z, double x_max, double *x, double lambda, double bg_ratio)
		{
			if (x_max <= 0)
			{
				printf("BPM::initialize() x_max must be > 0\n");
				exit(-1);
			} else if (lambda <= 0.0)
			{
				printf("BPM::initialize() lambda must be  > 0.0\n");
				exit(-1);
			} else if (boundary_guard_ratio < 0.0)
			{
				printf("BPM::initialize() boundary_guard_ratio must be >= 0\n");
				exit(-1);
			}

			double k0 = 2.0*Pi/lambda;
			double dk = 2.0*Pi/(x_max*(1.0 + 2.0*(double)ZERO_PADDING/(double)N));
			
			boundary_guard_ratio = bg_ratio;

/*		

			// Use this x grid to get optimal error
			double dx = x_max/((double)N);
			double x[N];
			for(int i = 0; i < N; i++)
			{
				double index = ((double)i)*((double)N-1.0)/((double)N-1.0);
				x[i] = dx*( (double)i - (double)N/2.0);
			}
*/

			if (N > 1)
			{
			
				for(int i = 0; i < (N+2*ZERO_PADDING)/2.0; i++)
				{
					k_r[i] = dk*i;
				}
				for(int i = (N+2*ZERO_PADDING)/2.0; i < N+2*ZERO_PADDING; i++)
				{
					k_r[i] = dk*(i - (N+2*ZERO_PADDING));
				}
				
				for(int i = 0; i < N + 2*ZERO_PADDING; i++)
				{
				//	prop[i] = exp(-I*k_r[i]*k_r[i]*target_Z/(2.0*k0)); // paraxial
					prop[i] = exp( I*target_Z*(sqrt(k0*k0 - (k_r[i]*k_r[i])) - k0)); // full propagator
				}

				// Add Boundary guard. Super gaussian
				if (boundary_guard_ratio > 0.0)
				{
					int degree = 16;

					double sigma = boundary_guard_ratio*x_max/(2.0*sqrt(2.0)*std::pow(log(2.0),1.0/(2.0*degree)));
					double b = std::pow(sqrt(2.0)*sigma,2.0*degree);
					for(int i = 0; i < N; i++)
					{
						double dx = std::pow(x[i],2.0*degree);
						boundary_guard_2[i] = exp(-dx/b);
					}

				} else {
					for(int i = 0; i < N; i++)
					{
						boundary_guard_2[i] = 1.0;
					}
					
				}
			} else {
				prop[0] = 1.0;
				boundary_guard_2[0] = 1.0;
			}
			
			
			// Add CARRIER frequency
			//carrier_freq = exp(I*k0*target_Z);
		}

		//! Construct Lens1+FreeSpace+Lens2 propagator
		/*! Equvalent cavity to any other using ABCD matrix
			\param focus1 Lens at start of propgation
			\param length air-gap to propagate through
		        \param focus2 Lens at end of propagation
		*/
		void initialize_focusing_lens_FFT_BPM(double focus1_new, double f_new, double focus2_new, double x_max, double *x, double lambda, double bg_ratio)
		{
			if (x_max <= 0)
			{
				printf("BPM::initialize() x_max must be > 0\n");
				exit(-1);
			} else if (lambda <= 0.0)
			{
				printf("BPM::initialize() lambda must be  > 0.0\n");
				exit(-1);
			} else if (boundary_guard_ratio < 0.0)
			{
				printf("BPM::initialize() boundary_guard_ratio must be >= 0\n");
				exit(-1);
			} else if ((focus1_new<0.0)||(focus2_new<0.0))
			{
				printf("BPM::initialize() focus1 and focus2 must be >= 0\n");
				exit(-1);
			} else if (f_new <= 0.0)
			{
				printf("BPM::initialize() air gap length must be > 0\n");
				exit(-1);
			} else if (focus1_new == focus2_new)
			{
				printf("BPM::initialize() focus1 != focus2 \n");
				exit(-1);
			}


	//		double waist0 = (1.0/sqrt(2.0))*FWHM0/sqrt(2.0*log(2.0)); // fwhm of |E(z)|
	//		double zR1 = Pi*waist0*waist0/lambda;
	//		double waist_z = waist0*sqrt(1.0+f1*f1/(zR1*zR1));
	//		cout << "LENS: f1 = " << f1 << " [m], w(0) = " << waist0/um << " [um], w(f1) = " << waist_z/um << " [um]" << endl;

			double x_max_new = x_max*(1.0 + 2.0*(double)ZERO_PADDING/(double)N);

			double k0 = 2.0*Pi/lambda;
			double dk = 2.0*Pi/x_max_new;
			
			boundary_guard_ratio = bg_ratio;


/*		
			// Use this x grid to get optimal error
			double dx = x_max/((double)N);
			double x[N];
			for(int i = 0; i < N; i++)
			{
				double index = ((double)i)*((double)N-1.0)/((double)N-1.0);
				x[i] = dx*( (double)i - (double)N/2.0);
			}
*/

			if (N > 1)
			{
				for(int i = 0; i < (N+2*ZERO_PADDING)/2.0; i++)
				{
					k_r[i] = dk*i;
				}
				for(int i = (N+2*ZERO_PADDING)/2.0; i < N+2*ZERO_PADDING; i++)
				{
					k_r[i] = dk*(i - (N+2*ZERO_PADDING));
				}
				
				for(int i = 0; i < N + 2*ZERO_PADDING; i++)
				{
					//	prop[i] = exp(-I*k[i]*k[i]*f_new/(2.0*k0)); // paraxial
						prop[i] = exp( I*f_new*(sqrt(k0*k0 - (k_r[i]*k_r[i])) - k0)); // full propagator
				}
				

				// STEP 2: Lens -> Plane
				// Apply Lens
				for(int i = 0; i < N + 2*ZERO_PADDING; i++)
				{
					double x_target;
					if (ZERO_PADDING > 0)
					{
						double dx = x[1]-x[0];
						x_target = dx*( i - ZERO_PADDING - (double)N/2.0);
					} else {
						x_target = x[i];
					}


					double r = sqrt(x_target*x_target);
					if(focus1_new>0)
					{
					lens[i]  = exp(-I*k0*r*r/(2.0*focus1_new));
					}else
					{
					lens[i]  = 1.0;
					}
					if(focus2_new>0)
					{
					lens2[i]  = exp(-I*k0*r*r/(2.0*focus2_new));
					}else
					{
					lens2[i]  = 1.0;
					}
				}


				// Add Boundary guard. Super gaussian
				if (boundary_guard_ratio > 0.0)
				{

					int degree = 16;
					
					// Small BG
					double sigma = boundary_guard_ratio*x_max/(2.0*sqrt(2.0)*std::pow(log(2.0),1.0/(2.0*degree)));
					double b = std::pow(sqrt(2.0)*sigma,2.0*degree);
					for(int i = 0; i < N; i++)
					{
						double dx = std::pow(x[i],2.0*degree);
						boundary_guard_2[i] = exp(-dx/b);
					}

					for(int i = 0; i < N; i++)
					{
						lens[i] = lens[i]*boundary_guard_2[i];
					}

					

				} else {
				
					for(int i = 0; i < N; i++)
					{
						boundary_guard_2[i] = 1.0;
					}
				}

			
	
			} else {
				prop[0] = 1.0;
				lens[0]  = 1.0;
				lens2[0]  = 1.0;
				boundary_guard_2[0] = 1.0;
			}
			
			// Add CARRIER frequency
			//carrier_freq = exp(I*k0*(f1_new));
		}

		//! Compute FreeSpace propagator on transverse grid x[]
		void compute_freeSpace_FFT_BPM(std::complex<double> *x)
		{


		#if defined(__ICC) || defined(__INTEL_COMPILER)
			memset(A_fft, 0, (2*N + 2*2*ZERO_PADDING)*sizeof(double)); // A_fft is longer than x
			memcpy(&A_fft[2*ZERO_PADDING], &(reinterpret_cast<double*>(x)[0]), 2*N*sizeof(double));

			// Step 0: Boundary Guard
			for(int i =0; i < N; i++)
			{
				A_fft[2*(ZERO_PADDING+i)]	= A_fft[2*(ZERO_PADDING+i)]*boundary_guard_2[i];
				A_fft[2*(ZERO_PADDING+i)+1]	= A_fft[2*(ZERO_PADDING+i)+1]*boundary_guard_2[i];
			}

			// Step 1: FFT
			status = DftiComputeForward( mkl_fft_handle, A_fft);


			// Step 2: Propagate
			double *prop_tmp = reinterpret_cast<double*>(prop);
			for(int i =0; i < N + 2*ZERO_PADDING; i++)
			{
				double re_tmp = A_fft[2*i]*prop_tmp[2*i  ] - A_fft[2*i+1]*prop_tmp[2*i+1];
				double im_tmp = A_fft[2*i]*prop_tmp[2*i+1] + A_fft[2*i+1]*prop_tmp[2*i];
				A_fft[2*i  ]	= re_tmp;
				A_fft[2*i+1]	= im_tmp;
			}

			// Step 3: Inverse spectrum
			status = DftiComputeBackward( mkl_fft_handle, A_fft);
			
/*
			// Step 4: Add CARRIER frequency and transfer into output array
			for(int i =0 ; i < N; i++)
			{
				b[i]	= (A_fft[2*i] + I*A_fft[2*i+1])*carrier_freq;
			}
			return b;
*/

			// Step 5: In place change, carrier freq taken outside
			memcpy(&(reinterpret_cast<double*>(x)[0]), &A_fft[2*ZERO_PADDING], 2*N*sizeof(double));

			// Step 6: Scale with length of FFT array for each inverse transform
			std::complex<double> scale = 1.0/(N + 2*ZERO_PADDING);
			cblas_zscal(N, &scale, x, 1);

			
		#elif defined(__GNUC__) || defined(__GNUG__)
			memset(dfft_x, 0, (N + 2*ZERO_PADDING)*sizeof(fftw_complex));
			memcpy(&(reinterpret_cast<double*>(dfft_x)[2*ZERO_PADDING]), &(reinterpret_cast<double*>(x)[0]), 2*N*sizeof(double));

			// Step 0: Boundary Guard
			double *prop_tmp = reinterpret_cast<double*>(dfft_x);
			for(int i =0; i < N; i++)
			{
				prop_tmp[2*(ZERO_PADDING+i)]	= prop_tmp[2*(ZERO_PADDING+i)]*boundary_guard_2[i];
				prop_tmp[2*(ZERO_PADDING+i)+1]	= prop_tmp[2*(ZERO_PADDING+i)+1]*boundary_guard_2[i];
			}

			// Step 1: FFT
			fftw_execute(plan_fft);

			// Step 2: Propagate
			double *data_tmp = reinterpret_cast<double*>(dfft_w);
			prop_tmp = reinterpret_cast<double*>(prop);
			for(int i =0; i < N + 2*ZERO_PADDING; i++)
			{
				double re_tmp = data_tmp[2*i]*prop_tmp[2*i  ] - data_tmp[2*i+1]*prop_tmp[2*i+1];
				double im_tmp = data_tmp[2*i]*prop_tmp[2*i+1] + data_tmp[2*i+1]*prop_tmp[2*i];
				data_tmp[2*i  ]	= re_tmp;
				data_tmp[2*i+1]	= im_tmp;
			}

			// Step 3: Inverse spectrum
			fftw_execute(plan_ifft);

/*
			// Step 4: Add CARRIER frequency and transfer into output array
			for(int i =0 ; i < N; i++)
			{
				b[i]	= (A_fft[2*i] + I*A_fft[2*i+1])*carrier_freq;
			}
			return b;
*/

			// Step 5: In place change, carrier freq taken outside
			memcpy(&(reinterpret_cast<double*>(x)[0]), &(reinterpret_cast<double*>(dfft_x)[2*ZERO_PADDING]), 2*N*sizeof(double));

			// Step 6: Scale with length of FFT array for each inverse transform
			std::complex<double> scale = 1.0/(N + 2*ZERO_PADDING);
			cblas_zscal(N, &scale, x, 1);
		#endif

			
		}

		//! Compute FreeSpace+Lens propagator on transverse grid x[]
		void compute_lens_FFT_BPM(std::complex<double> *x)
		{
		#if defined(__ICC) || defined(__INTEL_COMPILER)
			memset(A_fft, 0, (2*N + 2*2*ZERO_PADDING)*sizeof(double)); // A_fft is longer than x
			memcpy(&A_fft[2*ZERO_PADDING], &(reinterpret_cast<double*>(x)[0]), 2*N*sizeof(double));

			// Step 0: Apply first lens
			// Contains a boundary guard
			double *lens_tmp = reinterpret_cast<double*>(lens);
			for(int i =0; i < N+2*ZERO_PADDING; i++)
			{
				double re_tmp = A_fft[2*i]*lens_tmp[2*i  ] - A_fft[2*i+1]*lens_tmp[2*i+1];
				double im_tmp = A_fft[2*i]*lens_tmp[2*i+1] + A_fft[2*i+1]*lens_tmp[2*i];
				A_fft[2*i  ]	= re_tmp;
				A_fft[2*i+1]	= im_tmp;
			}

			// Step 1: FFT
			status = DftiComputeForward( mkl_fft_handle, A_fft);

			// Step 2: Propagate
			double *prop_tmp = reinterpret_cast<double*>(prop);
			for(int i =0; i < N + 2*ZERO_PADDING; i++)
			{
				double re_tmp = A_fft[2*i]*prop_tmp[2*i  ] - A_fft[2*i+1]*prop_tmp[2*i+1];
				double im_tmp = A_fft[2*i]*prop_tmp[2*i+1] + A_fft[2*i+1]*prop_tmp[2*i];
				A_fft[2*i  ]	= re_tmp;
				A_fft[2*i+1]	= im_tmp;
			}

			// Step 3: Inverse spectrum
			status = DftiComputeBackward( mkl_fft_handle, A_fft);
			
			// Step 4: Apply second Lens
			lens_tmp = reinterpret_cast<double*>(lens2);
			for(int i =0; i < N+2*ZERO_PADDING; i++)
			{
				double re_tmp = A_fft[2*i]*lens_tmp[2*i  ] - A_fft[2*i+1]*lens_tmp[2*i+1];
				double im_tmp = A_fft[2*i]*lens_tmp[2*i+1] + A_fft[2*i+1]*lens_tmp[2*i];
				A_fft[2*i  ]	= re_tmp;
				A_fft[2*i+1]	= im_tmp;
			}

/*
			// Step 8: Add CARRIER frequency and transfer into output array
			for(int i =0 ; i < N; i++)
			{
				b[i]	= (A_fft[2*i] + I*A_fft[2*i+1])*carrier_freq;
			}
			return b;
*/

			// Step 9: In place change, carrier freq taken outside
			memcpy(&(reinterpret_cast<double*>(x)[0]), &A_fft[2*ZERO_PADDING], 2*N*sizeof(double));

			// Step 10: Scale with length of FFT array for each inverse transform
			std::complex<double> scale = 1.0/(N + 2*ZERO_PADDING);
			cblas_zscal(N, &scale, x, 1);

		#elif defined(__GNUC__) || defined(__GNUG__)
			memset(dfft_x, 0, (N + 2*ZERO_PADDING)*sizeof(fftw_complex));
			memcpy(&(reinterpret_cast<double*>(dfft_x)[2*ZERO_PADDING]), &(reinterpret_cast<double*>(x)[0]), 2*N*sizeof(double));

			// Step 0: Apply Lens
			// Contains a boundary guard
			double *data_tmp = reinterpret_cast<double*>(dfft_x);
			double *lens_tmp = reinterpret_cast<double*>(lens);
			for(int i =0; i < N+2*ZERO_PADDING; i++)
			{
				double re_tmp = data_tmp[2*i]*lens_tmp[2*i  ] - data_tmp[2*i+1]*lens_tmp[2*i+1];
				double im_tmp = data_tmp[2*i]*lens_tmp[2*i+1] + data_tmp[2*i+1]*lens_tmp[2*i];
				data_tmp[2*i  ]	= re_tmp;
				data_tmp[2*i+1]	= im_tmp;
			}

			// Step 1: FFT
			fftw_execute(plan_fft);

			// Step 2: Propagate
			data_tmp = reinterpret_cast<double*>(dfft_w);
			double *prop_tmp = reinterpret_cast<double*>(prop);
			for(int i =0; i < N + 2*ZERO_PADDING; i++)
			{
				double re_tmp = data_tmp[2*i]*prop_tmp[2*i  ] - data_tmp[2*i+1]*prop_tmp[2*i+1];
				double im_tmp = data_tmp[2*i]*prop_tmp[2*i+1] + data_tmp[2*i+1]*prop_tmp[2*i];
				data_tmp[2*i  ]	= re_tmp;
				data_tmp[2*i+1]	= im_tmp;
			}

			// Step 3: Inverse spectrum
			fftw_execute(plan_ifft);

			// Step 4: Apply second Lens
			data_tmp = reinterpret_cast<double*>(dfft_x);
			lens_tmp = reinterpret_cast<double*>(lens2);
			for(int i =0; i < N+2*ZERO_PADDING; i++)
			{
				double re_tmp = data_tmp[2*i]*lens_tmp[2*i  ] - data_tmp[2*i+1]*lens_tmp[2*i+1];
				double im_tmp = data_tmp[2*i]*lens_tmp[2*i+1] + data_tmp[2*i+1]*lens_tmp[2*i];
				data_tmp[2*i  ]	= re_tmp;
				data_tmp[2*i+1]	= im_tmp;
			}
/*
			// Step 8: Add CARRIER frequency and transfer into output array
			for(int i =0 ; i < N; i++)
			{
				b[i]	= (A_fft[2*i] + I*A_fft[2*i+1])*carrier_freq;
			}
			return b;
*/

			// Step 9: In place change, carrier freq taken outside
			memcpy(&(reinterpret_cast<double*>(x)[0]), &(reinterpret_cast<double*>(dfft_x)[2*ZERO_PADDING]), 2*N*sizeof(double));

			// Step 10: Scale with length of FFT array for each inverse transform
			std::complex<double> scale = 1.0/(N + 2*ZERO_PADDING);
			cblas_zscal(N, &scale, x, 1);
		#endif

		}
	private:
		
		std::complex<double> *prop; // BPM propagator
		std::complex<double> *lens; // BPM propagator
		std::complex<double> *lens2; // BPM propagator
		double *boundary_guard_2; // BPM boundary guard in real space
		
		int N; // Number of transverse points
		double *k_r; // Momentum grid
		double boundary_guard_ratio; // Ratio of domain used for zero padding (ZERO_PADDING_DOMAIN_SIZE)
		int ZERO_PADDING; // Number of zero padded elements in the fourier transform (Increases domain size, reduces need for BG)
		double ZERO_PADDING_DOMAIN_SIZE; // minimal size that the zero padding extends out to
		
	#if defined(__ICC) || defined(__INTEL_COMPILER)
		double *A_fft; // Size is 2*N for complex numbers
		DFTI_DESCRIPTOR_HANDLE mkl_fft_handle;
		MKL_LONG status;
	#elif defined(__GNUC__) || defined(__GNUG__)
		fftw_complex *dfft_x;
		fftw_complex *dfft_w;
		fftw_plan plan_fft;
		fftw_plan plan_ifft;
	#endif
};

#endif






