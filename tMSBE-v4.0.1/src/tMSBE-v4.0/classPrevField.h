

#ifndef __PREVFIELD_H_INCLUDED__
#define __PREVFIELD_H_INCLUDED__


#include <stdio.h>
#include <iostream>
#include <complex>
#include <cstring> // memcpy
#include "fileIO.cpp"
#include "setup_simulation_variables.h"

#ifdef USE_OPENMP
	#include <omp.h>
#endif

using namespace std;

//! Storage for 'Filter()' and convolution of data with filter function
/*! This class contains stiores an E-filed and a filter function. The purpose is to compute their convolution at given timesteps.\n
    This class is similar to 'Cyclic()' but with additional filter coefficients/methods. They both rotate the front component of the field when the appropriate function is called.  \n
    \sa Cyclic 
    \sa Filter\n
    2018: Add more comments\n
    Isak Kilen
*/
class prev_field{

	public:
		//! A constructor
		prev_field()
		{
			array_length =0;
			array = NULL;
			head_index = 0;


			ind_start = 0;
			ind_stop = 0;
		}

		//! A constructor
		prev_field(int num)
		{
			if (num%2 == 0)
			{
				cout << "prev_field: Number of filter elements should be an ODD number " << endl;
				exit(-1);
			}


			array_length = num;
			array = new std::complex<double>[array_length];
			head_index = 0;


			for(int i = 0; i < array_length; i++)
			{
				array[i] = 0;
			}

			ind_start = 0;
			ind_stop = array_length-1;
		}
		
		//! A destructor
		~prev_field()
		{
			if (array!=NULL)
			{
				delete array;
			}
		}

		//! Copy an obj into current obj
		void copy(prev_field *obj)
		{
			array_length = obj->array_length;
			head_index   = obj->head_index;
			ind_start    = obj->ind_start;
			ind_stop     = obj->ind_stop;
			
			if (array==NULL)
			{
				array = new std::complex<double>[array_length];
			}	
			
			memcpy(array, obj->array, array_length*sizeof(std::complex<double>));
		}

		
		//! Set the head index of the array to a given number 
		void set_head_index(int n)
		{
			head_index = n;
		}

		//! Increment head index of array during iteration
		void update_head_index(void)
		{
			head_index = head_index+1;
			head_index = head_index%array_length;
		}

		//! Set first element
		void set_top_element(std::complex<double> *nw)
		{
			set_element(0,nw);
		}
		
		//! Set element number
		void set_element(int i, std::complex<double> *nw)
		{
			int index = (i-head_index);
			if (index < 0)
			{
				index += array_length;
			}
			array[index] = *nw;
		}

		//! Return element number
		std::complex<double> get_element(int i)
		{
			int index = (i-head_index);
			if (index < 0)
			{
				index += array_length;
			}

			return array[index];
		}

		//! Set all elements to zero
		void zero_all_elements(void)
		{
			for(int i = 0; i < array_length; i++)
			{
				array[i] = 0.0;
			}
		}

		//! Set convolution indices
		void c_array_conv_ind_set(int sta,int fin)
		{
			ind_start = sta;
			ind_stop = fin;
		}

		//! Set convolution indices. Attempt to restrict ind_start and ind_stop based on filter amplitude cutoff
		void c_array_conv_ind_set(double cutoff)
		{
			int ind0 = 0;
			while(abs(array[ind0]) < cutoff)
			{
				ind0 += 1;
			}

			ind_start = ind0;

			int ind1 = array_length-1;
			while(abs(array[ind1]) < cutoff)
			{
				ind1 -= 1;
			}
			ind_stop = ind1;

	//		cout << "[start,stop] = " << ind_start << " / " << ind_stop << " -> " << ind_stop-ind_start << endl;
		}

		//! Set convolution indices. Restrict ind_start and ind_stop based on rato of 'filter volume'
		void c_array_conv_ind_set_center_volume(double ratio)
		{
			// Make a symetic index set from M/2-ind to M/2+ind that contains the
			// "rato" volume of data.
			// rato is a number [0,1]
			double vol_t = 0;
			for(int i = 0; i<array_length; i++)
			{
				vol_t += abs(array[i]);
			}

			int M = (array_length-1.0)/2.0;
			int ind0 = 0;
			while((M-ind0 > 0)&&(M+ind0<array_length))
			{
				// Calculate volume of region
				double vol_r = 0;
				for(int i = M-ind0; i<M+ind0; i++)
				{
					vol_r += abs(array[i]);
				}

				// Termination condition
				if (vol_r/vol_t >= ratio)
				{
					break;
				} else {
					ind0 += 1;
				}
			}

			ind_start = M-ind0;
			ind_stop  = M+ind0;

	//		cout << "1-r = " << 1.0-ratio << " [start,stop] = " << ind_start << " / " << ind_stop << " -> " << ind_stop-ind_start << endl;
		}
		
		//! Calculate convolution of array with filter function
		/*!
			Calculate y = h*x\n

			Where x[n] has head element x[h] and previous element x[h-1]
 			head_index points at the current signal time\n
		*/
		std::complex<double> c_array_conv(std::complex<double> *coeff_array)
		{
			std::complex<double> sum = 0.0;
			double sum_re = 0.0;
			double sum_im = 0.0;
			int index;

			// CIRCULAR CONVOLUTION (TRANSPOSE E(n))
			#ifdef USE_OPENMP
			#pragma omp parallel shared(sum)
			{
				std::complex< double > priv_sum = 0.0;
				#pragma omp for
				for(int i = 0; i < array_length; i++)
				{
					// The data is stored in transposed order so this is the CORRECT sum based on the storage.
					priv_sum += coeff_array[i]*get_element(i); // c[0]*E[0]+c[1]*E[1]+
				}
				#pragma omp critical
				{
					sum += priv_sum;
				}
				

			}
			#endif
			return sum;
		}

		//! Print class information to screen
		void print()
		{
			cout << "array["<< head_index <<"] " << endl;
			for(int i = 0; i < array_length; i++)
			{
				cout << "index = " << i << " -> " << get_element(i) << endl;
			}
		}

		//! Print class information to screen
		void Print_array()
		{
			for(int i = 0; i < array_length; i++)
			{
				printf("%.0f",real(array[i]));
				if (i < array_length-1)
				{
					printf(", ");
				}
			}
			printf("\t|\t");
			for(int i = 0; i < array_length; i++)
			{
				printf("%.0f",real(get_element(i)));
				if (i < array_length-1)
				{
					printf(", ");
				}
			}
		}

		//! Save current array and head_index to binary file
		void saveArray(const std::string &fileName)
		{
			std::stringstream contName;

			contName.str("");
			contName << fileName << "_array.dat";
			saveBinary(contName.str(), array, array_length);

			contName.str("");
			contName << fileName << "_head.dat";
			saveBinary(contName.str(), &head_index, 1);
		}

		//! Load current array and head_index from binary file
		void loadArray(const std::string &fileName)
		{
			std::stringstream contName;

			contName.str("");
			contName << fileName << "_array.dat";
			loadBinary(contName.str(), array, array_length);
			
			contName.str("");
			contName << fileName << "_head.dat";
			loadBinary(contName.str(), &head_index, 1);
		}
		

	private:
		//! The propagating E-field in the cavity
		std::complex<double> *array;
		//! Length of the stored field
		int array_length;

		//! index of the first element in the array
		int head_index;

		//! Lower integration limit for truncated convolution
		int ind_start;

		//! Upper integration limit for truncated convolution
		int ind_stop;
};

#endif


