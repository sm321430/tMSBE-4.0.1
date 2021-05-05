#ifndef __SORTWEIGHT3_H_INCLUDED__
#define __SORTWEIGHT3_H_INCLUDED__


#include <iostream>

using namespace std;


//! Sort elements in three 1D arrays A1, A2, A3
/*! For the given arrays, the class can find an index set I such that (A1+A2+A3)(I) = B is sorted, as much as possible. \n
    Additional helper functions can apply I to OTHER arrays.\n
    This is a helper function for sorting indices for carrier scattering in the 'Device()' class. \n
    2018: Added more comments \n
    Isak Kilen
*/
class sortWeight3
{
	public:
		//! A constructor
		sortWeight3()
		{
			a_1 = NULL;
			a_2 = NULL;
			a_3 = NULL;
			index = NULL;
			array_length = 0;
		}
		
		//! A constructor
		sortWeight3(int *A_1, int *A_2, int *A_3, int N)
		{
			a_1 = A_1;
			a_2 = A_2;
			a_3 = A_3;
			index = NULL;
			array_length = N;
			

			index = new int[array_length];
			
			for(int i =0; i < array_length; i++)
			{
				index[i] = i;
			}
		}

		//! A destructor
		~sortWeight3()
		{
			if (index != NULL)
			{
				delete [] index;
			}
		}

		//! A copy constructor
		sortWeight3(const sortWeight3 &obj)
		{
			array_length = obj.array_length;
			if (index != NULL)
			{
				index = new int[array_length];
				a_1 = obj.a_1;
				a_2 = obj.a_2;
				a_3 = obj.a_3;
			}
		}
		
		//! Initializes sorting argument function
		/*! \param A_1 First 1D array to be sorted
		    \param A_2 Second 1D array to be sorted
		    \param A_3 Third 1D array to be sorted
	     	    \param N length of array */
		void initSort(int *A_1, int *A_2, int *A_3,int N)
		{
			array_length = N;
			
			a_1 = A_1;
			a_2 = A_2;
			a_3 = A_3;
			
			if(index==NULL)
			{
				index = new int[array_length];
			} else {
				
				delete [] index;
				index = new int[array_length];
			}
			
			for(int i =0; i < array_length; i++)
			{
				index[i] = i;
			}
		}
		
		//! Run sorting on arrays A1, A2 and A3, set with 'initSort()'
		/*! Does not alter array A1, A2, or A3
		   \return a pointer to a list of indices such that (A1+A2+A3)(I) = B is sorted as much as possible */
		int *runSort()
		{
			quick_sort(index, array_length);
			return index;
		}
		
		//! Sort (int) array using index set I and a dummy array for tmp storage
		void applyIndexSort_int(int *array, int *dummy)
		{
			for(int i = 0; i < array_length; i++)
			{
				dummy[i] = array[index[i]];
			}
			
			for(int i = 0; i < array_length; i++)
			{
				array[i] = dummy[i];
			}
		}
		
		//! Sort (double) array using index set I and a dummy array for tmp storage
		void applyIndexSort_double(double *array, double *dummy)
		{
			for(int i = 0; i < array_length; i++)
			{
				dummy[i] = array[index[i]];
			}
			
			for(int i = 0; i < array_length; i++)
			{
				array[i] = dummy[i];
			}
		}
		
		
	
	private:
	
		//! Argument function for quicksort: If F(i) <= F(j)
		bool my3Sort(int i, int j) 
		{
			int sum_i = a_1[i] + a_2[i] + a_3[i];
			int sum_j = a_1[j] + a_2[j] + a_3[j];
			if (sum_i==sum_j)
			{
				int min_i = a_1[i];
				if (a_2[i] < min_i)
				{
					min_i = a_2[i];
				}
				if (a_3[i] < min_i)
				{
					min_i = a_3[i];
				}

				int min_j = a_1[j];
				if (a_2[j] < min_j)
				{
					min_j = a_2[j];
				}
				if (a_3[j] < min_j)
				{
					min_j = a_3[j];
				}

				return (min_i <= min_j);
				
			} else{
				return (sum_i < sum_j);
			}
		}
		
		
		//! Quick sort helper function
		int partition(int list[], int p, int r)
		{
			int index, exchange_temp;
			int pivot = list[r];
			index = p - 1;
			for(int i = p; i < r; i++)
			{
				//if(list[i] <= pivot)
				if (my3Sort(list[i],pivot))
				{
					index++;
					exchange_temp = list[i];
					list[i] = list[index];
					list[index] = exchange_temp;
				}
			}
			exchange_temp = list[r];
			list[r] = list[index+1];
			list[index+1] = exchange_temp;
			return index+1;
		}

		//! Quick sort helper function
		void quicksort_aux(int list[], int p, int r)
		{
			int q;
			if(p<r)
			{
				q = partition(list, p, r);
				quicksort_aux(list, p, q-1);
				quicksort_aux(list, q+1, r);
			}
		}

		//! Quick sort main function
		void quick_sort(int list[], int size)
		{
			quicksort_aux(list,0, size-1);
		}
	
		//! The external array A1 that the index set I should be generated from
		int *a_1;
		//! The external array A2 that the index set I should be generated from
		int *a_2;
		//! The external array A3 that the index set I should be generated from
		int *a_3;
		//! The index set such that (A1+A2)(I) is sorted as much as possible
		int *index;
		//! Length of arrays
		int array_length;
		
};

#endif
