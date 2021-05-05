
#ifndef __SORTWEIGHT1_H_INCLUDED__
#define __SORTWEIGHT1_H_INCLUDED__


#include <iostream>

using namespace std;

//! Sort elements in 1D array A
/*! For a given array A, the class can find an index set I such that A(I) = B is sorted according to argument function. \n
    Additional helper functions can apply I to OTHER arrays.\n
    This is a helper function for sorting indices for carrier scattering in the 'Device()' class. \n
    2018: Added more comments \n
    Isak Kilen
*/
class sortWeight1
{
	public:
		//! A constructor
		sortWeight1()
		{
			a_1 = NULL;
			index = NULL;
			array_length = 0;
		}
		
		//! A constructor
		sortWeight1(int *A_1, int N)
		{
			a_1 = A_1;
			index = NULL;
			array_length = N;
			

			index = new int[array_length];
			
			for(int i =0; i < array_length; i++)
			{
				index[i] = i;
			}
		}

		//! A destructor
		~sortWeight1()
		{
			if (index != NULL)
			{
				delete [] index;
			}
		}

		//! A copy constructor
		sortWeight1(const sortWeight1 &obj)
		{
			array_length = obj.array_length;
			if (index != NULL)
			{
				index = new int[array_length];
				a_1 = obj.a_1;
			}
		}
		//! Initializes sorting argument function
		/*! \param A_1 1D array to be sorted
	     	    \param N length of array */
		void initSort(int *A_1,int N)
		{
			array_length = N;
			
			a_1 = A_1;
			
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
		//! Run sorting on array A, set with 'initSort()'
		/*! Does not alter array A 
		   \return a pointer to a list of indices such that A(I) = B is sorted */
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
		bool my1Sort(int i, int j) 
		{
			return a_1[i] <= a_1[j];
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
				if (my1Sort(list[i],pivot))
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
	
		//! The external array A that the index set I should be generated from
		int *a_1;

		//! The index set such that A(I) is sorted
		int *index;
		
		//! Length of arrays
		int array_length;
		
};

#endif
