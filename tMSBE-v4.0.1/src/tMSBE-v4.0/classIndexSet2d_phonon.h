/*
	To store a variable indexd by phonon scattering
	* 
	Isak Kilen @ 2017
*/

#ifndef __INDEXSET2D_PHONON_H_INCLUDED__
#define __INDEXSET2D_PHONON_H_INCLUDED__


#include <iostream>
#include "classSortWeight1.h"

using namespace std;


class indexSet2d_phonon
{
	public:
		indexSet2d_phonon()
		{
			weight_e= NULL; 
			i_kq 	= NULL; 
			b_kq 	= NULL; 
			index_length_row = 0;
			index_length_col = NULL;

			indexSort = NULL;
		}
		
		~indexSet2d_phonon()
		{
			if (weight_e != NULL)
			{
				for(int i = 0; i < index_length_row; i++)
				{
					delete [] weight_e[i];
					delete [] i_kq[i];
					delete [] b_kq[i];
				}
				delete [] weight_e;
				delete [] i_kq;
				delete [] b_kq;
				delete [] index_length_col;

				delete indexSort;
			}
		}

		indexSet2d_phonon(const indexSet2d_phonon &obj)
		{
			weight_e= NULL; 
			i_kq 	= NULL; 
			b_kq 	= NULL; 
			index_length_row = obj.index_length_row;
			index_length_col = NULL;

			if (obj.weight_e!= NULL)
			{
				index_length_col = new int[index_length_row];
				for(int i = 0; i < index_length_row; i++)
				{
					index_length_col[i] = obj.index_length_col[i];
				}

				weight_e= new double*[index_length_row];
				i_kq 	= new int*[index_length_row]; 
				b_kq 	= new double*[index_length_row]; 
				
				for(int i = 0; i < index_length_row; i++)
				{
					weight_e[i]	= new double[index_length_col[i]];
					i_kq[i]		= new int[index_length_col[i]];
					b_kq[i]		= new double[index_length_col[i]];
					
					for(int j = 0; j < index_length_col[i]; j++)
					{
						weight_e[i][j] 	= obj.weight_e[i][j];
						i_kq[i][j]		= obj.i_kq[i][j];
						b_kq[i][j]		= obj.b_kq[i][j];
					}
				}
				
				indexSort = new sortWeight1(*obj.indexSort);
			}
			
		}
		
		indexSet2d_phonon(int N, int *array_M)
		{
			index_length_row = N;
			index_length_col = new int[index_length_row];
			for(int i = 0; i < index_length_row; i++)
			{
				index_length_col[i] = array_M[i];
			}
			
			weight_e= new double*[index_length_row];
			i_kq 	= new int*[index_length_row]; 
			b_kq 	= new double*[index_length_row]; 

			indexSort = new sortWeight1();
			
			for(int i = 0; i < index_length_row; i++)
			{
				weight_e[i]	= new double[index_length_col[i]];
				i_kq[i]		= new int[index_length_col[i]];
				b_kq[i]		= new double[index_length_col[i]];
				
				for(int j = 0; j < index_length_col[i]; j++)
				{
					weight_e[i][j] 	= 0;
					i_kq[i][j]		= 0;
					b_kq[i][j]		= 0;
				}
			}
		}
		
		void updateWeight(int index, double Fe, int k, int q, int kq, double B_kq)
		{
			weight_e[k][index] 	= Fe;
			i_kq[k][index]		= kq;
			b_kq[k][index] 		= B_kq;
		}

		void sortIndices()
		{
			// Initialize temporary arrays
			double **weight_dummy	= new double*[index_length_row];
			int **i_dummy 		= new int*[index_length_row]; 
			double **b_dummy 	= new double*[index_length_row]; 

			for(int i = 0; i < index_length_row; i++)
			{
				weight_dummy[i]	= new double[index_length_col[i]];
				i_dummy[i]	= new int[index_length_col[i]];
				b_dummy[i]	= new double[index_length_col[i]];
				
				for(int j = 0; j < index_length_col[i]; j++)
				{
					weight_dummy[i][j]	= 0;
					i_dummy[i][j]		= 0;
					b_dummy[i][j]		= 0;
				}
			}
			for(int i = 0; i < index_length_row; i++)
			{
				// Sort indices, does not change arrays
				indexSort->initSort(i_kq[i], index_length_col[i]);
				int *index_list = indexSort->runSort();

				// Sort Q vectors
				for(int j = 0; j < index_length_col[i]; j++)
				{
					weight_dummy[i][j]	= weight_e[i][index_list[j]];
					i_dummy[i][j]		= i_kq[i][index_list[j]];
					b_dummy[i][j]		= b_kq[i][index_list[j]];
				}
				for(int j = 0; j < index_length_col[i]; j++)
				{
					weight_e[i][j]	= weight_dummy[i][j];
					i_kq[i][j]	= i_dummy[i][j];
					b_kq[i][j]	= b_dummy[i][j];
				}
			}

			// Delete dummy arrays
			for(int i = 0; i < index_length_row; i++)
			{
				delete [] weight_dummy[i];
				delete [] i_dummy[i];
				delete [] b_dummy[i];
			}
			delete [] weight_dummy;
			delete [] i_dummy;
			delete [] b_dummy;
			
		}
		
		inline double F(int i, int index)
		{
			return weight_e[i][index];
		}
		
		inline int I_kq(int i, int index)
		{
			return i_kq[i][index];
		}
		
		inline double beta_kq(int i, int index)
		{
			return b_kq[i][index];
		}
		
		inline int get_num_cols(int i)
		{
			return index_length_col[i];
		}
	
	private:
	
		double **weight_e; // Primary variable
		int **i_kq; // k-q or k+q
		double **b_kq; // Interpolation constant
		int index_length_row; // Number of row elements
		int *index_length_col; // Number of col elements

		sortWeight1 *indexSort;
};

#endif







