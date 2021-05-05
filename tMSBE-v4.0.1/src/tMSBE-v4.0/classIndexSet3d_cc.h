/*
	To store a variable indexd by carrier scattering
	* 
	Isak Kilen @ 2018
*/

#ifndef __INDEXSET3D_CC_H_INCLUDED__
#define __INDEXSET3D_CC_H_INCLUDED__


#include <iostream>
#include "classSortWeight3.h"

using namespace std;


class indexSet3d_cc
{
	public:
		indexSet3d_cc()
		{
			weight	= NULL; 
			i_kp 	= NULL; 
			i_k_q 	= NULL; 
			i_kp_q 	= NULL; 
			b_kp 	= NULL; 
			b_k_q 	= NULL; 
			b_kp_q 	= NULL; 
			index_length_row = 0;
			index_length_col = 0;
		
			index_length_K = NULL;
			
			indexSort = NULL;
		}
		
		~indexSet3d_cc()
		{
			if (weight != NULL)
			{
				for(int i = 0; i < index_length_row; i++)
				{
					for(int j = 0; j < index_length_col; j++)
					{
						delete [] weight[i][j];
						delete [] i_kp[i][j];
						delete [] i_k_q[i][j];
						delete [] i_kp_q[i][j];
						delete [] b_kp[i][j];
						delete [] b_k_q[i][j];
						delete [] b_kp_q[i][j];
					}
					delete [] weight[i];
					delete [] i_kp[i];
					delete [] i_k_q[i];
					delete [] i_kp_q[i];
					delete [] b_kp[i];
					delete [] b_k_q[i];
					delete [] b_kp_q[i];

					delete [] index_length_K[i];
				}
				delete [] weight;
				delete [] i_kp;
				delete [] i_k_q;
				delete [] i_kp_q;
				delete [] b_kp;
				delete [] b_k_q;
				delete [] b_kp_q;

				delete [] index_length_K;

				delete indexSort;
			}
		}

		indexSet3d_cc(const indexSet3d_cc &obj)
		{
			index_length_row = obj.index_length_row;
			index_length_col = obj.index_length_col;

			if (obj.weight!= NULL)
			{
				index_length_K = new int*[index_length_row];
				for(int i = 0; i < index_length_row; i++)
				{
					index_length_K[i] = new int[index_length_col];
					for(int j = 0; j < index_length_col; j++)
					{
						index_length_K[i][j] = obj.index_length_K[i][j];
					}
				}
				weight	= new double**[index_length_row];
				i_kp 	= new int**[index_length_row]; 
				i_k_q 	= new int**[index_length_row]; 
				i_kp_q 	= new int**[index_length_row]; 
				b_kp 	= new double**[index_length_row]; 
				b_k_q 	= new double**[index_length_row]; 
				b_kp_q 	= new double**[index_length_row]; 

				for(int i = 0; i < index_length_row; i++)
				{
					weight[i]	= new double*[index_length_col];
					i_kp[i]		= new int*[index_length_col];
					i_k_q[i]	= new int*[index_length_col];
					i_kp_q[i]	= new int*[index_length_col];
					b_kp[i]		= new double*[index_length_col];
					b_k_q[i]	= new double*[index_length_col];
					b_kp_q[i]	= new double*[index_length_col];
					
					for(int j = 0; j < index_length_col; j++)
					{
						weight[i][j] 	= new double[index_length_K[i][j]];
						i_kp[i][j]	= new int[index_length_K[i][j]];
						i_k_q[i][j]	= new int[index_length_K[i][j]];
						i_kp_q[i][j]	= new int[index_length_K[i][j]];
						b_kp[i][j]	= new double[index_length_K[i][j]];
						b_k_q[i][j]	= new double[index_length_K[i][j]];
						b_kp_q[i][j]	= new double[index_length_K[i][j]];

						for(int k = 0; k < index_length_K[i][j]; k++)
						{
							weight[i][j][k] = obj.weight[i][j][k];
							i_kp[i][j][k]   = obj.i_kp[i][j][k];
							i_k_q[i][j][k]  = obj.i_k_q[i][j][k];
							i_kp_q[i][j][k] = obj.i_kp_q[i][j][k];
							b_kp[i][j][k]   = obj.b_kp[i][j][k];
							b_k_q[i][j][k]  = obj.b_k_q[i][j][k];
							b_kp_q[i][j][k] = obj.b_kp_q[i][j][k];
						}
					}
				}


				indexSort = new sortWeight3(*obj.indexSort);
			}
		}
		
		indexSet3d_cc(int N, int M, int **K)
		{
			index_length_row = N;
			index_length_col = M;
			
			index_length_K = new int*[index_length_row];
			for(int i = 0; i < index_length_row; i++)
			{
				index_length_K[i] = new int[index_length_col];
				for(int j = 0; j < index_length_col; j++)
				{
					index_length_K[i][j] = K[i][j];
				}
			}
			
			weight	= new double**[index_length_row];
			i_kp 	= new int**[index_length_row]; 
			i_k_q 	= new int**[index_length_row]; 
			i_kp_q 	= new int**[index_length_row]; 
			b_kp 	= new double**[index_length_row]; 
			b_k_q 	= new double**[index_length_row]; 
			b_kp_q 	= new double**[index_length_row]; 

			indexSort = new sortWeight3();
			
			for(int i = 0; i < index_length_row; i++)
			{
				weight[i]	= new double*[index_length_col];
				i_kp[i]		= new int*[index_length_col];
				i_k_q[i]	= new int*[index_length_col];
				i_kp_q[i]	= new int*[index_length_col];
				b_kp[i]		= new double*[index_length_col];
				b_k_q[i]	= new double*[index_length_col];
				b_kp_q[i]	= new double*[index_length_col];
				
				for(int j = 0; j < index_length_col; j++)
				{
					weight[i][j] 	= new double[index_length_K[i][j]];
					i_kp[i][j]	= new int[index_length_K[i][j]];
					i_k_q[i][j]	= new int[index_length_K[i][j]];
					i_kp_q[i][j]	= new int[index_length_K[i][j]];
					b_kp[i][j]	= new double[index_length_K[i][j]];
					b_k_q[i][j]	= new double[index_length_K[i][j]];
					b_kp_q[i][j]	= new double[index_length_K[i][j]];
				}
			}
		}
		
		void updateWeight(int k, int kp, int index, double F, int ind_kp, double beta_kp, int ind_k_q, double beta_k_q, int ind_kp_q, double beta_kp_q)
		{
			weight[k][kp][index] 	= F;
			i_kp[k][kp][index]	= ind_kp;
			i_k_q[k][kp][index]	= ind_k_q;
			i_kp_q[k][kp][index]	= ind_kp_q;
			b_kp[k][kp][index] 	= beta_kp;
			b_k_q[k][kp][index] 	= beta_k_q;
			b_kp_q[k][kp][index] 	= beta_kp_q;
		}

		void sortIndices()
		{
			// Initialize temporary arrays
			double ***weight_dummy	= new double**[index_length_row];
			int ***i_dummy 		= new int**[index_length_row]; 
			double ***b_dummy 	= new double**[index_length_row]; 

			for(int i = 0; i < index_length_row; i++)
			{
				weight_dummy[i]	= new double*[index_length_col];
				i_dummy[i]	= new int*[index_length_col];
				b_dummy[i]	= new double*[index_length_col];
				
				for(int j = 0; j < index_length_col; j++)
				{
					weight_dummy[i][j]	= new double[index_length_K[i][j]];
					i_dummy[i][j]		= new int[index_length_K[i][j]];
					b_dummy[i][j]		= new double[index_length_K[i][j]];
				}
			}
			for(int i = 0; i < index_length_row; i++)
			{
				for(int j = 0; j < index_length_col; j++)
				{
					// Sort indices, does not change arrays
					indexSort->initSort(i_kp[i][j], i_k_q[i][j], i_kp_q[i][j], index_length_K[i][j]);
					int *index_list = indexSort->runSort();

					// Kp
					for(int k = 0; k < index_length_K[i][j]; k++)
					{
						weight_dummy[i][j][k]	= weight[i][j][index_list[k]];
						i_dummy[i][j][k]	= i_kp[i][j][index_list[k]];
						b_dummy[i][j][k]	= b_kp[i][j][index_list[k]];
					}
					for(int k = 0; k < index_length_K[i][j]; k++)
					{
						weight[i][j][k]	= weight_dummy[i][j][k];
						i_kp[i][j][k]	= i_dummy[i][j][k];
						b_kp[i][j][k]	= b_dummy[i][j][k];
					}

					// k-q
					for(int k = 0; k < index_length_K[i][j]; k++)
					{
						i_dummy[i][j][k]	= i_k_q[i][j][index_list[k]];
						b_dummy[i][j][k]	= b_k_q[i][j][index_list[k]];
					}
					for(int k = 0; k < index_length_K[i][j]; k++)
					{
						i_k_q[i][j][k]	= i_dummy[i][j][k];
						b_k_q[i][j][k]	= b_dummy[i][j][k];
					}

					// k'-q
					for(int k = 0; k < index_length_K[i][j]; k++)
					{
						i_dummy[i][j][k]	= i_kp_q[i][j][index_list[k]];
						b_dummy[i][j][k]	= b_kp_q[i][j][index_list[k]];
					}
					for(int k = 0; k < index_length_K[i][j]; k++)
					{
						i_kp_q[i][j][k]	= i_dummy[i][j][k];
						b_kp_q[i][j][k]	= b_dummy[i][j][k];
					}
				}
			}

			// Delete dummy arrays
			for(int i = 0; i < index_length_row; i++)
			{
				for(int j = 0; j < index_length_col; j++)
				{
					delete [] weight_dummy[i][j];
					delete [] i_dummy[i][j];
					delete [] b_dummy[i][j];
				}

				delete [] weight_dummy[i];
				delete [] i_dummy[i];
				delete [] b_dummy[i];
			}
			delete [] weight_dummy;
			delete [] i_dummy;
			delete [] b_dummy;
			
		}
		
		inline double F(int i, int j, int index)
		{
			return weight[i][j][index];
		}

		inline int I_kp(int i, int j, int index)
		{
			return i_kp[i][j][index];
		}
		
		inline int I_kp_q(int i, int j, int index)
		{
			return i_kp_q[i][j][index];
		}

		inline int I_k_q(int i, int j, int index)
		{
			return i_k_q[i][j][index];
		}

		inline double beta_kp(int i, int j, int index)
		{
			return b_kp[i][j][index];
		}

		inline double beta_k_q(int i, int j, int index)
		{
			return b_k_q[i][j][index];
		}

		inline double beta_kp_q(int i, int j, int index)
		{
			return b_kp_q[i][j][index];
		}
		
		inline int num_indices(int i, int j)
		{
			return index_length_K[i][j];
		}
	
	private:
	
		double 	***weight; // Primary variable
		int 	***i_kp; // k-q or k+q
		int	***i_k_q; // k-q or k+q
		int 	***i_kp_q; // k-q or k+q
		double ***b_kp; // Interpolation constant
		double ***b_k_q; // Interpolation constant
		double ***b_kp_q; // Interpolation constant
		int index_length_row; // Number of row elements
		int index_length_col; // Number of row elements
		int **index_length_K; // Number of element in (row,col)
		
		sortWeight3 *indexSort;
};

#endif







