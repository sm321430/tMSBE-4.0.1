#ifndef __PARSCHEDULE_H_INCLUDED__
#define __PARSCHEDULE_H_INCLUDED__

//! Create a balanced workload for a parallel schedulers
/*! Take an array of run-times (A) from multiple tasks and the number of parallel nodes (N(, and re-group the tasks into N sub-groups that each have a run-time that is very similar (~avg(A)).\n
    The result is an index array I and a vector T such that A(I) = B is a sorted array where T is the number of elements in each parallel grouping such that sum(T) = length(A)\n
2019: First version\n
Isak Kilen
*/

#include <cstring> //memset
#include <iostream>
#include <stdio.h>
#include <cstdlib>


using namespace std;

class parSchedule
{

	public:
		//! A constrcutor
		parSchedule()
		{
			number_of_elements = -1;
			number_of_sub_groups = -1;
			run_times = NULL;
			run_times_index = NULL;
			sub_group_info = NULL;
			sub_group_index = NULL;
			run_time_0_offset = 0.0;
			run_time_0_scale = 1.0;

		}

		//! A destructor
		~parSchedule()
		{
			if (run_times != NULL)
			{
				delete [] run_times;
				delete [] run_times_index;
				delete [] sub_group_info;
				for(int i = 0; i < number_of_sub_groups; i++)
				{
					delete [] sub_group_index[i];
				}
				delete [] sub_group_index;
			}
		}

		//! A copy-constructor
		parSchedule(const parSchedule &obj)
		{
			if (run_times != NULL)
			{
				run_time_0_offset = obj.run_time_0_offset;
				run_time_0_scale  = obj.run_time_0_scale;

				number_of_elements = obj.number_of_elements;
				number_of_sub_groups = obj.number_of_sub_groups;
	
				for(int i = 0; i < number_of_elements; i++)
				{
					run_times[i] = obj.run_times[i];
				}
				for(int i = 0; i < number_of_sub_groups; i++)
				{
					run_times_index[i] = obj.run_times_index[i];
					sub_group_info[2*i  ] = obj.sub_group_info[2*i  ];
					sub_group_info[2*i+1] = obj.sub_group_info[2*i+1];

					for(int j = 0; j < sub_group_info[2*i]; j++)
					{
						sub_group_index[i][j] = obj.sub_group_index[i][j];
					}
				}
			}
		}
		
		//! Run optimization on dataset
		void optimize_schedule(int length, double *data, int sub_groups, double group_0_offset, double group_0_scale)
		{
			run_time_0_offset = group_0_offset;
			run_time_0_scale  = group_0_scale;
			number_of_sub_groups = sub_groups;

			number_of_elements = length;
			run_times = new double[number_of_elements];
			run_times_index = new int[number_of_elements];

			if (number_of_sub_groups <= 0)
			{
				cout << "parSchedule::init_scheduler() Error, cannot divide data into " << number_of_sub_groups << " sub_groups.." << endl;
				exit(-1);
			} else if ( number_of_sub_groups == 0){

				double sum = 0.0;
				for(int i = 0; i < number_of_elements; i++)
				{
					run_times_index[i] = i;
					if(data[i] <= 0.0)
					{
						cout << "parSchedule::init_scheduler() Error, data array should only contain positive (data[i]>0) numbers " << endl;
						exit(-1);
					}
					run_times[i] = data[i];
					sum += data[i];
				}
				sub_group_info = new double[2*number_of_sub_groups];
				sub_group_info[0] = number_of_elements;
				sub_group_info[1] = sum;
				
				return;
							
			} else if (number_of_sub_groups > length)
			{
				cout << "parSchedule::init_scheduler() Error, cannot divide " << length << " into " << number_of_sub_groups << " sub_groups.." << endl;
				exit(-1);
			}
			sub_group_info = new double[2*number_of_sub_groups];
			for(int i = 0; i < 2*number_of_sub_groups; i++)
			{
				sub_group_info[i] = 0.0;
			}

			for(int i = 0; i < number_of_elements; i++)
			{
				run_times_index[i] = i;
				if(data[i] <= 0.0)
				{
					cout << "parSchedule::init_scheduler() Error, data array should only contain positive (data[i]>0) numbers " << endl;
					exit(-1);
				}
				run_times[i] = data[i];
			}


			quick_sort(run_times_index, number_of_elements);
			double dummy[number_of_elements];
			for(int i = 0; i < number_of_elements; i++)
			{
				dummy[i] = run_times[run_times_index[i]];
			}
			memcpy(run_times, dummy, number_of_elements*sizeof(double));

			
			balance_array_info(0, number_of_sub_groups);
/*
			cout << "sub_group_info:" << endl;
			int num_tot = 0;
			double sum_tot = 0.0;
			for(int i = 0; i < number_of_sub_groups; i++)
			{
				cout << "       N = " << sub_group_info[2*i] << ", sum = " << sub_group_info[2*i+1] << endl;
				num_tot += sub_group_info[2*i];
				sum_tot += sub_group_info[2*i+1];
			}
			cout << "Total: N = " << num_tot << ", sum = " << sum_tot << ", ideal = " << sum_tot/((double)number_of_sub_groups) << endl;

*/
			// 3. Allocate balanced array structure
			sub_group_index = new int*[number_of_sub_groups];
			for(int i = 0; i < number_of_sub_groups; i++)
			{
				sub_group_index[i] = new int[(int)sub_group_info[2*i]];
				memset(sub_group_index[i],  -1, sub_group_info[2*i]*sizeof(int));
			}

			// 4. Fill balanced array
			for(int i = 0; i < number_of_sub_groups; i++)
			{
				sub_group_info[2*i+1] = 0.0; // reset sums
			}

			balance_array(0, number_of_sub_groups);

			// Sort index groups
			for(int i = 0; i < number_of_sub_groups; i++)
			{
				quick_sort2(sub_group_index[i], sub_group_info[2*i]);
				
			}
/*
			for(int i = 0; i < number_of_sub_groups; i++)
			{
				
				for(int j = 0; j < sub_group_info[2*i]; j++)
				{
					cout << sub_group_index[i][j] << " ";
				}
				cout << endl;
			}
*/
			
			
			// Structured output
			int cnt = 0;
			for(int i = 0; i < number_of_sub_groups; i++)
			{
				for(int j = 0; j < sub_group_info[2*i]; j++)
				{
					run_times_index[cnt] = sub_group_index[i][j];
					run_times[cnt] = data[run_times_index[cnt]];
					cnt++;
				}
			}
/*
			cout << "new array:" << endl;
			for(int i = 0; i < number_of_elements; i++)
			{
				cout << run_times[i] << " ";
			}
			cout << endl;
			cout << "new array index:" << endl;
			for(int i = 0; i < number_of_elements; i++)
			{
				printf("%2d ",run_times_index[i]);
			}
			cout << endl;

			cout << "grouping:" << endl;
			for(int i = 0; i < number_of_sub_groups; i++)
			{
				cout << sub_group_info[2*i] << " ";
			}
			cout << endl;
			cout << "done.." << endl;
*/

			
		}

		//! Copy number of elements in each subgroup
		void get_sub_group_info(int * info)
		{
			for(int i = 0; i < number_of_sub_groups; i++)
			{
				info[i] = sub_group_info[2*i];
			}
		}

		//! Copy index set I such that A(I) = B is balanced
		void get_sub_group_index(int * info)
		{
			for(int i = 0; i < number_of_elements; i++)
			{
				info[i] = run_times_index[i];
			}
		}

		//! Save results to file
		/*! Format is:\n
			Line 1: N = #sub-groups
			Line 2->N+1: #size-of-each-balanced-subgroup, #weight-in-this-subgroup \n
			Line N+2: M = #total-elements
			Line N+3: #index-of-each-element-balanced\n
		*/
		void file_write_structure(const std::string & newName)
		{
			FILE *fid = fopen(newName.c_str(), "w+");
			fprintf(fid,"%d\n",number_of_sub_groups);
			for(int i = 0; i < number_of_sub_groups; i++)
			{
				fprintf(fid,"%d %.3f\n",(int)sub_group_info[2*i], sub_group_info[2*i+1]);
			}
			fprintf(fid, "\n");

			fprintf(fid,"%d\n", number_of_elements);
			for(int i = 0; i < number_of_elements; i++)
			{
				fprintf(fid,"%d ",run_times_index[i]);
			}
			fprintf(fid, "\n");

			fprintf(fid,"\n");
			for(int i = 0; i < number_of_elements; i++)
			{
				fprintf(fid,"[%4d]: %.6f\n",i,run_times[i]);
			}
			fprintf(fid, "\n");
			
			fclose(fid);
		}
		
	private:

		int number_of_elements;
		int number_of_sub_groups;
		double *run_times;
		int *run_times_index;
		double *sub_group_info;
		int **sub_group_index;
		double run_time_0_offset; // A potential offset to ONLY the first one
		double run_time_0_scale; // A potential scaling of elements added to first group

		
		//! Quick sort helper function
		int partition(int list[], int p, int r)
		{
			int index, exchange_temp;
			int pivot = list[r];
			index = p - 1;
			for(int i = p; i < r; i++)
			{
				if (run_times[list[i]] >= run_times[pivot])
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

		//! Quick sort helper function
		int partition2(int list[], int p, int r)
		{
			int index, exchange_temp;
			int pivot = list[r];
			index = p - 1;
			for(int i = p; i < r; i++)
			{
				if (list[i] < pivot)
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
		void quicksort_aux2(int list[], int p, int r)
		{
			int q;
			if(p<r)
			{
				q = partition2(list, p, r);
				quicksort_aux2(list, p, q-1);
				quicksort_aux2(list, q+1, r);
			}
		}

		void quick_sort2(int list[], int size)
		{
			quicksort_aux2(list,0, size-1);
		}

		void balance_array_info(int index0, int num_sub_array)
		{
			if (num_sub_array > 1)
			{
				double S = 0.0;
				for(int i = index0; i < number_of_elements; i++)
				{
					S += run_times[i];
				}
				if (index0 == 0)
				{
					// Add offset to first group only
					sub_group_info[1] = run_time_0_offset;
					S += run_time_0_offset;
				}
				S = S/((double)num_sub_array);

				int ind1 = index0;
				for(int i = 0; i < num_sub_array; i++)
				{
					// Find sub-array with the lowest score
					double min_score = sub_group_info[1];
					int target_sub_array = 0;
					for(int j = 1; j < number_of_sub_groups; j++)
					{
						if (sub_group_info[2*j+1] < min_score)
						{
							min_score   = sub_group_info[2*j+1];
							target_sub_array = j;
						}
					}
					
					// Populate sub-array
					double tmp_sum = 0.0;
					if (index0 == 0)
					{
						tmp_sum = sub_group_info[2*target_sub_array+1];
					}
					// Scale element of first group
					double potential_new_element = run_times[ind1];
					if (target_sub_array == 0)
					{
						potential_new_element *= run_time_0_scale;
					}
					while (true)
					{
						if (tmp_sum + potential_new_element < S)
						{
							tmp_sum += potential_new_element;
							sub_group_info[2*target_sub_array  ] += 1.0;
							sub_group_info[2*target_sub_array+1] += potential_new_element;
							ind1 += 1;
						} else {
							break;
						}
					}
				}
				balance_array_info(ind1, num_sub_array-1);

			} else {

				for(int k = index0; k < number_of_elements; k++)
				{
					// Find sub-array with the lowest score
					double min_score = sub_group_info[1];
					int target_sub_array = 0;
					for(int j = 1; j < number_of_sub_groups; j++)
					{
						if (sub_group_info[2*j+1] < min_score)
						{
							min_score   = sub_group_info[2*j+1];
							target_sub_array = j;
						}
					}

					// Scale element of first group
					double potential_new_element = run_times[k];
					if (target_sub_array == 0)
					{
						potential_new_element *= run_time_0_scale;
					}

					sub_group_info[2*target_sub_array  ] += 1;
					sub_group_info[2*target_sub_array+1] += potential_new_element;
				}
			}
		}

		
		void balance_array(int index0, int num_sub_array)
		{
			if (num_sub_array > 1)
			{
				double S = 0.0;
				for(int i = index0; i < number_of_elements; i++)
				{
					S += run_times[i];
				}
				if (index0 == 0)
				{
					// Add offset to first group only
					sub_group_info[1] = run_time_0_offset;
					S += run_time_0_offset;
				}
				S = S/((double)num_sub_array);

				int ind1 = index0;
				for(int i = 0; i < num_sub_array; i++)
				{
					// Find sub-array with the lowest score
					double min_score = sub_group_info[1];
					int target_sub_array = 0;
					for(int j = 1; j < number_of_sub_groups; j++)
					{
						if (sub_group_info[2*j+1] < min_score)
						{
							min_score   = sub_group_info[2*j+1];
							target_sub_array = j;
						}
					}
				
					// Find first non-zero array element
					int ind_sub = 0;
					for(int j = 0; j < sub_group_info[2*target_sub_array]; j++)
					{
						if (sub_group_index[target_sub_array][j] == -1)
						{
							ind_sub = j;
							break;
						}
					}
					
					// Populate sub-array
					double tmp_sum = 0.0;
					if (index0 == 0)
					{
						tmp_sum = sub_group_info[2*target_sub_array+1];
					}
					// Scale element of first group
					double potential_new_element = run_times[ind1];
					if (target_sub_array == 0)
					{
						potential_new_element *= run_time_0_scale;
					}
					while (true)
					{
						if (tmp_sum + potential_new_element < S)
						{
							tmp_sum += potential_new_element;
							sub_group_info[2*target_sub_array+1]  += potential_new_element;
							sub_group_index[target_sub_array][ind_sub] = run_times_index[ind1];
							ind1    += 1;
							ind_sub += 1;

						} else {
							break;
						}
					}
				}

				balance_array(ind1, num_sub_array-1);

			} else {

				for(int k = index0; k < number_of_elements; k++)
				{
					// Find sub-array with the lowest score
					double min_score = sub_group_info[1];
					int target_sub_array = 0;
					for(int j = 1; j < number_of_sub_groups; j++)
					{
						if (sub_group_info[2*j+1] < min_score)
						{
							min_score   = sub_group_info[2*j+1];
							target_sub_array = j;
						}
					}
				
					// Find first non-zero array element
					int ind_sub = 0;
					for(int j = 0; j < sub_group_info[2*target_sub_array]; j++)
					{
						if (sub_group_index[target_sub_array][j] == -1)
						{
							ind_sub = j;
							break;
						}
					}

					// Scale element of first group
					double potential_new_element = run_times[k];
					if (target_sub_array == 0)
					{
						potential_new_element *= run_time_0_scale;
					}

					sub_group_info[2*target_sub_array+1]   	   += potential_new_element;
					sub_group_index[target_sub_array][ind_sub]  = run_times_index[k];
				}
			}
		}
		
};
#endif















