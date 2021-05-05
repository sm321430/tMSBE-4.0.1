
/*
	The LARGEVOLUME class is a storage system used for fast easy storage of
	complex elements
	Isak Kilen @ 2016
*/

#ifndef __LARGEVOLSTORE_H_INCLUDED__
#define __LARGEVOLSTORE_H_INCLUDED__


#include <stdio.h>
#include <iostream>
#include <complex>
#include <cstring> // memcpy
#include "fileIO.cpp"
#include <sstream>
#include <fstream>


using namespace std;

//! A class to hold a very large amount of data
/*! This class is designed to split up datasets that are too big to store in a single block of memory.\n
    In particular, it maintains the cyclic nature of the 'Cyclic()' class storage.\n
    \sa Cyclic \n
    2018: Added more comments \n
    Isak Kilen
*/
class largeVolumeStorage{

	public:
		//! A Constructor
		largeVolumeStorage()
		{
			MAX_BUCKET_SIZE = 100000;
			head = NULL;
			iter_i = 0;
			iter_j = 0;
			
		}
		
		//! A Constructor
		/*! \param num The size of array to store */
		largeVolumeStorage(int num)
		{
			MAX_BUCKET_SIZE = 100000;
			iter_i = 0;
			iter_j = 0;
			
			number_of_buckets = floor(num/MAX_BUCKET_SIZE);

			if (num % MAX_BUCKET_SIZE > 0)
			{
				number_of_buckets += 1;
			}

			int counter = num;
			head = new container[number_of_buckets];
			for(int i = 0; i < number_of_buckets; i++)
			{
				if (counter > MAX_BUCKET_SIZE)
				{
					head[i].number_of_elements = MAX_BUCKET_SIZE;
				} else {
					head[i].number_of_elements = counter;
				}
				head[i].data = new std::complex<double>[head[i].number_of_elements];
				
				counter -= MAX_BUCKET_SIZE;
			}

			zero_all_containers();
		}
	
		//! A destructor
		~largeVolumeStorage()
		{
			for(int i  =0; i < number_of_buckets; i++)
			{
				delete [] head[i].data;
			}
			delete [] head;

		}
		//! Insert element into beginning of array and move head to next position
		void seq_set_iterate(std::complex<double> el)
		{
			head[iter_i].data[iter_j] = el;
			seq_iterate();
		}
		//! Return element at head position and move head to next position
		std::complex<double> seq_get_iterate()
		{
			std::complex<double> el = head[iter_i].data[iter_j];
			seq_iterate();
			return el;
		}

		//! Reset counters for head position
		void seq_iterate_reset()
		{
			iter_i = 0;
			iter_j = 0;
		}

		//! Increment position of head in array
		void seq_iterate()
		{
			iter_j += 1;
			if (iter_j == head[iter_i].number_of_elements)
			{
				iter_j = 0;

				iter_i += 1;
				iter_i = iter_i % number_of_buckets;
			}
		}

		//! Set element at a given position
		/*! \param num position in array\n
		    \param el element to insert*/
		void set_element(int num, std::complex<double> el)
		{
			int i  = floor(num/MAX_BUCKET_SIZE);
			int j  = num % MAX_BUCKET_SIZE;
			head[i].data[j] = el;
		}

		//! Return element at a given positin in array
		std::complex<double> get_element(int num)
		{
			int i  = floor(num/MAX_BUCKET_SIZE);
			int j  = num % MAX_BUCKET_SIZE;

			return head[i].data[j];
		}
		
		//! Set all containers to zero
		void zero_all_containers()
		{
			for(int i  =0; i < number_of_buckets; i++)
			{
				for(int j = 0; j < head[i].number_of_elements; j++)
				{
					head[i].data[j] = 0.0;
				}
			}
		}

		//! Print information about stored array to screen
		void Print()
		{
			cout << "Large volume container" << endl;
			cout << "# of buckets    = " << number_of_buckets << endl;
			cout << "MAX bucket size = " << MAX_BUCKET_SIZE << endl;
			for(int i  =0; i < number_of_buckets; i++)
			{
				cout << "bucket[" << i << "]" << endl;
				for(int j = 0; j < head[i].number_of_elements; j++)
				{
					cout << " i =  " << i*MAX_BUCKET_SIZE + j << " , " << head[i].data[j] << endl;
				}
				cout << endl;
			}
		}

		//! Write stored array to binary file with given filename prepended
		/*! Also writes information about head node position to file*/
		void saveArray(const std::string &fileName)
		{
			std::stringstream contName;
			for(int i = 0; i  < number_of_buckets; i++)
			{
				contName.str("");
				contName << fileName << "_container_" << i << ".dat";
				saveBinary(contName.str(), head[i].data, head[i].number_of_elements);
			}

			contName.str("");
			contName << fileName << "_iterI.dat";
			saveBinary(contName.str(), &iter_i, 1);

			contName.str("");
			contName << fileName << "_iterJ.dat";
			saveBinary(contName.str(), &iter_j, 1);
			
		}

		//! Load array from binary file with given filename into class
		/*! No information about the number of elements is stored, the class must have enough storage to recover from the file.\n
		    Also loads information about position of head node */
		void loadArray(const std::string &fileName)
		{
			std::stringstream contName;
			for(int i = 0; i  < number_of_buckets; i++)
			{
				// Edit filename to match container
				contName.str("");
				contName << fileName << "_container_" << i << ".dat";
				loadBinary(contName.str(), head[i].data, head[i].number_of_elements);
			}

			contName.str("");
			contName << fileName << "_iterI.dat";
			loadBinary(contName.str(), &iter_i, 1);

			contName.str("");
			contName << fileName << "_iterJ.dat";
			loadBinary(contName.str(), &iter_j, 1);
		}

	private:
		//! The data storage unit for each segment of an array
		struct container {
			int number_of_elements;
			std::complex<double> *data;
		};
		
		//! Compile time set limit for largest array size
		int MAX_BUCKET_SIZE;
		//! Number of buckets used to store data
		int number_of_buckets;
		//! Pointre to head container
		struct container *head;
		
		//! Bucket counter for first element
		int iter_i;
		//! Array counter for first element
		int iter_j;
};


#endif
