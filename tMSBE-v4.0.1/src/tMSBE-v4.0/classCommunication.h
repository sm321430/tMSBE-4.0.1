
#ifndef __CLASS_COMM_H_INCLUDED__
#define __CLASS_COMM_H_INCLUDED__

#include <cstdio>
#include <iostream>
#include "fileIO.cpp"

//! Enable simple communication between external MatLab program through filesystem
/*! This is a class that enables the program to take input through the modification of a spesific file in the run/ folder. \n
    The external program and this running program uses the filesystem to exchange commands. \n
    This program was used to allow a genetic algorithm running in MatLab to evaluate a function based on output from this program when running on two seperate computers.\n
    In order to ease communication, mirroring of the filesystem is very helpfull.\n
    2018: Additional comments\n
    Isak Kilen
*/

class MatlabCommunication {

	public:
		//! A constructor
		MatlabCommunication()
		{
			stack = NULL;
			NumVariables = 0;
			
			comm_initialize_folder_structure();
			comm_wait_for_handshake();
		}

		//! A destructor
		~MatlabCommunication()
		{
			if (stack != NULL)
			{
				delete [] stack;
			}
		}


		//! A constructor
		/*! \param n is the number of elements to take as input */
		MatlabCommunication(int n)
		{
			if (n<=0)
			{
				std::cout << "MatlabCommunication():: Initializing with <=0 elements not allowed" << std::endl;
				exit(-1);
			}

			stack = new double[n];
			NumVariables = n;
			
			comm_initialize_folder_structure();
			comm_wait_for_handshake();
		}
		
		
		//! Check if Matlab wants to terminate the program
		/*! \return true if should terminate else false
		*/
		bool comm_check_for_termination()
		{
			if (fileExists("optimization/comm_terminate_program.txt"))
			{
				return true;	
			}
			return false;
		}

		//! Check if Matlab wants to evaluate the current input
		/*! \return true if we should evaluate else false
		*/
		bool comm_check_for_input()
		{
			if (fileExists("optimization/comm_evaluate_input.txt"))
			{
				return true;
			}
			return false;
		}

		//! Collect input from input array file
		double *comm_get_input()
		{
			// Check for errors
			if (!dirExists("optimization"))
			{
				dirMake("optimization");
			}
			if (!fileExists("optimization/comm_input.txt"))
			{
				cout << "cannot find: optimization/comm_input.txt" << endl;
				exit(-1);
			}

			// Program is NOT free, and results are NOT pending
			if (fileExists("optimization/comm_free.txt"))
			{
				remove("optimization/comm_free.txt");
			}
			if (fileExists("optimization/comm_results_pending.txt"))
			{
				remove("optimization/comm_results_pending.txt");
			}

			// Start importing
			std::ifstream myfile("optimization/comm_input.txt");

			//
			// Jump to last line of input file
			//
			myfile.seekg(-1,ios_base::end);
			char ch;
			myfile.get(ch);
			if (ch == '\n')
			{
				myfile.seekg(-2,ios_base::end);
			}
			
			bool keepLooping = true;
			while(keepLooping)
			{
				myfile.get(ch);
				
				if((int)myfile.tellg() <= 1) {             // If the data was at or before the 0th byte
					myfile.seekg(0);                       // The first line is the last line
					keepLooping = false;                // So stop there
				}
				else if(ch == '\n') {                   // If the data was a newline
					keepLooping = false;                // Stop at the current position.
				}
				else {                                  // If the data was neither a newline nor at the 0 byte
					myfile.seekg(-2,ios_base::cur);        // Move to the front of that data, then to the front of the data before it
				}
			}

			// Read, will ignore seperators: Whitespace, tab,..
			for(int i = 0; i < NumVariables; i++)
			{
				myfile >> stack[i];
			}
			myfile.close();


			return stack;
		}

		//! Collect input density from input file
		double comm_get_input_density()
		{
			// Check for errors
			if (!dirExists("optimization"))
			{
				dirMake("optimization");
			}
			if (!fileExists("optimization/comm_input_density.txt"))
			{
				cout << "cannot find: optimization/comm_input_density.txt" << endl;
				exit(-1);
			}

			// Program is NOT free, and results are NOT pending
			/*
			if (fileExists("optimization/comm_free.txt"))
			{
				remove("optimization/comm_free.txt");
			}
			if (fileExists("optimization/comm_results_pending.txt"))
			{
				remove("optimization/comm_results_pending.txt");
			}
			*/

			// Start importing
			std::ifstream myfile("optimization/comm_input_density.txt");

			//
			// Jump to last line of input file
			//
			myfile.seekg(-1,ios_base::end);
			char ch;
			myfile.get(ch);
			if (ch == '\n')
			{
				myfile.seekg(-2,ios_base::end);
			}
			
			bool keepLooping = true;
			while(keepLooping)
			{
				myfile.get(ch);
				
				if((int)myfile.tellg() <= 1) {             // If the data was at or before the 0th byte
					myfile.seekg(0);                       // The first line is the last line
					keepLooping = false;                // So stop there
				}
				else if(ch == '\n') {                   // If the data was a newline
					keepLooping = false;                // Stop at the current position.
				}
				else {                                  // If the data was neither a newline nor at the 0 byte
					myfile.seekg(-2,ios_base::cur);        // Move to the front of that data, then to the front of the data before it
				}
			}

			// Read, will ignore seperators: Whitespace, tab,..
			double new_dens = 0;
			myfile >> new_dens;
			myfile.close();
			return new_dens;
		}
		
		//! Set up the flags for completion and write to the output file
		/*! This tells Matlab that the computation is done and it is safe to collect the output\n
		    File: comm_output.txt
		*/
		void comm_write_output(double res)
		{
			// Check for errors
			if (!dirExists("optimization"))
			{
				dirMake("optimization");
			}
			if (!fileExists("optimization/comm_output.txt"))
			{
				cout << "cannot find: optimization/comm_output.txt" << endl;
				exit(-1);
			}

			FILE *fid = fopen("optimization/comm_output.txt","a+");
			fprintf(fid,"%.16e\n",res);
			fclose(fid);
			
			// Tell Matlab that we are free, results are pending and remove the evaluation flag
			fclose(fopen("optimization/comm_results_pending.txt","a+"));
			fclose(fopen("optimization/comm_free.txt","a+"));
			remove("optimization/comm_evaluate_input.txt");
		}
		
	private:
		//! Storage for input array
		double *stack;
		//! Length of 'stack' array
		int NumVariables;


		//! Create folder structure if it does not exist
		void comm_initialize_folder_structure()
		{
			if (dirExists("optimization"))
			{
				// Clear old optimization data
				remove("optimization");
			}

			dirMake("optimization");
			
			// Create empty files
			fclose(fopen("optimization/comm_input.txt","a+"));
			fclose(fopen("optimization/comm_input_density.txt","a+"));
			fclose(fopen("optimization/comm_output.txt","a+"));
			fclose(fopen("optimization/comm_free.txt","a+"));
			
		}
		
		//! Wait for file comm_hand_extend.txt to exist, then create file comm_hand_shake.txt
		void comm_wait_for_handshake()
		{
			while(1==1)
			{
				if (fileExists("optimization/comm_hand_extend.txt"))
				{
					fclose(fopen("optimization/comm_hand_shake.txt","a+"));
					break;
				} else {
					sleep(1);
				}
			}
		}


};

#endif


