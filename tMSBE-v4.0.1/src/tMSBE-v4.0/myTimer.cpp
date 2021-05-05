
#ifndef __MYTIMER_H_INCLUDED__
#define __MYTIMER_H_INCLUDED__

#include <stdio.h>
#include <iostream>
#include <complex>
#include <cstring> // memcpy
#include <sys/time.h>
#include <stdlib.h> // For system()

using namespace std;

//! A timer with additional information. Used in the class 'myTimerCentral()'
/*!
   These timers can be stored/found with a string name search and run-time information of the computaional time and mean/M2-variance of execution time/ number of times used.
   \sa myTimerCentral
*/
class myTimer{

	public:
		//! A constructor
		myTimer()
		{
			my_name = "n/a";
			id = 0;
			counter = 0;
			real_time = 0;
			start_timer = false;
			real_time_mean = 0.0;
			real_time_M2 = 0.0;
		}

		//! A copy-constructor
		myTimer(const std::string & newName)
		{
			my_name = newName;
			id = hash(my_name.c_str());
			counter = 0;
			real_time = 0;
			start_timer = false;
			real_time_mean = 0.0;
			real_time_M2 = 0.0;
		}

		myTimer(const myTimer &obj)
		{
			my_name = obj.my_name;
			id = obj.id;
			counter = obj.counter;
			real_time = obj.real_time;
			start_timer = obj.start_timer;
			real_time_mean = obj.real_time_mean;
			real_time_M2 = obj.real_time_M2;
		}

		//! Return the timer name
		inline std::string getName()
		{
			return my_name;
		}

		//! Return the time stored in timer
		inline double get_real_time()
		{
			return real_time;
		}

		//! Return the number of times the timer has been used
		inline double get_counter()
		{
			return counter;
		}

		//! Return real_time/counter
		inline double get_real_time_mean()
		{
			return real_time_mean;
		}
		
		//! Return sample standard deviation
		/*! real_time = real_time_mean +/- sqrt(s^2)*/
		inline double get_real_time_deviation()
		{
			return sqrt(real_time_M2/(counter-1.0));
		}

		//! Return the ID number of the timer computed from the timer_name
		inline int getId()
		{
			return id;
		}

		//! Start timer using 'clock_gettime()'
		void start()
		{
		//	gettimeofday(&t1,NULL);
			clock_gettime(CLOCK_REALTIME,&t1);
			start_timer = true;
		}

		//! Stop a running timer and store computed time using 'clock_gettime()'
		/*! Also updates the number of times timer has been used, mean time computed, M2-variance*/
		void stop()
		{
			if (start_timer)
			{
			//	gettimeofday(&t2,NULL);
			//	double  diff = ((double)t2.tv_sec + (double)t2.tv_usec * .000001 - (double)t1.tv_sec - (double)t1.tv_usec * .000001);
				clock_gettime(CLOCK_REALTIME,&t2);
				double  diff = ((double)t2.tv_sec + (double)t2.tv_nsec * 1e-9 - (double)t1.tv_sec - (double)t1.tv_nsec * 1e-9);
				real_time += diff;
				

				// Welfords algorithm for running sample variance
				counter += 1;
				double delta = diff - real_time_mean;
				real_time_mean += delta/counter;
				double delta2 = diff - real_time_mean;
				real_time_M2 += delta*delta2;
			} else {
				Print();

				cout << " " << endl;
				cout << "ERROR" << endl;
				cout << "Asking for end() before calling start()" << endl;
				exit(-1);
			}
		}

		//! Print information about timer
		void Print()
		{
			cout << "\n" << endl;
			cout << "================================" << endl;
			cout << "name:         " << my_name << endl;
			cout << "id:           " << id << endl;
			cout << "iterations:   " << counter << endl;
			cout << "total_time:   " << real_time << " [s]" << endl;
			if (real_time > 3600)
			{
				double time_day = real_time/(60.0*60.0*24.0);
				double time_hr  = 24.0*(time_day - floor(time_day));
				double time_min = 60.0*(time_hr  - floor(time_hr));
				double time_sec = 60.0*(time_min - floor(time_min));
				cout << "[d:hr:min:s] " << floor(time_day) << ":" << floor(time_hr) << ":" << floor(time_min) << ":" << floor(time_sec) << endl;
			} else {
				double time_min = real_time/60.0;
				double time_sec = 60.0*(time_min  - floor(time_min));
				double time_ms  = 1000.0*(time_sec - floor(time_sec));
				double time_us  = 1000.0*(time_ms - floor(time_ms));
				cout << "[min:s:ms:us] " << floor(time_min) << ":" << floor(time_sec) << ":" << floor(time_ms) << ":" << floor(time_us) << endl;
			}
		}

		//! Print information about timer with additonal formating
		void Print(int step)
		{
			
			std::stringstream my_step;
			for(int i =0; i < step; i++)
			{
				my_step << "\t";
			}

			cout << my_step.str() << "\n" << endl;
			cout << my_step.str() << "name:         " << my_name << endl;
			cout << my_step.str() << "id:           " << id << endl;
			cout << my_step.str() << "iterations:   " << counter << endl;
			cout << my_step.str() << "total_time:   " << real_time << " [s]" << endl;
			if (real_time > 3600)
			{
				double time_day = real_time/(60.0*60.0*24.0);
				double time_hr  = 24.0*(time_day - floor(time_day));
				double time_min = 60.0*(time_hr  - floor(time_hr));
				double time_sec = 60.0*(time_min - floor(time_min));
				cout << my_step.str() << "[d:hr:min:s] " << floor(time_day) << ":" << floor(time_hr) << ":" << floor(time_min) << ":" << floor(time_sec) << endl;
			} else {
				double time_min = real_time/60.0;
				double time_sec = 60.0*(time_min  - floor(time_min));
				double time_ms  = 1000.0*(time_sec - floor(time_sec));
				double time_us  = 1000.0*(time_ms - floor(time_ms));
				cout << my_step.str() << "[min:s:ms:us] " << floor(time_min) << ":" << floor(time_sec) << ":" << floor(time_ms) << ":" << floor(time_us) << endl;
			}
			

		}
		//! Print information about timer with additonal formating
		void Print_table(int step)
		{
			//---- Print tabs/name
			std::stringstream my_step;
			for(int i =0; i < step; i++)
			{
				my_step << "      ";
			}

			cout << my_step.str() << my_name;

			my_step.str();
			for(int i =0; i < 60-12*step-my_name.length(); i++)
			{
				my_step << " ";
			}
			my_step << ": ";
			cout << my_step.str();

			//---- Print iterations
			cout << counter;
			my_step.str("");
			for(int i =0; i < 10 - GetNumberOfDigits(counter); i++)
			{
				my_step << " ";
			}
			my_step << ": ";
			cout << my_step.str();
			

			//---- Print nice format

			if (real_time > 3600)
			{
				double time_day = real_time/(60.0*60.0*24.0);
				double time_hr  = 24.0*(time_day - floor(time_day));
				double time_min = 60.0*(time_hr  - floor(time_hr));
				double time_sec = 60.0*(time_min - floor(time_min));
				cout << "[d:hr:min:s] " << floor(time_day) << ":" << floor(time_hr) << ":" << floor(time_min) << ":" << floor(time_sec) << " ";
			} else {
				double time_min = real_time/60.0;
				double time_sec = 60.0*(time_min  - floor(time_min));
				double time_ms  = 1000.0*(time_sec - floor(time_sec));
				double time_us  = 1000.0*(time_ms - floor(time_ms));
				cout << "[min:s:ms:us] " << floor(time_min) << ":" << floor(time_sec) << ":" << floor(time_ms) << ":" << floor(time_us) << " ";
			}

			// ---- Print seconds
			cout << "\t: " << real_time << " [s]" << endl;
			

		}

		//! Print information about timer with additonal formating
		void Print_minimal(int step, double parent_real_time, double parent_counter)
		{
			//---- Print tabs/name
			std::stringstream my_step;
			for(int i =0; i < step; i++)
			{
				my_step << "      ";
			}
			my_step << my_name;
			printf("%-60s: ",my_step.str().c_str());

			//---- Print iterations
			printf("%.0f",counter);
			my_step.str("");
			for(int i =0; i < 10 - GetNumberOfDigits(counter); i++)
			{
				my_step << " ";
			}
			my_step << ": ";
			cout << my_step.str();
			

			//---- Print nice format

			if (real_time > 3600)
			{
				double time_day = real_time/(60.0*60.0*24.0);
				double time_hr  = 24.0*(time_day - floor(time_day));
				double time_min = 60.0*(time_hr  - floor(time_hr));
				double time_sec = 60.0*(time_min - floor(time_min));
			//	cout << "[d:hr:min:s] " << floor(time_day) << ":" << floor(time_hr) << ":" << floor(time_min) << ":" << floor(time_sec) << " ";
				printf("[d:hr:min:s] %.0f:%2.0f:%2.0f:%2.0f ", floor(time_day), floor(time_hr), floor(time_min), floor(time_sec));
			} else {
				double time_min = real_time/60.0;
				double time_sec = 60.0*(time_min  - floor(time_min));
				double time_ms  = 1000.0*(time_sec - floor(time_sec));
				double time_us  = 1000.0*(time_ms - floor(time_ms));
				//cout << "[min:s:ms:us] " << floor(time_min) << ":" << floor(time_sec) << ":" << floor(time_ms) << ":" << floor(time_us) << " ";
				printf("[min:s:ms:\xC2\xB5s] %2.0f:%2.0f:%3.0f:%3.0f ",floor(time_min), floor(time_sec), floor(time_ms), floor(time_us));
	
			}

			// ---- Print percentage of parent utilization
			cout << "\t: ";
			if (real_time > 0)
			{
				double tmp = 100.0*(real_time)/parent_real_time;
				for(int i = 0; i < 3-GetNumberOfDigits(floor(tmp)); i++)
				{
					cout << " ";
				}
				printf("%.2f %%",tmp);
			} else {
				cout << "       %";
			}

			// ---- Print seconds
			//cout << "\t: " << real_time << " [s]";
			double tmp = real_time;
			if (tmp < 1e-6)
			{
				printf("\t: %10.2f [ns]",tmp/1.0e-9);
			} else if (tmp < 1e-3)
			{
				printf("\t: %10.2f [\xC2\xB5s]",tmp/1.0e-6);
			} else if (tmp < 1)
			{
				printf("\t: %10.2f [ms]",tmp/1.0e-3);
			} else {
				printf("\t: %10.2f [s]",tmp);
			}

			// ---- Print average time +/- Sample variance s^2
			// Real number [m - sqrt(s^2), m+sqrt(s^2)]
			if (counter > 1)
			{
				if (real_time_mean < 1e-6)
				{
					printf("\t: %6.2f [ns]",real_time_mean/1.0e-9);
				} else if (real_time_mean < 1e-3)
				{
					printf("\t: %6.2f [\xC2\xB5s]",real_time_mean/1.0e-6);
				} else if (real_time_mean < 1)
				{
					printf("\t: %6.2f [ms]",real_time_mean/1.0e-3);
				} else {
					printf("\t: %6.2f [s]",real_time_mean);
				}

				double tmp = sqrt(real_time_M2/(counter-1.0));
								
				if (tmp < 1e-6)
				{
					printf(" \u00B1 %6.2f [ns]\n",tmp/1.0e-9);
				} else if (tmp < 1e-3)
				{
					printf(" \u00B1 %6.2f [\xC2\xB5s]\n",tmp/1.0e-6);
				} else if (tmp < 1)
				{
					printf(" \u00B1 %6.2f [ms]\n",tmp/1.0e-3);
				} else {
					printf(" \u00B1 %6.2f [s]\n",tmp);
				}
			} else {
				cout << endl;
			}
			

		}

		
	private:
		//! numerical ID computed from 'my_name'
		int id;
		//! Given name of timer
		std::string my_name;
		//! Counter that keeps track of number of uses
		double counter;
		//! Time while timer has been in use
		double real_time;
		//! Mean time based of off 'total_time' and 'counter'
		double real_time_mean;
		//! Variance of timer computed using counter, and mean time
		double real_time_M2;
		//! Flag that keeps track of when the timer is running
		bool start_timer;
		//! Internal clocks that keep track of REAL TIME when timer was used
		struct timespec t1, t2;


		//! Create djb2 hash integer from given string
		int hash(const char *str)
		{
			unsigned long hash = 5381;
			int c;

			while (c = *str++)
			{
				hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
			}

			return hash;
		}

		//! Returns the number of digits in a number used for formating in 'print()'
		int GetNumberOfDigits (double x)
		{
			x = abs(x);  
			return (x < 10 ? 1 :   
				(x < 100 ? 2 :   
				(x < 1000 ? 3 :   
				(x < 10000 ? 4 :   
				(x < 100000 ? 5 :   
				(x < 1000000 ? 6 :   
				(x < 10000000 ? 7 :  
				(x < 100000000 ? 8 :  
				(x < 1000000000 ? 9 :  
				10)))))))));
		}

};

//! A class for timing the exec. speed of any code segment.
/*! This class creates an object that manages various timers and basic display/analysis of data.\n
    an object of this class can manage multiple object from the 'myTimer' class.\n
    The timers are divided into levels, where any timer can have a number of sub-timers included.\n
    Standalone usage:\n
    // setup \n
    myTimerCentral *a = new myTimerCentral();\n
    a->newTimer("1st main timer"); \\ create first main timer \n
    a->newSubTimer("1st main timer", "1st sub timer"); \\ create subtimer of first main timer \n
    a->newSubTimer("1st main timer", "2nd sub timer"); \\ create subtimer of first main timer \n
    a->newSubTimer("2nd sub timer", "3rd sub timer"); \\ create sub-sub-timer of first main timer i.e. a sub-timer of the 2nd sub timer \n
    a->newTimer("2nd main timer");\n
    a->newSubTimer("2nd main timer","another sub timer"); \n
    // runtime \n
    a->start("1st main timer");\n
    a->start("2nd main timer");\n
    for some number of timesteps... \n
     a->start("1st sub timer");\n
     \\do work.. \n
     a->stop("1st sub timer"); \n\n
     a->start("2nd sub timer");\n
     a->start("3rd sub timer");\n
     \\do more work... \n
     a->stop("2nd sub timer");\n
     a->stop("3rd sub timer");\n
    a->stop("1st main timer"); \n
    \\do even more work... \n
    a->stop("2nd main timer"); \n
    // analysis\n
    a->Print();
    a->Print_short();

    \sa myTimer
*/
class myTimerCentral{

	public:
		//! A constructor
		myTimerCentral()
		{
			timer_tree_root = NULL;
			MAX_SUB_TIMERS = 100; // A single node can have 10 sub timers
		}

		//! Create new top level timer with given name
		void newTimer(const std::string & newName)
		{
			if (timer_tree_root == NULL)
			{
				timer_tree_root = new timer_tree();

				timer_tree_root->next_same_level = NULL;
				timer_tree_root->next_sub_level = new timer_tree*[MAX_SUB_TIMERS];
				for(int i = 0; i < MAX_SUB_TIMERS; i++)
				{
					timer_tree_root->next_sub_level[i] = NULL;
				}

				timer_tree_root->p_timer = new myTimer(newName);

			} else {

				// Find next level
				struct timer_tree *timer_tree_tmp = timer_tree_root;
				while (timer_tree_tmp->next_same_level != NULL)
				{
					timer_tree_tmp = timer_tree_tmp->next_same_level;
				}
				timer_tree_tmp->next_same_level = new timer_tree();
				timer_tree_tmp = timer_tree_tmp->next_same_level;

				timer_tree_tmp->next_same_level = NULL;
				timer_tree_tmp->next_sub_level = new timer_tree*[MAX_SUB_TIMERS];
				for(int i = 0; i < MAX_SUB_TIMERS; i++)
				{
					timer_tree_tmp->next_sub_level[i] = NULL;
				}
			
				timer_tree_tmp->p_timer = new myTimer(newName);
			}

		}

		//! Create new sub timer under 'main_timer' with given name
		void newSubTimer(const std::string & main_timer, const std::string & newName)
		{
			int id0 = hash(main_timer.c_str());
			struct timer_tree *tmp = find_node(id0);
			
			if (tmp == NULL)
			{
				cout << "ERROR: myTimerCentral::newSubTimer() Could not locate correct timer.. spelling?" << endl;
				cout << "Tried to find timer = " << main_timer << endl;
				cout << " " << endl;
				cout << "also had input new_timer = " << newName << endl;
				exit(-1);
			}

			// Check if space for new timer
			int j = -1;
			for(int i = 0; i < MAX_SUB_TIMERS; i++)
			{
				if (tmp->next_sub_level[i] == NULL)
				{
					j = i;
					break;
				}
			}

			if (j == -1)
			{
				cout << "ERROR: myTimerCentral::newSubtimer() Don't have enough space for this sub-tumer. Increase cap.." << endl;
				cout << "Main timer name = " << main_timer << endl;
				cout << "Sub timer name  = " << newName << endl;
				exit(-1);
			}

			tmp->next_sub_level[j] = new timer_tree();
			tmp = tmp->next_sub_level[j];

			tmp->next_same_level = NULL;
			tmp->next_sub_level = new timer_tree*[MAX_SUB_TIMERS];
			for(int i = 0; i < MAX_SUB_TIMERS; i++)
			{
				tmp->next_sub_level[i] = NULL;
			}
			
			tmp->p_timer = new myTimer(newName);
		}


		//! Print tree with full information
		void Print()
		{
			cout << "\n\n" << endl;
			cout << "Print all timer information:" << endl;
			
			arg_Print(timer_tree_root,0);
		}

		//! Print only times in tree
		void Print_short()
		{
			cout << "\n\n" << endl;
			cout << "Print all timer information:" << endl;
			
			arg_Print_minimal(timer_tree_root,0, timer_tree_root->p_timer->get_real_time(), timer_tree_root->p_timer->get_counter());

		}

		//! Start counter in given timer
		void start(const std::string & main_timer)
		{
			int id0 = hash(main_timer.c_str());
			struct timer_tree *tmp = find_node(id0);
			if (tmp == NULL)
			{
				cout << "ERROR: myTimerCentral::start() Could not locate correct timer.. spelling?" << endl;
				cout << "Tried to find timer = " << main_timer << endl;
				exit(-1);
			}
			
			tmp->p_timer->start();
		}
		//! Start counter in given sub-timer
		void start(const std::string & main_timer, const std::string & sub_timer)
		{
			int id0 = hash(main_timer.c_str());
			struct timer_tree *tmp = find_node(id0);
			if (tmp == NULL)
			{
				cout << "ERROR: myTimerCentral::start(main,sub) Could not locate correct MAIN-timer.. spelling?" << endl;
				cout << "Tried to find MAIN timer  = " << main_timer << endl;
				cout << "with additional sub timer = " << sub_timer << endl;
				exit(-1);
			}

			int id1 = hash(sub_timer.c_str());
			struct timer_tree *tmp2 = find_node(id1,tmp);
			if (tmp2 == NULL)
			{
				cout << "ERROR: myTimerCentral::start(main,sub) Could not locate correct SUB-timer.. spelling?" << endl;
				cout << "Tried to find MAIN timer  = " << main_timer << endl;
				cout << "with additional sub timer = " << sub_timer << endl;
				exit(-1);
			}
			tmp2->p_timer->start();
		}
		//! Stop counter in given timre
		void stop(const std::string & main_timer)
		{
			int id0 = hash(main_timer.c_str());
			struct timer_tree *tmp = find_node(id0);
			if (tmp == NULL)
			{
				cout << "ERROR: myTimerCentral::stop() Could not locate correct timer.. spelling?" << endl;
				cout << "Tried to find timer = " << main_timer << endl;
				exit(-1);
			}
			
			tmp->p_timer->stop();
		}
		//! Stop counter in given sub-timer
		void stop(const std::string & main_timer, const std::string & sub_timer)
		{
			int id0 = hash(main_timer.c_str());
			struct timer_tree *tmp = find_node(id0);
			if (tmp == NULL)
			{
				cout << "ERROR: myTimerCentral::stop(main,sub) Could not locate correct MAIN-timer.. spelling?" << endl;
				cout << "Tried to find MAIN timer  = " << main_timer << endl;
				cout << "with additional sub timer = " << sub_timer << endl;
				exit(-1);
			}

			int id1 = hash(sub_timer.c_str());
			struct timer_tree *tmp2 = find_node(id1,tmp);
			if (tmp2 == NULL)
			{
				cout << "ERROR: myTimerCentral::stop(main,sub) Could not locate correct SUB-timer.. spelling?" << endl;
				cout << "Tried to find MAIN timer  = " << main_timer << endl;
				cout << "with additional sub timer = " << sub_timer << endl;
				exit(-1);
			}
			tmp2->p_timer->stop();
		}

	private:
		//! Timer object organized with pointer to other levels of timers
		struct timer_tree {

			class myTimer *p_timer;
			
			struct timer_tree *next_same_level;
			struct timer_tree **next_sub_level;
		};

		//! Main root of timer tree
		struct timer_tree *timer_tree_root;

		//! Maximal number of sub_timers (default = 55)
		int MAX_SUB_TIMERS;
		
		//! Print information helper function
		void arg_Print(struct timer_tree *tmp, int step)
		{
			if (tmp!= NULL)
			{
				tmp->p_timer->Print(step);

				for(int i = 0; i < MAX_SUB_TIMERS; i++)
				{
					arg_Print(tmp->next_sub_level[i],step+1);
				}
				arg_Print(tmp->next_same_level, step);
			}
			
		}
		//! Print short information helper function
		void arg_Print_minimal(struct timer_tree *tmp, int step, double parent_time, double parent_count)
		{
			if (tmp!= NULL)
			{
				tmp->p_timer->Print_minimal(step, parent_time, parent_count);

				for(int i = 0; i < MAX_SUB_TIMERS; i++)
				{
					arg_Print_minimal(tmp->next_sub_level[i],step+1, tmp->p_timer->get_real_time(), tmp->p_timer->get_counter());
				}
				arg_Print_minimal(tmp->next_same_level, step, tmp->p_timer->get_real_time(), tmp->p_timer->get_counter());
			}
			
		}

		//! Find and return a pointer to a sub-timer from the given 'target_id'
		/*! This is a helper function for 'find_node()' 
		    \param tmp structure to search through
		    \param target_id numerical id of search target
		    \return NULL if nothing is found matching target_id
		*/
		struct timer_tree * arg_find_node(struct timer_tree *tmp, int target_id)
		{
			if (tmp != NULL)
			{
				if (tmp->p_timer->getId()==target_id)
				{
					return tmp;
				}
				
				for(int i = 0; i < MAX_SUB_TIMERS; i++)
				{
					struct timer_tree *tmp2 =  arg_find_node(tmp->next_sub_level[i], target_id);
					if (tmp2 != NULL)
					{
						return tmp2;
					}
				}
				struct timer_tree *tmp3 =  arg_find_node(tmp->next_same_level, target_id);
				if (tmp3 != NULL)
				{
					return tmp3;
				}

			}
			return NULL;
		}
		//! Find and return a pointer to a timer for the given 'target_id'
		struct timer_tree * find_node(int target_id)
		{
			struct timer_tree *tmp = arg_find_node(timer_tree_root, target_id);

			if (tmp != NULL)
			{
				return tmp;
			}
			return NULL;
		}
		//! Find and return a pointer to a timer for the given 'target_id' starting from 'start_node'
		struct timer_tree * find_node(int target_id, struct timer_tree *start_node)
		{
			struct timer_tree *tmp = arg_find_node(start_node, target_id);

			if (tmp != NULL)
			{
				return tmp;
			}
			return NULL;
		}


		//! Create djb2 hash integer from given string
		int hash(const char *str)
		{
			unsigned long hash = 5381;
			int c;

			while (c = *str++)
			{
				hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
			}

			return hash;
		}
};



#endif
