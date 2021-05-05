
/*
	The CYCLIC class is a storage system used for fast easy storage of
	arbitrary elements.
	Isak Kilen @ 2016
*/

#ifndef __CYCLIC_H_INCLUDED__
#define __CYCLIC_H_INCLUDED__

#include <string>
#include <cstring>

#include <stdlib.h>
#include <iostream>

#include <fstream>
#include <sstream>

//#include <sys/stat.h>
#include "fileIO.cpp"

using namespace std;


//! A container to store propagating electromagnetic fields
/*! The storage is optimized for the transfer matrix method where one want to cycle the stored components forward in time.\n
    After initzialization the 'updateStorage()' is called to cycle the first element to the second storage, 2nd to 3rd, etc..\n
    Then the first element is free to be overwritten. \n

    The class creates containers of a given type that can store store a single number or an array of the given type.

    2018: Added more comments\n
    Isak Kilen

    2020: Added some additional features for looking at a prescribed additional container
    S. A. McLaren
*/
template <class A>
class Cyclic
{
	public:
		//! A constructor
		Cyclic()
		{
			head = NULL;
			tail = NULL;
			exter_iter = NULL;
			number_containers = 0;
			number_elements = 0;
		}
		
		//! A destructor
		~Cyclic()
		{
			deleteAllContainers();
		}
		
		//! A constructor
		Cyclic(int numContainers)
		{
			head = NULL;
			tail = NULL;
			exter_iter = NULL;
			number_containers = numContainers;
			number_elements = 1;
			
			if (numContainers <= 0)
			{
				cout << "Cyclic(a): Need number of containers >= 0" << endl;
				cout << "Cyclic(a): got numContainers = " << numContainers << endl;
				exit(-1);
			}
	
			for(unsigned i = 0; i < numContainers; i++)
			{
				addContainer();
			}
		}
		
		//! A constructor
		Cyclic(int numContainers,int numElements)
		{
			head = NULL;
			tail = NULL;
			exter_iter = NULL;
			number_containers = numContainers;
			number_elements = numElements;
			
			if (numElements <= 0)
			{
				cout << "Cyclic(): Need number of elements > 0" << endl;
				exit(-1);
			}
			
			if (numContainers <= 0)
			{
				cout << "Cyclic(a,b): Need number of containers >= 0" << endl;
				cout << "Cyclic(a,b): got numContainers = " << numContainers << endl;
				cout << "Cyclic(a,b): got numElements   = " << numElements << endl;
				exit(-1);
			}
	
			for(unsigned i = 0; i < numContainers; i++)
			{
				addContainer();
			}
		}

		//! A copy-constructor
		Cyclic(const Cyclic<A> &obj)
		{
			head = NULL;
			tail = NULL;
			exter_iter = NULL;

			number_containers = obj.number_containers;
			number_elements = obj.number_elements;
			for(unsigned i = 0; i < number_containers; i++)
			{
				addContainer();
			}

			// Copy content
			for(unsigned i = 0; i < number_containers; i++)
			{
				setContainerNr(i, obj.getContainerNr(i));
			}
		}

		//! Print container information to file
		void Print() const
		{
			cout << "Cyclic:" << endl;
			cout << "Size of container  = " << sizeof(struct container) << endl;
			cout << "Number of containers = " << number_containers << endl;
			cout << "Number of elements in each container = " << number_elements << endl;
			struct container *iter = head;
			while(iter != NULL)
			{
				cout << "id = " << iter->id << " | ";
				cout << endl;
				iter = iter->next;

			}
			
			// Check consistency of containers
			bool check = false;
			int code = -1;
			if ((head != NULL)&&(head == tail))
			{
				if (head->prev != NULL) {code = 1;}
				if (head->next != NULL) {code = 2;}
			}
			if ((head != NULL)&&(head != tail))
			{
				if (head->prev == NULL) {code = 3;}
				if (head->next != NULL) {code = 4;}
				if (tail->prev != NULL) {code = 5;}
				if (tail->next == NULL) {code = 6;}
				
				iter = head;
				while (iter->next != NULL)
				{
					if (iter->prev == NULL) {code = 7;}
					if (iter->next == NULL) {code = 8;}
					iter = iter->next;
				}
				
			}
			
			cout << "Container consistency: ";
			if (check)
			{
				cout << "Failed! Code = " << code << endl;
				exit(-1);
			} else {
				cout << "Ok!" << endl;
			}
		}
	
		//! Fill an empty cyclic object
		/*! Takes container and element info and fills a cyclic object.*/
		void initializeCyclic(int numContainers, int numElements)
		{
			head = NULL;
			tail = NULL;
			exter_iter = NULL;
			number_containers = numContainers;
			number_elements = numElements;
			
			if (numElements <= 0)
			{
				cout << "Cyclic(): Need number of elements > 0" << endl;
				exit(-1);
			}
			
			if (numContainers <= 0)
			{
				cout << "Cyclic(a,b): Need number of containers >= 0" << endl;
				cout << "Cyclic(a,b): got numContainers = " << numContainers << endl;
				cout << "Cyclic(a,b): got numElements   = " << numElements << endl;
				exit(-1);
			}
	
			for(unsigned i = 0; i < numContainers; i++)
			{
				addContainer();
			}

		}
	
		//! Cylcle storage one step forward
		/*! Moves the 1st container to the 2nd location, the 2nd container to the 3rd location, etc. The final container becomes available as the 'new first' container.*/
		void updateStorage(void)
		{
			struct container *tmp = tail;
			// Fix tail
			tail = tail->prev;
			tail->next = NULL;
	
			// Fix last node
			tmp->prev = NULL;
			tmp->next = head;
	
			// Fix head
			head->prev = tmp;
			head = tmp;	
		}

		//! Set all container elements to the value '0.0'
		void zeroAllFields(void)
		{
			struct container *iter = head;
			while (iter != NULL)
			{
				for (unsigned i = 0; i < number_elements; i++)
				{
					iter->field[i] = 0;
				}
				iter = iter->next;
			}
		}

		//! Remove all containers from this structure
		void deleteAllContainers()
		{
			while (head != NULL)
			{
				deleteContainer(0);
			}
	
			head = NULL;
			tail = NULL;
		}
		
		//! Set the contents of the first container to a new value
		inline void setFirstContainer(A *new_field)
		{
			memcpy(head->field,new_field,number_elements*sizeof(A));
		}

		//! Return a pointer to the first container value
		inline A *getFirstContainer(void) const
		{
			return head->field;
		}
		//! Return a pointer to the second container value
		inline A *get2ndContainer(void) const
		{
			return head->next->field;
		}

		//! Set the contents of the last container to a new value
		inline void setLastContainer(A *new_field)
		{
			memcpy(tail->field,new_field,number_elements*sizeof(A));
		}

		//! Return a pointer to the last container value
		inline A *getLastContainer(void) const
		{
			return tail->field;
		}
		
		//! Return a pointer to the value in a given container
		inline A * getContainerNr(int nr) const
		{
			int i = 0;
			struct container *iter = head;
			if (nr >= number_containers)
			{
				cout << "getContainerNr(): Requesting container # not in collection" << endl;
				exit(-1);
			}
			while (i != nr)
			{
				iter = iter->next;
				i++;
			}
			return iter->field;
		}
		
		//! Return a pointer to the value in a given container. Cycle direction is based on which end is closer
		inline A getContainerNr_point(int nr, int ncav) const
		{
			if (nr >= number_containers)
			{
				cout << "getContainerNr(): Requesting container # not in collection" << endl;
				exit(-1);
			}
			
			int i = 0;
			struct container *iter;
			if (nr < number_containers/2)
			{
				iter=head;
				while (i != nr)
				{
					iter = iter->next;
					i++;
				}
			} else
			{
				iter=tail;
				nr=number_containers-nr;
				while (i != nr)
				{
					iter = iter->prev;
					i++;
				}
			} 
			return iter->field[ncav];
		}
			
		//! Set the value in a given container
		inline void setContainerNr(int nr, A *new_field) const
		{
			A *field = getContainerNr(nr);
			memcpy(field,new_field,number_elements*sizeof(A));
		}
		
		//! Return a pointer to the value inside the 2nd last container	
		inline A * get2ndLastContainer()
		{
			return tail->prev->field;
		}

		//! Return a pointer to the value inside the 3rd last container	
		inline A * get3rdLastContainer()
		{
			return tail->prev->prev->field;
		}
		
		//! Save all container data to file with the given filename
		/*! Only the container data is saved, i.e. not any other configuration */
		void saveAllContainers(const std::string &fileName)
		{
			// Create block of memory for storage
			A tmp[number_elements];
			struct container *iter = head;
			while (iter != NULL)
			{
				memcpy(tmp,iter->field,number_elements*sizeof(A));

				iter = iter->next;
				saveAppendBinary(fileName, tmp, number_elements);
			}
		}

		//! Load all container data form file with given filename
		/*! Assumes that the current datastructure is the same that was saved. Only data is loaded*/
		void loadAllContainers(const std::string &fileName)
		{
			A tmp[number_elements];
			struct container *iter = head;
			unsigned i = 1;
			while (iter != NULL)
			{
				loadBinaryLine(fileName, tmp, number_elements, i);
				memcpy(iter->field,tmp,number_elements*sizeof(A));

				i++;
				iter = iter->next;
			}
		}

		//! Diagnostics/Interface:  For iterating through the nodes from the outside
		inline void externalIter_reset(void)
		{
			exter_iter = head;
		}

		//! Diagnostics/Interface:  For iterating through the nodes from the outside
		inline void externalIter_step(void)
		{
			exter_iter = exter_iter->next;
		}

		//! Diagnostics/Interface:  For iterating through the nodes from the outside
		inline A * externalIter_get(void)
		{
			return exter_iter->field;
		}
		
		//! Diagnostics/Interface:  For iterating through the nodes from the outside
		inline A externalIter_getField( int i )
		{
			return exter_iter->field[i];
		}

		//! Diagnostics/Interface: For setting exter_iter.
		inline void externalIter_setContainerNr( int i )
		{
			exter_iter->getContainerNr(i);
		} 
		
	private:
		//! A container for the stored data and pointers to the next elements in the list
		struct container
		{
			A *field;
			//int number_of_elements;
			int id;
			
			struct container *next;
			struct container *prev;
		};
		
		//! Number of containers in the class
		int number_containers;

		//! Number of elements in each container
		int number_elements;
		
		//! Pointer to the current head of the list
		struct container *head;

		//! Pointer to the current tail of the list
		struct container *tail;

		//! Iterator for running through the data
		struct container *exter_iter;
		
		//! Private function to add containers to the list
		/*! Containers are allways added at the end of the list */
		void addContainer()
		{
			if (head == NULL)
			{
				// No containers in list
				head = new container;
				head->id = number_containers;
				head->next = NULL;
				head->prev = NULL;
			
	
				head->field = new A[number_elements];
				for (unsigned i = 0; i < number_elements; i++)
				{
					head->field[i] = 0;
				}
		
				tail = head;
		
			} else {
				// Add at the end
				struct container *tmp = new container;
				tmp->id = number_containers;
				tmp->prev = tail;
				tmp->next = NULL;
				
				tmp->field = new A[number_elements];
				for (unsigned i = 0; i < number_elements; i++)
				{
					tmp->field[i] = 0;
				}
				
				tail->next = tmp;
				tail = tmp;
			}
		}
		
		/*
			Delete a container with number 'num'
			'num' is between [0,inf)
		*/

		//! Private function to delete a container from the list
		/*! The number is relative to the head node (num==0)*/
		void deleteContainer(int num)
		{
			if (head != NULL)
			{
				struct container *tmp = head;
				if (num == 0)
				{
					// Delete first node
					if (head->next != NULL)
					{
						head = head->next;
						head->prev = NULL;
					} else {
						head = NULL;
					}
				} else if (num > 0)
				{
			
					//	Find container number 'num'
			
					int counter = 0;
					while (tmp != NULL)
					{
						if (counter == num)
						{
							break;
						}
			
						counter++;	
						tmp = tmp->next;
					}
			
					if (tmp == NULL)
					{
						cout << "Cyclic::deleteContainer(): Could not find the correct container" << endl;
						cout << "Looking for number " << num << endl;
						exit(-1);
					}
			
					// Previous node points to next
					tmp->prev->next = tmp->next;
			
					// Next node is
					if (tmp != tail)
					{		
						// A middle node
						tmp->next->prev = tmp->prev;
					} else {
			
						// Move tail pointer
						tail = tmp->prev;
					}
				} else if (num < 0)
				{
					cout << "Cyclic::deleteContainer(): Cannot delete negative index" << endl;
					cout << "Looking for number " << num << endl;
					exit(-1);
				}
		
				if (tmp != NULL)
				{
					tmp->prev = NULL;
					tmp->next = NULL;
					delete [] tmp->field;
					delete  tmp;
				}
		
			} else {
				cout << "Cyclic::deleteContainer(): Nothing to delete" << endl;
				cout << "Looking for number " << num << endl;
				exit(-1);
			}
	
			number_containers--;
		}
		
};



#endif
