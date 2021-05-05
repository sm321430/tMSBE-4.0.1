

/*!
	Here are a few functions for input and output of the simualtino data
	The data can be compressed into binary files\n
	There are functions for checking if a directory/file exists\n
	The default behavior is to overwrite files with the same name.\n
	2018: Additional comments \n
	2016: Build and design \n
	Isak Kilen
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

//#include <errno.h>
#include <sys/stat.h>
#include <string>

#ifndef __FILEIO_CPP_INCLUDED__
#define __FILEIO_CPP_INCLUDED__

using namespace std;


static bool fileExists(const std::string& filename);
static bool  dirExists(const std::string& direname);
static bool  dirMake(const std::string& dirname);

template <typename A>
static int detectNumberOfLines(const std::string&, A*, int);

template <typename A>
static void saveBinary(const std::string&, A *, int);

template <typename A>
static void loadBinary(const std::string&, A *, int);

template <typename A>
static void saveAppendBinary(const std::string&, A *, int);

template <typename A>
static void loadBinaryLine(const std::string&, A *, int, int);

template <typename A>
static void loadBinaryIndex(const std::string&, A *, int, int);

static void openAppendBinary(std::ofstream *,const std::string&);




//! Check if a file exists
/*! \param[in] filename - the name of the file to check
    \return    true if the file exists, else false
*/
static bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}

//! Check if a directory exists -> First check if it exists, then if it is a directory
/*! \param[in] direname - the name of the directory to check
   \return    true if the dir exists, else false
*/
static bool dirExists(const std::string& direname)
{
    struct stat buf;
    if ((stat(direname.c_str(), &buf) != -1)&&(S_ISDIR(buf.st_mode)))
    {
        return true;
    }
    return false;
}

//! Create a directory in currect directory
/*! \param[in] direname - the name of the directory to create
    \return    true if sucsessfull, false if failed
*/

static bool dirMake(const std::string& dirname)
{
	if (mkdir(dirname.c_str(), 0777)!=-1) 
	{
		return true;
	}
    return false;
}

//! Load from binary file
/*! If file does not exist: Will return 0 else return 1 */

template <typename A>
static void loadBinary(const std::string& filename, A *store, int num)
{
	if (fileExists(filename))
	{
		std::ifstream loadOut;
		loadOut.open(filename.c_str(), std::ifstream::in|std::ifstream::binary);
		if (num > 0)
		{
			loadOut.read((char *)(store),num*sizeof(A));
		} else {
			cout << "LOAD BINARY: Cannot load <=0 number of numbers" << endl;
			exit(-1);
		}
		loadOut.close();
	} else {
		// Failed to find file to load
		cout << "LOAD BINARY: failed to find file: " << filename << endl;
		exit(-1);
	}
}

//! Load a single line 'num' from a file 'filename' and store it in 'store'
/*!
	Load 'num' into A from line in file\n
	\param num is length of line\n
	\param line is what line number [1,2,3,...]
*/

template <typename A>
static void loadBinaryLine(const std::string& filename, A *store, int num, int line)
{
	if (fileExists(filename))
	{
		std::ifstream loadOut;
		loadOut.open(filename.c_str(), std::ifstream::in|std::ifstream::binary);
		if (num > 0)
		{
			if (line > 0)
			{
				// Set to beginning
				for(unsigned i = 0; i < line; i++)
				{
					if (!loadOut.read((char *)(store),num*sizeof(A)))
					{
						cout << "LOAD BINARY LINE: line > # of lines in file" << endl;
						exit(-1);
					}
				}
				
			} else {
				cout << "LOAD BINARY LINE: Cannot load line <= 0" << endl;
				exit(-1);
			}
			
		} else {
			cout << "LOAD BINARY LINE: Cannot load <= 0 number of numbers" << endl;
			exit(-1);
		}
		loadOut.close();
	} else {
		// Failed to find file to load
		cout << "LOAD BINARY LINE: failed to find file: " << filename << endl;
		exit(-1);
	}
}

//! Load a single an array of length 'num' from a binary file 'filename' and store it in 'store'
/*!
	Load 'num' into A from line in file\n
	\param num is length of data to load\n
	\param index0 is the starting index of data to load\n
*/

template <typename A>
static void loadBinaryIndex(const std::string& filename, A *store, int num, int index0)
{
	if (fileExists(filename))
	{
		std::ifstream loadOut;
		loadOut.open(filename.c_str(), std::ifstream::in|std::ifstream::binary);
		if (num > 0)
		{
			if (index0 > 0)
			{
				A dummy;
				for(unsigned i = 0; i < index0; i++)
				{
					if (!loadOut.read((char *)(&dummy),sizeof(A)))
					{
						cout << "LOAD BINARY INDEX: index0 > # of lines in file" << endl;
						exit(-1);
					}
				}

				if (!loadOut.read((char *)(store),num*sizeof(A)))
				{
					cout << "LOAD BINARY INDEX: index0 > # of lines in file" << endl;
					exit(-1);
				}

				
			} else {
				cout << "LOAD BINARY INDEX: Cannot load index0 <= 0" << endl;
				exit(-1);
			}
			
		} else {
			cout << "LOAD BINARY INDEX: Cannot load <= 0 number of numbers" << endl;
			exit(-1);
		}
		loadOut.close();
	} else {
		// Failed to find file to load
		cout << "LOAD BINARY INDEX: failed to find file: " << filename << endl;
		exit(-1);
	}
}

//! Return number of lines in binary file with column length 'num'
/*! \param num is number of elements in each column of the file
    \param filename is the name of the given file
*/

template <typename A>
static int detectNumberOfLines(const std::string& filename, A*dummy, int num)
{
	if (fileExists(filename))
	{
		std::ifstream loadOut;
		loadOut.open(filename.c_str(), std::ifstream::in|std::ifstream::binary);
		if (num > 0)
		{
			A *tmpStore = new A[num];
			// Set to beginning
			int i = 0;
			while (loadOut.read((char *)(tmpStore),num*sizeof(A)))
			{
				i++;
			}
			
			return i;
			
		} else {
			cout << "detectNumberOfLines: Cannot load <= 0 number of numbers" << endl;
			exit(-1);
		}
		loadOut.close();
	} else {
		// Failed to find file to load
		cout << "detectNumberOfLines: failed to find file: " << filename << endl;
		exit(-1);
	}
}

//! Save information in a binary file with truncation
/*! If file exists: Will truncate file.\n
\param filename name of target file
\param store data to save to file
\param num number of elements in data to store
*/

template <typename A>
static void saveBinary(const std::string& filename, A *store, int num)
{	
	if (num >0)
	{
		std::ofstream saveOut;
		saveOut.open(filename.c_str(), std::ofstream::out|std::ofstream::binary|std::ofstream::trunc);
		saveOut.write(reinterpret_cast<const char*>(store),num*sizeof(A));
		saveOut.close();
	} else {
		cout << "SAVE BINARY: Cannot save <=0 number of numbers" << endl;
		cout << "filename = " << filename << endl;
		cout << "num      = " << num << endl;
		exit(-1);
	}
}

//! Save information in a binary file
/*! If file exists: Will append to the end.\n
\param filename name of target file
\param store data to save to file
\param num number of elements in data to store
*/

template <typename A>
static void saveAppendBinary(const std::string& filename, A *store, int num)
{	
	if (num >0)
	{
		std::ofstream saveOut;
		saveOut.open(filename.c_str(), std::ofstream::out|std::ofstream::binary|std::ofstream::app);
		saveOut.write(reinterpret_cast<const char*>(store),num*sizeof(A));
		saveOut.close();
	} else {
		cout << "SAVE APPEND BINARY: Cannot save <=0 number of numbers" << endl;
		cout << "filename = " << filename << endl;
		cout << "num      = " << num << endl;
		exit(-1);
	}
}

//! Open binary file and append at the end
/*! Open stream for appending
	If file exists: Will delete so as to truncate the file
	If not: Just open.
*/

static void openAppendBinary(std::ofstream *saveOut,const std::string& filename)
{
	if (fileExists(filename))
	{
		// File exists, delete and open
		cout << "OPEN APPEND BINARY: FOUND AND DELETING FILE: " << filename << endl;
		std::remove(filename.c_str());
		//exit(-1);
	}
	
	(*saveOut).open(filename.c_str(), std::ofstream::out|std::ofstream::binary|std::ofstream::app);
}



#endif
