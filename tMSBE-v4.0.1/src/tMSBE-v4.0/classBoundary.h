
#ifndef __BOUNDARY_H_INCLUDED__
#define __BOUNDARY_H_INCLUDED__

#include <string>
#include <complex>
#include <cstring>


//! The BOUNDARY class is placed inside a 'Module()' on the edges of the domain.
/*! 
	In order to serve as boundaries for the computational domain and can be used as frequency independent reflective mirrors with a given amount of loss.\n	

	The boundary is also used to store the incoming E+/E- field from outside of the domain.
	\sa Module

	2018: Added more comments \n
	Isak Kilen
*/
class Boundary
{
	public:
		//! A constructor
		Boundary();

		//! A constructor
		Boundary(double);

		//! A destructor
		~Boundary() 
		{
			if (boundary_Epluss!=NULL)
			{
				delete [] boundary_Epluss;
				delete [] boundary_Eminus;
			}
		};

		//! A copy-constructor
		Boundary(const Boundary &obj)
		{
			boundary_reflection_coefficient = obj.boundary_reflection_coefficient;
			boundary_Epluss = obj.boundary_Epluss;
			boundary_Eminus = obj.boundary_Eminus;
			boundary_next_cavity_index = obj.boundary_next_cavity_index;

			if (obj.boundary_Epluss != NULL)
			{
				number_transverse_points = obj.number_transverse_points;
				boundary_Epluss = new std::complex<double>[number_transverse_points];
				boundary_Eminus = new std::complex<double>[number_transverse_points];
			}
		}

		//! Print all information to screen
		void Print( void ) const;

		//! Return the reflection coefficient
		double getRefCoeff(void) const;

		//! Set the reflection coefficient
		void   setRefCoeff(double);

		//! Initialize storage for boundary
		void initializeZero(int);

		//! Return E+ stored at the boundary.
		void getEpluss(std::complex<double> *);

		//! Set E+ stored in this boundary
		void setEpluss(std::complex<double>*);
		std::complex<double> * setEpluss(void);

		//! Return E- stored at the boundary
		void getEminus(std::complex<double> *);

		//! Set E- stored at the boundary
		void setEminus(std::complex<double>*);
		std::complex<double> * setEminus(void);
		
		//! Set the external index of refraction. (default=1)
		void setNextCavityIndex(double);

		//! Return the external index of refraction
		double getNextCavityIndex(void);
		
	private:
		//! Reflection coefficient used in propagation code
		double boundary_reflection_coefficient;

		//! The refractive index of the material outside the domain.
		double boundary_next_cavity_index;

		//! number of transverse points
		int number_transverse_points;

		//! Input E+ stored at each transverse point of the boundary
		std::complex<double> *boundary_Epluss;

		//! Input E- stored at each transverse point of the boundary
		std::complex<double> *boundary_Eminus;
};

inline void Boundary::setRefCoeff(double newCoeff)
{
	boundary_reflection_coefficient = newCoeff;
}

inline double Boundary::getRefCoeff(void) const
{
	return boundary_reflection_coefficient;
}

inline void Boundary::setEpluss(std::complex<double> *Ep)
{
	memcpy(boundary_Epluss, Ep, number_transverse_points*sizeof(std::complex<double>));
/*
	for(int i = 0; i < number_transverse_points; i++)
	{
		boundary_Epluss[i] = Ep[i];
	}
*/
}

inline void Boundary::setEminus(std::complex<double> *Em)
{
	memcpy(boundary_Eminus, Em, number_transverse_points*sizeof(std::complex<double>));
/*
	for(int i = 0; i < number_transverse_points; i++)
	{
		boundary_Eminus[i] = Em[i];
	}
*/
}

inline std::complex<double> * Boundary::setEpluss()
{
	return boundary_Epluss;
}

inline std::complex<double> * Boundary::setEminus()
{
	return boundary_Eminus;
}

inline void Boundary::getEpluss(std::complex<double> * tmp)
{
	for(int i = 0; i < number_transverse_points; i++)
	{
		tmp[i] = boundary_Epluss[i];
	}
}

inline void Boundary::getEminus(std::complex<double> * tmp)
{
	for(int i = 0; i < number_transverse_points; i++)
	{
		tmp[i] = boundary_Eminus[i];
	}
}

inline void Boundary::setNextCavityIndex(double n)
{
	boundary_next_cavity_index = n;
}

inline double Boundary::getNextCavityIndex(void)
{
	return boundary_next_cavity_index;
}




#endif












