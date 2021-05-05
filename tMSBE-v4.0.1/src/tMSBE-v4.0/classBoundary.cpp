
#include <iostream>
#include "classBoundary.h"

using namespace std;

/*
	Construction function
*/

Boundary::Boundary()
{
	setRefCoeff(0);
	boundary_next_cavity_index = 1.0;
	boundary_Epluss = NULL;
	boundary_Eminus = NULL;
}

Boundary::Boundary(double new_ref_coeff)
{
	setRefCoeff(new_ref_coeff);
	boundary_Epluss = NULL;
	boundary_Eminus = NULL;
	boundary_next_cavity_index = 1.0;
}

void Boundary::Print() const
{
	cout << "Print boundary:" << endl;
	cout << " -> Reflection coeff  = " << getRefCoeff() << endl;
}

void Boundary::initializeZero(int numT)
{
	number_transverse_points = numT;
	boundary_Epluss = new std::complex<double>[number_transverse_points];
	boundary_Eminus = new std::complex<double>[number_transverse_points];
	for(int i = 0; i < number_transverse_points; i++)
	{
		boundary_Epluss[i] = 0.0;
		boundary_Eminus[i] = 0.0;
	}
}
