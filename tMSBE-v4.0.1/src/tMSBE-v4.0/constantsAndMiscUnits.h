
#ifndef __CONSTANTSANDMISCUNITS_H_INCLUDED__
#define __CONSTANTSANDMISCUNITS_H_INCLUDED__

#include "complex"

// Permanent constants
const double Pi   = acos(-1.0);                           //  Pi
const double c0   = 299792458;                           //  light velocity in m/s
const double e    = 1.6021773349e-19;                           //  electron charge in As 
const double m0   = 9.109389754e-31;                           //  free electron mass in kg
const double mu0  = 4.0e-7*Pi;                           //  permeability in Vs/Am
const double eps0 = 1.0/(mu0*c0*c0);                           //  vac.diel.const. in As/Vm
const double h    = 6.626075540e-34;                           //  Planck constant in Js
const double hbar = h/(2.0*Pi);                           //  h*2*Pi
const double kB   = 1.38065812e-23;                           //  Boltzmann constant in J/K
const std::complex<double> I(0,1);                          //  imaginary unit
//const double eps  = 10.970339;                           //  rel.diel.constant
//const double a0   = 1.062146e-08;                           //  Bohr radius in m


// Jorg New units April 13 2018
const double eps  = 11.0028;                           //  rel.diel.constant
const double a0   = 1.0650426e-08;                           //  Bohr radius in m


// Time scale
const double fs   	= 1.0e-15;
const double ps 	= 1.0e-12;
const double ns 	= 1.0e-9;

// Length scales
const double nm 	= 1.0e-9;
const double um 	= 1.0e-6;
const double cm 	= 1.0e-2;


#endif
