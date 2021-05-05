
#ifndef __LOSSELEMENT_H_INCLUDED__
#define __LOSSELEMENT_H_INCLUDED__

#include <string>
#include <complex>

//! The LOSS ELEMENT class is placed inside the 'Module()' container wherever needed.
/*! 
   The loss can be in one direction or in both. There a cavity on both sides.\n
   \sa Module \n

   2018: Additonal comments added\n
   2016: start\n
   Isak Kilen
*/

class LossElement
{
	public:
		//! A constructor
		LossElement();
		//! A constructor
		LossElement(double,double);
		//! Print information about the loss element to screen
		void Print( void ) const;
		//! Return loss coefficient in both directions
		void  getLossCoeff(double*, double*) const;
		//! Set the loss coefficients in both directions
		void  setLossCoeff(double,double);
		
	private:
		//! Loss for right moving wave
		double loss_plus;
		//! Loss for left moving wave
		double loss_minus;
};

inline void LossElement::setLossCoeff(double new_plus, double new_minus)
{
	loss_plus  = new_plus;
	loss_minus = new_minus;
}

inline void LossElement::getLossCoeff(double *data_plus, double *data_minus) const
{
	*data_plus  = loss_plus;
	*data_minus = loss_minus;
}

#endif












