#include <iostream>
#include "classLossElement.h"

using namespace std;

/*
	Construction function
*/

LossElement::LossElement()
{
	loss_plus = 0;
	loss_minus = 0;
}

LossElement::LossElement(double lm, double lp)
{
	loss_minus = lm;
	loss_plus = lp;
}

void LossElement::Print() const
{
	cout << "Print loss element:" << endl;
	cout << " -> Loss -  = " << loss_minus << endl;
	cout << " -> Loss +  = " << loss_plus << endl;
}
