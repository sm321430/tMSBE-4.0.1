
/*
	The MODULE class is a container for any of the other
	classes: Cavity, Device,...
	Isak Kilen @ 2016

	Modified to include TwoArm structures
	S. A. McLaren @ 2020	
*/

#include <iostream>
#include <stdlib.h>
#include "classModule.h"

using namespace std;

Module::Module()
{
	cavity = NULL;
	device = NULL;
	twoArmCavity = NULL;
	birefringentCrystal = NULL;
	twoArmDevice = NULL;
	twoArmInterface = NULL;
	boundary = NULL;
	lossElement = NULL;
	filter = NULL;
	
	module_position_x0 = 0.0;
	module_position_x1 = 0.0;
	module_width = 0.0;
	module_refractive_index = 0.0;
	module_refractive_index_im = 0.0;
	module_cos_th_left = 1.0;
	module_cos_th_right = 1.0;
	module_angle = 0.0;
	module_intAngle = 0.0;
	module_prev_cav = 0;
	module_post_cav = 0;	

	output_to_file = 0;
	output_to_file_level = 0;
}

void Module::setRefInd(double n1)
{
	if (n1 <= 0.0)
	{
		cout << "Module::setRefInd: Need refractive index >0 (= " << n1 << ")" << endl;
		exit(-1);
	} else if (isCavity())
	{
		module_refractive_index = n1;
		cavity->setRefInd(n1);
	} else if (isTwoArmCavity())
	{
		module_refractive_index = n1;
		twoArmCavity->setRefInd(n1);
	} else if (isBirefringentCrystal())
	{
		module_refractive_index = n1;
		birefringentCrystal->setRefInd(n1);
	} else if (isTwoArmInterface())
	{
		module_refractive_index = n1;
		twoArmInterface->setRefInd(n1);
	} else if (isBoundary())
	{
		module_refractive_index = 0.0;	
	} else if (isDevice())
	{
		module_refractive_index = 0.0;
	} else if (isTwoArmDevice())
	{
		module_refractive_index = 0.0;
	} else if (isLossElement())
	{
		module_refractive_index = 0.0;
	} else if (isFilter())
	{
		module_refractive_index = 0.0;
	} else {
		cout << "Module::setRefInd(): Unknown module" << endl;
		exit(-1);
	}
}


void Module::setRefInd_extraAxis(double n1)
{
	if (n1 <= 0.0)
	{
		module_refractive_index_extraAxis = 0.0;
	} else if (isCavity())
	{
		module_refractive_index_extraAxis = 0.0;
	} else if (isTwoArmCavity())
	{
		module_refractive_index_extraAxis = 0.0;
	} else if (isBirefringentCrystal())
	{
		module_refractive_index_extraAxis = n1;
		birefringentCrystal->setRefInd_extraAxis(n1);
	} else if (isTwoArmInterface())
	{
		module_refractive_index_extraAxis = 0.0;
	} else if (isBoundary())
	{
		module_refractive_index_extraAxis = 0.0;
	} else if (isDevice())
	{
		module_refractive_index_extraAxis = 0.0;
	} else if (isTwoArmDevice())
	{
		module_refractive_index_extraAxis = 0.0;
	} else if (isLossElement())
	{
		module_refractive_index_extraAxis = 0.0;
	} else if (isFilter())
	{
		module_refractive_index_extraAxis = 0.0;
	} else {
		module_refractive_index_extraAxis = 0.0;
	}
}

void Module::setRefInd_im(double n1)
{
	if (isCavity())
	{
		module_refractive_index_im = n1;
		cavity->setRefInd_im(n1);
	} else if (isTwoArmCavity())
	{
		module_refractive_index_im = n1;
		twoArmCavity->setRefInd_im(n1);
	} else if (isBirefringentCrystal())
	{
		module_refractive_index_im = n1;
		birefringentCrystal->setRefInd_im(n1);
	} else if (isBoundary())
	{
		module_refractive_index_im = 0.0;	
	} else if (isDevice())
	{
		module_refractive_index_im = 0.0;
	} else if (isTwoArmDevice())
	{
		module_refractive_index_im = 0.0;
	} else if (isTwoArmInterface())
	{
		module_refractive_index_im = n1;
		twoArmInterface->setRefInd_im(n1);
	} else if (isLossElement())
	{
		module_refractive_index_im = 0.0;
	} else if (isFilter())
	{
		module_refractive_index_im = 0.0;
	} else {
		cout << "Module::setRefInd_im(): Unknown module" << endl;
		exit(-1);
	}
}

void Module::setWidth(double newW)
{
	if (newW <= 0)
	{
		cout << "Module::setWidth: Must have width > 0 (= " << newW << ")" << endl;
	}
	
	if (isCavity())
	{
		cavity->setWidth(newW);
		module_width = newW;
	} else if (isBirefringentCrystal())
	{
		birefringentCrystal->setWidth(newW);
		module_width = newW;
	} else if (isTwoArmCavity())
	{
		twoArmCavity->setWidth(newW);
		module_width = newW;
	} else if (isTwoArmInterface())
	{
		twoArmInterface->setWidth(newW);
		module_width = newW;
	} else if (isBoundary())
	{
		module_width = 0.0;
	} else if (isDevice())
	{
		module_width = 0.0;
	} else if (isTwoArmDevice())
	{
		module_width = 0.0;
	} else if (isLossElement())
	{
		module_width = 0.0;
	} else if (isFilter())
	{
		module_width = 0.0;
	} else {
		cout << "Module::setWidth(): Unknown module" << endl;
		exit(-1);
	}
}
void Module::setAngle(double angle)
{
	if (angle>30||angle<0)
	{
		cout << "Module::setAngle: Angle="<<angle<<"must be between 0 and 30" << endl;
		cout << "Set angle to zero" << endl;
		angle = 0;
	}

	if (isCavity())
	{
		module_angle=angle;
	} else if (isBirefringentCrystal())
	{
		module_angle=angle;
	} else if (isTwoArmCavity())
	{
		module_angle=angle;
	}  else if (isTwoArmInterface())
	{
		twoArmInterface->setAngle(angle);
		module_angle=angle;
	}else if (isBoundary())
	{
		module_angle=angle;
	} else if (isDevice())
	{
		module_angle=angle;
	} else if (isTwoArmDevice())
	{
		module_angle=angle;
	} else if (isLossElement())
	{
		module_angle=angle;
	} else if (isFilter())
	{
		module_angle=angle;
	} else {
		cout << "Module::setAngle(): Unknown module" << endl;
		exit(-1);
	}
}

void Module::setIntAngle(double intAngle)
{
	if (intAngle>30||intAngle<0)
	{
		cout << "Module::setAngle: Angle="<<intAngle<<"must be between 0 and 30" << endl;
		cout << "Set intAngle to zero" << endl;
		intAngle = 0;
	}

	if (isCavity())
	{
		module_intAngle=intAngle;
	} else if (isBirefringentCrystal())
	{
		module_intAngle=intAngle;
	} else if (isTwoArmCavity())
	{
		module_intAngle=intAngle;
	}  else if (isTwoArmInterface())
	{
		twoArmInterface->setIntAngle(intAngle);
		module_intAngle=intAngle;
	}else if (isBoundary())
	{
		module_intAngle=intAngle;
	} else if (isDevice())
	{
		module_intAngle=intAngle;
	} else if (isTwoArmDevice())
	{
		module_intAngle=intAngle;
	} else if (isLossElement())
	{
		module_intAngle=intAngle;
	} else if (isFilter())
	{
		module_intAngle=intAngle;
	} else {
		cout << "Module::setIntAngle(): Unknown module" << endl;
		exit(-1);
	}
}

void Module::setPrevCav(double prev_cav)
{
	if (prev_cav<0)
	{
		cout << "Module::setPostCav: Must be positive"<< endl;
		exit(-1);
	}

	if (isCavity())
	{
		module_prev_cav=prev_cav;
	} else if (isBirefringentCrystal())
	{
		module_prev_cav=prev_cav;
	} else if (isTwoArmCavity())
	{
		module_prev_cav=prev_cav;
	}  else if (isTwoArmInterface())
	{
		twoArmInterface->setPrevCav(prev_cav);
		module_prev_cav=prev_cav;
	}else if (isBoundary())
	{
		module_prev_cav=prev_cav;
	} else if (isDevice())
	{
		module_prev_cav=prev_cav;
	} else if (isTwoArmDevice())
	{
		module_prev_cav=prev_cav;
	} else if (isLossElement())
	{
		module_prev_cav=prev_cav;
	} else if (isFilter())
	{
		module_prev_cav=prev_cav;
	} else {
		cout << "Module::setPreCav(): Unknown module" << endl;
		exit(-1);
	}
}

void Module::setPostCav(double post_cav)
{
	if (post_cav<0)
	{
		cout << "Module::setPostCav: Must be positive"<< endl;
		exit(-1);
	}

	if (isCavity())
	{
		module_post_cav=post_cav;
	} else if (isBirefringentCrystal())
	{
		module_post_cav=post_cav;
	} else if (isTwoArmCavity())
	{
		module_post_cav=post_cav;
	}  else if (isTwoArmInterface())
	{
		twoArmInterface->setPostCav(post_cav);
		module_post_cav=post_cav;
	}else if (isBoundary())
	{
		module_post_cav=post_cav;
	} else if (isDevice())
	{
		module_post_cav=post_cav;
	} else if (isTwoArmDevice())
	{
		module_post_cav=post_cav;
	} else if (isLossElement())
	{
		module_post_cav=post_cav;
	} else if (isFilter())
	{
		module_post_cav=post_cav;
	} else {
		cout << "Module::setPostCav(): Unknown module" << endl;
		exit(-1);
	}
}

void Module::setCosTh(double left, double right)
{
	if ((left == 0)||(right==0))
	{
		cout << "Module::setCosTh: Must have  cos(th) != 0 (= " << left << ", " << right << ")" << endl;
	}
	
	if (isCavity())
	{
		cavity->setCosTh(left, right);
		module_cos_th_left  = left;
		module_cos_th_right = right;
	} else if (isBirefringentCrystal())
	{
		birefringentCrystal->setCosTh(left, right);
		module_cos_th_left  = left;
		module_cos_th_right = right;
	} else if (isTwoArmCavity())
	{
		twoArmCavity->setCosTh(left, right);
		module_cos_th_left  = left;
		module_cos_th_right = right;
	}  else if (isTwoArmInterface())
	{
		twoArmInterface->setCosTh(left, right);
		module_cos_th_left  = left;
		module_cos_th_right = right;
	}else if (isBoundary())
	{
		module_cos_th_left = 1.0;
		module_cos_th_right = 1.0;
	} else if (isDevice())
	{
		module_cos_th_left = 1.0;
		module_cos_th_right = 1.0;
	} else if (isTwoArmDevice())
	{
		module_cos_th_left = 1.0;
		module_cos_th_right = 1.0;
	} else if (isLossElement())
	{
		module_cos_th_left = 1.0;
		module_cos_th_right = 1.0;
	} else if (isFilter())
	{
		module_cos_th_left = 1.0;
		module_cos_th_right = 1.0;
	} else {
		cout << "Module::setCosTh(): Unknown module" << endl;
		exit(-1);
	}
}

void Module::setToFileOutputKey(const std::string & newName)
{
	
	if (isCavity())
	{
		cavity->setToFileOutputKey(newName);
	} else if (isDevice())
	{
		device->setToFileOutputKey(newName);
	} else if (isBirefringentCrystal())
	{
		birefringentCrystal->setToFileOutputKey(newName);
	} else if (isTwoArmCavity())
	{
		twoArmCavity->setToFileOutputKey(newName);
	} else if (isTwoArmInterface())
	{
		twoArmInterface->setToFileOutputKey(newName);
	} else if (isTwoArmDevice())
	{
		twoArmDevice->setToFileOutputKey(newName);
	}

}

void Module::setPosition(double x0, double x1)
{
	if (x0 > x1)
	{
		cout << "Module::setPosition(): Must have x0 <= x1 (x1-x0= " << x1-x0 << ")" << endl;
		exit(-1);
	}
	
	if (isCavity())
	{
		if (x0 == x1)
		{
			cout << "Module::setPosition: Cannot have CAVITY of length 0 (x0 == x1)" << endl;
			exit(-1);
		}
		module_position_x0 = x0;
		module_position_x1 = x1;
		cavity->setPosition0(x0);
		cavity->setPosition1(x1);
		cavity->setWidth(x1-x0);
		
	} else if (isTwoArmCavity())
	{
		if (x0 == x1)
		{
			cout << "Module::setPosition: Cannot have CAVITY of length 0 (x0 == x1)" << endl;
			exit(-1);
		}
		module_position_x0 = x0;
		module_position_x1 = x1;
		twoArmCavity->setPosition0(x0);
		twoArmCavity->setPosition1(x1);
		twoArmCavity->setWidth(x1-x0);
		
	} else if (isBirefringentCrystal())
	{
		if (x0 == x1)
		{
			cout << "Module::setPosition: Cannot have CAVITY of length 0 (x0 == x1)" << endl;
			exit(-1);
		}
		module_position_x0 = x0;
		module_position_x1 = x1;
		birefringentCrystal->setPosition0(x0);
		birefringentCrystal->setPosition1(x1);
		birefringentCrystal->setWidth(x1-x0);
		
	} else if (isTwoArmInterface())
	{
		if (x0 == x1)
		{
			cout << "Module::setPosition: Cannot have CAVITY of length 0 (x0 == x1)" << endl;
			exit(-1);
		}
		module_position_x0 = x0;
		module_position_x1 = x1;
		twoArmInterface->setPosition0(x0);
		twoArmInterface->setPosition1(x1);
		twoArmInterface->setWidth(x1-x0);

	} else if (isBoundary())
	{
		if (x0 != x1)
		{
			cout << "Module::setPosition: Unclear where you want BOUNDARY (x0 != x1)" << endl;
			exit(-1);
		}
		
		module_position_x0 = x0;
		module_position_x1 = x0;
	} else if (isLossElement())
	{
		if (x0 != x1)
		{
			cout << "Module::setPosition: Unclear where you want LOSS ELEMENT (x0 != x1)" << endl;
			exit(-1);
		}
		
		module_position_x0 = x0;
		module_position_x1 = x0;
	} else if (isDevice())
	{
		if (x0 != x1)
		{
			cout << "Module::setPosition: Unclear where you want DEVICE (x0 != x1)" << endl;
			exit(-1);
		}
		
		module_position_x0 = x0;
		module_position_x1 = x0;
		device->setPosition(x0);
	} else if (isTwoArmDevice())
	{
		if (x0 != x1)
		{
			cout << "Module::setPosition: Unclear where you want DEVICE (x0 != x1)" << endl;
			exit(-1);
		}
		
		module_position_x0 = x0;
		module_position_x1 = x0;
		twoArmDevice->setPosition(x0);
	} else if (isFilter())
	{
		if (x0 != x1)
		{
			cout << "Module::setPosition: Unclear where you want FILTER (x0 != x1)" << endl;
			exit(-1);
		}
		
		module_position_x0 = x0;
		module_position_x1 = x0;
		//device->setPosition(x0);
	}  else {
		cout << "Module::setPosition(): Unknown module" << endl;
		exit(-1);
	}
}

Cavity * Module::addCavity()
{
	cavity = new Cavity();
	return cavity;
}

Cavity * Module::addCavity(double refInd, double width, double x0, double x1, double dt, int OUTPUT_TO_FILE)
{
	cavity = new Cavity(refInd,width,x0,x1,dt);
	setOutputToFile(OUTPUT_TO_FILE);
	setOutputToFileLevel(OUTPUT_TO_FILE);
	return cavity;
}

Device * Module::addDevice()
{
	device = new Device();
	return device;
}

TwoArmCavity * Module::addBirefringentCrystal()
{
	birefringentCrystal = new TwoArmCavity();
	return birefringentCrystal;
}

TwoArmCavity * Module::addTwoArmCavity()
{
	twoArmCavity = new TwoArmCavity();
	return twoArmCavity;
}

TwoArmInterface * Module::addTwoArmInterface()
{
	twoArmInterface = new TwoArmInterface();
	return twoArmInterface;
}

TwoArmCavity * Module::addTwoArmCavity(double refInd, double width, double x0, double x1, double dt, int OUTPUT_TO_FILE)
{
	twoArmCavity = new TwoArmCavity(refInd,width,x0,x1,dt);
	setOutputToFile(OUTPUT_TO_FILE);
	setOutputToFileLevel(OUTPUT_TO_FILE);
	return twoArmCavity;
}

TwoArmDevice * Module::addTwoArmDevice()
{
	twoArmDevice = new TwoArmDevice();
	return twoArmDevice;
}

Boundary * Module::addBoundary()
{
	boundary = new Boundary();
	return boundary;
}

Boundary * Module::addBoundary(double refCoeff)
{
	boundary = new Boundary(refCoeff);
	return boundary;
}

LossElement * Module::addLossElement()
{
	lossElement = new LossElement;
	return lossElement;
}

LossElement * Module::addLossElement(double loss_minus, double loss_plus)
{
	lossElement = new LossElement(loss_minus,loss_plus);
	return lossElement;
}

Filter * Module::addFilter(int numEl)
{
	filter = new Filter(numEl);
	return filter;
}

void Module::Print() const
{
	if (isCavity())
	{	
		cavity->Print();
	} else if (isDevice())
	{
		device->Print();
	} else if (isTwoArmCavity())
	{	
		twoArmCavity->Print();
	} else if (isTwoArmInterface())
	{	
		twoArmInterface->Print();
	} else if (isTwoArmDevice())
	{
		twoArmDevice->Print();
	}else if (isBoundary())
	{
		boundary->Print();
	} else if (isLossElement())
	{
		lossElement->Print();
	}else if (isFilter()) 
	{
		filter->Print();
	} else {
		cout << "Module::Print(): Unknown module !!!!" << endl;
		//exit(-1);
	}
	cout << " -> Output to file = " << getOutputToFile() << endl;
}


void Module::file_output_write()
{
	if (getOutputToFile() > 0)
	{
		if (isCavity())
		{
			if(cavity->getName()=="CAVL")
			{
			cavity->file_output_write_leftWall(getOutputToFileLevel());
			} else
			{
			cavity->file_output_write(getOutputToFileLevel());
			}
		} else if (isDevice())
		{
			#ifdef ITERATE_QW
			device->file_output_write(getOutputToFileLevel());
			#endif
		} else if (isTwoArmCavity())
		{
			if(twoArmCavity->getName()=="TACFRONT")
			{
			twoArmCavity->file_output_write_backWall(getOutputToFileLevel());
			} else
			{
			twoArmCavity->file_output_write(getOutputToFileLevel());
			}
		}  else if (isTwoArmInterface())
		{
			twoArmInterface->file_output_write(getOutputToFileLevel());
		}  else if (isBirefringentCrystal())
		{
			birefringentCrystal->file_output_write_backWall(getOutputToFileLevel());
		} else if (isTwoArmDevice())
		{
			#ifdef ITERATE_QW
			twoArmDevice->file_output_write(getOutputToFileLevel());
			#endif
		} 
 
	}
}

void Module::file_output_close()
{
	if (getOutputToFile() > 0)
	{
		if (isCavity())
		{
			cavity->file_output_close(getOutputToFileLevel());
		} else if (isDevice())
		{
			#ifdef ITERATE_QW
			device->file_output_close(getOutputToFileLevel());
			#endif
		} else if (isTwoArmCavity())
		{
			twoArmCavity->file_output_close(getOutputToFileLevel());
		} else if (isBirefringentCrystal())
		{
			birefringentCrystal->file_output_close(getOutputToFileLevel());
		}  else if (isTwoArmInterface())
		{
			twoArmInterface->file_output_close(getOutputToFileLevel());
		} else if (isTwoArmDevice())
		{
			#ifdef ITERATE_QW
			twoArmDevice->file_output_close(getOutputToFileLevel());
			#endif
		} 
	}
}

void Module::file_output_open(int out_count)
{
	if (getOutputToFile() > 0)
	{
		if (isCavity())
		{
			cavity->file_output_open(out_count,getOutputToFileLevel());
		} else if (isDevice())
		{
			#ifdef ITERATE_QW
			device->file_output_open(out_count,getOutputToFileLevel());
			#endif
		} else if (isBirefringentCrystal())
		{
			birefringentCrystal->file_output_open(out_count,getOutputToFileLevel());
		} else if (isTwoArmCavity())
		{
			twoArmCavity->file_output_open(out_count,getOutputToFileLevel());
		} else if (isTwoArmInterface())
		{
			twoArmCavity->file_output_open(out_count,getOutputToFileLevel());
		} else if (isTwoArmDevice())
		{
			#ifdef ITERATE_QW
			twoArmDevice->file_output_open(out_count,getOutputToFileLevel());
			#endif
		}
	}
}

void Module::Remove()
{
	//cout << "Module: Remove() -> ";
	if (isCavity())
	{
	//	cout << "Cavity" << endl;
		delete cavity;
		cavity = NULL;
	} else if (isDevice())
	{
	//	cout << "Device" << endl;
		delete device;
		device = NULL;
	} else if (isTwoArmCavity())
	{
	//	cout << "Cavity" << endl;
		delete twoArmCavity;
		twoArmCavity = NULL;
	} else if (isTwoArmInterface())
	{
	//	cout << "Cavity" << endl;
		delete twoArmInterface;
		twoArmInterface = NULL;
	} else if (isTwoArmDevice())
	{
	//	cout << "Device" << endl;
		delete twoArmDevice;
		twoArmDevice = NULL;
	}else if (isBoundary())
	{
	//	cout << "Boundary" << endl;
		delete boundary;
		boundary = NULL;
	} else if (isLossElement())
	{
	//	cout << "Loss Element" << endl;
		delete lossElement;
		lossElement = NULL;
	} else if (isFilter())
	{
	//	cout << "Filter" << endl;
		delete filter;
		filter = NULL;
	} else {
	//	cout << " nothing.. " << endl;
	}
}








