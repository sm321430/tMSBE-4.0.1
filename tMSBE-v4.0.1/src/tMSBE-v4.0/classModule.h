

#ifndef __MODULE_H_INCLUDED__
#define __MODULE_H_INCLUDED__

#include <string>
#include "classCavity.h"
#include "classDevice.h"
#include "classTwoArmCavity.h"
#include "classTwoArmInterface.h"
#include "classTwoArmDevice.h"
#include "classBoundary.h"
#include "classLossElement.h"
#include "classFilter.h"

//! The MODULE class is a container for other objects
/*! A VECSEL is buildt up by multiple modules (Cavity(), Boundary, Device(),etc..) and the module class is the container that can hold any one of those objects. \n
    Some basic processing can be delegated to the module class, but most advanced commands should go directly to the respective components \n
    ex. by calling the module a->getCavity()->"myCommand"\n
    ex. by asking what kind of component is in the module a->isCavity()\n

    It is possible to extend this class to NEW objects, but you must duplicate the data in mutliple locations. See how the boundary, or device object is implemented and mimic that implementation. \n

    warning: Putting more than 1 object in a module can create unexpected results.\n

   \sa TwoArmCavity
   \sa BirefringentCrystal
   \sa Cavity
   \sa Boundary
   \sa Device
   \sa TwoArmDevice
   \sa Filter
   \sa LossElement \n

    2018: Additional comments\n
    Isak Kilen

   2020: Additional comments\n
   Sam McLaren
*/

class Module
{
	public:
		//! A constructor
		Module();
		//! A destructor
		~Module()
		{
			Remove();
		}
		//! A copy-constructor
		Module(const Module &obj)
		{
			module_position_x0 = obj.module_position_x0;
			module_position_x1 = obj.module_position_x1;
			module_width = obj.module_width;
			module_refractive_index = obj.module_refractive_index;
			module_refractive_index_extraAxis = obj.module_refractive_index_extraAxis;
			module_refractive_index_im = obj.module_refractive_index_im;
			module_cos_th_left = obj.module_cos_th_left;
			module_cos_th_right = obj.module_cos_th_right;
			output_to_file = obj.output_to_file;
			output_to_file_level = obj.output_to_file_level;

			
			cavity = NULL;
			device = NULL;
			twoArmCavity = NULL;
			twoArmDevice = NULL;
			birefringentCrystal = NULL;
			twoArmInterface = NULL;
			boundary = NULL;
			lossElement = NULL;
			filter = NULL;

			if (obj.isCavity())
			{
				cavity = new Cavity(*obj.cavity);
			} else if (obj.isDevice())
			{
				device = new Device(*obj.device);
			} else if (obj.isTwoArmInterface())
			{
				twoArmInterface = new TwoArmInterface(*obj.twoArmInterface);
			} else if (obj.isTwoArmCavity())
			{
				twoArmCavity = new TwoArmCavity(*obj.twoArmCavity);
			} else if (obj.isBirefringentCrystal())
			{
				birefringentCrystal = new TwoArmCavity(*obj.birefringentCrystal);
			} else if (obj.isTwoArmDevice())
			{
				twoArmDevice = new TwoArmDevice(*obj.twoArmDevice);
			} else if (obj.isBoundary())
			{
				boundary = new Boundary(*obj.boundary);
			} else if (obj.isLossElement())
			{
				lossElement = new LossElement(*obj.lossElement);
			} else if (obj.isFilter())
			{
				filter = new Filter(*obj.filter);
			} else {
				//cout << "Module::Module copy, dont recognize the module type.." << endl;
				//exit(-1);
			}
			
		}
		//! Calls the 'Print()' function from the containing object
		void Print( void ) const;

		//! Deletes the containing object
		void Remove(void);
		
		//! Returns true if containing object is an cavity
		/*! \return true or false */
		bool isCavity() const;
		//! Returns true if containing object is a device (QW,ABS,..)
		/*! \return true or false */
		bool isDevice() const;
		//! Returns true if containing object is a two arm cavity
		/*! \return true or false */
		bool isTwoArmCavity() const;
		//! Returns true if containing object is a birefringent crystal
		/*! \return true or false */
		bool isBirefringentCrystal() const;
		//! Returns true if containing object is a two arm interface
		/*! \return true or false */
		bool isTwoArmInterface() const;
		//! Returns true if containing object is a two arm device (QW,ABS,..)
		/*! \return true or false */
		bool isTwoArmDevice() const;
		//! Returns true if containing object is a boundary
		/*! \return true or false */
		bool isBoundary() const;
		//! Returns true if containing object is a loss element
		/*! \return true or false */
		bool isLossElement() const;
		//! Returns true if containing object is a filter
		/*! \return true or false */
		bool isFilter() const;
		
		//! Returns a pointer to a possible cavity object
		/*! Is NULL unless this module contains an object */
		Cavity * getCavity() const;
		//! Returns a pointer to a possible device object
		/*! Is NULL unless this module contains an object */
		Device * getDevice() const;
		//! Returns a pointer to a possible two arm cavity object
		/*! Is NULL unless this module contains an object */
		TwoArmCavity * getTwoArmCavity() const;
		//! Returns a pointer to a possible birefringent crystal object
		/*! Is NULL unless this module contains an object */
		TwoArmCavity * getBirefringentCrystal() const;
		//! Returns a pointer to a possible two arm interface object
		/*! Is NULL unless this module contains an object */
		TwoArmInterface * getTwoArmInterface() const;
		//! Returns a pointer to a possible two arm device object
		/*! Is NULL unless this module contains an object */
		TwoArmDevice * getTwoArmDevice() const;
		//! Returns a pointer to a possible boundary object
		/*! Is NULL unless this module contains an object */
		Boundary * getBoundary() const;
		//! Returns a pointer to a possible loss element object
		/*! Is NULL unless this module contains an object */
		LossElement * getLossElement() const;
		//! Returns a pointer to a possible filter object
		/*! Is NULL unless this module contains an object */
		Filter * getFilter() const;
		
		//! Adds a cavity object to the module
		Cavity * addCavity();
		//! Adds a cavity object to the module
		Cavity * addCavity(double, double, double, double, double, int);
		//! Adds a device object to the module
		Device * addDevice();
		//! Adds a two arm cavity object to the module
		TwoArmCavity * addTwoArmCavity();
		//! Adds a birefringent crystal object to the module
		TwoArmCavity * addBirefringentCrystal();
		//! Adds a two arm cavity object to the module
		TwoArmCavity * addTwoArmCavity(double, double, double, double, double, int);
		//! Adds a two arm interface object to the module
		TwoArmInterface * addTwoArmInterface();
		//! Adds a two arm interface object to the module
		TwoArmInterface * addTwoArmInterface(double, double, double, double, double, int);
		//! Adds a two arm device object to the module
		TwoArmDevice * addTwoArmDevice();
		//! Adds a boundary object to the module
		Boundary * addBoundary();
		//! Adds a boundary object to the module
		Boundary * addBoundary(double);
		//! Adds a loss element object to the module
		LossElement * addLossElement();
		//! Adds a loss element object to the module
		LossElement * addLossElement(double, double);
		//! Adds a filter object to the module
		Filter * addFilter(int);
		
		//! Return position of left object-boundary
		/*! The left and right boundaries can be equal if module contains QW, boundary,..*/
		double getPosition0() const;
		//! Return position of right object-boundary
		/*! The left and right boundaries can be equal if module contains QW, boundary,..*/
		double getPosition1() const;
		//! Return width of object
		/*! The width can be zero if module contains QW, boundary,..*/
		double getWidth() const;
		//! Return REAL refractive index from cavity object extraordinary axis
		double getRefInd_extraAxis() const;
		//! Return REAL refractive index from cavity object
		double getRefInd() const;
		//! Return IMAG refractive index from cavity object
		double getRefInd_im() const;
		//! Return cos(th) from cavity object
		void getCosTh(double *, double *) const;
		//! Set coordinate of left and right object-boundary.
		/*! Input depends on object: A cavity cannot have equal boundaries, but a device must have equal boundaries.*/
		void setPosition(double, double);
		//! Set the width of the containing object
		/*! Input depends on object: A cavity cannot have zero width, but a device will get zero width.*/
		void setWidth(double);
		//! Set the angle for the containing object
		void setAngle(double);
		//! Set the interference angle for the containing object
		void setIntAngle(double);
		//! Set the previous cavity for two arm cavities
		void setPrevCav(double);
		//! Set the post cavity for two arm cavities
		void setPostCav(double);
		//! Set REAL refractive index
		void setRefInd(double);
		//! Set extraordinary axis refractive index
		void setRefInd_extraAxis(double);
		//! Set IMAG refractive index
		void setRefInd_im(double);
		//! Set cos(th) for cavity
		void setCosTh(double, double);
		//! Set output file name for cavity or device
		void setToFileOutputKey(const std::string &);
		//! Set output file LEVEL for cavity or device

		//! Set output to file LEVEL
		/*! Each object handes its own output, but shuold include a reference to this type of variable where =0 gives no/minimal outout from the contained object.\n Higher levels gives more and more output.*/
		void setOutputToFileLevel(int flag);

		//! Set the variable 'output_to_file', indicating if there should be file-output from this module
		void  setOutputToFile(int);
		//! Return the variable 'output_to_file', indicating if there should be file-output from this module
		int  getOutputToFile(void) const;
		//! Return the variable 'output_to_file_level', indicating what LEVEL of output is requested.
		int  getOutputToFileLevel(void) const;
		
		//! Calls the corresponding file-output write function if output_to_file==1
		void file_output_write(void);

		//! Calls the corresponding file-output close function if output_to_file==1
		void file_output_close(void);

		//! Calls the corresponding file-output open function if output_to_file==1
		void file_output_open(int);
		
	
	private:
		//! A pointer to a possible cavity object, NULL if none is present
		Cavity *cavity;
		//! A pointer to a possible device object, NULL if none is present
		Device *device;
		//! A pointer to a possible two arm cavity object, NULL if none is present
		TwoArmCavity *twoArmCavity;
		//! A pointer to a possible birefringent crystal cavity object, NULL if none is present
		TwoArmCavity *birefringentCrystal;
		//! A pointer to a possible two arm interface object, NULL if none is present
		TwoArmInterface *twoArmInterface;
		//! A pointer to a possible two arm device object, NULL if none is present
		TwoArmDevice *twoArmDevice;
		//! A pointer to a possible boundary object, NULL if none is present
		Boundary *boundary;
		//! A pointer to a possible loss element object, NULL if none is present
		LossElement *lossElement;
		//! A pointer to a possible filter object, NULL if none is present
		Filter *filter;
		
		//! The left position of the obejct-boundary
		double module_position_x0;
		//! The right position of the obejct-boundary
		double module_position_x1;
		//! The width of the obejct-boundary
		double module_width;
		//! The REAL refractive index of the cavity, =0 if object is not a birefringent crystal cavity
		double module_refractive_index_extraAxis;
		//! The REAL refractive index of the cavity, =0 if object is not a cavity
		double module_refractive_index;
		//! The IMAG refractive index of the cavity, =0 if object is not a cavity
		double module_refractive_index_im;
		//! The angle factor cos(th) on left of the cavity, =1 if object is not a cavity
		double module_cos_th_left;
		//! The angle factor cos(th) on right of the cavity, =1 if object is not a cavity
		double module_cos_th_right;
	
		//! The angle of incidence incoming the cavity, =0 if not a two arm cavity
		double module_angle;
		//! The angle of incidence for interference fringes, =0 if not a two arm cavity
		double module_intAngle;
		//! The previous cavity to this module, =0 if not a two arm interface
		double module_prev_cav;
		//! The post cavity to this module, =0 if not a two arm interface
		double module_post_cav;
		//! Flag for output to files.
		/*! =1 if the containing object should output to file*/
		int output_to_file;
		//! The output to file level
		/*! Used to limit the amount of output to file data. Each object should handle their own output level*/
		int output_to_file_level;
};

inline void Module::setOutputToFile(int flag)
{
	output_to_file = flag;
}

inline int Module::getOutputToFile(void) const
{
	return output_to_file;
}

inline void Module::setOutputToFileLevel(int flag)
{
	output_to_file_level = flag;
}

inline int Module::getOutputToFileLevel(void) const
{
	return output_to_file_level;
}

inline double Module::getPosition0() const
{
	return module_position_x0;
}

inline double Module::getPosition1() const
{
	return module_position_x1;
}

inline double Module::getWidth() const
{
	return module_width;
}

inline double Module::getRefInd_extraAxis() const
{
	return module_refractive_index_extraAxis;
}

inline double Module::getRefInd() const
{
	return module_refractive_index;
}

inline double Module::getRefInd_im() const
{
	return module_refractive_index_im;
}

inline void Module::getCosTh(double *left, double *right) const
{
	*left =  module_cos_th_left;
	*right =  module_cos_th_right;
}

inline Cavity * Module::getCavity() const
{
	return cavity;
}

inline Device * Module::getDevice() const
{
	return device;
}

inline TwoArmCavity * Module::getTwoArmCavity() const
{
	return twoArmCavity;
}

inline TwoArmCavity * Module::getBirefringentCrystal() const
{
	return birefringentCrystal;
}

inline TwoArmInterface * Module::getTwoArmInterface() const
{
	return twoArmInterface;
}

inline TwoArmDevice * Module::getTwoArmDevice() const
{
	return twoArmDevice;
}

inline Boundary * Module::getBoundary() const
{
	return boundary;
}

inline LossElement * Module::getLossElement() const
{
	return lossElement;
}

inline Filter* Module::getFilter() const
{
	return filter;
}
inline bool Module::isCavity() const
{
	if (cavity == NULL)
	{
		return false;
	} else {
		return true;
	}
}

inline bool Module::isDevice() const
{
	if (device == NULL)
	{
		return false;
	} else {
		return true;
	}
}

inline bool Module::isTwoArmCavity() const
{
	if (twoArmCavity == NULL)
	{
		return false;
	} else {
		return true;
	}
}

inline bool Module::isBirefringentCrystal() const
{
	if (birefringentCrystal == NULL)
	{
		return false;
	} else {
		return true;
	}
}

inline bool Module::isTwoArmInterface() const
{
	if (twoArmInterface == NULL)
	{
		return false;
	} else {
		return true;
	}
}

inline bool Module::isTwoArmDevice() const
{
	if (twoArmDevice == NULL)
	{
		return false;
	} else {
		return true;
	}
}

inline bool Module::isBoundary() const
{
	if (boundary == NULL)
	{
		return false;
	} else {
		return true;
	}
}

inline bool Module::isLossElement() const
{
	if (lossElement == NULL)
	{
		return false;
	} else {
		return true;
	}
}

inline bool Module::isFilter() const
{
	if (filter == NULL)
	{
		return false;
	} else {
		return true;
	}
}

#endif
