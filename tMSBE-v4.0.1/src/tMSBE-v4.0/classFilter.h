
#ifndef __FILTER_H_INCLUDED__
#define __FILTER_H_INCLUDED__

#include <string>
#include <complex>
//#include "classCyclic.h"
#include "classPrevField.h"
#include <cstdlib> // exit()
#include "largeVolumeStorage.h"

class Filter
{
	public:
		Filter();
		Filter(int);
		void Print( void ) const;		// Print all device data to screen
		void setFilterLength(int);
		int getFilterLength(void) const;
		void setFilterCoeff_pluss(std::complex<double>, int);
		void setFilterCoeff_minus(std::complex<double>, int);
		void setFilterCoeff_pluss_delay(std::complex<double>, int);
		std::complex<double> getFilterCoeff_pluss(int) const;
		std::complex<double> getFilterCoeff_minus(int) const;
		std::complex<double> getFilterCoeff_pluss_delay(int) const;
		void setEpluss(std::complex<double>*);
		void setEminus(std::complex<double>*);

		void setFilter_pluss_active_steps(int);
		void clearFields(void);

		void setFilter_pluss_diagnostic_square(double);
		void setFilter_pluss_doubleExp_flatTopWindow(double, double);
		void setFilter_pluss_doubleGauss_flatTopWindow(double, double, double,double);
		void setFilter_minus_PassThroughDelay(void);
		void setFilter_pluss_PassThroughDelay(void);
		std::complex<double> applyFilterPluss(std::complex<double> *);
		std::complex<double> applyFilterMinus(std::complex<double> *);
		std::complex<double> applyFilterPluss_cavity(std::complex<double> *);

		void setFilter_stage1_phase_diff(std::complex<double>);

		void Print_data(void);
		std::complex<double> getEminus(int);
		std::complex<double> getEpluss(int);

		void file_save_variables(int,int);
		void file_load_variables(int,int);
		
	private:
		int filter_length;
		std::complex<double> *filter_coeff_pluss;
		std::complex<double> *filter_coeff_minus;
		std::complex<double> *filter_coeff_pluss_delay;
		prev_field *electric_field_pluss;
		prev_field *electric_field_pluss_buffer; // Extra storage
		prev_field *electric_field_minus;

		largeVolumeStorage *cavity_field;
		int filter_coeff_pluss_cavity_field_max;
		int filter_coeff_pluss_cavity_field_counter;
		bool filter_coeff_pluss_cavity_field;

		bool filter_coeff_pluss_stage1; // Iterate through cavity signal
		bool filter_coeff_pluss_stage2; // Iterate through stored signal
		int  filter_coeff_pluss_rtt_steps; // How many steps for a full rtt
		int  filter_coeff_pluss_stage1_steps; // # iter to process cavity signal
		int  filter_coeff_pluss_stage1_counter;
		int  filter_coeff_pluss_stage2_counter;

		int  filter_coeff_pluss_2ndrun_counter;
		int  filter_coeff_pluss_2ndrun_max;

		bool filter_coeff_pluss_init_done;
		int filter_coeff_pluss_init_counter;
		int filter_coeff_pluss_init_max;
		
		int filter_counter;
		std::complex<double> filter_stage1_phase_diff;

		FILE *output_E;
		FILE *input_E;
};

inline void Filter::setFilter_stage1_phase_diff(std::complex<double> newP)
{
	filter_stage1_phase_diff = newP;
}


inline void Filter::setFilterLength(int numel)
{
	filter_length = numel;
}

inline int Filter::getFilterLength(void) const
{
	return filter_length;
}

inline void Filter::setFilterCoeff_pluss(std::complex<double> a,int num)
{
	filter_coeff_pluss[num] = a;
}
inline void Filter::setFilterCoeff_minus(std::complex<double> a,int num)
{
	filter_coeff_minus[num] = a;
}
inline void Filter::setFilterCoeff_pluss_delay(std::complex<double> a,int num)
{
	filter_coeff_pluss_delay[num] = a;
}

inline std::complex<double> Filter::getFilterCoeff_pluss(int num) const
{
	return filter_coeff_pluss[num];
}
inline std::complex<double> Filter::getFilterCoeff_minus(int num) const
{
	return filter_coeff_minus[num];
}
inline std::complex<double> Filter::getFilterCoeff_pluss_delay(int num) const
{
	return filter_coeff_pluss_delay[num];
}

inline void Filter::setEpluss(std::complex<double> *target)
{
	//electric_field_pluss->setFirstContainer(target);
	electric_field_pluss->set_top_element(target);
}

inline void Filter::setEminus(std::complex<double> *target)
{
	//electric_field_minus->setFirstContainer(target);
	electric_field_minus->set_top_element(target);
}

#endif












