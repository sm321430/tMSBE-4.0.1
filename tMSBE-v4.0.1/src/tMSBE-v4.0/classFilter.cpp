
/*
	The FILTER class is a frequency filter that should be put between two layers of indentical refractive index
*/

#include <stdlib.h>
#include <iostream>
#include "classFilter.h"
#include "constantsAndMiscUnits.h"
#include <stdio.h> // To write coeff to file

using namespace std;

void Filter::clearFields()
{
	if (electric_field_pluss != NULL)
	{
		electric_field_pluss->zero_all_elements();
	}
	if (electric_field_pluss_buffer != NULL)
	{
		electric_field_pluss_buffer->zero_all_elements();
	}
	if (electric_field_minus != NULL)
	{
		electric_field_minus->zero_all_elements();
	}

	filter_coeff_pluss_stage1 = false;
	filter_coeff_pluss_stage2 = false;
	filter_coeff_pluss_stage1_counter = 0;
	filter_coeff_pluss_stage2_counter = 0;

	filter_coeff_pluss_2ndrun_counter = 0;

	filter_coeff_pluss_init_done = false;
	filter_coeff_pluss_init_counter = 0;
	
	filter_counter = 0;
	filter_stage1_phase_diff = 1.0;


	
	filter_coeff_pluss_cavity_field_counter = 0;
	filter_coeff_pluss_cavity_field = false;
	if (cavity_field != NULL)
	{
		cavity_field->zero_all_containers();
		cavity_field->seq_iterate_reset();
	}
}

/*
	Construction function
*/

Filter::Filter()
{
	filter_length = 0;
	filter_coeff_pluss = NULL;
	filter_coeff_minus = NULL;
	filter_coeff_pluss_delay = NULL;
	electric_field_pluss = NULL;
	electric_field_minus = NULL;
	electric_field_pluss_buffer = NULL;

	cavity_field = NULL;
	filter_coeff_pluss_cavity_field_max = -1;

	filter_coeff_pluss_rtt_steps = -1;
	filter_coeff_pluss_stage1_steps = -1; // # iter to process cavity signal
	filter_coeff_pluss_2ndrun_max = -1;
	filter_coeff_pluss_init_max = -1;

	clearFields();

//	output_E = fopen("out__filter.dat","w+");
}

Filter::Filter(int numEl)
{
//	output_E = fopen("out__filter_out.dat","w+");
//	input_E  = fopen("out__filter_in.dat","w+");

	filter_coeff_pluss_rtt_steps = -1;
	filter_coeff_pluss_stage1_steps = -1; // # iter to process cavity signal
	filter_coeff_pluss_2ndrun_max = -1;
	filter_coeff_pluss_init_max = -1;
	filter_coeff_pluss_cavity_field_max = -1;


	cavity_field = NULL;


	// Force filter length to be ODD
	if (numEl%2 == 0)
	{
		cout << "Filter:: Forcing filter length to be ODD" << endl;
		cout << "Filter:: New length is = " << numEl+1 << endl;
		numEl += 1;
	}
	setFilterLength(numEl);

	filter_coeff_pluss = NULL;
	filter_coeff_pluss = new std::complex<double>[getFilterLength()];
	if (filter_coeff_pluss==NULL)
	{
		cout << "FILTER():: Failed to allocate memory for filter_coeff_pluss" << endl;
		exit(-1);
	}
	for(int i = 0; i < getFilterLength(); i++)
	{
		filter_coeff_pluss[i] = 0;
	}
	filter_coeff_pluss_delay = NULL;
	filter_coeff_pluss_delay = new std::complex<double>[getFilterLength()];
	if (filter_coeff_pluss_delay==NULL)
	{
		cout << "FILTER():: Failed to allocate memory for filter_coeff_pluss_delay" << endl;
		exit(-1);
	}
	for(int i = 0; i < getFilterLength(); i++)
	{
		filter_coeff_pluss_delay[i] = 0;
	}
	filter_coeff_minus = NULL;
	filter_coeff_minus = new std::complex<double>[getFilterLength()];
	if (filter_coeff_minus==NULL)
	{
		cout << "FILTER():: Failed to allocate memory for filter_coeff_minus" << endl;
		exit(-1);
	}
	for(int i = 0; i < getFilterLength(); i++)
	{
		filter_coeff_minus[i] = 0;
	}
	
	electric_field_pluss = NULL;
	electric_field_minus = NULL;

	electric_field_pluss = new prev_field(getFilterLength());
	electric_field_minus = new prev_field(getFilterLength());
	if ((electric_field_pluss==NULL)||(electric_field_minus==NULL))
	{
		cout << "FILTER():: Failed to allocate memory for electric field" << endl;
		exit(-1);
	}
	electric_field_pluss->zero_all_elements();
	electric_field_minus->zero_all_elements();


	electric_field_pluss_buffer = NULL;
	electric_field_pluss_buffer = new prev_field(getFilterLength());
	if (electric_field_pluss_buffer==NULL)
	{
		cout << "FILTER():: Failed to allocate memory for electric field buffer" << endl;
		exit(-1);
	}
	electric_field_pluss_buffer->zero_all_elements();


	clearFields();
}


void Filter::setFilter_pluss_diagnostic_square(double ratio)
{
	electric_field_pluss->zero_all_elements();
	if ((ratio < 0)||(ratio > 1))
	{
		cout << "setFilter_pluss_diagnostic_square(): Need ratio in [0,1]" << endl;
		exit(-1);
	}
	
	int numZeros = floor(getFilterLength()*(1-ratio)/2.0);
	for(int i = numZeros; i < getFilterLength()-numZeros; i++)
	{
		std::complex<double> coeff = 1.0;
		setFilterCoeff_pluss(coeff, i);
	}

	// Output coefficients
	FILE *fid = fopen("filter_coeff_pluss.dat","w+");
	for(int i =0; i < getFilterLength(); i++)
	{
		std::complex<double> cf = getFilterCoeff_pluss(i);
		fprintf(fid,"%.16e\t%.16e\n",real(cf),imag(cf));
	}
	fclose(fid);
}

/*
	Create the filter: 
		h(n) = exp(-i*wa_s*n) + exp(-i*wb_s*n)
	where wa_s and wb_s are sample frequencies found by
		wa_s = wa*dt
		wb_s = wb*dt
	where wa is the frequency [1/s] and dt is the sampling timestep [s]
*/
void Filter::setFilter_pluss_doubleExp_flatTopWindow(double wa_s, double wb_s)
{
	electric_field_pluss->zero_all_elements();
	// Flat Top parameters
	double a0 = 1.0;
	double a1 = 1.93;
	double a2 = 1.29;
	double a3 = 0.388;
	double a4 = 0.028;

	for(int i = 0; i < getFilterLength(); i++)
	{
		// Double Exp Filter
		double n = (double)i - (getFilterLength()-1.0)/2.0;

		std::complex<double> coeff = exp(-I*n*wa_s);
		if (wb_s != 0)
		{
			coeff += exp(-I*n*wb_s);
		}

		// Flat Top Window
		n = ((double)i )/((double)getFilterLength()-1.0);
		double wn = a0 - a1*cos(2.0*Pi*n) + a2*cos(4.0*Pi*n) - a3*cos(6.0*Pi*n) + a4*cos(8.0*Pi*n);

		coeff *= wn;
		setFilterCoeff_pluss(coeff, i);
	}

	// NORMALIZE for frequency respose at wa
	std::complex<double> norm = 0.0;
	for(int i = 0; i < getFilterLength(); i++)
	{
		norm += getFilterCoeff_pluss(i)*exp(I*wa_s*((double)i));
	}
	//norm *= sqrt(getFilterLength());

	for(int i = 0; i < getFilterLength(); i++)
	{
		std::complex<double> cf = getFilterCoeff_pluss(i)/norm;
		setFilterCoeff_pluss(cf,i);
	}
}

/*
	wa_s = wa*dt;
	wb_s = wb*dt;
	width_spec_s = width_spec*dt
	w0_s = w0*dt;
*/

void Filter::setFilter_pluss_doubleGauss_flatTopWindow(double wa_s, double wb_s, double width_spec_s, double w0_s)
{
	electric_field_pluss->zero_all_elements();
	// Flat Top parameters
	double a0 = 1.0;
	double a1 = 1.93;
	double a2 = 1.29;
	double a3 = 0.388;
	double a4 = 0.028;

	//double gauss_width = 2.355/(width_spec_s); // Width of abs(S)
	double gauss_width = 0.5*sqrt(2.0)*2.355/(width_spec_s); // Width of abs(S)

	for(int i = 0; i < getFilterLength(); i++)
	{
		// Double Exp Filter
		double n = ((double)i - ((double)getFilterLength()-1.0)/2.0);

		double den = 2.0*gauss_width*gauss_width;
		std::complex<double> coeff = exp(-I*n*wa_s -n*n/den);
		if (wb_s != 0)
		{
			coeff += exp(-I*n*wb_s -n*n/den);
		}

		// Flat Top Window
	//	n = ((double)i )/((double)getFilterLength()-1.0);
	//	double wn = a0 - a1*cos(2.0*Pi*n) + a2*cos(4.0*Pi*n) - a3*cos(6.0*Pi*n) + a4*cos(8.0*Pi*n);
	//	coeff *= wn;

		setFilterCoeff_pluss(coeff, i);
	}

	/*
		BUILD IN SIGNAL TIME DELAY INTO COEFF
	*/
/*	
	for(int i = 0; i < getFilterLength(); i++)
	{
		//double n = ((double)i((double)getFilterLength()-1.0)/2.0);
		std::complex<double> phase = exp(I*((double)i)*w0_s);
		std::complex<double> coff = getFilterCoeff_pluss(i)*phase;
		setFilterCoeff_pluss(coff,i);
	}
*/

	/*
		NORMALIZE for frequency respose at: wa + wb
	*/
	for(int i = 0; i < 2*getFilterLength(); i++)
	{
		std::complex<double> coff = exp(-I*wa_s*((double)(i)));
		electric_field_pluss->update_head_index();
		electric_field_pluss->set_top_element(&coff);
	//	electric_field_pluss->set_element(i,&coff);
	}
	std::complex<double> norm = electric_field_pluss->c_array_conv(filter_coeff_pluss);
	cout << "setFilter()::Gauss Norm = " << norm << " |Norm| = " << abs(norm) << endl;
	//exit(-1);
	
	for(int i = 0; i < getFilterLength(); i++)
	{
		std::complex<double> coff = getFilterCoeff_pluss(i)/norm;
		setFilterCoeff_pluss(coff,i);
	}
	electric_field_pluss->zero_all_elements();

	/* 
		Gaussian filter decays very fast, can exclude many points
	*/

	double max=0.0;

	for(int i = 0; i < getFilterLength(); i++)
	{
		double fc = abs(getFilterCoeff_pluss(i));
		if (fc>max)
		{
			max = fc;
		}
	}
	if (abs(getFilterCoeff_pluss(0))/max >1e-12)
	{
		cout << "Filter():: WARNING Filter coefficient at edges are too high, make filter longer" << endl;
		cout << "Filter():: c(0)/max = " << abs(getFilterCoeff_pluss(0))/max << endl;
	}

	double cutoff = (1.0e-16)*max;
	std::complex<double> cf;
/*
	for(int i = 0; i < getFilterLength(); i++)
	{
		cf = getFilterCoeff_pluss(i);
		electric_field_pluss->set_element(i,&cf);
	}
	electric_field_pluss->c_array_conv_ind_set(cutoff);
	//electric_field_pluss->c_array_conv_ind_set_center_volume(1.0); 
	for(int i = 0; i < getFilterLength(); i++)
	{
		cf = 0.0;
		electric_field_pluss->set_element(i,&cf);
	}
*/
	// Output coefficients
	FILE *fid = fopen("filter_coeff_pluss.dat","w+");
	for(int i =0; i < getFilterLength(); i++)
	{
		cf = getFilterCoeff_pluss(i);
		fprintf(fid,"%.16e\t%.16e\n",real(cf),imag(cf));
	}
	fclose(fid);

	// SET DELAY FILTER COEFFICIENTS
	for(int i = 0; i < getFilterLength(); i++)
	{
		// Double Exp Filter
		double n = ((double)i - ((double)getFilterLength()-1.0)/2.0);

		std::complex<double> coeff = 0.0;
		if (n==0)
		{
			coeff = 1.0;
		} else {
			coeff = sin(n)/n;
		}


		// Flat Top Window
		n = ((double)i )/((double)getFilterLength()-1.0);
		double wn = a0 - a1*cos(2.0*Pi*n) + a2*cos(4.0*Pi*n) - a3*cos(6.0*Pi*n) + a4*cos(8.0*Pi*n);
		coeff *= wn;

		setFilterCoeff_pluss_delay(coeff, i);
	}

	/*
		NORMALIZE DELAY for frequency respose at: w0
	*/
	norm = 0.0;
	for(int i = 0; i < getFilterLength(); i++)
	{
		norm += getFilterCoeff_pluss_delay(i);
	}
	for(int i = 0; i < getFilterLength(); i++)
	{
		std::complex<double> coff = getFilterCoeff_pluss_delay(i)/norm;
		setFilterCoeff_pluss_delay(coff,i);
	}
	// Output coefficients
	fid = fopen("filter_coeff_pluss_delay.dat","w+");
	for(int i =0; i < getFilterLength(); i++)
	{
		cf = getFilterCoeff_pluss_delay(i);
		fprintf(fid,"%.16e\t%.16e\n",real(cf),imag(cf));
	}
	fclose(fid);
}

void Filter::setFilter_minus_PassThroughDelay()
{
	electric_field_minus->zero_all_elements();
	std::complex<double> coeff = 1.0;
	setFilterCoeff_minus(coeff, (getFilterLength()-1)/2); // Half delay
	electric_field_minus->c_array_conv_ind_set((getFilterLength()-1)/2,(getFilterLength()-1)/2);

	// Output coefficients
	FILE *fid = fopen("filter_coeff_minus.dat","w+");
	for(int i =0; i < getFilterLength(); i++)
	{
		std::complex<double> cf = getFilterCoeff_minus(i);
		fprintf(fid,"%.16e\t%.16e\n",real(cf),imag(cf));
	}
	fclose(fid);
}

void Filter::setFilter_pluss_PassThroughDelay()
{
	electric_field_pluss->zero_all_elements();
	std::complex<double> coeff = 1.0;
	setFilterCoeff_pluss(coeff, (getFilterLength()-1)/2); // Half delay
	electric_field_pluss->c_array_conv_ind_set((getFilterLength()-1)/2,(getFilterLength()-1)/2);

	// Output coefficients
	FILE *fid = fopen("filter_coeff_pluss.dat","w+");
	for(int i =0; i < getFilterLength(); i++)
	{
		std::complex<double> cf = getFilterCoeff_pluss(i);
		fprintf(fid,"%.16e\t%.16e\n",real(cf),imag(cf));
	}
	fclose(fid);
}

// Apply filter to the fields, get filtered value
std::complex<double> Filter::applyFilterPluss(std::complex<double> *Ft)
{
//	filter_coeff_pluss_2ndrun_counter = filter_coeff_pluss_2ndrun_max; // To ship init phase
	std::complex<double> yn;

	// INIT PHASE
	if (filter_coeff_pluss_init_counter<=filter_coeff_pluss_init_max)
	{
		// Update filter
		electric_field_pluss->update_head_index(); // Cycle one timestep forward
		setEpluss(Ft); // Insert new data

		filter_coeff_pluss_init_counter += 1;
		//yn = electric_field_pluss->get_element(getFilterLength()-1); //Get only the oldest filter element
		yn = electric_field_pluss->get_element((getFilterLength()-1)/2); // Get middle filter element
		
		filter_counter++;
		return yn;
	} else if (!filter_coeff_pluss_init_done) {
		
		printf("[%d] FILTER SWITCH : INIT DONE: start waiting\n",filter_counter);
		filter_coeff_pluss_init_done = true;
	}

	// DELETE
/*
	electric_field_pluss->update_head_index(); // Cycle one timestep forward
	setEpluss(Ft); // Insert new data

	// Filter cavity signal
	yn = electric_field_pluss->c_array_conv(filter_coeff_pluss);
	fprintf(output_E,"%.16e\t%.16e\n",real(yn),imag(yn));
	yn = electric_field_pluss->get_element((getFilterLength()-1)/2);
	return yn;
*/
	// DELETE


	// RUN PHASE
	if (filter_coeff_pluss_stage1)
	{
		// Phase difference from propagation
/*
		if (filter_coeff_pluss_stage1_counter==filter_coeff_pluss_stage1_steps-1)
		{
			std::complex<double> newEl = electric_field_pluss->get_element(0);
			std::complex<double> oldEl = electric_field_pluss_buffer->get_element((getFilterLength()+1)/2);
			// Should be same element
			
			//filter_stage1_phase_diff = exp(I*(arg(newEl) - arg(oldEl)));
			filter_stage1_phase_diff = newEl/oldEl;
			cout << "filter():: phase = " << abs(filter_stage1_phase_diff) << " exp(I " << arg(filter_stage1_phase_diff) << " )" << endl;
	//		cout << "array[0]   = " << electric_field_pluss->get_element(0) << endl;
	//		cout << "buffer[mid] = " << electric_field_pluss_buffer->get_element((getFilterLength()+1)/2) << endl;
		}
*/

		// Update filter
		electric_field_pluss->update_head_index(); // Cycle one timestep forward
		setEpluss(Ft); // Insert new data

		// Filter cavity signal
		yn = electric_field_pluss->c_array_conv(filter_coeff_pluss);

		filter_coeff_pluss_stage1_counter+=1;	
		if (filter_coeff_pluss_stage1_counter >= filter_coeff_pluss_stage1_steps)
		{
			filter_coeff_pluss_stage1_counter = 0;
			filter_coeff_pluss_stage1 = false;
			filter_coeff_pluss_stage2 = true;

			printf("[%d] FILTER SWITCH : from stage1 to stage2\n",filter_counter);
		} 
		
	} else if (filter_coeff_pluss_stage2)
	{
		// Use buffered signal to process last part
		// Send new signal into buffer
		// Once done, switch buffer and regular signal

		// Update filter
		electric_field_pluss->update_head_index(); // Cycle one timestep forward
		electric_field_pluss_buffer->update_head_index(); // Cycle buffer one timestep forward
		std::complex<double> tmp1 = electric_field_pluss_buffer->get_element((getFilterLength()-1)/2)*filter_stage1_phase_diff;
		setEpluss(&tmp1); // Insert buffer data into filter
		electric_field_pluss_buffer->set_top_element(Ft); // Move cavity field into buffer
		
		// Filter cavity signal
		yn = electric_field_pluss->c_array_conv(filter_coeff_pluss);

		filter_coeff_pluss_stage2_counter+=1;
		if (filter_coeff_pluss_stage2_counter >= (getFilterLength()+1)/2)
		{
			filter_coeff_pluss_stage2_counter = 0;
			filter_coeff_pluss_stage1 = false;
			filter_coeff_pluss_stage2 = false;

			printf("[%d] FILTER SWITCH: from stage2 to waiting\n",filter_counter);
			

			electric_field_pluss->copy(electric_field_pluss_buffer); // Copy buffer into cavity field

/*
			FILE *fid = fopen("out_filter_Epluss_start_step3.dat","w+");
			for(int i =0; i < getFilterLength(); i++)
			{
				std::complex<double> cf = electric_field_pluss->get_element(i);
				fprintf(fid,"%.16e\t%.16e\n",real(cf),imag(cf));
			}
			fclose(fid);
*/
		}

	} else {
		// Update filter
		electric_field_pluss->update_head_index(); // Cycle one timestep forward
		setEpluss(Ft); // Insert new data

		// No filtering, wait
		//yn = electric_field_pluss->get_element(getFilterLength()-1);
		yn = electric_field_pluss->get_element((getFilterLength()-1)/2); // Half delay
		//yn = electric_field_pluss->c_array_conv(filter_coeff_pluss_delay);

		filter_coeff_pluss_2ndrun_counter+=1;
		if (filter_coeff_pluss_2ndrun_counter >= filter_coeff_pluss_2ndrun_max)
		{
			filter_coeff_pluss_2ndrun_counter = 0;
			filter_coeff_pluss_stage1 = true;
			filter_coeff_pluss_stage2 = false;

			printf("[%d] FILTER SWITCH: from waiting to stage1\n",filter_counter);

			// Copy field to buffer
			electric_field_pluss_buffer->copy(electric_field_pluss);
/*
			FILE *fid = fopen("out_filter_Epluss_start_step1.dat","w+");
			for(int i =0; i < getFilterLength(); i++)
			{
				std::complex<double> cf = electric_field_pluss->get_element(i);
				fprintf(fid,"%.16e\t%.16e\n",real(cf),imag(cf));
			}
			fclose(fid);
*/
		} 
		
	}
	
	filter_counter++;
	
	return yn;
}

std::complex<double> Filter::applyFilterPluss_cavity(std::complex<double> *Ft)
{
//	filter_coeff_pluss_2ndrun_counter = filter_coeff_pluss_2ndrun_max; // To ship init phase
	std::complex<double> yn;

	//--------------------------
	// Update filter
	electric_field_pluss->update_head_index(); // Cycle one timestep forward
	setEpluss(Ft); // Insert new data

	yn = electric_field_pluss->c_array_conv(filter_coeff_pluss);
//	cavity_field->seq_set_iterate(yn);
	
	return yn;
/*
	//-----------------------

	// INIT PHASE
	if (filter_coeff_pluss_init_counter<=filter_coeff_pluss_init_max)
	{
		// Update filter
		electric_field_pluss->update_head_index(); // Cycle one timestep forward
		setEpluss(Ft); // Insert new data

		filter_coeff_pluss_init_counter += 1;
		//yn = electric_field_pluss->get_element(getFilterLength()-1); //Get only the oldest filter element
		yn = electric_field_pluss->get_element((getFilterLength()-1)/2); // Get middle filter element
		
		filter_counter++;
		return yn;
	} else if (!filter_coeff_pluss_init_done) {
		
		printf("[%d] FILTER SWITCH : INIT DONE: start waiting\n",filter_counter);
		filter_coeff_pluss_init_done = true;
	}

	// RUN PHASE
	if (filter_coeff_pluss_stage1)
	{
		// Update filter
		electric_field_pluss->update_head_index(); // Cycle one timestep forward
		setEpluss(Ft); // Insert new data

		// Filter cavity signal
		yn = electric_field_pluss->c_array_conv(filter_coeff_pluss);
	//	fprintf(output_E,"%.16e\t%.16e\n",real(yn),imag(yn));
	//	fprintf(input_E,"%.16e\t%.16e\n",real(*Ft),imag(*Ft));
		cavity_field->seq_set_iterate(yn);
	//	yn = electric_field_pluss->get_element((getFilterLength()-1)/2);



		filter_coeff_pluss_stage1_counter+=1;	
		if (filter_coeff_pluss_stage1_counter == filter_coeff_pluss_cavity_field_max)
		{
			filter_coeff_pluss_stage1_counter = 0;
			filter_coeff_pluss_stage1 = false;
			filter_coeff_pluss_stage2 = true;

			std::complex<double> zn = cavity_field->get_element(0);
			yn = 0.5*(abs(yn)+abs(zn))*exp(I*0.5*(arg(yn)+arg(zn)));

	//		fclose(output_E);
	//		exit(-1);
			cavity_field->seq_iterate_reset();
//
//			for(int i = 0; i < (getFilterLength()-1)/4;i++)
//			{
//				cavity_field->seq_iterate();
//			}
//
	//		printf("[%d] FILTER SWITCH : from stage1 to stage2\n",filter_counter);
		} 
		
	} else if (filter_coeff_pluss_stage2)
	{
		// Update filter
		electric_field_pluss->update_head_index(); // Cycle one timestep forward
		setEpluss(Ft); // Insert new data

		// Filter cavity signal
		yn = cavity_field->seq_get_iterate();

//		if (filter_coeff_pluss_stage2_counter == 0)
//		{
//			std::complex<double> zn = electric_field_pluss->get_element((getFilterLength()-1)/2);
//			yn = 0.5*(abs(yn)+abs(zn))*exp(I*0.5*(arg(yn)+arg(zn)));
//		}



		filter_coeff_pluss_stage2_counter+=1;
		if (filter_coeff_pluss_stage2_counter == filter_coeff_pluss_cavity_field_max)
		{
//			std::complex<double> newP = 0.5*(abs(yn)+abs(*Ft))*exp(I*0.5*(arg(yn)+arg(*Ft)));
//			setEpluss(&newP); // Insert last point
			
			filter_coeff_pluss_stage2_counter = 0;
			filter_coeff_pluss_stage1 = false;
			filter_coeff_pluss_stage2 = false;

			cavity_field->seq_iterate_reset();

	//		printf("[%d] FILTER SWITCH: from stage2 to waiting\n",filter_counter);
		}

	} else {
		// Update filter
		electric_field_pluss->update_head_index(); // Cycle one timestep forward
		setEpluss(Ft); // Insert new data

		// No filtering, wait
		//yn = electric_field_pluss->get_element(getFilterLength()-1);
		yn = electric_field_pluss->get_element((getFilterLength()-1)/2); // Half delay
		//yn = electric_field_pluss->c_array_conv(filter_coeff_pluss_delay);

		filter_coeff_pluss_2ndrun_counter+=1;
		if (filter_coeff_pluss_2ndrun_counter >= filter_coeff_pluss_2ndrun_max)
		{
			filter_coeff_pluss_2ndrun_counter = 0;
			filter_coeff_pluss_stage1 = true;
			filter_coeff_pluss_stage2 = false;

	//		printf("[%d] FILTER SWITCH: from waiting to stage1\n",filter_counter);
		} 
		
	}
	
	filter_counter++;
	
	return yn;
*/
}


std::complex<double> Filter::applyFilterMinus(std::complex<double> *Ft)
{
	// Set
	//electric_field_pluss->updateStorage(); // Cycle one timestep forward
	electric_field_minus->update_head_index(); // Cycle one timestep forward
	setEminus(Ft); // Insert new data

	// Apply filter
	//std::complex<double> yn = electric_field_minus->c_array_conv(filter_coeff_minus);
	std::complex<double> yn = electric_field_minus->get_element((getFilterLength()-1)/2);
	return yn;
}

void Filter::Print() const
{
	cout << "Print filter:" << endl;
	cout << " -> filter length  = " << getFilterLength() << endl;
	cout << " -> step1 phase    = " << abs(filter_stage1_phase_diff) << " exp(I " << arg(filter_stage1_phase_diff) << " )" << endl;
}


void Filter::Print_data()
{
	printf("array[]:  ");
	electric_field_pluss->Print_array();
	printf("\n");
	printf("buffer[]: ");
	electric_field_pluss_buffer->Print_array();
}

std::complex<double> Filter::getEminus(int num)
{
	return electric_field_minus->get_element(num);
}
std::complex<double> Filter::getEpluss(int num)
{
	return electric_field_pluss->get_element(num);
}

/*
	Tell the filter how large of a cavity we have.
	ext_rtt_steps is the FULL round trip time of the cavity NOT including the 2*(filter_length)
		This is the number of iterations for a field component to propagate from a position in the caivity, make a full round trip and end back at the same point WITHOUT THE FILTER.

*/
void Filter::setFilter_pluss_active_steps(int ext_rtt_steps)
{
	if (ext_rtt_steps <= 0)
	{
		cout << "setFilter_pluss_active(): No cavity has negative or zero # steps.." << endl;
		cout << "Trying to set to " << ext_rtt_steps << endl;
		exit(-1);
	}
	filter_coeff_pluss_rtt_steps    = ext_rtt_steps +  getFilterLength()-1; // How many steps for a full rtt
	filter_coeff_pluss_stage1_steps = ext_rtt_steps + (getFilterLength()-1)/2; // Stage1Steps
	
	//filter_coeff_pluss_2ndrun_max 	= (ext_rtt_steps + getFilterLength()-1);    // Wait time between running the filter
	filter_coeff_pluss_2ndrun_max 	= 1*(ext_rtt_steps + getFilterLength()-1);    // Wait time between running the filter
	filter_coeff_pluss_init_max 	= (ext_rtt_steps + getFilterLength()-1); 	// Initial wait time before starting

	cout << "setFilter_pluss_active_steps:: Wait time = " << filter_coeff_pluss_2ndrun_max << " steps, " << filter_coeff_pluss_2ndrun_max/filter_coeff_pluss_rtt_steps << " [# roundtrips]" << endl;

	// STORAGE FOR CONTAINERS
	if (cavity_field != NULL)
	{
		delete [] cavity_field;
	}
	filter_coeff_pluss_cavity_field_max = filter_coeff_pluss_rtt_steps;
	cavity_field = new largeVolumeStorage(filter_coeff_pluss_cavity_field_max); // Extra storage space
	
}

void Filter::file_save_variables(int save_count, int filter_nr)
{
	// Save all dynamic variables
	std::stringstream fileName;
	
	fileName.str("");
	fileName << "save/save_" << save_count << "_FILTER_" << filter_nr << "_Ep";
	electric_field_pluss->saveArray(fileName.str());

	fileName.str("");
	fileName << "save/save_" << save_count << "_FILTER_" << filter_nr << "_Ep_buff";
	electric_field_pluss_buffer->saveArray(fileName.str());
	
	fileName.str("");
	fileName << "save/save_" << save_count << "_FILTER_" << filter_nr << "_Em";
	electric_field_minus->saveArray(fileName.str());

	fileName.str("");
	fileName << "save/save_" << save_count << "_FILTER_" << filter_nr << "_CavField";
	cavity_field->saveArray(fileName.str());
}

void Filter::file_load_variables(int load_count, int filter_nr)
{
	// Load all dynamic variables
	std::stringstream fileName;
	
	fileName.str("");
	fileName << "save/save_" << load_count << "_FILTER_" << filter_nr << "_Ep";
	electric_field_pluss->loadArray(fileName.str());

	fileName.str("");
	fileName << "save/save_" << load_count << "_FILTER_" << filter_nr << "_Ep_buff";
	electric_field_pluss_buffer->loadArray(fileName.str());
	
	fileName.str("");
	fileName << "save/save_" << load_count << "_FILTER_" << filter_nr << "_Em";
	electric_field_minus->loadArray(fileName.str());

	fileName.str("");
	fileName << "save/save_" << load_count << "_FILTER_" << filter_nr << "_CavField";
	cavity_field->loadArray(fileName.str());
}



