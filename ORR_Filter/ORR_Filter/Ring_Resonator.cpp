#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the ORR class
// R. Sheehan 18 - 2 - 2019

ORR::ORR()
{
	// Default Constructor
	params_defined = false;
	res_order = 0; 
	min_length = notch_lambda = notch_neff = notch_ngrp = kappa = gamma = gamma_conj = rho = eff_OPL = grp_OPL = 0.0;
	R = Lcoup = L = X = X_conj = Y = Y_conj = XY = XY_conj = XY_pconj = X_conj_Y_conj = 0.0;
	Tmax = Tmin = phi_FWHM = F = nu_FWHM = wl_FWHM = nu_FSR = wl_FSR = 0.0;
}

ORR::ORR(double &min_size, double &notch_wavelength, double &eff_indx, double &grp_index, double &coup_coeff, double &IL, double &BL, double &ring_coup_len)
{
	// Constructor

	compute_coefficients(min_size, notch_wavelength, eff_indx, grp_index, coup_coeff, IL, BL, ring_coup_len);
}

ORR::~ORR()
{
	params_defined = false;
}

void ORR::compute_coefficients(double &min_size, double &notch_wavelength, double &eff_indx, double &grp_index, double &coup_coeff, double &IL, double &BL, double &ring_coup_len)
{
	// assign parameters to the members of the ORR class, these should be updated anytime the parameters change
	// notch_wavelength is the reference wavelength for the calculation
	// coup_coeff is the waveguide to ring coupling coefficient
	// IL is the waveguide insertion loss
	// BL is the ring bend loss
	// ring_coup_len is the length of the ring coupling section (I have some queries about this parameter)
	// R. Sheehan 13 - 5 - 2019

	try {			
		bool c1 = min_size > 0.1 ? true : false; 
		bool c5 = notch_wavelength > 0.0 ? true : false;
		bool c7 = eff_indx > 1.0 ? true : false; 
		bool c9 = grp_index > 1.0 ? true : false;
		bool c2 = coup_coeff > 0.0 && coup_coeff < 1.0 ? true : false;
		bool c3 = IL > 0.0 ? true : false;
		bool c4 = BL > 0.0 ? true : false;		
		bool c6 = ring_coup_len > 0.0 ? true : false;		
		bool c8 = notch_wavelength > 0.0 ? true : false;
		bool c10 = c1 && c5 && c7 && c9 && c2 && c3 && c4 && c6 && c8 ? true : false;

		if (c10) {
			// Define member values
			min_length = PI * min_size;
			
			notch_lambda = notch_wavelength; // units of um
			notch_neff = eff_indx; 
			notch_ngrp = grp_index; 
			
			// this is the minimum resonator order required
			res_order = 1 + static_cast<int>(std::ceil( ( notch_neff / notch_lambda ) * min_length )  );

			kappa = coup_coeff; // units of dB / um
			gamma = IL; // units of dB / um
			gamma_conj = 1.0 - gamma;

			rho = BL; // bend loss, units of dB / um
			Lcoup = ring_coup_len; // coupling length (?) units of um
			L = 150 * res_order * (notch_lambda / notch_neff); // need to optimise how the value of m is chosen
			R = L / Two_PI; // ring radius value 

			eff_OPL = notch_neff * L; // optical path length based on effective index 
			grp_OPL = notch_ngrp * L; // optical path length based on group index

			nu_FSR = (300.0 / grp_OPL); // frequency FSR assuming c = 300 um.THz
			wl_FSR = template_funcs::DSQR(notch_lambda) / grp_OPL; // wavelength FSR units of um

			X = sqrt(gamma_conj)*exp(-0.5*rho*L);
			Y = cos(kappa*L);
			XY = X * Y;

			X_conj = 1.0 - template_funcs::DSQR(X); // 1 - x^{2}
			Y_conj = 1.0 - template_funcs::DSQR(Y); // 1 - y^{2}
			X_conj_Y_conj = X_conj * Y_conj; // (1 - x^{2}) (1 - y^{2})
			XY_conj = template_funcs::DSQR(1.0 - XY); // ( 1 - x y)^{2}
			XY_pconj = template_funcs::DSQR(1.0 + XY); // ( 1 + x y)^{2}

			Tmax = gamma_conj * (template_funcs::DSQR(X + Y) / XY_pconj);
			Tmin = gamma_conj * (template_funcs::DSQR(X - Y) / XY_conj);

			phi_FWHM = 2.0 * ( sqrt(XY_conj) / sqrt(XY) );  // 2.0*( 1 - x y) / ( x y )^{1/2}

			F = PI / phi_FWHM; // \pi * ( x y )^{1/2} / ( 1 - x y), had to exclude multiplicative factor of 2 to make it match example in Okamoto

			nu_FWHM = nu_FSR / F; // frequency FWHM units of THz
			wl_FWHM = wl_FSR / F; // wavelength FWHM units of um

			params_defined = true;
		}
		else {
			std::string reason;
			reason = "Error: void ORR::compute_coefficients(double &notch_wavelength, double &coup_coeff, double &IL, double &BL, double &ring_coup_len)\n";
			if (!c2) reason += "coup_coeff: " + template_funcs::toString(coup_coeff, 2) + " is not correct\n";
			if (!c3) reason += "IL: " + template_funcs::toString(IL, 2) + " is not correct\n";
			if (!c4) reason += "BL: " + template_funcs::toString(BL, 2) + " is not correct\n";
			if (!c6) reason += "ring_coup_len: " + template_funcs::toString(ring_coup_len, 2) + " is not correct\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void ORR::report()
{
	// Report on the computed parameters for the ORR device
	// R. Sheehan 19 - 2 - 2018

	try {
		if (params_defined) {
			std::cout << "Optical Ring Resonator Parameters\n\n";
			std::cout << "Reference wavelength: " << notch_lambda << " um\n";
			std::cout << "Waveguide effective index: " << notch_neff << "\n";
			std::cout << "Waveguide group index: " << notch_ngrp << "\n";
			std::cout << "Minimum Size: " << min_length << " um\n"; 
			std::cout << "Resonator Order: " << res_order << "\n\n";
			
			std::cout << "Insertion loss: "<<gamma<<" dB / um\n";
			std::cout << "Bend loss: " << rho << " dB / um\n";
			std::cout << "Coupling coefficient: " << kappa << " dB/um\n\n"; 

			std::cout << "Ring radius: " << R << " um\n";
			std::cout << "Coupler length: " << L << " um\n";
			std::cout << "Total ring length: " << L << " um\n"; 
			std::cout << "Effective optical path length: " << eff_OPL / 1000.0 << " mm\n";
			std::cout << "Group optical path length: " << grp_OPL / 1000.0 << " mm\n";
			std::cout << "kappa L: " << kappa * L << "\n\n";

			std::cout << "Attenuation term X: " << X << "\n"; 
			std::cout << "Coupling term Y: " << Y << "\n\n";

			std::cout << "Tmax: " << Tmax << "\n"; 
			std::cout << "Tmin: " << Tmin << "\n"; 
			std::cout << "Finesse: " << F << "\n"; 
			std::cout << "Frequency FSR: " << 1000.0*nu_FSR << " GHz\n"; 
			std::cout << "Frequency FWHM: " << 1.0e+6*nu_FWHM << " MHz\n"; 
			std::cout << "Wavelength FSR: " << 1000.0*wl_FSR << " nm\n"; 
			std::cout << "Wavelength FWHM: " << 1.0e+6*wl_FWHM << " pm\n\n"; 
		}
		else {
			std::string reason;
			reason = "Error: void ORR::report()\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double ORR::through_spctrm(double lambda)
{
	 /*Compute the value of the intensity transmission spectrum at a wavelength
	 In the transmission spectrum phi is the only term that is not constant
	 Code must be made to depend on wavelength dependent effective index data*/

	try {
		if (params_defined) {

			double phase = (PI * eff_OPL) / lambda; // actually \phi / 2 
			double denom = XY_conj + 4.0 * XY * template_funcs::DSQR( sin(phase) ); 

			if (abs(denom) > 0.0) {
				double t1 = X_conj_Y_conj / denom; 
				return gamma_conj * (1.0 - t1); 
			}
			else {
				return 0.0; 
				std::string reason;
				reason = "Error: double ORR::through_spctrm(int i)\n";
				reason += "Attempt to divide by zero\n";
				throw std::runtime_error(reason);
			}
		}
		else {
			std::string reason; 
			reason = "Error: double ORR::through_scptrm(int i)\n"; 
			reason += "Parameters not defined\n";
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what(); 
	}
}

double ORR::drop_scptrm(double lambda)
{
	// Compute the value of the intensity reflection spectrum at a wavelength

	try {
		if (params_defined) {
			return 1.0 - through_spctrm(lambda);
		}
		else {
			return 0.0; 
			std::string reason;
			reason = "Error: double ORR::drop_scptrm(double lambda)\n";
			reason += "Parameters not defined\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void ORR::spctrm_scan(double start, double end, double increment, std::string filename)
{
	// scan the wavelength spectrum to compute the through spectrum
	// R. Sheehan 19 - 2 - 2019

	try {

		bool c1 = start < end ? true : false;
		bool c2 = end > 0.0 ? true : false;
		bool c3 = increment > 0.0 ? true : false;
		bool c4 = abs(end - start) > increment ? true : false;
		bool c10 = params_defined && c1 && c2 && c3 && c4;
		
		if (c10) {

			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

			if ( write.is_open() ) {
				
				// loop over wavelength values
				double x; 

				int n_wl = static_cast<int>(1 + ((end - start) / increment));

				x = start; 
				for (int i = 0; i < n_wl; i++) {

					//double Tval = through_spctrm(i); // compute transmission value
					double Tval = through_spctrm(x); // compute transmission value

					// write the data to a file
					write << std::setprecision(10) << x << " , " << Tval << " , " << 1.0 - Tval << "\n";

					x += increment;
				}			
			}
			else {
				std::string reason;
				reason = "Error: void ORR::spctrm_scan(double start, double end, double increment, std::string filename)\n";
				reason += "Could not open file: " + filename + "\n";
				throw std::runtime_error(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void ORR::spctrm_scan(double start, double end, double increment, std::string filename)\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}