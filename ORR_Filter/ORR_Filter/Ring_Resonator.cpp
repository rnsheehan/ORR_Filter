#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the ORR class
// R. Sheehan 18 - 2 - 2019

ORR::ORR()
{
	// Default Constructor
	params_defined = false; 
	lambda = neff = ng = kappa = gamma = gamma_conj = rho = eff_OPL = grp_OPL = 0.0;
	R = Lcoup = L = X = X_conj = Y = Y_conj = XY = XY_conj = XY_pconj = X_conj_Y_conj = 0.0;
	Tmax = Tmin = FWHM = F = nu_FWHM = wl_FWHM = nu_FSR = wl_FSR = 0.0;
}

ORR::ORR(double wavelength, double eff_indx, double group_indx, double coup_coeff, double IL, double BL, double ring_rad, double ring_coup_len)
{
	// Constructor

	set_params(wavelength, eff_indx, group_indx, coup_coeff, IL, BL, ring_rad, ring_coup_len);
}

void ORR::set_params(double wavelength, double eff_indx, double group_indx, double coup_coeff, double IL, double BL, double ring_rad, double ring_coup_len)
{
	// assign parameters to the members of the ORR class
	// wavelength is the reference wavelength for the calculation
	// eff_indx is the waveguide effective index, not sure if it needs to be wavelength dependent yet
	// group indx is the waveguide group effective index, it definitely depends on wavelength
	// coup_coeff is the waveguide to ring coupling coefficient
	// IL is the waveguide insertion loss
	// BL is the ring bend loss
	// ring_rad is the ring bend radius
	// ring_coup_len is the length of the ring coupling section
	// R. Sheehan 18 - 2 - 2019
	
	try {
		bool c1 = eff_indx > 1.0 && eff_indx < 5.0 ? true : false; 
		bool c2 = coup_coeff > 0.0 && coup_coeff < 1.0 ? true : false;
		bool c3 = IL > 0.0 ? true : false; 
		bool c4 = BL > 0.0 ? true : false; 
		bool c5 = ring_rad > 0.0 ? true : false; 
		bool c6 = ring_coup_len > 0.0 ? true : false;
		bool c7 = abs(group_indx) > 0.0 ? true : false; 
		bool c8 = wavelength > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8 ? true : false; 

		if (c10) {
			// Define member values
			lambda = wavelength; 
			neff = eff_indx; 
			ng = group_indx; 
			kappa = coup_coeff; 
			gamma = IL; 
			gamma_conj = 1.0 - gamma; 

			rho = BL; 
			R = ring_rad; 
			Lcoup = ring_coup_len;
			L = 2.0*(Lcoup + PI * R); 

			eff_OPL = neff * L; // optical path length based on effective index 
			grp_OPL = ng * L; // optical path length based on group index

			nu_FSR = 300.0 / grp_OPL; // frequency FSR assuming c = 300 um.THz
			wl_FSR = template_funcs::DSQR(lambda) / grp_OPL; // wavelength FSR units of um

			X = sqrt(gamma_conj)*exp(-0.5*rho*L); 
			Y = cos(kappa*L); 
			XY = X * Y; 

			X_conj = 1.0 - template_funcs::DSQR(X); // 1 - x^{2}
			Y_conj = 1.0 - template_funcs::DSQR(Y); // 1 - y^{2}
			X_conj_Y_conj = X_conj * Y_conj; // (1 - x^{2}) (1 - y^{2})
			XY_conj = template_funcs::DSQR(1.0 - XY); // ( 1 - x y)^{2}
			XY_pconj = template_funcs::DSQR(1.0 + XY); // ( 1 + x y)^{2}

			Tmax = gamma_conj * ( template_funcs::DSQR(X + Y) / XY_pconj );
			Tmin = gamma_conj * ( template_funcs::DSQR(X - Y) / XY_conj );

			FWHM = 2.0 * sqrt( XY_conj / XY );  // 2.0*( 1 - x y) / ( x y )^{1/2}

			F = Two_PI / FWHM; // \pi * ( x y )^{1/2} / ( 1 - x y)

			nu_FWHM = nu_FSR / F; // frequency FWHM units of THz
			wl_FWHM = wl_FSR / F; // wavelength FWHM units of um

			params_defined = true; 
		}
		else {
			std::string reason; 
			reason = "Error: void ORR::set_params(double eff_indx, double kappa, double IL, double BL, double ring_rad, double ring_coup_len)\n"; 
			if (!c1) reason += "eff_indx: " + template_funcs::toString(eff_indx, 2) + " is not correct\n"; 
			if (!c2) reason += "coup_coeff: " + template_funcs::toString(coup_coeff, 2) + " is not correct\n";
			if (!c3) reason += "IL: " + template_funcs::toString(IL, 2) + " is not correct\n";
			if (!c4) reason += "BL: " + template_funcs::toString(BL, 2) + " is not correct\n";
			if (!c5) reason += "ring_rad: " + template_funcs::toString(ring_rad, 2) + " is not correct\n";
			if (!c6) reason += "ring_coup_len: " + template_funcs::toString(ring_coup_len, 2) + " is not correct\n";
			if (!c7) reason += "group_indx: " + template_funcs::toString(group_indx, 2) + " is not correct\n";
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

double ORR::T(double wavelength)
{
	// Compute the value of the intensity transmission spectrum at a wavelength

	try {
		if (params_defined) {

			double phase = (PI * eff_OPL) / wavelength; // actually \phi / 2 
			double denom = XY_conj + 4.0 * XY * template_funcs::DSQR(phase); 

			if (abs(denom) > 0.0) {
				double t1 = X_conj_Y_conj / denom; 
				return gamma_conj * (1.0 - t1); 
			}
			else {
				return 0.0; 
				std::string reason;
				reason = "Error: double ORR::T(double wavelength)\n";
				reason += "Attempt to divide by zero\n";
				throw std::runtime_error(reason);
			}
		}
		else {
			std::string reason; 
			reason = "Error: double ORR::T(double wavelength)\n"; 
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