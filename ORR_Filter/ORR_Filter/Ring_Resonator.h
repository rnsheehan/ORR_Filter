#ifndef RING_RESONATOR_H
#define RING_RESONATOR_H

// Code for defining the optical ring resonator
// Implementation follows that given in K. Okamoto, "Fundamentals of Optical Waveguides", sect. 4.5.2
// Code implements a single notch-filter type resonator with ring in racetrack configuration
// R. Sheehan 13 - 5 - 2019

// Make sure to add code to estimate Q-factor of the resonator
// R. Sheehan 13 - 5 - 2019

class ORR {
public:
	ORR();
	ORR(std::vector<double> &wl_vals, std::vector<double> &neff_vals, std::vector<double> &ngrp_vals);
	~ORR(); 

	void set_params(std::vector<double> &wl_vals, std::vector<double> &neff_vals, std::vector<double> &ngrp_vals);
	
	void compute_coefficients(double &notch_wavelength, double &coup_coeff, double &IL, double &BL, double &ring_coup_len);

	void report();

	void spctrm_scan(std::string filename);

private:
	double through_spctrm(int i);

	double drop_scptrm(int i);

private:
	bool inputs_defined; // switch that decides if inputs have been input correctly
	bool params_defined; // switch that decides if parameters have been computed correctly

	int n_wl; // no. wavelengths in input spectrum

	double notch_lambda; // resonance wavelength
	double notch_neff; // effective index at resonance wavelength
	double notch_ngrp; // group index at resonance wavelength
	double kappa; // coupling coefficient between ring waveguide and input waveguide, this will include waveguide ring separation, units of um^{-1}
	double gamma; // intensity insertion loss coefficient, units of um^{-1}
	double gamma_conj; // 1-gamma
	double rho; // intensity attenuation coefficient in the ring, bend loss, units of um^{-1}

	double R; // ring radius, units of um
	double Lcoup; // length of coupling section, units of um
	double L; // ring round trip length, this should include length of directional coupler section, units of um
	double eff_OPL; // effective optical path length
	double grp_OPL; // group optical path length
	
	double X; // attenuation term
	double Y; // coupling strength term
	double XY; // parameter used in calculation of T

	double X_conj; // parameter used in calculation of T
	double Y_conj; // parameter used in calculation of T
	double X_conj_Y_conj; // parameter used in calculation of T
	double XY_conj; // parameter used in calculation of T
	double XY_pconj; // parameter used in calculation of T

	double Tmax; // max value of transmission spectrum
	double Tmin; // min value of transmission spectrum
	double F; // ORR finesse
	double FWHM; // resonance peak FWHM in phase space
	double nu_FWHM; // resonance peak FWHM in frequency space
	double wl_FWHM; // resonance peak FWHM in wavelength space
	double nu_FSR; // frequency free spectral range
	double wl_FSR; // wavelength free spectral range

	std::vector<double> wavelengths; // must be in units of um
	std::vector<double> effective_index; 
	std::vector<double> group_index; 
};

#endif
