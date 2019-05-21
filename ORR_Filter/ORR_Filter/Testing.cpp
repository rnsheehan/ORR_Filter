#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the functions from the testing namespace

void testing::material_values()
{
	// How do you access the material classes in a generic way without having to have a switch statement of some kind? 
	// R. Sheehan 20 - 12 - 2018

	material *mat1;
	material *mat2;

	Si smpl1;

	SiO2 smpl2;

	InGaAs smpl3;

	double wavelength = 1.55;

	std::cout << "Wavelength: " << wavelength << " um\n"; 

	mat1 = &smpl1; // Si
	mat1->set_wavelength(wavelength);

	std::cout<<"Refractive index of Si: " << mat1->refractive_index() << "\n"; // this should return the RI of Si

	smpl2.set_wavelength(wavelength);
	std::cout<<"Refractive index of SiO2: " << smpl2.refractive_index() << "\n";// this should return the RI of SiO2

	mat2 = &smpl3; // InGaAs
	mat2->set_wavelength(wavelength);
	double infrac = 0.3; 

	std::cout <<"Refractive index of InGaAs: " << smpl3.refractive_index(infrac) << "\n";
	std::cout<<"Refractive index of InGaAs: " << mat2->refractive_index(infrac) << "\n"; // this should return the RI of InGaAs
}

void testing::disp_curve_wire()
{
	// Compute the dispersion curve for a wire waveguide object
	// R. Sheehan 7 - 3 - 2019

	bool pol = TM;

	int n_pts;
	double start, stop, W, H;

	sweep WL;
	wg_dims wire_dims;

	Air ri_air;
	SiN ri_si;
	SiO2 ri_sio2;

	n_pts = 200; start = 1.4; stop = 1.6;
	WL.set_vals(n_pts, start, stop);

	W = 2; H = 0.3;

	wire_dims.set_rect_wire(W, H);

	wire_dispersion disp_calc;

	std::string filename;

	filename = "Wire_WG_Dispersion_W_" + template_funcs::toString(W, 2) + "_H_" + template_funcs::toString(H, 2) + dottxt;

	disp_calc.compute_dispersion_data(pol, WL, wire_dims, &ri_si, &ri_sio2, &ri_air, filename);
}

void testing::disp_curve_rib()
{
	// Compute the dispersion curve for a rib waveguide object
	// R. Sheehan 8 - 3 - 2019

	bool pol = TM;

	int n_pts;
	double start, stop, W, E, T;

	sweep WL;
	wg_dims rib_dims;

	Air ri_air;
	Si ri_si;
	SiO2 ri_sio2;

	n_pts = 200; start = 1.4; stop = 1.6;
	WL.set_vals(n_pts, start, stop);


	W = 1; E = 0.5; T = 0.3;
	rib_dims.set_rib(W, E, T);

	rib_dispersion disp_calc;

	std::string filename = "Rib_WG_Dispersion.txt";

	disp_calc.compute_dispersion_data(pol, WL, rib_dims, &ri_si, &ri_sio2, &ri_air, filename);
}

void testing::orr_values()
{
	// test calculation for the ORR shown in Okamoto
	// R. Sheehan 19 - 2 - 2018

	// read the required data from a file
	std::string filename; 
		
	//filename = "Wire_WG_Dispersion_W_2.00_H_0.30.txt";
	filename = "Rib_WG_Dispersion.txt";

	int nrows, ncols; 
	std::vector<double> wavelengths, neff, ngroup; 
	std::vector<std::vector<double>> data; 

	vecut::read_into_matrix(filename, data, nrows, ncols); // read the data from the file

	wavelengths = vecut::get_col(data, 0); 
	neff = vecut::get_col(data, 1); 
	ngroup = vecut::get_col(data, 2);

	// Do the ORR calculation
	double min_feature = 0.5; 
	double notch_wl = 1.55; 
	double notch_neff, notch_ngrp; 
	double gamma = 2.3E-7; // dB / um
	double rho = gamma; // dB / um
	double kappa = 1.31E-5; // dB / um
	double Lcoup = 20000; // units of um
	double bendR = 3000; 

	// estimate the effective and group indices by interpolation from dispersion data
	double delta; 
	interpolation::polint(wavelengths, neff, notch_wl, notch_neff, delta); 
	interpolation::polint(wavelengths, ngroup, notch_wl, notch_ngrp, delta);

	ORR device; 

	device.compute_coefficients(min_feature, notch_wl, notch_neff, notch_ngrp, kappa, gamma, rho, Lcoup);

	device.report(); 

	double start = 1.4, end = 1.6, increment = 1.0/1000.0;
	filename = "Spectrum_Scan_Values.txt";
	device.spctrm_scan(start, end, increment, filename);
}

void testing::orr_design_wire()
{
	// Use this script to examine the properties of an ORR based on a certain waveguide
	// R. Sheehan 14 - 5 - 2019

	// parameters
	bool pol = TM;

	int n_pts;
	double start, stop, end, increment, W, H, separ, delta;
	double min_feature = 0.2; // minimum feature size possible with lithography
	double notch_wl = 1.5; // desired resonant wavelength, units of um
	double notch_neff, notch_ngrp; 
	double gamma = 2.3E-7; // insertion loss, dB / um
	double rho = gamma; // bend loss, dB / um
	double kappa = 1.31E-5; // dB / um
	double Lcoup = 20000; // units of um

	std::string filename;

	std::vector<double> separation_vals, coupling_vals, wavelengths, neff, ngroup;

	std::vector< std::vector< double > > disp_data;

	sweep WL;
	wg_dims wire_dims;

	Air ri_air;
	SiN ri_si;
	SiO2 ri_sio2;

	ri_vals wire_ri; 

	W = 1; H = 0.3;

	wire_dims.set_rect_wire(W, H);

	ri_air.set_wavelength(notch_wl); ri_si.set_wavelength(notch_wl); ri_sio2.set_wavelength(notch_wl); 

	wire_ri.set_rib_wire(ri_si.refractive_index(), ri_sio2.refractive_index(), ri_air.refractive_index(), notch_wl); 

	// 1. Compute the coupling coefficient as a function of waveguide separation
	// Save the coupling coefficient data and coupling length data
	Wire wire_coupling; 

	wire_coupling.set_params(pol, wire_dims, wire_ri); 

	wire_coupling.reduce_wg(); 

	wire_coupling.get_index(true); 

	separ = min_feature; 
	for (int i = 0; i < 50; i++) {
		double cc = wire_coupling.coupling_coefficient(separ); 
		
		//std::cout << separ << " , " << cc << "\n";

		separation_vals.push_back(separ); coupling_vals.push_back(cc); 

		separ += min_feature; 
	}

	// 2. Compute the dispersion curve for the waveguide so you can estimate the waveguide group index value
	// Having the dispersion data available means that you can examine if the structure resonates at other wavelengths
	n_pts = 200; start = 1.4; stop = 1.6;
	WL.set_vals(n_pts, start, stop);

	wire_dispersion disp_calc;

	filename = "Wire_WG_Dispersion_W_" + template_funcs::toString(W, 2) + "_H_" + template_funcs::toString(H, 2) + dottxt;

	disp_calc.compute_dispersion_data(pol, WL, wire_dims, &ri_si, &ri_sio2, &ri_air, filename);

	disp_data = disp_calc.get_computed_data(false); 

	wavelengths = vecut::get_col(disp_data, 0);
	neff = vecut::get_col(disp_data, 1);
	ngroup = vecut::get_col(disp_data, 2);	

	// 3. Compute the properties of the ORR using the data obtained from previous calculations
	
	// estimate the effective and group indices by interpolation from dispersion data
	interpolation::polint(wavelengths, neff, notch_wl, notch_neff, delta);

	interpolation::polint(wavelengths, ngroup, notch_wl, notch_ngrp, delta);

	// estimate the coupling coefficient based on waveguide separation
	separ = 0.3; // units of um
	double kappa_alt = kappa, l_alt; 
	interpolation::polint(separation_vals, coupling_vals, separ, kappa, delta);

	l_alt = 20000; // this assumes units of cm? 
	std::cout << "kappa estimation\n"; 
	std::cout << "kappa: " << kappa << ", kappa_alt: " << kappa_alt << "\n"; 
	std::cout << "kappa L: " << kappa*l_alt << ", kappa_alt L: " << kappa_alt*l_alt << "\n";
	std::cout << "cos(kappa L): " << cos(kappa*l_alt) << ", cos(kappa_alt L): " << cos(kappa_alt*l_alt) << "\n\n";

	// what's up with the weird scaling for kappa?
	//kappa = log10(kappa); 

	ORR device;

	//Lcoup = PI / (2.0*kappa); 

	Lcoup = 20e+3; // units of um!!

	device.compute_coefficients(min_feature, notch_wl, notch_neff, notch_ngrp, kappa, gamma, rho, Lcoup);

	device.report();

	start = 1.4; end = 1.6; increment = 1.0 / 1000.0;
	filename = "Spectrum_Scan_Values.txt";
	device.spctrm_scan(start, end, increment, filename);

	// Clear memory
	disp_data.clear(); wavelengths.clear(); neff.clear(); ngroup.clear(); separation_vals.clear(); coupling_vals.clear();
}