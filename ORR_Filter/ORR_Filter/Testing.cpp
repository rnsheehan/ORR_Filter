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