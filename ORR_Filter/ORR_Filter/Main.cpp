#ifndef ATTACH_H
#include "Attach.h"
#endif

// Projects implements the code needed to compute the transmission spectrum of an optical ring resonator
// Implementation follows theory given in K. Okamoto, "Fundamentals of Optical Waveguides", sect. 4.5.2
// R. Sheehan 18 - 2 - 2019

int main()
{
	//testing::material_values(); 

	//testing::disp_curve_wire(); 

	//testing::disp_curve_rib(); 

	testing::orr_values(); 	

	std::cout << "Press enter to close\n"; 
	std::cin.get(); 

	return 0; 
}