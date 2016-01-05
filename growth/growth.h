#include <stdlib.h>
#include <math.h>
#include <string.h>




				
/* Code Flags and Parameters */
	
#define choice_of_growth_function 1
	/* 1 = solve diff. eq. assuming GR. Eq. 7.77 in Dodelson 
	   2 = modified gravity parameterization */
				
/* Cosmological parameters */
//const char * cospar = COSMOLOGICAL_PARAMETERS_SECTION;


	/*
	The following two quantities characterise growth discrepancy w/ GR, so that
	actual growth function = (GR prediction)*exp(growth_G0)*(1.0+z)^growth_gamma.
	
	Right now it's implemented only in the camb power spectrum.
	*/

#define GROWTH_G0 0.0				
#define GROWTH_GAMMA 0.0

/* Support quantities */

#define TOL (1.0e-6)
#define PI 3.14159265358979324 /* 17-digits */
#define SQRT2 1.414213562


