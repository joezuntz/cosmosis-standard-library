#ifndef _H_LIMBER
#define _H_LIMBER

#include <gsl/gsl_spline.h>
#include <stdbool.h>
#include "interp2d.h"

// These are the options you can set for
// the Limber integrator.
typedef struct limber_config{
	bool xlog;  // The output spline will be in terms of log(ell) not ell
	bool ylog;  // The output spline will return log(C_ell) not C_ell
	int n_ell;  // Number of ell values you want in the spline
	double * ell;  // The chosen ell values you want
} limber_config;


// Do a flat universe Limber approximation integral, of the form:
// C^{XY}(\ell) 
// = \int_0^{\chi_{\mathrm{max}}} 
//       \frac{W^X(\chi) W^Y(\chi)}{\chi^2} P(\ell/\chi, \chi)

// The two splines and the matter power passed into the function
// be FUNCTIONS OF CHI NOT Z

// TODO - Fortran version

gsl_spline * limber_integral(limber_config * config, 
	gsl_spline * WX, gsl_spline * WY, Interpolator2D * P);


#endif