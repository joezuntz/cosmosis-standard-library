#include "gsl/gsl_spline.h"

gsl_spline * shear_shear_kernel(double chi_max, gsl_spline * n_of_z, 
	gsl_spline * a_of_chi);


gsl_spline * cmb_wl_kappa_kernel(double chi_max, double chi_star, gsl_spline * a_of_chi);
