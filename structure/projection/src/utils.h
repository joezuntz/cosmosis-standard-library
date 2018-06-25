#ifndef _H_LIMBER_UTILS
#define _H_LIMBER_UTILS

#include "gsl/gsl_spline.h"
#include "cosmosis/datablock/c_datablock.h"
#include "interp2d.h"



gsl_spline * spline_from_arrays(int n, double * x, double *y);

DATABLOCK_STATUS save_spline(c_datablock * block, const char * section, 
	const char * n_name, const char * x_name, const char * y_name,
	gsl_spline * s);

gsl_spline * load_spline(c_datablock * block, const char * section, 
	const char * x_name, const char * y_name);

double f_K(double K, double chi);


void reverse(double * x, int n);


Interpolator2D * 
load_interpolator(c_datablock * block,
	const char * section,
	const char * k_name, const char * z_name, const char * P_name);


Interpolator2D * 
load_interpolator_chi(c_datablock * block, gsl_spline * chi_of_z,
	const char * section,
	const char * k_name, const char * z_name, const char * P_name);

Interpolator2D * 
load_interpolator_function(c_datablock * block, 
	const char * section,
	const char * k_name, const char * z_name, const char * P_name,
	interp2d_modifier_function function, void * args);

Interpolator2D * 
load_interpolator_chi_function(c_datablock * block, 
	gsl_spline * chi_of_z_spline,
	const char * section, 
	const char * k_name, const char * z_name, const char * P_name,
	interp2d_modifier_function function, void * args);

#endif