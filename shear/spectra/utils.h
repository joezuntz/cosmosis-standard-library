#include "gsl/gsl_spline.h"
#include "cosmosis/datablock/c_datablock.h"


gsl_spline * spline_from_arrays(int n, double * x, double *y);

DATABLOCK_STATUS save_spline(c_datablock * block, const char * section, 
	const char * n_name, const char * x_name, const char * y_name,
	gsl_spline * s);

void reverse(double * x, int n);
