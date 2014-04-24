#include "utils.h"
#include "assert.h"

gsl_spline * spline_from_arrays(int n, double * x, double *y)
{
	gsl_spline * output = gsl_spline_alloc(gsl_interp_akima, n);
	assert (output!=NULL);
	gsl_spline_init(output,x,y,n);
	return output;
}


DATABLOCK_STATUS save_spline(c_datablock * block, const char * section, 
	const char * n_name, const char * x_name, const char * y_name,
	gsl_spline * s)
{

	DATABLOCK_STATUS status = 0;
	status |= c_datablock_put_int(block, section, n_name, s->size);
	status |= c_datablock_put_double_array_1d(block, section, x_name, 
		s->x, s->size);
	status |= c_datablock_put_double_array_1d(block, section, y_name, 
		s->y, s->size);
	return status;
}

void reverse(double * x, int n)
{
	for (int i=0; i<n/2; i++){
		double tmp = x[i];
		x[i] = x[n-1-i];
		x[n-1-i] = tmp;
	}
}