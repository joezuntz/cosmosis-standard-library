#include "utils.h"
#include <assert.h>
#include "stdio.h"
#include "string.h"
#include "cosmosis/datablock/c_datablock.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

double f_K(double K, double chi){
	if (K==0) return chi;
	if (K>0){
		double r = sqrt(K);
		return sin(r*chi)/r;
	}
	double r = sqrt(-K);
	return sinh(r*chi)/r;

}

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
	if (strlen(n_name)>0){
		status |= c_datablock_put_int(block, section, n_name, s->size);
	}
	if (strlen(x_name)>0){
		status |= c_datablock_put_double_array_1d(block, section, x_name, 
			s->x, s->size);
	}
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

static 
double identity_function(double x, double y, double z, void * a)
{
  return z;
}



Interpolator2D * 
load_interpolator_chi_function(c_datablock * block, gsl_spline * chi_of_z_spline, 
	const char * section,
	const char * k_name, const char * z_name, const char * P_name,
	interp2d_modifier_function function, void * args
	)
{
	int nk, nz, nP;
	double *k=NULL, *z=NULL;
	double **P=NULL;
	int status = 0;

	status = c_datablock_get_double_grid(block, section, 
		k_name, &nk, &k, 
		z_name, &nz, &z, 
		P_name, &P);

	if (status){
		fprintf(stderr, "Could not load interpolator for P(k).  Error %d\n",status);
		return NULL;
	}

	for (int j=0; j<nk; j++){
		for (int i=0; i<nz; i++){
			P[j][i] = function(k[j], z[i], P[j][i], args);
		}
	}


	// What we have now is P(k, z).
	// We can optionally convert to P(k, chi)
	// If so we loop, lookup, and replace
	if (chi_of_z_spline){
		for (int i=0; i<nz; i++){
			double zi = z[i];
			double chi_i = gsl_spline_eval(chi_of_z_spline, zi, NULL);
			z[i] = chi_i;
		}
	}

	if (status) return NULL;
	Interpolator2D * interp = init_interp_2d_akima_grid(k, z, P, nk, nz);
	deallocate_2d_double(&P, nk);
	free(k);
	free(z);
	return interp;
}


Interpolator2D * 
load_interpolator_chi(c_datablock * block, gsl_spline * chi_of_z_spline, 
	const char * section,
	const char * k_name, const char * z_name, const char * P_name)
{

	return load_interpolator_chi_function(block, chi_of_z_spline, section, k_name, z_name, P_name, identity_function, NULL);
}


Interpolator2D * 
load_interpolator_function(c_datablock * block, 
	const char * section,
	const char * k_name, const char * z_name, const char * P_name,
	interp2d_modifier_function function, void * args
	)
{
	return load_interpolator_chi_function(block, NULL, section, k_name, z_name, P_name, function, args);

}


Interpolator2D * 
load_interpolator(c_datablock * block, 
	const char * section,
	const char * k_name, const char * z_name, const char * P_name)
{
	return load_interpolator_chi_function(block, NULL, section, k_name, z_name, P_name, identity_function, NULL);
}

gsl_spline * load_spline(c_datablock * block, const char * section, 
	const char * x_name, const char * y_name)
{
	int status = 0 ;
	int nx, ny;
	double *x, *y;
	status |= c_datablock_get_double_array_1d(block, section, x_name, &x, &nx);
	if (status) return NULL;
	status |= c_datablock_get_double_array_1d(block, section, y_name, &y, &ny);
	if (status){
		free(x);
		return NULL;
	}

	gsl_spline * output = NULL;

	if (x[1]<x[0]) {reverse(x, nx); reverse(y, ny);}


	if (nx==ny) output=spline_from_arrays(nx, x, y);
	else {fprintf(stderr, "Non-matching array sizes %s=%d, %s=%d\n", x_name, nx, y_name, ny);}
	free(x);
	free(y);
	return output;
}

/*
gsl_spline2d *
load_gsl_interpolator_chi(c_datablock * block, gsl_spline * chi_of_z_spline, 
	const char * section,
	const char * k_name, const char * z_name, const char * P_name,
	)
{
	int nk, nz, nP;
	double *k=NULL, *z=NULL;
	double **P=NULL;
	int status = 0;

	status = c_datablock_get_double_grid(block, section, 
		k_name, &nk, &k, 
		z_name, &nz, &z, 
		P_name, &P);

	if (status){
		fprintf(stderr, "Could not load interpolator for P(k).  Error %d\n",status);
		return NULL;
	}

	// What we have now is P(k, z).
	// We can optionally convert to P(k, chi)
	// If so we loop, lookup, and replace
	if (chi_of_z_spline){
		for (int i=0; i<nz; i++){
			double zi = z[i];
			double chi_i = gsl_spline_eval(chi_of_z_spline, zi, NULL);
			z[i] = chi_i;
		}
	}

	if (status) return NULL;
        const gsl_interp2d_type *T = gsl_interp2d_bicubic;
	gsl_spline2d *spline = gsl_spline2d_alloc(T, nz, nk);
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
	gsl_interp_accel *yacc = gsl_interp_accel_alloc();
	gsl_spline_init

	Interpolator2D * interp = init_interp_2d_akima_grid(k, z, P, nk, nz);
	deallocate_2d_double(&P, nk);
	free(k);
	free(z);
	return interp;
}
*/
