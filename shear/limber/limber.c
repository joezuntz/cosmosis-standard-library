#include "stdio.h"
#include "gsl/gsl_integration.h"
#include "limber.h"

// This is a workspace size for the gsl integrator
#define LIMBER_FIXED_TABLE_SIZE 1000

// data that is passed into the integrator
// This is everything we need to compute the
// integrand
typedef struct IntegrandData{
	double chimin;
	double chimax;
	double ell;
	gsl_spline * WX;
	gsl_spline * WY;
	Interpolator2D * P;
	gsl_interp_accel * accelerator;
} IntegrandData;

// These two convenience functions
// peer into the internals of the gsl_spline.
// This is probably a bit naughty, since they could
// in theory change the internals.
static double inline limber_gsl_spline_min_x(gsl_spline * s)
{
	return s->x[0];
}
static double inline limber_gsl_spline_max_x(gsl_spline * s)
{
	return s->x[s->size-1];
}

// the integrand W_X(chi) * W_Y(chi) * P(ell/chi, chi) / chi^2
static double integrand(double chi, void * data_void)
{
	IntegrandData * data = (IntegrandData*) data_void;
	// Return 0 if outside range, for convenience.
	// Important to ensure that ranges are wide enough.
	if(chi < data->chimin || chi > data->chimax) return 0.0;

	// Get W^X(chi) and W^Y(chi)
	double wx = gsl_spline_eval(data->WX,chi,data->accelerator);
	double wy;
	// A short-cut - if the two splines are the same we do not need to 
	// do the interpolation twice
	if (data->WX==data->WY) wy = wx;
	else wy = gsl_spline_eval(data->WY,chi,data->accelerator);
	if (wx==0 || wy==0) return 0.0;

	// Get P(k,z) using k=ell/chi.
	// The interp_2d interpolator returns 0 if either 
	// parameter is outside its range
	double k = data->ell / chi;
	double p = interp_2d(k, chi, data->P);

	// Integrand result.
	double result = wx * wy * p / chi / chi;
	return result;

}

// The only function in this little library callable from the outside
// world.  The limber_config structure is defined in limber.h but is fairly
// obvious.  The splines and the interpolator need to be functions of 
// chi NOT z.
gsl_spline * limber_integral(limber_config * config, gsl_spline * WX, 
	                 gsl_spline * WY, Interpolator2D * P)
{
	// Get the appropriate ranges over which to integrate
	// It is assumed that (at least one of) the kernel
	// splines should go to zero in some kind of reasonable
	// place, so we just use the range they specify
	IntegrandData data;
	double chimin_x = limber_gsl_spline_min_x(WX);
	double chimin_y = limber_gsl_spline_min_x(WY);
	double chimax_x = limber_gsl_spline_max_x(WX);
	double chimax_y = limber_gsl_spline_max_x(WY);
	// Take the smallest range since we want both the
	// splines to be valid there.
	// This range as well as all the data needed to compute
	// the integrand is put into a struct to be passed
	// through the integrator to the function above.
	data.chimin = chimin_x>chimin_y ? chimin_x : chimin_y;
	data.chimax = chimax_x<chimax_y ? chimax_x : chimax_y;
	data.WX = WX;
	data.WY = WY;
	data.P = P;
	data.accelerator = gsl_interp_accel_alloc();

	// Set up the workspace and inputs to the integrator.
	// Not entirely sure what the table is.
	gsl_function F;
	F.function = integrand;
	F.params = &data;
	gsl_integration_glfixed_table *table = 
	    gsl_integration_glfixed_table_alloc((size_t) LIMBER_FIXED_TABLE_SIZE);

	// results of the integration go into these arrays.
	double c_ell_vector[config->n_ell];
	double ell_vector[config->n_ell];

	// loop through ell values according to the input configuration
	for (int i_ell = 0; i_ell<config->n_ell; i_ell++){
		double ell = config->ell[i_ell];
		data.ell=ell;

		// Perform the main integration.
		// This particular function is used because that's what Matt Becker 
		// found to work best.
		double c_ell = gsl_integration_glfixed(&F,data.chimin,data.chimax,table);

		//Include the prefactor scaling
		c_ell *= config->prefactor;
		printf("Limber: %d  %le  %le  %le  %le  %le\n", i_ell, ell, c_ell, config->prefactor, data.chimin, data.chimax);

		// It is often useful to interpolate into the logs of the functions
		// This is optional in the config
		if (config->xlog) ell = log(ell);
		if (config->ylog) c_ell = log(c_ell);

		// Record the results into arrays
		c_ell_vector[i_ell] = c_ell;
		ell_vector[i_ell] = ell;
	}

	// Create a spline of the arrays as the output
	gsl_spline * output = gsl_spline_alloc(gsl_interp_akima, (size_t) config->n_ell);
	gsl_spline_init(output, ell_vector, c_ell_vector, (size_t) config->n_ell);

	// Tidy up
	gsl_interp_accel_free(data.accelerator);
	gsl_integration_glfixed_table_free(table);	

	// And that's it
	return output;
}