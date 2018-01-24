#include "stdio.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "non_limber.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_sf_bessel.h>

// This is a workspace size for the gsl integrator
#define LIMBER_FIXED_TABLE_SIZE 4096

// data that is passed into the integrator
// This is everything we need to compute the
// integrand
struct ClIntegrandData {
	double ell;
	double chimin;
	double chimax;
	gsl_spline * WX;
	gsl_spline * WY;
    gsl_spline * D_chi;
	gsl_spline2d * P;
	gsl_interp_accel * accelerator_wx;
	gsl_interp_accel * accelerator_wy;
	gsl_interp_accel * accelerator_px;
	gsl_interp_accel * accelerator_py;
	gsl_interp_accel * accelerator_growth;
};

double cl_integrand(double *x, size_t dim, void *data_void)
{
	struct ClIntegrandData * data = (struct ClIntegrandData*) data_void;
	
	double chi1 = x[0];
	double log_delta = x[1];
	double log_nu = x[2];
	double nu = exp(log_nu);
	double delta = exp(log_delta);
	double chi2_p = chi1 * ( 1 + delta);
	double chi2_m = chi1 * (1 - delta);
	double logk = log_nu / chi1;

	if(chi2_m < data->chimin || chi2_p > data->chimax) return 0.0;

	double FX = gsl_spline_eval(data->WX, chi1, data->accelerator_wx);
	double FY_p = gsl_spline_eval(data->WY, chi2_p, data->accelerator_wy);
	double FY_m = gsl_spline_eval(data->WY, chi2_m, data->accelerator_wy);


	double growth_1 = gsl_spline_eval(data->D_chi, chi1, data->accelerator_growth);
	double growth_2_p = gsl_spline_eval(data->D_chi, chi2_p, data->accelerator_growth);
	double growth_2_m = gsl_spline_eval(data->D_chi, chi2_m, data->accelerator_growth);

	// Get P(k,chi)
	// Deactivate error handling 
	gsl_error_handler_t * old_handler = gsl_set_error_handler_off();

	double p;
	int status = gsl_spline2d_eval_e(data->P, chi1, logk, data->accelerator_px, data->accelerator_py, &p);

	// Restore the old error handler
	gsl_set_error_handler(old_handler); 
	if (status) 
	{
		//printf(stderr,'spline failed, for now, assume this is because k was out of range and return 0...should probably be more careful here though.');
		return 0.0;
	}
	double j1 = gsl_sf_bessel_jl(data->ell, nu);
	double j2_p = gsl_sf_bessel_jl(data->ell, nu*(1+delta));
	double j2_m = gsl_sf_bessel_jl(data->ell, nu*(1-delta));
	double prefac = 2 * delta * 2 * nu * nu * nu * p / growth_1 / M_PI / chi1 / chi1 ;
	return prefac * ( FY_p * growth_2_p * j2_p + FY_m * growth_2_m * j2_m );
}

static double inline gsl_spline_min_x(gsl_spline * s)
{
	return s->x[0];
}
static double inline gsl_spline_max_x(gsl_spline * s)
{
	return s->x[s->size-1];
}

int get_chi_mean_sigma(double chimin, double chimax, gsl_spline * WX, gsl_spline * WY,
				gsl_interp_accel * accelerator_wx,	gsl_interp_accel * accelerator_wy,
				double *chi_mean, double *chi_sigma)
{
    int nchi = 1000;
    double dchi = (chimax-chimin)/nchi;
    double chi = chimin;
    double wxwy_sum = 0.;
    double wxwy2_sum = 0.;
    double wx,wy;
    for (int i=0; i<nchi; i++){
    	wx = gsl_spline_eval(WX, chi, accelerator_wx);
    	wy = gsl_spline_eval(WY, chi, accelerator_wy);	
    	wxwy_sum += wx*wy*chi;
    	wxwy2_sum += wx*wy*chi*chi;
    	chi += dchi;
    }
    *chi_mean = wxwy_sum*dchi;
    *chi_sigma = sqrt(wxwy2_sum*dchi - *chi_mean * *chi_mean);
    return 0;
}

int cl_integral(cl_config * config, gsl_spline * WX, 
			    gsl_spline * WY, gsl_spline2d * P, gsl_spline * D_chi,
			    double * cl_out, double * cl_err_out)
{

    config->status = LIMBER_STATUS_ERROR;
    int any_parameter_error=0;
    if (WX==NULL){
        fprintf(stderr, "NULL WX parameter in limber_integral\n");
        any_parameter_error = 1;
    }
    if (WY==NULL){
        fprintf(stderr, "NULL WY parameter in limber_integral\n");
        any_parameter_error = 1;
    }
    if (P==NULL){
        fprintf(stderr, "NULL P parameter in limber_integral\n");
        any_parameter_error = 1;
    }
    if (D_chi==NULL){
        fprintf(stderr, "NULL D_chi parameter in limber_integral\n");
        any_parameter_error = 1;
    }
    if (config->n_ell<0){
        fprintf(stderr, "Negative n_ell parameter in limber_integral\n");
        any_parameter_error = 1;
    }
    if (any_parameter_error){
        return 1;
    }

	config->status = LIMBER_STATUS_OK;

    static int n_ell_zero_warning = 0;
    if (config->n_ell==0){
        if (n_ell_zero_warning==0){
            fprintf(stderr, "Warning: n_ell=0 in Limber. Will not be treated as an error. Warning once only per process.\n");
        }
        n_ell_zero_warning = 1;
        return 1;
    }

	// Get the appropriate ranges over which to integrate
	// It is assumed that (at least one of) the kernel
	// splines should go to zero in some kind of reasonable
	// place, so we just use the range they specify
	struct ClIntegrandData data;
	double chimin_x = gsl_spline_min_x(WX);
	double chimin_y = gsl_spline_min_x(WY);
	double chimax_x = gsl_spline_max_x(WX);
	double chimax_y = gsl_spline_max_x(WY);
	double c_ell, c_ell_error;

	// Workspaces
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    size_t calls = 5000000;

	double reltol = config->relative_tolerance;
	double abstol = config->absolute_tolerance;

	// Take the smallest range since we want both the
	// splines to be valid there.
	// This range as well as all the data needed to compute
	// the integrand is put into a struct to be passed
	// through the integrator to the function above.
	double chimin = chimin_x>chimin_y ? chimin_x : chimin_y;
    double chimax = chimax_x<chimax_y ? chimax_x : chimax_y;
	data.chimin = chimin_y;
	data.chimax = chimax_y;
	data.WX = WX;
	data.WY = WY;
	data.P = P;
	data.accelerator_wx = gsl_interp_accel_alloc();
	data.accelerator_wy = gsl_interp_accel_alloc();
	data.accelerator_px = gsl_interp_accel_alloc();
	data.accelerator_py = gsl_interp_accel_alloc();
	data.accelerator_growth = gsl_interp_accel_alloc();
	data.D_chi = D_chi;

	double chi_mean, chi_sigma;
	get_chi_mean_sigma(chimin, chimax, WX, WY,
				data.accelerator_wx, data.accelerator_wy,
				&chi_mean, &chi_sigma);

	// save ell values for return
	double ell_vector[config->n_ell];
	double nu_0;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (3);

	// loop through ell values according to the input configuration
	for (int i_ell = 0; i_ell<config->n_ell; i_ell++){
		double ell = config->ell[i_ell];
		double log_nu_0 = log(ell+0.5);
		double log_nu_min = log_nu_0 - config->log_nu_range;
		double log_nu_max = log_nu_0 + config->log_nu_range;
		printf("delta_range_factor %f\n", config->delta_range_factor);
		double delta_min = 1. / ell / config->delta_range_factor;
		double delta_max = 1. * config->delta_range_factor / ell;

		double log_delta_min = log(delta_min);
		double log_delta_max = log(delta_max);
		printf("log_delta_min,max %f %f\n", delta_min, delta_max);
		double xl[3] = { chimin_x, log_delta_min, log_nu_min };
		double xu[3] = { chimax_x, log_delta_max, log_nu_max };
		data.ell=ell;
		gsl_monte_function G = {&cl_integrand, 3, &data};

        gsl_monte_vegas_integrate (&G, xl, xu, 3, 10000, r, s,
                               &c_ell, &c_ell_error);
        gsl_monte_vegas_integrate (&G, xl, xu, 3, calls/5, r, s,
                               &c_ell, &c_ell_error);
        gsl_monte_vegas_init (s);

		//Include the prefactor scaling
		c_ell *= config->prefactor;
		// Record the results into arrays
		cl_out[i_ell] = c_ell;
		cl_err_out[i_ell] = c_ell_error;
		ell_vector[i_ell] = ell;
	}
	gsl_monte_vegas_free (s);


	// Tidy up
	gsl_interp_accel_free(data.accelerator_wx);
	gsl_interp_accel_free(data.accelerator_wy);
	gsl_interp_accel_free(data.accelerator_px);
	gsl_interp_accel_free(data.accelerator_py);
	gsl_interp_accel_free(data.accelerator_growth);

	// These two are not deallocated because they are static and only initialized once.
	// gsl_integration_glfixed_table_free(table);	
	// gsl_integration_workspace_free(W);

	// And that's it
	return 0;
}