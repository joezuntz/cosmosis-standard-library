#include "shear_shear.h"
#include "gsl/gsl_integration.h"
#include "utils.h"
#include "gsl/gsl_errno.h"

#define NGLT 2048


typedef struct w_integrand_data{
	gsl_spline * n_of_z;
	gsl_spline * a_of_chi;
	gsl_interp_accel * acc_chi;
	gsl_interp_accel * acc_z;
	double chi;
	double zmax;
	double chi_max;
} w_integrand_data;


static double n_of_chi(double chi, void *p)
{
	w_integrand_data * data = (w_integrand_data*) p;
	if (chi>data->chi_max) return 0.0;
	// Get a(chi), da/dchi(chi) and z(chi)
	double a = gsl_spline_eval(data->a_of_chi, chi, data->acc_chi);
	double da_dchi = gsl_spline_eval_deriv(data->a_of_chi, chi, data->acc_chi);
	double z = 1.0/a - 1.0;
	if (z>=data->zmax) {
		// short circuit this in future
		data->chi_max = chi;
		return 0.0;
	}
	// Get n(chi) from n(z)
	double nz = gsl_spline_eval(data->n_of_z, z, data->acc_z);
	return -(1+z)*(1+z)*nz*da_dchi;
}

static double w_of_chi_integrand(double chis, void *p)
{
	w_integrand_data * data = (w_integrand_data*) p;
	double nchi = n_of_chi(chis, data);
	return nchi*(chis - data->chi)/chis;
}

// this function does not include the factor
// 1.5 omega_m H0^2 in it.  That should be included in the
// prefactor later
gsl_spline * shear_shear_kernel(double chi_max, gsl_spline * n_of_z, 
	gsl_spline * a_of_chi)
{


	// Integrand data
	w_integrand_data data;
	data.n_of_z = n_of_z;
	data.a_of_chi = a_of_chi;
	data.acc_chi = gsl_interp_accel_alloc();
	data.acc_z = gsl_interp_accel_alloc();
	data.zmax = n_of_z->x[n_of_z->size-1];

	// Limit of our range.
	double chi_min = 0.0;
	data.chi_max = chi_max;


	// Integration workspace
	gsl_integration_glfixed_table *table = 
		gsl_integration_glfixed_table_alloc((size_t) NGLT);

	// First compute the normalization
	gsl_function F;
	F.params = &data;
	F.function = &n_of_chi;
	double norm = 1.0 / gsl_integration_glfixed(&F, chi_min, chi_max, table);

	// Now do the main integral, chi-by-chi.
	// First, space for the results
	int n_chi = NGLT;
	double Chi[n_chi];
	double W[n_chi];

	// The evaluation points separation, and set up the integrator
	double delta_chi = (chi_max-chi_min)/ (n_chi - 1);
	F.function = &w_of_chi_integrand;

	gsl_error_handler_t * old_error_handler = gsl_set_error_handler_off();
	int error_status=0;
	// Loop through samples to be evenly spaced in chi
	for(int i=0; i<n_chi; i++)
	{
		// Get chi at this evaluation point
		double chi = delta_chi * i;
		Chi[i] = chi;
		// We need it during the integration, so save.
		data.chi = chi;
		// ingredient for the prefactor chi/a
		double a;
		int err = gsl_spline_eval_e(a_of_chi, chi, data.acc_chi, &a);
		if (err) {error_status=1; break;} 
		// and calculate the integral
		W[i] = norm * (chi/a) * gsl_integration_glfixed(&F, chi, chi_max, table);
	}
	gsl_set_error_handler(old_error_handler);

	// Convert the static vectors into a spline and return
	gsl_spline * output;
	if (error_status) output = NULL;
	else output = spline_from_arrays(n_chi, Chi, W);

	// Tidy up
	gsl_integration_glfixed_table_free(table);
	gsl_interp_accel_free(data.acc_chi);
	gsl_interp_accel_free(data.acc_z);

	// Finish
	return output;
}




gsl_spline * cmb_wl_kappa_kernel(double chi_max, double chi_star, gsl_spline * a_of_chi)
{

	int n_chi = NGLT;
	double chi_min = 0.0;
	int error_status = 0;
	// The evaluation points separation, and set up the integrator
	double delta_chi = (chi_max-chi_min)/ (n_chi - 1);

	double W[n_chi];
	double Chi[n_chi];

	// Loop through samples to be evenly spaced in chi
	for(int i=0; i<n_chi; i++)
	{
		// Get chi at this evaluation point
		double chi = delta_chi * i;
		Chi[i] = chi;
		double a;
		int err = gsl_spline_eval_e(a_of_chi, chi, NULL, &a);
		if (err) {error_status=1; break;} 
		// and calculate the integral
		W[i] = (chi/a) * (chi_star-chi)/chi_star;
	}

	// Convert the static vectors into a spline and return
	gsl_spline * output;
	if (error_status) output = NULL;
	else output = spline_from_arrays(n_chi, Chi, W);

	// Finish
	return output;
}

