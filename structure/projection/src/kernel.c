#include "cosmosis/datablock/c_datablock.h"
#include "gsl/gsl_spline.h"
//#include "shear_shear.h"
#include "utils.h"
#include <stdio.h>
//#include "limber.h"
#include <math.h>
#include "assert.h"
#include "string.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"

// Short-hand names for the sections we will
// be looking at
const char * wl_nz = WL_NUMBER_DENSITY_SECTION;
const char * dist = DISTANCES_SECTION;
const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
const char * ia = INTRINSIC_ALIGNMENT_PARAMETERS_SECTION;
const char * lum = GALAXY_LUMINOSITY_FUNCTION_SECTION;

#define NGLT 1024


typedef struct w_integrand_data{
  gsl_spline * n_of_z;
  gsl_spline * a_of_chi;
  gsl_interp_accel * acc_chi;
  gsl_interp_accel * acc_z;
  double chi;
  double zmax;
  double chi_max;
  double K;
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
  return nchi*f_K(data->K, chis - data->chi) / f_K(data->K,chis);
}

gsl_spline * get_reduced_kernel(gsl_spline * orig_kernel, gsl_spline * growth_of_chi)
{
  int nchi = orig_kernel->size;
  double kernel_out_array[nchi];	      
  double Chi[nchi];
  double chi;
  double growth;
  double chi_max = growth_of_chi->x[growth_of_chi->size-1];
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  for (int i=0; i<nchi; i++)
    {
      chi = (orig_kernel->x[i]);
      Chi[i] = chi;
      if (chi>0. && chi<chi_max){
        growth = gsl_spline_eval(growth_of_chi, chi, acc);
        kernel_out_array[i] = orig_kernel->y[i] * growth / sqrt(chi);
      }
      else {
        kernel_out_array[i] = 0.;
      }
  }
  gsl_interp_accel_free(acc);
  //convert to spline
  return spline_from_arrays(nchi, Chi, kernel_out_array);
}


// this function does not include the factor
// 1.5 omega_m H0^2 in it.  That should be included in the
// prefactor later
gsl_spline * shear_shear_kernel(double chi_max, gsl_spline * n_of_z, 
				gsl_spline * a_of_chi, double K)
{
  // Now do the main integral, chi-by-chi.
  // First, space for the results
  int n_chi = NGLT;
  double Chi[n_chi];
  double W[n_chi];


  // Limit of our range.
  double chi_min = 0.0;

  // Integration workspace
  gsl_integration_glfixed_table *table = 
    gsl_integration_glfixed_table_alloc((size_t) NGLT);

  gsl_error_handler_t * old_error_handler = gsl_set_error_handler_off();

  volatile int error_status=0; //volatile means that one thread can set it for other threads

#pragma omp parallel shared(error_status)
  {
    // Integrand data
    w_integrand_data data;
    data.n_of_z = n_of_z;
    data.a_of_chi = a_of_chi;
    data.acc_chi = gsl_interp_accel_alloc();
    data.acc_z = gsl_interp_accel_alloc();
    data.zmax = n_of_z->x[n_of_z->size-1];
    data.chi_max = chi_max;
    data.K = K;

    // First compute the normalization
    gsl_function F;
    F.params = &data;
    F.function = &n_of_chi;
    double norm = 1.0 / gsl_integration_glfixed(&F, chi_min, chi_max, table);


    // The evaluation points separation, and set up the integrator
    double delta_chi = (chi_max-chi_min)/ (n_chi - 1);
    F.function = &w_of_chi_integrand;

    // Loop through samples to be evenly spaced in chi
#pragma omp for schedule(dynamic)
    for(int i=0; i<n_chi; i++)
      {
        if (error_status) continue;
        // Get chi at this evaluation point
        double chi = delta_chi * i;
        Chi[i] = chi;
        // We need it during the integration, so save.
        data.chi = chi;
        // ingredient for the prefactor chi/a
        double a;
        int err = gsl_spline_eval_e(a_of_chi, chi, data.acc_chi, &a);
        if (err) {
#pragma omp atomic          
          error_status++; 
          continue;
        } 
        // and calculate the integral
        W[i] = norm * (f_K(K,chi)/a) * gsl_integration_glfixed(&F, chi, chi_max, table);
      }

    // Tidy up
    gsl_interp_accel_free(data.acc_chi);
    gsl_interp_accel_free(data.acc_z);
  }
  gsl_integration_glfixed_table_free(table);
  gsl_set_error_handler(old_error_handler);

    // Convert the static vectors into a spline and return
    gsl_spline * output;
    if (error_status) output = NULL;
    else output = spline_from_arrays(n_chi, Chi, W);

  // Finish
  return output;
}


/*
  Read z, N(z) values from the datablock and initialize a spline to
  use in the Limber integral. 
*/
gsl_spline * get_named_nchi_spline(c_datablock * block, const char * section, int bin, double *z,
				   gsl_spline * a_of_chi, gsl_spline * chi_of_z)
{
  int status = 0;
  double * N;
  int n;
  char name[20];
  snprintf(name, 20, "bin_%d", bin);

  status |= c_datablock_get_double_array_1d(block, section, name, &N, &n);
  double Chi[n];

  for (int i=0; i<n; i++){
    double chi = gsl_spline_eval(chi_of_z, z[i], NULL);
    double da_dchi = gsl_spline_eval_deriv(a_of_chi, chi, NULL);
    Chi[i] = chi;
    N[i] *= -(1+z[i])*(1+z[i])*da_dchi;
  }

  gsl_spline * n_of_chi = spline_from_arrays(n, Chi, N);
  free(N);
  return n_of_chi;

}

gsl_spline * get_nchi_spline(c_datablock * block, int bin, double *z,
			     gsl_spline * a_of_chi, gsl_spline * chi_of_z)
{
  return get_named_nchi_spline(block, wl_nz, bin, z, a_of_chi, chi_of_z);
}

gsl_spline * get_named_w_spline(c_datablock * block, const char * section, int bin, double * z,
				double chi_max, gsl_spline * a_of_chi_spline)
{
  char name[20];
  double * n_of_z;
  int nz1;
  int status = 0;

  // Load n(z)
  snprintf(name, 20, "bin_%d", bin);
  status |= c_datablock_get_double_array_1d(block, section, name, &n_of_z,
					    &nz1);

  // Make a spline from n(z)
  gsl_spline * n_of_z_spline = spline_from_arrays(nz1, z, n_of_z);

  double K;
  status |= c_datablock_get_double_default(block, COSMOLOGICAL_PARAMETERS_SECTION, "K", 0.0, &K);

  // Do the main computation
  gsl_spline * W = shear_shear_kernel(chi_max, n_of_z_spline,
				      a_of_chi_spline, K);

  // tidy up bin-specific data
  gsl_spline_free(n_of_z_spline);
  free(n_of_z);
  return W;
}

/*
  Read z, N(z) values from the datablock and initialize a spline
  of the lensing kernel W to use in the Limber integral.
*/
gsl_spline * get_w_spline(c_datablock * block, int bin, double * z,
			  double chi_max, gsl_spline * a_of_chi_spline)
{
  return get_named_w_spline(block, wl_nz, bin, z, chi_max, a_of_chi_spline);
}

int get_wchi_array(c_datablock * block, const char * nz_section,
		   int bin, double * z, double chi_max, gsl_spline * a_of_chi_spline, double * chi_out, 
		   int nchi_out, double * wchi_out)
{
  // Do the main computation
  gsl_spline * W = get_named_w_spline(block, nz_section, bin, z, chi_max, a_of_chi_spline);

  //Now interpolate to output chi
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  for (int i_chi = 0; i_chi<nchi_out; i_chi++){
    wchi_out[i_chi] = gsl_spline_eval(W, chi_out[i_chi], acc);
  }

  // tidy up bin-specific data
  gsl_spline_free(W);
  return 0;
}



gsl_spline * cmb_wl_kappa_kernel(double chi_max, double chi_star, gsl_spline * a_of_chi, double K)
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
      W[i] = f_K(K,chi) / a * f_K(K,chi_star-chi) / f_K(K,chi_star);
    }

  // Convert the static vectors into a spline and return
  gsl_spline * output;
  if (error_status) output = NULL;
  else output = spline_from_arrays(n_chi, Chi, W);

  // Finish
  return output;
}
