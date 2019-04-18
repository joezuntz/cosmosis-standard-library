#include "stdio.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "limber.h"
#include "utils.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_spline.h>

#if defined(_OPENMP)
#include "omp.h"
#endif

// This is a workspace size for the gsl integrator
#define LIMBER_FIXED_TABLE_SIZE 4096
#define LIMBER_NTHREAD_MAX 128


// data that is passed into the integrator
// This is everything we need to compute the
// integrand
typedef struct IntegrandData{
	double chimin;
	double chimax;
	double ell;
	gsl_spline * WX;
	gsl_spline * WY;
	gsl_spline2d * P;
	gsl_interp_accel * accelerator_wx;
	gsl_interp_accel * accelerator_wy;
	gsl_interp_accel * accelerator_px;
	gsl_interp_accel * accelerator_py;	
	double K;
} IntegrandData;

// data that is passed into the integrator
// This is everything we need to compute the
// integrand
typedef struct SigmaCritIntegrandData{
	double chimin;
	double chimax;
	gsl_spline * WX;
	gsl_spline * WY;
	gsl_interp_accel * accelerator_wx;
	gsl_interp_accel * accelerator_wy;
} SigmaCritIntegrandData;


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

// the Limber integrand: 
// P(nu/chi, chi) / ( chi * D^2(chi) ) * 
// { Wr_X(chi) * Wr_Y(chi)  //normal Limber up to here
// - 0.5 * chi^2/nu^2 * [Wr_X''(chi) * Wr_Y(chi) + Wr_Y''(chi) * Wr_X(chi)] } //extended part
// where Wr_X(chi) is D(chi)*W_X(chi)/sqrt(chi).
static double limber_integrand(double chi, void * data_void)
{
	IntegrandData * data = (IntegrandData*) data_void;
	// Return 0 if outside range, for convenience.
	// Important to ensure that ranges are wide enough.
	if(chi < data->chimin || chi > data->chimax) return 0.0;

	// Get W^X(chi) and W^Y(chi)
    double wx, wy;

    int status = gsl_spline_eval_e(data->WX,chi,data->accelerator_wx, &wx);


	// A short-cut - if the two splines are the same we do not need to 
	// do the interpolation twice
	if (data->WX==data->WY){
        wy = wx;
    }
    else{
        status |= gsl_spline_eval_e(data->WY,chi,data->accelerator_wy, &wy);
    }

    if (status) return 0.0;

	if ((wx==0) || (wy==0)) return 0.0;

	double nu = (data->ell+0.5);
	double k = nu / f_K(data->K, chi);

	double p;
	status |= gsl_spline2d_eval_e(data->P, chi, log(k), data->accelerator_px, data->accelerator_py, &p);

	if (status) 
	{
		return 0.0;
	}

    double result = wx * wy * p / chi / chi;

    // This is a total hack because we really really need to rewrite things
    // cleanly.
    result *= pow( chi/f_K(data->K, chi), 2);
    return result;
}

double get_kernel_peak(gsl_spline * WX, gsl_spline * WY, int n_chi)
{
  double chimin_x = limber_gsl_spline_min_x(WX);
  double chimin_y = limber_gsl_spline_min_x(WY);
  double chimax_x = limber_gsl_spline_max_x(WX);
  double chimax_y = limber_gsl_spline_max_x(WY);
  double chimin = chimin_x>chimin_y ? chimin_x : chimin_y;
  double chimax = chimax_x<chimax_y ? chimax_x : chimax_y;
  double dchi = (chimax - chimin)/n_chi;
  double chi_peak = chimin;
  double chi;
  double kernel_val=0.;
  double kernel_peak=0.;
  for (int i_chi=0; i_chi<=n_chi; i_chi++){
    chi=chimin+i_chi*dchi;
    kernel_val = gsl_spline_eval(WX,chi,NULL) * gsl_spline_eval(WY,chi,NULL) / chi / chi;
    if (kernel_val>kernel_peak){
      kernel_peak = kernel_val;
      chi_peak = chi;
    }
  }
  // printf("chi_peak = %f\n",chi_peak);
  return chi_peak;
}


static double sigma_crit_integrand(double chi,  void * data_void)
{
  IntegrandData * data = (IntegrandData*) data_void;
  // Return 0 if outside range, for convenience.                                                                                                                            
  // Important to ensure that ranges are wide enough.                                                                                                                       
  if(chi < data->chimin || chi > data->chimax) return 0.0;

  // Get W^X(chi) and W^Y(chi)                                                                                                                                              
  double wx = gsl_spline_eval(data->WX,chi,data->accelerator_wx);
  double wy;
  // A short-cut - if the two splines are the same we do not need to                                                                                                        
  // do the interpolation twice                                                                                                                                             
  if (data->WX==data->WY) wy = wx;
  else wy = gsl_spline_eval(data->WY,chi,data->accelerator_wy);
  if (wx==0 || wy==0) return 0.0;

  double result = wx * wy / chi;
  return result;
}

static double sigma_crit_chi_integrand(double chi, void * data_void)
{
  IntegrandData * data = (IntegrandData*) data_void;
  // Return 0 if outside range, for convenience.                                                                                                                            
  // Important to ensure that ranges are wide enough.                                                                                                                       
  if(chi < data->chimin || chi > data->chimax) return 0.0;

  // Get W^X(chi) and W^Y(chi)                                                                                                                                              
  double wx = gsl_spline_eval(data->WX,chi,data->accelerator_wx);
  double wy;
  // A short-cut - if the two splines are the same we do not need to                                                                                                        
  // do the interpolation twice                                                                                                                                             
  if (data->WX==data->WY) wy = wx;
  else wy = gsl_spline_eval(data->WY,chi,data->accelerator_wy);
  if (wx==0 || wy==0) return 0.0;
  double result = wx * wy;
  return result;
    }

void sigma_crit(gsl_spline * WX, gsl_spline * WY, double * sigma_crit, double * chi_mean)
{
	SigmaCritIntegrandData data;


	double sigma_crit_value;
	double chimin_x = limber_gsl_spline_min_x(WX);
	double chimin_y = limber_gsl_spline_min_x(WY);
	double chimax_x = limber_gsl_spline_max_x(WX);
	double chimax_y = limber_gsl_spline_max_x(WY);

	// Set up the workspace and inputs to the integrator.
	// Not entirely sure what the table is.
	gsl_function F;
	F.function = sigma_crit_integrand;
	F.params = &data;

	// The workspace for the integrator.
	gsl_integration_glfixed_table *table = 
	   gsl_integration_glfixed_table_alloc((size_t) LIMBER_FIXED_TABLE_SIZE);

	data.chimin = chimin_x>chimin_y ? chimin_x : chimin_y;
	data.chimax = chimax_x<chimax_y ? chimax_x : chimax_y;
	data.WX = WX;
	data.WY = WY;
	data.accelerator_wx = gsl_interp_accel_alloc();
	data.accelerator_wy = gsl_interp_accel_alloc();

	printf("Doing sigma_crit integral\n");

	sigma_crit_value = gsl_integration_glfixed(&F, data.chimin, data.chimax, table);
	*sigma_crit = sigma_crit_value;

	printf("Sigma_crit: %f\n",sigma_crit_value);
	
	F.function = sigma_crit_chi_integrand;
	*chi_mean = gsl_integration_glfixed(&F, data.chimin, data.chimax, table);
	*chi_mean /= sigma_crit_value;

	printf("chi_mean: %f\n",*chi_mean);	

	// Tidy up
	gsl_interp_accel_free(data.accelerator_wx);
	gsl_interp_accel_free(data.accelerator_wy);
	gsl_integration_glfixed_table_free(table);	
}



static
void limber_gsl_fallback_integrator(gsl_integration_workspace * W, gsl_integration_glfixed_table *table, gsl_function * F, double chimin, double chimax, 
	double abstol, double reltol, double * c_ell, double * error){

	// Only one warning per process
	static int fallback_warning_given = 0;

	// Deactivate error handling - if this one fails we will fall back to a more reliable but slower integrator
	gsl_error_handler_t * old_handler = gsl_set_error_handler_off();

	// Try the fast but flaky integrator.
	int status = gsl_integration_qag(F, chimin, chimax, abstol, reltol, LIMBER_FIXED_TABLE_SIZE, GSL_INTEG_GAUSS61, W, c_ell, error);

	// Restore the old error handler
	gsl_set_error_handler(old_handler); 

	// If the fast integrator failed fall back to the old one.
	if (status){
		IntegrandData * data = (IntegrandData*) F->params;
		double ell = data->ell;
		if (fallback_warning_given==0){
			fprintf(stderr, "Falling back to the old integrator for ell=%lf (status=%d)\n", ell,status);
			fallback_warning_given=1;
		}
		*c_ell = gsl_integration_glfixed(F,chimin,chimax,table);
	}

}

// The only function in this little library callable from the outside
// world.  The limber_config structure is defined in limber.h but is fairly
// obvious.  The splines and the interpolator need to be functions of 
// chi NOT z.
int limber_integral(limber_config * config, gsl_spline * WX, 
			       gsl_spline * WY, gsl_spline2d * P,
			       double * cl_out)
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
	double chimin_x = limber_gsl_spline_min_x(WX);
	double chimin_y = limber_gsl_spline_min_x(WY);
	double chimax_x = limber_gsl_spline_max_x(WX);
	double chimax_y = limber_gsl_spline_max_x(WY);


	double reltol = config->relative_tolerance;
	double abstol = config->absolute_tolerance;


    static gsl_integration_workspace * workspace[LIMBER_NTHREAD_MAX];
    static gsl_integration_glfixed_table * table[LIMBER_NTHREAD_MAX];

#if defined(_OPENMP)
    int nthread_max = omp_get_max_threads();
#else
    int nthread_max = 1;
#endif
    
    if (nthread_max>LIMBER_NTHREAD_MAX){
        fprintf(stderr, "Tried to use more than %d OpenMP threads - modify limber.c", LIMBER_NTHREAD_MAX);
        exit(1);
    }

    // Get P(k,z) using k=nu/chi.
    // Deactivate error handling 
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();


#pragma omp parallel
{

    //Figure out which thread is running this
#if defined(_OPENMP)
        int thread = omp_get_thread_num();
#else
        int thread = 0;
#endif

    // We maintain an array of per-thread integration workspaces.
    // The downside of this is that we need a fixed size, but we just use
    // a fairly large one.
    if (workspace[thread]==NULL){
        workspace[thread] = gsl_integration_workspace_alloc(LIMBER_FIXED_TABLE_SIZE);
    }
    if (table[thread]==NULL){
        table[thread] = gsl_integration_glfixed_table_alloc((size_t) LIMBER_FIXED_TABLE_SIZE);
    }

	// Take the smallest range since we want both the
	// splines to be valid there.
	// This range as well as all the data needed to compute
	// the integrand is put into a struct to be passed
	// through the integrator to the function above.
    IntegrandData data;    
	data.chimin = chimin_x>chimin_y ? chimin_x : chimin_y;
	data.chimax = chimax_x<chimax_y ? chimax_x : chimax_y;
	data.WX = WX;
	data.WY = WY;
	data.P = P;
	data.accelerator_wx = gsl_interp_accel_alloc();
	data.accelerator_wy = gsl_interp_accel_alloc();
	data.accelerator_px = gsl_interp_accel_alloc();
	data.accelerator_py = gsl_interp_accel_alloc();
    data.K = config->K;

	// Set up the workspace and inputs to the integrator.
	// Not entirely sure what the table is.
	gsl_function F;
	F.function = limber_integrand;
	F.params = &data;



	// loop through ell values according to the input configuration
#pragma omp for schedule(dynamic) nowait
	for (int i_ell = 0; i_ell<config->n_ell; i_ell++){
		double ell = config->ell[i_ell];
		data.ell=ell;
        double c_ell, error;

		// Perform the main integration.
		limber_gsl_fallback_integrator(workspace[thread], table[thread], &F, data.chimin, data.chimax, 
			abstol, reltol, &c_ell, &error);
		//Include the prefactor scaling
		c_ell *= config->prefactor;
		// Record the results into arrays
		cl_out[i_ell] = c_ell;
	}

	// Tidy up
	gsl_interp_accel_free(data.accelerator_wx);
	gsl_interp_accel_free(data.accelerator_wy);
	gsl_interp_accel_free(data.accelerator_px);
	gsl_interp_accel_free(data.accelerator_py);

} // end of parallel region

    // Restore the old error handler
    gsl_set_error_handler(old_handler); 

	// And that's it
	return 0;
}

