#include "stdio.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "limber.h"

// This is a workspace size for the gsl integrator
#define LIMBER_FIXED_TABLE_SIZE 4096

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
	gsl_interp_accel * accelerator_x;
	gsl_interp_accel * accelerator_y;
} IntegrandData;

// data that is passed into the integrator
// This is everything we need to compute the
// integrand
typedef struct SigmaCritIntegrandData{
	double chimin;
	double chimax;
	gsl_spline * WX;
	gsl_spline * WY;
	gsl_interp_accel * accelerator_x;
	gsl_interp_accel * accelerator_y;
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


// the integrand W_X(chi) * W_Y(chi) * P(ell/chi, chi) / chi^2
static double integrand(double chi, void * data_void)
{
	IntegrandData * data = (IntegrandData*) data_void;
	// Return 0 if outside range, for convenience.
	// Important to ensure that ranges are wide enough.
	if(chi < data->chimin || chi > data->chimax) return 0.0;

	// Get W^X(chi) and W^Y(chi)
	double wx = gsl_spline_eval(data->WX,chi,data->accelerator_x);
	double wy;
	// A short-cut - if the two splines are the same we do not need to 
	// do the interpolation twice
	if (data->WX==data->WY) wy = wx;
	else wy = gsl_spline_eval(data->WY,chi,data->accelerator_y);
	if (wx==0 || wy==0) return 0.0;

	// Get P(k,z) using k=ell/chi.
	// The interp_2d interpolator returns 0 if either 
	// parameter is outside its range
	double k = (data->ell+0.5) / chi;
	double p = interp_2d(k, chi, data->P);
	// printf("p, wx, wy = %le  %le  %le\n", p, wx, wy);
	// Integrand result.
	double result = wx * wy * p / chi / chi;
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
  double wx = gsl_spline_eval(data->WX,chi,data->accelerator_x);
  double wy;
  // A short-cut - if the two splines are the same we do not need to                                                                                                        
  // do the interpolation twice                                                                                                                                             
  if (data->WX==data->WY) wy = wx;
  else wy = gsl_spline_eval(data->WY,chi,data->accelerator_y);
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
  double wx = gsl_spline_eval(data->WX,chi,data->accelerator_x);
  double wy;
  // A short-cut - if the two splines are the same we do not need to                                                                                                        
  // do the interpolation twice                                                                                                                                             
  if (data->WX==data->WY) wy = wx;
  else wy = gsl_spline_eval(data->WY,chi,data->accelerator_y);
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
	data.accelerator_x = gsl_interp_accel_alloc();
	data.accelerator_y = gsl_interp_accel_alloc();

	printf("Doing sigma_crit integral\n");

	sigma_crit_value = gsl_integration_glfixed(&F, data.chimin, data.chimax, table);
	*sigma_crit = sigma_crit_value;

	printf("Sigma_crit: %f\n",sigma_crit_value);
	
	F.function = sigma_crit_chi_integrand;
	*chi_mean = gsl_integration_glfixed(&F, data.chimin, data.chimax, table);
	*chi_mean /= sigma_crit_value;

	printf("chi_mean: %f\n",*chi_mean);	

	// Tidy up
	gsl_interp_accel_free(data.accelerator_x);
	gsl_interp_accel_free(data.accelerator_y);
	gsl_integration_glfixed_table_free(table);	
}


gsl_integration_workspace * W = NULL;
gsl_integration_glfixed_table *table = NULL;

void setup_integration_workspaces(){
	if (W==NULL){
		W = gsl_integration_workspace_alloc(LIMBER_FIXED_TABLE_SIZE);
	}
	if (table==NULL){
		table = gsl_integration_glfixed_table_alloc((size_t) LIMBER_FIXED_TABLE_SIZE);
	}
}


static
void limber_gsl_fallback_integrator(gsl_function * F, double chimin, double chimax, 
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
gsl_spline * limber_integral(limber_config * config, gsl_spline * WX, 
	                 gsl_spline * WY, Interpolator2D * P)
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
    if (config->n_ell==0){
        fprintf(stderr, "Error: n_ell=0 in limber calculation.\n");
        any_parameter_error = 1;
    }
    if (any_parameter_error){
        return NULL;
    }

	config->status = LIMBER_STATUS_OK;


	// Get the appropriate ranges over which to integrate
	// It is assumed that (at least one of) the kernel
	// splines should go to zero in some kind of reasonable
	// place, so we just use the range they specify
	IntegrandData data;
	double chimin_x = limber_gsl_spline_min_x(WX);
	double chimin_y = limber_gsl_spline_min_x(WY);
	double chimax_x = limber_gsl_spline_max_x(WX);
	double chimax_y = limber_gsl_spline_max_x(WY);
	double c_ell, error;

	// Workspaces for the main and falback integrators.
	// Static, so only allocated once as it has a fixed size.
	setup_integration_workspaces();

	double reltol = config->relative_tolerance;
	double abstol = config->absolute_tolerance;
	// double reltol = 0.001;
	// double abstol = 0.00001;
	// printf("TOLS: %le %le\n",reltol,abstol);


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
	data.accelerator_x = gsl_interp_accel_alloc();
	data.accelerator_y = gsl_interp_accel_alloc();

	// Set up the workspace and inputs to the integrator.
	// Not entirely sure what the table is.
	gsl_function F;
	F.function = integrand;
	F.params = &data;

	//gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(2048);

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
		//c_ell = gsl_integration_glfixed(&F,data.chimin,data.chimax,table);
		// New function still attributable to the legacy of Matt Becker's integrator wisdom.
		// gsl_integration_qag(&F, data.chimin, data.chimax, abstol, reltol, LIMBER_FIXED_TABLE_SIZE, GSL_INTEG_GAUSS61, W, table, &c_ell, &error);
		//printf("%d %f %f\n",i_ell,c_ell_old,c_ell);
		limber_gsl_fallback_integrator(&F, data.chimin, data.chimax, 
			abstol, reltol, &c_ell, &error);

		//Include the prefactor scaling
		c_ell *= config->prefactor;

		// Record the results into arrays
		c_ell_vector[i_ell] = c_ell;
		ell_vector[i_ell] = ell;
	}

		// It is often useful to interpolate into the logs of the functions
		// This is optional in the config. We move this outside the main loop
		// since we may have all zeros in the output
		if (config->xlog) {
			for (int i_ell = 0; i_ell<config->n_ell; i_ell++){
				ell_vector[i_ell] = log(ell_vector[i_ell]);
			}
		}
		if (config->ylog){
			for (int i_ell = 0; i_ell<config->n_ell; i_ell++){
				if (c_ell_vector[i_ell]<0){
					config->status = LIMBER_STATUS_NEGATIVE;
				}
					// negative is worse than zero so only set to zero it not already negative
				else if ((c_ell_vector[i_ell]==0) && (config->status<LIMBER_STATUS_ZERO)){
					config->status = LIMBER_STATUS_ZERO;
				}
			}
			// If none of the values are <= 0 then we are okay to go ahead and take the logs.
			if (config->status == LIMBER_STATUS_OK){
				for (int i_ell = 0; i_ell<config->n_ell; i_ell++) c_ell_vector[i_ell] = log(c_ell_vector[i_ell]);
			}

		}

	// Create a spline of the arrays as the output
	gsl_spline * output = gsl_spline_alloc(gsl_interp_akima, (size_t) config->n_ell);
	gsl_spline_init(output, ell_vector, c_ell_vector, (size_t) config->n_ell);

	// Tidy up
	gsl_interp_accel_free(data.accelerator_x);
	gsl_interp_accel_free(data.accelerator_y);

	// These two are not deallocated because they are static and only initialized once.
	// gsl_integration_glfixed_table_free(table);	
	// gsl_integration_workspace_free(W);

	// And that's it
	return output;
}
