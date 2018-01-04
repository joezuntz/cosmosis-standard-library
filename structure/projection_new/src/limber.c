#include "stdio.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include "limber.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_spline.h>

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
typedef struct ExtIntegrandData{
	double chimin;
	double chimax;
	double ell;
	gsl_spline * WX_red;
	gsl_spline * WY_red;
    gsl_spline * D_chi;
	gsl_spline2d * P;
	gsl_interp_accel * accelerator_wx;
	gsl_interp_accel * accelerator_wy;
	gsl_interp_accel * accelerator_px;
	gsl_interp_accel * accelerator_py;
	gsl_interp_accel * accelerator_growth;
    int ext_order;
} ExtIntegrandData;


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

// the limber integrand: 
// P(nu/chi, chi) * W_X(chi) * W_Y(chi) / chi^2
static double limber_integrand(double chi, void * data_void)
{
	ExtIntegrandData * data = (ExtIntegrandData*) data_void;
	// Return 0 if outside range, for convenience.
	// Important to ensure that ranges are wide enough.
	if(chi < data->chimin || chi > data->chimax) return 0.0;

	// Get W^X(chi) and W^Y(chi)
	double Wr_X = gsl_spline_eval(data->WX_red,chi,data->accelerator_wx);
	double Wr_Y;
	// A short-cut - if the two splines are the same we do not need to 
	// do the interpolation twice
	if (data->WX_red==data->WY_red) Wr_Y = Wr_X;
	else Wr_Y = gsl_spline_eval(data->WY_red,chi,data->accelerator_wy);
	if (Wr_X==0 || Wr_Y==0) return 0.0;

	double nu = (data->ell+0.5);
	double k = nu / chi;

	// Get P(k,z) using k=nu/chi.
	// Deactivate error handling 
	gsl_error_handler_t * old_handler = gsl_set_error_handler_off();

	double p;
	int status = gsl_spline2d_eval_e(data->P, chi, log(k), data->accelerator_px, data->accelerator_py, &p);

	// Restore the old error handler
	gsl_set_error_handler(old_handler); 
	if (status) 
	{
		//printf(stderr,'spline failed, for now, assume this is because k was out of range and return 0...should probably be more careful here though.');
		return 0.0;
	}

	// Integrand result.
	double result = p * Wr_X * Wr_Y / (chi*chi);
	return result;
}

// the extended limber integrand: 
// P(nu/chi, chi) / ( chi * D^2(chi) ) * 
// { Wr_X(chi) * Wr_Y(chi) * D^2(chi) } //normal Limber up to here
// - chi^2/nu^2 * [Wr_X''(chi) * Wr_Y(chi) + Wr_Y''(chi) * Wr_X(chi)] } //extended part
// where Wr_X(chi) is D(chi)*W_X(chi)/sqrt(chi).
static double ext_limber_integrand(double chi, void * data_void)
{
	ExtIntegrandData * data = (ExtIntegrandData*) data_void;
	// Return 0 if outside range, for convenience.
	// Important to ensure that ranges are wide enough.
	if(chi < data->chimin || chi > data->chimax) return 0.0;

	// Get W^X(chi) and W^Y(chi)
	double Wr_X = gsl_spline_eval(data->WX_red,chi,data->accelerator_wx);
	double Wr_Y;
	// A short-cut - if the two splines are the same we do not need to 
	// do the interpolation twice
	if (data->WX_red==data->WY_red) Wr_Y = Wr_X;
	else Wr_Y = gsl_spline_eval(data->WY_red,chi,data->accelerator_wy);
	if (Wr_X==0 || Wr_Y==0) return 0.0;

	double nu = (data->ell+0.5);
	double k = nu / chi;
	// 
	double growth = gsl_spline_eval(data->D_chi, chi, data->accelerator_growth);
	//if doing extended limber, get kernel 2nd derivatives
	double d2Wr_x, d2Wr_y;
	if (data->ext_order > 0){
	    double d2Wr_x = gsl_spline_eval_deriv2(data->WX_red, chi, data->accelerator_wx);
	    double d2Wr_y;
    	if (data->WX_red==data->WY_red) d2Wr_y = d2Wr_x;
	    else d2Wr_y = gsl_spline_eval_deriv2(data->WY_red, chi, data->accelerator_wy);
	}

	// Get P(k,z) using k=nu/chi.
	// Deactivate error handling 
	gsl_error_handler_t * old_handler = gsl_set_error_handler_off();

	double p;
	int status = gsl_spline2d_eval_e(data->P, chi, log(k), data->accelerator_px, data->accelerator_py, &p);

	// Restore the old error handler
	gsl_set_error_handler(old_handler); 
	if (status) 
	{
		//printf(stderr,'spline failed, for now, assume this is because k was out of range and return 0...should probably be more careful here though.');
		return 0.0;
	}

	// Integrand result.
	p /= (chi * growth * growth);
	double result = Wr_X * Wr_Y;  //zeroth order
	if (data->ext_order>0){
    	double x2y1 =  - d2Wr_x * Wr_Y; //second order
	    double x1y2 =  - d2Wr_y * Wr_X; //second order
	    result += ( x2y1 + x1y2 ) * chi * chi / (nu * nu) ;
	}
	result *= p;
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
int limber_integral(limber_config * config, gsl_spline * WX_red, 
			       gsl_spline * WY_red, gsl_spline2d * P, gsl_spline * D_chi,
			       int ext_order, double * cl_out)
{

    config->status = LIMBER_STATUS_ERROR;
    int any_parameter_error=0;
    if (WX_red==NULL){
        fprintf(stderr, "NULL WX parameter in limber_integral\n");
        any_parameter_error = 1;
    }
    if (WY_red==NULL){
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
	ExtIntegrandData data;
	double chimin_x = limber_gsl_spline_min_x(WX_red);
	double chimin_y = limber_gsl_spline_min_x(WY_red);
	double chimax_x = limber_gsl_spline_max_x(WX_red);
	double chimax_y = limber_gsl_spline_max_x(WY_red);
	double c_ell, error;

	// Workspaces for the main and falback integrators.
	// Static, so only allocated once as it has a fixed size.
	setup_integration_workspaces();

	double reltol = config->relative_tolerance;
	double abstol = config->absolute_tolerance;

	// Take the smallest range since we want both the
	// splines to be valid there.
	// This range as well as all the data needed to compute
	// the integrand is put into a struct to be passed
	// through the integrator to the function above.
	data.chimin = chimin_x>chimin_y ? chimin_x : chimin_y;
	data.chimax = chimax_x<chimax_y ? chimax_x : chimax_y;
	data.WX_red = WX_red;
	data.WY_red = WY_red;
	data.P = P;
	data.accelerator_wx = gsl_interp_accel_alloc();
	data.accelerator_wy = gsl_interp_accel_alloc();
	data.accelerator_px = gsl_interp_accel_alloc();
	data.accelerator_py = gsl_interp_accel_alloc();
	data.accelerator_growth = gsl_interp_accel_alloc();
	data.D_chi = D_chi;
	data.ext_order = ext_order;

	// Set up the workspace and inputs to the integrator.
	// Not entirely sure what the table is.
	gsl_function F;
	if (ext_order>0) {
		F.function = ext_limber_integrand;
	}
	else {
		F.function = limber_integrand;
	}
	F.params = &data;

	// save ell values for return
	double ell_vector[config->n_ell];

	// loop through ell values according to the input configuration
	for (int i_ell = 0; i_ell<config->n_ell; i_ell++){
		double ell = config->ell[i_ell];
		data.ell=ell;

		// Perform the main integration.
		limber_gsl_fallback_integrator(&F, data.chimin, data.chimax, 
			abstol, reltol, &c_ell, &error);
		//Include the prefactor scaling
		c_ell *= config->prefactor;
		// Record the results into arrays
		cl_out[i_ell] = c_ell;
		ell_vector[i_ell] = ell;
	}

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
