#include "stdio.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
//include "limber.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_sf_bessel.h>

// The integrated C_ell are in general allowed to be zero or negative if
// they describe cross-correlations. We use these statuses to describe
// errors where a log was requested also.
// LIMBER_STATUS_NEGATIVE is probably always an error
// but LIMBER_STATUS_ZERO is not necessarily
#define LIMBER_STATUS_OK 0
#define LIMBER_STATUS_ZERO 1
#define LIMBER_STATUS_NEGATIVE 2
#define LIMBER_STATUS_ERROR 3

// These are the options you can set for
// the Limber integrator.
typedef struct cl_config{
	int n_ell;  // Number of ell values you want in the spline
	int * ell;  // The chosen ell values you want
	double prefactor; //Scaling prefactor
    int status; // did everything go okay?
    double delta_range_factor;
    double log_nu_range;
    double absolute_tolerance;
    double relative_tolerance;
} cl_config;

int cl_integral(cl_config * config, gsl_spline * WX_red, 
			    gsl_spline * WY_red, gsl_spline2d * P, gsl_spline * D_chi,
			    double * cl_out, double * cl_err_out);