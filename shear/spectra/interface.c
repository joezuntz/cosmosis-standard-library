
#include "cosmosis/datablock/c_datablock.h"
#include "gsl/gsl_spline.h"
#include "shear_shear.h"
#include "utils.h"
#include <stdio.h>
#include "limber.h"
#include <math.h>

// Short-hand names for the sections we will 
// be looking at
const char * wl_nz = WL_NUMBER_DENSITY_SECTION;
const char * dist = DISTANCES_SECTION;
const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;

typedef struct shear_spectrum_config {
	int n_ell;
	double ell_min;
	double ell_max;
	bool intrinsic_alignments;
} shear_spectrum_config;


void * setup(c_datablock * options){

	shear_spectrum_config * config = malloc(sizeof(shear_spectrum_config));
	int status = 0;
	status |= c_datablock_get_int(options, OPTION_SECTION, "n_ell", &(config->n_ell));
	status |= c_datablock_get_double(options, OPTION_SECTION, "ell_min", &(config->ell_min));
	status |= c_datablock_get_double(options, OPTION_SECTION, "ell_max", &(config->ell_max));
	status |= c_datablock_get_bool(options, OPTION_SECTION, "intrinsic_alignments", &(config->intrinsic_alignments));

	if (status){
		fprintf(stderr, "Please specify n_ell, ell_min, and ell_max in the shear spectra module.\n");
		exit(status);
	}

	return config;
}

gsl_spline * get_nchi_spline(c_datablock * block, int bin, double *z, 
	gsl_spline * a_of_chi, gsl_spline * chi_of_z)
{
	int status = 0;
	double * N;
	int n;
	char name[20];
	snprintf(name, 20, "bin_%d", bin);

    status |= c_datablock_get_double_array_1d(block, wl_nz, name, &N, &n);
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

gsl_spline * get_w_spline(c_datablock * block, int bin, double * z, 
	double chi_max, gsl_spline * a_of_chi_spline)
{
	char name[20];
	double * n_of_z;
	int nz1;
	int status = 0;

	// Load n(z)
	snprintf(name, 20, "bin_%d", bin);
	status |= c_datablock_get_double_array_1d(block, wl_nz, name, &n_of_z,
		&nz1);

	// Make a spline from n(z)
	gsl_spline * n_of_z_spline = spline_from_arrays(nz1, z, n_of_z);

	// Do the main computation
	gsl_spline * W = shear_shear_kernel(chi_max, n_of_z_spline, 
		a_of_chi_spline);

	// tidy up bin-specific data
	gsl_spline_free(n_of_z_spline);
	free(n_of_z);
	return W;
}



int shear_shear_config(c_datablock * block, limber_config * lc, shear_spectrum_config * config){

	// Compute the prefactor, (1.5 Omega_M H0^2)^2
	// The units of this need to be consistent with 
	// what we 		
	double omega_m;
	int status = c_datablock_get_double(block, cosmo, "omega_m", &omega_m);
	// in (Mpc/h)^-2
	// We leave the factor of h in because our P(k) 
	// should have it in.  That's why there is a 100 in here.
	const double c_kms = 299792.4580; //km/s
	double scaling = 1.5 * (100.0*100.0)/(c_kms*c_kms) * omega_m; 
	// This scaling value is the bit that goes in front
	// of the W(chi) kernel.  So we square it for the 
	// shear-shear since that has two copies

	// Set up the integration rules.
	// We can interpolate in the logs because
	// everything should be positive
	lc->xlog = true;
	lc->ylog = true;
	lc->n_ell = config->n_ell;
	lc->prefactor = scaling*scaling;

	// Log spaced target ell based on input
	lc->ell = malloc(sizeof(double)*lc->n_ell);
	double alpha = (log(config->ell_max) - log(config->ell_min)) / (config->n_ell - 1);
	for (int i=0; i<lc->n_ell; i++) lc->ell[i] = config->ell_min*exp(alpha * i);

	return status;
}



int shear_intrinsic_config(c_datablock * block, limber_config * lc, shear_spectrum_config * config){

	// Compute the prefactor, (1.5 Omega_M H0^2)^2
	// The units of this need to be consistent with 
	// what we 		
	double omega_m;
	int status = c_datablock_get_double(block, cosmo, "omega_m", &omega_m);
	// in (Mpc/h)^-2
	// We leave the factor of h in because our P(k) 
	// should have it in.  That's why there is a 100 in here.
	const double c_kms = 299792.4580; //km/s
	double scaling = 1.5 * (100.0*100.0)/(c_kms*c_kms) * omega_m; 
	// This scaling value is the bit that goes in front
	// of the W(chi) kernel.  So we square it for the 
	// shear-shear since that has two copies

	// Set up the integration rules.
	lc->xlog = false;
	lc->ylog = false;
	lc->n_ell = config->n_ell;
	lc->prefactor = scaling;

	// Log spaced target ell based on input
	lc->ell = malloc(sizeof(double)*lc->n_ell);
	double alpha = (log(config->ell_max) - log(config->ell_min)) / (config->n_ell - 1);
	for (int i=0; i<lc->n_ell; i++) lc->ell[i] = config->ell_min*exp(alpha * i);

	return status;
}




int intrinsic_intrinsic_config(c_datablock * block, limber_config * lc, shear_spectrum_config * config){
	// Set up the integration rules.
	// We can interpolate in the logs because
	// everything should be positive, just as in shear-shear
	lc->xlog = true;
	lc->ylog = true;
	lc->n_ell = config->n_ell;
	lc->prefactor = 1.0;

	// Log spaced target ell based on input
	lc->ell = malloc(sizeof(double)*lc->n_ell);
	double alpha = (log(config->ell_max) - log(config->ell_min)) / (config->n_ell - 1);
	for (int i=0; i<lc->n_ell; i++) lc->ell[i] = config->ell_min*exp(alpha * i);

	return 0;
}

 

static
int save_c_ell(c_datablock * block, const char * section, 
	int bin1, int bin2, gsl_spline * s, limber_config * lc)
{
	char name[64];
	snprintf(name, 64, "bin_%d_%d",bin1,bin2);
	int n_ell = lc->n_ell;
	double c_ell[n_ell];
	for (int i=0; i<n_ell; i++){
		double ell = lc->ell[i];
		if (lc->xlog) ell = log(ell);
		c_ell[i] = gsl_spline_eval(s, ell, NULL);
		if (lc->ylog) c_ell[i] = exp(c_ell[i]);
	}

	int status = c_datablock_put_double_array_1d(block, section, name, c_ell, n_ell);
	status |= c_datablock_put_metadata(block, section, name, "note", "no ell^2 prefactor");
			return status;

}

int shear_intrinsic_spectra(c_datablock * block, int nbin, 
	gsl_spline * W[nbin], gsl_spline * Nchi[nbin], Interpolator2D * PK, shear_spectrum_config * config)
{

	// Get the prefactor
	const char * section = SHEAR_CL_GI_SECTION;

	limber_config lc;
	int status = shear_intrinsic_config(block, &lc, config);
	if (status) {free(lc.ell); return status;}
	status |= c_datablock_put_double_array_1d(block, section, "ell", lc.ell, lc.n_ell);
	status |= c_datablock_put_int(block, section, "nbin", nbin);

	for (int bin1=1; bin1<=nbin; bin1++){
		for (int bin2=1; bin2<=nbin; bin2++){
			gsl_spline * c_ell = limber_integral(&lc, Nchi[bin1-1], W[bin2-1], PK);
			int status = save_c_ell(block, section, bin1, bin2,  c_ell, &lc);
			gsl_spline_free(c_ell);
			if (status) return status;
		}
	}
	free(lc.ell);
	return status;
}	


int shear_shear_spectra(c_datablock * block, int nbin, 
	gsl_spline * W[nbin], Interpolator2D * PK, shear_spectrum_config * config, const char * section)
{

	// Get the prefactor

	limber_config lc;
	int status = shear_shear_config(block, &lc, config);
	if (status) {free(lc.ell); return status;}
	status |= c_datablock_put_double_array_1d(block, section, "ell", lc.ell, lc.n_ell);
	status |= c_datablock_put_int(block, section, "nbin", nbin);

	for (int bin1=1; bin1<=nbin; bin1++){
		for (int bin2=1; bin2<=bin1; bin2++){
			gsl_spline * c_ell = limber_integral(&lc, W[bin1-1], W[bin2-1], PK);
			int status = save_c_ell(block, section, bin1, bin2,  c_ell, &lc);
			gsl_spline_free(c_ell);
			if (status) return status;
		}
	}
	free(lc.ell);
	return status;
}	


int intrinsic_intrinsic_spectra(c_datablock * block, int nbin, 
	gsl_spline * Nchi[nbin], Interpolator2D * PK, shear_spectrum_config * config)
{

	// Get the prefactor
	const char * section = SHEAR_CL_II_SECTION;
	limber_config lc;
	int status = intrinsic_intrinsic_config(block, &lc, config);
	if (status) {free(lc.ell); return status;}
	status |= c_datablock_put_double_array_1d(block, section, "ell", lc.ell, lc.n_ell);
	status |= c_datablock_put_int(block, section, "nbin", nbin);

	for (int bin1=1; bin1<=nbin; bin1++){
		for (int bin2=1; bin2<=bin1; bin2++){
			gsl_spline * c_ell = limber_integral(&lc, Nchi[bin1-1], Nchi[bin2-1], PK);
			int status = save_c_ell(block, section, bin1, bin2,  c_ell, &lc);
			gsl_spline_free(c_ell);
			if (status) return status;
		}
	}
	free(lc.ell);
	return status;
}	



int execute(c_datablock * block, void * config_in)
{
	DATABLOCK_STATUS status=0;
	double * chi;
	double * a;
	double * z;
	int nz1, nz2;
	int nbin;

	shear_spectrum_config * config = (shear_spectrum_config*) config_in;


	// Load the number of bins we 
	status |= c_datablock_get_int(block, wl_nz, "nbin", &nbin);

	// Load z
	status |= c_datablock_get_double_array_1d(block, wl_nz, "z", &z, &nz1);

	// Load another, different z spacing.
	// I am putting this in an array called "a"
	// because I am about to convert it
	status |= c_datablock_get_double_array_1d(block, dist, "z", &a, &nz2);

	// Load chi(z)
	status |= c_datablock_get_double_array_1d(block, dist, "d_m", &chi, &nz2);

	// Reverse ordering so a is increasing - that is what
	// gsl_spline wants
	reverse(a, nz2);
	reverse(chi, nz2);

	// Convert chi from Mpc to Mpc/h
	double h0=0.0;
	status |= c_datablock_get_double(block, cosmo, "h0", &h0);
	for (int i=0; i<nz2; i++) chi[i]*=h0;


	// At the moment "a" is still actually redshift z
	gsl_spline * chi_of_z_spline = spline_from_arrays(nz2, a, chi);

	// Replace z->a
	for (int i=0; i<nz2; i++) a[i] = 1.0/(1+a[i]);

	double chi_max = chi[nz2-1];

	// Make spline of a(chi)
	gsl_spline * a_of_chi_spline = spline_from_arrays(nz2, chi, a);

	// Make the W()
	int error_status=0;
	gsl_spline * W_splines[nbin];
	gsl_spline * Nchi_splines[nbin];
	for (int bin=1; bin<=nbin; bin++){
		W_splines[bin-1] = get_w_spline(block, bin, z, chi_max, a_of_chi_spline);
		Nchi_splines[bin-1] = get_nchi_spline(block, bin, z, a_of_chi_spline, chi_of_z_spline);
		if (W_splines[bin-1]==NULL) error_status=1;
		if (Nchi_splines[bin-1]==NULL) error_status=1;
		
	}
	if (error_status){
		free(chi);
		free(a);
		free(z);
		gsl_spline_free(a_of_chi_spline);
		gsl_spline_free(chi_of_z_spline);
		return 1;		
	}


	// Get the P(k) we need
	Interpolator2D * PK = load_interpolator(
		block, chi_of_z_spline, MATTER_POWER_NL_SECTION, "k_h", "z", "P_k");

	if (PK==NULL) return 1;

	Interpolator2D * PK_GI = NULL;
	Interpolator2D * PK_II = NULL;
	if (config->intrinsic_alignments){
		PK_II = load_interpolator(
			block, chi_of_z_spline, "intrinsic_alignment_parameters", "k_h", "z", "P_II");

		PK_GI = load_interpolator(
			block, chi_of_z_spline, "intrinsic_alignment_parameters", "k_h", "z", "P_GI");
		
		if (PK_II==NULL) return 2;
		if (PK_GI==NULL) return 3;

	}



	// This is P(k,z)
	// We need P(k, chi)
	if (!PK) {
		free(chi);
		free(a);
		free(z);
		gsl_spline_free(a_of_chi_spline);
		gsl_spline_free(chi_of_z_spline);
		return 1;
	}

	// Make the C_ell and save them
	const char * section;
	if (config->intrinsic_alignments){
		section = SHEAR_CL_GG_SECTION;
	}
	else{
		section=SHEAR_CL_SECTION;
	}
	
	status |= shear_shear_spectra(block, nbin, W_splines, PK, config, section);

	if (config->intrinsic_alignments){
		status |= intrinsic_intrinsic_spectra(block, nbin, 
			Nchi_splines, PK_II, config);

		status |= shear_intrinsic_spectra(block, nbin, 
			W_splines, Nchi_splines, PK_GI, config);
	}


	// tidy up global data
	for (int bin=0; bin<nbin; bin++) gsl_spline_free(W_splines[bin]);
	gsl_spline_free(a_of_chi_spline);
	gsl_spline_free(chi_of_z_spline);
	destroy_interp_2d(PK);
	free(chi);
	free(a);
	free(z);

  return status;

}

int cleanup(void * config_in)
{
	// Free the memory that we allocated in the 
	// setup
	shear_spectrum_config * config = (shear_spectrum_config*) config_in;
	free(config);	
}