
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
const char * ppf_section = POST_FRIEDMANN_PARAMETERS_SECTION;

typedef enum spectrum_type_t {
	shear_shear,
	shear_intrinsic,
	intrinsic_intrinsic,
	matter_matter
} spectrum_type_t;

typedef struct shear_spectrum_config {
	int n_ell;
	double ell_min;
	double ell_max;
	bool intrinsic_alignments;
	bool matter_spectra;
} shear_spectrum_config;


void * setup(c_datablock * options){

	shear_spectrum_config * config = malloc(sizeof(shear_spectrum_config));
	int status = 0;
	status |= c_datablock_get_int(options, OPTION_SECTION, "n_ell", &(config->n_ell));
	status |= c_datablock_get_double(options, OPTION_SECTION, "ell_min", &(config->ell_min));
	status |= c_datablock_get_double(options, OPTION_SECTION, "ell_max", &(config->ell_max));

	status |= c_datablock_get_bool(options, OPTION_SECTION, "intrinsic_alignments", 
		&(config->intrinsic_alignments));

	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "matter_spectra", 
		false, &(config->matter_spectra));

	if (config->intrinsic_alignments) {
		fprintf(stderr, "\n\n***********\nIMPORTANT***********\n\n");
		fprintf(stderr, "No one has yet thought through how to do intrinsic");
		fprintf(stderr, "alignments in modified gravity.  So you will need");
		fprintf(stderr, "to set intrinsic_alignments=F in the ini file, sorry.\n");
		fprintf(stderr, "\n\n**********************************\n\n");
		exit(1);
	}

	if (status){
		fprintf(stderr, "Please specify n_ell, ell_min, and ell_max in the shear spectra module.\n");
		fprintf(stderr, "And set intrinsic_alignments=T or F\n");
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

const char * choose_output_section(spectrum_type_t spectrum_type, 
	shear_spectrum_config * config)
{
	switch(spectrum_type){
		case matter_matter:
			return "matter_cl";
			break;
		case shear_intrinsic:
			return 	SHEAR_CL_GI_SECTION;
			break;
		case intrinsic_intrinsic:
			return SHEAR_CL_II_SECTION;
			break;
		case shear_shear:
			return (config->intrinsic_alignments ? SHEAR_CL_GG_SECTION : SHEAR_CL_SECTION);
			break;
		default:
			return NULL;
	}
}


int choose_configuration(c_datablock * block, spectrum_type_t spectrum_type, 
	limber_config * lc, shear_spectrum_config * config)
{

	// First we need to decide whether to use logs in the splines.
	// This is more accurate, but we cannot use it for cross-spectra
	// since they can go negative.

	if (spectrum_type==shear_shear || spectrum_type==matter_matter || 
		spectrum_type==intrinsic_intrinsic){
		lc->xlog=true;
		lc->ylog=true;
	}
	else{
		lc->xlog=false;
		lc->ylog=false;
	}


	// Now we need to set up the ell range we wish to compute.
	// We need the number of ell values and the ell min and max.
	lc->n_ell = config->n_ell;

	// We log space the actual ell values
	lc->ell = malloc(sizeof(double)*lc->n_ell);
	double alpha = (log(config->ell_max) - log(config->ell_min)) / (config->n_ell - 1);
	for (int i=0; i<lc->n_ell; i++) lc->ell[i] = config->ell_min*exp(alpha * i);


	// We just fix these for now.
	lc->relative_tolerance = 1e-3;
	lc->absolute_tolerance = 1e-5;


	//Scalings

	// This scaling value is the bit that goes in front
	// of the kernels. For n(chi) this is 1.0 and for 
	// shear W(chi) it's this more complex number.

	// shear-shear picks up two copies of this scaling,
	// shear-intrinsic picks up one, 
	// and the other spectra just have a unity prefactor.

	double omega_m;
	int status = c_datablock_get_double(block, cosmo, "omega_m", &omega_m);
	// in (Mpc/h)^-2
	// We leave the factor of h in because our P(k) 
	// should have it in.  That's why there is a 100 in here.
	const double c_kms = 299792.4580; //km/s
	double shear_scaling = 1.5 * (100.0*100.0)/(c_kms*c_kms) * omega_m; 

	switch(spectrum_type){
		case matter_matter:
			lc->prefactor = 1.0;
			break;
		case shear_intrinsic:
			lc->prefactor = shear_scaling;
			break;
		case intrinsic_intrinsic:
			lc->prefactor = 1.0;
			break;
		case shear_shear:
			lc->prefactor = shear_scaling*shear_scaling;
			break;
		default:
			lc->prefactor = 1.0/0.0;
			return 1;
	}
	return 0;
}



int choose_kernels(spectrum_type_t spectrum_type, int nbin, gsl_spline * W[nbin], 
	gsl_spline * Nchi[nbin], gsl_spline ** K1, gsl_spline ** K2)
{
	// The different spectra use different kernels to integrate over:
	// Here we choose which to use (we do it like this so we can use
	// the same code for the different spectra)

	switch(spectrum_type){
		case matter_matter:
			for(int b=0; b<nbin; b++) {K1[b] = Nchi[b]; K2[b] = Nchi[b];};
			break;
		case shear_intrinsic:
			for(int b=0; b<nbin; b++) {K1[b] = Nchi[b]; K2[b] = W[b];};
			break;
		case intrinsic_intrinsic:
			for(int b=0; b<nbin; b++) {K1[b] = Nchi[b]; K2[b] = Nchi[b];};
			break;
		case shear_shear:
			for(int b=0; b<nbin; b++) {K1[b] = W[b]; K2[b] = W[b];};
			break;
		default:
			return 1;
	}
	return 0;

}

int choose_bin2_max(spectrum_type_t spectrum_type, int nbin, int bin1)
{
	if (spectrum_type==shear_shear || spectrum_type==matter_matter || 
		spectrum_type==intrinsic_intrinsic){
		return bin1;
	}
	else{
		return nbin;
	}

}



int compute_spectra(c_datablock * block, int nbin, spectrum_type_t spectrum_type,
	gsl_spline * W[nbin], gsl_spline * Nchi[nbin], Interpolator2D * PK, 
	shear_spectrum_config * config)
{
	// This is a unified interface for generating different 
	// spectra.

	// We need:
	//  - to choose a section to save the output into
	//  - to configure the ell ranges and other options
	//  - to choose the W/N kernels to integrate over.

	limber_config lc; 
	lc.ell = NULL;
	gsl_spline * K1[nbin];
	gsl_spline * K2[nbin];

	const char * section = choose_output_section(spectrum_type, config);
	int status = choose_configuration(block, spectrum_type, &lc, config);
	status |= choose_kernels(spectrum_type, nbin, W, Nchi, &K1[0], &K2[0]);

	// If any of these have gone wrong we quit
	// after cleaning up any allocated memory
	if (section==NULL) status = 1;
	if (status) {
		free(lc.ell);
		return status;
	}

	// First save the ell values and number of bins
	status |= c_datablock_put_double_array_1d(block, section, "ell", lc.ell, lc.n_ell);
	status |= c_datablock_put_int(block, section, "nbin", nbin);

	for (int bin1=1; bin1<=nbin; bin1++){
		// We also need to choose the max value of bin 2. 
		// For auto-spectra we don't compute both bin_i_j and bin_j_i,
		// whereas for cross-spectra we do
		int bin2_max = choose_bin2_max(spectrum_type, nbin, bin1);
		for (int bin2=1; bin2<=bin2_max; bin2++){
			gsl_spline * c_ell = limber_integral(&lc, K1[bin1-1], K2[bin2-1], PK);
			if (c_ell == NULL)  return 1;
			int status = save_c_ell(block, section, bin1, bin2,  c_ell, &lc);
			gsl_spline_free(c_ell);
			if (status) return status;
		}
	}
	free(lc.ell);
	return status;
}

double shear_shear_mg_scaling(double k,double z,double P, void* args)
{
	Interpolator2D * D = (Interpolator2D*) args;
	double d_kz = interp_2d(k, z, D);
	// If we are outside the range where d_kz is defined then
	// this function will return zero.
	// in this case just revert to the unmodified P(k,z)
	if (d_kz==0) return P;
	return P*d_kz*d_kz;
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
	Interpolator2D * MG_D = load_interpolator_chi(
		block, chi_of_z_spline, ppf_section, "k_h", "z", "D");

	if (MG_D==NULL) return 1;


	// Get the P(k) we need
	Interpolator2D * PK = load_interpolator_chi_function(
		block, chi_of_z_spline, MATTER_POWER_NL_SECTION, "k_h", "z", "P_k",
		shear_shear_mg_scaling, (void*) MG_D
		);

	destroy_interp_2d(MG_D);

	if (PK==NULL) return 1;

	Interpolator2D * PK_GI = NULL;
	Interpolator2D * PK_II = NULL;
	if (config->intrinsic_alignments){
		PK_II = load_interpolator_chi(
			block, chi_of_z_spline, "intrinsic_alignment_parameters", "k_h", "z", "P_II");

		PK_GI = load_interpolator_chi(
			block, chi_of_z_spline, "intrinsic_alignment_parameters", "k_h", "z", "P_GI");
		
		if (PK_II==NULL) return 2;
		if (PK_GI==NULL) return 3;

	}

	if (!PK) {
		free(chi);
		free(a);
		free(z);
		gsl_spline_free(a_of_chi_spline);
		gsl_spline_free(chi_of_z_spline);
		return 1;
	}


	status |= compute_spectra(block, nbin, shear_shear,
		W_splines, Nchi_splines, PK, config);
	

	if (config->matter_spectra){
		status |= compute_spectra(block, nbin, matter_matter,
			W_splines, Nchi_splines, PK, config);
	}

	if (config->intrinsic_alignments){

		status |= compute_spectra(block, nbin, intrinsic_intrinsic,
			W_splines, Nchi_splines, PK_II, config);

		status |= compute_spectra(block, nbin, shear_intrinsic,
			W_splines, Nchi_splines, PK_GI, config);
	
		destroy_interp_2d(PK_II);
		destroy_interp_2d(PK_GI);
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
        return 0;
}
