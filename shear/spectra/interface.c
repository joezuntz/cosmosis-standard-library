
#include "cosmosis/datablock/c_datablock.h"
#include "gsl/gsl_spline.h"
#include "shear_shear.h"
#include "utils.h"
#include <stdio.h>
#include "limber.h"
#include <math.h>
#include "assert.h"
#include "string.h"

// Short-hand names for the sections we will
// be looking at
const char * wl_nz = WL_NUMBER_DENSITY_SECTION;
const char * dist = DISTANCES_SECTION;
const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
const char * ia = INTRINSIC_ALIGNMENT_PARAMETERS_SECTION;
const char * lum = GALAXY_LUMINOSITY_FUNCTION_SECTION;

typedef enum spectrum_type_t {
	// Shape-shape
	shear_shear,
	shear_intrinsic,
	intrinsic_intrinsic,

	// position-position
	matter_matter,
	matter_magnification,
	magnification_magnification,

	// position-shape
	matter_shear,
	matter_intrinsic,
	magnification_intrinsic,
	magnification_shear
} spectrum_type_t;

typedef struct shear_spectrum_config {
	int n_ell;
	double ell_min;
	double ell_max;
	bool verbose;
	bool shear_shear;
	bool intrinsic_alignments;
	bool matter_spectra;
	bool ggl_spectra;
	bool gal_IA_cross_spectra;
	bool galaxy_magnification;
	bool magnification_magnification;
	bool magnification_intrinsic;
	bool magnification_shear;
	
} shear_spectrum_config;

int destroy_and_null(Interpolator2D ** P){
	if (*P) destroy_interp_2d(*P);
	*P = NULL;
        return 0;
}

void * setup(c_datablock * options){

	shear_spectrum_config * config = malloc(sizeof(shear_spectrum_config));
	int status = 0;
	status |= c_datablock_get_int(options, OPTION_SECTION, "n_ell", &(config->n_ell));
	status |= c_datablock_get_double(options, OPTION_SECTION, "ell_min", &(config->ell_min));
	status |= c_datablock_get_double(options, OPTION_SECTION, "ell_max", &(config->ell_max));
	// Selectors for different types of angular power spectra to compute
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "verbose", false,
		&(config->verbose));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "shear_shear", true,
		&(config->shear_shear));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "intrinsic_alignments", false,
		&(config->intrinsic_alignments));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "matter_spectra", false,
		&(config->matter_spectra));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "ggl_spectra", false,
		&(config->ggl_spectra));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "gal_IA_cross_spectra", false,
		&(config->gal_IA_cross_spectra));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "mag_gal_cross_spectra", false,
		&(config->galaxy_magnification));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "mag_mag", false,
			&(config->magnification_magnification));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "mag_IA", false,
			&(config->magnification_intrinsic));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "mag_shear", false,
			&(config->magnification_shear));

	if (status){
		fprintf(stderr, "Please specify n_ell, ell_min, and ell_max in the shear spectra module.\n");
		exit(status);
	}

	return config;
}

gsl_spline * linear_spline_from_arrays(int n, double * x, double *y)
{
	gsl_spline * output = gsl_spline_alloc(gsl_interp_linear, n);
	assert (output!=NULL);
	gsl_spline_init(output,x,y,n);
	return output;
}


#define NZ_FORMAT_HISTOGRAM 1
#define NZ_FORMAT_SAMPLE 2

void nz_string_error_message(){
	fprintf(stderr, "\nIf specifying n(z) format when loading, please use one of:\n");
	fprintf(stderr, "format=histogram\n");
	fprintf(stderr, "format=sample\n\n");	
}

// n(z) values come in two general forms.
// Most n(z) codes generate histograms in small z bins, with zlow, zhigh.
// This isn't a great model for n(z) of course - we don't really think that
// the number density is really blocky like that for real.  So we also support
// another mode where we take the z as sample values through which a spline should
// be drawn.
int detect_nz_format(c_datablock * block, const char * section)
{
	char *nz_format_string=NULL;
	int nz_format;
	int status = 0;
	bool has_format = c_datablock_has_value(block, section, "format");
	if (has_format){
		status |= c_datablock_get_string(block, section, "format", &nz_format_string);
		if (nz_format_string==NULL){
			nz_string_error_message();
			status = 1;
		}
		if (!strcmp(nz_format_string, "histogram")){
			nz_format = NZ_FORMAT_HISTOGRAM;
		}
		else if (!strcmp(nz_format_string, "sample")){
			nz_format = NZ_FORMAT_SAMPLE;
		}
		else{
			nz_string_error_message();
			status = 1;
		}

		free(nz_format_string);
		if (status){
			return -1;
		}
	}
	else{
		nz_format = NZ_FORMAT_SAMPLE;	
	}
	return nz_format;
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

int get_nchi_array(c_datablock * block, const char * nz_section,
	 int bin, double * z, double z_max, gsl_spline * a_of_chi, gsl_spline * chi_of_z, 
	  double * chi_out, int num_chi_out, double * nofchi_out)
{
	gsl_spline * nofchi_spline;
	int status = 0;
	double chi_max = gsl_spline_eval(chi_of_z, z_max, NULL);

	nofchi_spline = get_named_nchi_spline(block, nz_section, bin, z, a_of_chi, chi_of_z);
	//Now interpolate to output chi - if chi_out>chi, set to zero
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	for (int i_chi = 0; i_chi<num_chi_out; i_chi++){
		if (chi_out[i_chi] <= chi_max){
			nofchi_out[i_chi] = gsl_spline_eval(nofchi_spline, chi_out[i_chi], acc);
		}
		else {
			nofchi_out[i_chi] = 0.;
		}
	}
	printf("got nofchi_out\n");

	gsl_spline_free(nofchi_spline);
	return status;
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

	// Do the main computation
	gsl_spline * W = shear_shear_kernel(chi_max, n_of_z_spline,
		a_of_chi_spline);

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


/*
	Save an angular power spectrum to the datablock (for a single pair of z bins).
*/
static
int save_c_ell(c_datablock * block, const char * section,
	int bin1, int bin2, double coeff, gsl_spline * s, limber_config * lc)
{
	char name[64];
	snprintf(name, 64, "bin_%d_%d",bin1,bin2);
	int n_ell = lc->n_ell;
	double c_ell[n_ell];
	for (int i=0; i<n_ell; i++){
		double ell = lc->ell[i];
		if (lc->xlog) ell = log(ell);
		c_ell[i] = gsl_spline_eval(s, ell, NULL);
		if (lc->ylog && (lc->status==LIMBER_STATUS_OK)){
			c_ell[i] = exp(c_ell[i]);
		}
		c_ell[i]*=coeff;
		
	}

	int status = c_datablock_put_double_array_1d(block, section, name, c_ell, n_ell);
	status |= c_datablock_put_metadata(block, section, name, "note", "no ell^2 prefactor");
	return status;

}

/*
	Decide the name of the datablock section for saving the angular power spectrum.
*/
const char * choose_output_section(spectrum_type_t spectrum_type,
	shear_spectrum_config * config)
{
	switch(spectrum_type){
		case matter_matter: // config->matter_spectra
			return "matter_cl";
			break;
		case shear_intrinsic: // config->intrinsic_alignments
			return 	SHEAR_CL_GI_SECTION;
			break;
		case intrinsic_intrinsic: // config->intrinsic_alignments
			return SHEAR_CL_II_SECTION;
			break;
		case shear_shear: // config->shear_shear
		        //return (config->intrinsic_alignments ? SHEAR_CL_GG_SECTION : SHEAR_CL_SECTION);
		        return SHEAR_CL_SECTION;
			break;
		case matter_shear: // config->ggl_spectra
			return "ggl_cl";
			break;
		case matter_intrinsic: // config->gal_IA_cross_spectra
			return "gal_IA_cross_cl";
			break;
		case matter_magnification: // config->galaxy_magnification
			return "galaxy_magnification_cl";
			break;
		case magnification_magnification: // config->magnification_magnification
			return "magnification_magnification_cl";
			break;
		case magnification_intrinsic: // config->magnification_intrinsic
			return "magnification_intrinsic_cl";
			break;
		case magnification_shear: // config->magnification_shear
			return "magnification_shear_cl";
			break;
		default:
			return NULL;
	}
}

/*
	Choose the range of support and the amplitude scaling for the angular power spectrum model
*/
int choose_configuration(c_datablock * block, spectrum_type_t spectrum_type,
	limber_config * lc, shear_spectrum_config * config)
{
	int status = 0;
	// First we need to decide whether to use logs in the splines.
	// This is more accurate, but we cannot use it for cross-spectra
	// since they can go negative.
  	// Also, can't use for intrinsic_intrinsic if A=0 (i.e. zero signal)
	if (spectrum_type==intrinsic_intrinsic) {
	        double A;
	        status = c_datablock_get_double(block, "intrinsic_alignment_parameters", "A", &A);
		double eps=1e-9;
		if ( A>-1*eps && A<eps){
		  lc->xlog=false;
		  lc->ylog=false;
		}
		else {
		  lc->xlog=true;
		  lc->ylog=true;
		}
	}
	else if (spectrum_type==shear_shear || spectrum_type==matter_matter ||
		       spectrum_type==magnification_magnification){
		lc->xlog=true;
		lc->ylog=true;

	}
	else{
		lc->xlog=false;
		lc->ylog=false;
	}

	// We just fix these for now.
	lc->relative_tolerance = 1e-3;
	lc->absolute_tolerance = 1e-5;

	// Now we need to set up the ell range we wish to compute.
	// We need the number of ell values and the ell min and max.
	lc->n_ell = config->n_ell;

	// We log space the actual ell values
	lc->ell = malloc(sizeof(double)*lc->n_ell);
	double alpha = (log(config->ell_max) - log(config->ell_min)) / (config->n_ell - 1);
	for (int i=0; i<lc->n_ell; i++) lc->ell[i] = config->ell_min*exp(alpha * i);

	//Scalings

	// This scaling value is the bit that goes in front
	// of the kernels. For n(chi) this is 1.0 and for
	// shear W(chi) it's this more complex number.

	// shear-shear picks up two copies of this scaling,
	// shear-intrinsic picks up one,
	// and the other spectra just have a unity prefactor.

	// Note that this scaling is redshift independent.
	// The magnification spectra also have a luminosity
	// term for each bin pairing, but this is dealt with
	// separately.  

	double omega_m;
	status = c_datablock_get_double(block, cosmo, "omega_m", &omega_m);
	// in (Mpc/h)^-2
	// We leave the factor of h in because our P(k)
	// should have it in.  That's why there is a 100 in here.
	const double c_kms = 299792.4580; //km/s
	double shear_scaling = 1.5 * (100.0*100.0)/(c_kms*c_kms) * omega_m;

	switch(spectrum_type){
		case shear_shear:
			lc->prefactor = shear_scaling*shear_scaling;
			break;
		case shear_intrinsic:
			lc->prefactor = shear_scaling;
			break;
		case intrinsic_intrinsic:
			lc->prefactor = 1.0;
			break;
		case matter_matter:
			lc->prefactor = 1.0;
			break;
		case matter_shear:
			lc->prefactor = shear_scaling;
			break;
		case matter_intrinsic:
			lc->prefactor = 1.0;
			break;
		case matter_magnification:
			lc->prefactor = shear_scaling;
			break;
		case magnification_magnification:
			lc->prefactor = shear_scaling*shear_scaling;
			break;
		case magnification_intrinsic:
			lc->prefactor = shear_scaling;
			break;
		case magnification_shear:
			lc->prefactor = shear_scaling*shear_scaling;
			break;
		default:
			lc->prefactor = 1.0/0.0;
			return 1;
	}
	return status;
}

/*
	Select the combinations of line-of-sight kernels for the Limber integral.
*/
int choose_kernels(spectrum_type_t spectrum_type, int nbin, gsl_spline * W[nbin],
	gsl_spline * Nchi[nbin], gsl_spline ** K1, gsl_spline ** K2)
{
	// The different spectra use different kernels to integrate over:
	// Here we choose which to use (we do it like this so we can use
	// the same code for the different spectra)

	switch(spectrum_type){
		case shear_shear:
		case magnification_magnification:
		case magnification_shear:
			for(int b=0; b<nbin; b++) {K1[b] = W[b]; K2[b] = W[b];};
			break;
		case intrinsic_intrinsic:
		case matter_matter:
		case matter_intrinsic:
			for(int b=0; b<nbin; b++) {K1[b] = Nchi[b]; K2[b] = Nchi[b];};
			break;
		case shear_intrinsic:
		// This is actually IG rather than GI
		case matter_shear:
		case matter_magnification:
			for(int b=0; b<nbin; b++) {K1[b] = Nchi[b]; K2[b] = W[b];};
			break;
		case magnification_intrinsic:
			for(int b=0; b<nbin; b++) {K1[b] = W[b]; K2[b] = Nchi[b];};
			break;
		default:
			return 1;
	}
	return 0;

}

int choose_bin2_max(spectrum_type_t spectrum_type, int nbin, int bin1)
{
	if (spectrum_type==shear_shear || spectrum_type==matter_matter ||
		spectrum_type==intrinsic_intrinsic || spectrum_type==magnification_magnification){
		return bin1;
	}
	else{
		return nbin;
	}

}

/*
	Choose the appropriate prefactors for the magnification spectra
*/
double choose_limber_coefficient(spectrum_type_t spectrum_type, double alpha_1, double alpha_2)
{
	double coeff;

	switch(spectrum_type){
		case magnification_magnification:
			coeff = 4.0*(alpha_1-1.0)*(alpha_2-1.0);
			break;
		case magnification_intrinsic:
		case magnification_shear:
			coeff = 2.0*(alpha_1-1.0);
			break;
		case matter_magnification:
			coeff = 2.0*(alpha_2-1.0);
			break;
		case shear_shear:
		case intrinsic_intrinsic:
		case matter_matter:
		case matter_intrinsic:
		case shear_intrinsic:
		default:
			coeff = 1.0/0.0;
	}

	return coeff;
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

	// In addition to the Limber kernel, the magnification 
	// spectra also require the slope of the luminosity 
	// function at the limiting magnitude of the survey
	double *alpha, coeff;
	int Nzbin;
	int is_mag =   spectrum_type==magnification_magnification 
	            || spectrum_type==matter_magnification 
	            || spectrum_type==magnification_shear
	            || spectrum_type==magnification_intrinsic;

	if (is_mag){
		status |= c_datablock_get_double_array_1d(block, lum, "alpha_binned", &alpha, &Nzbin);
	}
	if (status) return status;

	for (int bin1=1; bin1<=nbin; bin1++){
		// We also need to choose the max value of bin 2.
		// For auto-spectra we don't compute both bin_i_j and bin_j_i,
		// whereas for cross-spectra we do
		int bin2_max = choose_bin2_max(spectrum_type, nbin, bin1);
		for (int bin2=1; bin2<=bin2_max; bin2++){
			// reset the status
			lc.status = LIMBER_STATUS_OK;
			gsl_spline * c_ell = limber_integral(&lc, K1[bin1-1], K2[bin2-1], PK);
			if (c_ell == NULL)  return 1;
			if (is_mag) coeff = choose_limber_coefficient(spectrum_type, alpha[bin1-1], alpha[bin2-1]);
			else coeff=1;
			int status = save_c_ell(block, section, bin1, bin2, coeff, c_ell, &lc);
			gsl_spline_free(c_ell);
			if (status) return status;
		}
	}
	free(lc.ell);
	if (is_mag) free(alpha);
	return status;
}

/*
	If the relevant shear spectra have already been calculated the magnification spectra do not need to be calculated a priori. 
	These other spectra can simply be rescaled using the logarithmic slope of the luminosity function.
	See e.g. Kirk, Rassat, Host and Bridle (2011)
*/
int rescale_spectrum(c_datablock * block, spectrum_type_t spectrum_type, int nbin, shear_spectrum_config * config)
{
	int status = 0;

	// Choose the appropriate sections to read from and save to
	const char *in_section;
	const char * out_section = choose_output_section(spectrum_type, config);

	switch(spectrum_type){
		case magnification_magnification:
			in_section = choose_output_section(shear_shear, config);
			break;
		case matter_magnification:
			in_section = choose_output_section(matter_shear, config);
			break;
		case magnification_intrinsic:
			in_section = choose_output_section(shear_intrinsic, config);
			break;
		case magnification_shear:
			in_section = choose_output_section(shear_shear, config);
			break;
                default:
                        return -1;
                
	}

	// The l sampling can be copied directly across to the new section
	double *l = NULL;
	int nl;
	status|=c_datablock_get_double_array_1d(block, in_section, "ell", &l, &nl);
	status|=c_datablock_put_double_array_1d(block, out_section, "ell", l, nl);
	if (status){
		fprintf(stderr, "Could not find ell value in %s\n", in_section);
		return status;
	}

	double *alpha = NULL;
	status|=c_datablock_get_double_array_1d(block, lum, "alpha_binned", &alpha, &nbin);
	if (status){
		fprintf(stderr, "Could not find alpha value in %s\n", lum);
		return status;
	}


	char spectrum_name[64];
	double coeff, *C_l_new, *C_l_old;
	C_l_new = (double*)malloc(sizeof(double)*nl);

	int bin2_max;
	
	// Loop first over the redshift bins
	// Load the relevant spectrum for each pair of tomographic bins
	// Then cycle through the angular frequency bins, reweighting by the luminosity function
	for (int bin1=1 ; bin1< nbin+1; ++bin1){
		bin2_max = choose_bin2_max(spectrum_type, nbin, bin1);
		for (int bin2=1 ; bin2< bin2_max+1 ; ++bin2){
			// This line deals with the case of the mG spectrum, which isn't symmetric but is 
			// obtained by rescaling the GG spectrum, which is (and of which we only have half the i,j pairs) 
			if (spectrum_type==magnification_shear && bin2>bin1){ snprintf(spectrum_name, 64, "bin_%d_%d",bin2,bin1); }
			else { snprintf(spectrum_name, 64, "bin_%d_%d", bin1, bin2); }

			status|=c_datablock_get_double_array_1d(block, in_section, spectrum_name, &C_l_old, &nl);
		
			coeff = choose_limber_coefficient(spectrum_type, alpha[bin1-1], alpha[bin2-1]);

			for (int bin3=0 ; bin3<nl ; ++bin3){ 
				C_l_new[bin3] = coeff*C_l_old[bin3];
			}

			snprintf(spectrum_name, 64, "bin_%d_%d",bin1,bin2);
			status |= c_datablock_put_double_array_1d(block, out_section, spectrum_name, C_l_new, nl);
		}	
	}

	free(C_l_old);
	free(C_l_new);
	free(l);
	free(alpha);

	return status;

}

#define printf_verbose(...) if (config->verbose) printf(__VA_ARGS__)


int execute(c_datablock * block, void * config_in)
{
	DATABLOCK_STATUS status=0;
	double * chi;
	double * a;
	double * z;
	int nz1, nz2;
	int nbin;

	shear_spectrum_config * config = (shear_spectrum_config*) config_in;

	// Load the number of bins we will use
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

	// ------------
	// Get the P(k)s we need
	printf_verbose("Loading interpolation splines\n");
	/*				Matter power spectrum			 		*/

	Interpolator2D * PK = load_interpolator_chi(
		block, chi_of_z_spline, MATTER_POWER_NL_SECTION, "k_h", "z", "P_k");

	/*				Galaxy power spectrum			 		*/

	// Get a galaxy power spectrum with a bias model that distinguishes from the
	// matter power spectrum if present in the data block. Otherwise, just use the
	// matter power spectrum for the 'galaxy power spectrum'

	// It seems useful to highlight this decision. If the code does this
	// silently, there's potential for confusion if an earlier module is omitted.
	Interpolator2D * PK_gg = NULL;
	int free_pk_gg;
	if (config->matter_spectra) {
		if (c_datablock_has_section(block, "matter_power_gal")) {
			PK_gg = load_interpolator_chi(block, chi_of_z_spline, "matter_power_gal", "k_h", "z", "P_k");
			free_pk_gg = 1;
		}
		if (PK_gg==NULL){
			printf_verbose("No galaxy power spectrum found in datablock. Using the matter power spectrum instead.\n"); 
			PK_gg=PK;
			free_pk_gg = 0;
		}
	}

	if (PK==NULL) return 1;
	if (!PK) {
		free(chi);
		free(a);
		free(z);
		gsl_spline_free(a_of_chi_spline);
		gsl_spline_free(chi_of_z_spline);
		return 1;
	}

	/*				Intrinsic alignment spectra			                   */

	Interpolator2D * PK_GI = NULL;
	Interpolator2D * PK_II = NULL;
	if (config->intrinsic_alignments){
		PK_II = load_interpolator_chi(
			block, chi_of_z_spline, "intrinsic_power", "k_h", "z", "p_k");

		PK_GI = load_interpolator_chi(
			block, chi_of_z_spline, "matter_intrinsic_power", "k_h", "z", "p_k");

		if (PK_II==NULL) return 2;
		if (PK_GI==NULL) return 3;

	}

	/*				galaxy position - mass/shear spectrum			           */
	// Get a galaxy-mass cross-power spectrum with a bias model that distinguishes
	// from the matter power spectrum if present in the data block. Otherwise,
	// just use the matter power spectrum for the 'galaxy-mass power spectrum'
	Interpolator2D * PK_gm = NULL;
	int free_pk_gm;
	if (config->ggl_spectra | config->galaxy_magnification) {
		if (c_datablock_has_section(block, "matter_power_gal_mass")) {
			PK_gm = load_interpolator_chi(block, chi_of_z_spline, "matter_power_gal_mass", "k_h", "z", "P_k");
			free_pk_gm = 1;
		}
		if (PK_gm==NULL){
			printf_verbose("No galaxy power spectrum found in datablock. Using the matter power spectrum instead.\n");
			PK_gm=PK;
			free_pk_gm = 0;
		}
	}

	/*				galaxy position - intrinsic shape spectrum			 		*/
	Interpolator2D * PK_gI = NULL;
	int free_pk_gi;
	if (config->gal_IA_cross_spectra) {
		if (c_datablock_has_section(block, "matter_power_gal_intrinsic")) {
			PK_gI = load_interpolator_chi(block, chi_of_z_spline, "matter_power_gal_intrinsic", "k_h", "z", "P_k");
			free_pk_gi = 1;
		}
		if (PK_gI==NULL){
			printf_verbose("No galaxy position-intrinsic power spectrum found in datablock. Using the II spectrum instead.\n");
			PK_gI=PK_II;
			free_pk_gi = 0;
		}
	}

	// ------------
	// Compute the angular power spectra
	if (config->shear_shear){
		status |= compute_spectra(block, nbin, shear_shear,
			W_splines, Nchi_splines, PK, config);
		printf_verbose("Saved shear-shear spectrum.\n");
	}
	if (config->intrinsic_alignments){

		status |= compute_spectra(block, nbin, intrinsic_intrinsic,
			W_splines, Nchi_splines, PK_II, config);

		status |= compute_spectra(block, nbin, shear_intrinsic,
			W_splines, Nchi_splines, PK_GI, config);

		printf_verbose("Saved IA spectra.\n");

	}
	if (config->matter_spectra){
		status |= compute_spectra(block, nbin, matter_matter,
			W_splines, Nchi_splines, PK_gg, config);
		printf_verbose("Saved galaxy-galaxy spectrum.\n");
		
	}
	if (config->ggl_spectra){
		status |= compute_spectra(block, nbin, matter_shear,
			W_splines, Nchi_splines, PK_gm, config);
		printf_verbose("Saved galaxy-shear spectrum.\n");
		
	}
	if (config->gal_IA_cross_spectra){
		status |= compute_spectra(block, nbin, matter_intrinsic,
			W_splines, Nchi_splines, PK_gI, config);
		printf_verbose("Saved galaxy-intrinsic shape spectrum.\n"); 
		if (status) return status;
	}
	// For the spectra involving magnification, check the spectra already saved to the datablock
	// It makes things quicker if the shear spectra can be rescaled rather than evaluating the
	// full Limber integral 
	if (config->galaxy_magnification){
		if (c_datablock_has_section(block, "ggl_cl")) {
			status |= rescale_spectrum(block, matter_magnification, nbin, config);
		}
		else { 
			status |= compute_spectra(block, nbin, matter_magnification,
				W_splines, Nchi_splines, PK, config);
		}
		if (status) return status;
		printf_verbose("Saved galaxy-magnification spectrum.\n");
	}

	if (config->magnification_magnification) {
		if ( c_datablock_has_section(block, SHEAR_CL_GG_SECTION) ||   
		     c_datablock_has_section(block, SHEAR_CL_SECTION) ) {
			status |= rescale_spectrum(block, magnification_magnification, nbin, config);
		}
		else { 
			status |= compute_spectra(block, nbin, magnification_magnification,
				W_splines, Nchi_splines, PK, config);
		}
		if (status) return status;
		printf_verbose("Saved magnification-magnification spectrum.\n");
	}

	if (config->magnification_intrinsic) {
		if (c_datablock_has_section(block, SHEAR_CL_GI_SECTION)) {
			status |= rescale_spectrum(block, magnification_intrinsic, nbin, config);
		}
		else { 
			status |= compute_spectra(block, nbin, matter_intrinsic,
				W_splines, Nchi_splines, PK, config);
		}
		if (status) return status;
		printf_verbose("Saved magnification-intrinsic shape spectrum.\n");
	}
	
	if (config->magnification_shear ) {
		if ( c_datablock_has_section(block, SHEAR_CL_GG_SECTION) ||
		     c_datablock_has_section(block, SHEAR_CL_SECTION) ) {
			status |= rescale_spectrum(block, magnification_shear, nbin, config);
		}
		else { 
			status |= compute_spectra(block, nbin, magnification_shear,
				W_splines, Nchi_splines, PK, config);
		}
		if (status) return status;
		printf_verbose("Saved shear-shear spectrum.\n");
	}
	
	// tidy up global data
	for (int bin=0; bin<nbin; bin++) gsl_spline_free(W_splines[bin]);
	gsl_spline_free(a_of_chi_spline);
	gsl_spline_free(chi_of_z_spline);

	destroy_and_null(&PK);
	destroy_and_null(&PK_GI);
	destroy_and_null(&PK_II);
	if (free_pk_gm) destroy_and_null(&PK_gm);
	if (free_pk_gg) destroy_and_null(&PK_gg);
	if (free_pk_gi) destroy_and_null(&PK_gI);
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
