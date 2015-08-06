#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include "utils.h"
#include <stdio.h>
#include <math.h>

const char * nl = MATTER_POWER_NL_SECTION;
const char * mps = MATTER_POWER_LIN_SECTION;

// Define some structs to hold the configuration parameters and biases

typedef enum spectrum_type_t {
	// Intrinsic aligments
	mass_intrinsic,
	intrinsic_intrinsic,

	// Galaxy clustering
	galaxy_galaxy,
	galaxy_mass,

	// Cross spectrum
	galaxy_intrinsic
	
} spectrum_type_t;

typedef struct bias_config {
	int verbosity;
	bool intrinsic_alignments;
	bool galaxy_bias;
	
} bias_config;

typedef struct biases {
	double ** b_I;
	double ** r_I;
	double ** b_g;
	double ** r_g;
	double ** P_nl;
	double *k_h;
	int nk;
	double *z;
	int nz;
} biases;

// Define functions called by the interface

const char * choose_output_section(spectrum_type_t spectrum_type){
	switch(spectrum_type){
		case intrinsic_intrinsic:
			return IA_SPECTRUM_II_SECTION;
			break;
		case mass_intrinsic:
			return IA_SPECTRUM_GI_SECTION;
			break;
		case galaxy_galaxy:
			return MATTER_POWER_GAL_SECTION;
			break;
		case galaxy_mass:
			return MATTER_POWER_GAL_MASS_SECTION;
			break;
		case galaxy_intrinsic:
			return MATTER_POWER_GAL_INTRINSIC_SECTION;
			break;
		default:
			return NULL;
	}
}

const char * choose_spectrum_name(spectrum_type_t spectrum_type){
	switch(spectrum_type){
		case intrinsic_intrinsic:
			return "P_II";
			break;
		case mass_intrinsic:
			return "P_GI";
			break;
		case galaxy_galaxy:
			return "P_k";
			break;
		case galaxy_mass:
			return "P_k";
		case galaxy_intrinsic:
			return "P_k";
			break;
		default:
			return NULL;
	}
}

double apply_bias_coefficients(spectrum_type_t spectrum_type, double b_I, double r_I, double b_g, double r_g, double P_nl){

	double Pk; 

	switch(spectrum_type){
		case intrinsic_intrinsic:
			Pk = b_I * b_I * P_nl;
			break;
		case galaxy_galaxy:
			Pk = b_g * b_g * P_nl; 
			break;
		case mass_intrinsic:
			Pk = r_I * b_I * P_nl; 
			break;
		case galaxy_mass:
			Pk = r_g * b_g * P_nl; 
			break;
		case galaxy_intrinsic:
			Pk = r_I * b_I * r_g * b_g * P_nl; 
			break;
	}

	return Pk; 

}

int get_power_spectrum(c_datablock * block, spectrum_type_t spectrum_type, biases * b){
	int status=0;


	const char * section = choose_output_section(spectrum_type);
	const char * spectrum_name = choose_spectrum_name(spectrum_type);

	double **Pk;
	Pk= (double **)malloc(b->nk*sizeof(double*));
		Pk[0] = (double *)malloc(b->nk*b->nz*sizeof(double));
		for(int i= 1 ; i< b->nk ; i++){
			Pk[i]= Pk[0]+i*b->nz;}

	// Multiply by the coefficients needed for the current spectrum
	// This gives a measurable spectrum, which represents the underlying mass distribution in a biased way
	for (int i=0 ; i<b->nz ; ++i){
		for (int j=0 ; j<b->nk ; ++j){
			Pk[j][i] = apply_bias_coefficients(spectrum_type, b->b_I[j][i], b->r_I[j][i], b->b_g[j][i], b->r_g[j][i], b->P_nl[j][i]);
		}
	}	

	status|=c_datablock_put_double_grid(block, section, "k_h", b->nk, b->k_h, "z", b->nz, b->z, spectrum_name, Pk);

	free(Pk);

	return status ;
}

int load_biases(c_datablock * block, biases * b, bias_config * config ){
	int status=0;

	// First load the nonlinear spectrum 
	status |= c_datablock_get_double_grid(block, MATTER_POWER_NL_SECTION, "k_h", &(b->nk), &(b->k_h), "z", &(b->nz), &(b->z), "P_k", &(b->P_nl));
	if (config->verbosity>1) printf("Loaded matter power spectrum.\n");

	// Then the biases
	int nk,nz;
	if (config->intrinsic_alignments){
		status |= c_datablock_get_double_grid(block, "intrinsic_alignment_parameters", "k_h", &(b->nk), &(b->k_h), "z", &(b->nz), &(b->z), "b_I", &(b->b_I) );
		status |= c_datablock_get_double_grid(block, "intrinsic_alignment_parameters", "k_h", &(b->nk), &(b->k_h), "z", &(b->nz), &(b->z), "r_I", &(b->r_I) );
		if(config->verbosity > 1){ printf("Loaded intrinsic alignment bias.\n"); }
	} 

	if (config->galaxy_bias){
		status |= c_datablock_get_double_grid(block, "bias_field", "k_h", &(b->nk), &(b->k_h), "z", &(b->nz), &(b->z), "b_g", &(b->b_g) );
		status |= c_datablock_get_double_grid(block, "bias_field", "k_h", &(b->nk), &(b->k_h), "z", &(b->nz), &(b->z), "r_g", &(b->r_g) );
		if(config->verbosity > 1){ printf("Loaded galaxy bias.\n"); }
	}

return status;

} 

int get_all_spectra(c_datablock * block, biases * b, bias_config * config){
	int status = 0;

	// Apply the bias arrays loaded earlier to the nonlinear power spectrum to get the required spectra
	if(config->intrinsic_alignments){
		status|=get_power_spectrum(block, intrinsic_intrinsic, b);
		if(config->verbosity>0){ printf("Saved II spectrum.\n");}
		status|=get_power_spectrum(block,  mass_intrinsic, b);	
		if(config->verbosity>0){ printf("Saved GI spectrum.\n");}	
	}
	if(config->galaxy_bias){
		status|=get_power_spectrum(block, galaxy_galaxy, b);
		if(config->verbosity>0){ printf("Saved gg spectrum.\n");}
	}
	if (config->intrinsic_alignments && config->galaxy_bias){
		status|=get_power_spectrum(block, galaxy_mass ,b);
		if(config->verbosity>0){ printf("Saved gG spectrum.\n");}
		status|=get_power_spectrum(block, galaxy_intrinsic ,b);
		if(config->verbosity>0){ printf("Saved gI spectrum.\n");}
	}

	return status;
}
