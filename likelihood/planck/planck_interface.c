#include "clik.h"
#include "cosmosis/datablock/c_datablock.h"
#include <math.h>

typedef struct configuration_data{
	int ready;
	clik_object * T_low_data;
	clik_object * P_low_data;
	clik_object * T_high_data;
	clik_object * lensing_data;


} configuration_data;


static int find_max_lmax(clik_object * like_obj)
{
	error * err = initError();
	int n_cl = 6;
	int lmax_by_type[n_cl];
	clik_get_lmax(like_obj, lmax_by_type,&err);
	int max_lmax = 0;
	int i;
	for (i=0; i<n_cl; i++) {
		if (lmax_by_type[i]>max_lmax){
			max_lmax = lmax_by_type[i];
		}
	}
	return max_lmax;
	endError(&err);
}



configuration_data * setup(c_datablock * options){
	error * err = initError();

	configuration_data * config = malloc(sizeof(configuration_data));
	char * T_low_file=""; 
	char * P_low_file="";
	char * T_high_file="";
	char * lensing_file="";

	int status = 0;

	status |= c_datablock_get_string_default(options, OPTION_SECTION, "t_low_file", "", &T_low_file);
	status |= c_datablock_get_string_default(options, OPTION_SECTION, "p_low_file", "", &P_low_file);
	status |= c_datablock_get_string_default(options, OPTION_SECTION, "t_high_file", "", &T_high_file);
	status |= c_datablock_get_string_default(options, OPTION_SECTION, "lensing_file", "", &lensing_file);


	printf("Likelihood files to be used by Planck:\n");
	printf("t_low_file = %s\np_low_file = %s\nt_high_file = %s\nlensing_file = %s\n\n", T_low_file, P_low_file, T_high_file, lensing_file);
	config->T_low_data = NULL;
	config->P_low_data = NULL;
	config->T_high_data = NULL;
	config->lensing_data = NULL;

	if (strlen(T_low_file))   config->T_low_data   = clik_init(T_low_file, &err);
	if (strlen(P_low_file))   config->P_low_data   = clik_init(P_low_file, &err);
	if (strlen(T_high_file))  config->T_high_data  = clik_init(T_high_file,&err);
	if (strlen(lensing_file)) config->lensing_data = clik_lensing_init(lensing_file,&err);

	if (status){
		fprintf(stderr,"There was an getting parameters to initialize the Planck likelihood.\n");
		fprintf(stderr,"The error code was %d\n\n", status);
		exit(1);
	}

	if (isError(err)){
		fprintf(stderr,"There was an error initializating the Planck likelihoods.\n");
		fprintf(stderr,"Here is the error message:\n");
		printError(stderr,err);
		fprintf(stderr,"\nThis probably means that a likelihood file was not found.\n\n");
		fprintf(stderr,"If you are running demo two or demo four then it just means that you\n");
		fprintf(stderr,"did not download the Planck likelihood when installing.\n");
		fprintf(stderr,"You can either get it manually from:\n");
		fprintf(stderr,"     http://pla.esac.esa.int/pla/aio/planckProducts.html\n");
		fprintf(stderr,"and edit demos/demo2.ini, or just skip demos two and four for now.\n\n");
		exit(1);
	}

	if (config->T_low_data==NULL && config->P_low_data==NULL && config->T_high_data==NULL && config->lensing_data==NULL){
		fprintf(stderr, "No data files were specified for Planck at all!\n");
		fprintf(stderr, "In the options file (not the values.ini file) you need to set in the section for Planck at least one of:\n");
		fprintf(stderr, "   t_low_file\n");
		fprintf(stderr, "   t_high_file\n");
		fprintf(stderr, "   p_low_file\n");
		fprintf(stderr, "   lensing_file\n");
		fprintf(stderr, "The Planck data files are available via the cosmosis boostrap script with the -d option,\n");
		fprintf(stderr, "or by downloading them manually at the following url:\n");
		fprintf(stderr, "http://pla.esac.esa.int/pla/aio/planckProducts.html\n");
		exit(1);
	}


	endError(&err);

	return config;
}


static 
int extract_c_ell(c_datablock * block, int n_cl, 
	char * package_names[n_cl], int lmax_by_type[n_cl], int n_nuisance, 
	parname nuisance_names[n_nuisance], double *p, double *cl)
{
	int i;
	int status=0;

	for (i=0;i<n_cl;i++) {
		// Check lmax for this spectrum.  If not needed, continue
		int lmax = lmax_by_type[i];

		if (lmax<=0) continue;
		// Otherwise fill in first the trivial cl and then the rest, loaded fromt the package
		cl[0]=0.0;
		cl[1]=0.0;
		int nread=0;
		double * ell_squared_cl = NULL;
		status |= c_datablock_get_double_array_1d(block, CMB_CL_SECTION, package_names[i], 
			&ell_squared_cl, &nread);

		if (status){
			fprintf(stderr, "Could not get C_ell values for %s (status=%d)\n", 
				package_names[i],status);
			status = 1;
			return status;
		}
		if (nread<lmax-1){
			fprintf(stderr, "Could not get enough C_ell values for %s (read lmax=%d but needed lmax=%d).  May need to adjust lmax.\n", 
				package_names[i], nread+1, lmax_by_type[i]);
			status = 1;
			return status;
		}


		// Transform to raw cl from l**2 cl
		int ell;
		for (ell=2; ell<=lmax; ell++) {
			cl[ell] = ell_squared_cl[ell-2] / (ell*(ell+1.0)/2.0/M_PI);
		}
		free(ell_squared_cl);

		// Move the spectrum pointer forward ready for the next spectrum (e.g. TT->EE)
		cl+=lmax+1;
	}

	// ... and now get the nuisance parameters
	for (i=0; i<n_nuisance; i++){
		int this_status=c_datablock_get_double(block, PLANCK_SECTION, nuisance_names[i], 	cl+i  ); //ok
		// Report any missing parameters
		if (this_status) fprintf(stderr, "Could not get required Planck parameter %s in section %s\n",nuisance_names[i], PLANCK_SECTION);
		status |= this_status;
	}
	return status;
}


static 
int run_clik_cosmosis(c_datablock * block, clik_object * like_obj, double * output_like)
{
	int i;
#define N_CL 6
	int param_count=0;
	int status = 0;
	char * package_names[N_CL] = {"TT","EE","BB","TE","TB","EB"};
	int lmax_by_type[N_CL];

	// Planck internal error system
	error * err = initError();

	// Get the lmax values for the different spectra and sum them
	clik_get_lmax(like_obj, lmax_by_type,&err);
	for (i=0;i<N_CL;i++) if (lmax_by_type[i]>0) param_count += lmax_by_type[i]+1; 

	// Check for any errors from the Planck library
	if (isError(err)){
		fprintf(stderr,"There was an error getting lmax from a Planck likelihood.\n");
		endError(&err);
		return 1;
	}

	// Get the number and names of extra parameters and add to total
	int n_nuisance;
	parname * nuisance_names;
	n_nuisance = clik_get_extra_parameter_names(like_obj, &nuisance_names, &err);
	param_count += n_nuisance;

	// Check for any errors from the Planck library
	if (isError(err)){
		fprintf(stderr,"There was an error getting nuisance params from a Planck likelihood.\n");
		endError(&err);
		return 2;
	}

	// p is the overall space and cl is the space within it for each spectra.
	double * p = malloc(sizeof(double)*param_count);
	double * cl = p;

	// Get the c_ell we need, and the nuisance parameters (if any)
	status |= extract_c_ell(block, N_CL, package_names, lmax_by_type, n_nuisance, 
		nuisance_names, p,cl);

	free(nuisance_names);
	// Check for any errors reading things in
	if (status){
		free(p);
		return status;
	}

	// Compute the actual likelihood
	double like = clik_compute(like_obj, p, &err);
	free(p);
	// ...  and check for errors
	if (isError(err)){
		status = 1;
		fprintf(stderr,"Planck likelihood error\n");
	}
	// Clean up Planck error system and 
	// return required outputs
	endError(&err);	
	*output_like = like;
	return 0;
}
#undef N_CL

static 
int run_clik_cosmosis_lensing(c_datablock * block, clik_object * like_obj, double * output_like)
{
	int i;
#define N_CL 2
	int param_count=0;
	int status = 0;
	char * package_names[N_CL] = {"PP","TT"};
	error * err = initError();

	// Get the lmax values for the different spectra and sum them
	int lmax = clik_lensing_get_lmax(like_obj,&err);
	int lmax_by_type[N_CL] = {lmax, lmax};
	param_count = N_CL*(lmax+1);

	// Check for any errors from the Planck library
	if (isError(err)){
		fprintf(stderr,"There was an error getting lmax from a Planck likelihood.\n");
		return 1;
	}

	// Get the number and names of extra parameters and add to total
	int n_nuisance;
	parname * nuisance_names;
	n_nuisance = clik_lensing_get_extra_parameter_names(like_obj, &nuisance_names, &err);
	param_count += n_nuisance;

	// Check for any errors from the Planck library
	if (isError(err)){
		fprintf(stderr,"There was an error getting nuisance params from a Planck likelihood.\n");
		return 2;
	}

	// p is the overall space and cl is the space within it for each spectra.
	double * p = malloc(sizeof(double)*param_count);
	double * cl = p;

	// Get the c_ell we need, and the nuisance parameters (if any)
	status |= extract_c_ell(block, N_CL, package_names, lmax_by_type, n_nuisance, 
		nuisance_names, p,cl);

	free(nuisance_names);

	// Check for any errors reading things in and quit early if problem
	if (status){
		free(p);
		return status;
	}


	// Compute the actual likelihood and check for errors
	double like = clik_lensing_compute(like_obj, p, &err);
	free(p);

	//Clean up
	if (isError(err)){
		status = 1;
		fprintf(stderr, "Error running Planck lensing\n");
	}

	// Clean up Planck error system and 
	// return required outputs
	endError(&err);	
	*output_like = like;
	return 0;

	return 0;
}
#undef N_CL

int execute(c_datablock * block, configuration_data * config){

	double tt_high_like = 0.0;
	double tt_low_like = 0.0;
	double p_low_like = 0.0;
	double lensing_like = 0.0;
	int status = 0;

	// Run the individual clik likelihoods
	if (config->T_high_data){
		status = run_clik_cosmosis(block, config->T_high_data, &tt_high_like);
		status |= c_datablock_put_double(block, LIKELIHOODS_SECTION, "PLANCK_LIKE_TT_HIGH", tt_high_like);
	}
	if (config->T_low_data){
		status |= run_clik_cosmosis(block, config->T_low_data, &tt_low_like );
		status |= c_datablock_put_double(block, LIKELIHOODS_SECTION, "PLANCK_LIKE_TT_LOW", tt_low_like);
	}
	if (config->P_low_data){
		status |= run_clik_cosmosis(block, config->P_low_data, &p_low_like  );	
		status |= c_datablock_put_double(block, LIKELIHOODS_SECTION, "PLANCK_LIKE_P_LOW", p_low_like);
	} 
	// The lensing one is a little different
	if (config->lensing_data) {
		status |= run_clik_cosmosis_lensing(block, config->lensing_data, &lensing_like  );		
		status |= c_datablock_put_double(block, LIKELIHOODS_SECTION, "PLANCK_LIKE_LENSING", lensing_like);
	}

	if (status) return status;

	// Save the total likelihood
	double like = tt_high_like + tt_low_like + p_low_like + lensing_like;
	status |= c_datablock_put_double(block, LIKELIHOODS_SECTION, "PLANCK_LIKE", like);

	return 0;
}
