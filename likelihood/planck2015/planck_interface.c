#include "clik.h"
#include "cosmosis/datablock/c_datablock.h"
#include <math.h>

// Hard-coded maximum number of likelihoods.
// This is far more than have been released.
// I'm just doing this to avoid some mallocs and stuff.
#define MAX_NUMBER_LIKES 128

typedef struct configuration_data{
	int ndata;
	clik_object * clik_data[MAX_NUMBER_LIKES];


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

	int status = 0;
	config->ndata=0;
	char param[512];
	while (1){
		snprintf(param, 128, "data_%d", config->ndata+1);
		char * value;
		status |= c_datablock_get_string_default(options, OPTION_SECTION, 
			param, "", &value);

		if (strlen(value)>0){
			printf("Looking for clik Planck likelihood file %d: %s\n", 
				config->ndata+1, value);
			config->clik_data[config->ndata] = clik_init(value, &err);
			config->ndata++;


			if (isError(err)){
				fprintf(stderr,"There was an error initializating one of the Planck likelihoods.\n");
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
			else{
			printf("%s successfully initialized.\n", value);

			}

		}		
		else{
			break;
		}
	}



	if (config->ndata==0){
		fprintf(stderr, "No data files were specified for Planck at all!\n");
		fprintf(stderr, "In the options file (not the values.ini file) you need to set in the section for Planck at least the first of:\n");
		fprintf(stderr, "   data_1 = ...\n");
		fprintf(stderr, "   data_2 = ...\n");
		fprintf(stderr, "   data_3 = ...\n");
		fprintf(stderr, "   data_4 = ...\n");
		fprintf(stderr, "   etc.\n");
		fprintf(stderr, "The Planck data files can be downloaded from:\n");
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
				package_names[i], nread+1, lmax_by_type[i], status);
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
	const int N_CL = 6;
	int param_count=0;
	int status = 0;
	char * package_names[] = {"TT","EE","BB","TE","TB","EB"};
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


int execute(c_datablock * block, configuration_data * config){

	double like = 0.0;
	int status = 0;

	for (int i=0; i<config->ndata; i++){
		// Compute the likelihood for this file.
		double like_i = 0.0;
		status = run_clik_cosmosis(block, config->clik_data[i], &like_i);

		//
		char name[64];
		snprintf(name, 64, "PLANCK_%d_LIKE", i+1);
		status |= c_datablock_put_double(block, LIKELIHOODS_SECTION, name, like_i);

		like += like_i;
	}

	if (status) return status;

	// Save the total likelihood
	status |= c_datablock_put_double(block, LIKELIHOODS_SECTION, "PLANCK2015_LIKE", like);

	return 0;
}
