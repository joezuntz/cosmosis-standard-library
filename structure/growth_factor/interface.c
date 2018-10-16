#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include "growthfactor.h"


//Module to calculate the linear growth factor D, and linear growth rate, f. Where D, f are defined by the growth of a
//linear perturbation, delta, with scale factor a: 
//delta(a') = delta(a)*(D(a')/D(a)) and f = dlnD/dlna
//Anyone using Komatsu's CRL library should note: growth_factor_crl = D *(1+z) and growth_rate_crl = f/(1+z)


const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
const char * like = LIKELIHOODS_SECTION;
const char * growthparameters = GROWTH_PARAMETERS_SECTION;

typedef struct growth_config {
	double zmin;
	double zmax;
	double dz;
	int nz_lin;
	double zmax_log;
	int nz_log;
} growth_config;

void reverse(double * x, int n)
{
	for (int i=0; i<n/2; i++){
		double tmp = x[i];
		x[i] = x[n-1-i];
		x[n-1-i] = tmp;
	}
}

growth_config * setup(c_datablock * options)
{
	int status = 0;
	growth_config * config = malloc(sizeof(growth_config));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "zmin", 0.0, &(config->zmin));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "zmax", 3.0, &(config->zmax));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "dz", 0.01, &(config->dz));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "zmax_log", 1100.0, &(config->zmax_log));
	status |= c_datablock_get_int_default(options, OPTION_SECTION, "nz_log", 0, &(config->nz_log));

	config->nz_lin = (int)((config->zmax-config->zmin)/config->dz)+1;

	printf("Will calculate f(z) and d(z) in %d bins (%lf:%lf:%lf)\n", config->nz_lin, config->zmin, config->zmax, config->dz);
	if (config->nz_log > 0){
	  printf("Will extend growth calculation using %d logarithmically spaced z-bins to zmax = %lf \n", config->nz_log, config->zmax_log);
	}

	if (config->zmax_log <= config->zmax){
		fprintf(stderr, "zmax_log parameter must be more than zmax in growth function module\n");
		exit(1);
	}
	
	// status |= c_datablock_get_double(options, OPTION_SECTION, "redshift", config);
	if (status){
		fprintf(stderr, "Please specify the redshift in the growth function module.\n");
		exit(status);
	}
	return config;

} 
    
int execute(c_datablock * block, growth_config * config)
{

	int i,status=0;
	double w,wa,omega_m,omega_v;
	int nz_lin = config->nz_lin;
	int nz_log = config->nz_log;
	int nz = nz_lin + nz_log;
	double zmin_log = config->zmax + config->dz;
	

	//read cosmological params from datablock
	status |= c_datablock_get_double_default(block, cosmo, "w", -1.0, &w);
	status |= c_datablock_get_double_default(block, cosmo, "wa", 0.0, &wa);
	status |= c_datablock_get_double(block, cosmo, "omega_m", &omega_m);
	status |= c_datablock_get_double_default(block, cosmo, "omega_lambda", 1-omega_m, &omega_v);

	if (status){
		fprintf(stderr, "Could not get required parameters for growth function (%d)\n", status);
		return status;
	}

	//allocate memory for single D, f and arrays as function of z
	double *a = malloc(nz*sizeof(double));
	double *dz = malloc(nz*sizeof(double));
	double *fz = malloc(nz*sizeof(double));
	double *z = malloc(nz*sizeof(double));


	// First specify the z bins in increasing z, for my own sanity
	for (i=0;i<nz_lin;i++){
		z[i] = config->zmin + i*config->dz;
	}	

	for (i=0;i<nz_log;i++){
		z[nz_lin + i] = zmin_log*exp(i*(log(config->zmax_log) - log(zmin_log))/(nz_log-1));
	}


	for (i=0;i<nz;i++){
		a[i] = 1.0/(1+z[i]);
	}

	// Now reverse them to increasing a, decreasing z.
	// Note that we don't need to reverse z as we don't use it
	// in the function below
	reverse(a,nz);

	// Compute D and f
	status = get_growthfactor(nz, a, omega_m, omega_v, w, wa, dz, fz);
	
	// Now reverse everything back to increasing z
	// Note that we do not unreverse z as we never reversed it in the first place.
	reverse(a,nz);
	reverse(dz,nz);
	reverse(fz,nz);


	status |= c_datablock_put_double_array_1d(block,growthparameters, "d_z", dz, nz);
	status |= c_datablock_put_double_array_1d(block,growthparameters, "f_z", fz, nz);
	status |= c_datablock_put_double_array_1d(block,growthparameters, "z", z, nz);
	status |= c_datablock_put_double_array_1d(block,growthparameters, "a", a, nz);

	free(fz);
	free(dz);
	free(z);
	free(a);

return status;
}


int cleanup(growth_config * config)
{
	free(config);
	return 0;
}
