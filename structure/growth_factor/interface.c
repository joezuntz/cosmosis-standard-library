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
	int nz;
        double zmax_log;
        int nz_log;
} growth_config;

growth_config * setup(c_datablock * options)
{
	int status = 0;
	growth_config * config = malloc(sizeof(growth_config));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "zmin", 0.0, &(config->zmin));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "zmax", 3.0, &(config->zmax));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "dz", 0.01, &(config->dz));
	status |= c_datablock_get_double_default(options, OPTION_SECTION, "zmax_log", 3.0, &(config->zmax_log));
	status |= c_datablock_get_int_default(options, OPTION_SECTION, "nz_log", 0, &(config->nz_log));

	config->nz = (int)((config->zmax-config->zmin)/config->dz)+1;
	// status |= c_datablock_get_double(options, OPTION_SECTION, "redshift", config);
	// status |= c_datablock_get_double(options, OPTION_SECTION, "redshift", config);

	printf("Will calculate f(z) and d(z) in %d bins (%lf:%lf:%lf)\n", config->nz, config->zmin, config->zmax, config->dz);
	if (config->nz_log > 0){
	  printf("Will extend growth calculation using %d logarithmically spaced z-bins to zmax = %lf \n", config->nz_log, config->zmax_log);
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
	double w,wa,omega_m;
	double *dz,*fz,*zbins;
	int nz_lin_bins = config->nz;
	int nz_log_bins = config->nz_log;
	int nzbins = nz_lin_bins + nz_log_bins;
	
	//allocate memory for single D, f and arrays as function of z
	double gf[2];
	dz = malloc(nzbins*sizeof(double));
	fz = malloc(nzbins*sizeof(double));
	zbins = malloc(nzbins*sizeof(double));
	//read cosmological params from datablock
        status |= c_datablock_get_double_default(block, cosmo, "w", -1.0, &w);
        status |= c_datablock_get_double_default(block, cosmo, "wa", 0.0, &wa);
        status |= c_datablock_get_double(block, cosmo, "omega_m", &omega_m);
	if (status){
		fprintf(stderr, "Could not get required parameters for growth function (%d)\n", status);
		return status;
	}

	//returns linear growth factor and growth function for flat cosmology with either const w or variable DE eos w(a) = w + (1-a)*wa	
	// status = get_growthfactor(1.0/(1.0+redshift),omega_m,w,wa,gf);
	//save to datablock
	// status |= c_datablock_put_double(block, growthparameters, "delta", gf[0]);
	// status |= c_datablock_put_double(block, growthparameters, "dln_delta_dlna", gf[1]);
	// status |= c_datablock_put_double(block, growthparameters, "growth_z", redshift);
	// z=0
	// status = get_growthfactor(1.0,omega_m,w,wa,gf);
	// status |= c_datablock_put_double(block, growthparameters, "delta_z0", gf[0]);

	//output D and f over a range of z
	for (i=0;i<nz_lin_bins;i++)
	{	
		double z = config->zmin + i*config->dz;
		status = get_growthfactor(1.0/(1.0+z),omega_m,w,wa,gf);
		dz[i] = gf[0];
		fz[i] = gf[1];
		zbins[i] = z;
	}	

	double zmin_log = config->zmax + config->dz;
	for (i=0;i<nz_log_bins;i++)
	{	
	  double z = zmin_log*exp(i*(log(config->zmax_log) - log(zmin_log))/(nz_log_bins-1));
		status = get_growthfactor(1.0/(1.0+z),omega_m,w,wa,gf);
		dz[nz_lin_bins + i] = gf[0];
		fz[nz_lin_bins + i] = gf[1];
		zbins[nz_lin_bins + i] = z;
	}	

	status |= c_datablock_put_double_array_1d(block,growthparameters, "d_z", dz, nzbins);
	status |= c_datablock_put_double_array_1d(block,growthparameters, "f_z", fz, nzbins);
	status |= c_datablock_put_double_array_1d(block,growthparameters, "z", zbins, nzbins);
	free(fz);
	free(dz);
	free(zbins);

return status;
}


int cleanup(growth_config * config)
{
// Config is whatever you returned from setup above
// Free it 	
	free(config);

	return 0;
}
