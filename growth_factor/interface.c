#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../../cosmosis/datablock/c_datablock.h"
#include "../../cosmosis/datablock/section_names.h"
#include "growthfactor.h"


//Module to calculate the linear growth factor D, and linear growth rate, f. Where D, f are defined by the growth of a
//linear perturbation, delta, with scale factor a: 
//delta(a') = delta(a)*(D(a')/D(a)) and f = dlnD/dlna
//Anyone using Komatsu's CRL library should note: growth_factor_crl = D *(1+z) and growth_rate_crl = f/(1+z)


const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
const char * like = LIKELIHOODS_SECTION;
const char * growthparameters = GROWTH_PARAMETERS_SECTION;

double * setup(c_datablock * options)
{
	int status = 0;
	double * config;
	config = malloc(sizeof(double));
	 status |= c_datablock_get_double(options, OPTION_SECTION, "redshift", config);
        if (status){
                fprintf(stderr, "Please specify the redshift in the growth function module.\n");
                exit(status);
        }
	return config;

} 
    
int execute(c_datablock * block,double * config)
{
	int i,status=0,nzbins = 1000;
	double w,wa,omega_m,z=0;
	double *gf,*dz,*fz,*zbins;
	double redshift = *config;
	
	//allocate memory for single D, f and arrays as function of z
	gf = malloc(2*sizeof(double));
	dz = malloc(nzbins*sizeof(double));
	fz = malloc(nzbins*sizeof(double));
	zbins = malloc(nzbins*sizeof(double));
	//read cosmological params from datablock
        status |= c_datablock_get_double(block, cosmo, "w", &w);
        status |= c_datablock_get_double(block, cosmo, "omega_m", &omega_m);
        status |= c_datablock_get_double(block, cosmo, "wa", &wa);
	if (status){
		wa = 0.0;
		status = 0.0;
	}

	//returns linear growth factor and growth function for flat cosmology with either const w or variable DE eos w(a) = w + (1-a)*wa	
	status = get_growthfactor(1.0/(1.0+redshift),omega_m,w,wa,gf);
	//save to datablock
	status |= c_datablock_put_double(block, growthparameters, "delta", gf[0]);
	status |= c_datablock_put_double(block, growthparameters, "dln_delta_dlna", gf[1]);
	status |= c_datablock_put_double(block, growthparameters, "growth_z", redshift);
	// z=0
	status = get_growthfactor(1.0,omega_m,w,wa,gf);
	status |= c_datablock_put_double(block, growthparameters, "delta_z0", gf[0]);

	//output D and f over a range of z
	for (i=0;i<nzbins;i++)
	{	
		status = get_growthfactor(1.0/(1.0+z),omega_m,w,wa,gf);
		dz[i] = gf[0];
		fz[i] = gf[1];
		zbins[i] = z;
		z = z + 0.1;
	}	
	status |= c_datablock_put_double_array_1d(block,growthparameters, "d_z", dz, nzbins);
	status |= c_datablock_put_double_array_1d(block,growthparameters, "f_z", fz, nzbins);
	status |= c_datablock_put_double_array_1d(block,growthparameters, "zbins", zbins, nzbins);
	free(gf);
	free(fz);
	free(dz);
	free(zbins);

return status;
}


int cleanup(double * config)
{
// Config is whatever you returned from setup above
// Free it 	
 	double * redshift = config;
	free(redshift);

	return 0;
}
