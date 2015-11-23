#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include <stdio.h>
#include <math.h>
#include "header.h"

const char * cospar_sec = COSMOLOGICAL_PARAMETERS_SECTION;
const char * growth = GROWTH_PARAMETERS_SECTION;

//Define a struct to contain the options set in the ini file
typedef struct growthconfig {
	int nz;
} growthconfig;

growthconfig * setup(c_datablock * options){
	growthconfig * config;
	return config; 
}

int execute(c_datablock * block, growthconfig * config){

	int status= 0; 
 
	int nz=4000;
	//nz= config->nz;
 
	double zmin,zmax,z,dz;
	double a;

	int status_r;

	//Read cosmological parameters required for growth calculation
	
	double *cospar;
	cospar= (double*)malloc(sizeof(double)*6);

	status |= c_datablock_get_double(block, cospar_sec, "omega_m",&cospar[0]);
	status |= c_datablock_get_double(block, cospar_sec, "omega_lambda",&cospar[1]);
	status_r |= c_datablock_get_double(block, cospar_sec, "omega_r",&cospar[2]);
	if (status_r != 0){
		cospar[2]=0.0;}
	status |= c_datablock_get_double(block, cospar_sec, "w",&cospar[3]);
	status |= c_datablock_get_double(block, cospar_sec, "wa",&cospar[4]);
	status |= c_datablock_get_double(block, cospar_sec, "omega_k",&cospar[5]);

	Initialise_cosmological_parameters(cospar);
	
	fprintf(stderr, "Calculating growth factor. \n");

	double *D_z;
	D_z= (double*)malloc(sizeof(double)*nz);
	double * z_arr;
	z_arr= (double*)malloc(sizeof(double)*nz);
   
	zmin=0.0;
	zmax=4.0;
	dz=(zmax-zmin)/((double) nz-1);   
	for(int i=0 ; i<nz ; i++){
		z=zmin+((double) i)*dz; 
		a=1.0/(1.0+z); 
		z_arr[i]= z;
		D_z[i]= growth_function(a)/growth_function(1.0);
	}

	status |= c_datablock_put_double_array_1d(block,growth,"d_z",D_z,nz);
	status |= c_datablock_put_double_array_1d(block,growth,"z",z_arr,nz);
	
	free(D_z);
	free(z_arr); 
	
	fprintf(stderr,"Done.\n ");
		 
    	return status;
}


int cleanup(void * config)
{
    // Config is whatever you returned from setup above
    // Free it 
//	pktest_config * config = (pktest_config*) config_in;
	free(config);
    return 0;
}
