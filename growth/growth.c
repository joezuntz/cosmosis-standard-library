#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include <stdio.h>
#include <math.h>
#include <growth.h>
#include "routines/growth_routines.c"

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
 
	double zmin,zmax,z,dz;
	double a;

	//Read cosmological parameters required for growth calculation
	printf("TESTING1\n");
	cosmology *cospar = initialise_cosmological_parameters(block);
	// Initialise support parameters
	printf("TESTING2\n");
	support *suppar = initialise_support();
	
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
		D_z[i]= growth_function(a, suppar, cospar)/growth_function(1.0, suppar, cospar);
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
