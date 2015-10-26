#include "cosmosis/datablock/c_datablock.h"
#include "apply_biases.h"
#include "utils.h"
#include <stdio.h>
#include <math.h>


void * setup(c_datablock * options){

	bias_config * config = malloc(sizeof(bias_config));
	int status = 0;
	status |= c_datablock_get_int(options, OPTION_SECTION, "verbosity", &(config->verbosity));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "galaxy_bias", true,
		&(config->galaxy_bias));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "intrinsic_alignments", true,
		&(config->intrinsic_alignments));

	if (status){
		fprintf(stderr, "Please specify intrinsic_alignments and galaxy_bias as true or false.\n");
		exit(status);
	}

	return config;
}

int execute(c_datablock * block, void * config_in)
{
	int status=0;

	bias_config * config = (bias_config*) config_in;
	biases b;

	
	status|=load_biases(block, &b, config);

	status|=get_all_spectra(block, &b, config);	
	free_memory(&b, config);	

	return status;

}

int cleanup(void * config_in)
{
	// Free the memory that we allocated in the
	// setup
	bias_config * config = (bias_config*) config_in;
	free(config);
}
