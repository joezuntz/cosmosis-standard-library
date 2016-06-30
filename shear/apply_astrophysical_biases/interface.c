#include "cosmosis/datablock/c_datablock.h"
#include "apply_biases.h"
#include "utils.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

void * setup(c_datablock * options){

	bias_config * config = malloc(sizeof(bias_config));
	int status = 0;
	int status_colour = 0;
	char *colour;
	status |= c_datablock_get_int(options, OPTION_SECTION, "verbosity", &(config->verbosity));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "galaxy_bias", true,
		&(config->galaxy_bias));
	status |= c_datablock_get_bool_default(options, OPTION_SECTION, "intrinsic_alignments", true,
		&(config->intrinsic_alignments));
	status_colour |= c_datablock_get_string(options, OPTION_SECTION, "colour", &colour);

	if (status){
		fprintf(stderr, "Please specify intrinsic_alignments and galaxy_bias as true or false.\n");
		exit(status);
	}

	// If a colour is specified a suffix is added to the power spectra and bias names
	// colour_switch is 1 if a colour variable is available and 0 otherwise
	if (!status_colour){
		snprintf(config->colour, 10, "_%s",colour); 
		config->colour_switch = 1;
	}
	else{
		config->colour_switch = 0;
		status_colour = 0; }

	printf("%d\n", status_colour);
	printf("%d \n", config->colour_switch);

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
