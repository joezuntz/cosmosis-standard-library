/*CosmoSIS interface file for going from shear C(l) to xi+/-
Uses function tpstat_via_hankel from Martin Kilbinger's nicaea 2.5
(http://cosmostat.org/nicaea). Adapted here by Niall MacCrann (niall.maccrann@gmail.com)
*/
#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include "maths.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

const char * shear_xi = SHEAR_XI_SECTION;
const char * shear_cl = SHEAR_CL_SECTION;
const char * wl_nz = WL_NUMBER_DENSITY_SECTION;


typedef struct cl_to_xi_config {
	int n_theta_user;
	double theta_min_user;
	double theta_max_user;

} cl_to_xi_config;


void * setup(c_datablock * options)
{
	cl_to_xi_config * config = malloc(sizeof(cl_to_xi_config));
	int status = 0;
	//status |= c_datablock_get_int(options, OPTION_SECTION, "n_theta", 10, &(config->n_theta_user));
	//status |= c_datablock_get_double(options, OPTION_SECTION, "theta_min", 0.0, &(config->theta_min_user));
	//status |= c_datablock_get_double(options, OPTION_SECTION, "theta_max", 0.0, &(config->theta_max_user));


	if (status){
		fprintf(stderr, "Please specify n_theta, theta_min, and theta_max in the cl_to_xi module.\n");
		exit(status);
	}

	return config;
}

typedef struct p_projected_data{
	interTable * cl_table;
	//FILE * f;
} p_projected_data;

static double P_projected(void *info, double ell, int bin_i, int bin_j, error ** err)
{
	p_projected_data * data = (p_projected_data*) info;
	interTable * cl_table = data->cl_table;
	//FILE * f = data->f;
	double cl;
	// if (ell<20.0) cl = 0.0;
	// else if (ell>3000.0) cl = 0.0;
	// else
	cl =  exp(interpol_wr(cl_table, log(ell), err));
	//fprintf(f,"%le  %le\n", ell, cl);
	return cl;
}

int execute(c_datablock * block, void * config_in)
{
	DATABLOCK_STATUS status=0;
	int num_z_bin;
	int count=2;
	double ** xi = malloc(count*sizeof(double*));
	for (int i=0; i<count; i++) xi[i] = malloc(sizeof(double)*N_thetaH);

	cl_to_xi_config * config = (cl_to_xi_config*) config_in;

	// Load the number of redshift bins
	status |= c_datablock_get_int(block, wl_nz, "nbin", &num_z_bin);
	//Also load ell array
	double * ell;
	
	int n_ell;
	status |= c_datablock_get_double_array_1d(block,shear_cl, "ell", &ell, &n_ell);
	double log_ell_min = log(ell[0]);
	double log_ell_max = log(ell[n_ell-1]);
	double dlog_ell=(log_ell_max-log_ell_min)/(n_ell-1);
	error * err = NULL;
	interTable* cl_table = init_interTable(n_ell, log_ell_min, log_ell_max, 
		dlog_ell, 1.0, -3.0, &err);

	char name_in[64],name_xip[64],name_xim[64];
    
	double log_theta_min, log_theta_max;

	//Now loop through bin combinations, loading ell,C_ell and computing xi+/-
   	for (int i_bin=1; i_bin<=num_z_bin; i_bin++) {
    	for (int j_bin=i_bin; j_bin<=num_z_bin; j_bin++) {
    		snprintf(name_in, 64, "bin_%d_%d",j_bin,i_bin);
    		snprintf(name_xip, 64, "xiplus_%d_%d",j_bin,i_bin);
    		snprintf(name_xim, 64, "ximinus_%d_%d",j_bin,i_bin);
    		/*FILE *f = fopen("cls.txt", "w");
    		if (f == NULL)
				{
				    printf("Error opening file!\n");
				    exit(1);
				}
    		*/
			//read in C(l)
			double * C_ell;
    		status |= c_datablock_get_double_array_1d(block,shear_cl, name_in, &C_ell, &n_ell);
    		//fill cl_table for P_projected
    		for (int i=0; i<n_ell; i++){
    			cl_table->table[i] = log(C_ell[i]);
    		}
    		tpstat_t tpstat = tp_xipm;
    		p_projected_data d;
    		d.cl_table = cl_table;
    		//d.f = f;
			
			tpstat_via_hankel(&d, xi, &log_theta_min, &log_theta_max,
				       tpstat, &P_projected, i_bin, j_bin, &err);
			//Now save to block
			c_datablock_put_double_array_1d(block, shear_xi, name_xip,
                  xi[0], N_thetaH);
			c_datablock_put_double_array_1d(block, shear_xi, name_xim,
                  xi[1], N_thetaH);
			free(C_ell);
			//fclose(f);
		}
	}
	//Now also save theta values...check these are being calculated correctly at some point
	double dlog_theta= (log_theta_max-log_theta_min)/((double)N_thetaH-1.0);
	double logtheta_center = 0.5*(log_theta_max+log_theta_min);
	int nc = N_thetaH/2+1;
	double theta_vals[N_thetaH]; 
	for (int i; i<N_thetaH; i++){
		theta_vals[i] = exp(log_theta_min+i*dlog_theta);
	}
	c_datablock_put_double_array_1d(block, shear_xi, "theta",
                  theta_vals, N_thetaH);

	//Clean up
	
	for (int i=0; i<count; i++) free(xi[i]);
	free(xi);

	return status;
}

int cleanup(void * config_in)
{
	// Free the memory that we allocated in the 
	// setup
	cl_to_xi_config * config = (cl_to_xi_config*) config_in;
	free(config);	
}





