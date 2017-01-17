/*
 *  interface.c
 *  
 *
 *  Created by Joe Zuntz for CosmoSIS on 7/01/15 based on main.c.  Original header follows:
 * 
 *  Created by Earl Lawrence on 9/17/09.
 *
 *  This program was prepared by Los Alamos National Security, LLC at Los Alamos National Laboratory (LANL) 
 *  under contract No. DE-AC52-06NA25396 with the U.S. Department of Energy (DOE). All rights in the program 
 *  are reserved by the DOE and Los Alamos National Security, LLC.  Permission is granted to the public to 
 *  copy and use this software without charge, provided that this Notice and any statement of authorship are 
 *  reproduced on all copies.  Neither the U.S. Government nor LANS makes any warranty, express or implied, 
 *  or assumes any liability or responsibility for the use of this software.
 *
 */
 
 
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "cosmosis/datablock/c_datablock.h"

#include "hubble.h"

int emu(double *xstar, double *ystar, int *outtype);

#define c_kms 299792.458


typedef struct emu_options {
    int output_type;
    int nz;
    double dz;
    bool do_distances;
} emu_options;

void reverse(double * x, int n)
{
    int i;
    for (i=0; i<n/2; i++){
        double tmp = x[i];
        x[i] = x[n-1-i];
        x[n-1-i] = tmp;
    }
}

void* setup(c_datablock * options){
    int status = 0;
    const int default_output_type = 2; // P(k).
    // Could also add mode for getting transfer functions

    emu_options * config = malloc(sizeof(emu_options));

    status |= c_datablock_get_int_default(options, OPTION_SECTION, "nz", 300, &(config->nz));
    status |= c_datablock_get_double_default(options, OPTION_SECTION, "dz", 0.01, &(config->dz));
    status |= c_datablock_get_bool_default(options, OPTION_SECTION, "do_distances", true, &(config->do_distances));

    if (status){
        fprintf(stderr,"Please set an integer nz and double dz in FrankenEmu.\n");
        exit(1);
    }
    config->output_type=default_output_type;
    return config;
}

int execute(c_datablock * block, emu_options * config) {
    int i,j;
    const int nk = 582;
    double xstar[7], ystar[2*nk];
    int status=0;

    // Load input parameters
    const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
    status |= c_datablock_get_double(block, cosmo, "ombh2",  &xstar[0]);
    status |= c_datablock_get_double(block, cosmo, "ommh2",  &xstar[1]);
    status |= c_datablock_get_double(block, cosmo, "n_s",    &xstar[2]);
    status |= c_datablock_get_double(block, cosmo, "h0",     &xstar[3]);
    status |= c_datablock_get_double(block, cosmo, "w",      &xstar[4]);
    status |= c_datablock_get_double(block, cosmo, "sigma_8",&xstar[5]);
    if (status) return status;

    // Convert h->H
    double h0 = xstar[3];
    xstar[3] = h0 * 100.0;


    // Generate z values to run at and space to save results to
    double ** PK = allocate_2d_double(config->nz, nk);
    double z[config->nz];
    double a[config->nz];
    double d_a[config->nz];
    double d_m[config->nz];
    double d_l[config->nz];
    double mu[config->nz];
    double h_z[config->nz];

    for (i=0; i<config->nz; i++) z[i] = i*config->dz;

    // Convert to the form used in the distance calculator
    struct cosmo mycosmo={ h0, xstar[4], 0.0, xstar[1], xstar[0]};

    double h3 = h0*h0*h0;
    // Run once for each redshift we want
    for (i=0; i<config->nz; i++){
        xstar[6] = z[i];
        a[i] = 1.0/(1.0+z[i]);
        // Run the emulator to get P(k,z)
        status |= emu(xstar, ystar, &(config->output_type));
        for (j=0; j<nk; j++) PK[i][j] = ystar[nk+j]*h3;

        // Also useful to save comoving and angular diameter distances
        h_z[i] = hubble(mycosmo, a[i]) / c_kms; // Convert to (Mpc/c)^{-1}
        d_m[i] = co_distance(mycosmo, a[i]) / h0; //Convert to Mpc
        d_a[i] = ang_dist(mycosmo, a[i]) / (1+z[i]) / h0; //convert to physical angular diameter in Mpc
        d_l[i] = d_m[i]*(1+z[i]);
        mu[i] = 5*log10(d_l[i])+25;
    }

    if (status) return status;






    // Convert output k to k/h
    double k_h[nk];
    for (i=0; i<nk; i++) k_h[i] = ystar[i]/h0;

    //Save output
    const char * nl = MATTER_POWER_NL_SECTION;
    status |= c_datablock_put_double_grid(block, nl, 
        "z", config->nz, z, 
        "k_h", nk, k_h, 
        "p_k", PK);

    if (config->do_distances){
        const char * dist = DISTANCES_SECTION;
        status |= c_datablock_put_double_array_1d(block, dist, "z", z, config->nz);
        status |= c_datablock_put_double_array_1d(block, dist, "a", a, config->nz);
        status |= c_datablock_put_double_array_1d(block, dist, "d_a", d_a, config->nz);
        status |= c_datablock_put_double_array_1d(block, dist, "d_m", d_m, config->nz);
        status |= c_datablock_put_double_array_1d(block, dist, "d_l", d_l, config->nz);
        status |= c_datablock_put_double_array_1d(block, dist, "mu", mu, config->nz);
        status |= c_datablock_put_double_array_1d(block, dist, "h", h_z, config->nz);
    }

    return status;
}

int cleanup(emu_options * config){
    free(config);
    return 0;
}