#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "cosmosis/datablock/c_datablock.h"

#include "params.h"

int emu(double *xstar, double *ystar);

typedef struct emu_options {
    int nz;
    double zmax;
    double * z;
} emu_options;


void* setup(c_datablock * options){
    int status = 0;

    emu_options * config = malloc(sizeof(emu_options));

    status |= c_datablock_get_int_default(options, OPTION_SECTION, "nz", 50, &(config->nz));
    status |= c_datablock_get_double_default(options, OPTION_SECTION, "zmax", 2.01, &(config->zmax));

    config->z = malloc(config->nz*sizeof(double));
    for (int i=0; i<config->nz; i++){
        config->z[i] = i*(config->zmax/(config->nz-1));
    }


    if (status){
        fprintf(stderr,"Please set an integer nz and double zmax in cosmic_emu.\n");
        exit(1);
    }

    return config;
}

int execute(c_datablock * block, emu_options * config){

    const int nparam = 9; // from emu.c
    const int nk = 351; // from emu.c

    const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
    double xstar[nparam];
    double ystar[nk];

    int status = 0;
    double wa;
    status |= c_datablock_get_double(block, cosmo, "ommh2",  &xstar[0]);
    status |= c_datablock_get_double(block, cosmo, "ombh2",  &xstar[1]);
    status |= c_datablock_get_double(block, cosmo, "sigma_8",&xstar[2]);
    status |= c_datablock_get_double(block, cosmo, "h0",     &xstar[3]);
    status |= c_datablock_get_double(block, cosmo, "n_s",    &xstar[4]);

    status |= c_datablock_get_double_default(block, cosmo, "w", -1.0,  &xstar[5]);

    // Need to do wa separately because the emu code messes
    // with its input.
    status |= c_datablock_get_double_default(block, cosmo, "wa", 0.0,  &wa);
    status |= c_datablock_get_double(block, cosmo, "omega_nu",  &xstar[7]);

    if (status) return status;

    double** PK = allocate_2d_double(config->nz, nk);

    for (int i=0; i<config->nz; i++){
        xstar[6] = wa; // Have to do this or wa is overwritten each time
        xstar[8] = config->z[i];
        status |= emu(xstar, PK[i]);
        if (status) break;
    }

    if (status) return status; 

    c_datablock_put_double_grid(block,
            MATTER_POWER_NL_SECTION,
            "z", config->nz, config->z,
            "k_h", nk, mode, // mode is defined in params.h
            "p_k", PK);

    return status;

}


int cleanup(emu_options * config){
    free(config->z);
    free(config);
    return 0;
}