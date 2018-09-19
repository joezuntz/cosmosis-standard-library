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

    if (config->zmax>2.02){
        fprintf(stderr, "The Cosmic Emulator can only go up to zmax=2.02.  You set zmax=%lf\n", config->zmax);
        fprintf(stderr, "This will not work - quitting\n");
        exit(1);
    }

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
    double k_h[nk];

    int status = 0;
    double wa, h0;

    status |= c_datablock_get_double(block, cosmo, "ommh2",  &xstar[0]);
    status |= c_datablock_get_double(block, cosmo, "ombh2",  &xstar[1]);
    status |= c_datablock_get_double(block, cosmo, "sigma_8",&xstar[2]);
    status |= c_datablock_get_double(block, cosmo, "n_s",    &xstar[4]);

    // We need h0 separately.
    status |= c_datablock_get_double(block, cosmo, "h0",    &h0);
    xstar[3] = h0;

    // w0 is optional - assume -1 if not set.
    status |= c_datablock_get_double_default(block, cosmo, "w", -1.0,  &xstar[5]);

    // Need to do wa separately because the emu code messes
    // with its input.
    status |= c_datablock_get_double_default(block, cosmo, "wa", 0.0,  &wa);
    status |= c_datablock_get_double_default(block, cosmo, "omnuh2",  0.0, &xstar[7]);

    double h3 = h0*h0*h0;

    if (status) return status;

    double** PK = allocate_2d_double(config->nz, nk);

    for (int i=0; i<config->nz; i++){
        xstar[6] = wa; // Have to do this or wa is overwritten each time
        xstar[8] = config->z[i];

        // Run emulator
        status |= emu(xstar, PK[i]);

        // Match camb normalization (h^3)
        for (int j=0; j<nk; j++){
            PK[i][j] *= h3;
        }

        if (status) break;
    }

    // Get the k*h array
    for (int j=0; j<nk; j++){
        k_h[j] = mode[j] / h0;
    }


    if (!status) status |= c_datablock_put_double_grid(block,
            MATTER_POWER_NL_SECTION,
            "z", config->nz, config->z,
            "k_h", nk, k_h, // mode is defined in params.h
            "p_k", PK);

    deallocate_2d_double(&PK, config->nz);

    return status;

}


int cleanup(emu_options * config){
    free(config->z);
    free(config);
    return 0;
}