/*
 *  emu_noh.c
 *  
 *
 *  Created by Earl Lawrence on 9/17/09.
 *  Update 11/30/2012
 *
 *  This program was prepared by Los Alamos National Security, LLC at Los Alamos National Laboratory (LANL) 
 *  under contract No. DE-AC52-06NA25396 with the U.S. Department of Energy (DOE). All rights in the program 
 *  are reserved by the DOE and Los Alamos National Security, LLC.  Permission is granted to the public to 
 *  copy and use this software without charge, provided that this Notice and any statement of authorship are 
 *  reproduced on all copies.  Neither the U.S. Government nor LANS makes any warranty, express or implied, 
 *  or assumes any liability or responsibility for the use of this software.  
 *
 *  
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#include "pcbasis_noh.h"
#include "design_noh.h"
#include "pcweights_noh.h"
#include "corrlengths_noh.h"
#include "precisions_noh.h"
#include "meansd.h"
#include "kemu.h"
#include "ksim.h"

// Sizes of stuff and number of redshifts
static int m=37, neta=5500, p=5, peta=6, rs=11, nsim=582;
// Kriging basis computed by emuInit, sizes should be peta and m
static double KrigBasis[6][37];


// Initialization function that computes the Kriging basis
void emuInit_noh() {
    int i,j,k,l;
    double cov;
    gsl_matrix *SigmaSim = gsl_matrix_alloc(m,m);
    gsl_vector *b = gsl_vector_alloc(m);
    
    // Do these one principal component at a time
    for(i=0; i<peta; i++) {
        // Fill in the covariance matrix for the principals components
        // Also make a gsl_vector with the weights.
        for(j=0; j<m; j++) {
            // Diagonal
            gsl_matrix_set(SigmaSim, j, j, (1.0/lamz[i]) + (1.0/lamws[i]));
            // Off-diagonals
            for(k=0; k<j; k++) {
                // Compute the covariance
                cov = 0.0;
                for(l=0; l<p; l++) {
                    cov -= beta[i][l]*pow(x[j][l]-x[k][l], 2.0);
                }
                cov = exp(cov)/lamz[i];
                gsl_matrix_set(SigmaSim, j, k, cov);
                gsl_matrix_set(SigmaSim, k, j, cov);
            } // for(k=0; k<j; k++)
            gsl_vector_set(b, j, w[i][j]);
        } // for(j=0; j<m; j++)
        
        // Cholesky and solve
        gsl_linalg_cholesky_decomp(SigmaSim);
        gsl_linalg_cholesky_svx(SigmaSim, b);
        
        // Copy into the Kriging Basis
        for(j=0; j<m; j++) {
            KrigBasis[i][j] = gsl_vector_get(b, j);
        }
    } // for(i=0; i<peta; i++)
    gsl_matrix_free(SigmaSim);
    gsl_vector_free(b);
}

// The actual emulation
// Cosmological parameters, placeholder for the output, type of output
void emu_noh(double *xstar, double *ystar, int *outtype) {
    static int inited=0;
    int i, j, k;
    double wstar[peta], Sigmastar[peta][m], ystaremu[neta], ystar_allz[rs*nsim], logc;
    double xstarstd[p];
    double zemu[rs], ybyz[rs];
    FILE *fp;
    // Interpolation stuff for k and then z
    gsl_spline *lininterp_k = gsl_spline_alloc(gsl_interp_linear, neta/rs);
    gsl_spline *lininterp_z = gsl_spline_alloc(gsl_interp_linear, rs);
    gsl_interp_accel *accel = gsl_interp_accel_alloc();
    
    // Iinitialize if necessary
    if(inited==0) {
        emuInit_noh();
        inited=1;
    }
    
    // Check the inputs to make sure we're interpolating.
    for(i=0; i<p; i++) {
        if((xstar[i] < xmin[i]) || (xstar[i] > xmin[i]+xrange[i])) {
            //printf("The inputs are outside the domain of the emulator.\n");
            switch(i) {
                case 0:
                    printf("omega_b must be between %f and %f.\n", xmin[i], xmin[i]+xrange[i]);
                    break;
                case 1:
                    printf("omega_m must be between %f and %f.\n", xmin[i], xmin[i]+xrange[i]);
                    break;
                case 2:
                    printf("n_s must be between %f and %f.\n", xmin[i], xmin[i]+xrange[i]);
                    break;
                case 3:
                    printf("w must be between %f and %f.\n", xmin[i], xmin[i]+xrange[i]);
                    break;
                case 4:
                    printf("sigma_8 must be between %f and %f.\n", xmin[i], xmin[i]+xrange[i]);
                    break;
            }
            exit(1);
        }
    } // for(i=0; i<p; i++)
    
    // Check redshift to make sure we're interpolating
    if((xstar[p] < 0) || (xstar[p] > 4)) {
        //printf("The inputs are outside the domain of the emulator.\n");
        printf("z must between 0 and 4\n");
        exit(1);
    }
    
    // Standardize the inputs
    for(i=0; i<p; i++) {
        xstarstd[i] = (xstar[i] - xmin[i]) / xrange[i];
    }
    
    // Compute the covariances between the new input and sims for all PCs
    for(i=0; i<peta; i++) {
        for(j=0; j<m; j++) {
            logc = 0.0;
            for(k=0; k<p; k++) {
                logc -= beta[i][k]*pow(x[j][k]-xstarstd[k], 2.0);
            }
            Sigmastar[i][j] = exp(logc)/lamz[i];
        }
    }
    
    // Compute wstar, the predicted PC weights for the new input
    for(i=0; i<peta; i++) {
        wstar[i]=0.0;
        for(j=0; j<m; j++) {
            wstar[i] += Sigmastar[i][j] * KrigBasis[i][j];
        }
    }
    
    // Compute ystar, the new output
    for(i=0; i<neta; i++) {
        ystaremu[i] = 0.0;
        for(j=0; j<peta; j++) {
            ystaremu[i] += K[i][j]*wstar[j];
        }
        ystaremu[i] = ystaremu[i]*sd + mean[i];
    }
    
    //printf("Emulated\n");
    
    
    // Interpolate the emulated output onto the original domain.
    for(i=0; i<rs; i++) {
        gsl_spline_init(lininterp_k, kemu, &ystaremu[i*neta/rs], neta/rs);
        for(j=0; j<nsim; j++) {
            ystar_allz[i*nsim+j] = gsl_spline_eval(lininterp_k, ksim[j], accel);
        }
        gsl_interp_accel_reset(accel);
    }
    
    //printf("k interped\n");
    
    // Fill in the k values for the final output
    for(i=0; i<nsim; i++) {
        ystar[i] = ksim[i];
    }
    
    // Interpolate on to the desired redshift
    // The order needs to be reversed here.
    for(i=0; i<rs; i++) {
        zemu[i] = (1/aemu[rs-i-1]) - 1.0;
    }
    for(i=0; i<nsim; i++) {
        // Build an array with the values of y for a given value of z
        // Reverse the order
        for(j=0; j<rs; j++) {
            ybyz[rs-j-1] = ystar_allz[j*nsim+i];
        }
        gsl_spline_init(lininterp_z, zemu, ybyz, rs);
        ystar[nsim+i] = gsl_spline_eval(lininterp_z, xstar[p], accel);
        gsl_interp_accel_reset(accel);
    }
    
    //printf("z interped\n");
    
    switch(outtype[0]) {
        default:
	   for(j=0; j<nsim; j++) {
	     ystar[nsim+j] = pow(10.0,ystar[nsim+j]);
	   }

            break;
        case 1:
            // Transform to Delta^2
            for(j=0; j<nsim; j++) {
                ystar[nsim+j] = ystar[nsim+j] + 1.5*log10(ksim[j]);
                ystar[nsim+j] = pow(10.0,ystar[nsim+j]);
            }
            break;
        case 2:
            // Transform to P(k)
            for(j=0; j<nsim; j++) {
                ystar[nsim+j] = ystar[nsim+j] - 1.5*log10(ksim[j]);
                ystar[nsim+j] = pow(10.0,ystar[nsim+j])*2.0*M_PI*M_PI;
            }
            break;
    }
    
    // Free some stuff.  I always forget this.  
    // Thanks to Tim Eifler for discovering it (apparently the hard way).
    gsl_spline_free(lininterp_z);
    gsl_spline_free(lininterp_k);
    gsl_interp_accel_free(accel);
}

// Linker function for use with Fortran
void emu_noh_(double *x, double *y, int *outtype) {
    void emu();
    emu(x,y,outtype);
}

