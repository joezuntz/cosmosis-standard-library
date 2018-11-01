#include <stdio.h>
#include <stdlib.h>

const int XI_PLUS = 0;
const int XI_MINUS = 1;
const int GAMMA_T = 2;
const int W_THETA = 3;

#define DEG2RAD 0.017453292519943295769
#define INV4PI 0.079577471545947667884


// This function is defined in wigner_d.c
void wigner_d(int l0, int l1, int n, int m, double theta, double* d);

/**
 * Transform C_ell to xi(theta)
 * NB: Note that the number ell_max is not the number of elements in c_ell: n_ell = ell_max+1
 * because the ell_max is inclusive.
 * Parameters:
 *     int calc_type: 0 (XI_PLUS), 1 (XI_MINUS), 2 (GAMMA_T), 3 (W_THETA)
 *     int ell_max: max ell value in input c_ell
 *     double * c_ell: input c_ell value, starting from zero and up to and including ell_max, length ell_max+1
 *     int n_theta: number of requested theta values
 *     double * theta: requested theta values in RADIANS, must be already allocated with length n_theta
 *     double * xi: : output xi values, must be already allocated with length n_theta
 * Returns:
 *     int status: 0 for OK, 1 for input parameter error
**/
int transform_cl_to_corr(int calc_type, int ell_max, double * c_ell, int n_theta, double * theta, double * xi){

    // Set the values of m1, m2, and ell_min based on the calculation type.
    int m1, m2, ell_min;
    if (calc_type==XI_PLUS){
        m1 = 2;
        m2 = 2;
        ell_min = 2;
    }
    else if (calc_type==XI_MINUS){
        m1 = 2;
        m2 = -2;
        ell_min = 2;
    }
    else if (calc_type==GAMMA_T){
        m1 = 2;
        m2 = 0;
        ell_min = 2;
    }
    else if (calc_type==W_THETA){
        m1 = 0;
        m2 = 0;
        ell_min = 0;
    }
    else{
        fprintf(stderr, "Error, calc_type parameter should be 0, 1, or 2, not %d\n",calc_type);
        return 1;
    }

    int n_ell = ell_max + 1;

    // Space for the Wigner components
    double * wigner = calloc(n_ell, sizeof(double));

    for (int i=0; i<n_theta; i++){
        // Compute the array of Wigner d elements
        // This could be made faster if needed by caching these,
        // at the cost of higher memory usage
        wigner_d(0, ell_max, m1, m2, theta[i], wigner);

        // Sum up the contributions to this value
        xi[i] = 0.0;
        for(int ell = ell_min; ell < n_ell; ell++){
            xi[i] += INV4PI*(2*ell+1)*wigner[ell]*c_ell[ell];
        }
    }

    free(wigner);
    return 0;
}