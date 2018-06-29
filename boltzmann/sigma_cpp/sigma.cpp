# include <time.h>
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <vector>
# include <omp.h>
# include "./files/sigma.h"
# include "./files/ConstMichel.h"

using namespace std;
using namespace Constants;


extern "C"
int executemain (     
        double  OmM             ,
        int*    int_config      ,
        double* Pk_k            ,
        double* Pk_z            ,
        double* Pk0             ,
        double* z_vec           ,
        double* m_vec           ,
        double* r_vec           ,
        double* sigma_m         
        ) {


    Omega_M = OmM; 
    int proc_num = int_config[0];
    int n_k = int_config[1];
    int nm_bin = int_config[2];
    int nz_bin = int_config[3];
    int nzk_bin = int_config[4];

    omp_set_num_threads(proc_num);

    kmin = Pk_k[0];
    kmax = Pk_k[n_k-1];
    printf("kmin=%.10f\tkmax=%.10f\n",kmin,kmax);

    double tstart, tstop, ttime;
    tstart = (double)clock()/CLOCKS_PER_SEC;
    printf("starting\n");


    if(prt_details==1)printf("Omega_M=%g\n",Omega_M);

    printf(" proc_num       : %d\n", proc_num       );
    printf(" n_k           : %d\n", n_k           );
    printf(" nm_bin         : %d\n", nm_bin         );
    printf(" nz_bin         : %d\n", nz_bin         );
    printf(" nzk_bin         : %d\n", nzk_bin         );

    for (int i=0;i<nm_bin;i++)
        r_vec[i] = radius(pow(10., m_vec[i]));

    VecDoub ps_k(n_k) ,ps(n_k);
    for (int i=0; i<n_k; i++)
            ps_k[i] = Pk_k[i];

    VecDoub ps_z(nzk_bin);
    for (int iz=0;iz<nzk_bin;iz++)
        ps_z[iz] = Pk_z[iz];

    printf("ok\n");

    VecDoub Pk_temp(nzk_bin);
    VecDoub Pk(n_k*nz_bin);
    for (int i=0; i<n_k; i++){

        for (int iz=0;iz<nzk_bin;iz++)
            Pk_temp[iz] = Pk0[iz*n_k+i];

        Spline_interp powerspec(ps_z, Pk_temp);

        for (int iz=0;iz<nz_bin;iz++)
            Pk[iz*n_k+i] = powerspec.interp(z_vec[iz]);
        }
            
    printf("ok\n");

    for (int iz=0;iz<nz_bin;iz++){

        // put power spectrum in spline
        for (int i=0; i<n_k; i++)
            ps[i] = Pk[i+iz*n_k];

        Spline_interp powerspec(ps_k,ps);
    
        // Compute sigma
        struct_sigma function_sigma(powerspec,kmin,kmax,eps_sigma);
    
        for (int i=0;i<nm_bin;i++)
            sigma_m[i+iz*nm_bin] = function_sigma( m_vec[i]*ln10 );
        }

    tstop = (double)clock()/CLOCKS_PER_SEC;
    ttime= tstop-tstart; /*ttime is how long your code run */
    return 0;}
