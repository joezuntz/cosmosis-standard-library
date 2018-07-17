#include <time.h>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <omp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include "InternalVars.h"                                 

double radius(const double M)
{
    return std::pow(M
            /(PARAMS::MATH::four_pi_o3
            *PARAMS::PHYS::rho_c
            *PARAMS::COSM::Omega_M),
            0.333333333333333333333333333333333333333333333333333333333333);	
}
//window function     
double w(const double k)
{
    return 3*((std::sin(k)/k)-std::cos(k))/(k*k);
}
double sig_arg(const double lnk,  void * params) {

    const double 
        k = exp(lnk),
        wkR = w(k*PARAMS::DOUB_V::r);

    return gsl_spline_eval(PARAMS::powerspec, k, NULL)
            *wkR*wkR*k*k*k;
}

gsl_spline* spline_from_vectors(std::vector<double> &x, std::vector<double> &y) {
    int n = x.size();
    gsl_spline *s = gsl_spline_alloc(gsl_interp_cspline, x.size());
    double* px = &x[0];
    double* py = &y[0];
    gsl_spline_init(s, px, py, n);
    return s;
  }

template <class T>
double gsl_integral(T &func, const double a, const double b, const double eps=1.0e-10)
{
  gsl_integration_romberg_workspace * w
    = gsl_integration_romberg_alloc(20);

  double result, error;
  size_t neval = 0;
  void* vd;

  gsl_function F;
  F.function = &func;

  gsl_integration_romberg(&F, a, b, 0, eps, &result, &neval, w);

  gsl_integration_romberg_free (w);

  return result;
}

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

    PARAMS::COSM::Omega_M = OmM; 
    int proc_num = int_config[0];
    int n_k = int_config[1];
    int nm_bin = int_config[2];
    int nz_bin = int_config[3];
    int nzk_bin = int_config[4];

    omp_set_num_threads(proc_num);

    PARAMS::DOUB_V::kmin = Pk_k[0];
    PARAMS::DOUB_V::kmax = Pk_k[n_k-1];
    printf("kmin=%.10f\tkmax=%.10f\n",PARAMS::DOUB_V::kmin,PARAMS::DOUB_V::kmax);

    double tstart, tstop, ttime;
    tstart = (double)clock()/CLOCKS_PER_SEC;
    printf("starting\n");

    //if(PARAMS::INT_V::prt_details==1)printf("Omega_M=%g\n",PARAMS::COSM::Omega_M);

    printf("proc_num:%d\n",proc_num);
    printf("n_k:%d\n",n_k);
    printf("nm_bin:%d\n",nm_bin);
    printf("nz_bin:%d\n",nz_bin);
    printf("nzk_bin:%d\n",nzk_bin);

    // Writes radius to vector

    for (int i=0;i<nm_bin;i++)
        r_vec[i] = radius(pow(10., m_vec[i]));

    // Converts k vec to std::vector

    std::vector<double> ps_k(n_k);
    for (int i=0; i<n_k; i++)
            ps_k[i] = Pk_k[i];

    // Converts k_z vec to std::vector

    std::vector<double> ps_z(nzk_bin);
    for (int iz=0;iz<nzk_bin;iz++)
        ps_z[iz] = Pk_z[iz];

    // Creates a table with Pk(kvals, z_bins)

    std::vector<double> Pk_temp(nzk_bin);
    std::vector<double> Pk(n_k*nz_bin);
    for (int i=0; i<n_k; i++){

        for (int iz=0;iz<nzk_bin;iz++)
            Pk_temp[iz] = Pk0[iz*n_k+i];

        gsl_spline * powerspec;
        powerspec = spline_from_vectors(ps_z, Pk_temp);

        for (int iz=0;iz<nz_bin;iz++)
            Pk[iz*n_k+i] = gsl_spline_eval(powerspec, z_vec[iz], NULL);
    }

    // Compute sigma

    double
        romb = 0;
    const double
        eps = PARAMS::EPS::eps_sigma,
        lnkmin = std::log(PARAMS::DOUB_V::kmin * 1.001),
        lnkmax = std::log(PARAMS::DOUB_V::kmax * 0.999);

    std::vector<double> ps(n_k);

    for (int iz=0;iz<nz_bin;iz++){

        // Gets Pk at right z -> ps
        for (int i=0; i<n_k; i++)
            ps[i] = Pk[i+iz*n_k];

        PARAMS::powerspec = spline_from_vectors(ps_k, ps);
    
        for (int i=0;i<nm_bin;i++){
            PARAMS::DOUB_V::r = r_vec[i];
            romb = gsl_integral(sig_arg, lnkmin, lnkmax, eps);
            sigma_m[i+iz*nm_bin] =  std::sqrt(.5*romb)/PARAMS::MATH::pi;
        }
    }

    tstop = (double)clock()/CLOCKS_PER_SEC;
    ttime= tstop-tstart; /*ttime is how long your code run */
    return 0;}
