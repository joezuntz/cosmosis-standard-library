#include <cmath>
#include <vector>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

namespace PARAMS
{
    const long double eps_sigma = 1e-6;
    const long double ln10 = log(10.0);
    const long double pi = 3.14159265358979323846264338327950;
    const long double pi2 = pi*pi;
    const long double four_pi_o3 = 4.*pi/3;
    const long double rho_c = 2.775e11; // rho_c/h^2
    const long double one_over_3pi2 = 1./(3.*pi*pi);
    double kmin;
    double kmax;
    double r;
    gsl_spline * powerspec;
};

double radius(const double M, const double Omega_M)
{
    return std::pow(M
            /(PARAMS::four_pi_o3
            *PARAMS::rho_c
            *Omega_M),
            0.333333333333333333333333333333333333333333333333333333333333);	
}
//window function     
double w(const double kR)
{
    return 3*((std::sin(kR)/kR)-std::cos(kR))/(kR*kR);
}
//derivate of window function     
double dwdlnkR(double kR){
    return 3*(((kR*kR-3)*std::sin(kR)/kR)+3*std::cos(kR))/(kR*kR);
}
double sig2_arg(const double lnk,  void * params) {
    const double k = exp(lnk);
    const double wkR = w(k*PARAMS::r);
    return gsl_spline_eval(PARAMS::powerspec, k, NULL)*wkR*wkR*k*k*k;
}
double dsig2dlnM_arg(const double lnk,  void * params) {
    const double k = exp(lnk);
    return gsl_spline_eval(PARAMS::powerspec, k, NULL)*w(k*PARAMS::r)*dwdlnkR(k*PARAMS::r)*k*k*k;
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
  double result, error;
  size_t limit = 100;
  int key = 1;

  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc(limit);

  gsl_function F;
  F.function = &func;

  gsl_integration_qag(&F, a, b, 0, eps, limit, key, w, &result, &error);

  gsl_integration_workspace_free (w);

  return result;
}

extern "C"
int executemain (     
        double  OmM,
        int*    int_config,
        double* Pk_k,
        double* Pk_z,
        double* Pk0,
        double* z_vec,
        double* m_vec,
        double* r_vec,
        double* sigma2_m, 
        double* dsigma2dlnM_m
        ) {

    int proc_num = int_config[0];
    int n_k = int_config[1];
    int nm_bin = int_config[2];
    int nz_bin = int_config[3];
    int nzk_bin = int_config[4];
    int prt_detail = int_config[5];
    int calc_dsig = int_config[6];

#ifdef _OPENMP
    omp_set_num_threads(proc_num);
#endif

    PARAMS::kmin = Pk_k[0];
    PARAMS::kmax = Pk_k[n_k-1];

    if (prt_detail == 1){
        printf("kmin=%.10f\tkmax=%.10f\n",PARAMS::kmin,PARAMS::kmax);
        printf("starting\n");
        printf("proc_num:%d\n",proc_num);
        printf("n_k:%d\n",n_k);
        printf("nm_bin:%d\n",nm_bin);
        printf("nz_bin:%d\n",nz_bin);
        printf("nzk_bin:%d\n",nzk_bin);
        }

    // Writes radius to vector
    for (int i=0;i<nm_bin;i++)
        r_vec[i] = radius(pow(10., m_vec[i]), OmM);

    // Converts arrays to std::vector
    std::vector<double> ps_k(Pk_k, Pk_k+n_k);
    std::vector<double> ps_z(Pk_z, Pk_z+nzk_bin);


    // Creates a table with Pk(kvals, z_bins)
    std::vector<double> Pk_temp(nzk_bin);
    std::vector<double> Pk(n_k*nz_bin);
    for (int i=0; i<n_k; i++){

        // Extract values for a fixed k values with varying z
        for (int iz=0;iz<nzk_bin;iz++){
            Pk_temp[iz] = Pk0[iz*n_k+i];
        }

        // Make a spline P(z) for this fixed k value
        gsl_spline * powerspec = spline_from_vectors(ps_z, Pk_temp);

        // Evalulate that spline at our new points
#pragma omp parallel for
        for (int iz=0;iz<nz_bin;iz++){
            Pk[iz*n_k+i] = gsl_spline_eval(powerspec, z_vec[iz], NULL);
        }

        // Free the spline
        gsl_spline_free(powerspec);

    }

    // Compute sigma2
    double romb = 0;
    const double eps = PARAMS::eps_sigma;
    const double lnkmin = std::log(PARAMS::kmin * 1.001);
    const double lnkmax = std::log(PARAMS::kmax * 0.999);

    std::vector<double> ps(n_k);

    for (int iz=0;iz<nz_bin;iz++){

        // Gets Pk at right z -> ps
        for (int i=0; i<n_k; i++){
            ps[i] = Pk[i+iz*n_k];
        }

        PARAMS::powerspec = spline_from_vectors(ps_k, ps);
    
        // Run the integration to calculate these sigma2_m values
        for (int i=0;i<nm_bin;i++){
            PARAMS::r = r_vec[i];
            romb = gsl_integral(sig2_arg, lnkmin, lnkmax, eps);
            sigma2_m[i+iz*nm_bin] =  0.5*romb/PARAMS::pi2;
        }

        if (calc_dsig==1){
            // Run the integration to calculate these dsigma2/dlnm values
            for (int i=0;i<nm_bin;i++){
                PARAMS::r = r_vec[i];
                romb = gsl_integral(dsig2dlnM_arg, lnkmin, lnkmax, eps);
                dsigma2dlnM_m[i+iz*nm_bin] =  romb*PARAMS::one_over_3pi2;
            }
        }

        gsl_spline_free(PARAMS::powerspec);
    }

    return 0;
}
