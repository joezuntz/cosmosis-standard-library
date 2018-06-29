# include <time.h>
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <vector>
# include <omp.h>
# include "ConstMichel.h"
# include "phys_funcs.h"
# define NDIV 1

using namespace std;

using namespace Constants;
using namespace params;


//structure for the argument of sigma that will be intgrated
struct sigma_integrando {
    double r ;
    Spline_interp powerspec1;    

    double operator () (const double & lnk) {
        double k = exp(lnk);
        double wkR = w(k*r);
        double F = power_k(powerspec1,k)*wkR*wkR*k*k*k;
        return F ;
    }
    sigma_integrando (Spline_interp powerspec,double R): powerspec1(powerspec), r(R) { }
};

//Sigma function
double Sigma_R (Spline_interp powerspec, double R, double kmin, double kmax, const Doub eps=1.0e-10){
            sigma_integrando integrando(powerspec,R);
            double romb=integral(integrando,log(kmin),log(kmax),eps);
            return sqrt(.5*romb)/pi;
            }


//structure for the argument of dsigma/dR that will be intgrated
struct dsigma_integrando {
    double r ;
    Spline_interp powerspec1;    

    double operator () (const double & lnk) {
        double k=exp(lnk);
        return power_k(powerspec1,k)*w(k*r)*RdwdR(k*r)*k*k*k ;
    }
    dsigma_integrando (Spline_interp powerspec,double R): powerspec1(powerspec), r(R) { }
};

// dln(sigma)/dlnM * Sigma^2 function
double dlnSigma_dlnM_times_sig2_R (Spline_interp powerspec, double R, double kmin, double kmax, const Doub eps=1.0e-10){
            dsigma_integrando integrando(powerspec,R);
            double romb=integral(integrando,log(kmin),log(kmax),eps);
            return romb/(6*pi*pi);
            }


// structs that allow sigma and dsigma to be writen in vectors 
struct struct_sigma {
    double EPS,KMIN,KMAX ;
    Spline_interp powerspec1;    
    double operator () (const double & lnM) const {
//        printf("part1_struc_sig\n");
        double R=radius(exp(lnM));
        return Sigma_R(powerspec1,R,KMIN,KMAX,EPS) ;
    }
    struct_sigma (Spline_interp powerspec, double kmin, double kmax, const Doub eps=1.0e-10):  powerspec1(powerspec), KMIN(kmin), KMAX(kmax), EPS(eps) { }
};
struct struct_dsigma {
    double EPS,KMIN,KMAX ;
    Spline_interp powerspec1;    
    double operator () (const double & lnM) const {
        double R=radius(exp(lnM));
        return dlnSigma_dlnM_times_sig2_R(powerspec1,R,KMIN,KMAX,EPS) ;
    }
    struct_dsigma (Spline_interp powerspec, double kmin, double kmax, const Doub eps=1.0e-10): powerspec1(powerspec), KMIN(kmin), KMAX(kmax), EPS(eps) { }
};
