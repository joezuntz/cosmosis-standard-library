/*
 *  hubble.c
 *  
 *  written by christian wagner
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

#define VC 299792.458
#define TWOPI 6.283185307179586
#define PI    3.141592653589793
#define Tcmb  2.725
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))

#include "hubble.h"

double hubble(struct cosmo mycosmo, double a){
  double H0=mycosmo.hub*100.;
  double Omega_m = mycosmo.wm/pow2(mycosmo.hub);
  double Omega_r  = 4.15e-5/pow2(mycosmo.hub);
  double Omega_X  = 1.-Omega_m;
  double Omega_k = 0.;

  return  H0*sqrt(Omega_r/pow4(a)+Omega_m/pow3(a)+Omega_k/pow2(a)+Omega_X*pow(a,-3*(1+mycosmo.w0+mycosmo.wa))*exp(-3.*mycosmo.wa*(1.-a)));

}


double co_distance_int(double a,void * params){
    double t;
    struct cosmo mycosmo = *(struct cosmo *) params;

    return  1./(a*a*hubble(mycosmo,a));

}

double co_distance(struct cosmo mycosmo, double a){

  gsl_function  F;
   double result;
   double epsabs=0.001;
   double epsrel=0.001;
   double abserr;
   size_t neval;

   F.function = &co_distance_int;
   F.params = &mycosmo;


   gsl_integration_qng (&F, a, 1., epsabs, epsrel, &result, &abserr, &neval);

   return result*VC*mycosmo.hub;

}

double ang_dist(struct cosmo mycosmo, double a){
  return co_distance(mycosmo,a); // assuming flatness
}


int
ode_growth (double a, const double y[], double f[], void *params)
{

  struct cosmo mycosmo = * (struct cosmo *) params;


  double w_a=mycosmo.w0+(1-a)*mycosmo.wa;
  double Omega_m=mycosmo.wm/pow2(mycosmo.hub);
  double Omega_X=1-Omega_m;
  double Xvar=(1-Omega_m)*pow(a,-3.*(1+mycosmo.w0+mycosmo.wa))*exp(-3.*mycosmo.wa*(1.-a))/(Omega_m/(a*a*a)+(1-Omega_m)*pow(a,-3.*(1+mycosmo.w0+mycosmo.wa))*exp(-3.*mycosmo.wa*(1.-a)));


    f[0] = y[1];
    f[1] = -1.5*(1.-w_a)*Xvar*y[0]/(a*a) - (3.5-1.5*w_a*Xvar)*y[1]/a ;

    return GSL_SUCCESS;
}

int
jac_growth (double a, const double y[], double *dfdy, 
     double dfdt[], void *params)
{
    return GSL_SUCCESS;
}


double growth(struct cosmo mycosmo, double a){

    const gsl_odeiv_step_type * T 
	= gsl_odeiv_step_rkf45;
    
    gsl_odeiv_step * s 
	= gsl_odeiv_step_alloc (T, 2);
    gsl_odeiv_control * c 
	= gsl_odeiv_control_y_new (1e-12, 0.0);
    gsl_odeiv_evolve * e 
	= gsl_odeiv_evolve_alloc (2);
    
    double mu = 10;
    gsl_odeiv_system sys = {ode_growth, jac_growth, 2, &mycosmo};
    
    

    double  t1=a;
    double  t = 1.E-12;

    double d1;

    double h = 1e-12;
    double y[2] = { 1.0, 0.0 };

    int status;
    
    while (t < t1){
      status = gsl_odeiv_evolve_apply (e, c, s,
					   &sys, 
					   &t, t1,
					   &h, y);
      if  (h>0.001)
	h=0.001;
      
      
      if (status != GSL_SUCCESS)
	break;
      
    }
	
    d1=y[0]*(t);

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);

    return d1;
}




double z_lastscattering(struct cosmo mycosmo){
  /*
  z_lastscattering(cosmo):
  Returns z_LS from Hu & White, DampingTail paper.
  */
  double wm = mycosmo.wm;
  double wb = mycosmo.wb;
  double b1 = 0.0783*pow(wb,(-0.238))/(1+39.5*pow(wb,0.763));
  double b2 = 0.560/(1+21.1*pow(wb,1.81));
  double zls= 1048.*(1+0.00124*pow(wb,(-0.738)))*(1+b1*pow(wm,b2));
  return zls;
}

double soundhorizon(struct cosmo mycosmo){
  /*
  soundhorizon(cosmo):
  A fit to the sound horizon, in Mpc, from Hu & Sugiyama (1995), Eq. B6
  */
  double wm = mycosmo.wm;
  double wb = mycosmo.wb;
  double wg = 2.4888e-5*pow4(Tcmb/2.73);
  double r0 = 0.75*wb/wg;
  double zeq= 5464.0*(wm/0.135)/pow4(Tcmb/2.725)/(1+0.6851)-1.0;
  double req= r0/(1.+zeq);
  double zLS= z_lastscattering(mycosmo);
  double rls= r0/(1+zLS);
  double tmp= (sqrt(1.+rls)+sqrt(rls+req))/(1.+sqrt(req));
  tmp= 3997.0*sqrt(wg/wm/wb)*log(tmp);
  return(tmp);
}

double distls(struct cosmo mycosmo){
  /*
  distls(cosmo):
  Returns the distance to LS, in Mpc.
  */
  double zLS = z_lastscattering(mycosmo);
  double dLS = ang_dist(mycosmo,1./(1.+zLS))/mycosmo.hub;
  return(dLS);
}


double solvehh(double dLS,struct cosmo mycosmo){
  /*
  solvehh(dLS,cosmo):
  Solves for h given the other cosmological parameters and the distance
  to last scattering (in Mpc).
  */
  struct cosmo hmin = mycosmo;
  hmin.hub = 0.3;
  // printf("min_cosmo: %g %g %g %g %g\n",hmin.hub,hmin.w0,hmin.wa,hmin.wm,hmin.wb);
  double ymin = distls(hmin)-dLS;
  //  printf("ymin: %g\n",ymin);
  struct cosmo hmax = mycosmo;
  hmax.hub = 1.0;
  double ymax = distls(hmax)-dLS;
  //  printf("ymax: %g\n",ymax);
  double ymid;
  struct cosmo hmid = mycosmo;
  while (fabs(hmax.hub-hmin.hub)>1e-3){
    hmid.hub = (hmax.hub+hmin.hub)/2.0;
    ymid = distls(hmid)-dLS;
    if (ymin*ymid<0){
	hmax = hmid;
	ymax = ymid;
      } else {
      hmin = hmid;
      ymin = ymid;
    }
  }
  //  printf("hubble from CMB: %g\n",hmid.hub);
  return hmid.hub;
}

//void getH0fromCMB(double *xstar, double myh, double *rs, double *z_lss, double *d_lss, double *h_CMB, int *writeout, char *fname){
void getH0fromCMB(double *xstar, double *stuff){

  FILE *fp;

  struct cosmo mycosmo={ 0.72, xstar[3], 0, xstar[1], xstar[0]};



  double rs=soundhorizon(mycosmo);
  stuff[0] = rs;

  double z_lss=z_lastscattering(mycosmo);
  stuff[1] = z_lss;
 
  double d_lss=302.4*rs/PI;
  stuff[2] = d_lss;

  double hubble_cmb= solvehh(d_lss,mycosmo);
  stuff[3] = hubble_cmb;
  
}

// Linker function for use with Fortran

void geth0fromcmb_(double *xstar, double *stuff){
  void getH0fromCMB();
  getH0fromCMB(xstar,stuff);
}
