#define NRANSI
#include "nrutil.h"

/*
 Modified to work on 2nd order differential equations, 
 output now on a vector that must be provided in input, 
 with memory properly allocated ahead of time.
 
 Note nvar *must* be 2.
 
*/


void drkdumb_2ode(double vstart[], int nvar, double x1, double x2, int nstep,
	void (*derivs)(double, double [], double [], cosmology *),double axis_vals[], double solution[], double firstder[], cosmology * cospar)
{
	void drk4(double y[], double dydx[], int n, double x, double h, double yout[],
		void (*derivs)(double, double [], double [], cosmology *), cosmology *cospar);
	int i,k;
	double x,h;
	double *v,*vout,*dv;

	v=dvector(1,nvar);
	vout=dvector(1,nvar);
	dv=dvector(1,nvar);
	for (i=1;i<=nvar;i++) v[i]=vstart[i];
	solution[1]=vstart[1];
	firstder[1]=vstart[2];
	axis_vals[1]=x1;
	x=x1;
	h=(x2-x1)/((double) nstep-1);
	for (k=1;k<nstep;k++) {
		(*derivs)(x,v,dv,cospar);
		drk4(v,dv,nvar,x,h,vout,derivs,cospar);
		if ((double)(x+h) == x) nrerror("Step size too small in routine drkdumb");
		x += h;
		axis_vals[k+1]=x;
		solution[k+1]=vout[1];
		firstder[k+1]=vout[2];
		for(i=1;i<=nvar;i++) v[i]=vout[i];
	}
	free_dvector(dv,1,nvar);
	free_dvector(vout,1,nvar);
	free_dvector(v,1,nvar);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 0!(5'R3. */
