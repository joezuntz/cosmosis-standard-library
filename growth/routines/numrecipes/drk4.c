#define NRANSI
#include "nrutil.h"

void drk4(double y[], double dydx[], int n, double x, double h, double yout[],
	void (*derivs)(double, double [], double [], cosmology *), cosmology * cospar)
{
	int i;
	double xh,hh,h6,*dym,*dyt,*yt;

	dym=dvector(1,n);
	dyt=dvector(1,n);
	yt=dvector(1,n);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt,cospar);
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym,cospar);
	for (i=1;i<=n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt,cospar);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	free_dvector(yt,1,n);
	free_dvector(dyt,1,n);
	free_dvector(dym,1,n);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 0!(5'R3. */
