# include <time.h>
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <vector>
# include <omp.h>
# include "ConstMichel.h"                                 
# include "base_funcs.h"

using namespace std;

using namespace Constants;
using namespace params;


// calculate mass given radius and vice-versa
double mass(double R){
	return (4*pi/3)*(R*R*R)*rho_c*Omega_M;	
	}
double radius(double M){
	double res= pow(3*M/(4*pi*rho_c*Omega_M),0.333333333333333333333333333333333333333333333333333333333333);	
	return res;
	}
// passes the power spectrum spline to a function
double power_k(Spline_interp powerspec,double k){
       return powerspec.interp(k);
       }      
   
//window function     
double w(double k){
                   return 3*((sin(k)/k)-cos(k))/(k*k);
                   }
//derivarive of the window function times R
double RdwdR(double x){
                   return 3*(((x*x-3)*sin(x)/x)+3*cos(x))/(x*x);
                   }
