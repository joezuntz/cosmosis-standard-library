# include <time.h>
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <vector>
# include <omp.h>
# include "ConstMichel.h"                                 
# include <vector>

using namespace std;

using namespace Constants;
using namespace params;

//integrates with romberg
template <class T>
double integral(T &func, Doub a, Doub b, const Doub eps=1.0e-10, const double N=1,const char* name="",const double z=0){
	double I=0,delta=(b-a)/N,a_int,b_int;int i=1;
	a_int=a; b_int=a_int+delta  ;
	while(i<=N){i++;
		I=I+qromb(func,a_int,b_int,eps);
		a_int=b_int; b_int=a_int+delta;
		}
	return I;
	}
