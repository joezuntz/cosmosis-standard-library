#include <math.h>
#define EPS 1.0e-5
#define JMAX 50
#define JMAXP (JMAX+1)
#define K 5
double qromo(double (*func)(double), double a, double b,
     double (*choose)(double(*)(double), double, double, int))
{
     void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
     void nrerror(char error_text[]);
     int j;
     double ss,dss,h[JMAXP+1],s[JMAXP];
     h[1]=1.0;
     for (j=1;j<=JMAX;j++) {
          s[j]=(*choose)(func,a,b,j);
          if (j >= K) {
              polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
              if (fabs(dss) <= EPS*fabs(ss)) return ss;
          }
          h[j+1]=h[j]/9.0;
     }
     nrerror("Too many steps in routing qromo");
     return 0.0;
}

