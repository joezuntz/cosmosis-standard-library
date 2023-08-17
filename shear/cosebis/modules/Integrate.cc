#include "Integrate.h"

int cache_table_size = -1;
gsl_integration_glfixed_table * cache_table = NULL;


//gsl Gaussian legendre integration
number gaussianIntegrate_gsl(function_cosebis& f, number a, number b, int N)
{
  if (a==b) return 0.0;

  number result = 0.0;

  if (cache_table == NULL){
    cache_table = gsl_integration_glfixed_table_alloc(N);
    cache_table_size = N;
  }
  else if (cache_table_size != N){
    gsl_integration_glfixed_table_free(cache_table);
    cache_table = gsl_integration_glfixed_table_alloc(N);
    cache_table_size = N;
  }

  number xi;
  number wi;

  for(int m=0;m<N;m++)
  {
    gsl_integration_glfixed_point(a, b, m, &xi, &wi, cache_table);
    result+=wi*f.integrant(xi);
  }
  
  return result;
}
