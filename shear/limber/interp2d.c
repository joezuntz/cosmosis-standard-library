#include <assert.h>
#include "stdio.h"
#include "interp2d.h"


Interpolator2D * load_interp_2d(const char * filename, int N1, int N2)
{
  FILE * infile = fopen(filename,"r");

  if (infile==NULL) {
    fprintf(stderr,"Could not open file %s\n",filename);
    return NULL;
  }
  int N = N1*N2;

  double * x1 = malloc(sizeof(double)*N1);
  double * x2 = malloc(sizeof(double)*N2);
  double * y = malloc(sizeof(double)*N);

  for(int i=0; i<N; i++){
    fscanf(infile, "%lf  %lf  %lf\n", x1+(i/N2), x2+(i%N2), y+i);
  }

  Interpolator2D * interp2d = init_interp_2d(x1, x2, y, N1, N2, gsl_interp_akima);
  free(x1);
  free(x2);
  free(y);
  return interp2d;

}

Interpolator2D * init_interp_2d(double *x1, double *x2, 
  double *y, int N1, int N2, const gsl_interp_type *T) {
  Interpolator2D * interp2d;
  int i,j;
  
  interp2d = (Interpolator2D*)malloc(sizeof(Interpolator2D));
  assert(interp2d != NULL);
  
  //copy in data
  interp2d->N1 = N1;
  interp2d->x1 = (double*)malloc(sizeof(double)*N1);
  assert(interp2d->x1 != NULL);
  for(i=0;i<N1;++i)
    interp2d->x1[i] = x1[i];
  
  interp2d->N2 = N2;
  interp2d->x2 = (double*)malloc(sizeof(double)*N2);
  assert(interp2d->x2 != NULL);
  for(i=0;i<N2;++i)
    interp2d->x2[i] = x2[i];
  
  interp2d->y = (double*)malloc(sizeof(double)*N1*N2);
  assert(interp2d->y != NULL);
  for(i=0;i<N1;++i)
    for(j=0;j<N2;++j)
      interp2d->y[i*N2+j] = y[i*N2+j];
  
  interp2d->y1 = (double*)malloc(sizeof(double)*N1);
  assert(interp2d->y1 != NULL);
  
  //init splines in 2 direction
  interp2d->x2dir_splines = (gsl_spline**)malloc(sizeof(gsl_spline*)*N1);
  assert(interp2d->x2dir_splines != NULL);
  interp2d->x2dir_accels = (gsl_interp_accel**)malloc(sizeof(gsl_interp_accel*)*N1);
  assert(interp2d->x2dir_accels != NULL);
  
  for(i=0;i<N1;++i)
    {
      interp2d->x2dir_splines[i] = gsl_spline_alloc(T,(size_t) N2);
      assert(interp2d->x2dir_splines[i] != NULL);
      
      gsl_spline_init(interp2d->x2dir_splines[i],x2,y+(i*N2),(size_t) N2);
      
      interp2d->x2dir_accels[i] = gsl_interp_accel_alloc();
      assert(interp2d->x2dir_accels[i] != NULL);
    }
  
  //alloc spline in 1 direction
  
	const gsl_interp_type *T_small = T;

  interp2d->x1dir_spline = gsl_spline_alloc(T_small,(size_t) 7);
  assert(interp2d->x1dir_spline != NULL);
  
  interp2d->x1dir_accel = gsl_interp_accel_alloc();
  assert(interp2d->x1dir_accel != NULL);

  return interp2d;
}

double interp_2d(double x1, double x2, Interpolator2D * interp2d)
{
  int i;
  int x1lo,x1hi,Nx1;
  size_t x1ind;
  
  
  if(x1 > interp2d->x1[interp2d->N1-1] || x1 < interp2d->x1[0] ||
     x2 > interp2d->x2[interp2d->N2-1] || x2 < interp2d->x2[0]){
    return 0.0;
  }
  
  x1ind = gsl_interp_bsearch(interp2d->x1,x1,0,interp2d->N1-1);
  x1lo = x1ind - 3;
  x1hi = x1ind + 3;
  if(x1lo < 0)
    {
      x1lo = 0;
      x1hi = 6;
    }
  else if(x1hi >= interp2d->N1)
    {
      x1lo = interp2d->N1-1 - 6;
      x1hi = interp2d->N1-1;
    }
  Nx1 = x1hi - x1lo + 1;
  
  for(i=0;i<Nx1;++i)
    interp2d->y1[i] = gsl_spline_eval(interp2d->x2dir_splines[i+x1lo],x2,interp2d->x2dir_accels[i+x1lo]);
      
  gsl_spline_init(interp2d->x1dir_spline,interp2d->x1+x1lo,interp2d->y1,(size_t) Nx1);
    
  return gsl_spline_eval(interp2d->x1dir_spline,x1,interp2d->x1dir_accel);
}

void destroy_interp_2d(Interpolator2D * interp2d)
{
  int i;
  
  free(interp2d->x1);
  free(interp2d->x2);
  free(interp2d->y);
  free(interp2d->y1);
  
  for(i=0;i<interp2d->N1;++i)
    {
      gsl_spline_free(interp2d->x2dir_splines[i]);
      gsl_interp_accel_free(interp2d->x2dir_accels[i]);
    }
  free(interp2d->x2dir_splines);
  free(interp2d->x2dir_accels);
  
  gsl_spline_free(interp2d->x1dir_spline);
  gsl_interp_accel_free(interp2d->x1dir_accel);
  
  free(interp2d);
}

