#include <assert.h>
#include "stdio.h"
#include "interp2d.h"



Interpolator2D * load_interp_2d_file(const char * filename, int N1, int N2)
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

  int status=0;

  for(int i=0; i<N; i++){
    int count = fscanf(infile, "%lf  %lf  %lf\n", x1+(i/N2), x2+(i%N2), y+i);
    if (count!=3) {
        status=i+1;
        break;
    }

  }

  Interpolator2D * interp2d = NULL;
  if (status==0){
    interp2d = init_interp_2d(x1, x2, y, N1, N2, gsl_interp_akima);
  }
  else{
    fprintf(stderr, "Error reading file at line %d\n", status);
  }
  free(x1);
  free(x2);
  free(y);
  return interp2d;

}

static
Interpolator2D * allocate_interp2d(int N1, int N2, const gsl_interp_type *T){
  int i;

  Interpolator2D * interp2d = (Interpolator2D*)malloc(sizeof(Interpolator2D));
  assert(interp2d != NULL);
  interp2d->N1 = N1;
  interp2d->x1 = (double*)malloc(sizeof(double)*N1);
  interp2d->N2 = N2;
  interp2d->x2 = (double*)malloc(sizeof(double)*N2);
  assert(interp2d->x2 != NULL);
  interp2d->y = (double*)malloc(sizeof(double)*N1*N2);
  assert(interp2d->y != NULL);

  interp2d->y1 = (double*)malloc(sizeof(double)*N1);
  assert(interp2d->y1 != NULL);


  interp2d->x2dir_splines = (gsl_spline**)malloc(sizeof(gsl_spline*)*N1);
  assert(interp2d->x2dir_splines != NULL);

  interp2d->x2dir_accels = (gsl_interp_accel**)malloc(sizeof(gsl_interp_accel*)*N1);
  assert(interp2d->x2dir_accels != NULL);

  for(i=0;i<N1;++i)
    {
      interp2d->x2dir_splines[i] = gsl_spline_alloc(T,(size_t) N2);
      assert(interp2d->x2dir_splines[i] != NULL);
      
      interp2d->x2dir_accels[i] = gsl_interp_accel_alloc();
      assert(interp2d->x2dir_accels[i] != NULL);      
    }

  const gsl_interp_type *T_small = T;

  interp2d->x1dir_spline = gsl_spline_alloc(T_small,(size_t) 7);
  assert(interp2d->x1dir_spline != NULL);
  
  interp2d->x1dir_accel = gsl_interp_accel_alloc();
  assert(interp2d->x1dir_accel != NULL);

  return interp2d;

}

static 
double identity_function(double x, double y, double z, void * a)
{
  return z;
}



static
Interpolator2D * init_interp_2d_core_function(double *x1, double *x2, 
  void *y, int N1, int N2, const gsl_interp_type *T, int from_2d,
  interp2d_modifier_function function, void * args) {

  int i,j;
  Interpolator2D * interp2d = allocate_interp2d(N1,N2, T);
  
  //copy in data
  for(i=0;i<N1;++i)
    interp2d->x1[i] = x1[i];
  
  for(i=0;i<N2;++i)
    interp2d->x2[i] = x2[i];
  
  if (from_2d){
    double **p = (double**) y;
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
        interp2d->y[i*N2+j] = function(x1[i], x2[j], p[i][j], args);
  }
  else {
    double *p = (double*) y;
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
        interp2d->y[i*N2+j] = p[i*N2+j]; 
  }
  

  for(i=0;i<N1;++i)
    {
      gsl_spline_init(interp2d->x2dir_splines[i],x2,interp2d->y+(i*N2),(size_t) N2);
    }
  

  return interp2d;
}

static
Interpolator2D * init_interp_2d_core(double *x1, double *x2, 
  void *y, int N1, int N2, const gsl_interp_type *T, int from_2d) {

  return init_interp_2d_core_function(x1, x2, y, N1, N2, T, from_2d, identity_function, NULL);
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



Interpolator2D * init_interp_2d(double *x1, double *x2, 
  double *y, int N1, int N2, const gsl_interp_type *T)
{
  int from_2d = 0;
  return init_interp_2d_core(x1, x2, y, N1, N2, T, from_2d);

}

Interpolator2D * init_interp_2d_grid(double *x1, double *x2, 
  double **y, int N1, int N2, const gsl_interp_type *T)
{
  int from_2d = 1;
  return init_interp_2d_core(x1, x2, y, N1, N2, T, from_2d);
}


Interpolator2D * init_interp_2d_akima(double *x1, double *x2, double *y, int N1, int N2)
{
  return init_interp_2d(x1, x2, y, N1, N2, gsl_interp_akima);
}

Interpolator2D * init_interp_2d_akima_grid(double *x1, double *x2, double **y, int N1, int N2)
{
  return init_interp_2d_grid(x1, x2, y, N1, N2, gsl_interp_akima);
}






Interpolator2D * init_interp_2d_function(double *x1, double *x2, 
  double *y, int N1, int N2, const gsl_interp_type *T, interp2d_modifier_function function,
  void * args

  )
{
  int from_2d = 0;
  return init_interp_2d_core_function(x1, x2, y, N1, N2, T, from_2d, function, args);

}

Interpolator2D * init_interp_2d_grid_function(double *x1, double *x2, 
  double **y, int N1, int N2, const gsl_interp_type *T, interp2d_modifier_function function,
  void * args)
{
  int from_2d = 1;
  return init_interp_2d_core_function(x1, x2, y, N1, N2, T, from_2d, function, args);
}


Interpolator2D * init_interp_2d_akima_function(double *x1, double *x2, double *y, int N1, int N2,
  interp2d_modifier_function function, void * args)
{
  return init_interp_2d_function(x1, x2, y, N1, N2, gsl_interp_akima, function, args);
}

Interpolator2D * init_interp_2d_akima_grid_function(double *x1, double *x2, double **y, int N1, int N2,
  interp2d_modifier_function function, void * args)
{
  return init_interp_2d_grid_function(x1, x2, y, N1, N2, gsl_interp_akima, function, args);
}



double interp_2d_xmin(Interpolator2D * interp2d){
  return interp2d->x1[0];
}

double interp_2d_xmax(Interpolator2D * interp2d){
  return interp2d->x1[interp2d->N1-1];
}

double interp_2d_ymin(Interpolator2D * interp2d){
  return interp2d->x2[0];
}

double interp_2d_ymax(Interpolator2D * interp2d){
  return interp2d->x2[interp2d->N2-1];
}
