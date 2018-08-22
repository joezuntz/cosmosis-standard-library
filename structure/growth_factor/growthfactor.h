int growth_de(int n, double *a, double *d_vec, double * f_vec);
double w (double a);
double w_int(double z, void *param);
double DarkEnergy_a( double a );
double Xde_int (double a ,void * params);
double Xde (double a);
int func (double t, const double y[], double f[], void *params);
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);
int get_growthfactor(int n, double *a,double om, double ov, double w, double w2, double *d, double *f);
