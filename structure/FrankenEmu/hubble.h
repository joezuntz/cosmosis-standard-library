
struct cosmo{
  double hub;
  double w0;
  double wa;
  double wm;
  double wb;
};


double hubble(struct cosmo mycosmo, double a);
double co_distance(struct cosmo mycosmo, double a);
double ang_dist(struct cosmo mycosmo, double a);
double growth(struct cosmo mycosmo, double a);
double z_lastscattering(struct cosmo mycosmo);
double soundhorizon(struct cosmo mycosmo);
double distls(struct cosmo mycosmo);