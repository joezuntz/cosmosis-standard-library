
typedef struct emu_params {
    double  x[111][8];
    double  xmin[8];
    double  xrange[8];
    double  xmax[8];
    double  z[8];
    double  z_asc[8];
    double  mode[351];
    double  K[2808][45];
    double  w[45][111];
    double  mean[2808];
    double  sd;
    double  beta[45][8];
    double  lamws[45];
    double  lamz[45];
    double KrigBasis[45][111]; // This bit gets generated internally not read from file
} emu_params;

void dump_params(const char * filename, emu_params * p);
void read_params(const char * filename, emu_params * p);
void emuInit(const char * filename, emu_params * par);
int emu(emu_params * par, double *xstar, double *ystar);