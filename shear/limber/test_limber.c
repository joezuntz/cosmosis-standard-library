#include "stdio.h"
#include "math.h"
#include "limber.h"

// This test program reads the two files supplied 
// in this directory, P(k,chi) and n(chi), and computes
// the galaxy power C_ell with bias=1
// The n(chi) is really an n(z) from the pipeline, not a real
// n(chi) (but the difference is small).  The P is actually a 
// linear LCDM power spectrum.

gsl_spline * load_nchi(const char * filename);

int main(){
	const int NCHI = 301;
	const int NK = 500;

	gsl_spline * n_of_chi = load_nchi("n_of_chi.txt");
	Interpolator2D * P = load_interp_2d_file("p_k_chi.txt", NK, NCHI);

	limber_config config;
	config.xlog = true;
	config.ylog = true;
	config.n_ell = 200;
	config.ell = malloc(sizeof(double)*config.n_ell);
	for (int i=0; i<config.n_ell; i++) config.ell[i] = 9.9 * pow(1.03, i);

	gsl_spline * limber = limber_integral(&config, n_of_chi, n_of_chi, P);
	if (limber == NULL)  exit (1);
    
	for(double ell=10.0; ell<2000; ell*=1.01){
		printf("%le   %le\n", ell, exp(gsl_spline_eval(limber, log(ell), NULL)));
	}
	destroy_interp_2d(P);

}



gsl_spline * load_nchi(const char * filename)
{
	const int NCHI = 301;
	FILE * infile = fopen(filename, "r");
	if (infile==NULL) {
		fprintf(stderr,"Could not open file %s\n",filename);
		exit(1);
	}

	double chi[NCHI];
	double nchi[NCHI];
	for(int i=0; i<NCHI; i++){
		int count = fscanf(infile, "%lf  %lf\n", chi+i, nchi+i);
		if (count!=2){
			fprintf(stderr, "wrong column(s) in file %s on line %d\n",filename, i+1);
			exit(1);
		}
	}
	gsl_spline * output = gsl_spline_alloc(gsl_interp_akima, NCHI);
	gsl_spline_init(output, chi, nchi, NCHI);
	return output;
}
