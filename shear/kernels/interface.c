// #include "shear_shear.h"
#include "cosmosis/datablock/c_datablock.h"
#include "gsl/gsl_spline.h"
#include "shear_shear.h"
#include "utils.h"
#include <stdio.h>
#include "limber.h"
#include <math.h>

// Short-hand names for the sections we will 
// be looking at
const char * wl_nz = WL_NUMBER_DENSITY_SECTION;
const char * dist = DISTANCES_SECTION;
const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;


void * setup(c_datablock * options){
	return NULL;
}


gsl_spline * get_w_spline(c_datablock * block, int bin, double * z, 
	double chi_max, gsl_spline * a_of_chi_spline)
{
	char name[20];
	double * n_of_z;
	int nz1;
	int status = 0;

	// Load n(z)
	snprintf(name, 20, "bin_%d", bin);
	status |= c_datablock_get_double_array_1d(block, wl_nz, name, &n_of_z,
		&nz1);

	// Make a spline from n(z)
	gsl_spline * n_of_z_spline = spline_from_arrays(nz1, z, n_of_z);

	// Do the main computation
	gsl_spline * W = shear_shear_kernel(chi_max, n_of_z_spline, 
		a_of_chi_spline);

	// tidy up bin-specific data
	gsl_spline_free(n_of_z_spline);
	free(n_of_z);
	return W;
}



int shear_shear_config(c_datablock * block, limber_config * lc){

	// Compute the prefactor, (1.5 Omega_M H0^2)^2
	// The units of this need to be consistent with 
	// what we 		
	double omega_m, H0, h;
	int status = c_datablock_get_double(block, cosmo, "omega_m", &omega_m);
	// in (Mpc/h)^-2
	// We leave the factor of h in because our P(k) 
	// should have it in.  That's why there is a 100 in here.
	const double c_kms = 299792.4580; //km/s
	double scaling = 1.5 * (100.0*100.0)/(c_kms*c_kms) * omega_m; 
	// This scaling value is the bit that goes in front
	// of the W(chi) kernel.  So we square it for the 
	// shear-shear since that has two copies

	// Set up the integration rules.
	// We can interpolate in the logs because
	// everything should be positive
	lc->xlog = true;
	lc->ylog = true;
	lc->n_ell = 100;// printf("only one ell\n");
	lc->prefactor = scaling*scaling;
	lc->ell = malloc(sizeof(double)*lc->n_ell);
	for (int i=0; i<lc->n_ell; i++) lc->ell[i] = 20*pow(1.05,i);
	// lc->ell[0] = 1500.0; printf("Here also\n");
	return status;
}

Interpolator2D * 
load_interpolator(c_datablock * block, gsl_spline * chi_of_z_spline, 
	const char * section,
	const char * k_name, const char * z_name, const char * P_name)
{
	int nk, nz, nP;
	double *k=NULL, *z=NULL;
	double **P=NULL;
	int status = 0;

	status = c_datablock_get_double_grid(block, section, 
		k_name, &nk, &k, 
		z_name, &nz, &z, 
		P_name, &P);

	// What we have now is P(k, z)
	// What we need is P(k, chi)
	// So we loop, lookup, and replace
	for (int i=0; i<nz; i++){
		double zi = z[i];
		double chi_i = gsl_spline_eval(chi_of_z_spline, zi, NULL);
		z[i] = chi_i;
	}

	if (status) return NULL;
	Interpolator2D * interp = init_interp_2d_akima_grid(k, z, P, nk, nz);
	return interp;

}
 




int shear_shear_spectra(c_datablock * block, int nbin, 
	gsl_spline * W[nbin], Interpolator2D * PK)
{

	// Get the prefactor

	limber_config lc;
	int status = shear_shear_config(block, &lc);
	if (status) {free(lc.ell); return status;}

	for (int bin1=1; bin1<=nbin; bin1++){
		for (int bin2=1; bin2<=bin1; bin2++){
			gsl_spline * c_ell = limber_integral(&lc, W[bin1-1], W[bin2-1], PK);
			for (double ell=50; ell<500; ell*=1.1) printf("%le  %le\n",ell, exp(gsl_spline_eval(c_ell, log(ell),NULL)));
		}
	}
	free(lc.ell);


}


int execute(c_datablock * block, void * config)
{
	DATABLOCK_STATUS status=0;
	double * chi;
	double * a;
	double * z;
	int nz1, nz2;
	int nbin;


	// Load the number of bins we 
	status |= c_datablock_get_int(block, wl_nz, "nbin", &nbin);

	// Load z
	status |= c_datablock_get_double_array_1d(block, wl_nz, "z", &z, &nz1);

	// Load another, different z spacing.
	// I am putting this in an array called "a"
	// because I am about to convert it
	status |= c_datablock_get_double_array_1d(block, dist, "z", &a, &nz2);

	// Load chi(z)
	status |= c_datablock_get_double_array_1d(block, dist, "d_m", &chi, &nz2);

	// Reverse ordering so a is increasing - that is what
	// gsl_spline wants
	reverse(a, nz2);
	reverse(chi, nz2);

	// At the moment "a" is still actually redshift z
	gsl_spline * chi_of_z_spline = spline_from_arrays(nz2, a, chi);

	// Replace z->a
	for (int i=0; i<nz2; i++) a[i] = 1.0/(1+a[i]);

	double chi_max = chi[nz2-1];

	// Make spline of a(chi)
	gsl_spline * a_of_chi_spline = spline_from_arrays(nz2, chi, a);

	// Make the W()
	gsl_spline * W_splines[nbin];
	for (int bin=1; bin<=nbin; bin++){
		W_splines[bin-1] = get_w_spline(block, bin, z, chi_max, a_of_chi_spline);
	}

	// Get the P(k) we need
	Interpolator2D * PK = load_interpolator(
		block, chi_of_z_spline, "grid_pk", "k_h", "z", "P_k");
	// This is P(k,z)
	// We need P(k, chi)


	// Make the C_ell and save them
	status |= shear_shear_spectra(block, nbin, W_splines, PK);



	// tidy up global data
	for (int bin=0; bin<nbin; bin++) gsl_spline_free(W_splines[bin]);
	gsl_spline_free(a_of_chi_spline);
	free(chi);
	free(a);
	free(z);

  return status;

}
