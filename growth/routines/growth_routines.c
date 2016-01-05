#define distance_eps 1e-5

//Cosmological parameters
typedef struct cosmology { 
	double Omega_m;
	double Omega_DE;
	double Omega_r;
	double w0;
	double wa;
	double Omega_k;
	int initialisation_switch;
} cosmology;

// Support variables
typedef struct support {
	int growth_steps;
	double growth_amin;  /* zmax=30.0 */
	double growth_amax;

	double *growth_atable;
	double *growth_gtable;
	double *growth_gder_table;

} support;

#include "routines/numrecipes/nrutil.h"
#include "routines/numrecipes/nrutil.c"

#include "routines/numrecipes/dspline.c"
#include "routines/numrecipes/dsplint.c"

#include "routines/numrecipes/drk4.c"
#include "routines/numrecipes/drkdumb_2ode.c"


const char * cospar_sec = COSMOLOGICAL_PARAMETERS_SECTION;

/* --------------------------------------------------------------------------
 -------------------------------------------------------------------------- 
REQUIRED PARAMETERS:

double Omega_m;       Matter density parameter today
double Omega_DE;      Dark energy density parameter today
double Omega_r;       Radiation density parameter today
double w0;	     First dark energy eos parameter
double wa;	     Dynamic dark energy eos parameter 

 --------------------------------------------------------------------------
 -------------------------------------------------------------------------- */

//Cosmology functions 
double dark_energy(double a, cosmology *cospar);   	/* Returns the dark energy density in units of the critical density today. */

double Evol(double a, cosmology *cospar);   		/* Returns evolution factor E(a)=H(a)/H_0. */

//Growth 
double growth_function(double a, support *suppar, cosmology *cospar);
void initialise_growth(cosmology *cospar, support *suppar);
void growth_support(double a, double *yvec, double *ydervec, cosmology *cospar);

/* ---------------------------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------------------------- */
/* Main functions. */

double dark_energy(double a, cosmology *cospar){
	double val;
	
	val=pow(a,-3.0*(1.0+cospar->w0+cospar->wa));
	val*=exp(-3.0*cospar->wa*(1.0-a));
	val*=cospar->Omega_DE;
 	
 	return(val);
}

double Evol(double a, cosmology *cospar){   

	double val;

	val=cospar->Omega_m/(pow(a,3.0));
	val+=dark_energy(a, cospar);
	val+=cospar->Omega_k/(a*a);      
	val+=cospar->Omega_r/(pow(a,4.0));

	val=sqrt(val);

	return(val);
}

void initialise_growth(cosmology *cospar, support *suppar){
	
	long j;
	double a,dela;
	double *vstart,*vout;
	
	if(cospar->initialisation_switch==0){
		suppar->growth_atable=dvector(1,suppar->growth_steps);
		suppar->growth_gtable=dvector(1,suppar->growth_steps);
		suppar->growth_gder_table=dvector(1,suppar->growth_steps);
	}
	
	vstart=dvector(1,2);
	vstart[1]=1.0;						/* Initialisation of growth is G(a)=1.0. */
	vstart[2]=0.0;						/* First derivative with respect to a is initialised to 0.0. */
	
	drkdumb_2ode(vstart,2,suppar->growth_amin,suppar->growth_amax,suppar->growth_steps, growth_support,
		suppar->growth_atable,suppar->growth_gtable,suppar->growth_gder_table, cospar);

	vstart[1]=suppar->growth_gder_table[1];
	vstart[2]=suppar->growth_gder_table[suppar->growth_steps];	
	dspline(suppar->growth_atable,suppar->growth_gtable,suppar->growth_steps,vstart[1],vstart[2],suppar->growth_gder_table);

	printf("Reset interpolation table for new cosmology.\n");
	cospar->initialisation_switch = 1;

	free_dvector(vstart,1,2);
		
}

double growth_function(double a, support *suppar, cosmology *cospar){
	double ainit,dela,val,err;
	long i;
	static double *vstart,*vout;


	if(choice_of_growth_function==1){
		if(cospar->initialisation_switch==0){initialise_growth(cospar, suppar);}

		dsplint(suppar->growth_atable,suppar->growth_gtable,suppar->growth_gder_table,suppar->growth_steps,a,&val);
		
		val=a*val;
		
		return(val);

	}
	else{
		fprintf(stderr,"\n\n");
		fprintf(stderr,"Error: called growth_function, but transfer function is not Eisenstein & Hu with Baryons.\n");
	}
}

/* Support functions */

void growth_support(double a, double *yvec, double *ydervec, cosmology * cospar){
	double weff;
	double val,val1,oka,odea;
	
	val1=Evol(a, cospar);
	val=a*val1;
	weff=(cospar->w0+cospar->wa)+cospar->wa*(1.0-a)/log(a);
	
	oka=cospar->Omega_k/(val*val);
	odea=cospar->Omega_DE/(val1*val1)/pow(a,3.0*(1.0+weff));

	
	ydervec[1]=yvec[2];
	ydervec[2]=-(2.0*oka+1.5*(1.0-weff)*odea)/(a*a)*yvec[1]-( 1.0+(2.5+0.5*(oka-3.0*weff*odea)) )/a*yvec[2];	
}

cosmology * initialise_cosmological_parameters(c_datablock * block){

	int status= 0;
	int status_r = 0;
	cosmology *cospar = malloc(sizeof(cosmology));

	status |= c_datablock_get_double(block, cospar_sec, "omega_m",&(cospar->Omega_m));
	
	status |= c_datablock_get_double(block, cospar_sec, "omega_lambda",&(cospar->Omega_DE));
	status_r |= c_datablock_get_double(block, cospar_sec, "omega_r",&(cospar->Omega_r));
	if (status_r != 0){
		cospar->Omega_r=0.0;}
	status |= c_datablock_get_double(block, cospar_sec, "w",&(cospar->w0));
	status |= c_datablock_get_double(block, cospar_sec, "wa",&(cospar->wa));
	status |= c_datablock_get_double(block, cospar_sec, "omega_k",&(cospar->Omega_k));

	// initialisation_switch is set to 0 for each new cosmology in
	// the pipeline
	// This tells the above functions that the interpolation table
	// for D(z) needs to be reinitialised for the new parameters
	cospar->initialisation_switch = 0;

	printf("Cosmological parameters:\n");
	printf("omega_m= %f\n", cospar->Omega_m);
	printf("omega_DE= %f\n", cospar->Omega_DE);
	printf("omega_r= %f\n", cospar->Omega_r);  
	printf("omega_K= %f\n", cospar->Omega_k);
	printf("w0= %f\n", cospar->w0); 
	printf("wa= %f\n", cospar->wa);

	return cospar;
}

support * initialise_support(){
	support *suppar  = malloc(sizeof(support));
	// Hard code some basic parameter values that don't need to vary. 
	suppar->growth_steps=100;
	suppar->growth_amin=0.0322581;  /* zmax=30.0 */
	suppar->growth_amax= 1.05;

	return suppar;
}
