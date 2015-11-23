#define distance_eps 1e-5

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

//Cosmology 
double dark_energy(double a);   	/* Returns the dark energy density in units of the critical density today. */

double Evol(double a);   		/* Returns evolution factor E(a)=H(a)/H_0. */

//Growth 
double growth_function(double a);
void initialise_growth();
void growth_support(double a, double *yvec, double *ydervec);

/* Support functions and variables */

#define growth_steps 100

double *growth_atable;
double *growth_gtable;
double *growth_gder_table;

#define growth_amin 0.0322581  /* zmax=30.0 */
#define growth_amax 1.05

int initialisation_switch;

//Cosmological parameters 
double Omega_m;       
double Omega_DE;      
double Omega_r;                  
double w0;	     
double wa;	      

double Omega_k; 	 // Curvature energy density

/* ---------------------------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------------------------- */
/* Main functions. */

double dark_energy(double a){
	double val;
	
	val=pow(a,-3.0*(1.0+w0+wa));
	val*=exp(-3.0*wa*(1.0-a));
	val*=Omega_DE;
 	
 	return(val);
}

double Evol(double a){   

  	double val;

  	val=Omega_m/(pow(a,3.0));
  	val+=dark_energy(a);
  	val+=Omega_k/(a*a);      
  	val+=Omega_r/(pow(a,4.0));

  	val=sqrt(val);

  	return(val);
}

void initialise_growth(){
	
	long j;
	double a,dela;
	double *vstart,*vout;
	
	if(growth_atable==NULL){
		growth_atable=dvector(1,growth_steps);
		growth_gtable=dvector(1,growth_steps);
		growth_gder_table=dvector(1,growth_steps);
	}
	
	vstart=dvector(1,2);
	vstart[1]=1.0;						/* Initialisation of growth is G(a)=1.0. */
	vstart[2]=0.0;						/* First derivative with respect to a is initialised to 0.0. */
	
	drkdumb_2ode(vstart,2,growth_amin,growth_amax,growth_steps,&growth_support,growth_atable,growth_gtable,growth_gder_table);

	vstart[1]=growth_gder_table[1];
	vstart[2]=growth_gder_table[growth_steps];	
	dspline(growth_atable,growth_gtable,growth_steps,vstart[1],vstart[2],growth_gder_table);

	printf("Reset interpolation table for new cosmology.\n");
	initialisation_switch = 1;
    
	free_dvector(vstart,1,2);
		
}

double growth_function(double a){
	double ainit,dela,val,err;
	long i;
	static double *vstart,*vout;


	if(choice_of_growth_function==1){
		if(growth_atable==NULL || initialisation_switch==0){initialise_growth();}

		dsplint(growth_atable,growth_gtable,growth_gder_table,growth_steps,a,&val);
		
		val=a*val;
		
		return(val);

	}
	else{
		fprintf(stderr,"\n\n");
		fprintf(stderr,"Error: called growth_function, but transfer function is not Eisenstein & Hu with Baryons.\n");
	}
}

/* Support functions */

void growth_support(double a, double *yvec, double *ydervec){
	double weff;
	double val,val1,oka,odea;
	
	val1=Evol(a);
	val=a*val1;
	weff=(w0+wa)+wa*(1.0-a)/log(a);
	
	oka=Omega_k/(val*val);
	odea=Omega_DE/(val1*val1)/pow(a,3.0*(1.0+weff));

	
	ydervec[1]=yvec[2];
	ydervec[2]=-(2.0*oka+1.5*(1.0-weff)*odea)/(a*a)*yvec[1]-( 1.0+(2.5+0.5*(oka-3.0*weff*odea)) )/a*yvec[2];	
}

void Initialise_cosmological_parameters(double *cospar){

	Omega_m= cospar[0];
	Omega_DE= cospar[1];
	Omega_r= cospar[2];
	w0= cospar[3];
	wa= cospar[4];
	Omega_k= cospar[5];

	// initialisation_switch is set to 0 for each new cosmology in
	// the pipeline
	// This tells the above functions that the interpolation table
	// for D(z) needs to be reinitialised for the new parameters
	initialisation_switch = 0;
		
	printf("Cosmological parameters:\n");
	printf("omega_m= %f\n", Omega_m);
	printf("omega_DE= %f\n", Omega_DE);
	printf("omega_r= %f\n", Omega_r);  
	printf("omega_K= %f\n", Omega_k);
	printf("w0= %f\n", w0); 
	printf("wa= %f\n", wa);
}

  
