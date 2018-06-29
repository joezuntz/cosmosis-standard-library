#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

# include "NR/romberg.h"

namespace Constants
{
//Michel's constants
    const double        ln10        = log(10.0);
    // em R=472 dn/dm=0! 
    // em R=27 a precisão de n é de 1e-6 em 30 é 1e-7

    const double        rho_c           = 2.775e11;//  rho_c/h^2
    

//Precision of integrals
    const double        eps_sigma       = 1e-7;
    const double        eps_dsigma      = 1e-5;
    const double        eps_G_z         = 1e-5;

//Values related to pi and e
    const long double   pi          = 3.14159265358979323846264338327950;
    const long double   two_pi      = 2*pi;
    const long double   four_pi     = 4*pi;
    const long double   four_pi_o3  = four_pi/3;
    const long double   pi_over_2   = pi/2;
    
    const long double   natural_e   = 2.71828182845904523536028747135266;

//Conversions for radians to degrees and degrees to radians
    const long double   degrees_to_radians = pi/180;
    const long double   radians_to_degrees = 180/pi;

//Physical constants
    const long double   lightspeed =    299792.458               ;// velocidade da luz em km/s
    const long double   c_H0       =    2997.92458               ;// velocidade da luz por H0 em Mpc/h
    const long double   c_H03      =    26944002417.373989539    ;// velocidade da luz por H0 em Mpc/h ao cubo

}
namespace params
{
/************************************ TO BE UPDATED/FILLED ****************************************/
    int	prt_details = 1 ;
// Cosm params
    double        Omega_M     = 0.274694;

// Other params
    double kmin,kmax;
    double delta_c = 1.686470199841145450157679146655379605036548423191285671119557;

}

#endif
