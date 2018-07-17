#ifndef _INTERNALVARS_H_
#define _INTERNALVARS_H_

#include <gsl/gsl_spline.h>

namespace PARAMS
{

namespace EPS {
//Precision of integrals
    const long double
    eps_sigma = 1e-6
    ;
};

namespace MATH {
//Mathematical constants
    const long double
    ln10 = log(10.0),
    pi = 3.14159265358979323846264338327950,
    two_pi = 2*pi,
    four_pi = 4*pi,
    four_pi_o3 = four_pi/3,
    pi_over_2 = pi/2,
    natural_e = 2.71828182845904523536028747135266,
    degrees_to_radians = pi/180,
    radians_to_degrees = 180/pi
    ;
};

namespace PHYS {
//Physical constants
    const long double
    rho_c = 2.775e11,// rho_c/h^2
    lightspeed = 299792.458,// in km/s
    c_H0 = 2997.92458,// c/H0
    c_H03 = 26944002417.373989539,// (c/H0)^3
    delta_c = 1.686470199841145450157679146655379605036548423191285671119557
    ;
};

enum class INT_V: int {
// int parameters
    prt_details = 1
};

namespace DOUB_V {
// double parameters
    double
    kmin,
    kmax,
    r
    ;
};

inline namespace COSM {
// cosmological parameters
    double
    Omega_M = 0.274694
    ;
};

gsl_spline * powerspec;

}

#endif
