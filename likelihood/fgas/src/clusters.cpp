// Module for X-ray cluster fgas data
// Copyright (c) 2013 Adam Mantz, amantz@slac.stanford.edu

// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:

// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.

#include "clusters.hpp"

#include <fstream>

#include <gsl/gsl_math.h>

using std::atan;
using std::log;
using std::pow;
using std::sqrt;

namespace Clusters {

  double NFWmodel::f_inverse(const double f) {
    const double a1 = 0.5116;
    const double a2 = -0.4283;
    const double a3 = -3.13e-3;
    const double a4 = -3.52e-5;
    double logf, p;
    logf = log(f);
    p = a2 + a3*logf + a4*logf*logf;
    return 1.0/sqrt( a1*pow(f, 2.0*p) + 0.5625) + 2.0*f;
  }

  double NFWmodel::unn_shear(const double x) {
    double x2, a, b, b2, b3, c;
    x2 = pow2(x);
    if (x < 1.0) {
      a = gsl_atanh(sqrt((1-x)/(1+x)));
      b = sqrt(1-x2);
      c = x2 - 1;
      return 8*a/(b*x2) + 4*log(x/2)/x2 - 2/c + 4*a/(b*c);
    } 
    else if (x > 1.0) {
      a = atan(sqrt((x-1)/(1+x)));
      b = sqrt(x2-1);
      b2 = pow2(b);
      b3 = b*b2;
      return 8*a/(b*x2) + 4*log(x/2)/x2 - 2/b2 + 4*a/b3;
    }
    // x == 1
    return 10.0/3.0 + 4.0*log(0.5);
  }

  double NFWmodel::unn_kappa(const double x) {
    double x2, a, b, b2, b3, c;
    x2 = pow2(x);
    if (x < 1.0) {
      a = gsl_atanh(sqrt((1-x)/(1+x)));
      b = sqrt(1-x2);
      c = 1./(x2 - 1);
      return c*(1 - 2.*a/b);
    } 
    else if (x > 1.0) {
      a = atan(sqrt((x-1)/(1+x)));
      b = sqrt(x2-1);
      c = 1./(x2 - 1);
      return c*(1 - 2.*a/b);
    } 
    // x == 1
    return 1./3.;
  }


  double ClusterModel::mHe_mp = 6.64465620e-24/1.67262158e-24; // ratio of alpha mass to proton mass
  
  bool ClusterModel::init(func2 &f, const string& filename) {
    AngularDiameterDistance = f;
    // default values
    U_baryon_0 = 0.824;
    U_baryon_1 = 0.0;
    U_coldgas_0 = 0.010;
    U_coldgas_1 = 0.0;
    fgas_scatter = 0.08;
    massScaling = 0.0;
    massPivot = 3.0e14;
    fb_cosmic = 0.16;
    set_from_YHe(0.24);
    return true;
  }

  bool ClusterModel::load_parameters(const double fb0, const double fb1, const double fc0, const double fc1, const double scat, const double rslope, const double Mslope, const double Mscaling, const double fbcos, const double YHe, const double nenh, const double ntne, const double myoo) {
    U_baryon_0 = fb0;
    U_baryon_1 = fb1;
    U_coldgas_0 = fc0;
    U_coldgas_1 = fc1;
    fgas_scatter = scat;
    fb_prof.index = rslope;
    mtot_prof.index = Mslope;
    massScaling = Mscaling;
    fb_cosmic = fbcos;
    if (YHe >= 0.0) set_from_YHe(YHe);
    else {
      ne_nH = nenh;
      ntot_ne = ntne;
      mu = myoo;
    }
    if (ne_nH < 0.0 || ntot_ne < 0.0 || mu < 0.0) return false;
    return true;
  }

  void ClusterModel::set_from_YHe(const double YHe) {
    // YHe is the primordial mass fraction of He (assume no other species apart from H)
    double nHe_nH = YHe / ((1.0 - YHe) * mHe_mp); // He/H number density
    ne_nH = 1.0 + 2.0 * nHe_nH; // electron/H number
    ntot_ne = 1.0 + (1.0+nHe_nH) / ne_nH; // (electron+baryon)/electron number
    mu = (1.0/ne_nH + mHe_mp*nHe_nH/ne_nH) / ntot_ne; // mean particle mass in units of mp
  }

}

std::istream& operator>>(std::istream &is, Clusters::NFWmodel &m) {
  return is >> m.nfwa >> m.nfwpot;
}

std::ostream& operator<<(std::ostream &os, Clusters::NFWmodel &m) {
  return os << m.nfwa <<' '<< m.nfwpot;
}
