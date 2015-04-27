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

#include "xray.hpp"
#include "util.hpp"

#include <fstream>
#include <iostream>

#include "ConfigFile.h"

namespace Clusters {

  double   XrayConversions::G  = 5.340e+43  // gravitational constant in keV Mpc/Msun^2
    ,      XrayConversions::mp = 8.4089e-58 // proton mass in Msun
    , XrayConversions::Mpc_cgs = 3.0857e24  // 1 Mpc in cm
    ;

  void XrayConversions::evaluate(const double z, const double rpp, const double ne_nH, const double ntot_ne, const double mu, const double dA, const double dL) {
    // measurement = conv * physical
    distance = 1.0 / (rpp*dA);
    density = G * mp * mu * pow2(rpp*dA);
    mtot = G * mp * mu / (rpp*dA);
    mgas = (1.0 + z) / ( 2.0e7 * sqrt(pow5(Mpc_cgs)*pi*ne_nH*pow3(rpp*dA)) * ntot_ne * mu * mp * dL );
  }

  bool XrayMeasurements::load(const string &file) {
    std::ifstream fin;
    fin.open(file.c_str());
    if ( !fin.good() ) return false;
    fin.close();

    ConfigFile config(file);
    string st;

    // for the moment, the redshift and pixel keys are redundant and will be ignored

    // read in NFW table
    if (!config.readInto(st, "nfw")) return false;
    else {
      std::istringstream iss(st);
      while (true) {
	mass.resize(mass.size() + 1);
	lens_cal_radius.resize(lens_cal_radius.size() + 1);
	iss >> mass.back() >> lens_cal_radius.back();
	if (!iss) break;
      }
      lens_cal_radius.pop_back();
      mass.pop_back();
    }
    if (mass.size() == 0) return false;

    return true;
  }

}

std::istream& operator>>(std::istream &is, Clusters::XrayMeasurements &x) {
  is >> x.rpp >> x.x1 >> x.x2 >> x.fgas >> x.ferr >> x.M >> x.Merr >> x.errcor >> x.dA >> x.G_mp_mu_rhocr;
  return is;
}

std::ostream& operator<<(std::ostream &os, Clusters::XrayMeasurements &x) {
  const char s = ' '
    ,        c = '#'
    ;
  return os << x.rpp <<s<< x.x1 <<s<< x.x2 <<s<< x.fgas <<s<< x.ferr <<s<< x.M <<s<< x.Merr <<s<< x.errcor <<s<< x.dA <<s<< x.G_mp_mu_rhocr;// <<s<<c<<s<< "LensCalRadius:"<<x.lens_cal_radius.size() <<s<< "XrayMass:"<<x.mass.size();
}
