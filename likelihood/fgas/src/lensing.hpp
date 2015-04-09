// Module for X-ray cluster fgas data
// Copyright (c) 2013 Anja von der Linden, Douglas Applegate,
//                    Adam Mantz (amantz@slac.stanford.edu)

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

#ifndef _LENSING_
#define _LENSING_

#include <iostream>
#include <string>
#include <vector>

#include "clusters.hpp"

using std::string;

typedef std::vector<double> Vector;

namespace Clusters {

  class ShearMeasurement {
  public:
    double r_arcmin
      ,    r_Mpc
      ,    g
      ,    sigma_g
      ,    sigma2_g
    ;
    inline void get_radius_Mpc(const double dA) {
      r_Mpc = r_arcmin * dA * (1./60.)*(pi/180.);
    }
  };
  typedef std::vector<ShearMeasurement> ShearProfile;

  class RedshiftPoint {
  public:
    double z, N;
  };
  typedef std::vector<RedshiftPoint> RedshiftHistogram;
  
  class LensingMeasurements {
  public:
    static double c2_over_G; // c^2/G in Msun/Mpc
    double betas, betas2, Sigmacrit;
    ShearProfile shear;
    RedshiftHistogram zhist;
    double chisq(NFWmodel&) const;
    inline bool has_data() {return shear.size()>0;}
    bool load(const string &file);
    void redshift_dependent_calculations(ClusterModel&, const double zcluster, const double dAcluster);
  };

}

std::istream& operator>>(std::istream &is, Clusters::ShearMeasurement&);
std::ostream& operator<<(std::ostream &os, Clusters::ShearMeasurement&);
std::istream& operator>>(std::istream &is, Clusters::RedshiftPoint&);
std::ostream& operator<<(std::ostream &os, Clusters::RedshiftPoint&);
std::ostream& operator<<(std::ostream &os, Clusters::LensingMeasurements&);

#endif

