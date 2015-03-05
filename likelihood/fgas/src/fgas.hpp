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

#ifndef _FGAS_
#define _FGAS_

#include <iostream>
#include <string>
#include <vector>

#include "clusters.hpp"
#include "lensing.hpp"
#include "xray.hpp"

using std::string;

typedef std::vector<double> Vector;

namespace Clusters {

  class NuisanceParameters {
  public:
    bool do_lensing_calibration, use_fgas;
    double nonthermal, calibration, calibration_evol, calibration_scatter, lensing_systematics;
    // no better place to keep these
    Vector lensing_concentrations, lensing_concentrations500;
    NuisanceParameters();
    inline double calibration_z(const double z) const {
      return lensing_systematics * calibration * (1.0 + z*calibration_evol);
    }
    inline void load_parameters(const double nth, const double cal, const double calev, const double calscat, const double lenssys) {
      nonthermal = nth;
      calibration = cal;
      calibration_evol = calev;
      calibration_scatter = calscat;
      lensing_systematics = lenssys;
    }
  };

  class Cluster {
  public:
    bool do_lensing_calibration, use_fgas;
    string name, data_path;
    double z // redshift
      ,    last_lnP
      ,    ln_lensing_mass
      ,    ln_lensing_c
    ;
    struct {
      NuisanceParameters *nuisance;
      double rho_c, lncalibration;
      NFWmodel mass;
    } for_integrals;
    XrayMeasurements xray;
    LensingMeasurements lensing;
    double lnP(ClusterModel&, NuisanceParameters&, const double dA, const double dL, const double rhocr);
    double lnP_lensing(ClusterModel&, NuisanceParameters&, const double dA, const double rhocr);
    double lnP_lensing_integrand(const double lnMwl);
    bool load_parameters(const double lnMwl, const double lnc);
  };
  typedef std::vector<Cluster*> ClusterVector;

  class Dataset {
  public:
    bool do_lensing_calibration, use_fgas;
    string data_path;
    Vector trial_dA, trial_dL, trial_rhocr;
    ClusterVector clusters;
    bool init(const string &file);
    double lnP(ClusterModel&, NuisanceParameters&);
    void load_trial(double *dAs, double *dLs, double *rhocrs);
    inline double redshift(const int i) const {return clusters[i]->z;}
    void simulate(ClusterModel&, NuisanceParameters&, const bool iscatter=false, const bool mscatter=false);
    inline size_t size() const {return clusters.size();}
  };

  double fgas_prediction(Cluster&, ClusterModel&, NuisanceParameters&, const double dA, const double dL, const double rhocr);

}

std::istream& operator>>(std::istream &is, Clusters::Cluster&);
std::ostream& operator<<(std::ostream &os, Clusters::Cluster&);
std::istream& operator>>(std::istream &is, Clusters::Dataset&);
std::ostream& operator<<(std::ostream &os, Clusters::Dataset&);

#endif
