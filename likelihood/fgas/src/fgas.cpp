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

#include "fgas.hpp"
#include "util.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_qrng.h>

using std::exp;
using std::log;
using std::log10;
using std::max;
using std::min;
using std::sqrt;

namespace Clusters {

  NuisanceParameters::NuisanceParameters() : do_lensing_calibration(true), nonthermal(1.0), calibration(1.0), calibration_evol(0.0), calibration_scatter(0.17) {
    lensing_concentrations.resize(1000);
    lensing_concentrations500.resize(lensing_concentrations.size());
    gsl_qrng *qrng = gsl_qrng_alloc(gsl_qrng_sobol, 1);
    double d;
    for (int i=0; i<lensing_concentrations.size(); ++i) {
      gsl_qrng_get(qrng, &d);
      lensing_concentrations[i] = exp( log(10.0)*(0.664 + 0.061*gsl_cdf_ugaussian_Pinv(d)) ); // Neto prior
      lensing_concentrations500[i] = 1.0 / NFWmodel::f_inverse( 500.0/200.0 * NFWmodel::f_fcn(1.0/lensing_concentrations[i]) ); // might as well convert to 500 once and be done with it
    }
    gsl_qrng_free(qrng);
  }

  double Cluster::lnP(ClusterModel &m, NuisanceParameters &n, const double dA, const double dL, const double rhocr) {
    xray.conv.evaluate(z, xray.rpp, m.ne_nH, m.ntot_ne, m.mu, dA, dL);
    double fgas_lnP = 0.0;
    if (use_fgas) {
      double var = pow2(m.fgas_scatter) + pow2(xray.ferr) + pow2(m.massScaling*xray.Merr) - 2.0*m.massScaling*xray.errcor*xray.ferr*xray.Merr;
      fgas_lnP = -0.5*( pow2(log(fgas_prediction(*this, m, n, dA, dL, rhocr) / xray.fgas)) / var ) - 0.5*log(var);
    }
    double lens_lnP = lnP_lensing(m, n, dA, rhocr);
    return fgas_lnP + lens_lnP;
  }

  double Cluster::lnP_lensing(ClusterModel &m, NuisanceParameters &n, const double dA, const double rhocr) {
    if (!(xray.has_mass_data() && lensing.has_data())) return 0.0;
    for_integrals.nuisance = &n;
    for_integrals.lncalibration = log( n.calibration_z(z) );
    for_integrals.rho_c = rhocr;
    lensing.redshift_dependent_calculations(m, z, dA);
    return lnP_lensing_integrand(ln_lensing_mass);
  }

  double Cluster::lnP_lensing_integrand(const double lnMwl) {
    double    P = 0.0
      ,  Mwl500 = exp(lnMwl)
      ;

    double r500_lens = pow(3.0*Mwl500 / (4.0*pi*500.0*for_integrals.rho_c), 0.33333);

    int i=0, j=0, n=0;
    const int xstep = 10, lstep=1; // stupid
    double Mx, c_lens, c500_lens, rs_lens, cal_radius_Mpc, Mwl, shear_chisq, cal_chisq, lnratio;

    while (i<xray.mass.size()) {
      ++n;
      // X-ray: calculate mass and convert to Msun
      Mx = xray.mass[i].mass( xray.lens_cal_radius[i] ) / xray.conv.mtot;
      // lensing: get everything in physical units and calculate mass
      c_lens = for_integrals.nuisance->lensing_concentrations[j];
      c500_lens = for_integrals.nuisance->lensing_concentrations500[j];
      rs_lens = r500_lens / c500_lens;
      for_integrals.mass.set_from_rs_c_rhocr(rs_lens, c_lens, for_integrals.rho_c);
      cal_radius_Mpc = xray.lens_cal_radius[i] / xray.conv.distance;
      Mwl = for_integrals.mass.mass(cal_radius_Mpc);
      // lensing data likelihood
      shear_chisq = lensing.chisq(for_integrals.mass);
      // X-ray/lensing ratio likelihood
      lnratio = log(Mwl/Mx);
      cal_chisq = pow2( (for_integrals.lncalibration - lnratio) / for_integrals.nuisance->calibration_scatter );
      // sum
      P += exp(-0.5*(cal_chisq + shear_chisq));
      i += xstep;
      j = (j + lstep) % for_integrals.nuisance->lensing_concentrations.size();
    }

    return log(P / (n*for_integrals.nuisance->calibration_scatter));
  }

  bool Cluster::load_parameters(const double lnMwl, const double lnc) {
    if (lensing.has_data()) {
      ln_lensing_mass = lnMwl;
      ln_lensing_c = lnc;
      return true;
    }
    return false;
  }

  bool Dataset::init(const string &file) {
    std::ifstream fin;
    bool ret;
    // get path to directory containing data file
    char *str = (char*)file.c_str()
      ,  *pch = std::strrchr(str, '/')
      ;
    if (pch != 0) {
      int len = (int)(pch - str) + 1;
      data_path = file.substr(0, len);
    }
    // read it
    fin.open(file.c_str());
    if ( fin.good() ) {
      fin >> *this;
      fin.close();
      return true;
    }
    else return false;
  }

  double Dataset::lnP(ClusterModel &m, NuisanceParameters &n) {
    double ll = 0.0;
    for (int j=0; j<clusters.size(); ++j) ll += ( clusters[j]->last_lnP = clusters[j]->lnP(m, n, trial_dA[j], trial_dL[j], trial_rhocr[j]) );
    return ll;
  }

  void Dataset::load_trial(double *dAs, double *dLs, double *rhocrs) {
    for (int i=0; i<size(); ++i) trial_dA[i] = dAs[i];
    for (int i=0; i<size(); ++i) trial_dL[i] = dLs[i];
    for (int i=0; i<size(); ++i) trial_rhocr[i] = rhocrs[i];
  }

  void Dataset::simulate(ClusterModel &m, NuisanceParameters &n, const bool iscatter, const bool mscatter) {
    for (int j=0; j<clusters.size(); ++j) {
      clusters[j]->xray.conv.evaluate(clusters[j]->z, clusters[j]->xray.rpp, m.ne_nH, m.ntot_ne, m.mu, trial_dA[j], trial_dL[j]);
      clusters[j]->xray.fgas = fgas_prediction(*clusters[j], m, n, trial_dA[j], trial_dL[j], trial_rhocr[j]);
      // not used: intrinsic and/or measurement scatter
    }
  }

  double fgas_prediction(Cluster &c, ClusterModel &m, NuisanceParameters &n, const double dA, const double dL, const double rhocr) {
    double R2ratio = pow2(dA)*XrayConversions::G*XrayConversions::mp*m.mu*rhocr / (pow2(c.xray.dA)*c.xray.G_mp_mu_rhocr) // ratio of square of measurement radius and trial rDelta
      ,      fbar = m.f_baryon(c.z, 1.0, R2ratio)
      ,     fcold = m.f_coldgas(c.z)
      ;
    // NB massPivot should be interpreted as biased mass
    // because of where nonthermal*calibration is placed here
    return max(fbar-fcold, 0.001) * n.nonthermal * n.calibration * c.xray.conv.mgas / c.xray.conv.mtot // magic number
      * pow(c.xray.M/(c.xray.conv.mtot*m.massPivot*m.mtot_prof(1.0, R2ratio)), m.massScaling)
      ;
  }

}

std::istream& operator>>(std::istream &is, Clusters::Cluster &c) {
  is >> c.name >> c.z >> c.xray;
  if (c.do_lensing_calibration) {
    c.xray.load(c.data_path + c.name + ".xray");
    c.lensing.load(c.data_path + c.name + ".lens");
  }
  return is;
}

std::ostream& operator<<(std::ostream &os, Clusters::Cluster &c) {
  const char s = ' ';
  return os << c.name <<s<< c.z <<s<< c.xray <<s<< c.lensing;
}


std::istream& operator>>(std::istream &is, Clusters::Dataset &d) {
  unsigned long Ncl;
  is >> Ncl;
  d.trial_dA.resize(Ncl);
  d.trial_dL.resize(Ncl);
  d.trial_rhocr.resize(Ncl);
  d.clusters.resize(Ncl);
  for (int i=0; i<Ncl; ++i) {
    d.clusters[i] = new Clusters::Cluster();
    d.clusters[i]->data_path = d.data_path;
    d.clusters[i]->do_lensing_calibration = d.do_lensing_calibration;
    d.clusters[i]->use_fgas = d.use_fgas;
    is >> *d.clusters[i];
  }
  return is;
}

std::ostream& operator<<(std::ostream &os, Clusters::Dataset &d) {
  const char eol = '\n';
  os << d.clusters.size();
  for (int i=0; i<d.clusters.size(); ++i) os << eol << *d.clusters[i];
  return os << std::endl;
}
