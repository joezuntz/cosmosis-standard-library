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

#include "wrapper.hpp"

#include <algorithm>
#include <iostream>

namespace Clusters {

  int Wrapper::init(func2 AngularDistance, const std::string &modelopt) {
    if (model.init(AngularDistance, modelopt)) return true;
    return false;
  }
  
  void Wrapper::load_datapars(IOArray &datapars) {
    int j=0;
    double *p;
    for (int i=0; i<data.size(); ++i) {
      p = &(datapars[j]);
      data[i].load_trial(p, p+data[i].size(), p+2*data[i].size());
      j += 3 * data[i].size(); // quasi-magic
    }
  }

  double Wrapper::lnP(IOArray &datapars) {
    load_datapars(datapars);
    double ll = 0.0;
    for (int j=0; j<data.size(); ++j) ll += data[j].lnP(model, nuisance);
    return ll;
  }

  int Wrapper::load_data(const std::string &file) {
    data.resize(data.size() + 1);
    data.back().do_lensing_calibration = nuisance.do_lensing_calibration;
    data.back().use_fgas = nuisance.use_fgas;
    if ( data.back().init(file) ) {
      nclusters += data.back().size();
      ndatapars += 3 * data.back().size(); // quasi-magic
      return true;
    }
    data.pop_back();
    return false;
  }

  int Wrapper::load_parameters(IOArray &modelpars) {
    // NB num_parameters() function!
    if (!model.load_parameters(modelpars[0], modelpars[1], modelpars[2], modelpars[3], modelpars[4], modelpars[5], modelpars[6], modelpars[7], modelpars[8], modelpars[9], modelpars[10], modelpars[11], modelpars[12])) return false;
    nuisance.load_parameters(modelpars[13], modelpars[14], modelpars[15], modelpars[16], modelpars[17]);
    int i = 18;
    for (int j=0; j<data.size() && i<num_parameters(); ++j) 
      for (int k=0; k<data[j].clusters.size(); ++k)
    	if (data[j].clusters[k]->load_parameters(modelpars[i], 1.4)) i += 1; // second is ignored
    return true;
  }

  void Wrapper::set_const(IOArray &physics) {
    if (physics[0] > 0.0) XrayConversions::G = physics[0];
    if (physics[1] > 0.0) XrayConversions::mp = physics[1];
    if (physics[2] > 0.0) XrayConversions::Mpc_cgs = physics[2];
    if (physics[3] > 0.0) ClusterModel::mHe_mp = physics[3];
    if (physics[4] > 0.0) LensingMeasurements::c2_over_G = physics[4];
  }

  void Wrapper::set_data_options(const int i, const int j) {
    nuisance.do_lensing_calibration = (i != 0);
    nuisance.use_fgas = (j != 0);
  }

  void Wrapper::set_massPivot(const double M) {
    model.massPivot = M;
  }

  void Wrapper::simulate(IOArray &datapars, const bool iscatter, const bool mscatter) {
    load_datapars(datapars);
    for (int i=0; i<data.size(); ++i) {
      data[i].simulate(model, nuisance, iscatter, mscatter);
      std::cout << data[i];
    }    
  }

}

extern "C" {

  Clusters::Wrapper* newClWrapper() {return new Clusters::Wrapper();}

  void freeClWrapper(Clusters::Wrapper *w) {delete w;}

  int ClWrapperDatasetSize(Clusters::Wrapper *w, const int i) {
    return w->dataset_size(i);
  }

  double ClWrapperGetRedshift(Clusters::Wrapper *w, const int i, const int j) {
    return w->get_redshift(i, j);
  }

  int ClWrapperInit(Clusters::Wrapper *w, func2 AngularDistance, const char *opt) {
    std::string s(opt);
    return w->init(AngularDistance, s);
  }

  double ClWrapperLnP(Clusters::Wrapper *w, double *datapars) {
    IOArray d(datapars, datapars + w->num_datapars());
    return w->lnP(d);
  }

  int ClWrapperLoadData(Clusters::Wrapper *w, const char *file) {
    std::string s(file);
    return w->load_data(s);
  }

  int ClWrapperLoadParameters(Clusters::Wrapper *w, double *modelpars) {
    IOArray d(modelpars, modelpars + w->num_parameters());
    return w->load_parameters(d);
  }

  int ClWrapperNumClusters(Clusters::Wrapper *w) {return w->num_clusters();}

  int ClWrapperNumDatapars(Clusters::Wrapper *w) {return w->num_datapars();}

  int ClWrapperNumParameters(Clusters::Wrapper *w) {return w->num_parameters();}

  int ClWrapperNumConst(Clusters::Wrapper *w) {return w->num_const();}

  void ClWrapperSetConst(Clusters::Wrapper *w, double *physics) {
    IOArray d(physics, physics + w->num_const());
    w->set_const(d);
  }

  void ClWrapperSetDataOptions(Clusters::Wrapper *w, const int i, const int j) {
    w->set_data_options(i, j);
  }

  void ClWrapperSetMassPivot(Clusters::Wrapper *w, const double M) {
    w->set_massPivot(M);
  }

  void ClWrapperSimulate(Clusters::Wrapper *w, double *datapars, const int iscatter, const int mscatter)  {
    IOArray d(datapars, datapars + w->num_datapars());
    w->simulate(d, iscatter!=0, mscatter!=0);
  }


}
