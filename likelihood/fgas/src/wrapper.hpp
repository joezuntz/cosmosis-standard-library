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

#ifndef _WRAPPER_
#define _WRAPPER_

#include <vector>
#include <string>

#include "clusters.hpp"
#include "fgas.hpp"

typedef std::vector<double> IOArray;

namespace Clusters {

  typedef std::vector<Dataset> DatasetList;

/*
  Workflow:
   init
   load_data
   get_redshifts for each dataset and cluster
   set_phys_const, optionally
   load_parameters
   simulate and/or lnP
 */

  class Wrapper {
  private:
    int nclusters, ndatapars;
    NuisanceParameters nuisance;
    ClusterModel model;
    DatasetList data;
    void load_datapars(IOArray &datapars);
  public:
    // functions that would normally return bools here return ints instead, to facilitate integration with other languages
    Wrapper() : nclusters(0), ndatapars(0) {}
    inline int dataset_size(const int i) const {return data[i].size();}
    inline double get_redshift(const int i, const int j) const {return data[i].clusters[j]->z;}
    int init(func2 AngularDistance, const std::string &modelopt);
    double lnP(IOArray &datapars);
    int load_data(const std::string &file);
    int load_parameters(IOArray &modelpars);
    inline int num_clusters() const {return nclusters;}
    inline int num_datapars() const {return ndatapars;}
    inline int num_parameters() const {return 18+1*12;} // magic
    inline int num_const() const {return 5;} // magic
    void set_data_options(const int, const int);
    void set_const(IOArray &physics);
    void set_massPivot(const double);
    void simulate(IOArray &datapars, const bool iscatter=false, const bool mscatter=false);
  };

}

extern "C" {
  Clusters::Wrapper* newClWrapper();
  void freeClWrapper(Clusters::Wrapper*);
  int ClWrapperDatasetSize(Clusters::Wrapper*, const int);
  double ClWrapperGetRedshift(Clusters::Wrapper*, const int, const int);
  int ClWrapperInit(Clusters::Wrapper*, func2 AngularDistance, const char *opt);
  double ClWrapperLnP(Clusters::Wrapper*, double *datapars);
  int ClWrapperLoadData(Clusters::Wrapper*, const char *file);
  int ClWrapperLoadParameters(Clusters::Wrapper*, double *modelpars);
  int ClWrapperNumClusters(Clusters::Wrapper*);
  int ClWrapperNumDatapars(Clusters::Wrapper*);
  int ClWrapperNumParameters(Clusters::Wrapper*);
  int ClWrapperNumConst(Clusters::Wrapper*);
  void ClWrapperSetConst(Clusters::Wrapper*, double *physics);
  void ClWrapperSetDataOptions(Clusters::Wrapper*, const int, const int);
  void ClWrapperSetMassPivot(Clusters::Wrapper*, const double);
  void ClWrapperSimulate(Clusters::Wrapper*, double *datapars, const int iscatter=0, const int mscatter=0);
}

#endif
