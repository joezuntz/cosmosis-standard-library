#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/section_names.h"
#include "src/wrapper.hpp"
#include <gsl/gsl_spline.h>
#include <algorithm>

namespace Clusters {
  const char *sectionName = "fgas";
  const char *paramsSection = "fgas_parameters";
  int Ndatasets;
  cosmosis::DataBlock *theDataBlock; // this will need to be assigned in execute
  double cal_mean, cal_sd, slp_mean, slp_sd, len_mean, len_sd;
  IOArray params, dataparams; // std::vector<double> (see wrapper.hpp)
  // copied from jla_v3 interface
  gsl_spline* spline_from_vectors(std::vector<double> &x, std::vector<double> &y) {
    int n = x.size();
    gsl_spline *s = gsl_spline_alloc(gsl_interp_akima, x.size());
    double* px = &x[0];
    double* py = &y[0];
    gsl_spline_init(s, px, py, n);
    return s;
  }
  class Distances {
  public:
    std::vector<double> vz, vdA, vdL, vrhocr;
    gsl_spline *sdA, *sdL, *srhocr;
    inline void clear() {
      vz.clear(); vdA.clear(); vdL.clear(); vrhocr.clear();
      gsl_spline_free(sdA); gsl_spline_free(sdL); gsl_spline_free(srhocr);
    }
    inline void setup() {
      sdA = spline_from_vectors(vz, vdA);
      sdL = spline_from_vectors(vz, vdL);
      srhocr = spline_from_vectors(vz, vrhocr);
    }
    inline bool backwards() const {
      return (vz[1] < vz[0]);
    }
    inline void reverse() {
      std::reverse(vz.begin(), vz.end());
      std::reverse(vdA.begin(), vdA.end());
      std::reverse(vdL.begin(), vdL.end());
      std::reverse(vrhocr.begin(), vrhocr.end());
    }
    // could use gsl "accelerator" for these
    inline double dA(const double z) const {return gsl_spline_eval(sdA, z, NULL);}
    inline double dL(const double z) const {return gsl_spline_eval(sdL, z, NULL);}
    inline double rhocr(const double z) const {return gsl_spline_eval(srhocr, z, NULL);}
  } distances;
  double angulardistance2(const double z1, const double z2) {
    // Angular diameter distance between two redshifts
    // This function is not needed for the moment (until we release lensing data).
    // Nevertheless, TODO and WARNING!!!
    return 0.0;
  }
  void tantrum(const std::string &s, const bool die=false) {
    std::cerr << s << std::endl;
    if (die) exit(1);
  }
}

extern "C" {

  void *setup(cosmosis::DataBlock *options) {
    // Read options from the CosmoSIS configuration ini file,
    // passed via the "options" argument

    //   initialize the module
    Clusters::Wrapper *wrapper = new Clusters::Wrapper();
    wrapper->init(Clusters::angulardistance2, "");

    // Set any global variables required
    Clusters::theDataBlock = 0;

    // Record any configuration information required

    // TODO: when lensing is added, remove default prior on cl_cal
    options->get_val<double>(Clusters::sectionName, "cl_cal_mean", 0.9, Clusters::cal_mean);
    options->get_val<double>(Clusters::sectionName, "cl_cal_sd", 0.09, Clusters::cal_sd);
    options->get_val<double>(Clusters::sectionName, "fgas_rslope_mean", 0.442, Clusters::slp_mean);
    options->get_val<double>(Clusters::sectionName, "fgas_rslope_sd", 0.035, Clusters::slp_sd);
    options->get_val<double>(Clusters::sectionName, "cl_lenssys_mean", 1.0, Clusters::len_mean);
    options->get_val<double>(Clusters::sectionName, "cl_lenssys_sd", 0.069, Clusters::len_sd);

    // use fgas but not lensing. TODO: make lensing optional when it's available.
    wrapper->set_data_options(0, 1);

    //   load in the data
    Clusters::Ndatasets = 0;
    options->get_val<int>(Clusters::sectionName, "numdatasets", Clusters::Ndatasets);
    if (Clusters::Ndatasets > 9) Clusters::tantrum("Error in fgas module: numdatasets >= 10 is not supported. This is easy to fix, but we're lazy. See fgas_cosmosis.cpp.");
    std::string datafile, key;
    char buffer[3]; // NB size
    for (int i=0; i<Clusters::Ndatasets; ++i) {
      sprintf(buffer, "%d", i+1);
      key = std::string("dataset") + std::string(buffer);
      options->get_val<std::string>(Clusters::sectionName, key, datafile);
      if (!wrapper->load_data(datafile)) Clusters::tantrum("Error in fgas module: load_data(" + datafile + ")", true);
    }
    Clusters::params.resize(wrapper->num_parameters());
    Clusters::dataparams.resize(wrapper->num_datapars());

    //    TODO: set physical constants from cosmosis for consistency

    // Pass back any object you like
    return wrapper;
  }

  DATABLOCK_STATUS execute(cosmosis::DataBlock *block, void *config) {
    // Config is whatever you returned from setup above
    // Block is the collection of parameters and calculations for
    // this set of cosmological parameters
    
    Clusters::Wrapper *wrapper = (Clusters::Wrapper*)config;
    Clusters::theDataBlock = block;
    DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
    const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;
    
    // load up model parameters
    double cl_cal, fgas_rslope, cl_lenssys;
    int i = 0;
    status = block->get_val(Clusters::paramsSection, "U_gas_0", Clusters::params[i++]); // Baryon depletion factor & evolution. We set these to the values for hot gas and fix the cold gas parameters to zero (next).
    if (status) return failure;
    status = block->get_val(Clusters::paramsSection, "U_gas_1", Clusters::params[i++]);
    if (status) return failure;
    Clusters::params[i++] = 0.0; // unused cold gas parameters
    Clusters::params[i++] = 0.0;
    status = block->get_val(Clusters::paramsSection, "fgas_scatter", Clusters::params[i++]);
    if (status) return failure;
    status = block->get_val(Clusters::paramsSection, "fgas_rslope", fgas_rslope);
    if (status) return failure;
    Clusters::params[i++] = fgas_rslope;
    Clusters::params[i++] = 0.0; // unused Mslope and Mscaling
    Clusters::params[i++] = 0.0;
    status = block->get_val(COSMOLOGICAL_PARAMETERS_SECTION, "baryon_fraction", Clusters::params[i++]); // fbaryon
    if (status) return failure;
    status = block->get_val(COSMOLOGICAL_PARAMETERS_SECTION, "yhe", Clusters::params[i++]); // helium mass fraction
    if (status) return failure;
    Clusters::params[i++] = -1.; // alternative to YHe is to specify ne_nH,
    Clusters::params[i++] = -1.; // ntot_ne,
    Clusters::params[i++] = -1.; // and mu
    Clusters::params[i++] = 1.0; // unused nonthermal pressure (separate from cal)
    status = block->get_val(Clusters::paramsSection, "cl_cal", cl_cal);
    if (status) return failure;
    Clusters::params[i++] = cl_cal;
    status = block->get_val(Clusters::paramsSection, "cl_calev", Clusters::params[i++]);
    if (status) return failure;
    status = block->get_val(Clusters::paramsSection, "cl_calscat", Clusters::params[i++]);
    if (status) return failure;
    status = block->get_val(Clusters::paramsSection, "cl_lenssys", cl_lenssys);
    Clusters::params[i++] = cl_lenssys;
    if (status) return failure;
    const int Nlensing = 12; // TODO: make this less stupid
    std::string key;
    char buffer[3]; // NB size
    for (int j=0; j<Nlensing; ++j) {
      sprintf(buffer, "%d", j+1);
      key = std::string("cl_lnMwl_") + std::string(buffer);
      status = block->get_val(Clusters::paramsSection, key, Clusters::params[i++]);
      if (status) return failure;
    }
    if (!wrapper->load_parameters(Clusters::params)) return failure;

    // get distances from cosmosis
    status = block->get_val(DISTANCES_SECTION, "Z", Clusters::distances.vz);
    if (status) return failure;
    status = block->get_val(DISTANCES_SECTION, "D_A", Clusters::distances.vdA); // Mpc
    if (status) return failure;
    status = block->get_val(DISTANCES_SECTION, "D_L", Clusters::distances.vdL); // Mpc
    if (status) return failure;
    status = block->get_val(DISTANCES_SECTION, "H", Clusters::distances.vrhocr);
    if (status) return failure;
    // convert H(z) to critical density
    const double c = 299792458e2; // cm/s
    const double Mpc = 3.0857e24; // cm
    const double pi = 3.14159265358979323846264338328;
    const double G = 6.672e-8; // cm^3/g/s^2
    const double Msun = 1.9891e33; // g
    for (int j=0; j<Clusters::distances.vrhocr.size(); ++j)
      Clusters::distances.vrhocr[j] *= Clusters::distances.vrhocr[j] * c*c * 3.0 * Mpc / (8.0 * pi * G * Msun); // units of Msun/Mpc^3
    // make sure redshift order is increasing
    if (Clusters::distances.backwards()) Clusters::distances.reverse();
    Clusters::distances.setup();

    // calculate data parameters
    i = 0;
    for (int j=0; j<Clusters::Ndatasets; ++j) {
      for (int k=0; k<wrapper->dataset_size(j); ++k) {
	Clusters::dataparams[i++] = Clusters::distances.dA(wrapper->get_redshift(j,k));
      }
      for (int k=0; k<wrapper->dataset_size(j); ++k) {
	Clusters::dataparams[i++] = Clusters::distances.dL(wrapper->get_redshift(j,k));
      }
      for (int k=0; k<wrapper->dataset_size(j); ++k) {
	Clusters::dataparams[i++] = Clusters::distances.rhocr(wrapper->get_redshift(j,k));
      }
    }

    // get the likelihood
    double lnP = wrapper->lnP(Clusters::dataparams);
    Clusters::distances.clear();

    // add Gaussian priors on nuisance parameters
    lnP += -0.5*Clusters::pow2( (cl_cal - Clusters::cal_mean) / Clusters::cal_sd );
    lnP += -0.5*Clusters::pow2( (fgas_rslope - Clusters::slp_mean) / Clusters::slp_sd );
    lnP += -0.5*Clusters::pow2( (cl_lenssys - Clusters::len_mean) / Clusters::len_sd );

    status = block->put_val<double>(LIKELIHOODS_SECTION, "FGAS_LIKE", lnP);
    return status;
  }


  int cleanup(void *config) {
    // Config is whatever you returned from setup above
    // Free it 
    Clusters::Wrapper *wrapper = (Clusters::Wrapper*)config;
    delete wrapper;
    return 0;
  }

} // end of extern C
