#ifndef __JLA_H__
#define __JLA_H__

#include <iostream>
#include <string>
#include <vector>

//--------------------------- Full Likelihood ---------------------------------------------

/*
 * Entries of the light-curve parameter file.
 */
struct LCPar{
  std::string name;
  double zcmb, zhel, dz, mb, dmb, x1, dx1, color, dcolor, thirdvar, dthirdvar, cov_m_s, cov_m_c, cov_s_c, set, ra, dec, biascor;
};

/*
 * Store values read from the configuration file.
 */
struct Configuration{
  int version;
  double scriptmcut;
  char* data_file;
  char* C00;
  char* C11;
  char* C22;
  char* C01;
  char* C02;
  char* C12;
};

/*
 * Read one entry of the LC parameter file
 */
std::istream & operator >> (std::istream & is, LCPar & SN);

/*
 * The full JLA likelihood as described in Betoule et al. (2014),
 * Sect. 6.1.
 *
 * See test.cc for a working use case.
 */
class JLALikelihood{
public:

  /*
   * Constructor for the JLA Likelihood object.
   *
   * The object remains empty until the read method is called.
   * 
   * Parameters:
   * -----------
   * verbosity: either 0: completely silent (no error reporting)
   *                   1: errors only
   *                   2: informative (default)
   *                  >2: annoying
   */
  JLALikelihood(int verbosity=2)
    :verbosity(verbosity), C00(NULL), C11(NULL),
    C22(NULL), C01(NULL), C02(NULL), C12(NULL){
  }
  
  ~JLALikelihood();

  /*
   * Load JLA data according to information provided in the provided
   * configuration file (typically data/jla.dataset)
   */
  void configure(Configuration &config);

  /*
   * Give the number of SNe Ia.
   */
  int size();
  
  /*
   * Return the redshift of all SN in the sample.
   */
  double * getZ();
  
  /*
   * Compute the negative log-likelihood of a set of distance modulii
   * given the JLA sample (see Betoule et al. 2014, Eq. 15)
   *
   * Parameters:
   * -----------
   * - distanceModulli: a size N vector of double containing the
   *   predicted distance modulus for each SN in the JLA sample.
   *
   * - nuisanceParameters: a size 4 vector of double containing the
   *   distance estimate nuisance parameters in the order: alpha,
   *   beta, M, deltaM.
   *
   * Return:
   * -------
   * - 2 ln (L) if the computation is sucessfull, NaN otherwise.
   */
  double computeLikelihood(double * distanceModulii, double * nuisanceParameters);

  /*
   * Compute the standardized residuals of the JLA sample to the
   * provided model of distance modulii.
   *
   * Minimisation algorithms specialised in quadratic criterion (such
   * as Levenberg-Marquardt) typically needs the output of this
   * method.
   *
   * Parameters:
   * -----------
   * - distanceModulli: a size N vector of double containing the
   *   predicted distance modulus for each SN in the JLA sample.
   *
   * - nuisanceParameters: a size 4 vector of double containing the
   *   distance estimate nuisance parameters in the order: alpha,
   *   beta, M, deltaM.
   *
   * - residuals: an allocated space for N double. Receive the
   *   standardized residuals r_i at the end of the execution.
   *   The minization criterion is $\chi^2 = \sum_i r_i^2$.
   *
   * Return:
   * -------
   * 0 if the computation is successful, -1 otherwise.
   */
  int computeResiduals(double * distance_moduli, double * nuisance_par, double * residuals);
  
  int verbosity;
  std::vector<LCPar> lcpars;
  double * C00, * C11, * C22, * C01, * C02, * C12;
  double scriptmcut;

};

#endif // __JLA_H__
