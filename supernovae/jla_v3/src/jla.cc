#include "jla.h"


#include <iostream>
#include <fstream>
// #include <cstring>
// #include <cstdlib>
#include <gsl/gsl_cblas.h>
#include <cmath>

#define sq(x) x*x


//---------------------------- LAPACK headers ------------------------------------

#define dpotrf dpotrf_
#define dtrtrs dtrtrs_

extern "C" {
  void dpotrf(const char * UPLO, int * N, double * A, int * LDA, int * INFO);
  void dtrtrs(const char * UPLO, const char * TRANS, const char * DIAG, int * N,
               int * NRHS, double * A, int * LDA, double * B, int * LDB, int *INFO);
}

using namespace std;


//-------------- Reading configuration and data from disk --------------------------

/*
 * Read one entry in the lightcurve file.
 */
istream & operator >> (istream & is, LCPar & sn)
{
  is >> sn.name;
  is >> sn.zcmb >> sn.zhel >> sn.dz >> sn.mb >> sn.dmb >> sn.x1 >> sn.dx1 >> sn.color
     >> sn.dcolor >> sn.thirdvar >> sn.dthirdvar >> sn.cov_m_s >> sn.cov_m_c
     >> sn.cov_s_c >> sn.set >> sn.ra >> sn.dec >> sn.biascor;
  return is;
}


/*
 * Read square matrices from text files.
 *
 * Matrix file must constain space separated values, starts with one
 * integer N (the matrix size) followed by at least N x N floating
 * point values. The rest of the file is ignored.
 *
 * Return a pointer on a N x N double region containing the matrix
 * values.
 */
double * read_matrix(const char * filename, int verbosity)
{
  int N;
  double * mat = NULL;

  std::ifstream fid(filename);

  if (!fid){
    fid.close();
    if (verbosity > 0)
      cerr << "Fail to read file " << filename << endl;
    exit(-1);
    return mat;
  }

  if (verbosity > 1)
    cout << "Reading file " << filename << endl;

  fid >> N;
  mat = (double *)malloc(sizeof(double) * N * N);

  for (int i = 0; i < N * N; i++)
    fid >> mat[i];

  if (verbosity > 1)
    cout << N << " x " << N << " values read from file " << filename << endl;

  fid.close();
  return mat;
}

/* 
 *  JAZ New read method based instead on block input
 * 
 *
 */
void JLALikelihood::configure(Configuration &config)
{
  C00 = read_matrix(config.C00, verbosity);
  C11 = read_matrix(config.C11, verbosity);
  C22 = read_matrix(config.C22, verbosity);
  C01 = read_matrix(config.C01, verbosity);
  C02 = read_matrix(config.C02, verbosity);
  C12 = read_matrix(config.C12, verbosity);


  char buffer[1024];

  if (verbosity > 1)
    cout << "Reading light-curve parameters from '" << config.data_file << endl;
  ifstream fid(config.data_file);
  if (!fid){
    if (verbosity > 0)
      cerr << "Can't load '" << config.data_file << "'\n";
    exit(-1);
  }
    
  while(fid.get() == '#')
    fid.getline(buffer, 1024);
  fid.unget();
  while (fid)
    {
      LCPar SN;
      fid >> SN;
      if (fid)
        lcpars.push_back(SN);
    }

  scriptmcut = config.scriptmcut;
   
  if (verbosity > 1)
    cout << "Read " << lcpars.size() << " SNe in file " << config.data_file << endl;
}



/*
 * Load JLA data according to information provided in the provided
 * configuration file (typically data/jla.dataset)
 */
/*
void JLALikelihood::read(const char * datasetfile)
{
  Configuration config;
  if (verbosity > 1)
    cout << "Reading config from " << datasetfile << endl;
  if (ini_parse(datasetfile, configuration_handler, &config) < 0)
    {
      if (verbosity > 0) 
        cerr << "Can't load '" << datasetfile << "'\n";
      exit(-1);
    }
  scriptmcut = config.scriptmcut;
  if (verbosity > 1)
    cout << "Config loaded from '" << datasetfile << "': scriptmcut=" << config.scriptmcut << endl;

  C00 = read_matrix(config.C00, verbosity);
  C11 = read_matrix(config.C11, verbosity);
  C22 = read_matrix(config.C22, verbosity);
  C01 = read_matrix(config.C01, verbosity);
  C02 = read_matrix(config.C02, verbosity);
  C12 = read_matrix(config.C12, verbosity);

  char buffer[1024];

  if (verbosity > 1)
    cout << "Reading light-curve parameters from '" << config.data_file << endl;
  ifstream fid(config.data_file);
  if (!fid){
    if (verbosity > 0)
      cerr << "Can't load '" << config.data_file << "'\n";
    exit(-1);
  }
    
  while(fid.get() == '#')
    fid.getline(buffer, 1024);
  fid.unget();
  while (fid)
    {
      LCPar SN;
      fid >> SN;
      if (fid)
        lcpars.push_back(SN);
    }
   
  if (verbosity > 1)
    cout << "Read " << lcpars.size() << " SNe in file " << config.data_file << endl;
}
*/

// ----------------------- Computations ---------------------------------


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
double JLALikelihood::computeLikelihood(double * distance_modulii,
                                        double * nuisance_parameters)
{
  double residuals[size()];
  int status;
  status = computeResiduals(distance_modulii, nuisance_parameters, residuals);
  
  double chi2 = 0;
  if (status == 0)
    {
      for (int i = 0; i < size(); i++){
        chi2 += sq(residuals[i]);
      }
    }
  else
    chi2 = NAN;
  if (verbosity > 2)
    cout << "JLA likelihood evaluation: " << chi2 << endl;
  return chi2;
}


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
int JLALikelihood::computeResiduals(double * distance_modulii, double * nuisance_parameters, double * residuals)
{
  int n = lcpars.size();

  double alpha = nuisance_parameters[0];
  double beta = nuisance_parameters[1];
  double M = nuisance_parameters[2];
  double DeltaM = nuisance_parameters[3];

  // Covariance matrix computation
  double * cov = new double[sq(n)];
  cblas_dcopy(sq(n), C00, 1, cov, 1);
  cblas_daxpy(sq(n), sq(alpha), C11, 1, cov, 1);
  cblas_daxpy(sq(n), sq(beta), C22, 1, cov, 1);
  cblas_daxpy(sq(n), 2.0 * alpha, C01, 1, cov, 1);
  cblas_daxpy(sq(n), -2.0 * beta, C02, 1, cov, 1);
  cblas_daxpy(sq(n), -2.0 * alpha * beta, C12, 1, cov, 1);

  
  for (int i = 0; i < n; i++)
    {
      LCPar sn = lcpars[i];
      // Compute distance modulus estimate
      residuals[i] = sn.mb - (M - alpha * sn.x1 + beta * sn.color + DeltaM * (sn.thirdvar > scriptmcut));
      // Compute residuals
      residuals[i] -= distance_modulii[i];
      // Update the diagonal terms of the covariance matrix with statistical error
      cov[i * n + i] += sq(sn.dmb) + sq(alpha * sn.dx1) + sq(beta * sn.dcolor)
        + 2.0 * alpha * sn.cov_m_s
        - 2.0 * beta * sn.cov_m_c
        - 2.0 * alpha * beta * sn.cov_s_c;
    }

  
  // Whiten the residuals
  int nhrs = 1, info = 0;
  dpotrf("U", &n, cov, &n, &info);
  
  if (info != 0){
    if (verbosity > 0)
      cerr << "Cholesky Error: " << info << endl;
    delete[] cov;
    return -1;
  }
  dtrtrs("U", "T", "N", &n, &nhrs, cov, &n, residuals, &n, &info);  
  if (info != 0){
    if (verbosity > 0)
      cerr << "Solve Error: " << info << endl;
    delete[] cov;
    return -1;
  }

  delete[] cov;
  return 0;
}


// -------------------------- Misc  ---------------------------------------

/*
 * Return the redshift of all SN in the sample.
 */
double * JLALikelihood::getZ()
{
  double * z = (double *)malloc(sizeof(double) * size());
  for (int i=0; i < size(); i++)
    {
      z[i] = lcpars[i].zcmb;
    }
  return z;
}

/*
 * Return the number of SNe
 */
int JLALikelihood::size()
{
  return lcpars.size();
}

/*
 * Free the allocated memory
 */
JLALikelihood::~JLALikelihood()
{
  if (C00 != NULL)
    {
      free(C00);
      free(C11);
      free(C22);
      free(C01);
      free(C02);
      free(C12);
    }
}

