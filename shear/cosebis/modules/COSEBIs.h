
//Written by Marika Asgari

///this class calculates En and Bn and their covariance matrix 
//It needs a convergence power spectrum or two point correlation functions as input
//for the covariance matrix it needs the noise value for the Gaussian term 
//and a non-Gaussian Cl covariance if you wish to calculate the non-Gaussian terms. 


#ifndef COSEBIS_H
#define COSEBIS_H

//This one calculates logarithmic Wn functions which are used in making En from theory:
//En=int_0^inf dl l/(2pi) Wn(ell) P_E(ell)
//Bn=int_0^inf dl l/(2pi) Wn(ell) P_B(ell)
#include "WnLog.h"
//This one calculates T_+- which are used to calculate COSEBIs from 2PCFs
//En=int_{theta_min}^{theta_max} d theta theta/2[T_+(theta)xi_+(theta)+T_-(theta)xi_-(theta)]
//Bn=int_{theta_min}^{theta_max} d theta theta/2[T_+(theta)xi_+(theta)-T_-(theta)xi_-(theta)]
#include "TpnLogRadian.h"
//This makes matrices and does matrix stuff (most of the functions are disables because they used NR)
#include "matrix.h"
#include "errors.h"

/*


CAUTION: The covariance calculation with the number of pairs depends 
on the definition of the number of pairs. For example in Athena_1.7 for the autocorrelations
when both gal_cats are given a factor of 0.5 is needed to correct the number of pairs. 
This is in the COSEBIs code now but this convention might be different if using different codes.
*/

class COSEBIs : public function_cosebis
{
public:
	//Constructors
	COSEBIs();
	COSEBIs(int nMax,number thetamin, number thetamax, int nPairs,
		string WnFolderName,string TnFolderName,string OutputTnFolderName,string WnFileName = "WnLog",int precision = 20);
	//destructor
	~COSEBIs();
	//Constructor calls initialize
	void initialize(int nMaximum,number thetamin, number thetamax, int nPairs
		,string WnFolderName1=COSEBIS_DIR "/WnLog/"
		,string TnFolderName1=COSEBIS_DIR "/TLogsRootsAndNorms/"
		,string OutputTnFolderName="./cosebis/TpnLog/"
		,string WnFileName = "WnLog"
		,int precision = 20);
	//sets the number of redshift bins
	void setZbins(int nPairs1);
	///sets the COSEBIs parameters
	void setEparam(int nMax1,number thetamin1, number thetamax1);
	///initializes Wn_vec
	void setWns(int nMax, int Nlbins = 1000000);
	///initializes Tp_n_vec and Tm_n_vec
	void setTs(int nMaximum);
	///sets the noise parameters for calculating the covariance 
	void setNoise(number A1,vector<number> sigma_e_vec,vector<number> nBar_vec);
	///sets neff_vec
	void set_neff(vector<number> neff_vec);
	///sets noise value for each bin separately 
	void setNoise_vec(vector<number> noise_vec1);
	///sets the already exsiting noise_vec to zero.
	void setNoiseToZero();
	///calculates power for a given redshift bin combination 
	void setPower(vector<number> log_ell_vec,vector<vector<number> > InputPower);
	void setPower_single(vector<number> log_x, vector<number> Input);
	void setPower_single_withExtrapolation(vector<number> log_x, vector<number> Input);
	///Return power
	number ReturnPower(number ell,int rPair);
	///sets 2PCFs from input KsiP and KsiM. Makes tables of them for interpolation. NOTE: the Ksi+- here need to
	///be from theory and smooth
	void setKsi(vector<number> theta,vector<vector<number> > InputKsiP,vector<vector<number> > InputKsiM);
	///integrand for En or Cov depends if Cov_on is true or not
	number integrant(number l);
	///return the integrand for the power spectrum case
	matrix returnIntegrandForPowerCase(vector <number> ell_vec, int n, int redshift);
	///min max finder for En setWns uses this and determines the integration limits once and for all
	void determine_integration_limits_En();
	///min max finder for the plus part of the En from 2PCFs integral
	void determine_integration_limits_En_Tp();
	///min max finder for the minus part of the En from 2PCFs integral
	void determine_integration_limits_En_Tm();
	///the integrand used to find En from the Trapezoidal integration rutine. if noisyKsi is true
	///noise will be added to Ksi_+ and Ksi_- before integration
	number trapIntegrant(number theta);
	///finds the value E_n from 2PCFs
	matrix valueEn2PCFs(int n);
	///calculates En for all n from 2PCFs
	matrix calEn2PCFs();
	///calculates En for a given value of n
	number valueEn(int n);
	///calculate En frm P(ell)
	matrix calEn();
	//only works for single binning
	matrix calEn(vector<int> index);
	//calculates the covariance for En Em, n=n1 m=m1
	number valueCov(int n1,int m1);
	//checks if input power soectrum is set. Returns false if not.
	bool checkPower();
	///calculates Covariance matrix and returns it
	matrix calCov();
	//takes in the value of sigma_m and returns a covariance that is the equal to
	// 4*sigma_m^2* E^{ij}_m * E^{kl}_n= C^{ij kl}_mn
	matrix calCovForSigma_m(number sigma_m);
	// given an input matrix with the same number of rows and columns as nPairs, calculates the covariance for m_bias uncertainity
	matrix calCovForSigma_m_from_m_cov(matrix m_cov);
	///reads Npairs from an input ksi file from Athena
	//added a check linear/log binning
	void readNpairs(vector<string> FileName_vec,int nCol_theta=0,int nCol_nPair=7,bool Athena1=true);
	///this calculates the Noise covariance assuming noise is Gaussian, it returns
	///\sigma_e^4/ 4* (d \theta)^2 \sum_i \theta_i^2/ N(\theta_i) [T_+m(\theta_i) T_+n(\theta_i)+ T_-m(\theta_i) T_-n(\theta_i)]
	number valueNoiseCov_fromInputNpair(int n1,int m1, int np,int bin1,int bin2);
	///From input napir
	matrix calNoiseCov_fromInputNpair();
	///this calculates the B-mode covariance assuming that there is no P_B(l) and noise is Gaussian, it returns
	///1/(2piA)*\int dl l W_n W_m
	number valueBCov(int n1,int m1);
	matrix calBCov();
	///returns 1 if i==j, zero otherwise
	number delta(int i, int j);
	///Determines the min max values for integration for the Covariance
	void determine_integration_limitsCov();
	///takes n for En and the number of bins for the trapezoidal integral and returns the En values
	///using 2PCFs
	//matrix valueEn2PCFsTrap(int n,int nBinsTrap,number h);
	//matrix calEn2PCFsTrap(int nBinsTrap,bool noisyKsi1);
	///calculates the En from the input Ksi_mat with trapezoidal integration resturns a matrix:
	///0:(Int_p+Int_m)/2 1:Int_p 2:Int_m
	matrix valueEn2PCFsKsiInput(matrix& Ksi_mat,int m,number h);
	///reads a Ksi file with theta ksi+ and ksi- returns the values in a 2D vector 
	vector<vector<number> > readKsi(string FileName);
	//checks what type of binning is used for theta
	string check_binning(vector<number> theta_vec);
	///finds the minTheta and maxTheta which are the closest to thetamin and thetamax
	///return the index of the MinTheta=index[0] and MaxTheta=index[1]
	vector<int> FindMinMaxTheta(vector<vector<number> > ksi_vec,int nCol_theta=0);
	///evaluates Tpm for the given theta's in theta_mat and saves them in Tpm_mat_vec, does this only once
	void FindTnKsiInput(matrix theta_mat);
	///calculates En from an input Ksi file using trapezoidal integration
	matrix calEn2PCFsFromInputKsi(vector<string> FileName,int nColumns=8);	
	///calculates En from an input Ksi file using trapezoidal integration
	matrix calEn2PCFsFromInputKsi(vector<string> FileName, vector<string> corrFile,int nColumns=8);
	///calculates En from an input Ksi file using trapezoidal integration
	matrix calEn2PCFsFromInputKsi(vector<string> FileName, vector<number> Corr_vec,int nColumns=8);
	///sets the value of the index from the bin
	int calP(int nBins,int fbin,int sbin);
	///calculates the covariance matrix for COSEBIs given an input covariance matrix for power spectra
	///mainly used for the non-Gaussian terms: connected from trispectrum and the super sample covariance (SSC)
	number valueCovFromInputPowerCov(int n1,int m1);
	void setTableCovIn(matrix& InputCov);
	matrix calCovFromInputPowerCov(matrix &InputCov,vector<number> ell_vec_in1);


private: //private variables
  string WnFolderName,TnFolderName,OutputTnFolderName,WnFileName;
  vector<WnLog> Wn_vec;
  vector<TpnLog> Tpn_vec; 
  vector<function_cosebis> powerspectrum_vec;
  vector<function_cosebis> Ksi_p_vec;
  vector<function_cosebis> Ksi_m_vec;

  function_cosebis IntegInner_table;


  vector <number> integ_limits,integ_limitsTp,integ_limitsTm,logell_vec_in;
  vector <vector<number> > integ_limits_vec;
  vector<vector<vector<function_cosebis*> > > CovIn_table_vecvecvec;
  vector <number> noise_vec,sigma_e_vec;
  vector <matrix> Ksi_mat_vec;
  vector <matrix> Tpm_mat_vec,pofz_mat_vec;
  vector <matrix> Npair_mat_vec;
  vector <number> neff_vec;
  
  number delta1noise,delta2noise,delta3noise,delta4noise;  
  number A,sigmaE,begin,end,nBar;
  number thetamin,thetamax,lthresh,LHIGH,NoiseKsi;
  
  int param,nPairs,powerswitch,nBins,nMaximum,derivative;
  int redshiftPair,rp1,rp2,rp3,rp4;
  int nW,mW,counter,nT,iRand;  
  int ell; 
  int precision;

  ///default is false
  bool TpmNotDone,noisyKsi, BCov_on,DEn_calculated,
		OneParam,Cov_on,Real,realPlus,WnSet,TnSet,EnInteglimitSet,
		Cov_Input,Inner_integ,
		linear_binning,log_binning,
		Athena;

  matrix En_data, Cov_mat;
};

#endif
