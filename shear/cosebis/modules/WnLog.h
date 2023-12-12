#ifndef  WNLOG_H
#define  WNLOG_H

#include "function_cosebis.h"
#include "Integrate.h"
#include "besselzero04.h"
#include <gsl/gsl_sf_bessel.h>

/// the lowest ell that is calculated
const number LLOW = 1.;
const number high =500.;
///the highest ell that is calculated
//	LHIGH=high*20./thetamax;

// The number of points in the Gauss-Legendre integration
const number accuracyG=40;

/// number of log-bins for Wn-table
//const int    NLBINS  = 1000000;
// const int    NLBINS  = 10000;

//code by Marika Asgari and Patrick Simon

class WnLog : public function_cosebis
{
public:

	WnLog();
	WnLog(number thetamin1,number thetamax1,int nMax
		,string TnFolderName=COSEBIS_DIR "/TLogsRootsAndNorms/"
		,string WnFolderName=COSEBIS_DIR "/WnLog/"
		,string WnFileName="WnLog"
		,int precision = 20
		,int Nlbins    = 1000000);
	~WnLog();

	///sets the number of digits for the saved table
	void setPrecision(int precision1);
	///sets the number of l bins in Wnlog table
	void setNlbins(int Nlbins1);
	/// integrant for Wn(l), depends on internal parameters (thetamin,thetamax,n)
	number integrant(number x);
	///sets thetamin thetamax and nMax to read the roots and the normalization from disk
	void setTheta(number thetamin1,number thetamax1,int nMax);
	///sets the Folder for WnLog and the root and normalisation and the start of the WnLog file name. 
	void setWnLogName(string TnFolderName,string WnFolderName,string WnFileName);
	/// sets the internal parameters (n) to open Wn file or make it
	void set(int order);
	///returns n 
	int show_n();
	/// returns Wn(l) 
	number get(number x);
	///writes TnLog on a matrix and on disk
	void writeTnLog(int n1);
	///calculates the value of WnLog using an step by step integration between zeros of bessel function and TnLogs
	number valueStep(number l);
	void StepFinder();
	void writeStep(int n1);



private:

	/// TLog+n(x)
	number TnLog(number theta);
	/// internal parameters
	function_cosebis Wn_table;
	vector<number> j0minmax;
	vector< vector<number> > root;
	vector< vector<number> > rootTheta;
	vector <number> norm;
	vector <number> StepN;
	int      n;
	number   l;
	number   lthresh;
	number   thetamin;
	number   thetamax;
	number   B;
	number LHIGH;
	string TnFolderName,WnFolderName,WnFileName;
	int precision,Nlbins;
};

#endif
