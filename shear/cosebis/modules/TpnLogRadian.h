#ifndef TPNLOGRADIAN_H
#define TPNLOGRADIAN_H

#include "function_cosebis.h"
#include "Integrate.h"

const int nThetaBins=10000;

// Code by Marika Asgari and Patrick Simon

class TpnLog : public function_cosebis
{
public:

	TpnLog();
	TpnLog(number thetamin1,number thetamax1,int nMax,string TnFolderName,string foldername="TpnLog");
	~TpnLog();

	/// integrant for T_{-n}(\theta), depends on internal parameters (thetamin,thetamax,n)
	number integrant(number y);
	///sets thetamin thetamax and nMax to read the roots and the normalization from disk
	void setTheta(number thetamin1,number thetamax1,int nMax,string TnFolderName,string foldername="TpnLog");
	///returns n 
	int show_n();
	/// Peter's TLog+n(x) z=theta/thetamin
	number TpLog(number z);
	///writes TnLog_+ on a matrix and on disk
	void writeTpLog(int n1);
	/// TLog-n(x)
	number TnLog(number theta);
	///writes TnLog_- on a matrix and on disk
	void writeTnLog(int n1);
	//number LHIGH;
	///prepares look up tables for Tp and Tn for a given n, if the table exists on the disk it will read from the disk and if not it will write to the disk
	void PrepareLookupTables(int n);
	///returns the value of Tp for a given theta, first must prepare look up table
	number TpValue(number theta);
	///returns the value of Tn for a given theta, first must prepare look up table
	number TnValue(number theta);

private:
	/// determines integration limits (based on max/min of J0(x))
	vector< vector<number> > root;
	vector< vector<number> > rootTheta;
	vector <number> norm;
	int      n;
	number   thetamin;
	number   thetamax;
	number zInput;
	function_cosebis Tp_table;
	function_cosebis Tn_table;
	bool lookupTableDone,writeInteg;
	string foldername;

};

#endif
