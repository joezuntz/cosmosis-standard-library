#include "WnLog.h"

WnLog::WnLog(): Wn_table()
{
	l=thetamin=thetamax=1.;
	n=0;
	lthresh=0.;
	TnFolderName= COSEBIS_DIR "/TLogsRootsAndNorms";
	WnFolderName= COSEBIS_DIR "/WnLog/";
	WnFileName="WnLog";
}

WnLog::WnLog(number thetamin1,number thetamax1,int nMax,string TnFolderName,string WnFolderName
		,string WnFileName, int precision,int Nlbins)
{

  setWnLogName(TnFolderName,WnFolderName,WnFileName);
  setTheta(thetamin1,thetamax1,nMax);
  setPrecision(precision);
  setNlbins(Nlbins);
}

WnLog::~WnLog(){}

number WnLog::integrant(number x)
{
	number integ=x*TnLog(x/l)*gsl_sf_bessel_J0(x);
	return integ;
}

void WnLog::setWnLogName(string TnFolderName1,string WnFolderName1,string WnFileName1)
{
	TnFolderName=TnFolderName1;
	WnFolderName=WnFolderName1;
	WnFileName=WnFileName1;
}

void WnLog::setPrecision(int precision1)
{
	precision=precision1;
}

void WnLog::setNlbins(int Nlbins1)
{
	Nlbins=Nlbins1;
}

void WnLog::setTheta(number thetamin1,number thetamax1,int nMax)
{
	//clog<<"setting thetamin,thetamax,nMax"<<endl;
	//clog<<"nMax="<<nMax<<endl;
	thetamin=thetamin1;
	thetamax=thetamax1;
	LHIGH=high*20./thetamax;
	root.clear();
	norm.clear();
	rootTheta.clear();
	string fRoot_name=TnFolderName+string("/Root_")+toString(thetamin/arcmin,2)+string("-")
	    +toString(thetamax/arcmin,2)+string(".table");
	string fNorm_name=TnFolderName+string("/Normalization_")+toString(thetamin/arcmin,2)+string("-")
	    +toString(thetamax/arcmin,2)+string(".table");
	ifstream fRoot((fRoot_name).c_str());
	ifstream fNorm((fNorm_name).c_str());
	if(fRoot.fail())
	{
		clog<<"in WnLog.cc: error occured during opening Root file:"<<fRoot_name<<endl;
		clog<<"exiting now ..."<<endl;
		exit(1);
	}


	if(fNorm.fail())
	{
		clog<<"in WnLog.cc: error occured during opening Normalization file:"<<fNorm_name<<endl;
		clog<<"exiting now ..."<<endl;
		exit(1);
	}

	for(n=1; n<=nMax; n++)
	{
		if(fRoot.eof())
		{
			clog<<"nMax is too big, nMax="<<nMax<<", exiting now ..."<<endl;
			exit(1);
		}
		int r=0;
		fRoot >>r;
		clog<<"reading T"<<r<<"  roots"<<"\t";
		vector<number> row; // Create an empty row
		vector<number> rowTheta;
		for (int i =1 ; i <=n+1 ; i++)
		{
			number temp;
			fRoot>>temp;
			row.push_back(temp); // Add an element (column) to the row
			rowTheta.push_back(thetamin*exp(temp));
		}
		root.push_back(row);
		rootTheta.push_back(rowTheta);
		number normal;
		fNorm >> normal;
		clog<<"reading T"<<normal<<"  normalization"<<endl;
		fNorm>>normal;
		norm.push_back(normal);
	}
	fRoot.close();
	fNorm.close();
}

void WnLog::set(int order)
{
	//clog<<"set order of WnLog"<<endl;
	n= order;
	// is there a table on disk?
	string myname =WnFolderName+string("/")+WnFileName+toString(n)+string("-")
		 +toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2);
	setName(myname.c_str(),function_cosebis::NONAMECOUNTER);
	ifstream fhandler((myname+string(".table")).c_str());
	// NO?
	if (fhandler.fail())
	{
		if(!CheckFolderExist(WnFolderName))
		{
			clog<<"Making the folder for Wn:"<<WnFolderName<<endl;
			mkdir((WnFolderName).c_str(), 0777);
		}
		if(!CheckFolderExist(WnFolderName))
		{
			clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			clog<<"!!! Can't make the Folder for Wn functions !!!!!!!"<<endl;
			clog<<"!!!!!!!!!!! WILL NOT SAVE FILES !!!!!!!!!!!!!!!!!"<<endl;
			clog<<"!!!!!!!!!!! YOU'VE BEEN WARNNED !!!!!!!!!!!!!!!!!"<<endl;
			clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		}

		clog<<"writing table:"<<myname<<endl;
		StepFinder();
		// transition point between the different orders for the Gauss-Legendere integration
		lthresh =  2.*pi/thetamax*order/10.;
		//clog<<"lthresh="<<lthresh<<endl;
		//cout << LLOW << " -- " << LHIGH << " -- " << Nlbins << endl;

		clog<<"I am calculating the Wn window function.  This takes a few minutes per node."<<endl;
		clog<<"You only need to calculate each node's table once."<<endl;
		clog<<"It is then stored for future calculations."<<endl;
	}
	fhandler.close();

	// make table of myself and save or load existing one
	//clog<<"LLOW="<<log(LLOW)<<" LHIGH="<<log(LHIGH)<<" Nlbins="<<Nlbins<<endl;
	loadTable(LLOW,LHIGH,Nlbins,true,precision);
	extrapolationOff();
}

int WnLog::show_n()
{
	return n;
}

number WnLog::get(number x)
{
	l = x;
	if (l<lthresh)
		// This is not accurate enough for low ell and W1. 
		// return gaussianIntegrate_gsl(*this,thetamin*l,thetamax*l,10*n)/l/l;
		return gaussianIntegrate_gsl(*this,thetamin*l,thetamax*l,200)/l/l;
	return valueStep(l);
}

number WnLog::TnLog(number theta)
{
	number z=log(theta/thetamin);
	number result=1.;

	for(int i=1; i<=(n+1); i++)
		result*=(z-root[n-1][i-1]);
	result*=norm[n-1];
	return result;
}

void WnLog::writeTnLog(int n1)
{
	n=n1;
	matrix Tn_mat(2,10000);

	clog<<"T"<<n<<" Log="<<norm[n-1]<<"  ";
	for(int i=1; i<=(n+1);i++)
		clog<<"(z-"<<root[n-1][i-1]<<")";
	clog<<endl;
	for(int i=0; i<10000; i++)
	{
		number theta=exp(log(thetamin)+log(thetamax/thetamin)/(10000-1)*i);
		Tn_mat.load(0,i,theta/arcmin);
		Tn_mat.load(1,i,TnLog(theta));
	}
	Tn_mat.printOut((string("TLog_")+toString(n)+string("_")+toString(thetamin/arcmin,2)
		+string("-")+toString(thetamax/arcmin,2)+string(".ascii")).c_str(),8);
}


//integrates between zeros of bessel and filter
number WnLog::valueStep(number l)
{
	number xmin = (l*thetamin);
	number xmax = (l*thetamax);
	number resultG=0.;
	int i=0;
	if (StepN[0]>xmax)
	{
		resultG=gaussianIntegrate_gsl(*this, xmin, xmax,accuracyG);
 		cout<<"# in if (StepN[0]>xmax) resultG="<<resultG<<endl;
	}
	else
	{
		for(i=0; xmin>(StepN[i]); i++);

		number result1=gaussianIntegrate_gsl(*this, xmin, StepN[i],accuracyG);
		resultG+=result1;
		while(StepN[i+1]<xmax)
		{
			number result1=gaussianIntegrate_gsl(*this, StepN[i],StepN[i+1],accuracyG);
			resultG+=result1;
			i++;
		}
		resultG+=gaussianIntegrate_gsl(*this, StepN[i],xmax,accuracyG);
	}
	return resultG/l/l;
}

void WnLog::StepFinder()
{
	//StepN is a vector of consequitive zeros of integrant for each n

	clog<<"StepFinder begins"<<endl;
	int i0=0;
	int iR=0;
	number arg=0.;
	StepN.clear();
	while((arg<LHIGH*thetamax)&&(i0<100000)&&(iR<100000))
	{
		if(iR>n)
		{
			arg=BesselJ0_Zeros[i0];
			i0++;
		}
		else
		{
			if(BesselJ0_Zeros[i0]<rootTheta[n-1][iR])
			{
				arg=BesselJ0_Zeros[i0];
				i0++;
			}
			else
			{
				arg=rootTheta[n-1][iR];
				iR++;
				//clog<<"iR="<<iR<<endl;
			}
		}

		StepN.push_back(arg);

	}
	clog<<"StepFinder ended"<<endl;
	//clog<<"i0="<<i0<<"iR="<<iR<<endl;

}

void WnLog::writeStep(int n1)
{
	n=n1;
	StepFinder();
	for(unsigned int i=0; i<StepN.size(); i++)
		cout<<i<<"\t"<<StepN[i]<<endl;
}

