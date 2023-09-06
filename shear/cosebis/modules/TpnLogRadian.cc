#include "TpnLogRadian.h"

TpnLog::TpnLog()
{
	thetamin=thetamax=1.; n=0;
}

TpnLog::TpnLog(number thetamin1,number thetamax1,int nMax,string TnFolderName,string foldername)
{
	setTheta(thetamin1,thetamax1,nMax,TnFolderName,foldername);
}

TpnLog::~TpnLog(){}

number TpnLog::integrant(number y)
{
	number integ=4*TpLog(y)*(exp(2.*(y-zInput))-3.*exp(4.*(y-zInput)));
	if(writeInteg)
		cout<<y<<"\t"<<integ<<endl;
	return integ;
}

void TpnLog::setTheta(number thetamin1,number thetamax1,int nMax,string TnFolderName,string foldername1)
{
	clog<<"setting thetamin,thetamax,nMax"<<endl;
	clog<<"nMax="<<nMax<<endl;
	foldername=foldername1;
	lookupTableDone=false;
	thetamin=thetamin1;
	thetamax=thetamax1;
	root.clear();
	norm.clear();
	rootTheta.clear();
	ifstream fRoot((TnFolderName+string("/Root_")+toString(thetamin/arcmin,2)+string("-")
	    +toString(thetamax/arcmin,2)+string(".table")).c_str());
	ifstream fNorm((TnFolderName+string("/Normalization_")+toString(thetamin/arcmin,2)+string("-")
	    +toString(thetamax/arcmin,2)+string(".table")).c_str());
	if(fRoot.fail())
	{
		clog<<"error occured during opening Root file"<<endl;
		exit(1);
	}

	if(fNorm.fail())
	{
		clog<<"error occured during opening Normalization file"<<endl;
		exit(1);
	}

	for(n=1; n<=nMax; n++)
	{
		if(fRoot.eof())
		{
			clog<<"nMax is too big"<<endl;
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

int TpnLog::show_n()
{
	return n;
}

number TpnLog::TpLog(number z)
{
	//number z=log(theta/thetamin);
	number result=1.;

	for(int i=1; i<=(n+1); i++)
		result*=(z-root[n-1][i-1]);
	result*=norm[n-1];
	return result;
}

void TpnLog::writeTpLog(int n1)
{
	n=n1;
	matrix Tpn_mat(nThetaBins);
	matrix theta_mat(nThetaBins);
	number eps=1e-15;

	clog<<"Tp"<<n<<" Log="<<norm[n-1]<<"  ";
	for(int i=1; i<=(n+1);i++)
		clog<<"(z-"<<root[n-1][i-1]<<")";
	clog<<endl;
	//number arcminlog=log(arcmin);
	number thetaminlog=log(thetamin);
	for(int i=0; i<nThetaBins; i++)
	{
		number thetalog=(thetaminlog-eps
			+log((thetamax+eps)/(thetamin-eps))/(nThetaBins-1)*i);
		theta_mat.load(i,thetalog);
		Tpn_mat.load(i,TpLog(thetalog-thetaminlog));
	}
	string myname =foldername+string("/TpRadian")+toString(n)+string("_")
		+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2);
	Tp_table.setName(myname.c_str(),function_cosebis::NONAMECOUNTER);
	Tp_table.loadWithValues(theta_mat, Tpn_mat,true);
	Tp_table.saveTable();
}

number TpnLog::TnLog(number zIn)
{
	zInput=zIn;
	writeInteg=false;
	number result=0.;
	vector<number> integ_steps;
	integ_steps.push_back(0.);
	int counter=0;
	while(zInput>root[n-1][counter] && counter<root[n-1].size())
	{
		integ_steps.push_back(root[n-1][counter]);
		counter++;
	}
	integ_steps.push_back(zInput);
	//number thetaminlog=log(thetamin);
	cout.precision(5);
	//cout<<"#n="<<n<<" size root[n]="<<root[n-1].size()<<endl;
	//cout<<"#zInput="<<zInput<<"\t"<<zInput/log(thetamax/thetamin)<<endl;
	
	for(int i=0; i<integ_steps.size()-1;i++)
	{
		number res=gaussianIntegrate_gsl(*this, integ_steps[i], integ_steps[i+1],20);
		//cout<<"#res="<<res<<"  integ_steps["<<i<<"]="
			//<<integ_steps[i]<<"  integ_steps["<<i+1<<"]="<<integ_steps[i+1]<<endl;
		result+=res;
	}
	result+=TpLog(zInput);

	return result;
}

void TpnLog::writeTnLog(int n1)
{
	n=n1;
	matrix Tnn_mat(nThetaBins);
	matrix theta_mat(nThetaBins);
	number eps=1e-15;

	//number arcminlog=log(arcmin);
	number thetaminlog=log(thetamin);
	for(int i=0; i<nThetaBins; i++)
	{
		number thetalog=(thetaminlog-eps
			+log((thetamax+eps)/(thetamin-eps))/(nThetaBins-1)*i);
		theta_mat.load(i,(thetalog));
		Tnn_mat.load(i,TnLog(thetalog-thetaminlog));
	}
	string myname =foldername+string("/TmRadian")+toString(n)+string("_")
		+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2);
	Tn_table.setName(myname.c_str(),function_cosebis::NONAMECOUNTER);
	Tn_table.loadWithValues(theta_mat, Tnn_mat,true);
	Tn_table.saveTable();
}

void TpnLog::PrepareLookupTables(int n1)
{
	lookupTableDone=true;
	clog<<"in PrepareLookupTables"<<endl;
	n=n1;
	string myname =foldername+
		string("/TpRadian")+toString(n)+string("_")
		+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2);
	clog<<myname<<endl;
	ifstream fTp((myname+string(".table")).c_str());
	clog<<myname<<endl;
	if(fTp.fail())
	{
		clog<<"Tp file doesn't exist:"<<myname+string(".table")<<endl;
		clog<<"making Tp table"<<endl;
		writeTpLog(n);
	}
	else
	{
		Tp_table.setName(myname.c_str(),function_cosebis::NONAMECOUNTER);
		clog<<"loading Tp"<<endl;
		Tp_table.loadTable(thetamin,thetamax,nThetaBins,true);
	}
	Tp_table.extrapolationOff();
	fTp.close();
	clog<<"Tp loaded"<<endl;
	myname =foldername+
		string("/TmRadian")+toString(n)+string("_")
		+toString(thetamin/arcmin,2)+string("-")+toString(thetamax/arcmin,2);
	ifstream fTn((myname+string(".table")).c_str());

	if(fTn.fail())
	{
		clog<<"Tn file doesn't exist:"<<myname+string(".table")<<endl;
		clog<<"making Tn table"<<endl;
		writeTnLog(n);
	}
	else
	{
		Tn_table.setName(myname.c_str(),function_cosebis::NONAMECOUNTER);
		Tn_table.loadTable(thetamin,thetamax,nThetaBins,true);
	}
	Tn_table.extrapolationOff();
	fTn.close();
}

number TpnLog::TpValue(number theta)
{
	if(lookupTableDone)
		return Tp_table.value(theta);
	else
		clog<<"lookupTable not Done, make it first!"<<endl;
	exit(1);
}

number TpnLog::TnValue(number theta)
{
	if(lookupTableDone)
		return Tn_table.value(theta);
	else
		clog<<"lookupTable not Done, make it first!"<<endl;
	exit(1);
}
