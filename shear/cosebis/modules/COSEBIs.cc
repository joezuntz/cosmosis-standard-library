#include "COSEBIs.h"

COSEBIs::COSEBIs(){}

COSEBIs::~COSEBIs(){}

COSEBIs::COSEBIs(int nMaximum,number thetamin, number thetamax, int nPairs
	,string WnFolderName,string TnFolderName,string OutputTnFolderName,string WnFileName,int precision)
{
	initialize(nMaximum,thetamin, thetamax, nPairs, WnFolderName, TnFolderName,OutputTnFolderName,WnFileName,precision);
}

void COSEBIs::initialize(int nMaximum,number thetamin, number thetamax, int nPairs
	,string WnFolderName1,string TnFolderName1,string OutputTnFolderName1
	,string WnFileName1,int precision1)
{
	clog<<"COSEBIs initialized"<<endl;
	DEn_calculated=false;
	Cov_on=false;
	BCov_on=false;
	WnSet=false;
	TnSet=false;
    EnInteglimitSet=false;
	TpmNotDone=true;
	Cov_Input=false;
	Inner_integ=false;
	WnFolderName=WnFolderName1;
	WnFileName = WnFileName1;
	TnFolderName=TnFolderName1;
	OutputTnFolderName=OutputTnFolderName1;
	WnFileName = WnFileName1;
	precision  = precision1;


	setEparam(nMaximum,thetamin*arcmin,thetamax*arcmin);
	setZbins(nPairs);

}

void COSEBIs::setZbins(int nPairs1)
{
	nPairs=nPairs1;
	nBins=int((-1+sqrt(1+8*nPairs))/2);
}

void COSEBIs::setEparam(int nMaximum1,number thetamin1, number thetamax1)
{
	clog<<"Setting the En parameters in COSEBIs"<<endl;
	DEn_calculated=false;
	Cov_on=false;
	WnSet=false;
	TnSet=false;
	TpmNotDone=true;
    EnInteglimitSet=false;
	nMaximum=nMaximum1;
	clog<<"nMaximum="<<nMaximum<<endl;
	thetamin=thetamin1;
	clog<<"In radians: theta_min="<<thetamin<<endl;
	thetamax=thetamax1;
	clog<<"In radians: theta_max="<<thetamax<<endl;
	LHIGH=high*20./thetamax;
}

void COSEBIs::setWns(int nMaximum, int Nlbins)
{
	if(!WnSet)
	{
        EnInteglimitSet=false;
		clog<<"Setting Wn now:"<<endl;
		WnLog WnL(thetamin,thetamax,nMaximum,TnFolderName,WnFolderName,WnFileName,precision,Nlbins);
		Wn_vec.clear();
		//these have to be separate otherwise it doesn't work
		for(int n=0; n<nMaximum; n++)
			Wn_vec.push_back(WnL);
		for(int n=0; n<nMaximum; n++)
			Wn_vec[n].set(n+1);
		WnSet=true;
	}
}

void COSEBIs::setTs(int nMaximum)
{
	if(!TnSet)
	{
		clog<<"setting Tpm in COSEBIs"<<endl;
		TpnLog Tpn(thetamin,thetamax,nMaximum,TnFolderName,OutputTnFolderName);
		//these have to be separate otherwsie it doesn't work
		for(int n=0; n<nMaximum; n++)
			Tpn_vec.push_back(Tpn);
		for(int n=0; n<nMaximum; n++)
		{
			Tpn_vec[n].PrepareLookupTables(n+1);

		}
		TnSet=true;
	}
}



void COSEBIs::setNoise(number A1,vector<number> sigma_e_vec1,vector<number> nBar_vec)
{
	clog<<"setting noise in COSEBIs"<<endl;
	A=A1;
	noise_vec.clear();
	sigma_e_vec=sigma_e_vec1;
	if(sigma_e_vec.size()==nBar_vec.size())
	{
		for(int bin=0; bin<sigma_e_vec.size(); bin++)
		{
			noise_vec.push_back(sigma_e_vec[bin]*sigma_e_vec[bin]/(2.*nBar_vec[bin]));
			clog<<"noise_vec["<<bin<<"]="<<noise_vec[bin]<<endl;
		}
	}
	else
	{
		clog<<"sigma_e_vec.size()="<<sigma_e_vec.size()<<" != nBar_vec.size()=="<<nBar_vec.size()<<endl;
		clog<<"You need to give the same number of simga_e and ngal values, equal to the number of tomographic bins."<<endl;
		clog<<"Exiting now..."<<endl;
		exit(1);
	}
}


void COSEBIs::setNoiseToZero()
{
	clog<<"setting noise to zero in COSEBIs, noise_vec.size()="<<noise_vec.size()<<endl;
	for(int bin=0; bin<noise_vec.size(); bin++)
	{
		noise_vec[bin]=0.;
		clog<<"noise_vec["<<bin<<"]="<<noise_vec[bin]<<endl;
	}
}

// void COSEBIs::setNoise_vec(vector<number> noise_vec1)
// {
// 	noise_vec=noise_vec1;
// 	clog<<"noise_vec is read: "<<noise_vec.size()<<endl;
// 	// exit(1);
// 	for(int bin=0; bin<noise_vec.size(); bin++)
// 		clog<<"noise_vec["<<bin<<"]="<<noise_vec[bin]<<endl;
// 	clog<<"test"<<endl;
// }

void COSEBIs::setKsi(vector<number> theta,vector<vector<number> > InputKsiP,vector<vector<number> > InputKsiM)
{
	clog<<"setting Ksi in COSEBIs"<<endl;
	Ksi_m_vec.clear();
	Ksi_p_vec.clear();
	if(InputKsiM.size()== nPairs && InputKsiP.size()== nPairs)
	{
		clog<<"Don't panic, number of redshift bin pairs matches the InputKsi vectors"<<endl;
	}
	else
	{
		clog<<"Panic! number of redshift bin pairs DOES NOT match the InputKsi vectors, exiting now ..."<<endl;
		exit(1);
	}
	//these have to be separate otherwise it doesn't work
	for(int r=0;r<nPairs;r++)
	{
		Ksi_p_vec.push_back(function_cosebis());
		Ksi_m_vec.push_back(function_cosebis());
	}
	for(int r=0;r<nPairs;r++)
	{
		Ksi_p_vec[r].loadWithValues(theta,InputKsiP[r],true);
		Ksi_m_vec[r].loadWithValues(theta,InputKsiM[r],true);
		Ksi_p_vec[r].extrapolationOff();
		Ksi_m_vec[r].extrapolationOff();
	}
}

void COSEBIs::setPower(vector<number> log_ell,vector<vector<number> > InputPower)
{
	if(InputPower.size()== nPairs)
	{
		//clog<<"Don't panic, number of redshift bin pairs matches the Input Power vector"<<endl;
	}
	else
	{
		clog<<"Panic! number of redshift bin pairs DOES NOT match the Input Power vector, exiting now ..."<<endl;
		exit(1);
	}
	
	powerspectrum_vec.clear();
	
	//these have to be separate otherwise does not work
	for(int r=0; r<InputPower.size(); r++)
		powerspectrum_vec.push_back(function_cosebis());
	for(int r=0; r<InputPower.size(); r++)
	{
		powerspectrum_vec[r].loadWithValues(log_ell,InputPower[r],true);
		powerspectrum_vec[r].extrapolationOff();
	}
}


void COSEBIs::setPower_single(vector<number> log_x, vector<number> Input)
{
	//clog<<"in set input, nPairs=";
	//nPairs=Input.size();
	//clog<<nPairs<<endl;
	//clog<<"log_x.size()="<<log_x.size()<<" Input.size()="<<Input.size()<<endl;
	nPairs=1;
	powerspectrum_vec.clear();
	powerspectrum_vec.push_back(function_cosebis());
	//clog<<"log_x.size()="<<log_x.size()<<" Input.size()="<<Input.size()<<endl;
	powerspectrum_vec[0].loadWithValues(log_x,Input,true);
	powerspectrum_vec[0].extrapolationOff();
}


void COSEBIs::setPower_single_withExtrapolation(vector<number> log_x, vector<number> Input)
{
	//clog<<"in set input, nPairs=";
	//nPairs=Input.size();
	//clog<<nPairs<<endl;
	//clog<<"log_x.size()="<<log_x.size()<<" Input.size()="<<Input.size()<<endl;
	nPairs=1;
	powerspectrum_vec.clear();
	powerspectrum_vec.push_back(function_cosebis());
	//clog<<"log_x.size()="<<log_x.size()<<" Input.size()="<<Input.size()<<endl;
	powerspectrum_vec[0].loadWithValues(log_x,Input,true);
	powerspectrum_vec[0].extrapolationOn();
}



number COSEBIs::ReturnPower(number ell,int rPair)
{
	return powerspectrum_vec[rPair].value(ell);
}


number COSEBIs::integrant(number l)
{
	if(BCov_on)
	{
		number integ=l*Wn_vec[nW].value(l)*Wn_vec[mW].value(l);
		return integ;
	}
	else if(Cov_Input)
	{
		//calculate covariance for COSEBIs from input Cl cov
		if(Inner_integ)
		{
			number integ=l*Wn_vec[nW].value(l)*CovIn_table_vecvecvec[rp1][rp2][ell]->value(l);
			return integ;
		}
		else
		{
			number integ=l*Wn_vec[mW].value(l)*IntegInner_table.value(l);
			return integ;
		}
	}
	else if(Cov_on)
	{
		counter++;
		number integ;
		number power1=powerspectrum_vec[rp1].value(l)+delta1noise;
		number power2=powerspectrum_vec[rp2].value(l)+delta2noise;
		number power3=powerspectrum_vec[rp3].value(l)+delta3noise;
		number power4=powerspectrum_vec[rp4].value(l)+delta4noise;
		number powers=power1*power2+power3*power4;
		integ=l*Wn_vec[nW].value(l)*Wn_vec[mW].value(l)*powers;
		return integ;
	}
	else if(Real)
	{
		if (realPlus)
		{
			number theta=l;
			number Tp=Tpn_vec[nT].TpValue(theta);
			number Ksip=Ksi_p_vec[redshiftPair].value(theta);
			number integ=theta*Tp*Ksip;
// 			clog<<"theta="<<l<<"  Ksip="<<Ksip<<endl;
			return integ;
		}
		else
		{
			number theta=l;
			number Tm=Tpn_vec[nT].TnValue(theta);
 			number integ=theta*Tm*Ksi_m_vec[redshiftPair].value(theta);
// 			clog<<"theta="<<l<<"  Ksim="<<Ksi_m_vec[redshiftPair].value(theta)<<endl;
			return integ;
		}
	}
	else
	{
		number power=powerspectrum_vec[redshiftPair].value(l);
// 		clog<<"power="<<power<<endl;
		number Wn=Wn_vec[nW].value(l);
// 		clog<<"Wn="<<Wn<<endl;
		number integ=power*Wn*l;
// 		clog<<"integ="<<integ<<endl;
// 		cout<<l<<"    "<<power<<"   "<<Wn<<endl;
		return integ;
	}

}



matrix COSEBIs::returnIntegrandForPowerCase(vector <number> ell_vec, int n, int redshift)
{
	matrix integ_mat(2, ell_vec.size());
	BCov_on=false;
	Cov_on=false;
	Real=false;
	realPlus=false;
	setWns(nMaximum);
	nW=n-1;
	redshiftPair=redshift;
	for(int i=0; i<ell_vec.size();i++)
	{
		integ_mat.load(0,i,ell_vec[i]);
		integ_mat.load(1,i,integrant(ell_vec[i]));
	}
	return integ_mat;
}

number COSEBIs::delta(int i1, int j1)
{
	return i1==j1? 1.: 0.;
}

void COSEBIs::determine_integration_limits_En()
{
    if(!EnInteglimitSet)
    {
        const int Nbins = 1000000;
        // free old list
        integ_limits_vec.clear();
        for(int n=0; n<nMaximum; n++)
        {
            integ_limits.clear();
            // make table of integrant values (Wn's only) on a very fine grid
            matrix table(2,Nbins);
            lthresh =  2.*pi/thetamax*(n+1.)/12.;
     		
            for(int i=0;i<Nbins;i++)
            {
                table.load(0,i,exp(log(lthresh)+log(LHIGH/lthresh)/(Nbins-1.)*i));
                table.load(1,i,Wn_vec[n].value(table.get(0,i)));
            }
            integ_limits.push_back(lthresh);
            
            for(int i=1;i<Nbins-1;i++)
                if ((table.get(1,i-1)<table.get(1,i) && table.get(1,i+1)<table.get(1,i))
                    || (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
                integ_limits.push_back(table.get(0,i));
            integ_limits.push_back(LHIGH);
            integ_limits_vec.push_back(integ_limits);
        }
        EnInteglimitSet=true;
    }
}


number COSEBIs::valueEn(int n)
{
	Cov_on=false;
	Real=false;
	lthresh =  2.*pi/thetamax*n/12.;
	number result= 0.;
	nW=n-1;
	if(LLOW<lthresh)
		result= gaussianIntegrate_gsl(*this,LLOW,lthresh,100);	

	for(unsigned int i=0; (i+1)<integ_limits_vec[nW].size(); i++)
	{
		number res=gaussianIntegrate_gsl(*this,integ_limits_vec[nW][i],integ_limits_vec[nW][i+1],20);
		result+=res;
	}

	return result/2./pi;
}


matrix COSEBIs::calEn()
{
	matrix En(nMaximum*nPairs);
	setWns(nMaximum);
    if(!EnInteglimitSet)
        determine_integration_limits_En();

	for(int r=0; r<nPairs; r++)
	{
		redshiftPair=r;
		for(int n1=nMaximum*r,m=1 ;n1<nMaximum*(r+1) ;n1++,m++)
		{
				En.load(n1,valueEn(m));
		}
	}
	return En;
}

matrix COSEBIs::calEn(vector<int> index)
{
	matrix En(index.size());
	setWns(nMaximum);
    if(!EnInteglimitSet)
        determine_integration_limits_En();

	redshiftPair=0;
	for(int b=0; b<index.size(); b++)
	{
		En.load(b,valueEn(index[b]));
	}

	return En;
}

void COSEBIs::determine_integration_limits_En_Tp()
{
	const int Nbins = 1000000;
	// free old list
	integ_limitsTp.clear();
	// make table of integrant values (Wn's only) on a very fine grid
	matrix table(2,Nbins);
	for(int i=0;i<Nbins;i++)
	{
		table.load(0,i,exp(log(thetamin)+log(thetamax/thetamin)/(Nbins-1.)*i));
		table.load(1,i,Tpn_vec[nT].TpValue(table.get(0,i)));
	}
	integ_limitsTp.push_back(thetamin);
	for(int i=1;i<Nbins-1;i++)
	{
		if ((table.get(1,i-1)<table.get(1,i) && table.get(1,i+1)<table.get(1,i))
			|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
		{
			integ_limitsTp.push_back(table.get(0,i));
		}
	}
	integ_limitsTp.push_back(thetamax);
}

/// min max finder
void COSEBIs::determine_integration_limits_En_Tm()
{
	const int Nbins = 100000;
	// free old list
	integ_limitsTm.clear();
	// make table of integrant values (Wn's only) on a very fine grid
	matrix table(2,Nbins);
	for(int i=0;i<Nbins;i++)
	{
		table.load(0,i,exp(log(thetamin)+log(thetamax/thetamin)/(Nbins-1.)*i));
		table.load(1,i,Tpn_vec[nT].TnValue(table.get(0,i)));
	}
	integ_limitsTm.push_back(thetamin);
	for(int i=1;i<Nbins-1;i++)//2->1
		if ((table.get(1,i-1)<table.get(1,i) && table.get(1,i+1)<table.get(1,i))
			|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
		integ_limitsTm.push_back(table.get(0,i));
	integ_limitsTm.push_back(thetamax);
}


matrix COSEBIs::valueEn2PCFs(int n)
{
	Cov_on=false;
	Real=true;
	number IntPlus= 0.;
	number IntMinus= 0.;
	nT=n-1;
	determine_integration_limits_En_Tp();
	determine_integration_limits_En_Tm();
	for(unsigned int i=0;(i+1)<integ_limitsTp.size();i++)
	{
		realPlus=true;
		number IntP=gaussianIntegrate_gsl(*this,integ_limitsTp[i],integ_limitsTp[i+1],10);	
		IntPlus+=IntP;		
	}
	for(unsigned int i=0;(i+1)<integ_limitsTm.size();i++)
	{
		realPlus=false;
		number IntM=gaussianIntegrate_gsl(*this,integ_limitsTm[i],integ_limitsTm[i+1],10);
		IntMinus+=IntM;
	}
	matrix E3(3);
	E3.load(0,(IntPlus+IntMinus)/2.);
	E3.load(1,IntPlus);
	E3.load(2,IntMinus);
	return E3;
}


matrix COSEBIs::calEn2PCFs()
{
	matrix En(5,nMaximum*nPairs);
	setTs(nMaximum);
	for(int bin1=0;bin1<nBins;bin1++)
	{
		for(int bin2=bin1;bin2<nBins;bin2++)
		{
			int m=0;
			int p1=calP(nBins,bin1,bin2);
			redshiftPair=p1;
			for(int n1=nMaximum*p1,m=1 ;n1<nMaximum*(p1+1) ;n1++,m++)
			{
				clog<<"loading E n="<<n1+1<<endl;
				matrix E3=valueEn2PCFs(m);
				En.load(0,n1,m);
				///E3: En, the plus integral, the minus integral
				En.load(1,n1,E3.get(0));
				En.load(2,n1,E3.get(1));
				En.load(3,n1,E3.get(2));
 				En.load(4,n1,0.5*(E3.get(1)-E3.get(2)));
			}
		}
	}
	return En;
}

matrix COSEBIs::valueEn2PCFsKsiInput(matrix& Ksi_mat,int m,number h)
{
	number IntPlus= 0.;
	number IntMinus= 0.;
	int nTheta=Ksi_mat.rows;
	number intp,intm;
	intp=(Ksi_mat.get(0,0)*Ksi_mat.get(1,0)*Tpm_mat_vec[m].get(0,0))/2.;
	IntPlus=intp;
	intp=(Ksi_mat.get(0,nTheta-1)*Ksi_mat.get(1,nTheta-1)*Tpm_mat_vec[m].get(0,nTheta-1))/2.;
	IntPlus+=intp;
	intm=(Ksi_mat.get(0,0)*Ksi_mat.get(2,0)*Tpm_mat_vec[m].get(1,0))/2.;
	IntMinus=intm;
	IntMinus+=(Ksi_mat.get(0,nTheta-1)*Ksi_mat.get(2,nTheta-1)*Tpm_mat_vec[m].get(1,nTheta-1))/2.;
	for(int itheta=1; itheta<(nTheta-1); itheta++)
	{
		intp=Ksi_mat.get(0,itheta)*Ksi_mat.get(1,itheta)*Tpm_mat_vec[m].get(0,itheta);
		IntPlus +=intp;
		intm=Ksi_mat.get(0,itheta)*Ksi_mat.get(2,itheta)*Tpm_mat_vec[m].get(1,itheta);
		IntMinus+=intm;
	}
	
	IntPlus*=h;
	IntMinus*=h;
	number IntTotal=(IntPlus+IntMinus)/2.;
	
	matrix E3(3);
	E3.load(0,IntTotal);
 	E3.load(1,IntPlus);
 	E3.load(2,IntMinus);
	return E3;
}

vector<vector<number> > COSEBIs::readKsi(string FileName)
{
	int nRow=0;
	string str;
	number temp;
	vector<vector<number> > Ksi_vecvec;
	ifstream KsiFile((FileName).c_str());
  	clog<<"reading, "<<FileName<<endl;
	if(KsiFile.fail())
	{
		clog<<"error occured during opening: "<<FileName<<endl;
		exit(1);
	}

	string line;
	while(getline(KsiFile,line))
	{
		stringstream stream(line);
		str=stream.peek();
		if (str=="#")
		{
			//nColumns++;
// 			string line;
// 			getline(KsiFile,line);
// 			clog<<line<<endl;
		}
		else
		{
			int cols_temp=0;
      		number temp;
      		vector<number> Ksi_vec;
      		while(stream>>temp)
      		{
        		Ksi_vec.push_back(temp);
        		cols_temp++;
      		}
      		if(cols_temp>0)
      		{
				Ksi_vecvec.push_back(Ksi_vec);
				nRow++;
			}
		}
	}

	KsiFile.close();
	return Ksi_vecvec;
}



string COSEBIs::check_binning(vector<number> theta_vec)
{

	linear_binning=true;
	log_binning=true;
	//tolorance for the difference between the theta bins
	number tolorance=0.001;

	clog<<"checking binning, theta_vec.size()="<<theta_vec.size()<<endl;

	number deltaTheta=theta_vec[1]-theta_vec[0];
	number deltaLogTheta=log(theta_vec[1])-log(theta_vec[0]);
	clog<<"deltaTheta="<<deltaTheta<<endl;
	clog<<"deltaLogTheta="<<deltaLogTheta<<endl;

	for(int i=0;i<theta_vec.size()-1; i++)
	{
		number deltaTheta_temp=theta_vec[i+1]-theta_vec[i];
		//clog<<"deltaTheta_temp="<<deltaTheta_temp<<endl;
		if(abs(deltaTheta-deltaTheta_temp)>tolorance)
		{
			linear_binning=false;
		}
	}

	for(int i=0;i<theta_vec.size()-1; i++)
	{
		number deltaLogTheta_temp=log(theta_vec[i+1])-log(theta_vec[i]);
		//clog<<"deltaLogTheta_temp="<<deltaLogTheta_temp<<endl;
		if(abs(deltaLogTheta-deltaLogTheta_temp)>tolorance)
		{
			log_binning=false;
			//clog<<"log binning is set to false"<<endl;
		}
	}

	if(linear_binning)
	{
		clog<<"linear binning for Npair"<<endl;
		return "linear";
	}
	else if(log_binning)
	{
		clog<<"log binning for Npair"<<endl;
		return "log";
	}
	else
	{
		clog<<"binning is neither log or linear"<<endl;
		return "unknown";
	}

}

vector<int> COSEBIs::FindMinMaxTheta(vector<vector<number> > Ksi_vec,int nCol_theta)
{
	int row=0;
	vector<int> index(2);
	clog.precision(10);
	vector<number> theta_vec;
	clog<<"Ksi_vec.size="<<Ksi_vec.size()<<endl;
	for(int i=0;i<Ksi_vec.size();i++)
	{
		theta_vec.push_back(Ksi_vec[i][nCol_theta]);
	}
	string binning=check_binning(theta_vec);
	if(binning=="linear")
	{
		number h=Ksi_vec[1][nCol_theta]-Ksi_vec[0][nCol_theta];

		number tmin=thetamin/arcmin;
		number tmax=thetamax/arcmin;
		if(Ksi_vec[0][nCol_theta]-h>tmin || Ksi_vec[Ksi_vec.size()-1][nCol_theta]+h<tmax)
		{
			clog<<"the theta range is not sufficient:"
				<<Ksi_vec[0][nCol_theta]-h<<"   "<<Ksi_vec[Ksi_vec.size()-1][nCol_theta]+h<<endl;
			exit(1);
		}

		int ntmin=ceil((tmin-Ksi_vec[0][nCol_theta])/h);
		int ntmax=Ksi_vec.size()-1;

		while(tmin>Ksi_vec[ntmin][nCol_theta])
			ntmin++;

		while(tmax<Ksi_vec[ntmax][nCol_theta])
			ntmax--;

		index[0]=ntmin;
		index[1]=ntmax;
		return index;
	}
	else if(binning=="log")
	{
		number h=log(Ksi_vec[1][nCol_theta])-log(Ksi_vec[0][nCol_theta]);

		number tmin=thetamin/arcmin;
		number tmax=thetamax/arcmin;
		if(log(Ksi_vec[0][nCol_theta])-h>log(tmin) || log(Ksi_vec[Ksi_vec.size()-1][nCol_theta])+h<log(tmax))
		{
			clog<<"the theta range is not sufficient:"
				<<Ksi_vec[0][nCol_theta]-exp(h)<<"   "<<Ksi_vec[Ksi_vec.size()-1][nCol_theta]+exp(h)<<endl;
			exit(1);
		}

		int ntmin=ceil((log(tmin)-log(Ksi_vec[0][nCol_theta]))/h);
		int ntmax=Ksi_vec.size()-1;

		while(tmin>Ksi_vec[ntmin][nCol_theta])
			ntmin++;

		while(tmax<Ksi_vec[ntmax][nCol_theta])
			ntmax--;

		index[0]=ntmin;
		index[1]=ntmax;
		return index;
	}
	else
	{
		clog<<"unknown binning"<<endl;
		exit(1);
	}
}


void COSEBIs::FindTnKsiInput(matrix theta_mat)
{

	if(TpmNotDone)
	{

		Tpm_mat_vec.clear();
		int nTheta=theta_mat.size();

		for(int n=0; n<nMaximum; n++)
		{
			matrix Tpm_mat(2,nTheta);
			for(int t=0; t<nTheta; t++)
			{
				number theta=theta_mat.get(t);
				number Tp=Tpn_vec[n].TpValue(theta);
				number Tm=Tpn_vec[n].TnValue(theta);

				Tpm_mat.load(0,t,Tp);
				Tpm_mat.load(1,t,Tm);
			}

			Tpm_mat_vec.push_back(Tpm_mat);
		}
	}
	TpmNotDone=false;
}


matrix COSEBIs::calEn2PCFsFromInputKsi(vector<string> FileName,int nColumns)
{
 	clog<<"in calE2PCFsFromInputKsi in COSEBIs"<<endl;

 	//number of redshift pairs
	nPairs=FileName.size();
	clog<<"nPairs="<<nPairs<<endl;
	nBins=(sqrt(8*nPairs+1)-1)/2;
	matrix En(5,nMaximum*nPairs);
	setTs(nMaximum);
	for(int r=0;r<nPairs;r++)
	{
		vector<vector<number> > Ksi_vec=readKsi(FileName[r]);
		
		vector<int> index=FindMinMaxTheta(Ksi_vec);
		matrix theta_mat(index[1]-index[0]+1);
		matrix Ksi_mat(3,index[1]-index[0]+1);
		for(int i=index[0];i<=index[1]; i++)
		{
			theta_mat.load(i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(1,i-index[0],Ksi_vec[i][1]);
			Ksi_mat.load(2,i-index[0],Ksi_vec[i][2]);
		}
		number h=Ksi_mat.get(0,1)-Ksi_mat.get(0,0);
		FindTnKsiInput(theta_mat);
		for(int n1=nMaximum*r,m=0 ;n1<nMaximum*(r+1) ;n1++,m++)
		{
			matrix E3=valueEn2PCFsKsiInput(Ksi_mat,m,h);
			En.load(0,n1,m+1);
			En.load(1,n1,E3.get(0));
			En.load(2,n1,E3.get(1));
			En.load(3,n1,E3.get(2));
			En.load(4,n1,(E3.get(1)-E3.get(2))/2.);
		}
	}
	return En;
}

matrix COSEBIs::calEn2PCFsFromInputKsi(vector<string> FileName,vector<string> corrFile,int nColumns)
{
	clog<<"in calE2PCFsFromInputKsi in COSEBIs with corrections"<<endl;
	for(unsigned int i=0; i<FileName.size(); i++)
	{
		clog<<"FileName is:"<<FileName[i]<<endl;
		clog<<"Correction FileName is:"<<corrFile[i]<<endl;
	}


	nPairs=FileName.size();
	nBins=(sqrt(8*nPairs+1)-1)/2;
	matrix En(5,nMaximum*nPairs);
	setTs(nMaximum);
	for(int r=0;r<nPairs;r++)
	{
		clog<<"r="<<r<<endl;
		vector<vector<number> > Ksi_vec=readKsi(FileName[r]);
		vector<vector<number> > Corr_vec=readKsi(corrFile[r]);
		
		vector<int> index=FindMinMaxTheta(Ksi_vec);
		matrix theta_mat(index[1]-index[0]+1);
		matrix Ksi_mat(3,index[1]-index[0]+1);
		for(int i=index[0];i<=index[1]; i++)
		{
			theta_mat.load(i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
 			Ksi_mat.load(1,i-index[0],Ksi_vec[i][1]/Corr_vec[i][1]);//xi_+
 			Ksi_mat.load(2,i-index[0],Ksi_vec[i][2]/Corr_vec[i][1]);//xi_-
		}
		number h=Ksi_mat.get(0,1)-Ksi_mat.get(0,0);

		FindTnKsiInput(theta_mat);
		for(int n1=nMaximum*r,m=0 ;n1<nMaximum*(r+1) ;n1++,m++)
		{
			matrix E3=valueEn2PCFsKsiInput(Ksi_mat,m,h);
			En.load(0,n1,m+1);
			En.load(1,n1,E3.get(0));
			En.load(2,n1,E3.get(1));
			En.load(3,n1,E3.get(2));
			En.load(4,n1,(E3.get(1)-E3.get(2))/2.);
		}
	}
	return En;
}

matrix COSEBIs::calEn2PCFsFromInputKsi(vector<string> FileName,vector<number> Corr_vec,int nColumns)
{
	clog<<"in calE2PCFsFromInputKsi in COSEBIs with correction_vector"<<endl;
	for(unsigned int i=0; i<FileName.size(); i++)
	{
		clog<<"FileName is:"<<FileName[i]<<endl;
	}
	nPairs=FileName.size();
	nBins=(sqrt(8*nPairs+1)-1)/2;
	matrix En(5,nMaximum*nPairs);
	setTs(nMaximum);
	for(int r=0;r<nPairs;r++)
	{
		clog<<"r="<<r<<endl;
		vector<vector<number> > Ksi_vec=readKsi(FileName[r]);
		
		vector<int> index=FindMinMaxTheta(Ksi_vec);
		matrix theta_mat(index[1]-index[0]+1);
		matrix Ksi_mat(3,index[1]-index[0]+1);
		for(int i=index[0];i<=index[1]; i++)
		{
			theta_mat.load(i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(0,i-index[0],Ksi_vec[i][0]*arcmin);
			Ksi_mat.load(1,i-index[0],Ksi_vec[i][1]/Corr_vec[r]);
			Ksi_mat.load(2,i-index[0],Ksi_vec[i][2]/Corr_vec[r]);

		}
		number h=Ksi_mat.get(0,1)-Ksi_mat.get(0,0);
		FindTnKsiInput(theta_mat);
		for(int n1=nMaximum*r,m=0 ;n1<nMaximum*(r+1) ;n1++,m++)
		{
			matrix E3=valueEn2PCFsKsiInput(Ksi_mat,m,h);
			En.load(0,n1,m+1);
			En.load(1,n1,E3.get(0));
			En.load(2,n1,E3.get(1));
			En.load(3,n1,E3.get(2));
			En.load(4,n1,(E3.get(1)-E3.get(2))/2.);
		}
	}
	return En;
}


//NOTE: the value of LHIGH is important for the off-diagonals. 
// For better precision use a bigger high
void COSEBIs::determine_integration_limitsCov()
{
/* Idea: find a possibly complete list of consecutive local 
minima/maxima of oscillating integrant and integrate between them
*/
	const int Nbins = 1000000;
	// free old list
	integ_limits.clear();
	// make table of integrant values (Wn's only) on a very fine grid
	matrix table(2,Nbins);
	number lthresh=pi/thetamax/2.;
	for(int i=0;i<Nbins;i++)
	{
		table.load(0,i,exp(log(lthresh)+log(LHIGH/lthresh)/(Nbins-1.)*i));
		table.load(1,i,Wn_vec[nW].value(table.get(0,i))*Wn_vec[mW].value(table.get(0,i)));
	}
// go through list and pick minima/maxima (sort-of; does not need to be awfully exact)
	integ_limits.push_back(lthresh);
	for(int i=1;i<Nbins-1;i++)//2->1
		if ((table.get(1,i-1)<table.get(1,i)&& table.get(1,i+1)<table.get(1,i))
		|| (table.get(1,i-1)>table.get(1,i)&& table.get(1,i+1)>table.get(1,i)))
		      integ_limits.push_back(table.get(0,i));
	integ_limits.push_back(LHIGH);
}

number COSEBIs::valueCov(int n1,int m1)
{
	Cov_on=true;
	nW=n1-1;
	mW=m1-1;
	number lthresh=pi/thetamax/2.;
	determine_integration_limitsCov();
	number result= gaussianIntegrate_gsl(*this,LLOW,lthresh,20);	
	for(unsigned int i=0;(i+1)<integ_limits.size();i++)
	{
		number res=gaussianIntegrate_gsl(*this,integ_limits[i],integ_limits[i+1],20);
		result+=res;
	}
	Cov_on=false;
	return result/A/2./pi;
}

//lthresh integration doesn't need to be separated, this was done as a test. The results change very little.
number COSEBIs::valueBCov(int n1,int m1)
{
	BCov_on=true;
	nW=n1-1;
	mW=m1-1;
	number lthresh=pi/thetamax/2.;
	determine_integration_limitsCov();
	number result= gaussianIntegrate_gsl(*this,LLOW,lthresh,20);	
	for(unsigned int i=0;(i+1)<integ_limits.size();i++)
	{
		number res=gaussianIntegrate_gsl(*this,integ_limits[i],integ_limits[i+1],20);
		result+=res;
	}
	BCov_on=false;
	return result/A/2./pi;
}

//check here if the input is linear or log binned
void COSEBIs::readNpairs(vector<string> FileName_vec,int nCol_theta,int nCol_nPair,bool Athena1)
{
	clog<<"in readNpairs in COSEBIs"<<endl;
	clog<<"FileName_vec.size is:"<<FileName_vec.size()<<endl;
	clog<<"nCol_theta="<<nCol_theta<<", nCol_nPair="<<nCol_nPair<<endl;
	nPairs=FileName_vec.size();
	Npair_mat_vec.clear();

	linear_binning= true;
	log_binning=true;
	Athena=Athena1;
	//tolorance for the difference between the theta bins
	number tolorance=0.001;
	for(int r=0; r<nPairs; r++)
	{
		vector<vector<number> > Ksi_vec=readKsi(FileName_vec[r]);
		vector<int> index=FindMinMaxTheta(Ksi_vec,nCol_theta);
		matrix Npair_mat(2,index[1]-index[0]+1);
		number deltaTheta=Ksi_vec[1][nCol_theta]*arcmin-Ksi_vec[0][nCol_theta]*arcmin;
		number deltaLogTheta=log(Ksi_vec[1][nCol_theta]*arcmin)-log(Ksi_vec[0][nCol_theta]*arcmin);
		//clog<<"deltaTheta="<<deltaTheta<<endl;
		//clog<<"deltaLogTheta="<<deltaLogTheta<<endl;
		for(int i=index[0];i<=index[1]; i++)
		{

			Npair_mat.load(0,i-index[0],Ksi_vec[i][nCol_theta]*arcmin);
			Npair_mat.load(1,i-index[0],Ksi_vec[i][nCol_nPair]);
		}

		for(int i=index[0];i<index[1]; i++)
		{
			number deltaTheta_temp=Ksi_vec[i+1][nCol_theta]*arcmin-Ksi_vec[i][nCol_theta]*arcmin;
			if(abs(deltaTheta-deltaTheta_temp)>tolorance)
			{
				linear_binning=false;
			}
		}

		for(int i=index[0];i<index[1]; i++)
		{
			number deltaLogTheta_temp=log(Ksi_vec[i+1][nCol_theta]*arcmin)-log(Ksi_vec[i][nCol_theta]*arcmin);
			//clog<<"deltaLogTheta_temp="<<deltaLogTheta_temp<<endl;
			if(abs(deltaLogTheta-deltaLogTheta_temp)>tolorance)
			{
				log_binning=false;
				//clog<<"log binning is set to false"<<endl;
			}
		}

		Npair_mat_vec.push_back(Npair_mat);
		clog<<"theta="<<Npair_mat.get(0,0)/arcmin<<"   nPairs="<<Npair_mat.get(1,0)<<endl;

		if(linear_binning)
		{
			clog<<"linear binning for Npair"<<endl;
		}
		else if(log_binning)
		{
			clog<<"log binning for Npair"<<endl;
		}
		else
		{
			clog<<"binning is neither log or linear"<<endl;
		}
	}
}


number COSEBIs::valueNoiseCov_fromInputNpair(int n1,int m1, int np,int bin1,int bin2)
{
	number result=0.;
	int nT=n1-1;
	int mT=m1-1;

	if(linear_binning)
	{
		number deltaTheta=(Npair_mat_vec[np].get(0,2)-Npair_mat_vec[np].get(0,1));

		for(int i=0; i<Npair_mat_vec[np].rows; i++)
		{
			number theta=Npair_mat_vec[np].get(0,i);
			number Npairs=Npair_mat_vec[np].get(1,i);
			number Tp_m=Tpn_vec[mT].TpValue(theta);
			number Tm_m=Tpn_vec[mT].TnValue(theta);
			number Tp_n=Tpn_vec[nT].TpValue(theta);
			number Tm_n=Tpn_vec[nT].TnValue(theta);
			number Npairs_th=2.*deltaTheta*theta*pi*A*neff_vec[bin1]*neff_vec[bin2];
			
			number res=theta*theta*(Tp_m*Tp_n+Tm_m*Tm_n)/Npairs;
			
	// 		number res_B=theta*theta*(Tp_m*Tp_n+Tm_m*Tm_n)/Npair_B;
			// 		number res_th=theta*theta*(Tp_m*Tp_n+Tm_m*Tm_n)/Npairs_th;
			//    		cout<<theta<<"   "<<Npairs<<"    "<<Npairs_th<<endl;
			// 		if( (i % 1000)==0. )
			// 		{
			// 			clog<<"theta="<<theta<<"  Npairs="<<Npairs<<",  Npairs_th="<<Npairs_th<<endl;
			// 			clog<<"res="<<res<<endl;
			// 		}
			//  		result+=res;
			result+=res;
			//  		result+=res_th;
		}

		return result*deltaTheta*deltaTheta;
	}
	else if(log_binning)
	{
		number deltaLogTheta=log(Npair_mat_vec[np].get(0,2))-log(Npair_mat_vec[np].get(0,1));
		//clog<<"deltaLogTheta="<<endl;
		deltaLogTheta=log(Npair_mat_vec[np].get(0,10))-log(Npair_mat_vec[np].get(0,9));

		for(int i=0; i<Npair_mat_vec[np].rows; i++)
		{
			number theta=Npair_mat_vec[np].get(0,i);
			number Npairs=Npair_mat_vec[np].get(1,i);
			number Tp_m=Tpn_vec[mT].TpValue(theta);
			number Tm_m=Tpn_vec[mT].TnValue(theta);
			number Tp_n=Tpn_vec[nT].TpValue(theta);
			number Tm_n=Tpn_vec[nT].TnValue(theta);
			number Npairs_th=2.*deltaLogTheta*theta*theta*pi*A*neff_vec[bin1]*neff_vec[bin2];
			
			number res=theta*theta*theta*theta*(Tp_m*Tp_n+Tm_m*Tm_n)/Npairs;
			
	// 		number res_B=theta*theta*(Tp_m*Tp_n+Tm_m*Tm_n)/Npair_B;
			// 		number res_th=theta*theta*(Tp_m*Tp_n+Tm_m*Tm_n)/Npairs_th;
			//    		cout<<theta<<"   "<<Npairs<<"    "<<Npairs_th<<endl;
			// 		if( (i % 1000)==0. )
			// 		{
			// 			clog<<"theta="<<theta<<"  Npairs="<<Npairs<<",  Npairs_th="<<Npairs_th<<endl;
			// 			clog<<"res="<<res<<endl;
			// 		}
			//  		result+=res;
			result+=res;
			//  		result+=res_th;
		}
		return result*deltaLogTheta*deltaLogTheta;
	}
	else
	{
		clog<<"not a recognised option for binning of Npair - returning nan"<<endl;
		clog<<endl;
		return NAN;
	// 	for(int i=0; i<Npair_mat_vec[np].rows; i++)
	// 	{
	// 		number deltaTheta=(Npair_mat_vec[np].get(0,i+1)-Npair_mat_vec[np].get(0,i));

	// 		number theta=Npair_mat_vec[np].get(0,i);
	// 		number Npairs=Npair_mat_vec[np].get(1,i);
	// 		number Tp_m=Tpn_vec[mT].TpValue(theta);
	// 		number Tm_m=Tpn_vec[mT].TnValue(theta);
	// 		number Tp_n=Tpn_vec[nT].TpValue(theta);
	// 		number Tm_n=Tpn_vec[nT].TnValue(theta);
	// 		number Npairs_th=2.*deltaTheta*theta*pi*A*neff_vec[bin1]*neff_vec[bin2];
			
	// 		number res=theta*theta*(Tp_m*Tp_n+Tm_m*Tm_n)/Npairs;
			
	// // 		number res_B=theta*theta*(Tp_m*Tp_n+Tm_m*Tm_n)/Npair_B;
	// 		// 		number res_th=theta*theta*(Tp_m*Tp_n+Tm_m*Tm_n)/Npairs_th;
	// 		//    		cout<<theta<<"   "<<Npairs<<"    "<<Npairs_th<<endl;
	// 		// 		if( (i % 1000)==0. )
	// 		// 		{
	// 		// 			clog<<"theta="<<theta<<"  Npairs="<<Npairs<<",  Npairs_th="<<Npairs_th<<endl;
	// 		// 			clog<<"res="<<res<<endl;
	// 		// 		}
	// 		//  		result+=res;
	// 		result+=res*deltaTheta*deltaTheta;
	// 		//  		result+=res_th;
	// 	}

	// 	return result;
	}
}


matrix COSEBIs::calNoiseCov_fromInputNpair()
{
	clog<<"calculating the noise covariance in COSEBIs from input npair"<<endl;
	matrix CMT(nMaximum*nPairs,nMaximum*nPairs);

	clog<<"Noise Cov file failed"<<endl;
	setTs(nMaximum);
	for(int n=1; n<nMaximum+1; n++)
	{
		for(int m=n; m<nMaximum+1;m++)
		{
			clog.precision(10);

			for(int bin1=0; bin1<nBins; bin1++)
			{
				for(int bin2=bin1; bin2<nBins; bin2++)
				{
					for(int bin3=bin1; bin3<nBins; bin3++)
					{
						for(int bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
						{
							rp1=calP(nBins,bin1,bin3);
							rp2=calP(nBins,bin2,bin4);
							rp3=calP(nBins,bin1,bin4);
							rp4=calP(nBins,bin2,bin3);
							// 							delta1noise=delta(bin1,bin3)*SigmaE_vec[bin1]*SigmaE_vec[bin1];
							// 							delta2noise=delta(bin2,bin4)*SigmaE_vec[bin2]*SigmaE_vec[bin2];
							// 							delta3noise=delta(bin1,bin4)*SigmaE_vec[bin1]*SigmaE_vec[bin1];
							// 							delta4noise=delta(bin2,bin3)*SigmaE_vec[bin2]*SigmaE_vec[bin2];
							// 							cout<<"delta1noise="<<delta1noise<<endl;
							// 							cout<<"delta2noise="<<delta2noise<<endl;
							// 							cout<<"delta3noise="<<delta3noise<<endl;
							// 							cout<<"delta4noise="<<delta4noise<<endl;
							//the pair considered for the Eparam_vec
							int p1=calP(nBins,bin1,bin2);
							int p2=calP(nBins,bin3,bin4);
							int i=nMaximum*p1+n-1;
							int j=nMaximum*p2+m-1;
							number deltas=delta(bin1,bin3)*delta(bin2,bin4)+delta(bin1,bin4)*delta(bin2,bin3);
							///Why this? This is because of the way Npair is accounted for in Athena. They are counted twice when i=j
							if(Athena)
							{
								if(bin1==bin2)
									deltas=deltas*0.5;
							}
							
							number CovNM=0.;
							if (deltas>0)
								CovNM=valueNoiseCov_fromInputNpair(n,m,p1,bin1,bin2);

							number CovNoise=CovNM*deltas*
								sigma_e_vec[bin1]*sigma_e_vec[bin1]*sigma_e_vec[bin1]*sigma_e_vec[bin1]/8.;
							CMT.load(i,j,CovNoise);
							CMT.load(j,i,CMT.get(i,j));
							CMT.load(i+(m-1)-(n-1),j+(n-1)-(m-1),CMT.get(i,j));
							CMT.load(j+(n-1)-(m-1),i+(m-1)-(n-1),CMT.get(i,j));
						}
					}
				}
			}
		}
	}
	return CMT;
}


bool COSEBIs::checkPower()
{
	if(powerspectrum_vec.size())
		return true;
	return false;
}

matrix COSEBIs::calCov()
{
	clog<<"calculating the covariance in COSEBIs"<<endl;
	matrix CMT(nMaximum*nPairs,nMaximum*nPairs);
	setWns(nMaximum);
	//clog<<"calculating the covariance in COSEBIs"<<endl;
	if(checkPower())
	{
		clog<<"power spectra are set, nBins="<<nBins<<endl;
		for(int bin1=0; bin1<nBins; bin1++)
		{
			for(int bin2=bin1; bin2<nBins; bin2++)
			{
				for(int bin3=bin1; bin3<nBins; bin3++)
				{
					for(int bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
					{
						//clog<<"noise_vec="<<noise_vec[bin1]<<"  "<<noise_vec[bin2]<<endl;
						int n,m=0;
						rp1=calP(nBins,bin1,bin3);
						rp2=calP(nBins,bin2,bin4);
						rp3=calP(nBins,bin1,bin4);
						rp4=calP(nBins,bin2,bin3);
						delta1noise=delta(bin1,bin3)*noise_vec[bin1];
						delta2noise=delta(bin2,bin4)*noise_vec[bin2];
						delta3noise=delta(bin1,bin4)*noise_vec[bin1];
						delta4noise=delta(bin2,bin3)*noise_vec[bin2];
						//the pair considered for the Eparam_vec
						int p1=calP(nBins,bin1,bin2);
						int p2=calP(nBins,bin3,bin4);
						for(int i=nMaximum*p1,n=1;i<nMaximum*(p1+1);i++,n++)
						{
							for(int j=nMaximum*p2+n-1,m=n; j<nMaximum*(p2+1); j++,m++)
							{
								CMT.load(i,j,valueCov(n,m));
								CMT.load(j,i,CMT.get(i,j));
								CMT.load(i+(m-1)-(n-1),j+(n-1)-(m-1),CMT.get(i,j));
								CMT.load(j+(n-1)-(m-1),i+(m-1)-(n-1),CMT.get(i,j));
							}
						}
					}
				}
			}
		}
	}
	else
	{
		clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
		clog<<"!!!!!WARNING power spectrum is not set. The results are not reliable!!!!!!!!!!!!!!"<<endl;
		clog<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	}
	return CMT;
}


matrix COSEBIs::calCovForSigma_m(number sigma_m)
{
	clog<<"calculating the covariance for sigma_m in COSEBIs"<<endl;
	matrix CMT(nMaximum*nPairs,nMaximum*nPairs);
	setWns(nMaximum);
	matrix En_mat=calEn();
	for(int i=0;i<nMaximum*(nPairs);i++)
	{
		for(int j=0; j<nMaximum*(nPairs); j++)
		{
			CMT.load(i,j,En_mat.get(i)*En_mat.get(j));
		}
	}

	return CMT*4.*sigma_m*sigma_m;
}


///test this
matrix COSEBIs::calCovForSigma_m_from_m_cov(matrix m_cov)
{
	clog<<"calculating the covariance for sigma_m with an input m_cov in COSEBIs"<<endl;
	matrix CMT(nMaximum*nPairs,nMaximum*nPairs);
	setWns(nMaximum);
	matrix En_mat=calEn();
	// C^ij,kl_nm=(C_m^ik+C_m^jl+C_m^il+C_m^jk)*E_n^ij*E_n^kl
	for(int n=1; n<nMaximum+1; n++)
	{
		for(int m=n; m<nMaximum+1;m++)
		{
			// clog.precision(10);
			for(int bin1=0; bin1<nBins; bin1++)
			{
				for(int bin2=bin1; bin2<nBins; bin2++)
				{
					for(int bin3=bin1; bin3<nBins; bin3++)
					{
						for(int bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
						{
							number Cm_1 = m_cov.get(bin1,bin3);
							number Cm_2 = m_cov.get(bin2,bin4);
							number Cm_3 = m_cov.get(bin1,bin4);
							number Cm_4 = m_cov.get(bin2,bin3);

							//the pair considered for the Eparam_vec
							int p1=calP(nBins,bin1,bin2);
							int p2=calP(nBins,bin3,bin4);
							int i=nMaximum*p1+n-1;
							int j=nMaximum*p2+m-1;
							number CovNM=En_mat.get(i)*En_mat.get(j);
							number Cov_sigma_m=CovNM*(Cm_1+Cm_2+Cm_3+Cm_4);
							CMT.load(i,j,Cov_sigma_m);
							CMT.load(j,i,CMT.get(i,j));
							CMT.load(i+(m-1)-(n-1),j+(n-1)-(m-1),CMT.get(i,j));
							CMT.load(j+(n-1)-(m-1),i+(m-1)-(n-1),CMT.get(i,j));
						}
					}
				}
			}
		}
	}

	return CMT;
}



matrix COSEBIs::calBCov()
{
	clog<<"calculating the B-mode covariance in COSEBIs"<<endl;
	matrix CMT(nMaximum*nPairs,nMaximum*nPairs);
	setWns(nMaximum);
	for(int n=1; n<nMaximum+1; n++)
	{
		for(int m=n; m<nMaximum+1;m++)
		{
			number CovNM=valueBCov(n,m);
			clog.precision(10);
			for(int bin1=0; bin1<nBins; bin1++)
			{
				for(int bin2=bin1; bin2<nBins; bin2++)
				{
					for(int bin3=bin1; bin3<nBins; bin3++)
					{
						for(int bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
						{
							rp1=calP(nBins,bin1,bin3);
							rp2=calP(nBins,bin2,bin4);
							rp3=calP(nBins,bin1,bin4);
							rp4=calP(nBins,bin2,bin3);
							delta1noise=delta(bin1,bin3)*noise_vec[bin1];
							delta2noise=delta(bin2,bin4)*noise_vec[bin2];
							delta3noise=delta(bin1,bin4)*noise_vec[bin1];
							delta4noise=delta(bin2,bin3)*noise_vec[bin2];

							//the pair considered for the Eparam_vec
							int p1=calP(nBins,bin1,bin2);
							int p2=calP(nBins,bin3,bin4);
							int i=nMaximum*p1+n-1;
							int j=nMaximum*p2+m-1;
							number CovB=CovNM*(delta1noise*delta2noise+delta3noise*delta4noise);
							CMT.load(i,j,CovB);
							CMT.load(j,i,CMT.get(i,j));
							CMT.load(i+(m-1)-(n-1),j+(n-1)-(m-1),CMT.get(i,j));
							CMT.load(j+(n-1)-(m-1),i+(m-1)-(n-1),CMT.get(i,j));
						}
					}
				}
			}
		}
	}
	return CMT;
}



///Need to be tested
number COSEBIs::valueCovFromInputPowerCov(int n1,int m1)
{
	//clog<<"In COSEBIs::valueCovFromInputPowerCov"<<endl;
	Cov_Input=true;
	nW=n1-1;
	mW=m1-1;
	// number lthresh=pi/thetamax/2.;
	
	int n_ell_cov=logell_vec_in.size();
	//starting inner integrals
	vector<number> Inner_integ_vec(n_ell_cov);
	///go through the ells and take integral over ell'
	for(ell=0; ell<n_ell_cov;ell++)
	{
		//Inner integral
		//clog<<"In inner integral"<<endl;
		Inner_integ=true;
		number result_Inner=gaussianIntegrate_gsl(*this,LLOW,integ_limits_vec[nW][0],100);
        int nInteg=integ_limits_vec[nW].size();
//         matrix Inner_res(3,nInteg);
//         Inner_res.load(0,0,LLOW);
//         Inner_res.load(1,0,result_Inner);
//         Inner_res.load(2,0,result_Inner);
		//clog<<LLOW<<"\t"<<integ_limits_vec[nW][0]<<"\t"<<result_Inner<<endl;
// 		clog<<"nInteg="<<nInteg<<endl;
		for(unsigned int i=0;(i+1)<nInteg;i++)
		{
           
			number res=gaussianIntegrate_gsl(*this,integ_limits_vec[nW][i],integ_limits_vec[nW][i+1],20);
			result_Inner+=res;
//             Inner_res.load(0,i+1,integ_limits_vec[nW][i]);
//             Inner_res.load(1,i+1,res);
//             Inner_res.load(2,i+1,result_Inner);
			// clog<<"integ_limits_vec["<<nW<<"]["<<i<<"]="<<integ_limits_vec[nW][i]
			// 	<<"  integ_limits_vec["<<nW<<"]["<<i+1<<"]="<<integ_limits_vec[nW][i+1]
			//  	<<"  inner intergal is: "<<result_Inner<<endl;
			//exit(1);
			// cout<<integ_limits_vec[nW][i]<<"\t"<<result_Inner<<"\t"<<res<<endl;
		}
// 		Inner_res.printOut((string("Inner_result_")+toString(ell)+string("ascii")).c_str(),10);
		Inner_integ_vec[ell]=result_Inner;
 		//clog<<"ell="<<ell<<" inner="<<Inner_integ_vec[ell]<<endl;
	}
	///Make a table of the inner integrals to be interpolated for the outer integral
	IntegInner_table.loadWithValues(logell_vec_in,Inner_integ_vec,true);
// 	IntegInner_table.setName((string("IntegTable")).c_str(),function::NONAMECOUNTER);
// 	IntegInner_table.saveTable();
	IntegInner_table.extrapolationOff();
    
// 	int n_ell=10000;
// 	number ellmin=1.0;
// 	number ellmax=1e7;
// 	matrix testInner(2,n_ell);
// 	
// 	for(int i=0; i<n_ell;i++)
// 	{
// 		number elltest=exp(log(ellmin)+log(ellmax/ellmin)/(n_ell)*(i+0.5));
// 		testInner.load(0,i,elltest);
// 		number Inner=IntegInner_table.value(elltest);
// 		testInner.load(1,i,Inner);
// 	}
// 	testInner.printOut("testInner.ascii",10);
    
	//outer integral
	Inner_integ=false;
	number result_Outer=gaussianIntegrate_gsl(*this,LLOW,integ_limits_vec[mW][0],100);
	for(unsigned int i=0;(i+1)<integ_limits_vec[mW].size();i++)
	{
		number res=gaussianIntegrate_gsl(*this,integ_limits_vec[mW][i],integ_limits_vec[mW][i+1],20);
		result_Outer+=res;
	}
	//clog<<" outer="<<result_Outer/2./pi/2./pi<<endl;
	Cov_Input=false;
	return result_Outer/2./pi/2./pi;
}


void COSEBIs::setTableCovIn(matrix& InputCov)
{
	//clog<<"in COSEBIs::setTableCovIn"<<endl;

	int n_ell_cov=logell_vec_in.size();
    
	CovIn_table_vecvecvec.clear();
    for(int np1=0; np1<nPairs; np1++)
    {
       	vector<vector<function_cosebis*> > CovIn_table_vecvec;
        for(int np2=0; np2<nPairs; np2++)
        {
            vector<function_cosebis*> CovIn_table_vec;
        	for(int i=n_ell_cov*np1,iell=0;i<n_ell_cov*(np1+1);i++, iell++)
        	{
        		//each row for each Covariance block with n_ell entries
       			CovIn_table_vec.push_back(new function_cosebis());

                vector<number> InputCov_rows(n_ell_cov);
                for(int j=n_ell_cov*np2, jell=0; j<n_ell_cov*(np2+1); j++,jell++)
                {
                    InputCov_rows[jell]=InputCov.get(i,j);
                }
                CovIn_table_vec[iell]->loadWithValues(logell_vec_in,InputCov_rows,true);   
                CovIn_table_vec[iell]->extrapolationOff();
                // CovIn_table_vec[iell].setName((string("CovIn_table")+toString(i)).c_str(),function_cosebis::NONAMECOUNTER);
                // CovIn_table_vec[iell].saveTable(); 
            }
            CovIn_table_vecvec.push_back(CovIn_table_vec);
        }
        CovIn_table_vecvecvec.push_back(CovIn_table_vecvec);
    }
	clog<<"I've set the Input covariance tables for interpolation and integration"<<endl;
}

matrix COSEBIs::calCovFromInputPowerCov(matrix& InputCov,vector<number> ell_vec_in)
{
	logell_vec_in.clear();
    clog.precision(5);
	for(int i=0; i<ell_vec_in.size(); i++)
	{
		logell_vec_in.push_back(log(ell_vec_in[i]));
		//clog<<"Logell["<<i<<"]="<<logell_vec_in[i]<<endl;
	}
	clog<<"calculating the covariance from input Cl covariance in COSEBIs"<<endl;
	matrix CMT(nMaximum*nPairs,nMaximum*nPairs);

	setWns(nMaximum);
	///just for Wn not to get the Inner integral first and then the outer integral
	///For the Gaussian terms the integration limits is based on Wn(ell)*Wm(ell)
	determine_integration_limits_En();
	//This takes the input Cov and devides it into block and rows to prepare for the inner integrals
	setTableCovIn(InputCov);
	clog<<"nBins="<<nBins<<endl;

	for(int bin1=0; bin1<nBins; bin1++)
	{
		for(int bin2=bin1; bin2<nBins; bin2++)
		{
			for(int bin3=bin1; bin3<nBins; bin3++)
			{
				for(int bin4=(bin3==bin1) ? bin2:bin3; bin4<nBins;bin4++)
				{
					//clog<<"Calculating Covariance now"<<endl;
					int n,m=0;
					// rp3=calP(nBins,bin1,bin4);
					// rp4=calP(nBins,bin2,bin3);
					//the pair considered for the Eparam_vec
					int p1=calP(nBins,bin1,bin2);
					rp1=p1;
                    redshiftPair=p1;
					int p2=calP(nBins,bin3,bin4);
					rp2=p2;
					for(int i=nMaximum*p1,n=1;i<nMaximum*(p1+1);i++,n++)
					{
						for(int j=nMaximum*p2+n-1,m=n; j<nMaximum*(p2+1); j++,m++)
						{
							//clog<<bin1+1<<bin2+1<<" "<<bin3+1<<bin4+1<<endl;
// 								bool writeInteg=false;
							CMT.load(i,j,valueCovFromInputPowerCov(n,m));
							//CMT.load(i,j,(bin1+1)*10000+(bin2+1)*1000+(bin3+1)*10+bin4+1);
							CMT.load(j,i,CMT.get(i,j));
							CMT.load(i+(m-1)-(n-1),j+(n-1)-(m-1),CMT.get(i,j));
							CMT.load(j+(n-1)-(m-1),i+(m-1)-(n-1),CMT.get(i,j));
						}
					}
				}
			}
		}
	}

	return CMT;
}

///////////////////////////covariance functions ended//////////////////////////////////////////


int COSEBIs::calP(int nBins,int fbin,int sbin)
{
	if(fbin>sbin)
	{
		int swap=fbin;
		fbin=sbin;
		sbin=swap;
	}
	int p=fbin*nBins;

	for(int i=0; i<fbin; i++)
	  	p-=i;
	return p+sbin-fbin;
}


