#include <iostream>
#include <tgmath.h>
#include <ctime>

using namespace std;

double dBinomFunc(double dn, double dm);
double dBinomFunc(int nn, int nm);
int tIncreasement(int *pt, int nr, int ns);//set pt to the next value. pt has length nr, sum to ns, monotone decreasing
int tSetValue(int *pt, int nb, int nr, int ns);//length nr, each value<=nb, sum to ns

int kBound(int nn, int nd, int nr);
int kBoundWeakForm(int nn, int nd, int nr);

void printUsage()
{
	cout<<"Usage: "<<endl;
        cout<<"./upperbound n d r"<<endl;
	cout<<"n is the code length, n>=6"<<endl;
	cout<<"d is the minimum distance, d>=5"<<endl;
	cout<<"r is the localtiy parameter, r>=2"<<endl;
}

int parameterCheck(int nn, int nd, int nr)
{
	if(nn<6 || nd<5 || nr<2)
	{
		cout<<"Please input valid code parameters."<<endl;
		return -1;
	}

	if(nn>255 && nd>8)
	{
		cout<<"Warning: the inputed code parameters are large, it may take several hours to compute the bound. Try '-f' mode to accelerate the computation."<<endl;
		return 1;
	}

	if(nn<=nd || nn<=nr)
	{
		cout<<"Please input valid code parameters. The codelength can not be smaller than the distance or the locality parameters"<<endl;
		return -1;
	}

	return 0;
}

int main(int argc, char** argv)
{

	if(argc != 4 && argc != 5)
	{
		printUsage();	
		return 0;
	}

	if(argc == 4)
	{
		int nn = atoi(argv[1]);
		int nd = atoi(argv[2]);
		int nr = atoi(argv[3]);

		switch(parameterCheck(nn,nd,nr))
		{
			case 1:
				cout<<"Computing in fast mode."<<endl;
				break;
			case -1:
				return 0;
				break;
			default:
				cout<<"Computing in fast mode."<<endl;
		}
	
		int kb = kBoundWeakForm(nn,nd,nr);	
		cout<<"The upper bound on the dimension is: k <= "<<kb<<endl;

	}	
	
	if(argc == 5)
	{
		int nn = atoi(argv[2]);
		int nd = atoi(argv[3]);
		int nr = atoi(argv[4]);

		switch(parameterCheck(nn,nd,nr))
		{
			case 1:
				cout<<"Type 'Y' to continue."<<endl;
				if(getchar()!='Y')
					return 0;
				cout<<"Continuing computation."<<endl;
				break;
			case -1:
				return 0;
				break;
			default:
				cout<<"Computing."<<endl;
		}
	
		int kb = kBound(nn,nd,nr);	
		cout<<"The upper bound on the dimension is: k <= "<<kb<<endl;

	}
	
	return 0;
}

double dBinomFunc(double dn, double dm)
{
	if(dn<0 || dn<dm)
		return 0;
	if(dm==0 || dn==dm)
		return 1;
	double dtmp = lgamma(dn+1)-lgamma(dm+1)-lgamma(dn-dm+1);
	dtmp = round(exp(dtmp));
	
	return dtmp;
}

double dBinomFunc(int nn, int nm)
{
	double dn = (double)nn;
	double dm = (double)nm;

	return dBinomFunc(dn, dm);
}

int tIncreasement(int *pt, int nr, int ns)
{
	int index = nr-1;
	int partialSum = 0;
	
	//find the last nonzero coordinates(other than the last)
	while(index>0 && pt[index-1]==0)
		index--;
	//if all coordinates is zero, cannot increase
	if(index==0)
		return 0;

	int setSuccess = 0;
	while(setSuccess==0)
	{
		pt[index-1]--;
		partialSum = 0;
		for(int ni=0;ni<index;ni++)
			partialSum = partialSum+pt[ni];
		setSuccess = tSetValue(pt+index, pt[index-1], nr-index, ns-partialSum);

		if(setSuccess==0)
			index--;
		if(index<1)
			break;
	}

	if(index<1)
		return 0;
	else
		return 1;

}

int tSetValue(int *pt, int nb, int nr, int ns)
{
	if(nb==0 && ns>0)
		return 0;

	if(nb==0 && ns==0)
	{
		for(int ni=0;ni<nr;ni++)
			pt[ni] = 0;
	}
	
	const int ntmp = (int)ceil((1.0*ns)/nb);
	if(ntmp>nr)
		return 0;

	for(int ni=0;ni<ntmp-1;ni++)
		pt[ni] = nb;
	pt[ntmp-1] = ns-(ntmp-1)*nb;
	for(int ni=ntmp;ni<nr;ni++)
		pt[ni] = 0;
	return 1;
}

int kBound(int nn, int nd, int nr)
{
	time_t timer = time(NULL);
	
	const int ntmp1 = (int)ceil((1.0*nn)/(nr+1));
	const int ntmp2 = (int)floor((2.0*nn)/(nr+2));
	const int ntmp3 = (int)floor((1.0*(nd-1))/4);

	double *psi = new double[ntmp2-ntmp1+1];

	for(int nl=ntmp1; nl<=ntmp2; nl++)
		psi[nl-ntmp1] = -1;

	int *nt = new int[nr];//t1+t2+...+tr = 2n-l(r+2); 0<=tr<=...<=t1<=l;
	
	int indexT = ntmp2;

	for(int nl=ntmp1; nl<=ntmp2; nl++)
	{
		timer = time(NULL);

		/****test print****/	
		//cout<<nl<<endl;
		
		int nrsum = 2*nn-nl*(nr+2);
		int *npi = new int[nl];//the index
		int *npr = new int[nl];//the ri

		int nIncreased = tSetValue(nt, nl, nr, nrsum);
		if(nIncreased!=1)
			cerr<<"Initialize error!"<<endl;

		while(nIncreased==1)
		{
			//compute ri
			int ntmpIndex = nr;
			for(int nj=0;nj<nl;nj++)
			{
				while(nj>=nt[ntmpIndex-1] && ntmpIndex>0)
					ntmpIndex--;
				npr[nj] = ntmpIndex;
			}

			//compute psi
			double tmpPsiValue = 0;
			double tmpProd = 1;
			int indexSum = 0;
			int updateIndex = 1;
			for(int nj=0;nj<nl;nj++)
				npi[nj] = 0;
			
			while(indexSum<=ntmp3)
			{
				tmpProd = 1;
				for(int nj=0;nj<nl;nj++)
				{
					if(tmpProd==0)
						break;
					tmpProd = tmpProd * (int)dBinomFunc(npr[nj]+1,2*npi[nj]);
				}	

				tmpPsiValue = tmpPsiValue + tmpProd;

				updateIndex = 1;
				while(updateIndex<=nl)
				{
					npi[updateIndex-1]++;
					indexSum++;
					
					if(indexSum<=ntmp3)
						break;

					indexSum = indexSum - npi[updateIndex-1];
					npi[updateIndex-1] = 0;
					updateIndex++;
				}
				if(updateIndex>nl)
					break;
			}			
			
			if(psi[nl-ntmp1]==-1 || psi[nl-ntmp1]>tmpPsiValue)
				psi[nl-ntmp1] = tmpPsiValue;	

			/****test print****/
			/*
			for(int nj=0;nj<nr;nj++)
				cout<<nt[nj]<<", ";
			cout<<endl;			
			for(int nj=0;nj<nl;nj++)
				cout<<npr[nj]<<", ";
			cout<<endl;
			*/

			//increase the index
			nIncreased = tIncreasement(nt, nr, nrsum);

		}
		
		delete[] npr;
		delete[] npi;

		if(time(NULL)-timer>100)
		{
			indexT = nl;
			break;
		}
	
	}	

	for(int nl=ntmp2; nl>indexT; nl--)
	{
		timer = time(NULL);

		/****test print****/	
		//cout<<nl<<endl;
		
		int nrsum = 2*nn-nl*(nr+2);
		int *npi = new int[nl];//the index
		int *npr = new int[nl];//the ri

		int nIncreased = tSetValue(nt, nl, nr, nrsum);
		if(nIncreased!=1)
			cerr<<"Initialize error!"<<endl;

		while(nIncreased==1)
		{
			//compute ri
			int ntmpIndex = nr;
			for(int nj=0;nj<nl;nj++)
			{
				while(nj>=nt[ntmpIndex-1] && ntmpIndex>0)
					ntmpIndex--;
				npr[nj] = ntmpIndex;
			}

			//compute psi
			double tmpPsiValue = 0;
			double tmpProd = 1;
			int indexSum = 0;
			int updateIndex = 1;
			for(int nj=0;nj<nl;nj++)
				npi[nj] = 0;
			
			while(indexSum<=ntmp3)
			{
				tmpProd = 1;
				for(int nj=0;nj<nl;nj++)
				{
					if(tmpProd==0)
						break;
					tmpProd = tmpProd * (int)dBinomFunc(npr[nj]+1,2*npi[nj]);
				}	

				tmpPsiValue = tmpPsiValue + tmpProd;

				updateIndex = 1;
				while(updateIndex<=nl)
				{
					npi[updateIndex-1]++;
					indexSum++;
					
					if(indexSum<=ntmp3)
						break;

					indexSum = indexSum - npi[updateIndex-1];
					npi[updateIndex-1] = 0;
					updateIndex++;
				}
				if(updateIndex>nl)
					break;
			}			
			
			if(psi[nl-ntmp1]==-1 || psi[nl-ntmp1]>tmpPsiValue)
				psi[nl-ntmp1] = tmpPsiValue;	

			/****test print****/
			/*
			for(int nj=0;nj<nr;nj++)
				cout<<nt[nj]<<", ";
			cout<<endl;			
			for(int nj=0;nj<nl;nj++)
				cout<<npr[nj]<<", ";
			cout<<endl;
			*/

			//increase the index
			nIncreased = tIncreasement(nt, nr, nrsum);

		}
		
		delete[] npr;
		delete[] npi;

		if(time(NULL)-timer>100)
			break;
	
	}	

	int iMin = -1;
	for(int nl=ntmp1; nl<=ntmp2; nl++)
	{
		if(psi[nl-ntmp1]==-1)
			continue;

		if(iMin==-1 || iMin> nl + (int)ceil(log2((double)psi[nl-ntmp1])))
			iMin =  nl + (int)ceil(log2((double)psi[nl-ntmp1]));
	}
	
	//delete[] nt;
	delete[] psi;

	return (nn-iMin);
}

int kBoundWeakForm(int nn, int nd, int nr)
{
	time_t timer = time(NULL);
	
	const int ntmp1 = (int)ceil((1.0*nn)/(nr+1));
	const int ntmp2 = (int)floor((2.0*nn)/(nr+2));
	const int ntmp3 = (int)floor((1.0*(nd-1))/4);

	double *psi = new double[ntmp2-ntmp1+1];

	for(int nl=ntmp1; nl<=ntmp2; nl++)
		psi[nl-ntmp1] = -1;

	int *nt = new int[nr];//t1+t2+...+tr = 2n-l(r+2); 0<=tr<=...<=t1<=l;
	
	int indexT = ntmp2;

	for(int nl=ntmp1; nl<=ntmp2; nl++)
	{
		timer = time(NULL);

		/****test print****/	
		//cout<<nl<<endl;
		
		int nrsum = 2*nn-nl*(nr+2);
		int *npi = new int[nl];//the index
		int *npr = new int[nl];//the ri

		int nIncreased = tSetValue(nt, nl, nr, nrsum);
		if(nIncreased!=1)
			cerr<<"Initialize error!"<<endl;

		while(nIncreased==1)
		{
			//compute ri
			int ntmpIndex = nr;
			for(int nj=0;nj<nl;nj++)
			{
				while(nj>=nt[ntmpIndex-1] && ntmpIndex>0)
					ntmpIndex--;
				npr[nj] = ntmpIndex;
			}

			//compute psi
			double tmpPsiValue = 0;
			double tmpProd = 1;
			int indexSum = 0;
			int updateIndex = 1;
			for(int nj=0;nj<nl;nj++)
				npi[nj] = 0;
			
			while(indexSum<=ntmp3)
			{
				tmpProd = 1;
				for(int nj=0;nj<nl;nj++)
				{
					if(tmpProd==0)
						break;
					tmpProd = tmpProd * (int)dBinomFunc(npr[nj]+1,2*npi[nj]);
				}	

				tmpPsiValue = tmpPsiValue + tmpProd;

				updateIndex = 1;
				while(updateIndex<=nl)
				{
					npi[updateIndex-1]++;
					indexSum++;
					
					if(indexSum<=ntmp3)
						break;

					indexSum = indexSum - npi[updateIndex-1];
					npi[updateIndex-1] = 0;
					updateIndex++;
				}
				if(updateIndex>nl)
					break;
			}			
			
			if(psi[nl-ntmp1]==-1 || psi[nl-ntmp1]>tmpPsiValue)
				psi[nl-ntmp1] = tmpPsiValue;	

			/****test print****/
			/*
			for(int nj=0;nj<nr;nj++)
				cout<<nt[nj]<<", ";
			cout<<endl;			
			for(int nj=0;nj<nl;nj++)
				cout<<npr[nj]<<", ";
			cout<<endl;
			*/

			//increase the index
			nIncreased = tIncreasement(nt, nr, nrsum);

		}
		
		delete[] npr;
		delete[] npi;

		if(time(NULL)-timer>100)
		{
			indexT = nl;
			break;
		}
	
	}	

	for(int nl=ntmp2; nl>indexT; nl--)
	{
		timer = time(NULL);

		/****test print****/	
		//cout<<nl<<endl;
		
		int nrsum = 2*nn-nl*(nr+2);
		int *npi = new int[nl];//the index
		int *npr = new int[nl];//the ri

		int nIncreased = tSetValue(nt, nl, nr, nrsum);
		if(nIncreased!=1)
			cerr<<"Initialize error!"<<endl;

		while(nIncreased==1)
		{
			//compute ri
			int ntmpIndex = nr;
			for(int nj=0;nj<nl;nj++)
			{
				while(nj>=nt[ntmpIndex-1] && ntmpIndex>0)
					ntmpIndex--;
				npr[nj] = ntmpIndex;
			}

			//compute psi
			double tmpPsiValue = 0;
			double tmpProd = 1;
			int indexSum = 0;
			int updateIndex = 1;
			for(int nj=0;nj<nl;nj++)
				npi[nj] = 0;
			
			while(indexSum<=ntmp3)
			{
				tmpProd = 1;
				for(int nj=0;nj<nl;nj++)
				{
					if(tmpProd==0)
						break;
					tmpProd = tmpProd * (int)dBinomFunc(npr[nj]+1,2*npi[nj]);
				}	

				tmpPsiValue = tmpPsiValue + tmpProd;

				updateIndex = 1;
				while(updateIndex<=nl)
				{
					npi[updateIndex-1]++;
					indexSum++;
					
					if(indexSum<=ntmp3)
						break;

					indexSum = indexSum - npi[updateIndex-1];
					npi[updateIndex-1] = 0;
					updateIndex++;
				}
				if(updateIndex>nl)
					break;
			}			
			
			if(psi[nl-ntmp1]==-1 || psi[nl-ntmp1]>tmpPsiValue)
				psi[nl-ntmp1] = tmpPsiValue;	

			/****test print****/
			/*
			for(int nj=0;nj<nr;nj++)
				cout<<nt[nj]<<", ";
			cout<<endl;			
			for(int nj=0;nj<nl;nj++)
				cout<<npr[nj]<<", ";
			cout<<endl;
			*/

			//increase the index
			nIncreased = tIncreasement(nt, nr, nrsum);

		}
		
		delete[] npr;
		delete[] npi;

		if(time(NULL)-timer>100)
			break;
	
	}	

	int iMin = -1;
	for(int nl=ntmp1; nl<=ntmp2; nl++)
	{
		if(psi[nl-ntmp1]==-1)
			continue;

		if(iMin==-1 || iMin> nl + (int)ceil(log2((double)psi[nl-ntmp1])))
			iMin =  nl + (int)ceil(log2((double)psi[nl-ntmp1]));
	}
	
	//delete[] nt;
	delete[] psi;

	return (nn-iMin);
}



