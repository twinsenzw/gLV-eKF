#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <algorithm>
#include <random>
#include <time.h>
#include <map>
#include <unordered_map>
#include <thread>
#include <armadillo>

using namespace std;
using namespace arma;

void vnorm(vector<double> &v)
{
	double sum=0;
	for(auto vi : v)
	{
		sum+=vi;
	}
	if(sum!=0)
	{
		for(int i=0; i<v.size(); ++i)
		{
			v[i]=v[i]/sum;
		}
	}
}

//=================The generalized lotka-voltera model=====================

struct gLV
{
	//suppose there are Nv microbes and Lv external stimuli, a total of Tv generations measured

	//parameters

	int Nv;
	int Tv;
	int Lv;

	vector<double> b; //intrinsic birth rate (replication rate) of each taxa, size=N
	vector<double> d; //intrinsic death rate

	vector<vector<double>> Bm; //influence of microbe i on j's birth rate, size=N*N
	vector<vector<double>> Dm; //influence of microbe i on j's death rate

	vector<vector<double>> Bs; //influence of stimuli i on microbe j's birth rate size=L*N
	vector<vector<double>> Ds; 


	//observations

	vector<vector<double>> RA; //relative abundances of microbes at time t; size=T*N
	vector<vector<double>> Stim; //stimuli measurement at time t; size=T*L

	//simulations
	vector<vector<double>> RAsim; //simulated relative abundances, size=T*N


	//methods
	void initialize(int N, int L, int T); //initialize a model

	void inputStim(vector<double> &stim, int t); //input stimuli abundances (size=L) to time t in Stim
	void inputRA(vector<double> &ra, int t); //input relative abundances (size=N) to time t in RA
	void inputRAsim(vector<double> &ra, int t); //input relative abundances (size=N) to time t in RAsim

	void inputInitialCondition(ifstream &input); //input initial conditions to get simulation started

	void simNextStep(int i, int start); //simulate abundance of microbe i at t=start+1, based on RAsim of t=start
	void simNextStepGroup(vector<int> &v, int start); //simulate a group of microbial abundance at t=start+1, one by one
	void simNextStepPara(int p, int start, int normalize); //simulate all microbial abundance at t=start+1, with multithreading

	void sim(int p, int start, int end, int normalize); //simulate all microbial abundance from t=start+1 to t=end, with multithreading

	void printPara(); //print parameters
	void printRAsim(int t); //print simulated results at time t
};

void gLV::initialize(int N, int L, int T)
{
	cout<<"===Initializing gLV...";

	Nv=N;
	Lv=L;
	Tv=T;

	vector<double> tmp1;
	vector<double> tmp2;

	for(int i=0; i<N; i++)
	{
		tmp1.push_back(0);
		d.push_back(0);
		b.push_back(0);
	}
	
	for(int i=0; i<N; i++)
	{
		Bm.push_back(tmp1);
		Dm.push_back(tmp1);
	}
	
	for(int l=0; l<L; l++)
	{
		tmp2.push_back(0);
		Bs.push_back(tmp1);
		Ds.push_back(tmp1);
	}

	for(int t=0; t<T; t++)
	{
		Stim.push_back(tmp2);
		RAsim.push_back(tmp1);
		RA.push_back(tmp1);
	}

	cout<<"Finished."<<endl;
}

void gLV::inputStim(vector<double> &stim, int t)
{
	cout<<"===Inputting stim profiles "<<t<<"...";

	Stim[t]=stim;
	
	cout<<"Finished."<<endl;		
}

void gLV::inputRA(vector<double> &ra, int t)
{
	cout<<"===Inputting relative abundance profiles to time "<<t<<" in observation...";

	RA[t]=ra;
	
	cout<<"Finished."<<endl;		
}

void gLV::inputRAsim(vector<double> &ra, int t)
{
	cout<<"===Inputting relative abundance profiles to time "<<t<<" in simulation...";

	RAsim[t]=ra;
	
	cout<<"Finished."<<endl;		
}	

void gLV::inputInitialCondition(ifstream &input)
{
	string line;
	//N, T, L
	vector<int> NLT;

	for(int cnt=0; (cnt<3)&&(getline(input,line)) ; ++cnt)
	{
		istringstream linestream(line);
		string f;
		linestream >> f >> f;
		NLT.push_back(stod(f));
	}

	initialize(NLT[0],NLT[1],NLT[2]);

	//RA
	getline(input,line);
	istringstream RAstream(line);
	string ra;
	getline(RAstream,ra,'\t');
	vector<double> ras;

	for(int i=0; i<Nv; i++)
	{
		getline(RAstream,ra,'\t');
		ras.push_back(stod(ra));
	}

	inputRAsim(ras, 0);

	//Stim
	getline(input,line);
	istringstream Stimstream(line);
	string st;
	getline(Stimstream,st,'\t');
	
	for(int t=0; t<Tv; t++)
	{	
		vector<double> sts;
		for(int l=0; l<Lv; l++)
		{
			getline(Stimstream,st,'\t');
			sts.push_back(stod(st));
		}
		inputStim(sts,t);
	}

	//b
	getline(input,line);
	istringstream bstream(line);
	string ib;
	getline(bstream,ib,'\t');
	
	for(int i=0; i<Nv; i++)
	{
		getline(bstream,ib,'\t');
		b[i]=stod(ib);
	}

	//d
	getline(input,line);
	istringstream dstream(line);
	string id;
	getline(dstream,id,'\t');
	
	for(int i=0; i<Nv; i++)
	{
		getline(dstream,id,'\t');
		d[i]=stod(id);
	}

	//bm
	getline(input,line);
	istringstream bmstream(line);
	string ibm;
	getline(bmstream,ibm,'\t');
	
	for(int j=0; j<Nv; j++)
	{
		for(int i=0; i<Nv; i++)
		{
			getline(bmstream,ibm,'\t');
			Bm[j][i]=stod(ibm);
		}
	}

	//dm
	getline(input,line);
	istringstream dmstream(line);
	string idm;
	getline(dmstream,idm,'\t');
	
	for(int j=0; j<Nv; j++)
	{
		for(int i=0; i<Nv; i++)
		{
			getline(dmstream,idm,'\t');
			Dm[j][i]=stod(idm);
		}
	}

	//bs
	getline(input,line);
	istringstream bsstream(line);
	string ibs;
	getline(bsstream,ibs,'\t');
	
	for(int l=0; l<Lv; l++)
	{
		for(int i=0; i<Nv; i++)
		{
			getline(bsstream,ibs,'\t');
			Bs[l][i]=stod(ibs);
		}
	}

	//ds
	getline(input,line);
	istringstream dsstream(line);
	string ids;
	getline(dsstream,ids,'\t');
	
	for(int l=0; l<Lv; l++)
	{
		for(int i=0; i<Nv; i++)
		{
			getline(dsstream,ids,'\t');
			Ds[l][i]=stod(ids);
		}
	}
}

void gLV::simNextStep(int i, int start)
{
	double BmTotal=0;
	double DmTotal=0;
	double BsTotal=0;
	double DsTotal=0;

	for(int j=0; j<Nv; j++)
	{
		BmTotal+=Bm[j][i]*RAsim[start][j];
		DmTotal+=Dm[j][i]*RAsim[start][j];
	}

	for(int l=0; l<Lv; l++)
	{
		BsTotal+=Bs[l][i]*Stim[start][l];
		DsTotal+=Ds[l][i]*Stim[start][l];
	}

	RAsim[start+1][i]=RAsim[start][i]*(1+b[i]-d[i]+BmTotal-DmTotal+BsTotal-DsTotal);
}


void gLV::simNextStepGroup(vector<int> &v, int start)
{
	for(auto i : v)
	{
		simNextStep(i, start);
	}
}

void gLV::simNextStepPara(int p, int start, int normalize)
{
	if(Nv<=p) p=Nv;

	vector<vector<int>> buck(p);

	//assign microbes to p buckets, run each bucket on a thread
	int bi=0;
	for(int i=0; i<Nv; i++)
	{
		if(bi>=p) bi=0;
		buck[bi].push_back(i);
		++bi;
	}
	
	thread th[p];
	for(int ti=0; ti<p; ti++)
	{
		th[ti]=thread(&gLV::simNextStepGroup, this, std::ref(buck[ti]), start);
	}
	for(int ti=0; ti<p; ti++)
	{
		th[ti].join();
	}


	//normalize or not
	if(normalize)
	{
		vnorm(std::ref(RAsim[start+1]));
	}
}

void gLV::sim(int p, int start, int end, int normalize)
{
	for(int i=start; i<end; ++i)
	{
		simNextStepPara(p, i, normalize);
	}
}

void gLV::printPara()
{
	cout<<"RAsimT0:";
	for(auto ra : RAsim[0])
	{
		cout<<"\t"<<ra;
	}
	cout<<endl;

	cout<<"Stim:";
	for(auto st : Stim)
	{
		for(auto sti : st)
			cout<<"\t"<<sti;
	}
	cout<<endl;

	cout<<"b:";
	for(auto bi : b)
	{
		cout<<"\t"<<bi;
	}
	cout<<endl;

	cout<<"d:";
	for(auto di : d)
	{
		cout<<"\t"<<di;
	}
	cout<<endl;

	cout<<"Bm:";
	for(auto bmi : Bm)
	{
		for(auto bmii : bmi)
			cout<<"\t"<<bmii;
	}
	cout<<endl;

	cout<<"Dm:";
	for(auto dmi : Dm)
	{
		for(auto dmii : dmi)
			cout<<"\t"<<dmii;
	}
	cout<<endl;

	cout<<"Bs:";
	for(auto bsi : Bs)
	{
		for(auto bsii : bsi)
			cout<<"\t"<<bsii;
	}
	cout<<endl;

	cout<<"Ds:";
	for(auto dsi : Ds)
	{
		for(auto dsii : dsi)
			cout<<"\t"<<dsii;
	}
	cout<<endl;
}


void gLV::printRAsim(int t)
{
	int i=0;
	for(auto rai : RAsim[t])
	{
		cout<<"species "<<i<<":\t"<<rai<<endl;
		++i;
	}
	cout<<"Done."<<endl;
}


//==================Extended Kalman filter=====================

struct eKF
{
	gLV model;
	
	//The state vector z at time t
	//first N elements - predicted relative abundances of the N microbes at time t
	//next N elements - birth rates
	//next N elements - death rates
	//next N*N - Bm
	//next N*N - Dm
	//next L*N - Bs
	//next L*N - Ds

	vector<vec> z; //T * (3N+2N^2+2LN)
	void zToModel(int t);
	void modelToz(int t);
	void predictz(int start, int end); //predict z[end] based on z[start], usually start=t-1, end=t
	void updatez(int t); //update based on measurements

	//error covariance matrix P at time t

	vector<mat> P; 
	void predictP(int start, int end); //predict P[end] based on P[start], usually start=t-1, end=t
	void updateP(int t); //update P based on measurements

	//other matrixes that needs calculation
	vector<mat> H; //at time k, H is (N+1) * (3N+2N^2+2LN)
	vector<mat> K; //Kalman gain matrix
	vector<mat> Omega; //at time k, Omega is (3N+2N^2+2LN) * (3N+2N^2+2LN)

	void computeH(int t);
	void computeK(int t);
	void computeOmega(int t);

	//operations
	void initialize(gLV &glv); //initialize z0 and P0 by the dimension of N,L,T in model
	
void eKF::modelToz(int t)
{
	vector<double> tmp;

	for(auto ra : model.RAsim[t])
		tmp.push_back(ra);

	for(auto bi : model.b)
		tmp.push_back(bi);

	for(auto di : model.d)
		tmp.push_back(di);

	for(auto bm : model.Bm)
	{
		for(auto bmi : bm)
			tmp.push_back(bmi);
	}

	for(auto dm : model.Dm)
	{
		for(auto dmi : dm)
			tmp.push_back(dmi);
	}	

	for(auto bs : model.Bs)
	{
		for(auto bsi : bs)
			tmp.push_back(bsi);
	}		

	for(auto ds : model.Ds)
	{
		for(auto dsi : ds)
			tmp.push_back(dsi);
	}	

	z[t]=vec(tmp);
}

void eKF::zToModel(int t)
{
	//first N elements - predicted relative abundances of the N microbes at time t
	//next N elements - birth rates
	//next N elements - death rates
	//next N*N - Bm
	//next N*N - Dm
	//next L*N - Bs
	//next L*N - Ds

	int zi=0;
	for(int i=0; i<model.N; ++i,++zi)
		model.RAsim[t][i]=z[t](zi);

	for(int i=0; i<model.N; ++i,++zi)
		model.b[i]=z[t](zi);

	for(int i=0; i<model.N; ++i,++zi)
		model.d[i]=z[t](zi);

	for(int i=0; i<model.N; ++i,++zi)
	{
		for(int j=0; j<model.N; ++j, ++zi)
			model.Bm[i][j]=z[t](zi);
	}

	for(int i=0; i<model.N; ++i,++zi)
	{
		for(int j=0; j<model.N; ++j, ++zi)
			model.Dm[i][j]=z[t](zi);
	}		

	for(int l=0; l<model.L; ++l,++zi)
	{
		for(int j=0; j<model.N; ++j, ++zi)
			model.Bs[l][j]=z[t](zi);
	}

	for(int l=0; l<model.L; ++l,++zi)
	{
		for(int j=0; j<model.N; ++j, ++zi)
			model.Ds[l][j]=z[t](zi);
	}

}


void eKF::initialize()
{
	model=glv;
	//TODO some initialization of the dimension of z and P

	//compute P0 and z0
	modelToz(0);
	P[0]=cov(z[0]);
}


int main()
{
	gLV glv;
	ifstream input("config");
	glv.inputInitialCondition(input);
	glv.printPara();
	glv.simNextStepPara(4, 0, 1);
	glv.printRAsim(1);
	glv.sim(4, 0, 3, 1);
	glv.printRAsim(3);
}
