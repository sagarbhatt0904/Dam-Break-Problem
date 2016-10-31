// McCormack's Method for solution to Dam Break Problem

#include<iostream>
#include<fstream>
#include<math.h>

using namespace std;

int main()
{	
	ofstream fout("McCormack_plot.csv");
	int n=100;
	double H[n], H1[n], H_old[n], v_old[n], v[n], v1[n], Q[n], Q1[n], Q_old[n], E_old[n], E1[n]; 
	for(int i=0;i<n;i++)
	{
		if(i<=50)
		{
			H_old[i]=10;
			H[i]=10;
			H1[i]=0;
		}
		else
		{
			H_old[i]=1;
			H[i]=1;
			H1[i]=0;
		}
		v_old[i]=0;
		v[i]=0;
		v1[i]=0;
		Q[i]=0;
		Q1[i]=0;
		Q_old[i]=0;
		E_old[i]=0;
		E1[i]=0;
	}
	int L=100;
	double dx=1;
	double dt=0.01;
	double k=(dt/dx);
	double g=9.81;
	for(double t=0;t<4;t+=dt)
	{
		for(int i=1;i<n;i++)
		{	Q[i]=(H[i]*v[i]);
			Q_old[i]=(H_old[i]*v_old[i]);
			E_old[i]=(H_old[i]*(pow(v_old[i],2)))+(0.5*g*(pow(H_old[i],2)));
			E1[i]=(H1[i]*(pow(v1[i],2)))+(0.5*g*(pow(H1[i],2)));
			H1[i]=H_old[i]-k*(Q_old[i+1]-Q_old[i]);
			Q1[i]=Q_old[i]-k*(E_old[i+1]-E_old[i]);
			H[i]=0.5*(H_old[i]+H1[i]-k*(Q1[i]-Q1[i-1]));
			Q[i]=0.5*(Q_old[i]+Q1[i]-k*(E1[i]-E1[i-1]));
		}
		H[0]=H[1];
		H[n-1]=H[n-2];
		Q[0]=Q[1];
		Q[n-1]=Q[n-2];
		H1[0]=H1[1];
		H1[n-1]=H1[n-2];
		Q1[0]=Q1[1];
		Q1[n-1]=Q1[n-2];
		for (int i=0;i<n; i++)
		{
			v[i]=Q[i]/H[i];
			v_old[i]=Q_old[i]/H_old[i];
			H_old[i]=H[i];
			Q_old[i]=Q[i];
		}
	}
	for (int i=0; i<n; i++)
	{
		fout<<i+1<<","<<H[i]<<endl;
	}

	fout.close(); 
	return 0;
}
