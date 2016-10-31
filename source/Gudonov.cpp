// Gudonov's Method for solution to Dam Break Problem

#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>

using namespace std;

int main()
{
	ofstream fout("gudonov_plot.csv");
	int n=100;
	double H[n], H_old[n],v_old[n],v[n],Q[n],Q_old[n],E_old[n],a[n],F1[n],F2[n];
	for( int i=0; i<n; i++)
	{	if (i<=50)
		{
			H_old[i]=10;
			H[i]=10;
		}
		else
		{
			H_old[i]=1;
			H[i]=1;
		}
		v_old[i]=0;
		v[i]=0;
		Q[i]=0;
		Q_old[i]=0;
		E_old[i]=0;
		F1[i]=0;
		F2[i]=0;
	}
	
	int L=100;
	double dx=1;
	double dt=0.01;
	double k=(dt/dx);
	double g=9.81;
	//int count=0;
	for(double t=0; t<=4; t+=dt)
	{	for(int i=1; i<n-1; i++)
		{	
			// Q[i]=(H[i]*v[i]);
			// Q_old[i]=(H_old[i]*v_old[i]);
			E_old[i]=(H_old[i]*(pow(v_old[i],2)))+(0.5*g*(pow(H_old[i],2)));
			//a[i]=fmax(abs((v_old[i]+sqrt(g*H_old[i]))),abs((v_old[i-1]+sqrt(g*H_old[i-1]))));
			
			F1[i]=0.5*(Q_old[i]+Q_old[i-1])-0.5*(H_old[i+1]-H_old[i]);
			F2[i]=0.5*(E_old[i]+E_old[i-1])-0.5*(Q_old[i+1]-Q_old[i]);
		}
		for(int i=1;i<n-1;i++)
		{	
			H[i]=H_old[i]-k*(F1[i]-F1[i-1]);
			Q[i]=Q_old[i]-k*(F2[i]-F2[i-1]);
		}
		
		H[0]=H[1];
		H[n-1]=H[n-2];
		Q[0]=Q[1];
		Q[n-1]=Q[n-2];
		for (int i=0; i<n; i++)
		{
			H_old[i]=H[i];
			Q_old[i]=Q[i];
			v_old[i]=v[i];			
			v[i]=Q[i]/H[i];
		}
		//cout<<count<<endl;
		//count+=1;
	}
	for (int i=0; i<n; i++)
	{
		fout<<i+1<<","<<H[i]<<endl;
	}

	fout.close(); 
	return 0;
}
