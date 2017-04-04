// Analytical Solution to Dam Break Problem

#include<iostream>
#include<fstream>
#include<math.h>

using namespace std;

int main()
{
	int L = 100;
	double dx = 1;
	double g=9.81;
	double t = 4;
	int n=100;
	double H1=10;
	double H2 = 1;
	double Hg = 3.962;
	double x1, x2,x3, H[n], p;
	ofstream fout("analytical_plot.csv");
	x1=-sqrt(g*H1)*t;
	x2= -sqrt(g*Hg)*t;
	x3=g*t;
	for(int i=0; i<n; i++)
	{
		if ((i-L/2)<=x1)
		{
			H[i] = H1;
		}	
	
		if ((i-L/2)<=x2 && (i-L/2)>x1)
		{
			p=(2*sqrt(g*H1)/3)-(i-L/2)/(3*t);
			H[i] = (pow(p,2))/g;
		}
		if ((i-L/2)<=x3 & (i-L/2)>x2)
		{
			H[i] = Hg;
		}
		if ((i-L/2)> x3)
		{
			H[i] = H2;
		}
		
		fout<<i+1<<","<<H[i]<<endl;
	}
	fout.close(); 
	return 0;
}


