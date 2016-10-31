
// Runge-Kutta Method for solution to Dam Break Problem

#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>

using namespace std;

int sign(double sgn)
{
	if (sgn==0)
	{return 0;}
	else
	{	if (sgn<0)
			{return -1;}
		 else
		{
			if (sgn>0)
			{return 1;}
		}
	}
}
int main()
{
	
	int n=100;
	double H_old[n],H[n],H1[n],H2[n],H3[n],H4[n],v_old[n],v[n],Q[n],Q_old[n],E_old[n],Q1[n],Q2[n],Q3[n],Q4[n],E1[n],E2[n],E3[n],E4[n],a_old_H[n],a_H[n],theta_old_H[n],theta_H[n],theta_old_Q[n],theta_Q[n],G_old_H[n],G_H[n],G_old_Q[n],G_Q[n],psi_old_H[n],psi_H[n],psi_old_Q[n],psi_Q[n],a_old_Q[n],a_Q[n];
	int sigma_old_H[n],sigma_H[n],sigma_old_Q[n],sigma_Q[n];
	for(int i=0; i<n;i++)
	{    if( i<=50)
		{	
			H_old[i]=10;
			H[i]=10;
			H1[i]=10;
			H2[i]=10;
			H3[i]=10;
			H4[i]=10;
		}
	    else
	    {
			H_old[i]=1;
			H[i]=1;
			H1[i]=1;
			H2[i]=1;
			H3[i]=1;
			H4[i]=1;
	    }
	    v_old[i]=0;
	    v[i]=0;
	    Q[i]=0;
	    Q_old[i]=0;
	    E_old[i]=0;
	    Q1[i]=0;
	    Q2[i]=0;
	    Q3[i]=0;
	    Q4[i]=0;
	    E1[i]=0;
	    E2[i]=0;
	    E3[i]=0;
	    E4[i]=0;
	    a_old_H[i]=0;
	    a_H[i]=0;
	    sigma_old_H[i]=0;
	    sigma_H[i]=0;
	    sigma_old_Q[i]=0;
	    sigma_Q[i]=0;
	    theta_old_H[i]=0;
	    theta_H[i]=0;
	    theta_old_Q[i]=0;
	    theta_Q[i]=0;
	    G_old_H[i]=0;
	    G_H[i]=0;
	    G_old_Q[i]=0;
	    G_Q[i]=0;
	    psi_old_H[i]=0;
	    psi_H[i]=0;
	    psi_old_Q[i]=0;
	    psi_Q[i]=0;
	}
	int L=100;
	double dx=1;
	double dt=0.001;
	double k=(dt/dx);
	double g=9.81;
	for(double t=0; t<=4; t+=dt)
	{    for (int i=1; i<n-1; i++)
	    {
		if( (H_old[i]-H_old[i-1])<=0.0005)
		{    
			a_old_H[i]=H_old[i];
		}
		else
		{    
			a_old_H[i]=((Q_old[i]-Q_old[i-1])/(H_old[i]-H_old[i-1]));
		}
		
		if( (H_old[i+1]-H_old[i])<=0.0005)
		{    
			a_H[i]=H_old[i];
		}
		else
		{    
			a_H[i]=((Q_old[i+1]-Q_old[i])/(H_old[i+1]-H_old[i]));
		}
		
		if((Q_old[i]-Q_old[i-1])<=0.0005)
		{    
			a_old_Q[i]=Q_old[i];
		}
		else
		{    
			a_old_Q[i]=((E_old[i]-E_old[i-1])/(Q_old[i]-Q_old[i-1]));
		}
		
		if(( Q_old[i+1]-Q_old[i])<=0.0005)
		{
		    a_Q[i]=Q_old[i];
		}
		else
		{
		   a_Q[i]=((E_old[i+1]-E_old[i])/(Q_old[i+1]-Q_old[i]));
		 }
		
		sigma_old_H[i]=sign(a_old_H[i]);
		sigma_H[i]=sign(a_H[i]);
		sigma_old_Q[i]=sign(a_old_Q[i]);
		sigma_Q[i]=sign(a_Q[i]);
		if( (H_old[i]-H_old[i-1])==0)
		{
		    theta_old_H[i]=0;
		 }
		else
		{
		    theta_old_H[i]=(H_old[(i+sigma_old_H[i])]-H_old[(i-1+sigma_old_H[i])])/(H_old[i]-H_old[i-1]);
		 }
		
		if(( H_old[i+1]-H_old[i])==0)
		 {
		    theta_H[i]=0;
		  }
		else
		 {
		    theta_H[i]=(H_old[(i+1+sigma_H[i])]-H_old[(i+sigma_H[i])])/(H_old[i+1]-H_old[i]);
		  }
		
		if(( Q_old[i]-Q_old[i-1])==0)
		 {
		    theta_old_Q[i]=0;
		  }
		else
		{
		    theta_old_Q[i]=(Q_old[(i+sigma_old_Q[i])]-Q_old[(i-1+sigma_old_Q[i])])/(Q_old[i]-Q_old[i-1]);
		 }
		
		if(( Q_old[i+1]-Q_old[i])==0)
		 {
		    theta_Q[i]=0;
		  }
		else
		 {
		    theta_Q[i]=(Q_old[(i+1+sigma_Q[i])]-Q_old[(i+sigma_Q[i])])/(Q_old[i+1]-Q_old[i]);
		   }
		
		G_old_H[i]=(theta_old_H[i]+abs(theta_old_H[i]))/(1+theta_old_H[i]);
		G_H[i]=(theta_H[i]+abs(theta_H[i]))/(1+theta_H[i]);
		G_old_Q[i]=(theta_old_Q[i]+abs(theta_old_Q[i]))/(1+theta_old_Q[i]);
		G_Q[i]=(theta_Q[i]+abs(theta_Q[i]))/(1+theta_Q[i]);
		psi_old_H[i]=(0.5*G_old_H[i]*(abs(a_old_H[i])+k*(pow(a_old_H[i],2)))-abs(a_old_H[i]))*(H_old[i]-H_old[i-1]);
		psi_H[i]=(0.5*G_H[i]*(abs(a_H[i])+k*(pow(a_H[i],2)))-abs(a_H[i]))*(H_old[i+1]-H_old[i]);
		psi_old_Q[i]=(0.5*G_old_Q[i]*(abs(a_old_Q[i])+k*(pow(a_old_Q[i],2)))-abs(a_old_Q[i]))*(Q_old[i]-Q_old[i-1]);
		psi_Q[i]=(0.5*G_Q[i]*(abs(a_Q[i])+k*(pow(a_Q[i],2)))-abs(a_Q[i]))*(Q_old[i+1]-Q_old[i]);
	    }
	    for(int i=0;i<n;i++)
	    {	Q[i]=(H[i]*v[i]);
		Q_old[i]=(H_old[i]*v_old[i]);
		E_old[i]=(H_old[i]*(pow(v_old[i],2)))+(0.5*g*(pow(H_old[i],2)));
		Q1[i]=(H1[i]*v_old[i]);
		Q2[i]=(H2[i]*v_old[i]);
		Q3[i]=(H3[i]*v_old[i]);
		Q4[i]=(H4[i]*v_old[i]);
		E1[i]=(H1[i]*(pow(v_old[i],2)))+(0.5*g*(pow(H1[i],2)));
		E2[i]=(H2[i]*(pow(v_old[i],2)))+(0.5*g*(pow(H2[i],2)));
		E3[i]=(H3[i]*(pow(v_old[i],2)))+(0.5*g*(pow(H3[i],2)));
		E4[i]=(H4[i]*(pow(v_old[i],2)))+(0.5*g*(pow(H4[i],2)));
	    }
	    for(int i=1;i<n-2;i++)
	    {	H1[i]=H_old[i];
		H2[i]=H_old[i]-0.25*0.5*k*(Q1[i+1]-Q1[i-1]);
		H3[i]=H_old[i]-0.33*0.5*k*(Q2[i+1]-Q2[i-1]);
		H4[i]=H_old[i]-0.5*0.5*k*(Q3[i+1]-Q3[i-1]);
		H[i]=H_old[i]-0.5*k*(Q4[i+1]-Q4[i-1]);
		Q1[i]=Q_old[i];
		Q2[i]=Q_old[i]-0.25*0.5*k*(E1[i+1]-E1[i-1]);
		Q3[i]=Q_old[i]-0.33*0.5*k*(E2[i+1]-E2[i-1]);
		Q4[i]=Q_old[i]-0.5*0.5*k*(E3[i+1]-E3[i-1]);
		Q[i]=Q_old[i]-0.5*k*(E4[i+1]-E4[i-1]);
	    }
	    for(int i=0;i<n;i++)
	    {	H[i]=H[i]-0.5*k*(psi_H[i]-psi_old_H[i]);
		Q[i]=Q[i]-0.5*k*(psi_Q[i]-psi_old_Q[i]);
	    }
	    H[0]=H[1];
	    H[n-1]=H[n-2];
	    Q[0]=Q[1];
	    Q[n-1]=Q[n-2];
	    for (int i=0;i<n;i++)
	    {
		    v[i]=Q[i]/H[i];
		    H_old[i]=H[i];
		    Q_old[i]=Q[i];
		    v_old[i]=v[i];
	    }
	
	}
	ofstream fout("rk_plot.csv");
	for (int i=0; i<n; i++)
	{
		fout<<i+1<<","<<H[i]<<endl;
	}

	fout.close(); 
	return 0;
}

