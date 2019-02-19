/*
   *** HEADER FOR J's.cpp ***
	Spectrally resolved tunnelling probability for a 
   perfectly stacked gr-hBN-gr device.
   Full-RPA for the conductivity.
   T = 0K
   Different doping combinations for the graphene sheets.
   Calculates Jeh(w) according to notes 11/feb/2019.
   I1(kp,w) is tabulated in the folder ../IntPhiWo 
   Tabulated frequencies are in w2.dat
      0.01 < w(eV) < 2.99667           nw=900
      0.01 < kp(1/nm) < 11.9891        nk=1100 
   Both of these are equally spaced
   Combinations of Ef are: 
      Ef2=0.3eV    Ef1=0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0eV
      Ef2=0.5eV    Ef1=1.0eV
*/ 
# pragma once

#include "sve/sve.h"
#include "sve/qtraps.h"
//#include "sve/interp_linear.h"

// w2.dat and IntPhiWo
int    nk = 1100;
int    nw = 900;

// w3.dat and IntW
//int nk = 510;
//int nw = 500;

double Ef1,Ef2,Vb,g0;
double w,kp;
double En1,En2;
double Qi01,Qi02,Qi03,Qi04;
double Qf01,Qf02,Qf03,Qf04;
double A1,A2,A3,A4;
double eta0;

int j,nev=60;
double sum;
double var1 = 0.001;
double var2 = 0.999*pi;
double *phi=new double[nev];

// --- Fermi-Dirac distribution at T=0K (Step function)
double FD(double En){
	if(En>0)       return 1.0;
	else if(En==0) return 0.5;
	else           return 0.0;
}

// --- Qi0: poles of the deltas
double Qi001(double kp,double vph,double g0){ 
	double num = kp*kp - g0*g0;
	double den = 2.0*(kp*cos(vph)-g0);

   double res;
   if(den==0) return 0.0;
   else res = num/den;
   
   if(res<0)  return 0.0;
	else       return res;
}

double Qi002(double kp,double vph,double g0){ 
	double num = kp*kp - g0*g0;
	double den = 2.0*(kp*cos(vph)+g0);
	//if(den==0) den = 1e-30;

   double res;
   if(den==0) return 0.0;
   else res = num/den;
   
   if(res<0) return 0.0;
	else      return res;
}

double Qi003(double kp,double vph,double g0){ 
	return Qi001(kp,vph,g0);
}


double Qi004(double kp,double vph,double g0){ 	
	return Qi002(kp,vph,g0);
}



// --- A = [1+Qbi*Qbf/Qi/Qf]/|F'(Qi0)| 
// ---     delta(F(Qi))=delta(Qi-Qi0)/|F'(Qi0)|
double A01(double kp,double vph,double g0){ 
	double Qi01 = Qi001(kp,vph,g0);
   double num = 2.0*Qi01 - g0 - kp*cos(vph);
   double den = fabs(-g0 + kp*cos(vph));
  	if(den==0) return 0.0;
   else return num/den;
}

double A02(double kp,double vph,double g0){ 
	double res = -g0 - kp*cos(vph);
   
   if(res<0)       return -1.0; 
   else if(res>0)  return  1.0; 
   else            return  0.0;
}

double A03(double kp,double vph,double g0){ 
	double res = g0 - kp*cos(vph);
   
   if(res<0)       return -1.0; 
   else if(res>0)  return  1.0; 
   else            return  0.0;
}

double A04(double kp,double vph,double g0){ 
	double Qi02 = Qi002(kp,vph,g0);
   double num = 2.0*Qi02 + g0 - kp*cos(vph);
   double den = fabs(g0 + kp*cos(vph));
   if(den==0) return 0.0;
   else return num/den;
}



// --- Qf0: Qf's evaluated at Qi0
double Qf001(double kp,double vph,double g0){
   double Qi01 = Qi001(kp,vph,g0);
   
   double res = Qi01 - g0;
   if(res<0) return 0.0;
	else      return res;
} 

double Qf002(double kp,double vph,double g0){
   double Qi02 = Qi002(kp,vph,g0);
   
   double res = - Qi02 - g0;
   if(res<0) return 0.0;
	else      return res;
} 

double Qf003(double kp,double vph,double g0){
   double Qi03m = Qi003(kp,vph,g0);
   
   double res = g0 - Qi03m;
   if(res<0) return 0.0;
	else      return res;
} 

double Qf004(double kp,double vph,double g0){
   double Qf04 = - Qf002(kp,vph,g0);
   
   if(Qf04<0) return 0.0;
	else       return Qf04;
} 