/*
   *** HEADER FOR Jhe.cpp ***
	Spectrally resolved tunnelling probability for a 
   perfectly stacked gr-hBN-gr device.
   Full-RPA for the conductivity.
   T = 0K
   Hole (h) doped graphene to electron (e) doped graphene.
   Calculates Jhe(w) according to notes 11/feb/2019.
   I1(kp,w) is tabulated in the folder ../2_EF.../IntPhiW 
   Tabulated frequencies are in w3.dat
      0.001 < w(eV) < 2.192            nw=500
      1e-5 < kp(1/nm) < 83             nk=510  (not eq. spaced)
   Combinations of Ef are: 
      Ef2=0.3eV    Ef1=0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0eV
      Ef2=0.5eV    Ef1=1.0eV
	Notation: 
		c: conduction
		v: valence
		Iif_x: integrand for integral over x that corresponds to the tunnel from
		         i (c or v) to f (c or v)
*/ 

#include "J.h"

//double err=1.0e-21;  int nmax=6000;   // integration parameters

// --- Integral 1: val-cond
double Ivc_vph(double vph){ // Integrand for  Int d_vph
   Qi02 = Qi002(kp,vph,g0);
   Qf02 = Qf002(kp,vph,g0);
	A2   = A02(kp,vph,g0); 
   
   En1  = vf*Qi02 - Ef1;
   En2  = vf*Qf02 - Ef2;

   // allowed energy region
   if     (w<vf*(-2.0*Qi02-kp-eta0))  return 0.0;
   else if(w>vf*(-2.0*Qi02+kp-eta0))  return 0.0;
   else  return Qi02*A2*FD(En1)*FD(En2);
}

double Ivc_kp(double kp_,double itpI1){ // Integrand for  Int d_kp
   kp = kp_;
   
   sum = 0.0;
   for(j=0;j<(nev-1);j++){
      sum += 0.5*(Ivc_vph(phi[j+1])+Ivc_vph(phi[j]))*(phi[j+1]-phi[j]);
   }

	return kp*itpI1*sum;
}
