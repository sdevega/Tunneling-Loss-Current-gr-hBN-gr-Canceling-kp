/*
   *** HEADER FOR Jee.cpp ***
	Spectrally resolved tunnelling probability for a 
   perfectly stacked gr-hBN-gr device.
   ELIMINATING Qi!! Final integral over kp
   Full-RPA for the conductivity.
   T = 0K
   Electron (e) doped graphene to electron (e) doped graphene.
   Calculates Jeh(w) according to notes 15/feb/2019.
   I1(kp,w) is tabulated.   
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


// --- Integral 1: val-cond
double Ivc_vph(double vph){ // Integrand for  Int d_vph
   Qi02 = Qi002(kp,vph,g0);
   Qf02 = Qf002(kp,vph,g0);
	A2   = B02(kp,vph,g0); 
   En2  = Ef2 - vf*Qf02;

   if     (w<vf*(-2.0*Qi02-kp-eta0)) return 0.0;
   else if(w>vf*(-2.0*Qi02+kp-eta0)) return 0.0;
   else  return  Qi02*A2*FD(En2);
}

double Ivc_kp(double kp_,double itpI1){ // Integrand for  Int d_kp
   kp = kp_;
   return kp*itpI1*apt.integrate(Ivc_vph,var1,var2);
}



// --- Integral 2: cond-cond
double Icc_vph(double vph){ // Integrand for  Int d_vph
   Qi01 = Qi001(kp,vph,g0);
   Qf01 = Qf001(kp,vph,g0);
	A1   = A01(kp,vph,g0); 
	En1  = Ef1 - vf*Qi01;
   En2  = vf*Qf01 - Ef2;
		
   return Qi01*A1*FD(En1)*FD(En2);
}

double Icc_kp(double kp_,double itpI1){ // Integrand for  Int d_kp
   kp = kp_;
   if     (w<vf*(-kp-eta0)) return 0.0;
   else if(w>vf*( kp-eta0)) return 0.0;
   else return kp*itpI1*apt.integrate(Icc_vph,var1,var2);
}
