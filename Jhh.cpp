/*
   Spectrally resolved tunnelling probability for a 
   perfectly stacked gr-hBN-gr device.
   Full-RPA for the conductivity.
   T = 0K
   Hole (h) doped graphene to hole (h) doped graphene.
   Calculates Jhh(w) according to notes 11/feb/2019.
   I1(kp,w) is tabulated
   Combinations of Ef are: 
      Ef2=0.3eV    Ef1=0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0eV
      Ef2=0.5eV    Ef1=1.0eV
*/ 
#include "Jhh.h"

int main(int argc, char **argv)
{		
	if(argc<4){printf("./jhh.x  Ef1(eV)  Ef2(eV)  Vb(eV)  \n\n"); exit(1);}
	
	int    i,j,l,m0;
	 Ef1 = atof(argv[1])/eV;      // Fermi energy of layer 1
	 Ef2 = atof(argv[2])/eV;      // Fermi energy of layer 2
	 Vb  = atof(argv[3])/eV;	    // bias voltage
 	
  double Jw_cv,Jw_vc,Jw_cc,Jw_vv,Jw;
  char   ninw[150], nin[200];         // input files
	char   nout[90];                    // output files


	//--- Output for the spectral tunnelling probability
  sprintf(nout,"Jhh_Ef1-%g_Ef2-%g_d1_Vb%g.dat",Ef1*eV,Ef2*eV,Vb*eV);  
	FILE *fout=fopen(nout,"w");
	   	
  // --- Tabulated frequencies     
	//sprintf(ninw,"./IntPhiWo/w2.dat"); // sonic
  sprintf(ninw,"../w2.dat");
  //sprintf(ninw,"../2_EFdep_gr-hBN-gr_no-rotat/IntPhiW/w3.dat");
	int     nw = norow(ninw);
	FILE *finw = fopen(ninw,"r");
	double *ww = new double[nw];
	for(i=0;i<nw;++i) fscanf(finw,"%lf",ww+i);
	fclose(finw);   

  double  *kpt = new double[nk]; 
	double  *I1t = new double[nk]; 
  double Ivc_tot,Ivv_tot;
  
	for(l=0;l<nw;l++){
//    l = 50;
		w  = ww[l]/eV;
    g0 = (-Ef1+Ef2-Vb+w)/vf; 
    eta0 = (-Ef1+Ef2-Vb)/vf;

    //sprintf(nin,"./IntPhiWo/IntW_Ef1-%g_Ef2-%g_d1_w%g.dat",Ef1*eV,Ef2*eV,ww[l]); // sonic
    sprintf(nin,"../IntPhiWo/IntW_Ef1-%g_Ef2-%g_d1_w%g.dat",Ef1*eV,Ef2*eV,ww[l]);
    //sprintf(nin,"../2_EFdep_gr-hBN-gr_no-rotat/IntPhiW/IntW_Ef1-%g_Ef2-%g_d1_w%g.dat",Ef1*eV,Ef2*eV,ww[l]);
			FILE   *fin = fopen(nin,"r");
			double  *kpf = new double[nk]; // kp from file
			double  *I1f = new double[nk]; // I1 from file
			for(i=0;i<nk;++i) fscanf(fin,"%lf %lf",kpf+i,I1f+i);
			fclose(fin);

			for(i=0;i<nk;i++){
        kpt[i]=kpf[i]/nm;		
        I1t[i]=-I1f[i];		
      }
      delete [] kpf; kpf = NULL;
      delete [] I1f; I1f = NULL;
    
    Ivc_tot=Ivv_tot=0.0;
    for(i=0;i<(nk-1);i++){
      Ivc_tot+=0.5*(Ivc_kp(kpt[i+1],I1t[i+1])+Ivc_kp(kpt[i],I1t[i]))*(kpt[i+1]-kpt[i]);		      
      Ivv_tot+=0.5*(Ivv_kp(kpt[i+1],I1t[i+1])+Ivv_kp(kpt[i],I1t[i]))*(kpt[i+1]-kpt[i]);		
    }

    Jw_vc = Ivc_tot/pi/pi/pi/vf;  
    Jw_vv = Ivv_tot/pi/pi/pi/vf;  
    Jw = Jw_vc+Jw_vv;
	  fprintf(fout,"%g %g ",w*eV,Jw_vc*nm*nm);
    fprintf(fout,"%g %g \n",Jw_vv*nm*nm,Jw*nm*nm);		
    fflush(fout);
    //printf("-->  w(eV)=%g   J=%g\n",ww[l],Jw*nm*nm);
	} 
	delete [] ww;  ww  = NULL;
  delete [] kpt; kpt = NULL;
  delete [] I1t; I1t = NULL;
  
  fclose(fout);

  return 0;
}

