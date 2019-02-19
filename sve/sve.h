#pragma once

#include <string>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdio>     
#include <algorithm>  
#include <vector>     

//typedef std::complex<double> complex;
//const complex I(0,1);

const double 
	pi = 3.141592653589793238462642,
	eV = 27.2113834,   
	a0 = 0.529177208,
	c  = 137.03599971;
const double  nm    = 10.0/a0;
const double  kb_au = 3.1668151e-06;      // kB in 1/(au K)
const double  au_s  = 2.4188843e-17;    	// 1au of time in sec
const double  vf    = c/300.0;            // Fermi velocity ingraphene

// ---- number of rows in a file ----//
int norow(char *nin){
	int number_of_lines = 0;
	    FILE *fin = fopen(nin, "r");
	    int ch;
	    while (EOF != (ch=getc(fin)))
	        if ('\n' == ch)
	            ++number_of_lines;
		 fclose(fin);
	    return number_of_lines;
}

double sqr(double x){return x*x;}
//complex sqr(complex x){return x*x;}

double fact(int n){
  if(n<0) {printf("fact: negative number"); return 1;}
  else if(n==0) return 1;
  else{
  	 double   p=1;
  	 while(n) p*=n--;  
  	 return   p;
  }
}

