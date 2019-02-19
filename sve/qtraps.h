#pragma once

// Integration routine by TRAPEZOIDAL RULE
// a: lower integral limit
// b: upper integral limit 
// f(x): integrand function
// n: number of evaluations

#include "sve.h"

double trapzd(double a, double b, double f(double x), int n)
{
  double x,tnm,sum,del;
  static double s;
  int it,j;
  if (n == 1) {
    return (s=0.5*(b-a)*(f(a)+f(b)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm; // This is the spacing of the points to be added.
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += f(x);
    s=0.5*(s+(b-a)*sum/tnm); // This replaces s by its refined value.
    return s;
  }
}

/*
Returns the integral of the function func from a to b. The parameters
EPS can be set to the desired fractional accuracy and JMAX so that 2
to the power JMAX-1 is the maximum allowed number of steps.
Integration is performed by the trapezoidal rule.
*/

#define EPS 1.0e-5
#define JMAX 25
double qtrap(double a, double b,double f(double x))
{
  double trapzd(double a, double b, double f(double x), int n);
  int j;
  double s,olds;
  olds = -1.0e30;         // Any number that is unlikely to be the average
  for (j=1;j<=JMAX;j++) { //  of the function at its endpoints will do here.
    s=trapzd(a,b,f,j);
    if (j > 5)            // Avoid spurious early convergence.
    if (fabs(s-olds) < EPS*fabs(olds) || (s == 0.0 && olds == 0.0)) return s;
    olds=s;
    //printf("%d\n",j);
  }
  fprintf(stderr,"Too many steps in routine qtrap");
  printf("\n");  
  exit(1);
  return 0.0; // Never get here.
}