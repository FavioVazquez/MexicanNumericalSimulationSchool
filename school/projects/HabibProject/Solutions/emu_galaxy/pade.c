#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "extra_logk.h"
#include "extra_pk.h"

double func(double x, double a0, double b1, double b2); 

int pade(int nk, double *k, double *Pk, int npts, double *k_out, double *Pk_out)
{
  int i; 
  double x0, y0, x1, y1, x2, y2; 

  int extra_pts = 51; 

  x0 = k[nk-100]; y0 = Pk[nk-100]; 
  x1 = k[nk-50]; y1 = Pk[nk-50]; 
  x2 = k[nk-1]; y2 = Pk[nk-1]; 

  double a0 = 0; 
  double b1 = 0; 
  double b2 = 0;
 
  b2 = (y1 - y0) + (y1 - y2)*(y0*x0-y1*x1)/(y1*x1-y2*x2); 
  b2 /= ((y0*x0*x0-y1*x1*x1)+(y2*x2*x2-y1*x1*x1)*(y0*x0-y1*x1)/(y1*x1-y2*x2)); 


  b1 = (y2-y1)+b2*(y2*x2*x2-y1*x1*x1); 
  b1 /= (y1*x1-y2*x2); 


  a0 = y0*(1.+b1*x0+b2*x0*x0); 


  float kmin = -7; 
  float kmax = log10(k[0]); 
  int nk_out = 0; 
  float binwidth = (kmax-kmin)/(float)(extra_pts-1); 

  float n = 0.963;
  float A = Pk[0]/pow(k[0],n);

  // extend to larger scales using power law
  for (i = 0; i < extra_pts; i++)
    {
      double new_k = pow(10.,kmin+binwidth*i);
      if (i != (extra_pts-1))
  	{
  	  k_out[nk_out] = new_k;
  	  Pk_out[nk_out] = A*pow(new_k, n);
  	  nk_out++;
  	}
      //      fprintf(stderr,"%lf %lf %lf %d\n", binwidth, new_k, Pk_out[i], nk_out);
    }

  //  fprintf(stderr, "%d\n", nk_out); 


  //original power spectrum
  for (i = 0; i < nk; i++)
    {
      k_out[nk_out] = k[i]; 
      Pk_out[nk_out] = Pk[i]; 
      nk_out++;
    }
  //  fprintf(stderr, "%d %d\n", nk_out, npts); 

  kmin = log10(k[nk-1]); 
  kmax = 4.;
  //  npts = 51; 
  binwidth = (kmax-kmin)/(float)(extra_pts-1); 


  /* float K = log10(-1*derivs[nk-2])/(A*log10(k[nk-2])); // some factors bundled into a single constant */
  /* float determinant = 1.- 4*K;  */
  /* float r1 = (1 + sqrt(determinant))*0.5/K; */
  /* float r2 = (1 - sqrt(determinant))*0.5/K; */
  /* fprintf(stderr, "%f %f %f %f\n", derivs[nk-2], sqrt(determinant), r1, r2);  */


  // extend to smaller scales using Pade approximation
  for (i = 0; i < 51; i++)
    {
      double new_k = pow(10.,kmin+binwidth*i); 
      double pade = func(new_k, a0, b1, b2); 
      if (i !=0)
      	{
       	  k_out[nk_out] = new_k;
      	  Pk_out[nk_out] = pade;
      	  nk_out++;
      	}
      //      fprintf(stderr,"%lf %lf %lf %d\n", kmin+binwidth*i, new_k, Pk_out[661+i], nk_out);
    }

  return(0);
}

double func(double x, double a0, double b1, double b2)
{

  double f = a0/(1+b1*x+b2*x*x); 

  return(f); 
}
