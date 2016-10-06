#include <math.h>
#include "gsl/gsl_sf_bessel.h"
#include"gsl/gsl_spline.h"

int pk_to_xi(int nk, double *k, double *Pk, double *r_out, double *xi_out)
{
  int n; 
  int max_sum = 300; 

  double *phi     = malloc(max_sum*sizeof(double)); 
  double *phidash = malloc(max_sum*sizeof(double)); 
  double *weights = malloc(max_sum*sizeof(double)); 
  double *besselterm = malloc(max_sum*sizeof(double)); 

  double *Pkint = malloc(nk*sizeof(double)); 
  double sigma_smooth = 1.;
  int i;
  for (i =0; i < nk; i++)
    Pkint[i] = pow(k[i],1.5)*Pk[i]*exp(-k[i]*k[i]*sigma_smooth*sigma_smooth);

  gsl_interp_accel *acc 
    = gsl_interp_accel_alloc ();
  gsl_spline *spline 
    = gsl_spline_alloc (gsl_interp_cspline, nk);

  gsl_spline_init (spline, k, Pkint, nk);

  int m; 
  double r_min = 1; 
  double r_max = 200; 

  for (n = 1; n < max_sum; n++)
    {
      int nzero = n;
      double zero = (double)n; // zero of Bessel/PI
      double h = 1./150.;
      double hr = h*zero;

      phi[n] = hr *tanh(0.5*M_PI * sinh(hr));
      phidash[n] = tanh(0.5*M_PI*sinh(hr))+0.5*M_PI*hr*cosh(hr)/cosh(M_PI*0.5*sinh(hr))/cosh(M_PI*0.5*sinh(hr));
      weights[n] = gsl_sf_bessel_Ynu(0.5,M_PI*zero)/gsl_sf_bessel_Jnu(1.5,M_PI*zero);
      besselterm[n] = gsl_sf_bessel_Jnu(0.5, M_PI*phi[n]/h);

      //      weights[n] = -1.*sqrt(2./(M_PI*M_PI*zero))*cos(M_PI*zero)/(gsl_sf_bessel_Jnu(1.5,M_PI*zero));
    }

  for (m = 0; m < 200; m++)
    {
      double r = log10(r_min) + m*((log10(r_max) - log10(r_min))/200.); // r = 80 is the limit for k = 1; 
      r = pow(10.,r); 
      //      r *= da; 
      double sum = 0; 
      /* for (n=0; n < nk; n++) */
      /* 	Pk[n] /= r;  */

      for (n = 1, sum=0; n < max_sum; n++)
	{
	  double nu = 0.5; 
	  int nzero = n; 
	  double zero = (double)n; //gsl_sf_bessel_zero_Jnu (2,nzero); // zero of Bessel/PI
 
	  double h = 1./150; 
	  double hr = h*zero;
	  /* phi[n] = hr *tanh(0.5*M_PI * sinh(hr)); */
	  /* phidash[n] = tanh(0.5*M_PI*sinh(hr))+0.5*M_PI*hr*cosh(hr)/cosh(M_PI*0.5*sinh(hr))/cosh(M_PI*0.5*sinh(hr)); */
	  /* double weights = gsl_sf_bessel_Ynu(2,M_PI*zero)/gsl_sf_bessel_Jnu(3,M_PI*zero); */
	  /* double weights = -1.*sqrt(2./(M_PI*M_PI*zero))*cos(M_PI*zero)/(gsl_sf_bessel_Jnu(1.5,M_PI*zero)); */
	  double x = M_PI*phi[n]/h; 
	  //double pkterm   = pow(x,1.5)*gsl_spline_eval (spline, x/r, acc);
	  double pkterm   = gsl_spline_eval (spline, x/r, acc);
	  /* double besselterm = gsl_sf_bessel_Jnu(0.5,x); */
	  /* double hr_p = h*(zero+1.);  */
	  /* double hr_m = h*(zero-1.);  */
	  /* double hr_pp = h*(zero+2.);  */
	  /* double hr_mm = h*(zero-2.);  */

	  /* double phi_p = hr_p *tanh(0.5*M_PI * sinh(hr_p)); */
	  /* double phi_m = hr_m *tanh(0.5*M_PI * sinh(hr_m)); */
	  /* double phi_pp = hr_pp *tanh(0.5*M_PI * sinh(hr_pp)); */
	  /* double phi_mm = hr_mm *tanh(0.5*M_PI * sinh(hr_mm)); */
	  /* //	  double derivs = (phi[n]-phi[n-1])/(h);  */
	  /* double derivs = (-phi_pp+8*phi_p-8*phi_m+phi_mm)/(12.*h);  */
	  //	  sum += M_PI*weights*pkterm*besselterm*derivs;
	  sum += M_PI*weights[n]*pkterm*besselterm[n]*phidash[n];

	  /* fprintf(stderr, "%d %f %f %f %f %f %f \n", n, sum, pkterm, besselterm, weights,  phi[n], derivs); */
	  /* fprintf(stderr, "%d %f %f %f %f \n", n, x, x/r, phi[n],  phidash[n] ); */
	}
      //double constant = sqrt(0.5*M_PI)/(2.*M_PI*M_PI)/(r*r*r);
      double constant = sqrt(0.5*M_PI)/(2.*M_PI*M_PI)/pow(r,1.5);
      r_out[m] = r; 
      xi_out[m] = sum*constant;
    }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);


  free(phi); free(phidash); free(weights); free(Pkint); free(besselterm);

  return(0); 
}
