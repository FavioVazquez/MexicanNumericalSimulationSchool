#include <stdlib.h>
#include <stdio.h>
#include "emu.h"
#include "params.h"
#include "design.h"
#include "dist_mat.h"
#include <math.h>


int nk = 661;   // Number of k bins
int numPC = 10;    // Number of pricipal components
int nmodels = 100; // How many models to cover parameter space
int nmodels_0_4985 = 149; 
int nparams = 5; // Number of cosmological parameters

int make_sigma_w(int nmodels, float designparams[][nparams], double *newparams, gsl_matrix *V11, gsl_matrix *V21, int PCnow, int numz);
int invert_matrix(int size, gsl_matrix *A, gsl_matrix *A_inv);

int emu(double *newparams, double outputredshift, double *output_pk)
{
  int n; 
  int offset = 0; 
  char inputs[256]; 

    // xstar contains the 5 emulator parameters plus the red shift.
  double xstar[6], stuff[4]; 
  xstar[0] = newparams[0];
  xstar[1] = newparams[1];
  xstar[2] = newparams[2];
  xstar[3] = newparams[4];
  xstar[4] = newparams[3];
  xstar[5] = outputredshift;

  //identify output redshift amongst input snapshots
  double scalefactor[6] = {0.4985, 0.6086, 0.6974, 0.8051, 0.9083, 1.0};
  double redshifts[6];

  for (n = 0; n < 6; n++)
    redshifts[n] = 1./scalefactor[n]-1; 

  double output_scalefactor = 1./(outputredshift+1.); 

  int numz_want = 0; 
  for (n=0; n < 6; n++)
    {
      if (output_scalefactor > scalefactor[n])
	{
	  numz_want = n;
	}
    }
  //  fprintf(stderr, "%d %lf %lf %d\n", numz_want, scalefactor[numz_want], scalefactor[numz_want+1], nk);

  double wpred[numPC]; //predicted weights at new P(k) relation
  int numz = 0; // redshift 
  double *pkpred = calloc(nk*6,sizeof(double));
  
  for (numz = numz_want; numz < numz_want+2; numz++)
  /* for (numz = 0; numz < 6; numz++) */
    {
      if (numz == 0)
	{
	  gsl_matrix * V11 = gsl_matrix_alloc (nmodels_0_4985, nmodels_0_4985);
	  gsl_matrix * V11_inv = gsl_matrix_alloc (nmodels_0_4985, nmodels_0_4985);
	  gsl_matrix * V21 = gsl_matrix_alloc (1, nmodels_0_4985); // leave this as a matrix in case several interpolations are wanted at once; 

	  for (n = 0; n < numPC; n++)
	    {
	      wpred[n] = 0;
	      
	      make_sigma_w(nmodels_0_4985, designparams_0_4985, newparams, V11, V21,n, numz);

	      invert_matrix(nmodels_0_4985, V11, V11_inv);

	      int i, j;

	      double *dummyreslt  = calloc(nmodels_0_4985,sizeof(double));
	      for (i = 0; i < nmodels_0_4985; i++)
		{
		  for(j = 0; j < nmodels_0_4985; j++)
		    {
		      dummyreslt[i] += gsl_matrix_get(V11_inv, i, j)*what_0_4985[j+n*nmodels_0_4985];
		    }
		  //		  fprintf(stderr, "%f\n", dummyreslt[i]);
		}
	      for (i = 0; i < nmodels_0_4985; i++)
		wpred[n] += gsl_matrix_get(V21, 0, i)*dummyreslt[i];

	      free(dummyreslt);

	      //	      fprintf(stderr, "%f\n", wpred[n]);

	    }
	  gsl_matrix_free(V11);
	  gsl_matrix_free(V11_inv);
	  gsl_matrix_free(V21);

	}
      else
	{
	  gsl_matrix * V11 = gsl_matrix_alloc (nmodels, nmodels);
	  gsl_matrix * V11_inv = gsl_matrix_alloc (nmodels, nmodels);
	  gsl_matrix * V21 = gsl_matrix_alloc (1, nmodels); // leave this as a matrix in case several interpolations are wanted at once; 

	  for (n = 0; n < numPC; n++)
	    {
	      wpred[n] = 0;
	      
	      make_sigma_w(nmodels, designparams, newparams, V11, V21,n, numz);

	      invert_matrix(nmodels, V11, V11_inv); 

	      int i, j;
	      double *dummyreslt  = calloc(nmodels,sizeof(double));
	      for (i = 0; i < nmodels; i++)
		{
		  for(j = 0; j < nmodels; j++)
		    dummyreslt[i] += gsl_matrix_get(V11_inv, i, j)*what[j+n*nmodels][numz-1];
		}
	      for (i = 0; i < nmodels; i++)
		wpred[n] += gsl_matrix_get(V21, 0, i)*dummyreslt[i];

	      free(dummyreslt);

	      /* if(numz==5) */
	      /* 	fprintf(stderr, "%f\n", wpred[n]); */
	    }

	  gsl_matrix_free(V11);
	  gsl_matrix_free(V11_inv);
	  gsl_matrix_free(V21);
	}

      int m; 


      for (n = 0; n < nk; n++)
	{
	  if (logk[n] > -2)       /* // P(k) emulator  */
	    {

	      if (offset==0)
		offset=n;

	      pkpred[6*n+numz] = ymean[n-offset][numz];
	      for (m = 0; m < numPC; m++)
		{
		  pkpred[6*n+numz] += phi[n-offset][numz*numPC+m]*wpred[m]*ysimstd[n-offset][numz]; 
		}
	    }
	  else
	    {
      // large scale modes for xi calculation only
	      pkpred[6*n+numz] = ymean[0][numz];
	      for (m = 0; m < numPC; m++)
		{
		  pkpred[6*n+numz] += phi[0][numz*numPC+m]*wpred[m]*ysimstd[0][numz];
		}
	    }

	}


    }

  gsl_spline *spline_bias = gsl_spline_alloc(gsl_interp_linear, 2); 
  gsl_interp_accel *acc_bias  = gsl_interp_accel_alloc ();

  gsl_spline *spline_pk = gsl_spline_alloc(gsl_interp_linear, 2); 
  gsl_interp_accel *acc_pk  = gsl_interp_accel_alloc ();




  //  Convert back to P(k) from log10(bias/k)
  for (n = 0; n < nk; n++)
    {
      gsl_spline_init (spline_bias, &(scalefactor[numz_want]), &(pkpred[6*n+numz_want]), 2);
      gsl_spline_init (spline_pk, &(scalefactor[numz_want]), &(pk_m[n][numz_want]), 2);
      /* gsl_spline_init (spline_bias, scalefactor, &(pkpred[6*n]), 6); */
      /* gsl_spline_init (spline_pk, scalefactor, &(pk_m[n][0]), 6); */

      double bias = gsl_spline_eval(spline_bias, output_scalefactor, acc_bias);
      double k_unlogged;
      if (logk[n] > -2)
	k_unlogged = pow(10.,logk[n]);
      else
	k_unlogged = pow(10.,logk[offset]); // these modes are extrapolated with the linear bias (for xi only)
      bias = pow(10.,bias);
      output_pk[n] = bias*k_unlogged*gsl_spline_eval(spline_pk, output_scalefactor, acc_pk);
           /* fprintf(stderr, "%f %f %f %f %f %f\n", pkpred[6*n+0], pkpred[6*n+1], pkpred[6*n+2], pkpred[6*n+3], pkpred[6*n+4], pkpred[6*n+5]); */
      gsl_interp_accel_reset(acc_bias);
      gsl_interp_accel_reset(acc_pk);
    }

  gsl_spline_free(spline_bias);   gsl_spline_free(spline_pk);
  gsl_interp_accel_free(acc_bias);   gsl_interp_accel_free(acc_pk);
  free(pkpred);

  return(0);
}


int make_sigma_w(int nmodels, float designparams[][nparams], double *newparams, gsl_matrix *V11, gsl_matrix *V21, int PCnow, int numz)
{
  // Do each PC individually. 

  int i, j, k; 
  for (i = 0; i < nmodels; i++)
    {
    for (j = 0; j < nmodels; j++)
    /* for (j = 0; j < i; j++) */
      {
	double V11dummy = 1.; 
	for (k = 0; k < nparams; k++)
	  {
	    double distance = designparams[i][k]-designparams[j][k]; // always substract the same type of parameter from the same type i.e. Omegam_1 - Omegam_2 not Omegam_1-w 
	    distance = 4.*distance*distance; 
	    V11dummy *= pow(rho_w[numz*nparams+k][PCnow], distance); 
	  }
	V11dummy /= lambda_w[numz*numPC+PCnow]; 
	if (i==j)
	  V11dummy += 1./lambdaP[numz];
	  //	V11dummy += invphi[i+PCnow*nmodels][j+PCnow*nmodels]/lambdaP; //check this part..
	gsl_matrix_set (V11, i, j, V11dummy);
      }
    }


  for (j = 0; j < nmodels; j++)
    {
      double V21dummy = 1.; 
      for (k = 0; k < nparams; k++)
	{
	  double distance = designparams[j][k]-newparams[k];
	  distance = 4.*distance*distance; 
	  V21dummy *= pow(rho_w[numz*nparams+k][PCnow], distance); 
	}
      V21dummy /= lambda_w[numz*numPC+PCnow]; 
      gsl_matrix_set(V21, 0, j, V21dummy); 
    }

  return(0); 
}

int invert_matrix(int size, gsl_matrix *A, gsl_matrix *A_inv)
{

  gsl_matrix *Adummy = gsl_matrix_alloc(size,size);
  gsl_matrix *Adummy2 = gsl_matrix_alloc(size,size);
  
  gsl_matrix_memcpy (Adummy, A); 
  /* gsl_linalg_cholesky_decomp (Adummy);  // Adummy will be overwritten  */
  /* gsl_linalg_cholesky_invert (Adummy);  // This is just the inverse of the lower half of the decomposed matrix */
  gsl_vector *work = gsl_vector_alloc(size);
  gsl_matrix *V = gsl_matrix_alloc(size,size); // V is untransposed
  gsl_vector *Sdiag = gsl_vector_alloc(size); 
  gsl_matrix *S = gsl_matrix_alloc(size,size); 

  gsl_linalg_SV_decomp(Adummy, V, Sdiag, work); //A_inv is actually U: A = USV', Adummy becomes U on output;
  int i,j;
  double reslt = 0;

  gsl_matrix_set_zero(S);
  for (i=0; i < size; i++)
    {
      double reslt = gsl_vector_get(Sdiag,i);
      gsl_matrix_set(S,i,i, 1./reslt);
    }

  gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0,S, Adummy, 0.0, Adummy2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0,V, Adummy2, 0.0, A_inv);
	  
  return(0); 
}


