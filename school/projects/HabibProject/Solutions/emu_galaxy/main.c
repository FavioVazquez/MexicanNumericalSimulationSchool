#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "logk.h"
#include "params.h"
#include <math.h>

int main(int argc, char**argv)
{
  float min_design[5] = {12.9, 13.5, 0.5, 0.5, 0.5};
  float max_design[5] = {14.0, 15.0, 1.2, 1.5, 1.5};
  double newparams[nparams]; 
  char inputs[256]; 
  char paramnames[5][20]; 
  double outputredshift; 

  if (argc < 3)
    {
      fprintf(stderr, "Some input files are missing.\n The correct arguments are ./emu.out params.ini output.txt\n"); 
      exit(1); 
    }

  sprintf(paramnames[0], "M_cut"); 
  sprintf(paramnames[1], "M1"); 
  sprintf(paramnames[2], "sigma"); 
  sprintf(paramnames[3], "kappa"); 
  sprintf(paramnames[4], "alpha"); 


  //read .ini file
  sprintf(inputs, "%s", argv[1]); 
  FILE *fpinputs = fopen(inputs,"r");
  if (fpinputs==NULL)
    {
      fprintf(stderr, "I can't find this parameter file: %s.\nExiting now. \n", argv[1]); 
      exit(1);
    }
  float buf; 
  char line[1000];
  int n=0;
  int c; 
  int which_2pt; 
  while (fgets(line, sizeof(line), fpinputs) || n < 7) {
    if (*line == '#') continue; // ignore comment line 
    if (sscanf(line, "%f", &buf) != 1) 
      {
	// one parameter per line only! 
	if (strlen(line)!=1) // ignore blank lines
	  {
	    fprintf(stderr, "I don't know how to read your file format.\nPlease modify it to contain one parameter per line only.\nFor comments, please use a '#' symbol at the start of the line.\n"); 
	    exit(1); 
	  }

      }
    else 
      {
	if (n==5)
	  outputredshift = buf; 
	if (n==6)
	  which_2pt = buf; 
	else 
	  newparams[n] = buf; 
	n++; 
      }
  }

  if (n < 7)
    {
      fprintf(stderr, "I didn't read enough input parameters. \nSome lines may be missing from your input file.\n");  
      exit(1);
    }

  fclose(fpinputs);

  // check if parameters are within emulation range: 
  for (n = 0; n < 5; n++)
    {
      if (newparams[n] < min_design[n] || newparams[n] > max_design[n])
  	{
	    fprintf(stderr, "%s = %lf is outside of the emulation range: %f -- %f. \nPlease adjust your parameters accordingly.\n", paramnames[n], newparams[n], min_design[n], max_design[n]);
  	  fflush;
  	  exit(1);
  	}
    }
  if (outputredshift > 1 || outputredshift < 0)
    {
      fprintf(stderr, "%s = %f is outside of the emulation range: %f -- %f. \nPlease adjust your parameters accordingly.\n", paramnames[n], newparams[n], min_design[n], max_design[n]);
      fflush;
      exit(1);
    }
  

  char outputfile[256]; 
  sprintf(outputfile, "%s", argv[2]); 
  FILE *fp = fopen(outputfile,"w");
  if (fp==NULL)
    {
      fprintf(stderr, "I can't open this file for writing: %s.\nExiting now. \n", argv[2]); 
      exit(1);
    }


  double *output_pk = malloc(nk*sizeof(double)); // this must be at least of length nk elements
  double h; 

  //write parameters to output file as a record of used parameters
  for (n = 0; n < 5; n++)
    fprintf(fp, "# %s = %f\n", paramnames[n], newparams[n]);

  fprintf(fp, "# z = %f\n", outputredshift);

  //Normalise inputs:
  for (n = 0; n < 5; n++)
    {
      newparams[n]-=min_design[n];
      newparams[n]/=(max_design[n]-min_design[n]);

    }

  // Now do the emulation! 
  emu(newparams, outputredshift, output_pk); 

  double *k_unlogged = malloc(nk*sizeof(double)); 
  for (n = 0; n < nk; n++)
    k_unlogged[n] = pow(10.,logk[n]);

  if (which_2pt ==1)
    {
      fprintf(fp, "# k [Mpc^-1]    P(k) [Mpc^3] \n");
      for (n = 0; n < nk; n++)
	{
	  if (k_unlogged[n] > 0.01)  // only print out the k modes that are emulated
	    fprintf(fp,"%f %f\n", k_unlogged[n], output_pk[n]);  
	}
    }
  else if (which_2pt == 2)
    {
      int extra_pts = 51; // no of points to extend in either direction                                                                                                           
      int npts_ext = nk+2*extra_pts-2;

      double  *k_ext = malloc(npts_ext*sizeof(double));
      double *Pk_ext = malloc(npts_ext*sizeof(double));

      pade(nk, k_unlogged, output_pk, npts_ext, k_ext, Pk_ext);

      // Convert to xi(r) using Ogata (2005)                                                                                                                
      double *r_out = malloc(200*sizeof(double));
      double *xi_out = malloc(200*sizeof(double));

      pk_to_xi(npts_ext, k_ext, Pk_ext, r_out, xi_out);

      fprintf(fp, "# r [Mpc]    xi(r) \n");
      for (n = 0; n < 200; n++)
	  fprintf(fp,"%f %f\n", r_out[n], xi_out[n]);

      free(k_ext); free(Pk_ext);
      free(r_out); free(xi_out);

    }
  else
    {
      fprintf(stderr,"Sorry, I didn't understand your parameters. Please select 1 for P(k) or 2 for xi(r).\n"); 
    }

  fclose(fp);
  // undo normalisation (to avoid confusion)
  for (n = 0; n < 5; n++)
    {
      newparams[n]*=(max_design[n]-min_design[n]);
      newparams[n]+=min_design[n];
    }

  free(k_unlogged); 
  free(output_pk); 
  return(0);
}



