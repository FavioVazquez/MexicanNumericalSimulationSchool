/*
 *  main.c
 *  
 *
 *  Created by Earl Lawrence on 9/17/09.
 *
 *  This program was prepared by Los Alamos National Security, LLC at Los Alamos National Laboratory (LANL) 
 *  under contract No. DE-AC52-06NA25396 with the U.S. Department of Energy (DOE). All rights in the program 
 *  are reserved by the DOE and Los Alamos National Security, LLC.  Permission is granted to the public to 
 *  copy and use this software without charge, provided that this Notice and any statement of authorship are 
 *  reproduced on all copies.  Neither the U.S. Government nor LANS makes any warranty, express or implied, 
 *  or assumes any liability or responsibility for the use of this software.
 *
 */
 
 
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"

main(int argc, char **argv) {
    int i,j,type=2, writeout=0;
    double xstar[7], ystar[2*582], stuff[4], xstarcmb[6];
    FILE *fp;

    char fname_nl[256];
    int cmbh=1;

// My code 

    int model = 0;
    
    FIRST_OPTION:
    printf("Enter model (1 = wmap7, 2=wmap9, 3=planck)\n");
    scanf("%i",&model);
    printf("You chose model %i\n", model);
    
    switch(model) {
	case 1:	
	    xstar[0]=0.02258;
	    xstar[1]=0.1334;
	    xstar[2]=0.963;
	    xstar[3]=70.8;
	    xstar[4]=-1;
	    xstar[5]=0.801;
	    xstar[6]=0;
            break;
        case 2:
	    xstar[0]=0.02264;
	    xstar[1]=0.1364;
	    xstar[2]=0.972;
	    xstar[3]=70.8;
	    xstar[4]=-1;
	    xstar[5]=0.821;
	    xstar[6]=0;
            break;
        case 3:
	    xstar[0]=0.02227;
	    xstar[1]=0.1413;
	    xstar[2]=0.9681;
	    xstar[3]=67.90;
	    xstar[4]=-1;
	    xstar[5]=0.8154;
	    xstar[6]=0;
            break;
        default:
	    printf("Choose a model between 1 and 3\n\n\n");
	    goto FIRST_OPTION;
    }

    printf("Enter filename for output : ");
    scanf("%s",fname_nl);
    printf("Output will be written to: %s.\n", fname_nl);
    
    printf("Will you be using h as derived from CMB constraints (0 for no, 1 for yes)?\n");
    scanf("%d", &cmbh);
/*
    printf("Enter omega_b (= Omega_b*h^2): ");
    scanf("%lf",&xstar[0]);
    printf("%g\n",xstar[0]);
    printf("Enter omega_m (= Omega_m*h^2): ");
    scanf("%lf",&xstar[1]);
    printf("%g\n",xstar[1]);
    printf("Enter n_s: ");
    scanf("%lf",&xstar[2]);
    printf("%g\n",xstar[2]);
    if(cmbh == 0) {
        printf("Enter H0: ");
        scanf("%lf",&xstar[3]);
        printf("%g\n",xstar[3]);
    }
    printf("Enter w: ");
    scanf("%lf",&xstar[4]);
    printf("%g\n",xstar[4]);
    printf("Enter sigma_8: ");
    scanf("%lf",&xstar[5]);
    printf("%g\n",xstar[5]);
    printf("Enter z: ");
    scanf("%lf",&xstar[6]);
    printf("%g\n",xstar[6]);
 */   

    int param = 0;

SECOND_OPTION:
    printf("What parameter do you want to change? (1 = Om_b, 2=Om_m,3=n_s,4=H0, 5=w, 6=sigma_8, 7=z)\n");
    scanf("%i",&param);
    printf("You chose parameter number %i\n", param);

    switch(param) {
        case 1:
              printf("Enter omega_b (= Omega_b*h^2): ");
              scanf("%lf",&xstar[0]);
              printf("%g\n",xstar[0]);
            break;
        case 2:
              printf("Enter omega_m (= Omega_m*h^2): ");
              scanf("%lf",&xstar[1]);
              printf("%g\n",xstar[1]);
            break;
        case 3:
              printf("Enter n_s: ");
              scanf("%lf",&xstar[2]);
              printf("%g\n",xstar[2]);
            break;
        case 4:
              printf("Enter H0: ");
              scanf("%lf",&xstar[3]);
              printf("%g\n",xstar[3]);
            break;
        case 5:
              printf("Enter w: ");
              scanf("%lf",&xstar[4]);
              printf("%g\n",xstar[4]);
            break;
        case 6:
              printf("Enter sigma_8: ");
              scanf("%lf",&xstar[5]);
              printf("%g\n",xstar[5]);
            break;
        case 7:
              printf("Enter z: ");
              scanf("%lf",&xstar[6]);
              printf("%g\n",xstar[6]);
            break;
        default:
	    printf("Choose a parameter between 1 and 7\n\n\n");
	    goto SECOND_OPTION;
    }

    printf("Enter output type (0: Delta^2/k^1.5; 1: Delta^2; 2: P(k)): ");
    scanf("%i",&type);
    printf("%i\n", type);
    
    if(cmbh == 1) {
        xstarcmb[0] = xstar[0];
        xstarcmb[1] = xstar[1];
        xstarcmb[2] = xstar[2];
        xstarcmb[3] = xstar[4];
        xstarcmb[4] = xstar[5];
        xstarcmb[5] = xstar[6];
        emu_noh(xstarcmb, ystar, &type);
        getH0fromCMB(xstarcmb, stuff);
        xstar[3] = 100.*stuff[3];
    } else {
        emu(xstar, ystar, &type);
    }

    // Write the nonlinear file
    if ((fp = fopen(fname_nl,"w"))==NULL) {
        printf("cannot open %s \n",fname_nl);
        exit(1);
    }
    
    fprintf(fp, "# Parameters:\n");
    fprintf(fp, "# omega_b = %f, omega_m = %f, n_s = %f, h = %f, w = %f, sigma_8 = %f\n", xstar[0], xstar[1], xstar[2], xstar[3], xstar[4], xstar[5]);
    fprintf(fp, "# z = %f\n", xstar[6]);
    fprintf(fp, "#\n");
    fprintf(fp, "# k[1/Mpc] ");
    
    switch(type) {
        default:
            fprintf(fp, "# Delta^2 / k^1.5:\n");
            break;
        case 1:
            fprintf(fp, "# Delta^2:\n");
            break;
        case 2:
            fprintf(fp, "# P(k):\n");
            break;
    }
    
    for(j=0; j<582; j++) {
        fprintf(fp ,"%e %e \n", ystar[j], ystar[582+j]);
    }
    fclose(fp);
}
