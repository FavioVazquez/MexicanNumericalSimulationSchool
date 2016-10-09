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
 *
 *  Modified by Irvin Martínez and Favio Vazquez on 10/07/2016
 *  
 *  This program was modified during the Mexican Numerical Simulation Schools on Salman Habib (LANL) project. 
 *
 *
 */
 
 
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"

main(int argc, char **argv) {
    int i,j,type=2, writeout=0, var=0, mod=0;
    double xstar[7], ystar[2*582], stuff[4], xstarcmb[6];
    FILE *fp;

    //char fname_nl[256], fname_nl1[256], fname_nl2[256], fname_nl3[256], fname_nl4[256], fname_nl5[256], fname_nl6[256], fname_nl7[256],fname_nl8[256],fname_nl9[256];
    char fname_nl[256], fname_nl1[256] = "output1.txt", fname_nl2[256] = "output2.txt", fname_nl3[256] = "output3.txt", fname_nl4[256] = "output4.txt", fname_nl5[256] = "output5.txt", fname_nl6[256] = "output6.txt", fname_nl7[256] = "output7.txt", fname_nl8[256] = "output8.txt", fname_nl9[256] = "output9.txt", fname_nl10[256] = "output10.txt";
    int cmbh=0;
    int numvar, count;

   

    printf("Which model do you want to use?\n");
    printf("(1) WMAP7, (2) WMAP9, (3) Planck\n");
    scanf("%i",&mod);
    printf("\n");


/*
    printf("Enter the general filename for output : ");
    scanf("%s",gfname_nl);
    printf("Output will be written to: %s.\n", gfname_nl);
*/
    //xstar[0] - Omega_m*h^2, xstar[1] - Omega_b*h^2, xstar[2] - n_s, xstar[3] - sigma_8 , xstar[4] - w, xstar[5] - H_0, xstar[6] - z
    
    
    //WMAP7 data
    if (mod == 1)
    {
        printf("You've choosen number WMAP7");
        printf("\n");
        xstar[0] = 0.02258;
        xstar[1] = 0.1334;
        xstar[2] = 0.963;
        xstar[3] = 70.8;
        xstar[4] = -1;
        xstar[5] = 0.801;
        xstar[6] = 0;
        printf("Your output files are output$.txt\n");
    }
    //WMAP9 data
    else if (mod == 2)
    {
        printf("You've choosen number WMAP9");
        printf("\n");
        xstar[0] = 0.02264;
        xstar[1] = 0.1364;
        xstar[2] = 0.972;
        xstar[3] = 70;
        xstar[4] = -1;
        xstar[5] = 0.821;
        xstar[6] = 0;  
        printf("Your output files are output$.txt\n");
    }
    //PLANCK data
    else if (mod == 3)
    {
        printf("You've choosen number Planck");
        printf("\n");
        xstar[0] = 0.02227;
        xstar[1] = 0.1413;
        xstar[2] = 0.9681;
        xstar[3] = 67.9;
        xstar[4] = -1;
        xstar[5] = 0.8154;
        xstar[6] = 0;  
        printf("Your output files are output$.txt\n");
    }
     else
    {
        printf("This is not valid\n");
        exit(0);
    }
    
    printf("Those are the parameters\n");
    printf("Omega_b*h^2 = %g\n", xstar[0]);
    printf("Omega_m*h^2 = %g\n", xstar[1]);
    printf("n_s = %g\n", xstar[2]);
    printf("H_0 = %g\n", xstar[3]);
    printf("w = %g\n", xstar[4]);
    printf("sigma_8 = %g\n", xstar[5]);
    printf("z = %g\n", xstar[6]);
    printf("\n");

    printf("Wich parameters do you want to vary?\n");
    printf("(1) Omega_b*h^2, (2) Omega_m*h^2, (3) n_s, (4)sigma_8, (5)w\n");
    scanf("%i",&var);
    printf("You've selected %i\n", var);

    printf("Enter the number of variations you want to do (up to ten): \n");
    scanf("%i", &numvar);
    printf("You're going to vary %i times\n", numvar);

    for ( i = 1; i <= numvar; ++i)
    {
        //
       if (i == 1)
       {
           /*  
            printf("Enter filename for output : ");
            scanf("%s",fname_nl);
            printf("Output will be written to: %s.\n", fname_nl);
            */

            strcpy(fname_nl, fname_nl1);
            printf("Output will be written to: %s.\n", fname_nl);


            if (var == 1)
            {   
                printf("0.0215 < omega_b*h^2 < 0.0235\n");
                printf("Enter omega_b (= Omega_b*h^2): \n");
                scanf("%lf",&xstar[0]);
                printf("%g\n",xstar[0]);
            }

            else if (var == 2)
            {
                printf("0.120 < omega_m*h^2 < 0.155\n");
                printf("Enter omega_m (= Omega_m*h^2): \n");
                scanf("%lf",&xstar[1]);
                printf("%g\n",xstar[1]);
            }

            else if (var == 3)
            {
                printf("0.85 < n_s < 1.05\n");
                printf("Enter n_s: \n");
                scanf("%lf",&xstar[2]);
                printf("%g\n",xstar[2]);
            }

            else if (var == 4)
            {
                printf("0.61 < sigma_8 < 0.9\n");
                printf("Enter sigma_8: \n");
                scanf("%lf",&xstar[5]);
                printf("%g\n",xstar[5]);
            }

            else if (var == 5)
            {
                printf("-1.30 < w < -0.70\n");
                printf("Enter w: \n");
                scanf("%lf",&xstar[4]);
                printf("%g\n",xstar[4]);
            }

            else
            {
                printf("This is not valid\n");
                exit(0);
            }


            type = 2;
    
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

        else if (i == 2)
       {
           /*  
            printf("Enter filename for output : ");
            scanf("%s",fname_nl);
            printf("Output will be written to: %s.\n", fname_nl);
            */

            strcpy(fname_nl, fname_nl2);
            printf("Output will be written to: %s.\n", fname_nl);


            if (var == 1)
            {   
                printf("Enter omega_b (= Omega_b*h^2): \n");
                scanf("%lf",&xstar[0]);
                printf("%g\n",xstar[0]);
            }

            else if (var == 2)
            {
                printf("Enter omega_m (= Omega_m*h^2): \n");
                scanf("%lf",&xstar[1]);
                printf("%g\n",xstar[1]);
            }

            else if (var == 3)
            {
                printf("Enter n_s: \n");
                scanf("%lf",&xstar[2]);
                printf("%g\n",xstar[2]);
            }

            else if (var == 4)
            {
                printf("Enter sigma_8: \n");
                scanf("%lf",&xstar[5]);
                printf("%g\n",xstar[5]);
            }

            else if (var == 5)
            {
                printf("Enter w: \n");
                scanf("%lf",&xstar[4]);
                printf("%g\n",xstar[4]);
            }

            else
            {
                printf("This is not valid\n");
                exit(0);
            }


            type = 2;
    
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

        else if (i == 3)
       {
           /*  
            printf("Enter filename for output : ");
            scanf("%s",fname_nl);
            printf("Output will be written to: %s.\n", fname_nl);
            */

            strcpy(fname_nl, fname_nl3);
            printf("Output will be written to: %s.\n", fname_nl);



            if (var == 1)
            {   
                printf("Enter omega_b (= Omega_b*h^2): \n");
                scanf("%lf",&xstar[0]);
                printf("%g\n",xstar[0]);
            }

            else if (var == 2)
            {
                printf("Enter omega_m (= Omega_m*h^2): \n");
                scanf("%lf",&xstar[1]);
                printf("%g\n",xstar[1]);
            }

            else if (var == 3)
            {
                printf("Enter n_s: \n");
                scanf("%lf",&xstar[2]);
                printf("%g\n",xstar[2]);
            }

            else if (var == 4)
            {
                printf("Enter sigma_8: \n");
                scanf("%lf",&xstar[5]);
                printf("%g\n",xstar[5]);
            }

            else if (var == 5)
            {
                printf("Enter w: \n");
                scanf("%lf",&xstar[4]);
                printf("%g\n",xstar[4]);
            }

            else
            {
                printf("This is not valid\n");
                exit(0);
            }


            type = 2;
    
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

        else if (i == 4)
       {
           /*  
            printf("Enter filename for output : ");
            scanf("%s",fname_nl);
            printf("Output will be written to: %s.\n", fname_nl);
            */

            strcpy(fname_nl, fname_nl4);
            printf("Output will be written to: %s.\n", fname_nl);



            if (var == 1)
            {   
                printf("Enter omega_b (= Omega_b*h^2): \n");
                scanf("%lf",&xstar[0]);
                printf("%g\n",xstar[0]);
            }

            else if (var == 2)
            {
                printf("Enter omega_m (= Omega_m*h^2): \n");
                scanf("%lf",&xstar[1]);
                printf("%g\n",xstar[1]);
            }

            else if (var == 3)
            {
                printf("Enter n_s: \n");
                scanf("%lf",&xstar[2]);
                printf("%g\n",xstar[2]);
            }

            else if (var == 4)
            {
                printf("Enter sigma_8: \n");
                scanf("%lf",&xstar[5]);
                printf("%g\n",xstar[5]);
            }

            else if (var == 5)
            {
                printf("Enter w: \n");
                scanf("%lf",&xstar[4]);
                printf("%g\n",xstar[4]);
            }

            else
            {
                printf("This is not valid\n");
                exit(0);
            }


            type = 2;
    
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

        else if (i == 5)
       {
           /*  
            printf("Enter filename for output : ");
            scanf("%s",fname_nl);
            printf("Output will be written to: %s.\n", fname_nl);
            */

            strcpy(fname_nl, fname_nl5);
            printf("Output will be written to: %s.\n", fname_nl);



            if (var == 1)
            {   
                printf("Enter omega_b (= Omega_b*h^2): \n");
                scanf("%lf",&xstar[0]);
                printf("%g\n",xstar[0]);
            }

            else if (var == 2)
            {
                printf("Enter omega_m (= Omega_m*h^2): \n");
                scanf("%lf",&xstar[1]);
                printf("%g\n",xstar[1]);
            }

            else if (var == 3)
            {
                printf("Enter n_s: \n");
                scanf("%lf",&xstar[2]);
                printf("%g\n",xstar[2]);
            }

            else if (var == 4)
            {
                printf("Enter sigma_8: \n");
                scanf("%lf",&xstar[5]);
                printf("%g\n",xstar[5]);
            }
            else if (var == 5)
            {
                printf("Enter w: \n");
                scanf("%lf",&xstar[4]);
                printf("%g\n",xstar[4]);
            }
            else
            {
                printf("This is not valid\n");
                exit(0);
            }


            type = 2;
    
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

        else if (i == 6)
       {
           /*  
            printf("Enter filename for output : ");
            scanf("%s",fname_nl);
            printf("Output will be written to: %s.\n", fname_nl);
            */

            strcpy(fname_nl, fname_nl6);
            printf("Output will be written to: %s.\n", fname_nl);



            if (var == 1)
            {   
                printf("Enter omega_b (= Omega_b*h^2): \n");
                scanf("%lf",&xstar[0]);
                printf("%g\n",xstar[0]);
            }

            else if (var == 2)
            {
                printf("Enter omega_m (= Omega_m*h^2): \n");
                scanf("%lf",&xstar[1]);
                printf("%g\n",xstar[1]);
            }

            else if (var == 3)
            {
                printf("Enter n_s: \n");
                scanf("%lf",&xstar[2]);
                printf("%g\n",xstar[2]);
            }

            else if (var == 4)
            {
                printf("Enter sigma_8: \n");
                scanf("%lf",&xstar[5]);
                printf("%g\n",xstar[5]);
            }

            else if (var == 5)
            {
                printf("Enter w: \n");
                scanf("%lf",&xstar[4]);
                printf("%g\n",xstar[4]);
            }

            else
            {
                printf("This is not valid\n");
                exit(0);
            }


            type = 2;
    
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

        else if (i == 7)
       {
           /*  
            printf("Enter filename for output : ");
            scanf("%s",fname_nl);
            printf("Output will be written to: %s.\n", fname_nl);
            */

            strcpy(fname_nl, fname_nl7);
            printf("Output will be written to: %s.\n", fname_nl);



            if (var == 1)
            {   
                printf("Enter omega_b (= Omega_b*h^2): \n");
                scanf("%lf",&xstar[0]);
                printf("%g\n",xstar[0]);
            }

            else if (var == 2)
            {
                printf("Enter omega_m (= Omega_m*h^2): \n");
                scanf("%lf",&xstar[1]);
                printf("%g\n",xstar[1]);
            }

            else if (var == 3)
            {
                printf("Enter n_s: \n");
                scanf("%lf",&xstar[2]);
                printf("%g\n",xstar[2]);
            }

            else if (var == 4)
            {
                printf("Enter sigma_8: \n");
                scanf("%lf",&xstar[5]);
                printf("%g\n",xstar[5]);
            }

            else if (var == 5)
            {
                printf("Enter w: \n");
                scanf("%lf",&xstar[4]);
                printf("%g\n",xstar[4]);
            }

            else
            {
                printf("This is not valid\n");
                exit(0);
            }


            type = 2;
    
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

        else if (i == 8)
       {
           /*  
            printf("Enter filename for output : ");
            scanf("%s",fname_nl);
            printf("Output will be written to: %s.\n", fname_nl);
            */

            strcpy(fname_nl, fname_nl8);
            printf("Output will be written to: %s.\n", fname_nl);



            if (var == 1)
            {   
                printf("Enter omega_b (= Omega_b*h^2): \n");
                scanf("%lf",&xstar[0]);
                printf("%g\n",xstar[0]);
            }

            else if (var == 2)
            {
                printf("Enter omega_m (= Omega_m*h^2): \n");
                scanf("%lf",&xstar[1]);
                printf("%g\n",xstar[1]);
            }

            else if (var == 3)
            {
                printf("Enter n_s: \n");
                scanf("%lf",&xstar[2]);
                printf("%g\n",xstar[2]);
            }

            else if (var == 4)
            {
                printf("Enter sigma_8: \n");
                scanf("%lf",&xstar[5]);
                printf("%g\n",xstar[5]);
            }

            else if (var == 5)
            {
                printf("Enter w: \n");
                scanf("%lf",&xstar[4]);
                printf("%g\n",xstar[4]);
            }

            else
            {
                printf("This is not valid\n");
                exit(0);
            }


            type = 2;
    
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

        else if (i == 9)
       {
           /*  
            printf("Enter filename for output : ");
            scanf("%s",fname_nl);
            printf("Output will be written to: %s.\n", fname_nl);
            */

            strcpy(fname_nl, fname_nl9);
            printf("Output will be written to: %s.\n", fname_nl);



            if (var == 1)
            {   
                printf("Enter omega_b (= Omega_b*h^2): \n");
                scanf("%lf",&xstar[0]);
                printf("%g\n",xstar[0]);
            }

            else if (var == 2)
            {
                printf("Enter omega_m (= Omega_m*h^2): \n");
                scanf("%lf",&xstar[1]);
                printf("%g\n",xstar[1]);
            }

            else if (var == 3)
            {
                printf("Enter n_s: \n");
                scanf("%lf",&xstar[2]);
                printf("%g\n",xstar[2]);
            }

            else if (var == 4)
            {
                printf("Enter sigma_8: \n");
                scanf("%lf",&xstar[5]);
                printf("%g\n",xstar[5]);
            }

            else if (var == 5)
            {
                printf("Enter w: \n");
                scanf("%lf",&xstar[4]);
                printf("%g\n",xstar[4]);
            }

            else
            {
                printf("This is not valid\n");
                exit(0);
            }


            type = 2;
    
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

        else if (i == 10)
       {
           /*  
            printf("Enter filename for output : ");
            scanf("%s",fname_nl);
            printf("Output will be written to: %s.\n", fname_nl);
            */

            strcpy(fname_nl, fname_nl10);
            printf("Output will be written to: %s.\n", fname_nl);



            if (var == 1)
            {   
                printf("Enter omega_b (= Omega_b*h^2): \n");
                scanf("%lf",&xstar[0]);
                printf("%g\n",xstar[0]);
            }

            else if (var == 2)
            {
                printf("Enter omega_m (= Omega_m*h^2): \n");
                scanf("%lf",&xstar[1]);
                printf("%g\n",xstar[1]);
            }

            else if (var == 3)
            {
                printf("Enter n_s: \n");
                scanf("%lf",&xstar[2]);
                printf("%g\n",xstar[2]);
            }

            else if (var == 4)
            {
                printf("Enter sigma_8: \n");
                scanf("%lf",&xstar[5]);
                printf("%g\n",xstar[5]);
            }

            else if (var == 5)
            {
                printf("Enter w: \n");
                scanf("%lf",&xstar[4]);
                printf("%g\n",xstar[4]);
            }

            else
            {
                printf("This is not valid\n");
                exit(0);
            }


            type = 2;
    
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

    }

    

    //printf("Enter filename for output : ");
    //scanf("%s",fname_nl);
    

    /*
    if (var == 1)
    {
        printf("\n");
    }
    */
    /*
    if (var == 2)
    {
        printf("Enter omega_b (= Omega_b*h^2): ");
        scanf("%lf",&xstar[0]);
        printf("%g\n",xstar[0]);
    }

    else if (var == 3)
    {
        printf("Enter omega_m (= Omega_m*h^2): ");
        scanf("%lf",&xstar[1]);
        printf("%g\n",xstar[1]);
    }

    else if (var == 4)
    {
        printf("Enter n_s: ");
        scanf("%lf",&xstar[2]);
        printf("%g\n",xstar[2]);

    }

    else if (var == 5)
    {
        printf("Enter sigma_8: ");
        scanf("%lf",&xstar[5]);
        printf("%g\n",xstar[5]);
    }
    else
    {
        printf("This is not valid\n");
        exit(0);
    }


    type = 2;
    
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
    */
}
