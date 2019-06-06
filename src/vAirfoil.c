//
//  vAirfoil.c
//  
//
//
//
//

#include "vAirfoil.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void readarfsize(char *WngArf,int ArfXUNumber[],int ArfXLNumber[])
{
    /* Read size of data in airfoil file */
    FILE *fp2;
    int cond,npoint;
    char line[110],code[8];
    char* arfname=malloc(strlen(WngArf) + strlen(".arf") + 2);
    
    /* Append .arf at the end of the airfoil name to obtain the filename */
    sscanf(WngArf,"%s",arfname);
    strcat(arfname, ".arf");
    
    /* Open airfoil file */
    fp2 = fopen(arfname,"r");
    if(fp2 == NULL) {
        fprintf(stderr,"Error:  Could not open file %s\n",arfname);
        exit(1);
    }
    // Deallocate memory
    free(arfname);
    /* Scan the file until lines GMD401 and GMD402 to read the number of points on the two surfaces */
    cond=0;
    while (cond == 0){
        fgets(line, 110, fp2);
        /* Number of points on upper surface */
        if ( strncmp("GMD401",line,6) == 0 ){
            sscanf(line,"%s %i",code,ArfXUNumber);
        }
        /* Number of points on lower surface */
        if ( strncmp("GMD402",line,6) == 0 ){
            sscanf(line,"%s %i",code,ArfXLNumber);
            cond=1;
        }
    }
    /* Close airfoil file */
    fclose(fp2);
}

void readarf(char *WngArf,double *ArfXU, double *ArfYU, double *ArfXL,double *ArfYL,int ArfXUNumber[],int ArfXLNumber[])
{
    /* Read airfoil file data */
    FILE *fp2;
    int cond,npoint;
    char line[110],code[8];
    char* arfname=malloc(strlen(WngArf) + strlen(".arf") + 2);
    
    /* Append .arf at the end of the airfoil name to obtain the filename */
    sscanf(WngArf,"%s",arfname);
    strcat(arfname, ".arf");
    
    /* Open airfoil file */
    fp2 = fopen(arfname,"r");
    if(fp2 == NULL) {
        fprintf(stderr,"Error:  Could not open file %s\n",arfname);
        exit(1);
    }
    // Deallocate memory
    free(arfname);
    /* Scan the file until lines GMD15 to read the upper and lower surface data */
    fgets(line, 110, fp2);
    cond=0;
    npoint=-1; /* Set up counter to count how many datapoints have been read */
    while (cond == 0){
        fgets(line, 110, fp2);
        if ( strncmp("GM15",line,4) == 0 ){
            npoint++;
            /* If there are more points on the upper surface than the lower surface */
            if (ArfXUNumber[0] > ArfXLNumber[0]){
                if (npoint < ArfXLNumber[0]){
                    sscanf(line,"%s %lf %lf %lf %lf",code,(ArfXU+npoint),(ArfYU+npoint),(ArfXL+npoint),(ArfYL+npoint));
                }
                if (npoint >= ArfXLNumber[0] && npoint < ArfXUNumber[0]){
                    sscanf(line,"%s %lf %lf",code,(ArfXU+npoint),(ArfYU+npoint));
                }
                /* If there are more points on the lower surface than the upper surface */
            }else if (ArfXUNumber[0] < ArfXLNumber[0]){
                if (npoint < ArfXUNumber[0]){
                    sscanf(line,"%s %lf %lf %lf %lf",code,(ArfXU+npoint),(ArfYU+npoint),(ArfXL+npoint),(ArfYL+npoint));
                }
                if (npoint >= ArfXUNumber[0] && npoint < ArfXLNumber[0]){
                    sscanf(line,"%s %lf %lf",code,(ArfXL+npoint),(ArfYL+npoint));
                }
                /* If the number of points on the upper and lower surface is the same */
            }else{
                if (npoint < ArfXUNumber[0]){
                    sscanf(line,"%s %lf %lf %lf %lf",code,(ArfXU+npoint),(ArfYU+npoint),(ArfXL+npoint),(ArfYL+npoint));
                }
            }
        }
        /* Stop condition for the while loop */
        if ( strncmp("GMD6",line,4) == 0 ){
            cond=1;
        }
    }
    /* Close airfoil file */
    fclose(fp2);
}

void interparf(double *ycamber, double *ArfXU, double *ArfYU, double *ArfXL, double *ArfYL, int ArfXUNumber[], int ArfXLNumber[], double *xpline, int mp1)
{
    /* Use interpolation to obtain airfoil camber data at the desired x-coordinates */
    int i,j;
    double slope;
    
    /* Carry out linear interpolation to find ycamber at every element of xpline */
    for (i=0;i<mp1;i++){
        j=0;
        /* Upper surface interpolation */
        while (*(ArfXU+j)/100.0 <= *(xpline+i) && j <= ArfXUNumber[0]-1){
            j++;
        }
        slope=(*(ArfYU+j)- *(ArfYU+j-1))/(*(ArfXU+j)- *(ArfXU+j-1));
        /* Store the upper surface y-coordinate in ycamber */
        *(ycamber+i)=(slope*(*(xpline+i)- *(ArfXU+j-1)/100)+ *(ArfYU+j-1)/100.0);
        
        /* Lower surface interpolation */
        j=0;
        while (*(ArfXL+j)/100.0 <= *(xpline+i) && j <= ArfXLNumber[0]-1){
            j++;
        }
        slope=(*(ArfYL+j)- *(ArfYL+j-1))/(*(ArfXL+j)- *(ArfXL+j-1));
        /* Add the lower surface y-coordinate in ycamber and divide by two */
        *(ycamber+i)=(*(ycamber+i)+(slope*(*(xpline+i)- *(ArfXL+j-1)/100)+ *(ArfYL+j-1)/100.0))/2.0;
    }
}

void nacafourfivedigit(double *xpline, double *ycamber, int mp1, char WngArf[])
{
    /* Calculate the camber line of a 4- or 5-digit NACA series airfoil */
    int i;
    double xhere,m,p,k1,meanline;
    char *end;
    long number = strtol(WngArf, &end, 10); /* Get airfoil name */
    
    /* Check if it's a 4- or 5-digit series */
    if (number > 9999){
        /* NACA 5-digit series */
        meanline=floor(((double)number)/100.0);
        /* Choose a mean line */
        if (meanline == 210.0){
            p=0.05;
            m=0.0580;
            k1=361.4;
        }else if (meanline == 220.0){
            p=0.1;
            m=0.1260;
            k1=51.64;
        }else if (meanline == 230.0){
            p=0.15;
            m=0.2025;
            k1=15.957;
        }else if (meanline == 240.0){
            p=0.2;
            m=0.29;
            k1=6.643;
        }else if (meanline == 250.0){
            p=0.25;
            m=0.391;
            k1=3.23;
        }
        /* Calculate camber line */
        for (i=0;i<mp1;i++){
            xhere=*(xpline+i);
            if (xhere <= m){
                *(ycamber+i)=1/6.0*k1*(pow(xhere,3)-3*m*pow(xhere,2)+m*m*(3-m)*xhere);
            }else{
                *(ycamber+i)=1/6.0*k1*m*m*m*(1-xhere);
            }
        }
    }else{
        /* NACA 4-digit series */
        /* Choose maximum camber m and maximum camber position p */
        m=floor(((double)number)/1000.0)/100;
        p=floor((((double)number)-m*100000.0)/100.0)/10.0;
        /* Calculate camber line */
        for (i=0;i<mp1;i++){
            if (m == 0){
                *(ycamber+i)=0;
            }else{
                xhere=*(xpline+i);
                if (xhere <= p){
                    *(ycamber+i)=m/(p*p)*(2.0*p*xhere-xhere*xhere);
                }else{
                    *(ycamber+i)=m/(pow((1.0-p),2))*((1-2.0*p)+2.0*p*xhere-xhere*xhere);
                }
            }
        }
    }
}

void treatarfwgl(char line[], double *xpline, double *ycamber, int mp1)
{
    /* Carries out the treatment of the airfoil on the winglet */
    int ArfXUNumber[1],ArfXLNumber[1],j;
    double *ArfXU, *ArfYU, *ArfXL, *ArfYL;
    char code[8], *WglArf, *nacanumber;
    char *token = NULL;
    
    token = strtok(line, "\t,\n");
    WglArf = strtok(NULL, "\t,\n");
    
    /* Inboard airfoil */
    if ( strncmp("NACA",WglArf,4) == 0 ){
        token = strtok(WglArf, " ,\n");
        nacanumber=strtok(NULL, " ,\n");
        nacafourfivedigit(xpline,ycamber,mp1,nacanumber);
    }else{
        readarfsize(WglArf,ArfXUNumber,ArfXLNumber); /* Read airfoil data size */
        /* Set up vectors to store airfoil data */
        ArfXU = (double *)malloc(sizeof(double)*ArfXUNumber[0]);
        ArfYU = (double *)malloc(sizeof(double)*ArfXUNumber[0]);
        ArfXL = (double *)malloc(sizeof(double)*ArfXLNumber[0]);
        ArfYL = (double *)malloc(sizeof(double)*ArfXLNumber[0]);
        /* Read airfoil data */
        readarf(WglArf,ArfXU,ArfYU,ArfXL,ArfYL,ArfXUNumber,ArfXLNumber);
        /* Interpolate airfoil at the desired coordinates xpline */
        interparf(ycamber,ArfXU,ArfYU,ArfXL,ArfYL,ArfXUNumber,ArfXLNumber,xpline,mp1);
        // Deallocate
        free(ArfXU);
        free(ArfYU);
        free(ArfXL);
        free(ArfYL);
    }
}

void treatarf(char line[], double *xpline, double *ycamber, double *ycamberall, int mp1, int iTS)
{
    /* Carries out the treatment of the airfoils on each end of each trapezoidal section */
    int ArfXUNumber[1],ArfXLNumber[1],j;
    double *ArfXU, *ArfYU, *ArfXL, *ArfYL;
    char code[8], *WngTSRtArf, *WngTSTpArf, *nacanumber;
    char *token = NULL;
    
    token = strtok(line, "\t,\n");
    WngTSRtArf = strtok(NULL, "\t,\n");
    WngTSTpArf = strtok(NULL, "\t,\n");
    
    /* Inboard airfoil */
    if ( strncmp("NACA",WngTSRtArf,4) == 0 ){
        token = strtok(WngTSRtArf, " ,\n");
        nacanumber=strtok(NULL, " ,\n");
        nacafourfivedigit(xpline,ycamber,mp1,nacanumber);
        for (j=0;j<mp1;j++){
            *(ycamberall+j+2*iTS*mp1)=*(ycamber+j);
        }
    }else{
        readarfsize(WngTSRtArf,ArfXUNumber,ArfXLNumber); /* Read airfoil data size */
        /* Set up vectors to store airfoil data */
        ArfXU = (double *)malloc(sizeof(double)*ArfXUNumber[0]);
        ArfYU = (double *)malloc(sizeof(double)*ArfXUNumber[0]);
        ArfXL = (double *)malloc(sizeof(double)*ArfXLNumber[0]);
        ArfYL = (double *)malloc(sizeof(double)*ArfXLNumber[0]);
        /* Read airfoil data */
        readarf(WngTSRtArf,ArfXU,ArfYU,ArfXL,ArfYL,ArfXUNumber,ArfXLNumber);
        /* Interpolate airfoil at the desired coordinates xpline */
        interparf(ycamber,ArfXU,ArfYU,ArfXL,ArfYL,ArfXUNumber,ArfXLNumber,xpline,mp1);
        // Deallocate
        free(ArfXU);
        free(ArfYU);
        free(ArfXL);
        free(ArfYL);
        for (j=0;j<mp1;j++){
            *(ycamberall+j+2*iTS*mp1)=*(ycamber+j);
        }
    }
    
    /* Outboard airfoil */
    if ( strncmp("NACA",WngTSTpArf,4) == 0 ){
        token = strtok(WngTSTpArf, " ,\n");
        nacanumber=strtok(NULL, " ,\n");
        nacafourfivedigit(xpline,ycamber,mp1,nacanumber);
        for (j=0;j<mp1;j++){
            *(ycamberall+j+(2*iTS+1)*mp1)=*(ycamber+j);
        }
    }else{
        readarfsize(WngTSTpArf,ArfXUNumber,ArfXLNumber); /* Read airfoil data size */
        /* Set up vectors to store airfoil data */
        ArfXU = (double *)malloc(sizeof(double)*ArfXUNumber[0]);
        ArfYU = (double *)malloc(sizeof(double)*ArfXUNumber[0]);
        ArfXL = (double *)malloc(sizeof(double)*ArfXLNumber[0]);
        ArfYL = (double *)malloc(sizeof(double)*ArfXLNumber[0]);
        /* Read airfoil data */
        readarf(WngTSTpArf,ArfXU,ArfYU,ArfXL,ArfYL,ArfXUNumber,ArfXLNumber);
        /* Interpolate airfoil at the desired coordinates xpline */
        interparf(ycamber,ArfXU,ArfYU,ArfXL,ArfYL,ArfXUNumber,ArfXLNumber,xpline,mp1);
        // Deallocate
        free(ArfXU);
        free(ArfYU);
        free(ArfXL);
        free(ArfYL);
        for (j=0;j<mp1;j++){
            *(ycamberall+j+(2*iTS+1)*mp1)=*(ycamber+j);
        }
    }
}