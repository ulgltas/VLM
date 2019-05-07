//
//  vInput.c
//
//  Data input
//
//

#include "vInput.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void importInputFile(char *Infile, double *UVW, double *rho, double *aoa, double *yaw, int *m, int *mht, int *mvt, int *n, int *nht, int *nvt,
                    int *ntimes, double *timestep_denom, int *freewake, double *delta, double *beta, double *eta, double *zeta,
                    char *Wngfile, char *HTailfile, char *VTailfile)
{
    char line[110], code[8];
    double Q, aoadegrees, yawdegrees, deltadegrees, betadegrees, etadegrees, zetadegrees;

    printf("Reading parameters from %s\n",Infile);

    FILE *ifp = fopen(Infile,"r");
    if(ifp == NULL) {
        fprintf(stderr,"Error:  Could not open %s\n",Infile);
        exit(1);
    }


    while (fgets(line, 110, ifp) != NULL)
    {
        if (strncmp("INP1", line, 4) == 0)
        {
            sscanf(line, "%s\t%lf\t%lf\t%lf\t%lf", code, &Q, rho, &aoadegrees, &yawdegrees);
            printf("%f %f %f %f\n", Q, *rho, aoadegrees, yawdegrees);
        }
        if (strncmp("INP2", line, 4) == 0)
        {
            sscanf(line, "%s\t%i\t%i\t%i", code, m, mht, mvt);
            printf("%i %i %i\n", *m, *mht, *mvt);
        }
        if (strncmp("INP3", line, 4) == 0)
        {
            sscanf(line, "%s\t%i\t%i\t%i", code, n, nht, nvt);
            printf("%i %i %i\n", *n, *nht, *nvt);
        }
        if (strncmp("INP4", line, 4) == 0)
        {
            sscanf(line, "%s\t%i\t%lf\t%i", code, ntimes, timestep_denom, freewake);
            printf("%i %f %i\n", *ntimes, *timestep_denom, *freewake);
        }
        if (strncmp("INP5", line, 4) == 0)
        {
            sscanf(line, "%s\t%lf\t%lf\t%lf\t%lf", code, &deltadegrees, &betadegrees, &etadegrees, &zetadegrees);
            printf("%f %f %f %f\n", deltadegrees, betadegrees, etadegrees, zetadegrees);
        }
        if (strncmp("INP6", line, 4) == 0)
        {
            sscanf(line, "%s\t%s", code, Wngfile);
            printf("%s\n", Wngfile);
        }
        if (strncmp("INP7", line, 4) == 0)
        {
            sscanf(line, "%s\t%s", code, HTailfile);
            printf("%s\n", HTailfile);
        }
        if (strncmp("INP8", line, 4) == 0)
        {
            sscanf(line, "%s\t%s", code, VTailfile);
            printf("%s\n", VTailfile);
        }
    }

    fclose(ifp);

    double pi=atan(1.0)*4.0;
    *aoa=aoadegrees*pi/180.0;
    *yaw=yawdegrees*pi/180.0;
    delta[0]=(-deltadegrees)*pi/180.0;    /* left aileron angle; negative means up*/
    delta[1]=deltadegrees*pi/180.0;     /* right aileron angle; negative means up */
    beta[0]=betadegrees*pi/180.0;      /* left flap angle; negative means up*/
    beta[1]=beta[0];            /* right flap angle; equal to left */
    eta[0]=etadegrees*pi/180.0;    /* left elevator angle; negative means up*/
    eta[1]=eta[0];              /* right elevator angle; equal to left */
    zeta[0]=zetadegrees*pi/180.0;      /* Real rudder */
    zeta[1]=zeta[0];            /* Rudder mirror image */
    /* Calculate free stream flow components */
    UVW[0]=Q*cos(*aoa)*cos(*yaw);
    UVW[1]=-Q*sin(*yaw);
    UVW[2]=Q*sin(*aoa)*cos(*yaw);
}