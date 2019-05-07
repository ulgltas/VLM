//
//  vInput.h
//
//  Data Input
//
//
#ifndef vInput_h
#define vInput_h
#include <stdio.h> 
void importInputFile(char *Infile, double *UVW, double *rho, double *aoa, double *yaw, int *m, int *mht, int *mvt, int *n, int *nht, int *nvt,
                    int *ntimes, double *timestep_denom, int *freewake, double *delta, double *beta, double *eta, double *zeta,
                    char *Wngfile, char *HTailfile, char *VTailfile);

#endif /* vInput_h */