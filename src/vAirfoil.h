//
//  vAirfoil.h
//  
//
//  Airfoil processing, reading and interpolating functions
//
//

#ifndef vAirfoil_h
#define vAirfoil_h

#include <stdio.h>

void readarfsize(char *WngArf,int ArfXUNumber[],int ArfXLNumber[]);

void readarf(char *WngArf,double *ArfXU, double *ArfYU, double *ArfXL,double *ArfYL,int ArfXUNumber[],int ArfXLNumber[]);

void interparf(double *ycamber, double *ArfXU, double *ArfYU, double *ArfXL, double *ArfYL, int ArfXUNumber[], int ArfXLNumber[], double *xpline, int mp1);

void nacafourfivedigit(double *xpline, double *ycamber, int mp1, char WngArf[]);

void treatarfwgl(char line[], double *xpline, double *ycamber, int mp1);

void treatarf(char line[], double *xpline, double *ycamber, double *ycamberall, int mp1, int iTS);

#endif /* vAirfoil_h */
