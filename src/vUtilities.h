//
//  vUtilities.h
//  
//
//  General things used throughout the code
//
//

#ifndef vUtilities_h
#define vUtilities_h

#include <stdio.h>

int compare_function(const void *a,const void *b);

int checkdoubles(double *ypos, int m);

int finddoubles(double *ypos, int m);

int findindex(double *ypline, int np1, double value);

#endif /* vUtilities_h */
