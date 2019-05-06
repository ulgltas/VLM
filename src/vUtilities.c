//
//  vUtilities.c
//  
//
//  General things used throughout the code
//
//

#include "vUtilities.h"
#include <math.h>

int compare_function(const void *a,const void *b) {
    /* Comparison function for use with qsort */
    double *x = (double *) a;
    double *y = (double *) b;
    if (*x < *y) return -1;
    else if (*x > *y) return 1; return 0;
}

int checkdoubles(double *ypos, int m)
{
    /* Check if double entries exist in vectors */
    int n,i,j,nj;
    
    n=0;
    for (i=0;i<m;i++){
        for (j=i+1;j<m;j++){
            if (fabs(*(ypos+i)- *(ypos+j)) < 0.01){
                n++;
            }
        }
    }
    return n;
}

int finddoubles(double *ypos, int m)
{
    /* Find double entries in vectors */
    /* If there are double entries, they are pushed towards the end of the vector and overwritten */
    int n,i,j,nj;
    
    n=0;
    for (i=0;i<m;i++){
        nj=0;
        for (j=i+1;j<m;j++){
            if (fabs(*(ypos+i)- *(ypos+j)) < 0.01){
                nj=j;
                n++;
                for (j=nj+1;j<m;j++){
                    *(ypos+j-1)=*(ypos+j); /* Move all values from double to end of vector one element earlier */
                }
                *(ypos+m-1)=-1.0*n; /* Overwrite last element of vector with a non-double value */
                j--;
            }
        }
    }
    return m-n;
}

int findindex(double *ypline, int np1, double value)
{
    /* Find the index of the element of a vector that is closest to a number */
    int i,index;
    double dist,mindist;
    
    mindist=10000.0;
    for (i=0;i<np1;i++){
        dist=fabs(value-*(ypline+i));
        if (dist < mindist){
            mindist=dist;
            index=i;
        }
    }
    return index;
}