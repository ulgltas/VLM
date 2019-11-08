//
//  vPanel.c
//  
//
//  Functions related to panel creation
//
//

#include "vPanel.h"
#include <math.h>
#include "vLiftsurf.h"

int countfaces(int *ijflap, int nflap, int mp1, int np1)
{
    /* Count the number of panels in lifting surface
     * ijflap contains the indices from the complete wing grid that
     * lie in the current lifting surface.
     * nflap is the number of vertices in the current lifting surface. */
    int nfaces, i, i1, j1, k, condij[3];
    
    nfaces=0; /* number of panels */
    for (i=0;i<nflap;i++){
        i1=*(ijflap+i); /* chordwise index */
        j1=*(ijflap+i+mp1*np1); /* Spanwise index */
        condij[0]=0;condij[1]=0;condij[2]=0;
        for (k=i+1;k<nflap;k++){
            /* Check if we can create a complete rectangle from point i1,j1 */
            if (*(ijflap+k)==i1+1 && *(ijflap+k+mp1*np1)==j1){
                condij[0]=k;
            }
            if (*(ijflap+k)==i1+1 && *(ijflap+k+mp1*np1)==j1+1){
                condij[1]=k;
            }
            if (*(ijflap+k)==i1 && *(ijflap+k+mp1*np1)==j1+1){
                condij[2]=k;
            }
        }
        /* If we can create the rectangle, increment nfaces */
        if (condij[0] != 0 && condij[1] != 0 && condij[2] != 0){
            nfaces++;
        }
    }
    return nfaces;
}

void arrangefaces(int *ijflap, int mp1, int np1, struct liftsurf *pflap)
{
    /* Create the connections between the vertices of the current lifting surface
     * that form the panels
     * pflap is the lifting surface structure
     * ijflap contains the indices from the complete wing grid that
     * lie in the current lifting surface.*/
    
    int nfaces, nfaced2, nvertd2, i, i1, j1, k, condij[3];
    
    nfaced2=pflap->nface/2;
    nvertd2=pflap->nvert/2;
    nfaces=0; /* number of panels */
    pflap->nshed=0; /* Initialize the number of shedding panels to 0 */
    for (i=0;i<nvertd2;i++){
        i1=*(ijflap+i); /* chordwise index */
        j1=*(ijflap+i+mp1*np1); /* Spanwise index */
        condij[0]=0;condij[1]=0;condij[2]=0;
        for (k=i+1;k<nvertd2;k++){
            /* Check if we can create a complete rectangle from point i1,j1
             * and store in condij the indices of each panel vertex */
            if (*(ijflap+k)==i1+1 && *(ijflap+k+mp1*np1)==j1){
                condij[0]=k;
            }
            if (*(ijflap+k)==i1+1 && *(ijflap+k+mp1*np1)==j1+1){
                condij[1]=k;
            }
            if (*(ijflap+k)==i1 && *(ijflap+k+mp1*np1)==j1+1){
                condij[2]=k;
            }
        }
        /* If we can create the rectangle, increment nfaces
         * and store the elements of condij in pflap->faces */
        if (condij[0] != 0 && condij[1] != 0 && condij[2] != 0){
            nfaces++;
            /* Right lifting surface */
            *(pflap->faces+nfaces-1)=i;
            *(pflap->faces+nfaces-1+pflap->nface)=condij[0];
            *(pflap->faces+nfaces-1+2*pflap->nface)=condij[1];
            *(pflap->faces+nfaces-1+3*pflap->nface)=condij[2];
            *(pflap->shedding+nfaces-1)=0;
            *(pflap->shedding+nfaces-1)=0;
            /* Left lifting surface */
            *(pflap->faces+nfaces-1+nfaced2)=i+nvertd2;
            *(pflap->faces+nfaces-1+nfaced2+pflap->nface)=condij[0]+nvertd2;
            *(pflap->faces+nfaces-1+nfaced2+2*pflap->nface)=condij[1]+nvertd2;
            *(pflap->faces+nfaces-1+nfaced2+3*pflap->nface)=condij[2]+nvertd2;
            *(pflap->shedding+nfaced2+nfaces-1)=0;
            if (i1 == mp1-2){ /* Check if this a trailing edge panel that sheds wake */
                /* Right lifting surface */
                *(pflap->shedding+nfaces-1)=1;
                pflap->nshed++;
                /* Left lifting surface */
                *(pflap->shedding+nfaces-1+nfaced2)=1;
                pflap->nshed++;
            }
        }
    }
}

void vortexpanel(double *xv, double *yv, double *zv, double *xp, double *yp, double *zp, double dxw, int m, int n)
{
    /* Calculate the collocation points and vector rings on a panel surface xp,yp,zp */
    int i,j,mp1,np1;
    double angp1;
    
    mp1=m+1;
    np1=n+1;
    /* Calculate vortex positions */
    for (j=0;j<np1;j++){
        for (i=0;i<m;i++){
            *(xv+i+j*mp1)=3.0* *(xp+i+j*mp1)/4.0+ *(xp+(i+1)+j*mp1)/4.0;
            *(yv+i+j*mp1)= *(yp+i+j*mp1);
            *(zv+i+j*mp1)=3.0* *(zp+i+j*mp1)/4.0+ *(zp+(i+1)+j*mp1)/4.0;
        }
        /* Trailing edge vortex panels */
        angp1=atan2(*(zp+m+j*mp1)- *(zp+(m-1)+j*mp1),*(xp+m+j*mp1)- *(xp+(m-1)+j*mp1));
        *(xv+m+j*mp1)= *(xp+m+j*mp1)+dxw* cos(angp1);
        *(yv+m+j*mp1)= *(yp+m+j*mp1);
        *(zv+m+j*mp1)= *(zp+m+j*mp1)+dxw* sin(angp1);
    }
}