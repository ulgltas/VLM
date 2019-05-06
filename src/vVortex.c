//
//  vVortex.c
//
//  Vortex related functions
//
//

#include "vVortex.h"
#include "vLiftsurf.h"
#include "vGeometry.h"
#include <math.h>
#include <stdlib.h>

void colvec(struct liftsurf *pflap)
{
    /* Calculate the collocation points and the vortex segment lengths on a lifting surface */
    int i,j;
    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    
    /* Calculate collocation points */  
    pflap->control=(double *)malloc(sizeof(double)*pflap->nface*3); 
    pflap->dxy=(double *)malloc(sizeof(double)*pflap->nface*2); 
    for (i=0;i<pflap->nface;i++){
        /* Vertices of each geometric panel */
        x1=*(pflap->vertices+*(pflap->faces+i));
        x2=*(pflap->vertices+*(pflap->faces+i+pflap->nface));
        x3=*(pflap->vertices+*(pflap->faces+i+2*pflap->nface));
        x4=*(pflap->vertices+*(pflap->faces+i+3*pflap->nface));
        y1=*(pflap->vertices+pflap->nvert+*(pflap->faces+i));
        y2=*(pflap->vertices+pflap->nvert+*(pflap->faces+i+pflap->nface));
        y3=*(pflap->vertices+pflap->nvert+*(pflap->faces+i+2*pflap->nface));
        y4=*(pflap->vertices+pflap->nvert+*(pflap->faces+i+3*pflap->nface));
        z1=*(pflap->vertices+2*pflap->nvert+*(pflap->faces+i));
        z2=*(pflap->vertices+2*pflap->nvert+*(pflap->faces+i+pflap->nface));
        z3=*(pflap->vertices+2*pflap->nvert+*(pflap->faces+i+2*pflap->nface));
        z4=*(pflap->vertices+2*pflap->nvert+*(pflap->faces+i+3*pflap->nface));
        /* Calculate collocation points */ 
        *(pflap->control+i)=(x1+x4)/8.0+3.0*(x2+x3)/8.0; /* x */
        *(pflap->control+i+pflap->nface)=(y1+y2+y3+y4)/4.0; /* y */
        *(pflap->control+i+2*pflap->nface)=(z1+z4)/8.0+3.0*(z2+z3)/8.0; /* z */
        /* Calculate vortex segment lengths */
        *(pflap->dxy+i+pflap->nface)=sqrt((x4-x1)*(x4-x1)+(y4-y1)*(y4-y1)+(z4-z1)*(z4-z1)); /* vortex segment length in y */
        *(pflap->dxy+i)=sqrt((x3-x4)*(x3-x4)+(y3-y4)*(y3-y4)+(z3-z4)*(z3-z4)); /* vortex segment length in x */
    }
}

void vortexblob(double uvw[], double x, double y, double z, double x1, double y1, double z1,
        double x2, double y2, double z2, double gama)
{
    /* Calculate the flow induced at point P(x,y,z) by a vortex blob segment with endpoints P1(x1,y1,z1) and P2(x2,y2,z2) */
    /* Inspired by the code in Katz and Plotkin */
    double rsquared,pi,rcut,r0[3],r1[3],r2[3],ir1,ir2,r1r2[3],square,coeff,blob,epsilon;
    
    pi=3.14159265358979;
    epsilon=0.1; /* Blob radius. This value seems to work very well */
    rcut=1.0e-11; /* Cut-off radius; if distance between point and segment smaller than this, flow is zero */
    /* Calculate vectors r0=P1-P2, r1=P-P1 and r2=P-P2 */
    r0[0]=x2-x1;
    r0[1]=y2-y1;
    r0[2]=z2-z1;
    r1[0]=x-x1;
    r1[1]=y-y1;
    r1[2]=z-z1;
    r2[0]=x-x2;
    r2[1]=y-y2;
    r2[2]=z-z2;
    /* Calculate lengths of r1 and r2 */
    ir1=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]);
    ir2=sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]);
    
    /* Calculate cross product between r1 and r2 */
    cross(r1r2,r1,r2);
    /* Calculate square of length of r1 x r2 */
    square=r1r2[0]*r1r2[0]+r1r2[1]*r1r2[1]+r1r2[2]*r1r2[2];
    /* Calculate blob factor */
    blob=(1.0-exp(-square/(epsilon*epsilon)));
    /* Check if P is too close to the vortex segment */
    if ((ir1 < rcut) || (ir2 < rcut) || (sqrt(square) < rcut)){
        uvw[0]=0.0;
        uvw[1]=0.0;
        uvw[2]=0.0;
    }else{
        /* Calculate flow induced by vortex segment on P */
        coeff=blob*gama/4/pi/square*(r0[0]*(r1[0]/ir1-r2[0]/ir2)+r0[1]*(r1[1]/ir1-r2[1]/ir2)+r0[2]*(r1[2]/ir1-r2[2]/ir2));
        uvw[0]=coeff*r1r2[0];
        uvw[1]=coeff*r1r2[1];
        uvw[2]=coeff*r1r2[2];
    }
}

void vortex(double uvw[], double x, double y, double z, double x1, double y1, double z1,
        double x2, double y2, double z2, double gama)
{
    /* Calculate the flow induced at point P(x,y,z) by a vortex segment with endpoints P1(x1,y1,z1) and P2(x2,y2,z2) */
    /* Inspired by the code in Katz and Plotkin */
    double rsquared,pi,rcut,r0[3],r1[3],r2[3],ir1,ir2,r1r2[3],square,coeff;
    
    pi=3.14159265358979;
    rcut=1.0e-11; /* Cut-off radius; if distance between point and segment smaller than this, flow is zero */
    /* Calculate vectors r0=P1-P2, r1=P-P1 and r2=P-P2 */
    r0[0]=x2-x1;
    r0[1]=y2-y1;
    r0[2]=z2-z1;
    r1[0]=x-x1;
    r1[1]=y-y1;
    r1[2]=z-z1;
    r2[0]=x-x2;
    r2[1]=y-y2;
    r2[2]=z-z2;
    /* Calculate lengths of r1 and r2 */
    ir1=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]);
    ir2=sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]);
    
    /* Calculate cross product between r1 and r2 */
    cross(r1r2,r1,r2);
    /* Calculate square of length of r1 x r2 */
    square=r1r2[0]*r1r2[0]+r1r2[1]*r1r2[1]+r1r2[2]*r1r2[2];
    /* Check if P is too close to the vortex segment */
    if ((ir1 < rcut) || (ir2 < rcut) || (square < rcut)){
        uvw[0]=0.0;
        uvw[1]=0.0;
        uvw[2]=0.0;
    }else{
        /* Calculate flow induced by vortex segment on P */
        coeff=gama/4/pi/square*(r0[0]*(r1[0]/ir1-r2[0]/ir2)+r0[1]*(r1[1]/ir1-r2[1]/ir2)+r0[2]*(r1[2]/ir1-r2[2]/ir2));
        uvw[0]=coeff*r1r2[0];
        uvw[1]=coeff*r1r2[1];
        uvw[2]=coeff*r1r2[2];
    }
}