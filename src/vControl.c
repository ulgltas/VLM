//
//  vControl.c
//
//  Flight control stuff, such as control surface rotation
//
//

#include "vControl.h"
#include "vLiftsurf.h"
#include <math.h>

void rotateail(struct liftsurf *paileron,double delta[2])
{
    /* Rotate ailerons by angle delta */
    int i,j,k,nfaced2,nvertd2,te,nchordpanels,maxj;
    double x1[3],x2[3],uvw[3];
    double L,sindelta[2],cosdelta[2],sqrtL;
    double u,v,w,a,b,c,u2,v2,w2;
    double x,y,z,xnew,ynew,znew,dx,dy,dz;
    
    nfaced2=paileron->nface/2;
    nvertd2=paileron->nvert/2;
    nchordpanels=paileron->nface/paileron->nshed;
    
    sindelta[0]=sin(delta[0]);
    cosdelta[0]=cos(delta[0]);
    sindelta[1]=sin(-delta[1]);
    cosdelta[1]=cos(-delta[1]);         
    
    te=1; /* 1 if the previous panel was a trailing edge panel */
    for (k=0;k<2;k++){ /* k=0 for left aileron, k=1 for right aileron */        
        for (i=0;i<nfaced2;i++){
            if (te == 1 ){
                /* Find axis of rotation */
                /* The axis is given by the leading edge points of the first panel */
                for (j=0;j<3;j++){
                    x1[j]=*(paileron->vertices+j*paileron->nvert+*(paileron->faces+k*nfaced2+i));
                    x2[j]=*(paileron->vertices+j*paileron->nvert+*(paileron->faces+k*nfaced2+i+3*paileron->nface));
                    uvw[j]=x2[j]-x1[j]; /* Calculate direction vector of rotation axis */
                }
                u=uvw[0];
                v=uvw[1];
                w=uvw[2];
                u2=u*u;
                v2=v*v;
                w2=w*w;
                L=u2+v2+w2;
                sqrtL=sqrt(L);
                a=x1[0];
                b=x1[1];
                c=x1[2];
                /* Calculate how much the leading edge would have moved if it had been rotated */
                x=*(paileron->vertices+0*paileron->nvert+*(paileron->faces+k*nfaced2+i));
                y=*(paileron->vertices+1*paileron->nvert+*(paileron->faces+k*nfaced2+i));
                z=*(paileron->vertices+2*paileron->nvert+*(paileron->faces+k*nfaced2+i));
                xnew=(a*(v2+w2)-u*(b*v+c*w-u*x-v*y-w*z)*(1-cosdelta[k])+L*x*cosdelta[k]+sqrtL*(-c*v+b*w-w*y+v*z)*sindelta[k])/L;
                ynew=(b*(u2+w2)-v*(a*u+c*w-u*x-v*y-w*z)*(1-cosdelta[k])+L*y*cosdelta[k]+sqrtL*(c*u-a*w+w*x-u*z)*sindelta[k])/L;
                znew=(c*(u2+v2)-w*(a*u+b*v-u*x-v*y-w*z)*(1-cosdelta[k])+L*z*cosdelta[k]+sqrtL*(-b*u+a*v-v*x+u*y)*sindelta[k])/L;
                dx=xnew-a;
                dy=ynew-b;
                dz=znew-c;                
            }   
            /* For all panels, only rotate first trailing vertex */
            /* For tip panes also rotate second trailing vertex */
            maxj=2;
            if (i>=nfaced2-nchordpanels)
                maxj=3;
            for (j=1;j<maxj;j++){
                /* Rotate geometric panels */
                x=*(paileron->vertices+0*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface));
                y=*(paileron->vertices+1*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface));
                z=*(paileron->vertices+2*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface));
                /* Equations below from http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ */
                xnew=(a*(v2+w2)-u*(b*v+c*w-u*x-v*y-w*z)*(1-cosdelta[k])+L*x*cosdelta[k]+sqrtL*(-c*v+b*w-w*y+v*z)*sindelta[k])/L;
                ynew=(b*(u2+w2)-v*(a*u+c*w-u*x-v*y-w*z)*(1-cosdelta[k])+L*y*cosdelta[k]+sqrtL*(c*u-a*w+w*x-u*z)*sindelta[k])/L;
                znew=(c*(u2+v2)-w*(a*u+b*v-u*x-v*y-w*z)*(1-cosdelta[k])+L*z*cosdelta[k]+sqrtL*(-b*u+a*v-v*x+u*y)*sindelta[k])/L;
                /* Assign rotated values to paileron->vertices */
                *(paileron->vertices+0*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface))=xnew-dx;
                *(paileron->vertices+1*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface))=ynew-dy;
                *(paileron->vertices+2*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface))=znew-dz;
                /* Rotate vortex panels */
                x=*(paileron->vortex+0*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface));
                y=*(paileron->vortex+1*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface));
                z=*(paileron->vortex+2*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface));
                /* Equations below from http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ */
                xnew=(a*(v2+w2)-u*(b*v+c*w-u*x-v*y-w*z)*(1-cosdelta[k])+L*x*cosdelta[k]+sqrtL*(-c*v+b*w-w*y+v*z)*sindelta[k])/L;
                ynew=(b*(u2+w2)-v*(a*u+c*w-u*x-v*y-w*z)*(1-cosdelta[k])+L*y*cosdelta[k]+sqrtL*(c*u-a*w+w*x-u*z)*sindelta[k])/L;
                znew=(c*(u2+v2)-w*(a*u+b*v-u*x-v*y-w*z)*(1-cosdelta[k])+L*z*cosdelta[k]+sqrtL*(-b*u+a*v-v*x+u*y)*sindelta[k])/L;
                /* Assign rotated values to paileron->vertices */
                *(paileron->vortex+0*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface))=xnew-dx;
                *(paileron->vortex+1*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface))=ynew-dy;
                *(paileron->vortex+2*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface))=znew-dz;
            }
            /* Is this a trailing edge panel? */
            te=*(paileron->shedding+k*nfaced2+i); 
        }
    }
}