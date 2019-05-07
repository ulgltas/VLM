#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vInput.h"
#include "vOutput.h"
#include "vForce.h"
#include "vGeometry.h"
#include "vLiftsurf.h"
#include "vVortex.h"
#include "vWake.h"
#include "vHTail.h"
#include "vVTail.h"
#include "vWing.h"


void adaptycoord(double *ypline, double *ypos, int np1, int nypos)
{
    /* Adapt values of vector ypline so that they include the values of vector ypos */
    int i,j,nmin;
    double mindist,dist;
    
    nmin=-1;
    for (i=0;i<nypos;i++){
        mindist=1000.0;
        for (j=nmin+1;j<np1;j++){
            dist=fabs(*(ypos+i)- *(ypline+j)); /* Look for the value of ypline nearest to this element of ypos */
            if (dist <= mindist){
                mindist=dist;
                nmin=j;
            }
        }
        *(ypline+nmin)=*(ypos+i); /* Change the nearest value of ypos */
    }
}

void findneighbours(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, int cwing, int cflap, int caileron)
{
    /* Find all the neighbours of all the wing panels in all wing, flap and aileron panels */
    int i,j,k,l,repeat[4];
    double vertx[4],verty[4],vertz[4],dist[3],totaldist;
    
    pwing->neighbours=(int *)malloc(sizeof(int)*pwing->nface*4);
    for (i=0;i<pwing->nface;i++){
        for (k=0;k<4;k++){
            vertx[k]=*(pwing->vertices+ *(pwing->faces+i+k*pwing->nface));
            verty[k]=*(pwing->vertices+pwing->nvert+ *(pwing->faces+i+k*pwing->nface));
            vertz[k]=*(pwing->vertices+2*pwing->nvert+ *(pwing->faces+i+k*pwing->nface));
            *(pwing->neighbours+i+k*pwing->nface)=-1;
        }
        /* Check if wing panels have neighbours on the wing */
        for (j=0;j<pwing->nface;j++){
            if (i != j){
                for (k=0;k<4;k++){
                    repeat[k]=0;
                    for (l=0;l<4;l++){
                        dist[0]=vertx[k]-*(pwing->vertices+ *(pwing->faces+j+l*pwing->nface));
                        dist[1]=verty[k]-*(pwing->vertices+pwing->nvert+ *(pwing->faces+j+l*pwing->nface));
                        dist[2]=vertz[k]-*(pwing->vertices+2*pwing->nvert+ *(pwing->faces+j+l*pwing->nface));
                        totaldist=sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
                        if (totaldist < 0.001)
                            repeat[k]=j+1; /* +1 so that j=0 does not results in repeat[k]=0 */
                    }
                }
                if (repeat[0] != 0 && repeat[1] !=0){
                    if (i<pwing->nface/2){
                        *(pwing->neighbours+i+pwing->nface)=cwing+j; /* right neighbour looking from leading edge to trailing edge*/
                    }else
                        *(pwing->neighbours+i)=cwing+j; /* left neighbour looking from leading edge to trailing edge*/
                }
                if (repeat[1] != 0 && repeat[2] !=0){
                    *(pwing->neighbours+i+3*pwing->nface)=cwing+j; /* downstream neighbour looking from leading edge to trailing edge*/
                }
                if (repeat[2] != 0 && repeat[3] !=0){
                    if (i<pwing->nface/2){
                        *(pwing->neighbours+i)=cwing+j; /* left neighbour looking from leading edge to trailing edge*/
                    }else
                        *(pwing->neighbours+i+pwing->nface)=cwing+j; /* right neighbour looking from leading edge to trailing edge*/
                }
                if (repeat[3] != 0 && repeat[0] !=0){
                    *(pwing->neighbours+i+2*pwing->nface)=cwing+j; /* upstream neighbour looking from leading edge to trailing edge*/
                }
            }
        }
        /* Now check if wing panels have neighbours on the flap */
        for (j=0;j<pflap->nface;j++){
            for (k=0;k<4;k++){
                repeat[k]=0;
                for (l=0;l<4;l++){
                    dist[0]=vertx[k]-*(pflap->vertices+ *(pflap->faces+j+l*pflap->nface));
                    dist[1]=verty[k]-*(pflap->vertices+pflap->nvert+ *(pflap->faces+j+l*pflap->nface));
                    dist[2]=vertz[k]-*(pflap->vertices+2*pflap->nvert+ *(pflap->faces+j+l*pflap->nface));
                    totaldist=sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
                    if (totaldist < 0.001)
                        repeat[k]=j+1; /* +1 so that j=0 does not result in repeat[k]=0 */
                }
            }
            if (repeat[0] != 0 && repeat[1] !=0){
                if (i<pwing->nface/2){
                    *(pwing->neighbours+i+pwing->nface)=cflap+j; /* right neighbour looking from leading edge to trailing edge*/
                }else
                    *(pwing->neighbours+i)=cflap+j; /* left neighbour looking from leading edge to trailing edge*/
            }
            if (repeat[1] != 0 && repeat[2] !=0){
                *(pwing->neighbours+i+3*pwing->nface)=cflap+j; /* downstream neighbour looking from leading edge to trailing edge*/
            }
            if (repeat[2] != 0 && repeat[3] !=0){
                if (i<pwing->nface/2){
                    *(pwing->neighbours+i)=cflap+j; /* left neighbour looking from leading edge to trailing edge*/
                }else
                    *(pwing->neighbours+i+pwing->nface)=cflap+j; /* right neighbour looking from leading edge to trailing edge*/
            }
            if (repeat[3] != 0 && repeat[0] !=0){
                *(pwing->neighbours+i+2*pwing->nface)=cflap+j; /* upstream neighbour looking from leading edge to trailing edge*/
            }
        }
        /* Now check if wing panels have neighbours on the aileron */
        for (j=0;j<paileron->nface;j++){
            for (k=0;k<4;k++){
                repeat[k]=0;
                for (l=0;l<4;l++){
                    dist[0]=vertx[k]-*(paileron->vertices+ *(paileron->faces+j+l*paileron->nface));
                    dist[1]=verty[k]-*(paileron->vertices+paileron->nvert+ *(paileron->faces+j+l*paileron->nface));
                    dist[2]=vertz[k]-*(paileron->vertices+2*paileron->nvert+ *(paileron->faces+j+l*paileron->nface));
                    totaldist=sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
                    if (totaldist < 0.001)
                        repeat[k]=j+1; /* +1 so that j=0 does not results in repeat[k]=0 */
                }
            }
            if (repeat[0] != 0 && repeat[1] !=0){
                if (i<pwing->nface/2){
                    *(pwing->neighbours+i+pwing->nface)=caileron+j; /* right neighbour looking from leading edge to trailing edge*/
                }else
                    *(pwing->neighbours+i)=caileron+j; /* left neighbour looking from leading edge to trailing edge*/
            }
            if (repeat[1] != 0 && repeat[2] !=0){
                *(pwing->neighbours+i+3*pwing->nface)=caileron+j; /* downstream neighbour looking from leading edge to trailing edge*/
            }
            if (repeat[2] != 0 && repeat[3] !=0){
                if (i<pwing->nface/2){
                    *(pwing->neighbours+i)=caileron+j; /* left neighbour looking from leading edge to trailing edge*/
                }else
                    *(pwing->neighbours+i+pwing->nface)=caileron+j; /* right neighbour looking from leading edge to trailing edge*/
            }
            if (repeat[3] != 0 && repeat[0] !=0){
                *(pwing->neighbours+i+2*pwing->nface)=caileron+j; /* upstream neighbour looking from leading edge to trailing edge*/
            }
        }
    }
}

void rotateail(struct liftsurf *paileron,double delta[])
{
    /* Rotate ailerons by angle delta */
    int i,j,k,nfaced2,nvertd2,te,nchordpanels,maxj;
    double x1[3],x2[3],uvw[3],T[16];
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

void infcoeff(struct liftsurf *plift1, struct liftsurf *plift2, double *AN, double *BN, int istart, int jstart,int mtn)
{
    /* Calculate the matrices of flow influence coefficients */
    /* This is the influence of plift2 on plift1 */
    /* AN is the matrix of normal flow influence coefficients */
    /* BN is the matrix of downwash influence coefficients */
    int i,j;
    double xc,yc,zc;
    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    double uvw1[3],uvw2[3],uvw3[3],uvw4[3],uvw[3],uvwstar[3];
    
    for (i=0; i < plift1->nface; i++) {
        /* Control point */
        xc=*(plift1->control+i);
        yc=*(plift1->control+i+plift1->nface);
        zc=*(plift1->control+i+2*plift1->nface);
        /* Vortex rings are quadrilaterals with four vortex segments each */
        for (j=0; j < plift2->nface; j++) {
            /* Vortex panel vertices */
            x1=*(plift2->vortex+*(plift2->faces+j));
            x2=*(plift2->vortex+*(plift2->faces+j+plift2->nface));
            x3=*(plift2->vortex+*(plift2->faces+j+2*plift2->nface));
            x4=*(plift2->vortex+*(plift2->faces+j+3*plift2->nface));
            y1=*(plift2->vortex+plift2->nvert+*(plift2->faces+j));
            y2=*(plift2->vortex+plift2->nvert+*(plift2->faces+j+plift2->nface));
            y3=*(plift2->vortex+plift2->nvert+*(plift2->faces+j+2*plift2->nface));
            y4=*(plift2->vortex+plift2->nvert+*(plift2->faces+j+3*plift2->nface));
            z1=*(plift2->vortex+2*plift2->nvert+*(plift2->faces+j));
            z2=*(plift2->vortex+2*plift2->nvert+*(plift2->faces+j+plift2->nface));
            z3=*(plift2->vortex+2*plift2->nvert+*(plift2->faces+j+2*plift2->nface));
            z4=*(plift2->vortex+2*plift2->nvert+*(plift2->faces+j+3*plift2->nface));
            if (j<plift2->nface/2){
                /* Influence of vortex segment 1 */
                vortex(uvw1,xc,yc,zc,x1,y1,z1,x4,y4,z4,1.0);
                /* Influence of vortex segment 2 */
                vortex(uvw2,xc,yc,zc,x4,y4,z4,x3,y3,z3,1.0);
                /* Influence of vortex segment 3 */
                vortex(uvw3,xc,yc,zc,x3,y3,z3,x2,y2,z2,1.0);
                /* Influence of vortex segment 4 */
                vortex(uvw4,xc,yc,zc,x2,y2,z2,x1,y1,z1,1.0);
            }else{
                /* Influence of vortex segment 1 */
                vortex(uvw1,xc,yc,zc,x4,y4,z4,x1,y1,z1,1.0);
                /* Influence of vortex segment 2 */
                vortex(uvw2,xc,yc,zc,x1,y1,z1,x2,y2,z2,1.0);
                /* Influence of vortex segment 3 */
                vortex(uvw3,xc,yc,zc,x2,y2,z2,x3,y3,z3,1.0);
                /* Influence of vortex segment 4 */
                vortex(uvw4,xc,yc,zc,x3,y3,z3,x4,y4,z4,1.0);
            }
            /* Addd the four influcences to the normal flow */
            uvw[0]=uvw1[0]+uvw2[0]+uvw3[0]+uvw4[0];
            uvw[1]=uvw1[1]+uvw2[1]+uvw3[1]+uvw4[1];
            uvw[2]=uvw1[2]+uvw2[2]+uvw3[2]+uvw4[2];
            /* Addd the two influcences to the downwash */
            uvwstar[0]=uvw2[0]+uvw4[0];
            uvwstar[1]=uvw2[1]+uvw4[1];
            uvwstar[2]=uvw2[2]+uvw4[2];
            /* Create the AN and BN matrices */
            *(AN+istart+i+(jstart+j)*mtn)=uvw[0]* *(plift1->normal+i)+uvw[1]* *(plift1->normal+i+plift1->nface) +
                    uvw[2]* *(plift1->normal+i+2*plift1->nface);
            *(BN+istart+i+(jstart+j)*mtn)=uvwstar[0]* *(plift1->normal+i)+uvwstar[1]* *(plift1->normal+i+plift1->nface) +
                    uvwstar[2]* *(plift1->normal+i+2*plift1->nface);
        }
    }
}

void cycleliftsurf(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder, double *AN, double *BN, int mtn)
{
    int m,n;
    
    /* Cycle through the lifting surfaces in the order wing, flap, aileron, htail, elevator, vtail, rudder */
    /* Treat the vtal and rudder problem as uncoupled */
    m=0;
    n=0;
    /* Influence on the wing from the wing */
    infcoeff(pwing,pwing,AN,BN,m,n,mtn);
    /* Influence on the wing from the flap */
    n+=pwing->nface;
    infcoeff(pwing,pflap,AN,BN,m,n,mtn);
    /* Influence on the wing from the aileron */
    n+=pflap->nface;
    infcoeff(pwing,paileron,AN,BN,m,n,mtn);
    /* Influence on the wing from the horizontal tail */
    n+=paileron->nface;
    infcoeff(pwing,phtail,AN,BN,m,n,mtn);
    /* Influence on the wing from the elevator */
    n+=phtail->nface;
    infcoeff(pwing,pelevator,AN,BN,m,n,mtn);
//     /* Influence on the wing from the vertical tail */
//     n+=pelevator->nface;
//     infcoeff(pwing,pvtail,AN,BN,m,n,mtn);
//     /* Influence on the wing from the rudder */
//     n+=pvtail->nface;
//     infcoeff(pwing,prudder,AN,BN,m,n,mtn);

    m+=pwing->nface;
    n=0;
    /* Influence on the flap from the wing */
    infcoeff(pflap,pwing,AN,BN,m,n,mtn);
    /* Influence on the flap from the flap */
    n+=pwing->nface;
    infcoeff(pflap,pflap,AN,BN,m,n,mtn);
    /* Influence on the flap from the aileron */
    n+=pflap->nface;
    infcoeff(pflap,paileron,AN,BN,m,n,mtn);
    /* Influence on the flap from the horizontal tail */
    n+=paileron->nface;
    infcoeff(pflap,phtail,AN,BN,m,n,mtn);
    /* Influence on the flap from the elevator */
    n+=phtail->nface;
    infcoeff(pflap,pelevator,AN,BN,m,n,mtn);
//     /* Influence on the flap from the vertical tail */
//     n+=pelevator->nface;
//     infcoeff(pflap,pvtail,AN,BN,m,n,mtn);
//     /* Influence on the flap from the rudder */
//     n+=pvtail->nface;
//     infcoeff(pflap,prudder,AN,BN,m,n,mtn);
    
    m+=pflap->nface;
    n=0;
    /* Influence on the aileron from the wing */
    infcoeff(paileron,pwing,AN,BN,m,n,mtn);
    /* Influence on the aileron from the flap */
    n+=pwing->nface;
    infcoeff(paileron,pflap,AN,BN,m,n,mtn);
    /* Influence on the aileron from the aileron */
    n+=pflap->nface;
    infcoeff(paileron,paileron,AN,BN,m,n,mtn);
    /* Influence on the aileron from the horizontal tail */
    n+=paileron->nface;
    infcoeff(paileron,phtail,AN,BN,m,n,mtn);
    /* Influence on the aileron from the elevator */
    n+=phtail->nface;
    infcoeff(paileron,pelevator,AN,BN,m,n,mtn);
//     /* Influence on the aileron from the vertical tail */
//     n+=pelevator->nface;
//     infcoeff(paileron,pvtail,AN,BN,m,n,mtn);
//     /* Influence on the aileron from the rudder */
//     n+=pvtail->nface;
//     infcoeff(paileron,prudder,AN,BN,m,n,mtn);
    
    m+=paileron->nface;
    n=0;
    /* Influence on the horizontal tail from the wing */
    infcoeff(phtail,pwing,AN,BN,m,n,mtn);
    /* Influence on the horizontal tail from the flap */
    n+=pwing->nface;
    infcoeff(phtail,pflap,AN,BN,m,n,mtn);
    /* Influence on the horizontal tail from the aileron */
    n+=pflap->nface;
    infcoeff(phtail,paileron,AN,BN,m,n,mtn);
    /* Influence on the horizontal tail from the horizontal tail */
    n+=paileron->nface;
    infcoeff(phtail,phtail,AN,BN,m,n,mtn);
    /* Influence on the horizontal tail from the elevator */
    n+=phtail->nface;
    infcoeff(phtail,pelevator,AN,BN,m,n,mtn);
//     /* Influence on the horizontal tail from the vertical tail */
//     n+=pelevator->nface;
//     infcoeff(phtail,pvtail,AN,BN,m,n,mtn);
//     /* Influence on the horizontal tail from the rudder */
//     n+=pvtail->nface;
//     infcoeff(phtail,prudder,AN,BN,m,n,mtn);
    
    m+=phtail->nface;
    n=0;
    /* Influence on the elevator from the wing */
    infcoeff(pelevator,pwing,AN,BN,m,n,mtn);
    /* Influence on the elevator from the flap */
    n+=pwing->nface;
    infcoeff(pelevator,pflap,AN,BN,m,n,mtn);
    /* Influence on the elevator from the aileron */
    n+=pflap->nface;
    infcoeff(pelevator,paileron,AN,BN,m,n,mtn);
    /* Influence on the elevator from the horizontal tail */
    n+=paileron->nface;
    infcoeff(pelevator,phtail,AN,BN,m,n,mtn);
    /* Influence on the elevator from the elevator */
    n+=phtail->nface;
    infcoeff(pelevator,pelevator,AN,BN,m,n,mtn);
//     /* Influence on the elevator from the vertical tail */
//     n+=pelevator->nface;
//     infcoeff(pelevator,pvtail,AN,BN,m,n,mtn);
//     /* Influence on the elevator from the rudder */
//     n+=pvtail->nface;
//     infcoeff(pelevator,prudder,AN,BN,m,n,mtn);

    m+=pelevator->nface;
    n=0;
    /* Influence on the vertical tail from the wing */
//    infcoeff(pvtail,pwing,AN,BN,m,n,mtn);
    /* Influence on the vertical tail from the flap */
    n+=pwing->nface;
//    infcoeff(pvtail,pflap,AN,BN,m,n,mtn);
    /* Influence on the vertical tail from the aileron */
    n+=pflap->nface;
//    infcoeff(pvtail,paileron,AN,BN,m,n,mtn);
    /* Influence on the vertical tail from the horizontal tail */
    n+=paileron->nface;
//    infcoeff(pvtail,phtail,AN,BN,m,n,mtn);
    /* Influence on the vertical tail from the elevator */
    n+=phtail->nface;
//    infcoeff(pvtail,pelevator,AN,BN,m,n,mtn);
    /* Influence on the vertical tail from the vertical tail */
    n+=pelevator->nface;
    infcoeff(pvtail,pvtail,AN,BN,m,n,mtn);
    /* Influence on the vertical tail from the rudder */
    n+=pvtail->nface;
    infcoeff(pvtail,prudder,AN,BN,m,n,mtn);
    
    m+=pvtail->nface;
    n=0;
    /* Influence on the rudder from the wing */
//    infcoeff(prudder,pwing,AN,BN,m,n,mtn);
    /* Influence on the rudder from the flap */
    n+=pwing->nface;
//    infcoeff(prudder,pflap,AN,BN,m,n,mtn);
    /* Influence on the rudder from the aileron */
    n+=pflap->nface;
//    infcoeff(prudder,paileron,AN,BN,m,n,mtn);
    /* Influence on the rudder from the horizontal tail */
    n+=paileron->nface;
//    infcoeff(prudder,phtail,AN,BN,m,n,mtn);
    /* Influence on the rudder from the elevator */
    n+=phtail->nface;
//    infcoeff(prudder,pelevator,AN,BN,m,n,mtn);
    /* Influence on the rudder from the vertical tail */
    n+=pelevator->nface;
    infcoeff(prudder,pvtail,AN,BN,m,n,mtn);
    /* Influence on the rudder from the rudder */
    n+=pvtail->nface;
    infcoeff(prudder,prudder,AN,BN,m,n,mtn);
}

void gauss(double *a, double *inva, int nrows)
{
    /* Carry out Gaussian elimination on matrix a to calculate its inverse, inva */
    /* Inspired by the algorithm in Numerical Recipes in C */
    int i,j,k;
    double *MI, factor,dummy;
    
    MI = (double *)malloc(sizeof(double)*(2*nrows)*(nrows));
    
    /* Create augmented matrix */
    for (j=0; j<nrows; j++) {
        for (i=0;i<nrows; i++){
            *(MI+i+j*nrows)=*(a+j*nrows+i);
            if (i==j){
                *(MI+i+(j+nrows)*nrows)=1.0;
            } else {
                *(MI+i+(j+nrows)*nrows)=0.0;
            }
        }
    }
    
    /* Create inverse matrix */
    for (i=0; i<nrows; i++) {
        dummy=*(MI+i+i*nrows);
        for (j=0; j<2*nrows; j++) {
            *(MI+i+j*nrows)=*(MI+i+j*nrows)/dummy;
        }
        for (k=0; k<nrows; k++) {
            if (k!=i){
                factor=- *(MI+k+i*nrows)/ *(MI+i+i*nrows);
                for (j=0; j<2*nrows; j++) {
                    *(MI+k+j*nrows) = *(MI+k+j*nrows) + factor* *(MI+i+j*nrows);
                }
            }
        }
    }
    
    /* Extract inverse matrix */
    for (j=0; j<nrows; j++) {
        for(i=0; i<nrows; i++) {
            *(inva+i+nrows*j)=*(MI+i+(j+nrows)*nrows);
        }
    }
    free(MI);
}

void calcRHS(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder, double *RHS, double UVW[])
{   
    /* Calculate the Right Hand Side vector */
    /* Cycle through the lifting surfaces in the order wing, flap, aileron, horizontal tail, elevator */
    
    int i,m;

    m=0;
    for (i=0; i < pwing->nface; i++){
        *(RHS+i)=-(*(pwing->uvw+i)+UVW[0])* *(pwing->normal+i)-(*(pwing->uvw+i+pwing->nface)+UVW[1])* *(pwing->normal+i+pwing->nface) -
                (*(pwing->uvw+i+2*pwing->nface)+UVW[2])* *(pwing->normal+i+2*pwing->nface);
    }
    m+=pwing->nface;
    for (i=0; i < pflap->nface; i++){
        *(RHS+i+m)=-(*(pflap->uvw+i)+UVW[0])* *(pflap->normal+i)-(*(pflap->uvw+i+pflap->nface)+UVW[1])* *(pflap->normal+i+pflap->nface) -
                (*(pflap->uvw+i+2*pflap->nface)+UVW[2])* *(pflap->normal+i+2*pflap->nface);
    }
    m+=pflap->nface;
    for (i=0; i < paileron->nface; i++){
        *(RHS+i+m)=-(*(paileron->uvw+i)+UVW[0])* *(paileron->normal+i)-(*(paileron->uvw+i+paileron->nface)+UVW[1])* *(paileron->normal+i+paileron->nface) -
                (*(paileron->uvw+i+2*paileron->nface)+UVW[2])* *(paileron->normal+i+2*paileron->nface);
    }
    m+=paileron->nface;
    for (i=0; i < phtail->nface; i++){
        *(RHS+i+m)=-(*(phtail->uvw+i)+UVW[0])* *(phtail->normal+i)-(*(phtail->uvw+i+phtail->nface)+UVW[1])* *(phtail->normal+i+phtail->nface) -
                (*(phtail->uvw+i+2*phtail->nface)+UVW[2])* *(phtail->normal+i+2*phtail->nface);
    }
    m+=phtail->nface;
    for (i=0; i < pelevator->nface; i++){
        *(RHS+i+m)=-(*(pelevator->uvw+i)+UVW[0])* *(pelevator->normal+i)-(*(pelevator->uvw+i+pelevator->nface)+UVW[1])* *(pelevator->normal+i+pelevator->nface) -
                (*(pelevator->uvw+i+2*pelevator->nface)+UVW[2])* *(pelevator->normal+i+2*pelevator->nface);
    }
    m+=pelevator->nface;
    for (i=0; i < pvtail->nface; i++){
        *(RHS+i+m)=-(*(pvtail->uvw+i)+UVW[0])* *(pvtail->normal+i)-(*(pvtail->uvw+i+pvtail->nface)+UVW[1])* *(pvtail->normal+i+pvtail->nface) -
                (*(pvtail->uvw+i+2*pvtail->nface)+UVW[2])* *(pvtail->normal+i+2*pvtail->nface);
    }
    m+=pvtail->nface;
    for (i=0; i < prudder->nface; i++){
        *(RHS+i+m)=-(*(prudder->uvw+i)+UVW[0])* *(prudder->normal+i)-(*(prudder->uvw+i+prudder->nface)+UVW[1])* *(prudder->normal+i+prudder->nface) -
                (*(prudder->uvw+i+2*prudder->nface)+UVW[2])* *(prudder->normal+i+2*prudder->nface);
    }
}

void assignGammas(double *Gammas, struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder, int it)
{
    /* Assign the vortex strengths in vector Gammas back onto the panels of all the lifting surfaces */
    int i,m;
        
    /* Assign wing vortex strengths */
    m=0;
    for (i=0;i<pwing->nface;i++){
        *(pwing->gamma+i+it*pwing->nface)=*(Gammas+i);
    }
    /* Assign flap vortex strengths */
    m+=pwing->nface;
    for (i=0;i<pflap->nface;i++){
        *(pflap->gamma+i+it*pflap->nface)=*(Gammas+i+m);
    }
    /* Assign aileron vortex strengths */
    m+=pflap->nface;
    for (i=0;i<paileron->nface;i++){
        *(paileron->gamma+i+it*paileron->nface)=*(Gammas+i+m);
    }
    m+=paileron->nface;
    /* Assign horizontal tail vortex strengths */
    for (i=0;i<phtail->nface;i++){
        *(phtail->gamma+i+it*phtail->nface)=*(Gammas+i+m);
    }  
    m+=phtail->nface;
    /* Assign elevator vortex strengths */
    for (i=0;i<pelevator->nface;i++){
        *(pelevator->gamma+i+it*pelevator->nface)=*(Gammas+i+m);
    }
    m+=pelevator->nface;
    /* Assign vertical tail vortex strengths */
    for (i=0;i<pvtail->nface;i++){
        *(pvtail->gamma+i+it*pvtail->nface)=*(Gammas+i+m);
    }
    m+=pvtail->nface;
    /* Assign rudder vortex strengths */
    for (i=0;i<prudder->nface;i++){
        *(prudder->gamma+i+it*prudder->nface)=*(Gammas+i+m);
    }    
}

void assignwind(double *wind, struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder)
{
    /* Assign the vortex strengths in vector Gammas back onto the panels of all the lifting surfaces */
    int i,m;
        
    /* Assign wing vortex strengths */
    m=0;
    for (i=0;i<pwing->nface;i++){
        *(pwing->wind+i)=*(wind+i);
    }
    /* Assign flap vortex strengths */
    m+=pwing->nface;
    for (i=0;i<pflap->nface;i++){
        *(pflap->wind+i)=*(wind+i+m);
    }
    /* Assign aileron vortex strengths */
    m+=pflap->nface;
    for (i=0;i<paileron->nface;i++){
        *(paileron->wind+i)=*(wind+i+m);
    }
    m+=paileron->nface;
    /* Assign horizontal tail vortex strengths */
    for (i=0;i<phtail->nface;i++){
        *(phtail->wind+i)=*(wind+i+m);
    }  
    m+=phtail->nface;
    /* Assign elevator vortex strengths */
    for (i=0;i<pelevator->nface;i++){
        *(pelevator->wind+i)=*(wind+i+m);
    }
    m+=pelevator->nface;
    /* Assign vertical tail vortex strengths */
    for (i=0;i<pvtail->nface;i++){
        *(pvtail->wind+i)=*(wind+i+m);
    }  
    m+=pvtail->nface;
    /* Assign rudder vortex strengths */
    for (i=0;i<prudder->nface;i++){
        *(prudder->wind+i)=*(wind+i+m);
    }      
}

int main(int argc, char *argv[])
{
    char *Infile, Wngfile[60], HTailfile[60], VTailfile[60], *Outfile;
    struct liftsurf flap,aileron,wing,htail,elevator,vtail,rudder;
    struct liftsurf *pflap = &flap;
    struct liftsurf *paileron = &aileron;
    struct liftsurf *pwing = &wing;
    struct liftsurf *phtail = &htail;
    struct liftsurf *pelevator = &elevator;
    struct liftsurf *pvtail = &vtail;
    struct liftsurf *prudder = &rudder;
    double *invAN,*AN,*BN,*RHS, *Gammas, *wind;
    double UVW[3],aoa,yaw,dt,timestep_denom,MAC,rho,*totalforce;
    double delta[2],beta[2],eta[2],zeta[2];
    int i,j,mtn,it,ntimes,m,n,freewake,mht,nht,mvt,nvt;
    
    if (argc>1) {
        Infile=argv[1];
    } else {
        Infile="infile.arp";
    }
    
    Outfile="outfile.m";
    
    importInputFile(Infile, UVW, &rho, &aoa, &yaw, &m, &mht, &mvt, &n, &nht, &nvt, &ntimes, &timestep_denom, &freewake,
                    delta, beta, eta, zeta, Wngfile, HTailfile, VTailfile);
        
    totalforce=(double *)calloc(ntimes*4, sizeof(double));
    
    /* Create the lifting surfaces that make up the wing */
    MAC=wingsetup(pflap,paileron,pwing,Wngfile,m,n);
    /* Create the lifting surfaces that make up the horizontal tail */
    htailsetup(phtail,pelevator,HTailfile,mht,nht);
    /* Create the lifting surfaces that make up the vertical tail */
    vtailsetup(pvtail,prudder,VTailfile,mht,nht);
    dt=MAC/UVW[0]/timestep_denom; /* dt is based on the length of the wing's Mean Aerodynamic Chord */
    printf("MAC=%f, dt=%f\n",MAC,dt);
    
    /* Rotate ailerons */
    rotateail(paileron,delta);
    /* Rotate flaps */
    rotateail(pflap,beta);
    /* Rotate elevator */
    rotateail(pelevator,eta);
    /* Rotate rudder */
    rotateail(prudder,zeta);
    
    findneighbours(pwing,pflap,paileron,0,10000,20000);
    findneighbours(pflap,pwing,paileron,10000,0,20000);
    findneighbours(paileron,pwing,pflap,20000,0,10000);
    findneighbours(phtail,pelevator,paileron,30000,40000,20000); /* We don't care about paileron being called here */
    findneighbours(pelevator,phtail,paileron,40000,30000,20000); /* We don't care about paileron being called here */
    findneighbours(pvtail,prudder,paileron,50000,60000,20000); /* We don't care about paileron being called here */
    findneighbours(prudder,pvtail,paileron,60000,50000,20000); /* We don't care about paileron being called here */

    /* Create the wake */
    createwake(pwing,0,ntimes);
    createwake(pflap,10000,ntimes);
    createwake(paileron,20000,ntimes); 
    createwake(phtail,30000,ntimes); 
    createwake(pelevator,40000,ntimes); 
    createwake(pvtail,50000,ntimes); 
    createwake(prudder,60000,ntimes); 
    
    /* Find which bound vortex panels correspond to which wake panel vertices */
    correspshedwake(pwing);
    correspshedwake(pflap);
    correspshedwake(paileron);
    correspshedwake(phtail);
    correspshedwake(pelevator);
    correspshedwake(pvtail);
    correspshedwake(prudder);
    
    /* Shed first wake element */
    shedwake(pwing);
    shedwake(pflap);
    shedwake(paileron);
    shedwake(phtail);
    shedwake(pelevator);
    shedwake(pvtail);
    shedwake(prudder);
    findallwakeneighbours(pwing,pflap,paileron,0,10000,20000);
    findallwakeneighbours(phtail,pelevator,paileron,30000,40000,20000); /* We don't care about paileron being called here */
    findallwakeneighbours(pvtail,prudder,paileron,50000,60000,20000); /* We don't care about paileron being called here */
    
    /* Calculate collocation points and vortex segment lengths */
    colvec(pflap);
    colvec(paileron);
    colvec(pwing);
    colvec(phtail);
    colvec(pelevator);
    colvec(pvtail);
    colvec(prudder);
    /* Calculate normal vectors and surfaces */
    normals(pflap);
    normals(paileron);
    normals(pwing);
    normals(phtail);
    normals(pelevator);
    normals(pvtail);
    normals(prudder);
    /* Calculate tangential vectors */
    tangentials(pflap);
    tangentials(paileron);
    tangentials(pwing);
    tangentials(phtail);
    tangentials(pelevator);
    tangentials(pvtail);
    tangentials(prudder);
    /* Calculate total number of panels */
    mtn=flap.nface+wing.nface+aileron.nface+htail.nface+elevator.nface+vtail.nface+rudder.nface;
    
    /* Assign memory for global matrices */
    AN=(double *)calloc(mtn*mtn, sizeof(double)); /* Normal flow coefficient matrix */
    invAN=(double *)malloc(sizeof(double)*mtn*mtn); /* Inverse of normal flow coefficient matrix */
    BN=(double *)calloc(mtn*mtn, sizeof(double)); /* Downwash coefficient matrix */
    RHS=(double *)malloc(sizeof(double)*mtn); /* Right Hand Side vector */
    Gammas = (double *)malloc(sizeof(double)*mtn); /* Vorticity vector history */
    wind=(double *)malloc(sizeof(double)*mtn); /* Induced windspeeds */

    /* Calculate influence coefficient matrices */
    cycleliftsurf(pwing,pflap,paileron,phtail,pelevator,pvtail,prudder,AN,BN,mtn);    

    /* Invert the influence coefficient matrix */
    gauss(AN,invAN,mtn);

    /* Assign memory for gamma and wind values in every liftsurf */
    pwing->gamma=(double *)calloc(pwing->nface*ntimes, sizeof(double));
    pwing->wind=(double *)calloc(pwing->nface, sizeof(double));
    pflap->gamma=(double *)calloc(pflap->nface*ntimes, sizeof(double));
    pflap->wind=(double *)calloc(pflap->nface, sizeof(double));
    paileron->gamma=(double *)calloc(paileron->nface*ntimes, sizeof(double));
    paileron->wind=(double *)calloc(paileron->nface, sizeof(double));
    phtail->gamma=(double *)calloc(phtail->nface*ntimes, sizeof(double));
    phtail->wind=(double *)calloc(phtail->nface, sizeof(double));
    pelevator->gamma=(double *)calloc(pelevator->nface*ntimes, sizeof(double));
    pelevator->wind=(double *)calloc(pelevator->nface, sizeof(double));
    pvtail->gamma=(double *)calloc(pvtail->nface*ntimes, sizeof(double));
    pvtail->wind=(double *)calloc(pvtail->nface, sizeof(double));
    prudder->gamma=(double *)calloc(prudder->nface*ntimes, sizeof(double));
    prudder->wind=(double *)calloc(prudder->nface, sizeof(double));
    
    /* Assign memory AND INITIALISE uvw values in every liftsurf */
    pwing->uvw=(double *)calloc(pwing->nface * 3, sizeof(double));
    pflap->uvw=(double *)calloc(pflap->nface * 3, sizeof(double));
    paileron->uvw=(double *)calloc(paileron->nface * 3, sizeof(double));
    phtail->uvw = (double *)calloc(phtail->nface * 3, sizeof(double));
    pelevator->uvw=(double *)calloc(pelevator->nface * 3, sizeof(double));
    pvtail->uvw = (double *)calloc(pvtail->nface * 3, sizeof(double));
    prudder->uvw=(double *)calloc(prudder->nface * 3, sizeof(double));
    
    /* Assign memory for Deltap, Deltad and aeroforce on every liftsurf */
    pwing->Deltap=(double *)malloc(sizeof(double)*pwing->nface);
    pflap->Deltap=(double *)malloc(sizeof(double)*pflap->nface);
    paileron->Deltap=(double *)malloc(sizeof(double)*paileron->nface);
    phtail->Deltap=(double *)malloc(sizeof(double)*phtail->nface);
    pelevator->Deltap=(double *)malloc(sizeof(double)*pelevator->nface);
    pvtail->Deltap=(double *)malloc(sizeof(double)*pvtail->nface);
    prudder->Deltap=(double *)malloc(sizeof(double)*prudder->nface);
    pwing->Deltad=(double *)malloc(sizeof(double)*pwing->nface);
    pflap->Deltad=(double *)malloc(sizeof(double)*pflap->nface);
    paileron->Deltad=(double *)malloc(sizeof(double)*paileron->nface);
    phtail->Deltad=(double *)malloc(sizeof(double)*phtail->nface);
    pelevator->Deltad=(double *)malloc(sizeof(double)*pelevator->nface);
    pvtail->Deltad=(double *)malloc(sizeof(double)*pvtail->nface);
    prudder->Deltad=(double *)malloc(sizeof(double)*prudder->nface);
    pwing->aeroforce=(double *)malloc(sizeof(double)*4);
    pflap->aeroforce=(double *)malloc(sizeof(double)*4);
    paileron->aeroforce=(double *)malloc(sizeof(double)*4);
    phtail->aeroforce=(double *)malloc(sizeof(double)*4);
    pelevator->aeroforce=(double *)malloc(sizeof(double)*4);
    pvtail->aeroforce=(double *)malloc(sizeof(double)*4);
    prudder->aeroforce=(double *)malloc(sizeof(double)*4);
    
    for (it=0;it<ntimes-1;it++){
        printf("it=%i\n",it);
        if (it > 0){
            /* Calculate influences of all the wings on all the liftsurfs */
            calcwakeinf(pwing,pflap,paileron,phtail,pelevator,pvtail,prudder,it);
        }
        /* Calculate Right Hand Side vector */
        calcRHS(pwing,pflap,paileron,phtail,pelevator,pvtail,prudder,RHS,UVW);
        
        /* Solve for panel vorticity */
        for (i=0; i < mtn; i++) {
            *(Gammas+i)=0.0;
            for (j=0; j < mtn; j++) {
                *(Gammas+i) += *(invAN+i+j*mtn)* *(RHS+j);
            }
        }

        /* Calculate induced airspeeds */
        for (i=0; i < mtn; i++) {
            *(wind+i)=0.0;
            for (j=0; j < mtn; j++) {
                *(wind+i) += *(BN+i+j*mtn)* *(Gammas+j); /* Downwash by vortices on surface */
            }
        }
        
        /* Assign Gammas and wind onto respective lifting surfaces */
        assignGammas(Gammas,pwing,pflap,paileron,phtail,pelevator,pvtail,prudder,it);
        assignwind(wind,pwing,pflap,paileron,phtail,pelevator,pvtail,prudder);       
        
        /* Propagate wake vortex strength */
        propwakevort(pwing,dt,it,UVW);
        propwakevort(pflap,dt,it,UVW);
        propwakevort(paileron,dt,it,UVW);
        propwakevort(phtail,dt,it,UVW);
        propwakevort(pelevator,dt,it,UVW);
        propwakevort(pvtail,dt,it,UVW);
        propwakevort(prudder,dt,it,UVW);
        
        /* Assign vortex strength of trailing edge panels to the first wake panels */
        wakegamma(pwing,it);
        wakegamma(pflap,it);
        wakegamma(paileron,it);         
        wakegamma(phtail,it);         
        wakegamma(pelevator,it);         
        wakegamma(pvtail,it);         
        wakegamma(prudder,it);         
        
        /* Add influence of latest wake element on wind */
        addlastwind(pwing,pflap,paileron,it,3);
        addlastwind(phtail,pelevator,paileron,it,2); /* paileron is irrelevant */
        addlastwind(pvtail,prudder,paileron,it,2); /* paileron is irrelevant */
        
        /* Calculate the pressure difference across the panels in every liftsurf */
        calcforces(pwing,pwing,pflap,paileron,it,dt,UVW,rho,0,10000,20000);
        calcforces(pflap,pwing,pflap,paileron,it,dt,UVW,rho,0,10000,20000);
        calcforces(paileron,pwing,pflap,paileron,it,dt,UVW,rho,0,10000,20000);
        calcforceshtail(phtail,phtail,pelevator,it,dt,UVW,rho,30000,40000);
        calcforceshtail(pelevator,phtail,pelevator,it,dt,UVW,rho,30000,40000);
        calcforcesvtail(pvtail,pvtail,prudder,it,dt,UVW,rho,50000,60000);
        calcforcesvtail(prudder,pvtail,prudder,it,dt,UVW,rho,50000,60000);
       
        /* Calculate local airspeeds on wake */
        if (it>0 && freewake!=0){
            infonwake(pwing,pwing,it,0);
            infonwake(pwing,pflap,it,1);
            infonwake(pwing,paileron,it,1);
            infonwake(pwing,phtail,it,1);
            infonwake(pwing,pelevator,it,1);
            infonwake(pflap,pwing,it,0);
            infonwake(pflap,pflap,it,1);
            infonwake(pflap,paileron,it,1);
            infonwake(pflap,phtail,it,1);
            infonwake(pflap,pelevator,it,1);
            infonwake(paileron,pwing,it,0);
            infonwake(paileron,pflap,it,1);
            infonwake(paileron,paileron,it,1);
            infonwake(paileron,phtail,it,1);
            infonwake(paileron,pelevator,it,1);
            infonwake(phtail,pwing,it,0);
            infonwake(phtail,pflap,it,1);
            infonwake(phtail,paileron,it,1);
            infonwake(phtail,phtail,it,1);
            infonwake(phtail,pelevator,it,1);
            infonwake(pelevator,pwing,it,0);
            infonwake(pelevator,pflap,it,1);
            infonwake(pelevator,paileron,it,1);
            infonwake(pelevator,phtail,it,1);
            infonwake(pelevator,pelevator,it,1);
            /* vtail and rudder are treated as uncoupled */
            infonwake(pvtail,pvtail,it,0);
            infonwake(pvtail,prudder,it,1);
            infonwake(prudder,pvtail,it,0);
            infonwake(prudder,prudder,it,1);
       }

        /* Propagate wake position */
        propwakexyz(pwing,dt,it,UVW);
        propwakexyz(pflap,dt,it,UVW);
        propwakexyz(paileron,dt,it,UVW);
        propwakexyz(phtail,dt,it,UVW);
        propwakexyz(pelevator,dt,it,UVW);
        propwakexyz(pvtail,dt,it,UVW);
        propwakexyz(prudder,dt,it,UVW);
        
        /* Calculate the total forces on all the liftsurfs */
        for (i=0;i<4;i++){
            *(totalforce+i+it*4)=*(pwing->aeroforce+i)+*(pflap->aeroforce+i)+*(paileron->aeroforce+i)+*(phtail->aeroforce+i)+*(pelevator->aeroforce+i)+*(pvtail->aeroforce+i)+*(prudder->aeroforce+i);
        }
    }   
    
    it--;
    printf("forcex=%f, forcey=%f, forcez=%f, draginduced=%f\n",*(totalforce+0+it*4),*(totalforce+1+it*4),*(totalforce+2+it*4),*(totalforce+3+it*4));
    
    /* Print output */
    exportTextOutput(Outfile, ntimes, it, dt, totalforce, flap, pflap, aileron, paileron, wing, pwing, htail, phtail, vtail, pvtail,
                    elevator, pelevator, rudder, prudder, mtn, AN, invAN, BN, RHS, Gammas, wind);
}