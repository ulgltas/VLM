//
//  vInfluence.c
//
//  Calculation of influence coefficients of different lifting surfaces
//
//

#include "vLiftsurf.h"
#include "vVortex.h"

void infcoeff(struct liftsurf *plift1, struct liftsurf *plift2, double *AN, double *BN, int istart, int jstart, int mtn)
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