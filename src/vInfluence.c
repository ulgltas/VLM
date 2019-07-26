//
//  vInfluence.c
//
//  Calculation of influence coefficients of different lifting surfaces
//
//

#include "vLiftsurf.h"
#include "vVLMData.h"
#include "vInfluence.h"
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

void cycleliftsurf(struct VLMData *data)
{
    int m,n;
    
    /* Cycle through the lifting surfaces in the order wing, flap, aileron, htail, elevator, vtail, rudder */
    /* Treat the vtal and rudder problem as uncoupled */
    m=0;
    n=0;
    /* Influence on the wing from the wing */
    infcoeff(&(data->wing),&(data->wing),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the wing from the flap */
    n+=(data->wing).nface;
    infcoeff(&(data->wing),&(data->flap),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the wing from the aileron */
    n+=(data->flap).nface;
    infcoeff(&(data->wing),&(data->aileron),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the wing from the horizontal tail */
    n+=(data->aileron).nface;
    infcoeff(&(data->wing),&(data->htail),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the wing from the elevator */
    n+=(data->htail).nface;
    infcoeff(&(data->wing),&(data->elevator),data->AN,data->BN,m,n,data->mtn);
//     /* Influence on the wing from the vertical tail */
//     n+=(data->elevator).nface;
//     infcoeff(&(data->wing),&(data->vtail),data->AN,data->BN,m,n,data->mtn);
//     /* Influence on the wing from the rudder */
//     n+=(data->vtail).nface;
//     infcoeff(&(data->wing),&(data->rudder),data->AN,data->BN,m,n,data->mtn);

    m+=(data->wing).nface;
    n=0;
    /* Influence on the flap from the wing */
    infcoeff(&(data->flap),&(data->wing),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the flap from the flap */
    n+=(data->wing).nface;
    infcoeff(&(data->flap),&(data->flap),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the flap from the aileron */
    n+=(data->flap).nface;
    infcoeff(&(data->flap),&(data->aileron),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the flap from the horizontal tail */
    n+=(data->aileron).nface;
    infcoeff(&(data->flap),&(data->htail),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the flap from the elevator */
    n+=(data->htail).nface;
    infcoeff(&(data->flap),&(data->elevator),data->AN,data->BN,m,n,data->mtn);
//     /* Influence on the flap from the vertical tail */
//     n+=(data->elevator).nface;
//     infcoeff(&(data->flap),&(data->vtail),data->AN,data->BN,m,n,data->mtn);
//     /* Influence on the flap from the rudder */
//     n+=(data->vtail).nface;
//     infcoeff(&(data->flap),&(data->rudder),data->AN,data->BN,m,n,data->mtn);
    
    m+=(data->flap).nface;
    n=0;
    /* Influence on the aileron from the wing */
    infcoeff(&(data->aileron),&(data->wing),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the aileron from the flap */
    n+=(data->wing).nface;
    infcoeff(&(data->aileron),&(data->flap),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the aileron from the aileron */
    n+=(data->flap).nface;
    infcoeff(&(data->aileron),&(data->aileron),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the aileron from the horizontal tail */
    n+=(data->aileron).nface;
    infcoeff(&(data->aileron),&(data->htail),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the aileron from the elevator */
    n+=(data->htail).nface;
    infcoeff(&(data->aileron),&(data->elevator),data->AN,data->BN,m,n,data->mtn);
//     /* Influence on the aileron from the vertical tail */
//     n+=(data->elevator).nface;
//     infcoeff(&(data->aileron),&(data->vtail),data->AN,data->BN,m,n,data->mtn);
//     /* Influence on the aileron from the rudder */
//     n+=(data->vtail).nface;
//     infcoeff(&(data->aileron),&(data->rudder),data->AN,data->BN,m,n,data->mtn);
    
    m+=(data->aileron).nface;
    n=0;
    /* Influence on the horizontal tail from the wing */
    infcoeff(&(data->htail),&(data->wing),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the horizontal tail from the flap */
    n+=(data->wing).nface;
    infcoeff(&(data->htail),&(data->flap),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the horizontal tail from the aileron */
    n+=(data->flap).nface;
    infcoeff(&(data->htail),&(data->aileron),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the horizontal tail from the horizontal tail */
    n+=(data->aileron).nface;
    infcoeff(&(data->htail),&(data->htail),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the horizontal tail from the elevator */
    n+=(data->htail).nface;
    infcoeff(&(data->htail),&(data->elevator),data->AN,data->BN,m,n,data->mtn);
//     /* Influence on the horizontal tail from the vertical tail */
//     n+=(data->elevator).nface;
//     infcoeff(&(data->htail),&(data->vtail),data->AN,data->BN,m,n,data->mtn);
//     /* Influence on the horizontal tail from the rudder */
//     n+=(data->vtail).nface;
//     infcoeff(&(data->htail),&(data->rudder),data->AN,data->BN,m,n,data->mtn);
    
    m+=(data->htail).nface;
    n=0;
    /* Influence on the elevator from the wing */
    infcoeff(&(data->elevator),&(data->wing),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the elevator from the flap */
    n+=(data->wing).nface;
    infcoeff(&(data->elevator),&(data->flap),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the elevator from the aileron */
    n+=(data->flap).nface;
    infcoeff(&(data->elevator),&(data->aileron),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the elevator from the horizontal tail */
    n+=(data->aileron).nface;
    infcoeff(&(data->elevator),&(data->htail),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the elevator from the elevator */
    n+=(data->htail).nface;
    infcoeff(&(data->elevator),&(data->elevator),data->AN,data->BN,m,n,data->mtn);
//     /* Influence on the elevator from the vertical tail */
//     n+=(data->elevator).nface;
//     infcoeff(&(data->elevator),&(data->vtail),data->AN,data->BN,m,n,data->mtn);
//     /* Influence on the elevator from the rudder */
//     n+=(data->vtail).nface;
//     infcoeff(&(data->elevator),&(data->rudder),data->AN,data->BN,m,n,data->mtn);

    m+=(data->elevator).nface;
    n=0;
    /* Influence on the vertical tail from the wing */
//    infcoeff(&(data->vtail),&(data->wing),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the vertical tail from the flap */
    n+=(data->wing).nface;
//    infcoeff(&(data->vtail),&(data->flap),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the vertical tail from the aileron */
    n+=(data->flap).nface;
//    infcoeff(&(data->vtail),&(data->aileron),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the vertical tail from the horizontal tail */
    n+=(data->aileron).nface;
//    infcoeff(&(data->vtail),&(data->htail),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the vertical tail from the elevator */
    n+=(data->htail).nface;
//    infcoeff(&(data->vtail),&(data->elevator),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the vertical tail from the vertical tail */
    n+=(data->elevator).nface;
    infcoeff(&(data->vtail),&(data->vtail),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the vertical tail from the rudder */
    n+=(data->vtail).nface;
    infcoeff(&(data->vtail),&(data->rudder),data->AN,data->BN,m,n,data->mtn);
    
    m+=(data->vtail).nface;
    n=0;
    /* Influence on the rudder from the wing */
//    infcoeff(&(data->rudder),&(data->wing),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the rudder from the flap */
    n+=(data->wing).nface;
//    infcoeff(&(data->rudder),&(data->flap),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the rudder from the aileron */
    n+=(data->flap).nface;
//    infcoeff(&(data->rudder),&(data->aileron),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the rudder from the horizontal tail */
    n+=(data->aileron).nface;
//    infcoeff(&(data->rudder),&(data->htail),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the rudder from the elevator */
    n+=(data->htail).nface;
//    infcoeff(&(data->rudder),&(data->elevator),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the rudder from the vertical tail */
    n+=(data->elevator).nface;
    infcoeff(&(data->rudder),&(data->vtail),data->AN,data->BN,m,n,data->mtn);
    /* Influence on the rudder from the rudder */
    n+=(data->vtail).nface;
    infcoeff(&(data->rudder),&(data->rudder),data->AN,data->BN,m,n,data->mtn);
}