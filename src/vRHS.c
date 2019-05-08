//
//  vRHS.c
//
//  Calculate the values in the Right Hand Side vector
//
//

#include "vRHS.h"
#include "vLiftsurf.h"

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