//
//  vAssignValues.c
//
//  Assign values back to the Liftsurf struct
//
//

#include "vLiftsurf.h"
#include "vAssignValues.h"

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