//
//  vSetup.c
//
//  Setup functions
//
//

#include "vLiftsurf.h"
#include "vSetup.h"
#include "vControl.h"
#include "vWing.h"
#include "vHTail.h"
#include "vVTail.h"
#include "vWake.h"
#include "vInput.h"
#include "vNeighbours.h"
#include <stdio.h>

double setup(char* Infile, double UVW[3], double *rho, int *ntimes, int *freewake, struct liftsurf *pwing, struct liftsurf *pflap,
             struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder)
{
    int m, n, mht, nht, mvt, nvt;
    double delta[2],beta[2],eta[2],zeta[2];
    double aoa, yaw, timestep_denom, MAC, dt;
    char Wngfile[60], HTailfile[60], VTailfile[60];
    int ChkAil, ChkWTED, ChkElev, ChkRdr;
    importInputFile(Infile, UVW, rho, &aoa, &yaw, &m, &mht, &mvt, &n, &nht, &nvt, ntimes, &timestep_denom, freewake,
                    delta, beta, eta, zeta, Wngfile, HTailfile, VTailfile);
        
    
    /* Create the lifting surfaces that make up the wing */
    MAC=wingsetup(pflap,paileron,pwing,Wngfile,m,n, &ChkAil, &ChkWTED);
    /* Create the lifting surfaces that make up the horizontal tail */
    htailsetup(phtail,pelevator,HTailfile,mht,nht, &ChkElev);
    /* Create the lifting surfaces that make up the vertical tail */
    vtailsetup(pvtail,prudder,VTailfile,mht,nht, &ChkRdr);
    dt=MAC/UVW[0]/timestep_denom; /* dt is based on the length of the wing's Mean Aerodynamic Chord */
    printf("MAC=%f, dt=%f\n",MAC,dt);
    
    /* Rotate ailerons */
    if (ChkAil == 1)
    {
        rotateail(paileron, delta);
    }
    
    /* Rotate flaps */
    if (ChkWTED == 1)
    {
        rotateail(pflap, beta);
    }
    
    /* Rotate elevator */
    if (ChkElev == 1)
    {
        rotateail(pelevator,eta);
    }
    /* Rotate rudder */
    if (ChkRdr == 1)
    {
        rotateail(prudder,zeta);
    }

    findneighbours(pwing,pflap,paileron,0,10000,20000);
    findneighbours(pflap,pwing,paileron,10000,0,20000);
    findneighbours(paileron,pwing,pflap,20000,0,10000);
    findneighbours(phtail,pelevator,paileron,30000,40000,20000); /* We don't care about paileron being called here */
    findneighbours(pelevator,phtail,paileron,40000,30000,20000); /* We don't care about paileron being called here */
    findneighbours(pvtail,prudder,paileron,50000,60000,20000); /* We don't care about paileron being called here */
    findneighbours(prudder,pvtail,paileron,60000,50000,20000); /* We don't care about paileron being called here */

    /* Create the wake */
    createwake(pwing,0,*ntimes);
    createwake(pflap,10000,*ntimes);
    createwake(paileron,20000,*ntimes);
    createwake(phtail,30000,*ntimes);
    createwake(pelevator,40000,*ntimes);
    createwake(pvtail,50000,*ntimes);
    createwake(prudder,60000,*ntimes);

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
    printf("%i\n",pwing->nface);
    printf("%i\n",pflap->nface);
    printf("%i\n",paileron->nface);
    printf("%i\n",phtail->nface);
    printf("%i\n",pelevator->nface);
    printf("%i\n",pvtail->nface);
    printf("%i\n",prudder->nface);
    return dt;
}