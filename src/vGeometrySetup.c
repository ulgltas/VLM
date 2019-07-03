//
//  vGeometrySetup.c
//
//  Geometry setup function
//
//

#include "vLiftsurf.h"
#include "vGeometrySetup.h"
#include "vGeometry.h"
#include "vVortex.h"
#include <stdio.h>

int geometry_setup(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail,
                    struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder)
{
    printf("%p\n",pwing);
    printf("%p\n",pflap);
    printf("%p\n",paileron);
    printf("%i\n",pwing->nface);
    printf("%i\n",pflap->nface);
    printf("%i\n",paileron->nface);
    printf("%i\n",phtail->nface);
    printf("%i\n",pelevator->nface);
    printf("%i\n",pvtail->nface);
    printf("%i\n",prudder->nface);
    /* Calculate collocation points and vortex segment lengths */
    colvec(pflap);
    colvec(paileron);
    colvec(pwing);
    colvec(phtail);
    colvec(pelevator);
    colvec(pvtail);
    colvec(prudder);
    printf("Colvec\n");
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
    return pflap->nface+pwing->nface+paileron->nface+phtail->nface+pelevator->nface+pvtail->nface+prudder->nface;
}