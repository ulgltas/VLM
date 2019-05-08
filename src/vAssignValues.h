//
//  vAssignValues.h
//
//  Assign values back to the Liftsurf struct
//
//
#ifndef vAssignValues_h
#define vAssignValues_h

#include "vLiftsurf.h"

void assignGammas(double *Gammas, struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder, int it);

void assignwind(double *wind, struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder);

#endif /* vAssignValues_h */