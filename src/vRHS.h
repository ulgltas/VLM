//
//  vRHS.h
//
//  Calculate the values in the Right Hand Side vector
//
//
#ifndef vRHS_h
#define vRHS_h

#include "vLiftsurf.h"

void calcRHS(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder, double *RHS, double UVW[]);


#endif /* vRHS_h */