//
//  vInfluence.h
//
//  Calculation of influence coefficients of different lifting surfaces
//
//
#ifndef vInfluence_h
#define vInfluence_h

#include "vLiftsurf.h"

void infcoeff(struct liftsurf *plift1, struct liftsurf *plift2, double *AN, double *BN, int istart, int jstart,int mtn);

void cycleliftsurf(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder, double *AN, double *BN, int mtn);

#endif /* vInfluence_h */