//
//  vInfluence.h
//
//  Calculation of influence coefficients of different lifting surfaces
//
//
#ifndef vInfluence_h
#define vInfluence_h

#include "vLiftsurf.h"
#include "vVLMData.h"

void infcoeff(struct liftsurf *plift1, struct liftsurf *plift2, double *AN, double *BN, int istart, int jstart,int mtn);

void cycleliftsurf(struct VLMData *data);

#endif /* vInfluence_h */