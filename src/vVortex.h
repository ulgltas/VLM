//
//  vVortex.h
//
//  Vortex related functions
//
//

#ifndef vVortex_h
#define vVortex_h

#include "vLiftsurf.h"

void colvec(struct liftsurf *pflap);

void vortexblob(double uvw[], double x, double y, double z, double x1, double y1, double z1,
        double x2, double y2, double z2, double gama);

void vortex(double uvw[], double x, double y, double z, double x1, double y1, double z1,
        double x2, double y2, double z2, double gama);

#endif /* vVortex_h */