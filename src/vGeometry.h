//
//  vGeometry.h
//
//  Geometrical functions for the panels
//
//

#ifndef vGeometry_h
#define vGeometry_h

#include "vLiftsurf.h"

void cross(double xcrossy[], double x[], double y[]);

void normals(struct liftsurf *pflap);

void tangentials(struct liftsurf *pflap);

#endif /* vGeometry_h */