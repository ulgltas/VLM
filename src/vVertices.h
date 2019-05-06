//
//  vVertices.h
//
//  Vertex treatment functions
//
//

#ifndef vVertices_h
#define vVertices_h

#include "vLiftsurf.h"

void assignvertices(struct liftsurf *pflap, double *xflap, double *yflap, double *zflap, double *xvflap, double *yvflap, double *zvflap);

void assignvertices_fin(struct liftsurf *pflap, double *xflap, double *yflap, double *zflap, double *xvflap, double *yvflap, double *zvflap, int OptVTFusMounted, int OptVTTailMounted, int OptVTWingMounted);

#endif /* vVertices_h */