//
//  vForce.h
//
//  Calculation of forces on different lifting surfaces
//
//

#ifndef vForce_h
#define vForce_h

#include "vLiftsurf.h"

void calcforces(struct liftsurf *plift, struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, int it, double dt, double UVW[], double rho, int cwing, int cflap, int caileron);

void calcforceshtail(struct liftsurf *plift, struct liftsurf *phtail, struct liftsurf *pelevator, int it, double dt, double UVW[], double rho, int chtail, int celevator);

void calcforcesvtail(struct liftsurf *plift, struct liftsurf *pvtail, struct liftsurf *prudder, int it, double dt, double UVW[], double rho, int cvtail, int crudder);

#endif /* vForce_h */