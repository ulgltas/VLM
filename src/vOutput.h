//
//  vOutput.h
//
//  Data output
//
//
#ifndef vOutput_h
#define vOutput_h

#include "vLiftsurf.h"
#include <stdio.h>

void exportTextOutput(char *Outfile, int ntimes, int it, double dt, double *totalforce, struct liftsurf flap, struct liftsurf *pflap,
                    struct liftsurf aileron, struct liftsurf *paileron, struct liftsurf wing, struct liftsurf *pwing,
                    struct liftsurf htail, struct liftsurf *phtail, struct liftsurf vtail, struct liftsurf *pvtail,
                    struct liftsurf elevator, struct liftsurf *pelevator, struct liftsurf rudder, struct liftsurf *prudder,
                    int mtn, double *AN, double *invAN, double *BN, double *RHS, double *Gammas, double *wind);

#endif /* vOutput_h */