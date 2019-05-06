//
//  vWake.h
//
//  Wake related functions
//
//

#ifndef vWake_h
#define vWake_h

#include "vLiftsurf.h"

void storewakeinds(struct liftsurf *pwing);

int findwakeneighbours(struct liftsurf *pwing, int ipanel, int lrneighs[],int cwing);

void setupwakes(struct liftsurf *pwing, int cwing);

void createwake(struct liftsurf *pwing, double cwing,int ntimes);

void correspshedwake(struct liftsurf *pwing);

void shedwake(struct liftsurf *pwing);

void calcdistwake(struct liftsurf *pwing, struct liftsurf *pflap, int i, int cflap);

void findallwakeneighbours(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, int cwing, int cflap, int caileron);

void wakegamma(struct liftsurf *pwing, int it);

void propwakexyz(struct liftsurf *pwing, double dt, int it, double UVW[]);

void propwakevort(struct liftsurf *pwing, double dt, int it, double UVW[]);

void infonwake(struct liftsurf *pwing, struct liftsurf *pflap, int it, int summode);

void wakeinf(struct liftsurf *pwing, struct liftsurf *pflap, int nw, int it, int summode);

void calcwakeinf(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder, int it);

void addlastwind(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, int it, int narg);

#endif /* vWake_h */