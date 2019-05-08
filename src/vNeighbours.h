//
//  vNeighbours.h
//
//  Functions for finding neighbours in wake+lifting surface contexts
//
//
#ifndef vNeighbours_h
#define vNeighbours_h

#include "vLiftsurf.h"

void findneighbours(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, int cwing, int cflap, int caileron);

void calcdistwake(struct liftsurf *pwing, struct liftsurf *pflap, int i, int cflap);

void findallwakeneighbours(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, int cwing, int cflap, int caileron);

#endif /* vNeighbours_h */