//
//  vPanel.h
//  
//
//  Functions related to panel creation
//
//

#ifndef vPanel_h
#define vPanel_h

#include <stdio.h>
#include "vLiftsurf.h"

int countfaces(int *ijflap, int nflap, int mp1, int np1);

void arrangefaces(int *ijflap, int mp1, int np1, struct liftsurf *pflap);

void vortexpanel(double *xv, double *yv, double *zv, double *xp, double *yp, double *zp, double dxw, int m, int n);

#endif /* vPanel_h */
