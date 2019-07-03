//
//  vSetup.h
//
//  Setup functions
//
//
#ifndef vSetup_h
#define vSetup_h
double setup(char* Infile, double UVW[3], double *rho, int *ntimes, int *freewake, struct liftsurf *pwing, struct liftsurf *pflap,
             struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder);
#endif /* vSetup_h */