//
//  vVLMData.h
//  
//
//  VLMData struct definition
//
//

#ifndef vVLMData_h
#define vVLMData_h

#include "vLiftsurf.h"
struct VLMData
{
    struct liftsurf flap;
    struct liftsurf aileron;
    struct liftsurf wing;
    struct liftsurf htail;
    struct liftsurf elevator;
    struct liftsurf vtail;
    struct liftsurf rudder;
    double UVW[3];
    double aoa; // Angle of attack in radians
    double rho;
    double MAC;
    double dt;
    int it;
    int ntimes;
    int freewake;
    int mtn;
    // Control surfaces' angles (radians)
    double delta[2];
    double beta[2];
    double eta[2];
    double zeta[2];
    // Data that change with the movement (should be SWIG immutable)
    double *invAN,*AN,*BN,*RHS;
    double *Gammas, *wind;
    double *totalforce; // Forces over the whole thing
};

#endif /* vVLMData_h */
