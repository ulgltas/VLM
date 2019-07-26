//
//  vIteration.c
//
//  Run VLM iterations
//
//

#include "vLiftsurf.h"
#include "vVLMData.h"
#include "vIteration.h"
#include "vAssignValues.h"
#include "vForce.h"
#include "vRHS.h"
#include "vWake.h"
#include <stdio.h>

void iteration(struct VLMData *data, int it)
{
    int i, j;
    printf("it=%i\n",it);
    if (it > 0){
        /* Calculate influences of all the wings on all the liftsurfs */
        calcwakeinf(data, it);
    }
    /* Calculate Right Hand Side vector */
    calcRHS(data);
    
    /* Solve for panel vorticity */
    for (i=0; i < (data->mtn); i++) {
        *((data->Gammas)+i)=0.0;
        for (j=0; j < (data->mtn); j++) {
            *((data->Gammas)+i) += *((data->invAN)+i+j*(data->mtn))* *((data->RHS)+j);
        }
    }

    /* Calculate induced airspeeds */
    for (i=0; i < (data->mtn); i++) {
        *((data->wind)+i)=0.0;
        for (j=0; j < (data->mtn); j++) {
            *((data->wind)+i) += *((data->BN)+i+j*(data->mtn))* *((data->Gammas)+j); /* Downwash by vortices on surface */
        }
    }
    
    /* Assign Gammas and wind onto respective lifting surfaces */
    assignGammas(data, it);
    assignwind(data);
    
    /* Propagate wake vortex strength */
    propwakevort(&(data->wing),(data->dt),it,(data->UVW));
    propwakevort(&(data->flap),(data->dt),it,(data->UVW));
    propwakevort(&(data->aileron),(data->dt),it,(data->UVW));
    propwakevort(&(data->htail),(data->dt),it,(data->UVW));
    propwakevort(&(data->elevator),(data->dt),it,(data->UVW));
    propwakevort(&(data->vtail),(data->dt),it,(data->UVW));
    propwakevort(&(data->rudder),(data->dt),it,(data->UVW));
    
    /* Assign vortex strength of trailing edge panels to the first wake panels */
    wakegamma(&(data->wing),it);
    wakegamma(&(data->flap),it);
    wakegamma(&(data->aileron),it);
    wakegamma(&(data->htail),it);
    wakegamma(&(data->elevator),it);
    wakegamma(&(data->vtail),it);
    wakegamma(&(data->rudder),it);
    
    /* Add influence of latest wake element on wind */
    addlastwind(&(data->wing),&(data->flap),&(data->aileron),it,3);
    addlastwind(&(data->htail),&(data->elevator),&(data->aileron),it,2); /* &(data->aileron) is irrelevant */
    addlastwind(&(data->vtail),&(data->rudder),&(data->aileron),it,2); /* &(data->aileron) is irrelevant */
    
    /* Calculate the pressure difference across the panels in every liftsurf */
    calcforces(&(data->wing),&(data->wing),&(data->flap),&(data->aileron),it,(data->dt),(data->UVW),(data->rho),0,10000,20000);
    calcforces(&(data->flap),&(data->wing),&(data->flap),&(data->aileron),it,(data->dt), (data->UVW), (data->rho),0,10000,20000);
    calcforces(&(data->aileron),&(data->wing),&(data->flap),&(data->aileron),it,(data->dt), (data->UVW), (data->rho),0,10000,20000);
    calcforceshtail(&(data->htail),&(data->htail),&(data->elevator),it,(data->dt), (data->UVW), (data->rho),30000,40000);
    calcforceshtail(&(data->elevator),&(data->htail),&(data->elevator),it,(data->dt), (data->UVW), (data->rho),30000,40000);
    calcforcesvtail(&(data->vtail),&(data->vtail),&(data->rudder),it,(data->dt), (data->UVW), (data->rho),50000,60000);
    calcforcesvtail(&(data->rudder),&(data->vtail),&(data->rudder),it,(data->dt), (data->UVW), (data->rho),50000,60000);
    
    /* Calculate local airspeeds on wake */
    if (it>0 && data->freewake!=0){
        infonwake(&(data->wing),&(data->wing),it,0);
        infonwake(&(data->wing),&(data->flap),it,1);
        infonwake(&(data->wing),&(data->aileron),it,1);
        infonwake(&(data->wing),&(data->htail),it,1);
        infonwake(&(data->wing),&(data->elevator),it,1);
        infonwake(&(data->flap),&(data->wing),it,0);
        infonwake(&(data->flap),&(data->flap),it,1);
        infonwake(&(data->flap),&(data->aileron),it,1);
        infonwake(&(data->flap),&(data->htail),it,1);
        infonwake(&(data->flap),&(data->elevator),it,1);
        infonwake(&(data->aileron),&(data->wing),it,0);
        infonwake(&(data->aileron),&(data->flap),it,1);
        infonwake(&(data->aileron),&(data->aileron),it,1);
        infonwake(&(data->aileron),&(data->htail),it,1);
        infonwake(&(data->aileron),&(data->elevator),it,1);
        infonwake(&(data->htail),&(data->wing),it,0);
        infonwake(&(data->htail),&(data->flap),it,1);
        infonwake(&(data->htail),&(data->aileron),it,1);
        infonwake(&(data->htail),&(data->htail),it,1);
        infonwake(&(data->htail),&(data->elevator),it,1);
        infonwake(&(data->elevator),&(data->wing),it,0);
        infonwake(&(data->elevator),&(data->flap),it,1);
        infonwake(&(data->elevator),&(data->aileron),it,1);
        infonwake(&(data->elevator),&(data->htail),it,1);
        infonwake(&(data->elevator),&(data->elevator),it,1);
        /* vtail and rudder are treated as uncoupled */
        infonwake(&(data->vtail),&(data->vtail),it,0);
        infonwake(&(data->vtail),&(data->rudder),it,1);
        infonwake(&(data->rudder),&(data->vtail),it,0);
        infonwake(&(data->rudder),&(data->rudder),it,1);
    }

    /* Propagate wake position */
    propwakexyz(&(data->wing), (data->dt), it, (data->UVW));
    propwakexyz(&(data->flap), (data->dt), it, (data->UVW));
    propwakexyz(&(data->aileron), (data->dt), it, (data->UVW));
    propwakexyz(&(data->htail), (data->dt), it, (data->UVW));
    propwakexyz(&(data->elevator), (data->dt), it, (data->UVW));
    propwakexyz(&(data->vtail), (data->dt), it, (data->UVW));
    propwakexyz(&(data->rudder), (data->dt), it, (data->UVW));
    
    /* Calculate the total forces on all the liftsurfs */
    for (i=0;i<4;i++){
        *((data->totalforce)+i+it*4)=*((data->wing).aeroforce+i)+*((data->flap).aeroforce+i)+*((data->aileron).aeroforce+i)+*((data->htail).aeroforce+i)+*((data->elevator).aeroforce+i)+*((data->vtail).aeroforce+i)+*((data->rudder).aeroforce+i);
    }
    printf("forcex=%f, forcey=%f, forcez=%f, draginduced=%f\n",*((data->totalforce)+0+it*4),*((data->totalforce)+1+it*4),*((data->totalforce)+2+it*4),*((data->totalforce)+3+it*4));
}