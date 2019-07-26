//
//  vSetup.c
//
//  Setup functions
//
//

#include "vLiftsurf.h"
#include "vVLMData.h"
#include "vSetup.h"
#include "vControl.h"
#include "vWing.h"
#include "vHTail.h"
#include "vVTail.h"
#include "vWake.h"
#include "vInput.h"
#include "vNeighbours.h"
#include <stdio.h>
#include <stdlib.h>

void setup(char* Infile, struct VLMData *data)
{
    int m, n, mht, nht, mvt, nvt;
    double yaw, timestep_denom;
    char Wngfile[70], HTailfile[70], VTailfile[70];
    int ChkAil, ChkWTED, ChkElev, ChkRdr;
    importInputFile(Infile, (data->UVW), &(data->rho), &(data->aoa), &yaw, &m, &mht, &mvt, &n, &nht, &nvt, &(data->ntimes), &timestep_denom, &(data->freewake),
                    (data->delta), (data->beta), (data->eta), (data->zeta), Wngfile, HTailfile, VTailfile);
        
    
    /* Create the lifting surfaces that make up the wing */
    data->MAC=wingsetup(&(data->flap),&(data->aileron),&(data->wing),Wngfile,m,n, &ChkAil, &ChkWTED);
    /* Create the lifting surfaces that make up the horizontal tail */
    htailsetup(&(data->htail),&(data->elevator),HTailfile,mht,nht, &ChkElev);
    /* Create the lifting surfaces that make up the vertical tail */
    vtailsetup(&(data->vtail),&(data->rudder),VTailfile,mht,nht, &ChkRdr);
    data->dt=data->MAC/data->UVW[0]/timestep_denom; /* dt is based on the length of the wing's Mean Aerodynamic Chord */
    printf("MAC=%f, dt=%f\n",data->MAC,data->dt);
    
    /* Rotate ailerons */
    if (ChkAil == 1)
    {
        rotateail(&(data->aileron), (data->delta));
    }
    
    /* Rotate flaps */
    if (ChkWTED == 1)
    {
        rotateail(&(data->flap), (data->beta));
    }
    
    /* Rotate elevator */
    if (ChkElev == 1)
    {
        rotateail(&(data->elevator),(data->eta));
    }
    /* Rotate rudder */
    if (ChkRdr == 1)
    {
        rotateail(&(data->rudder),(data->zeta));
    }

    findneighbours(&(data->wing),&(data->flap),&(data->aileron),0,10000,20000);
    findneighbours(&(data->flap),&(data->wing),&(data->aileron),10000,0,20000);
    findneighbours(&(data->aileron),&(data->wing),&(data->flap),20000,0,10000);
    findneighbours(&(data->htail),&(data->elevator),&(data->aileron),30000,40000,20000); /* We don't care about &(data->aileron) being called here */
    findneighbours(&(data->elevator),&(data->htail),&(data->aileron),40000,30000,20000); /* We don't care about &(data->aileron) being called here */
    findneighbours(&(data->vtail),&(data->rudder),&(data->aileron),50000,60000,20000); /* We don't care about &(data->aileron) being called here */
    findneighbours(&(data->rudder),&(data->vtail),&(data->aileron),60000,50000,20000); /* We don't care about &(data->aileron) being called here */

    /* Create the wake */
    createwake(&(data->wing),0,(data->ntimes));
    createwake(&(data->flap),10000,(data->ntimes));
    createwake(&(data->aileron),20000,(data->ntimes));
    createwake(&(data->htail),30000,(data->ntimes));
    createwake(&(data->elevator),40000,(data->ntimes));
    createwake(&(data->vtail),50000,(data->ntimes));
    createwake(&(data->rudder),60000,(data->ntimes));

    /* Find which bound vortex panels correspond to which wake panel vertices */
    correspshedwake(&(data->wing));
    correspshedwake(&(data->flap));
    correspshedwake(&(data->aileron));
    correspshedwake(&(data->htail));
    correspshedwake(&(data->elevator));
    correspshedwake(&(data->vtail));
    correspshedwake(&(data->rudder));

    /* Shed first wake element */
    shedwake(&(data->wing));
    shedwake(&(data->flap));
    shedwake(&(data->aileron));
    shedwake(&(data->htail));
    shedwake(&(data->elevator));
    shedwake(&(data->vtail));
    shedwake(&(data->rudder));
    findallwakeneighbours(&(data->wing),&(data->flap),&(data->aileron),0,10000,20000);
    findallwakeneighbours(&(data->htail),&(data->elevator),&(data->aileron),30000,40000,20000); /* We don't care about &(data->aileron) being called here */
    findallwakeneighbours(&(data->vtail),&(data->rudder),&(data->aileron),50000,60000,20000); /* We don't care about &(data->aileron) being called here */
    /* Calculate total number of panels */
    data->mtn = (data->flap).nface+(data->wing).nface+(data->aileron).nface+(data->htail).nface+(data->elevator).nface+(data->vtail).nface+(data->rudder).nface;
    /* Assign memory for global matrices */
    data->AN=(double *)calloc(data->mtn*data->mtn, sizeof(double)); /* Normal flow coefficient matrix */
    data->invAN=(double *)malloc(sizeof(double)*data->mtn*data->mtn); /* Inverse of normal flow coefficient matrix */
    data->BN=(double *)calloc(data->mtn*data->mtn, sizeof(double)); /* Downwash coefficient matrix */
    data->RHS=(double *)malloc(sizeof(double)*data->mtn); /* Right Hand Side vector */
    data->Gammas = (double *)malloc(sizeof(double)*data->mtn); /* Vorticity vector history */
    data->wind=(double *)malloc(sizeof(double)*data->mtn); /* Induced windspeeds */

    data->totalforce=(double *)calloc(data->ntimes*4, sizeof(double));
}

void memory_setup(struct VLMData *data)
{
    /* Assign memory for gamma and wind values in every liftsurf */
    (data->wing).gamma=(double *)calloc((data->wing).nface*(data->ntimes), sizeof(double));
    (data->wing).wind=(double *)calloc((data->wing).nface, sizeof(double));
    (data->flap).gamma=(double *)calloc((data->flap).nface*(data->ntimes), sizeof(double));
    (data->flap).wind=(double *)calloc((data->flap).nface, sizeof(double));
    (data->aileron).gamma=(double *)calloc((data->aileron).nface*(data->ntimes), sizeof(double));
    (data->aileron).wind=(double *)calloc((data->aileron).nface, sizeof(double));
    (data->htail).gamma=(double *)calloc((data->htail).nface*(data->ntimes), sizeof(double));
    (data->htail).wind=(double *)calloc((data->htail).nface, sizeof(double));
    (data->elevator).gamma=(double *)calloc((data->elevator).nface*(data->ntimes), sizeof(double));
    (data->elevator).wind=(double *)calloc((data->elevator).nface, sizeof(double));
    (data->vtail).gamma=(double *)calloc((data->vtail).nface*(data->ntimes), sizeof(double));
    (data->vtail).wind=(double *)calloc((data->vtail).nface, sizeof(double));
    (data->rudder).gamma=(double *)calloc((data->rudder).nface*(data->ntimes), sizeof(double));
    (data->rudder).wind=(double *)calloc((data->rudder).nface, sizeof(double));
    
    /* Assign memory AND INITIALISE uvw values in every liftsurf */
    (data->wing).uvw=(double *)calloc((data->wing).nface * 3, sizeof(double));
    (data->flap).uvw=(double *)calloc((data->flap).nface * 3, sizeof(double));
    (data->aileron).uvw=(double *)calloc((data->aileron).nface * 3, sizeof(double));
    (data->htail).uvw = (double *)calloc((data->htail).nface * 3, sizeof(double));
    (data->elevator).uvw=(double *)calloc((data->elevator).nface * 3, sizeof(double));
    (data->vtail).uvw = (double *)calloc((data->vtail).nface * 3, sizeof(double));
    (data->rudder).uvw=(double *)calloc((data->rudder).nface * 3, sizeof(double));
    
    /* Assign memory for Deltap, Deltad and aeroforce on every liftsurf */
    (data->wing).Deltap=(double *)malloc(sizeof(double)*(data->wing).nface);
    (data->flap).Deltap=(double *)malloc(sizeof(double)*(data->flap).nface);
    (data->aileron).Deltap=(double *)malloc(sizeof(double)*(data->aileron).nface);
    (data->htail).Deltap=(double *)malloc(sizeof(double)*(data->htail).nface);
    (data->elevator).Deltap=(double *)malloc(sizeof(double)*(data->elevator).nface);
    (data->vtail).Deltap=(double *)malloc(sizeof(double)*(data->vtail).nface);
    (data->rudder).Deltap=(double *)malloc(sizeof(double)*(data->rudder).nface);
    (data->wing).Deltad=(double *)malloc(sizeof(double)*(data->wing).nface);
    (data->flap).Deltad=(double *)malloc(sizeof(double)*(data->flap).nface);
    (data->aileron).Deltad=(double *)malloc(sizeof(double)*(data->aileron).nface);
    (data->htail).Deltad=(double *)malloc(sizeof(double)*(data->htail).nface);
    (data->elevator).Deltad=(double *)malloc(sizeof(double)*(data->elevator).nface);
    (data->vtail).Deltad=(double *)malloc(sizeof(double)*(data->vtail).nface);
    (data->rudder).Deltad=(double *)malloc(sizeof(double)*(data->rudder).nface);
    (data->wing).aeroforce=(double *)malloc(sizeof(double)*4);
    (data->flap).aeroforce=(double *)malloc(sizeof(double)*4);
    (data->aileron).aeroforce=(double *)malloc(sizeof(double)*4);
    (data->htail).aeroforce=(double *)malloc(sizeof(double)*4);
    (data->elevator).aeroforce=(double *)malloc(sizeof(double)*4);
    (data->vtail).aeroforce=(double *)malloc(sizeof(double)*4);
    (data->rudder).aeroforce=(double *)malloc(sizeof(double)*4);
}