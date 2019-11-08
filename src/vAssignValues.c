//
//  vAssignValues.c
//
//  Assign values back to the Liftsurf struct
//
//

#include "vLiftsurf.h"
#include "vAssignValues.h"
#include "vVLMData.h"

void assignGammas(struct VLMData *data, int it)
{
    /* Assign the vortex strengths in vector Gammas back onto the panels of all the lifting surfaces */
    int i,m;
        
    /* Assign wing vortex strengths */
    m=0;
    for (i=0;i<(data->wing).nface;i++){
        *((data->wing).gamma+i+it*(data->wing).nface)=*((data->Gammas)+i);
    }
    /* Assign flap vortex strengths */
    m+=(data->wing).nface;
    for (i=0;i<(data->flap).nface;i++){
        *((data->flap).gamma+i+it*(data->flap).nface)=*((data->Gammas)+i+m);
    }
    /* Assign aileron vortex strengths */
    m+=(data->flap).nface;
    for (i=0;i<(data->aileron).nface;i++){
        *((data->aileron).gamma+i+it*(data->aileron).nface)=*((data->Gammas)+i+m);
    }
    m+=(data->aileron).nface;
    /* Assign horizontal tail vortex strengths */
    for (i=0;i<(data->htail).nface;i++){
        *((data->htail).gamma+i+it*(data->htail).nface)=*((data->Gammas)+i+m);
    }  
    m+=(data->htail).nface;
    /* Assign elevator vortex strengths */
    for (i=0;i<(data->elevator).nface;i++){
        *((data->elevator).gamma+i+it*(data->elevator).nface)=*((data->Gammas)+i+m);
    }
    m+=(data->elevator).nface;
    /* Assign vertical tail vortex strengths */
    for (i=0;i<(data->vtail).nface;i++){
        *((data->vtail).gamma+i+it*(data->vtail).nface)=*((data->Gammas)+i+m);
    }
    m+=(data->vtail).nface;
    /* Assign rudder vortex strengths */
    for (i=0;i<(data->rudder).nface;i++){
        *((data->rudder).gamma+i+it*(data->rudder).nface)=*((data->Gammas)+i+m);
    }
}

void assignwind(struct VLMData *data)
{
    /* Assign the vortex strengths in vector Gammas back onto the panels of all the lifting surfaces */
    int i,m;
        
    /* Assign wing vortex strengths */
    m=0;
    for (i=0;i<(data->wing).nface;i++){
        *((data->wing).wind+i)=*((data->wind)+i);
    }
    /* Assign flap vortex strengths */
    m+=(data->wing).nface;
    for (i=0;i<(data->flap).nface;i++){
        *((data->flap).wind+i)=*((data->wind)+i+m);
    }
    /* Assign aileron vortex strengths */
    m+=(data->flap).nface;
    for (i=0;i<(data->aileron).nface;i++){
        *((data->aileron).wind+i)=*((data->wind)+i+m);
    }
    m+=(data->aileron).nface;
    /* Assign horizontal tail vortex strengths */
    for (i=0;i<(data->htail).nface;i++){
        *((data->htail).wind+i)=*((data->wind)+i+m);
    }  
    m+=(data->htail).nface;
    /* Assign elevator vortex strengths */
    for (i=0;i<(data->elevator).nface;i++){
        *((data->elevator).wind+i)=*((data->wind)+i+m);
    }
    m+=(data->elevator).nface;
    /* Assign vertical tail vortex strengths */
    for (i=0;i<(data->vtail).nface;i++){
        *((data->vtail).wind+i)=*((data->wind)+i+m);
    }  
    m+=(data->vtail).nface;
    /* Assign rudder vortex strengths */
    for (i=0;i<(data->rudder).nface;i++){
        *((data->rudder).wind+i)=*((data->wind)+i+m);
    }
}