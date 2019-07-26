//
//  vGeometrySetup.c
//
//  Geometry setup function
//
//

#include "vLiftsurf.h"
#include "vVLMData.h"
#include "vGeometrySetup.h"
#include "vGeometry.h"
#include "vPanel.h"
#include "vVortex.h"
#include "vWake.h"
#include <stdio.h>
#include <stdlib.h>

void geometry_setup(struct VLMData *data)
{
    /* Calculate collocation points and vortex segment lengths */
    colvec(&(data->flap));
    colvec(&(data->aileron));
    colvec(&(data->wing));
    colvec(&(data->htail));
    colvec(&(data->elevator));
    colvec(&(data->vtail));
    colvec(&(data->rudder));
    /* Calculate normal vectors and surfaces */
    normals(&(data->flap));
    normals(&(data->aileron));
    normals(&(data->wing));
    normals(&(data->htail));
    normals(&(data->elevator));
    normals(&(data->vtail));
    normals(&(data->rudder));
    /* Calculate tangential vectors */
    tangentials(&(data->flap));
    tangentials(&(data->aileron));
    tangentials(&(data->wing));
    tangentials(&(data->htail));
    tangentials(&(data->elevator));
    tangentials(&(data->vtail));
    tangentials(&(data->rudder));
}
}