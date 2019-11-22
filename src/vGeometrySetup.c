/**
 * Copyright 2019 Université de Liège
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
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

void reset_wake(struct VLMData *data)
{
    /* Reset wake to initial values */
    for (int i=0;i<((data->wing).nshed+(data->wing).nwakes)*data->ntimes;i++)
    {
        *((data->wing).xw+i)=0.0;
        *((data->wing).yw+i)=0.0;
        *((data->wing).zw+i)=0.0;
    }
    /* Set contents of pwing->gw to zero */
    for (int i=0;i<(data->wing).nshed*data->ntimes;i++)
    {
        *((data->wing).gw+i)=0.0;
    }
    shedwake(&data->wing);
}

void reset_geometry(struct VLMData *data)
{
    // De-allocate normal vectors
    free((data->flap).normal);
    free((data->aileron).normal);
    free((data->wing).normal);
    free((data->htail).normal);
    free((data->elevator).normal);
    free((data->vtail).normal);
    free((data->rudder).normal);
    free((data->flap).nsurf);
    free((data->aileron).nsurf);
    free((data->wing).nsurf);
    free((data->htail).nsurf);
    free((data->elevator).nsurf);
    free((data->vtail).nsurf);
    free((data->rudder).nsurf);
    // De-allocate tangential vectors
    free((data->flap).tangx);
    free((data->aileron).tangx);
    free((data->wing).tangx);
    free((data->htail).tangx);
    free((data->elevator).tangx);
    free((data->vtail).tangx);
    free((data->rudder).tangx);
    free((data->flap).tangy);
    free((data->aileron).tangy);
    free((data->wing).tangy);
    free((data->htail).tangy);
    free((data->elevator).tangy);
    free((data->vtail).tangy);
    free((data->rudder).tangy);
    // De-allocate vector collocation points
    free((data->flap).control);
    free((data->aileron).control);
    free((data->wing).control);
    free((data->htail).control);
    free((data->elevator).control);
    free((data->vtail).control);
    free((data->rudder).control);
    free((data->flap).dxy);
    free((data->aileron).dxy);
    free((data->wing).dxy);
    free((data->htail).dxy);
    free((data->elevator).dxy);
    free((data->vtail).dxy);
    free((data->rudder).dxy);
}