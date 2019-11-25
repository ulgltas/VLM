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