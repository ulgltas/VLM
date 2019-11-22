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
//  vRHS.c
//
//  Calculate the values in the Right Hand Side vector
//
//

#include "vRHS.h"
#include "vLiftsurf.h"
#include "vVLMData.h"

void calcRHS(struct VLMData *data)
{   
    /* Calculate the Right Hand Side vector */
    /* Cycle through the lifting surfaces in the order wing, flap, aileron, horizontal tail, elevator */
    
    int i,m;

    m=0;
    for (i=0; i < (data->wing).nface; i++){
        *((data->RHS)+i)=-(*((data->wing).uvw+i)+data->UVW[0])* *((data->wing).normal+i)-(*((data->wing).uvw+i+(data->wing).nface)+data->UVW[1])* *((data->wing).normal+i+(data->wing).nface) -
                (*((data->wing).uvw+i+2*(data->wing).nface)+data->UVW[2])* *((data->wing).normal+i+2*(data->wing).nface);
    }
    m+=(data->wing).nface;
    for (i=0; i < (data->flap).nface; i++){
        *((data->RHS)+i+m)=-(*((data->flap).uvw+i)+data->UVW[0])* *((data->flap).normal+i)-(*((data->flap).uvw+i+(data->flap).nface)+data->UVW[1])* *((data->flap).normal+i+(data->flap).nface) -
                (*((data->flap).uvw+i+2*(data->flap).nface)+data->UVW[2])* *((data->flap).normal+i+2*(data->flap).nface);
    }
    m+=(data->flap).nface;
    for (i=0; i < (data->aileron).nface; i++){
        *((data->RHS)+i+m)=-(*((data->aileron).uvw+i)+data->UVW[0])* *((data->aileron).normal+i)-(*((data->aileron).uvw+i+(data->aileron).nface)+data->UVW[1])* *((data->aileron).normal+i+(data->aileron).nface) -
                (*((data->aileron).uvw+i+2*(data->aileron).nface)+data->UVW[2])* *((data->aileron).normal+i+2*(data->aileron).nface);
    }
    m+=(data->aileron).nface;
    for (i=0; i < (data->htail).nface; i++){
        *((data->RHS)+i+m)=-(*((data->htail).uvw+i)+data->UVW[0])* *((data->htail).normal+i)-(*((data->htail).uvw+i+(data->htail).nface)+data->UVW[1])* *((data->htail).normal+i+(data->htail).nface) -
                (*((data->htail).uvw+i+2*(data->htail).nface)+data->UVW[2])* *((data->htail).normal+i+2*(data->htail).nface);
    }
    m+=(data->htail).nface;
    for (i=0; i < (data->elevator).nface; i++){
        *((data->RHS)+i+m)=-(*((data->elevator).uvw+i)+data->UVW[0])* *((data->elevator).normal+i)-(*((data->elevator).uvw+i+(data->elevator).nface)+data->UVW[1])* *((data->elevator).normal+i+(data->elevator).nface) -
                (*((data->elevator).uvw+i+2*(data->elevator).nface)+data->UVW[2])* *((data->elevator).normal+i+2*(data->elevator).nface);
    }
    m+=(data->elevator).nface;
    for (i=0; i < (data->vtail).nface; i++){
        *((data->RHS)+i+m)=-(*((data->vtail).uvw+i)+data->UVW[0])* *((data->vtail).normal+i)-(*((data->vtail).uvw+i+(data->vtail).nface)+data->UVW[1])* *((data->vtail).normal+i+(data->vtail).nface) -
                (*((data->vtail).uvw+i+2*(data->vtail).nface)+data->UVW[2])* *((data->vtail).normal+i+2*(data->vtail).nface);
    }
    m+=(data->vtail).nface;
    for (i=0; i < (data->rudder).nface; i++){
        *((data->RHS)+i+m)=-(*((data->rudder).uvw+i)+data->UVW[0])* *((data->rudder).normal+i)-(*((data->rudder).uvw+i+(data->rudder).nface)+data->UVW[1])* *((data->rudder).normal+i+(data->rudder).nface) -
                (*((data->rudder).uvw+i+2*(data->rudder).nface)+data->UVW[2])* *((data->rudder).normal+i+2*(data->rudder).nface);
    }
}