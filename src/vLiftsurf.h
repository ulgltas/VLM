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
//  vLiftsurf.h
//  
//
//  Liftsurf struct definition
//
//

#ifndef vLiftsurf_h
#define vLiftsurf_h

#include <stdio.h>
struct liftsurf
{
    int    nface; /* number of panels */
    int    nvert; /* number of vertices */
    int *faces; /* panel definitions */
    int *neighbours; /* neighbouring panels to each panel arranged in the order left, right, upstream, downstream */
    double *vertices; /* panel vertices */
    int nshed; /* Number of trailing edge panels that shed wake */
    int *shedding; /* Trailing edge panel that sheds wake (0 or 1) */
    double *control; /* control point positions */
    double *vortex; /* Vortex panel vertices */
    double *tangx; /* Tangential vectors in x direction */
    double *tangy; /* Tangential vectors in y direction */
    double *normal; /* Normal vectors */
    double *nsurf; /* Normal surfaces of panels */
    double *dxy; /* Vortex segment lengths in x and y directions */
    double *Deltap; /* Pressure difference across each panel */
    double *Deltad; /* Induced drag contribution of each panel */
    double *gamma; /* Vortex strength on each panel */
    double *wind; /* Induced airspeed on each panel */
    double *uvw; /* Local airspeed induced by wakes */
    double *aeroforce; /* Aerodynamic force components in x, y and z directions */
    int nwakes; /* Number of contiguous wakes shed from the lifting surface */
    int *wakeinds; /* Indices of panels that shed wake */
    int *wakelengths; /* Number of wake panels shed at each time step in each wake */
    int *wakecorrespx; /* Correspondence between bound vortex panels and wake panel vertices */
    int *wakecorrespy; /* Correspondence between bound vortex panels and wake panel vertices */
    int *wakecorrespz; /* Correspondence between bound vortex panels and wake panel vertices */
    int *wakecorrespwake; /* Correspondence between wake panel vertices in different wakes */
    double *xw; /* x-coordinates of wake vertices */
    double *yw; /* y-coordinates of wake vertices */
    double *zw; /* z-coordinates of wake vertices */
    double *uw; /* velocity in x-direction of wake vertices */
    double *vw; /* velocity in y-direction of wake vertices */
    double *ww; /* velocity in z-direction of wake vertices */
    double *gw; /* Vortex strength on each wake panel */
    double dxw; // Distance to collocation point
};

#endif /* vLiftsurf_h */
