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
//  vVortex.h
//
//  Vortex related functions
//
//

#ifndef vVortex_h
#define vVortex_h

#include "vLiftsurf.h"

void colvec(struct liftsurf *pflap);

void vortexblob(double uvw[], double x, double y, double z, double x1, double y1, double z1,
        double x2, double y2, double z2, double gama);

void vortex(double uvw[], double x, double y, double z, double x1, double y1, double z1,
        double x2, double y2, double z2, double gama);

#endif /* vVortex_h */