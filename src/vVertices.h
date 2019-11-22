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
//  vVertices.h
//
//  Vertex treatment functions
//
//

#ifndef vVertices_h
#define vVertices_h

#include "vLiftsurf.h"

void assignvertices(struct liftsurf *pflap, double *xflap, double *yflap, double *zflap, double *xvflap, double *yvflap, double *zvflap);

void assignvertices_fin(struct liftsurf *pflap, double *xflap, double *yflap, double *zflap, double *xvflap, double *yvflap, double *zvflap, int OptVTFusMounted, int OptVTTailMounted, int OptVTWingMounted);

#endif /* vVertices_h */