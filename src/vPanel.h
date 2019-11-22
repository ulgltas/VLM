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
//  vPanel.h
//  
//
//  Functions related to panel creation
//
//

#ifndef vPanel_h
#define vPanel_h

#include <stdio.h>
#include "vLiftsurf.h"

int countfaces(int *ijflap, int nflap, int mp1, int np1);

void arrangefaces(int *ijflap, int mp1, int np1, struct liftsurf *pflap);

void vortexpanel(double *xv, double *yv, double *zv, double *xp, double *yp, double *zp, double dxw, int m, int n);

#endif /* vPanel_h */
