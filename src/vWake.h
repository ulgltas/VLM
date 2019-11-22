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
//  vWake.h
//
//  Wake related functions
//
//

#ifndef vWake_h
#define vWake_h

#include "vLiftsurf.h"
#include "vVLMData.h"

void storewakeinds(struct liftsurf *pwing);

int findwakeneighbours(struct liftsurf *pwing, int ipanel, int lrneighs[],int cwing);

void setupwakes(struct liftsurf *pwing, int cwing);

void createwake(struct liftsurf *pwing, double cwing,int ntimes);

void correspshedwake(struct liftsurf *pwing);

void shedwake(struct liftsurf *pwing);

void wakegamma(struct liftsurf *pwing, int it);

void propwakexyz(struct liftsurf *pwing, double dt, int it, double UVW[]);

void propwakevort(struct liftsurf *pwing, double dt, int it, double UVW[]);

void infonwake(struct liftsurf *pwing, struct liftsurf *pflap, int it, int summode);

void wakeinf(struct liftsurf *pwing, struct liftsurf *pflap, int nw, int it, int summode);

void calcwakeinf(struct VLMData *data, int it);

void addlastwind(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, int it, int narg);

#endif /* vWake_h */