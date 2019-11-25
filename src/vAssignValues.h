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
//  vAssignValues.h
//
//  Assign values back to the Liftsurf struct
//
//
#ifndef vAssignValues_h
#define vAssignValues_h

#include "vLiftsurf.h"
#include "vVLMData.h"

void assignGammas(struct VLMData *data, int it);

void assignwind(struct VLMData *data);

#endif /* vAssignValues_h */