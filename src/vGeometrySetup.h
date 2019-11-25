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
//  vGeometrySetup.h
//
//  Geometry setup functions
//
//
#ifndef vGeometrySetup_h
#define vGeometrySetup_h
#include "vVLMData.h"
void geometry_setup(struct VLMData *data);
void reset_wake(struct VLMData *data);
void reset_geometry(struct VLMData *data);
#endif /* vGeometrySetup_h */