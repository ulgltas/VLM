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
//  vInput.h
//
//  Data Input
//
//
#ifndef vInput_h
#define vInput_h
#include <stdio.h> 
void importInputFile(char *Infile, double *UVW, double *rho, double *aoa, double *yaw, int *m, int *mht, int *mvt, int *n, int *nht, int *nvt,
                    int *ntimes, double *timestep_denom, int *freewake, double *delta, double *beta, double *eta, double *zeta,
                    char *Wngfile, char *HTailfile, char *VTailfile);

#endif /* vInput_h */