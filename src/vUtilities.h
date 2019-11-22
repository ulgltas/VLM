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
//  vUtilities.h
//  
//
//  General things used throughout the code
//
//

#ifndef vUtilities_h
#define vUtilities_h

#include <stdio.h>

int compare_function(const void *a,const void *b);

int checkdoubles(double *ypos, int m);

int finddoubles(double *ypos, int m);

int findindex(double *ypline, int np1, double value);

#endif /* vUtilities_h */
