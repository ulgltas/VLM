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
//  vVertices.c
//
//  Vertex treatment functions
//
//

#include "vVertices.h"
#include "vLiftsurf.h"
#include <stdlib.h>

void assignvertices(struct liftsurf *pflap, double *xflap, double *yflap, double *zflap, double *xvflap, double *yvflap, double *zvflap)
{
    /* Store the x, y and z coordinates of the current lifting surface in the vertices field
     * of the corresponding lifting surface array */
    int i,nvertd2;
    
    nvertd2=pflap->nvert/2;
    pflap->vertices= (double *)malloc(sizeof(double)*pflap->nvert*3);
    pflap->vortex= (double *)malloc(sizeof(double)*pflap->nvert*3);
    for (i=0;i<nvertd2;i++){
        /* Store in pflap->vertices */
        /* Right lifting surface */
        *(pflap->vertices+i)=*(xflap+i);
        *(pflap->vertices+i+pflap->nvert)=*(yflap+i);
        *(pflap->vertices+i+2*pflap->nvert)=*(zflap+i);
        *(pflap->vortex+i)=*(xvflap+i);
        *(pflap->vortex+i+pflap->nvert)=*(yvflap+i);
        *(pflap->vortex+i+2*pflap->nvert)=*(zvflap+i);
        /* Left lifting surface */
        *(pflap->vertices+i+nvertd2)=*(xflap+i);
        *(pflap->vertices+i+nvertd2+pflap->nvert)=-*(yflap+i);
        *(pflap->vertices+i+nvertd2+2*pflap->nvert)=*(zflap+i);
        *(pflap->vortex+i+nvertd2)=*(xvflap+i);
        *(pflap->vortex+i+nvertd2+pflap->nvert)=-*(yvflap+i);
        *(pflap->vortex+i+nvertd2+2*pflap->nvert)=*(zvflap+i);
    }
}

void assignvertices_fin(struct liftsurf *pflap, double *xflap, double *yflap, double *zflap, double *xvflap, double *yvflap, double *zvflap, int OptVTFusMounted, int OptVTTailMounted, int OptVTWingMounted)
{
    /* Store the x, y and z coordinates of the current lifting surface in the vertices field
     * of the corresponding lifting surface array */
    int i,nvertd2;
    
    if (OptVTFusMounted == 1){
        nvertd2=pflap->nvert/2;
        pflap->vertices= (double *)malloc(sizeof(double)*pflap->nvert*3);
        pflap->vortex= (double *)malloc(sizeof(double)*pflap->nvert*3);
        for (i=0;i<nvertd2;i++){
            /* Store in pflap->vertices */
            /* Real fin */
            *(pflap->vertices+i)=*(xflap+i);
            *(pflap->vertices+i+pflap->nvert)=*(yflap+i);
            *(pflap->vertices+i+2*pflap->nvert)=*(zflap+i);
            *(pflap->vortex+i)=*(xvflap+i);
            *(pflap->vortex+i+pflap->nvert)=*(yvflap+i);
            *(pflap->vortex+i+2*pflap->nvert)=*(zvflap+i);
            /* Mirror image */
            *(pflap->vertices+i+nvertd2)=*(xflap+i);
            *(pflap->vertices+i+nvertd2+pflap->nvert)=*(yflap+i);
            *(pflap->vertices+i+nvertd2+2*pflap->nvert)=-*(zflap+i);
            *(pflap->vortex+i+nvertd2)=*(xvflap+i);
            *(pflap->vortex+i+nvertd2+pflap->nvert)=*(yvflap+i);
            *(pflap->vortex+i+nvertd2+2*pflap->nvert)=-*(zvflap+i);
        }
    }
}