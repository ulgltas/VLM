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
//  vGeometry.c
//
//  Geometrical functions for the panels
//
//

#include "vGeometry.h"
#include "vLiftsurf.h"
#include <stdlib.h>
#include <math.h>

void cross(double xcrossy[], double x[], double y[])
{
    /* Calculate the cross produce of two vectors x and y */
    xcrossy[0]=x[1]*y[2]-x[2]*y[1];
    xcrossy[1]=x[2]*y[0]-x[0]*y[2];
    xcrossy[2]=x[0]*y[1]-x[1]*y[0];
}

void normals(struct liftsurf *pflap)
{
    /* Calculate the unit normal vector and normal surface to the surface panels xp,yp,zp */
    int i,j,k,l,lp1;
    double x[4],y[4],z[4],a;
    double r1[3],r2[3],r1r2[3],csum[3];
    
    pflap->normal=(double *)malloc(sizeof(double)*pflap->nface*3);
    pflap->nsurf=(double *)malloc(sizeof(double)*pflap->nface);
    /* Cycle through all the panels */
    for (i=0;i<pflap->nface;i++){
        /* Select panel endpoints */
        x[0]=*(pflap->vertices+*(pflap->faces+i));
        y[0]=*(pflap->vertices+pflap->nvert+*(pflap->faces+i));
        z[0]=*(pflap->vertices+2*pflap->nvert+*(pflap->faces+i));
        x[1]=*(pflap->vertices+*(pflap->faces+i+pflap->nface));
        y[1]=*(pflap->vertices+pflap->nvert+*(pflap->faces+i+pflap->nface));
        z[1]=*(pflap->vertices+2*pflap->nvert+*(pflap->faces+i+pflap->nface));
        x[2]=*(pflap->vertices+*(pflap->faces+i+2*pflap->nface));
        y[2]=*(pflap->vertices+pflap->nvert+*(pflap->faces+i+2*pflap->nface));
        z[2]=*(pflap->vertices+2*pflap->nvert+*(pflap->faces+i+2*pflap->nface));
        x[3]=*(pflap->vertices+*(pflap->faces+i+3*pflap->nface));
        y[3]=*(pflap->vertices+pflap->nvert+*(pflap->faces+i+3*pflap->nface));
        z[3]=*(pflap->vertices+2*pflap->nvert+*(pflap->faces+i+3*pflap->nface));
        
        /* Calculate distances between endpoints */
        r1[0]=x[2]-x[0];
        r1[1]=y[2]-y[0];
        r1[2]=z[2]-z[0];
        r2[0]=x[1]-x[3];
        r2[1]=y[1]-y[3];
        r2[2]=z[1]-z[3];
        
        /* Calculate normal vector */
        cross(r1r2,r1,r2);
        a=sqrt(r1r2[0]*r1r2[0]+r1r2[1]*r1r2[1]+r1r2[2]*r1r2[2]);
        
        /* Assign normal vector */
        if (i <pflap->nface/2){
            *(pflap->normal+i)=-r1r2[0]/a;
            *(pflap->normal+i+pflap->nface)=-r1r2[1]/a;
            *(pflap->normal+i+2*pflap->nface)=-r1r2[2]/a;
        }else{
            *(pflap->normal+i)=r1r2[0]/a;
            *(pflap->normal+i+pflap->nface)=r1r2[1]/a;
            *(pflap->normal+i+2*pflap->nface)=r1r2[2]/a;
        }
        
        /* Calculate normal area */
        for (k=0;k<3;k++){
            csum[k]=0.0;
        }
        for (l=0;l<4;l++){
            if (l<3){
                lp1=l+1;
            }else{
                lp1=0;
            }
            r1[0]=x[l];
            r1[1]=y[l];
            r1[2]=z[l];
            r2[0]=x[lp1];
            r2[1]=y[lp1];
            r2[2]=z[lp1];
            cross(r1r2,r1,r2);
            for (k=0;k<3;k++){
                csum[k]=csum[k]+r1r2[k];
            }
        }
        *(pflap->nsurf+i)=0.5*fabs(*(pflap->normal+i)*csum[0]+*(pflap->normal+i+pflap->nface)*csum[1]+*(pflap->normal+i+2*pflap->nface)*csum[2]);
    }
}

void tangentials(struct liftsurf *pflap)
{
    /* Calculate the two unit vectors tangent to surface panels xp,yp,zp */
    int i;
    double tx1,tx2,tx3,txall;
    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;

    pflap->tangx=(double *)malloc(sizeof(double)*pflap->nface*3);
    pflap->tangy=(double *)malloc(sizeof(double)*pflap->nface*3);
    /* Cycle through all the panels */
    for (i=0;i<pflap->nface;i++){
        /* Vertices of each geometric panel */
        x1=*(pflap->vertices+*(pflap->faces+i));
        x2=*(pflap->vertices+*(pflap->faces+i+pflap->nface));
        x3=*(pflap->vertices+*(pflap->faces+i+2*pflap->nface));
        x4=*(pflap->vertices+*(pflap->faces+i+3*pflap->nface));
        y1=*(pflap->vertices+pflap->nvert+*(pflap->faces+i));
        y2=*(pflap->vertices+pflap->nvert+*(pflap->faces+i+pflap->nface));
        y3=*(pflap->vertices+pflap->nvert+*(pflap->faces+i+2*pflap->nface));
        y4=*(pflap->vertices+pflap->nvert+*(pflap->faces+i+3*pflap->nface));
        z1=*(pflap->vertices+2*pflap->nvert+*(pflap->faces+i));
        z2=*(pflap->vertices+2*pflap->nvert+*(pflap->faces+i+pflap->nface));
        z3=*(pflap->vertices+2*pflap->nvert+*(pflap->faces+i+2*pflap->nface));
        z4=*(pflap->vertices+2*pflap->nvert+*(pflap->faces+i+3*pflap->nface));
        /* Tangential vector in the x-direction */
        tx1=((x2-x1)+(x3-x4))/2.0;
        tx2=((y2-y1)+(y3-y4))/2.0;
        tx3=((z2-z1)+(z3-z4))/2.0;
        /* Normalize it */
        txall=sqrt(tx1*tx1+tx2*tx2+tx3*tx3);
        *(pflap->tangx+i)=tx1/txall;
        *(pflap->tangx+i+pflap->nface)=tx2/txall;
        *(pflap->tangx+i+2*pflap->nface)=tx3/txall;
        /* Tangential vector in the y-direction */
        tx1=0.0;
        tx2=3.0*(y4-y1)/4.0+(y3-y2)/4.0;
        tx3=3.0*(z4-z1)/4.0+(z3-z2)/4.0;
        /* Normalize it */
        txall=sqrt(tx1*tx1+tx2*tx2+tx3*tx3);
        /* All the y-tangential vectors face in the positive y direction */
        if (i <pflap->nface/2){
            *(pflap->tangy+i)=tx1/txall;
            *(pflap->tangy+i+pflap->nface)=tx2/txall;
            *(pflap->tangy+i+2*pflap->nface)=tx3/txall;
        }else{
            *(pflap->tangy+i)=-tx1/txall;
            *(pflap->tangy+i+pflap->nface)=-tx2/txall;
            *(pflap->tangy+i+2*pflap->nface)=-tx3/txall;
        }
    }
}