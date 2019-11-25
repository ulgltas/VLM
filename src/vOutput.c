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
//  vOutput.c
//
//  Data output
// mtn: total number of panels
//

#include "vOutput.h"
#include "vLiftsurf.h"
#include "vVLMData.h"
#include <stdio.h>
#include <stdlib.h>

void exportTextOutput(char *Outfile, int it, struct VLMData data)
{
    int i, j;
    struct liftsurf *pwing = &data.wing;
    struct liftsurf *pflap = &data.flap;
    struct liftsurf *paileron = &data.aileron;
    struct liftsurf *phtail = &data.htail;
    struct liftsurf *pelevator = &data.elevator;
    struct liftsurf *pvtail = &data.vtail;
    struct liftsurf *prudder = &data.rudder;
    FILE *ofp = fopen(Outfile, "w");
    if (ofp == NULL) {
        fprintf(stderr, "Can't open output file %s!\n",Outfile);
        exit(1);
    }

    fprintf(ofp,"it=%i;\n",it);
    fprintf(ofp,"dt=%f;\n",data.dt);
    
    fprintf(ofp,"totalforce=[");
    for (i=0;i<4;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(data.totalforce+i+j*4));
            }
            else{
                fprintf(ofp,"%f   ",*(data.totalforce+i+j*4));
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"flap.faces=[");
    for (i=0;i<data.flap.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pflap->faces+i+j*data.flap.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pflap->faces+i+j*data.flap.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"aileron.faces=[");
    for (i=0;i<data.aileron.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(paileron->faces+i+j*data.aileron.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(paileron->faces+i+j*data.aileron.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"wing.faces=[");
    for (i=0;i<data.wing.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pwing->faces+i+j*data.wing.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pwing->faces+i+j*data.wing.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"htail.faces=[");
    for (i=0;i<data.htail.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(phtail->faces+i+j*data.htail.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(phtail->faces+i+j*data.htail.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"elevator.faces=[");
    for (i=0;i<data.elevator.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pelevator->faces+i+j*data.elevator.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pelevator->faces+i+j*data.elevator.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"vtail.faces=[");
    for (i=0;i<data.vtail.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pvtail->faces+i+j*data.vtail.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pvtail->faces+i+j*data.vtail.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"rudder.faces=[");
    for (i=0;i<data.rudder.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(prudder->faces+i+j*data.rudder.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(prudder->faces+i+j*data.rudder.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");


    fprintf(ofp,"flap.vertices=[");
    for (j=0;j<data.flap.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->vertices+j),*(pflap->vertices+j+data.flap.nvert),*(pflap->vertices+j+2*data.flap.nvert));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"aileron.vertices=[");
    for (j=0;j<data.aileron.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->vertices+j),*(paileron->vertices+j+data.aileron.nvert),*(paileron->vertices+j+2*data.aileron.nvert));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"wing.vertices=[");
    for (j=0;j<data.wing.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->vertices+j),*(pwing->vertices+j+data.wing.nvert),*(pwing->vertices+j+2*data.wing.nvert));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"htail.vertices=[");
    for (j=0;j<data.htail.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->vertices+j),*(phtail->vertices+j+data.htail.nvert),*(phtail->vertices+j+2*data.htail.nvert));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"elevator.vertices=[");
    for (j=0;j<data.elevator.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->vertices+j),*(pelevator->vertices+j+data.elevator.nvert),*(pelevator->vertices+j+2*data.elevator.nvert));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"vtail.vertices=[");
    for (j=0;j<data.vtail.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->vertices+j),*(pvtail->vertices+j+data.vtail.nvert),*(pvtail->vertices+j+2*data.vtail.nvert));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"rudder.vertices=[");
    for (j=0;j<data.rudder.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->vertices+j),*(prudder->vertices+j+data.rudder.nvert),*(prudder->vertices+j+2*data.rudder.nvert));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"pwing.uvw=[");
    for (j=0;j<data.wing.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->uvw+j),*(pwing->uvw+j+data.wing.nface),*(pwing->uvw+j+2*data.wing.nface));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"pflap.uvw=[");
    for (j=0;j<data.flap.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->uvw+j),*(pflap->uvw+j+data.flap.nface),*(pflap->uvw+j+2*data.flap.nface));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"paileron.uvw=[");
    for (j=0;j<data.aileron.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->uvw+j),*(paileron->uvw+j+data.aileron.nface),*(paileron->uvw+j+2*data.aileron.nface));
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"phtail.uvw=[");
    for (j=0;j<data.htail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->uvw+j),*(phtail->uvw+j+data.htail.nface),*(phtail->uvw+j+2*data.htail.nface));
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"pelevator.uvw=[");
    for (j=0;j<data.elevator.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->uvw+j),*(pelevator->uvw+j+data.elevator.nface),*(pelevator->uvw+j+2*data.elevator.nface));
    }
    fprintf(ofp,"];\n");     

    fprintf(ofp,"pvtail.uvw=[");
    for (j=0;j<data.vtail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->uvw+j),*(pvtail->uvw+j+data.vtail.nface),*(pvtail->uvw+j+2*data.vtail.nface));
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"prudder.uvw=[");
    for (j=0;j<data.rudder.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->uvw+j),*(prudder->uvw+j+data.rudder.nface),*(prudder->uvw+j+2*data.rudder.nface));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"flap_shedding=[");
    for (j=0;j<data.flap.nface;j++){
        fprintf(ofp,"%i\n",*(pflap->shedding+j));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"aileron_shedding=[");
    for (j=0;j<data.aileron.nface;j++){
        fprintf(ofp,"%i\n",*(paileron->shedding+j));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"wing_shedding=[");
    for (j=0;j<data.wing.nface;j++){
        fprintf(ofp,"%i\n",*(pwing->shedding+j));
    }    
    fprintf(ofp,"];\n");

    fprintf(ofp,"flapcp=[");
    for (j=0;j<data.flap.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->control+j),*(pflap->control+j+data.flap.nface),*(pflap->control+j+2*data.flap.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"flapcv=[");
    for (j=0;j<data.flap.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->vortex+j),*(pflap->vortex+j+data.flap.nvert),*(pflap->vortex+j+2*data.flap.nvert));
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"flapnorm=[");
    for (j=0;j<data.flap.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->normal+j),*(pflap->normal+j+data.flap.nface),*(pflap->normal+j+2*data.flap.nface));
    }
    fprintf(ofp,"];\n");     

    fprintf(ofp,"flaptangx=[");
    for (j=0;j<data.flap.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->tangx+j),*(pflap->tangx+j+data.flap.nface),*(pflap->tangx+j+2*data.flap.nface));
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"flaptangy=[");
    for (j=0;j<data.flap.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->tangy+j),*(pflap->tangy+j+data.flap.nface),*(pflap->tangy+j+2*data.flap.nface));
    }
    fprintf(ofp,"];\n");      
    
    fprintf(ofp,"flapdxy=[");
    for (j=0;j<data.flap.nface;j++){
        fprintf(ofp,"%f %f\n",*(pflap->dxy+j),*(pflap->dxy+j+data.flap.nface));
    }
    fprintf(ofp,"];\n");   
    
    fprintf(ofp,"flapnsurf=[");
    for (j=0;j<data.flap.nface;j++){
        fprintf(ofp,"%f\n",*(pflap->nsurf+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"flapDeltap=[");
    for (j=0;j<data.flap.nface;j++){
        fprintf(ofp,"%f\n",*(pflap->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"flapDeltad=[");
    for (j=0;j<data.flap.nface;j++){
        fprintf(ofp,"%f\n",*(pflap->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"flapwind=[");
    for (j=0;j<data.flap.nface;j++){
        fprintf(ofp,"%f\n",*(pflap->wind+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"flapgamma=[");
    for (i=0;i<data.flap.nface;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pflap->gamma+i+j*data.flap.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(pflap->gamma+i+j*data.flap.nface));
            }
        }
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"pflap.xw=[");
    for (i=0;i<(pflap->nshed+pflap->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pflap->xw+i+j*(pflap->nshed+pflap->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pflap->xw+i+j*(pflap->nshed+pflap->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");        

    fprintf(ofp,"pflap.yw=[");
    for (i=0;i<(pflap->nshed+pflap->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pflap->yw+i+j*(pflap->nshed+pflap->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pflap->yw+i+j*(pflap->nshed+pflap->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"pflap.zw=[");
    for (i=0;i<(pflap->nshed+pflap->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pflap->zw+i+j*(pflap->nshed+pflap->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pflap->zw+i+j*(pflap->nshed+pflap->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");


    fprintf(ofp,"pflap.uw=[");
    for (i=0;i<(pflap->nshed+pflap->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pflap->uw+i+j*(pflap->nshed+pflap->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pflap->uw+i+j*(pflap->nshed+pflap->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");        

    fprintf(ofp,"pflap.vw=[");
    for (i=0;i<(pflap->nshed+pflap->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pflap->vw+i+j*(pflap->nshed+pflap->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pflap->vw+i+j*(pflap->nshed+pflap->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"pflap.ww=[");
    for (i=0;i<(pflap->nshed+pflap->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pflap->ww+i+j*(pflap->nshed+pflap->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pflap->ww+i+j*(pflap->nshed+pflap->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");

    
    fprintf(ofp,"pflap.gw=[");
    for (i=0;i<pflap->nshed;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pflap->gw+i+j*pflap->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(pflap->gw+i+j*pflap->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"flapneighbours=[");
    for (i=0;i<data.flap.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pflap->neighbours+i+j*data.flap.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pflap->neighbours+i+j*data.flap.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"pflap.nwakes=%i ;",pflap->nwakes);
    fprintf(ofp,"pflap.nshed=%i ;",pflap->nshed);
    fprintf(ofp,"pflap.wakeinds=[");
    for (i=0;i<data.flap.nshed;i++){
        fprintf(ofp,"%i\n",*(pflap->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"pflap.wakelengths=[");
    for (i=0;i<data.flap.nwakes;i++){
        fprintf(ofp,"%i\n",*(pflap->wakelengths+i));
    }
    fprintf(ofp,"];\n");       
    
    fprintf(ofp,"aileroncp=[");
    for (j=0;j<data.aileron.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->control+j),*(paileron->control+j+data.aileron.nface),*(paileron->control+j+2*data.aileron.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"aileroncv=[");
    for (j=0;j<data.aileron.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->vortex+j),*(paileron->vortex+j+data.aileron.nvert),*(paileron->vortex+j+2*data.aileron.nvert));
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"aileronnorm=[");
    for (j=0;j<data.aileron.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->normal+j),*(paileron->normal+j+data.aileron.nface),*(paileron->normal+j+2*data.aileron.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"ailerontangx=[");
    for (j=0;j<data.aileron.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->tangx+j),*(paileron->tangx+j+data.aileron.nface),*(paileron->tangx+j+2*data.aileron.nface));
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"ailerontangy=[");
    for (j=0;j<data.aileron.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->tangy+j),*(paileron->tangy+j+data.aileron.nface),*(paileron->tangy+j+2*data.aileron.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"ailerondxy=[");
    for (j=0;j<data.aileron.nface;j++){
        fprintf(ofp,"%f %f\n",*(paileron->dxy+j),*(paileron->dxy+j+data.aileron.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"aileronnsurf=[");
    for (j=0;j<data.aileron.nface;j++){
        fprintf(ofp,"%f\n",*(paileron->nsurf+j));
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"aileronDeltap=[");
    for (j=0;j<data.aileron.nface;j++){
        fprintf(ofp,"%f\n",*(paileron->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"aileronDeltad=[");
    for (j=0;j<data.aileron.nface;j++){
        fprintf(ofp,"%f\n",*(paileron->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"aileronwind=[");
    for (j=0;j<data.aileron.nface;j++){
        fprintf(ofp,"%f\n",*(paileron->wind+j));
    }
    fprintf(ofp,"];\n");   
    
    fprintf(ofp,"ailerongamma=[");
    for (i=0;i<data.aileron.nface;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(paileron->gamma+i+j*data.aileron.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(paileron->gamma+i+j*data.aileron.nface));
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"paileron.xw=[");
    for (i=0;i<(paileron->nshed+paileron->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(paileron->xw+i+j*(paileron->nshed+paileron->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(paileron->xw+i+j*(paileron->nshed+paileron->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     

    fprintf(ofp,"paileron.yw=[");
    for (i=0;i<(paileron->nshed+paileron->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(paileron->yw+i+j*(paileron->nshed+paileron->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(paileron->yw+i+j*(paileron->nshed+paileron->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"paileron.zw=[");
    for (i=0;i<(paileron->nshed+paileron->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(paileron->zw+i+j*(paileron->nshed+paileron->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(paileron->zw+i+j*(paileron->nshed+paileron->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     

    fprintf(ofp,"paileron.uw=[");
    for (i=0;i<(paileron->nshed+paileron->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(paileron->uw+i+j*(paileron->nshed+paileron->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(paileron->uw+i+j*(paileron->nshed+paileron->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     

    fprintf(ofp,"paileron.vw=[");
    for (i=0;i<(paileron->nshed+paileron->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(paileron->vw+i+j*(paileron->nshed+paileron->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(paileron->vw+i+j*(paileron->nshed+paileron->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"paileron.ww=[");
    for (i=0;i<(paileron->nshed+paileron->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(paileron->ww+i+j*(paileron->nshed+paileron->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(paileron->ww+i+j*(paileron->nshed+paileron->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    
    fprintf(ofp,"paileron.gw=[");
    for (i=0;i<paileron->nshed;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(paileron->gw+i+j*paileron->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(paileron->gw+i+j*paileron->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"aileronneighbours=[");
    for (i=0;i<data.aileron.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(paileron->neighbours+i+j*data.aileron.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(paileron->neighbours+i+j*data.aileron.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"paileron.nwakes=%i ;",paileron->nwakes);
    fprintf(ofp,"paileron.nshed=%i ;",paileron->nshed);
    fprintf(ofp,"paileron.wakeinds=[");
    for (i=0;i<data.aileron.nshed;i++){
        fprintf(ofp,"%i\n",*(paileron->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"paileron.wakelengths=[");
    for (i=0;i<data.aileron.nwakes;i++){
        fprintf(ofp,"%i\n",*(paileron->wakelengths+i));
    }
    fprintf(ofp,"];\n");    
    
    
    fprintf(ofp,"wingcp=[");
    for (j=0;j<data.wing.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->control+j),*(pwing->control+j+data.wing.nface),*(pwing->control+j+2*data.wing.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"wingcv=[");
    for (j=0;j<data.wing.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->vortex+j),*(pwing->vortex+j+data.wing.nvert),*(pwing->vortex+j+2*data.wing.nvert));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"wingnorm=[");
    for (j=0;j<data.wing.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->normal+j),*(pwing->normal+j+data.wing.nface),*(pwing->normal+j+2*data.wing.nface));
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"wingtangx=[");
    for (j=0;j<data.wing.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->tangx+j),*(pwing->tangx+j+data.wing.nface),*(pwing->tangx+j+2*data.wing.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"wingtangy=[");
    for (j=0;j<data.wing.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->tangy+j),*(pwing->tangy+j+data.wing.nface),*(pwing->tangy+j+2*data.wing.nface));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"wingdxy=[");
    for (j=0;j<data.wing.nface;j++){
        fprintf(ofp,"%f %f\n",*(pwing->dxy+j),*(pwing->dxy+j+data.wing.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"wingnsurf=[");
    for (j=0;j<data.wing.nface;j++){
        fprintf(ofp,"%f\n",*(pwing->nsurf+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"wingDeltap=[");
    for (j=0;j<data.wing.nface;j++){
        fprintf(ofp,"%f\n",*(pwing->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"wingDeltad=[");
    for (j=0;j<data.wing.nface;j++){
        fprintf(ofp,"%f\n",*(pwing->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"wingwind=[");
    for (j=0;j<data.wing.nface;j++){
        fprintf(ofp,"%f\n",*(pwing->wind+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"winggamma=[");
    for (i=0;i<data.wing.nface;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pwing->gamma+i+j*data.wing.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(pwing->gamma+i+j*data.wing.nface));
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"pwing.xw=[");
    for (i=0;i<(pwing->nshed+pwing->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pwing->xw+i+j*(pwing->nshed+pwing->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pwing->xw+i+j*(pwing->nshed+pwing->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"pwing.yw=[");
    for (i=0;i<(pwing->nshed+pwing->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pwing->yw+i+j*(pwing->nshed+pwing->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pwing->yw+i+j*(pwing->nshed+pwing->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"pwing.zw=[");
    for (i=0;i<(pwing->nshed+pwing->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pwing->zw+i+j*(pwing->nshed+pwing->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pwing->zw+i+j*(pwing->nshed+pwing->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"pwing.uw=[");
    for (i=0;i<(pwing->nshed+pwing->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pwing->uw+i+j*(pwing->nshed+pwing->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pwing->uw+i+j*(pwing->nshed+pwing->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"pwing.vw=[");
    for (i=0;i<(pwing->nshed+pwing->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pwing->vw+i+j*(pwing->nshed+pwing->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pwing->vw+i+j*(pwing->nshed+pwing->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"pwing.ww=[");
    for (i=0;i<(pwing->nshed+pwing->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pwing->ww+i+j*(pwing->nshed+pwing->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pwing->ww+i+j*(pwing->nshed+pwing->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"pwing.gw=[");
    for (i=0;i<pwing->nshed;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pwing->gw+i+j*pwing->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(pwing->gw+i+j*pwing->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"wingneighbours=[");
    for (i=0;i<data.wing.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pwing->neighbours+i+j*data.wing.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pwing->neighbours+i+j*data.wing.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"pwing.nwakes=%i ;",pwing->nwakes);
    fprintf(ofp,"pwing.nshed=%i ;",pwing->nshed);
    fprintf(ofp,"pwing.wakeinds=[");
    for (i=0;i<data.wing.nshed;i++){
        fprintf(ofp,"%i\n",*(pwing->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"pwing.wakelengths=[");
    for (i=0;i<data.wing.nwakes;i++){
        fprintf(ofp,"%i\n",*(pwing->wakelengths+i));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"htailcp=[");
    for (j=0;j<data.htail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->control+j),*(phtail->control+j+data.htail.nface),*(phtail->control+j+2*data.htail.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"htailcv=[");
    for (j=0;j<data.htail.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->vortex+j),*(phtail->vortex+j+data.htail.nvert),*(phtail->vortex+j+2*data.htail.nvert));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"htailnorm=[");
    for (j=0;j<data.htail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->normal+j),*(phtail->normal+j+data.htail.nface),*(phtail->normal+j+2*data.htail.nface));
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"htailtangx=[");
    for (j=0;j<data.htail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->tangx+j),*(phtail->tangx+j+data.htail.nface),*(phtail->tangx+j+2*data.htail.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"htailtangy=[");
    for (j=0;j<data.htail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->tangy+j),*(phtail->tangy+j+data.htail.nface),*(phtail->tangy+j+2*data.htail.nface));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"htaildxy=[");
    for (j=0;j<data.htail.nface;j++){
        fprintf(ofp,"%f %f\n",*(phtail->dxy+j),*(phtail->dxy+j+data.htail.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"htailnsurf=[");
    for (j=0;j<data.htail.nface;j++){
        fprintf(ofp,"%f\n",*(phtail->nsurf+j));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"htailDeltap=[");
    for (j=0;j<data.htail.nface;j++){
        fprintf(ofp,"%f\n",*(phtail->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"htailDeltad=[");
    for (j=0;j<data.htail.nface;j++){
        fprintf(ofp,"%f\n",*(phtail->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"htailwind=[");
    for (j=0;j<data.htail.nface;j++){
        fprintf(ofp,"%f\n",*(phtail->wind+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"htailgamma=[");
    for (i=0;i<data.htail.nface;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(phtail->gamma+i+j*data.htail.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(phtail->gamma+i+j*data.htail.nface));
            }
        }
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"phtail.xw=[");
    for (i=0;i<(phtail->nshed+phtail->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(phtail->xw+i+j*(phtail->nshed+phtail->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(phtail->xw+i+j*(phtail->nshed+phtail->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"phtail.yw=[");
    for (i=0;i<(phtail->nshed+phtail->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(phtail->yw+i+j*(phtail->nshed+phtail->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(phtail->yw+i+j*(phtail->nshed+phtail->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"phtail.zw=[");
    for (i=0;i<(phtail->nshed+phtail->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(phtail->zw+i+j*(phtail->nshed+phtail->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(phtail->zw+i+j*(phtail->nshed+phtail->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"phtail.uw=[");
    for (i=0;i<(phtail->nshed+phtail->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(phtail->uw+i+j*(phtail->nshed+phtail->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(phtail->uw+i+j*(phtail->nshed+phtail->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"phtail.vw=[");
    for (i=0;i<(phtail->nshed+phtail->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(phtail->vw+i+j*(phtail->nshed+phtail->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(phtail->vw+i+j*(phtail->nshed+phtail->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"phtail.ww=[");
    for (i=0;i<(phtail->nshed+phtail->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(phtail->ww+i+j*(phtail->nshed+phtail->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(phtail->ww+i+j*(phtail->nshed+phtail->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"phtail.gw=[");
    for (i=0;i<phtail->nshed;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(phtail->gw+i+j*phtail->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(phtail->gw+i+j*phtail->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");       
    
    fprintf(ofp,"htailneighbours=[");
    for (i=0;i<data.htail.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(phtail->neighbours+i+j*data.htail.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(phtail->neighbours+i+j*data.htail.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"phtail.nwakes=%i ;",phtail->nwakes);
    fprintf(ofp,"phtail.nshed=%i ;",phtail->nshed);
    fprintf(ofp,"phtail.wakeinds=[");
    for (i=0;i<data.htail.nshed;i++){
        fprintf(ofp,"%i\n",*(phtail->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"phtail.wakelengths=[");
    for (i=0;i<data.htail.nwakes;i++){
        fprintf(ofp,"%i\n",*(phtail->wakelengths+i));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"elevatorcp=[");
    for (j=0;j<data.elevator.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->control+j),*(pelevator->control+j+data.elevator.nface),*(pelevator->control+j+2*data.elevator.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"elevatorcv=[");
    for (j=0;j<data.elevator.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->vortex+j),*(pelevator->vortex+j+data.elevator.nvert),*(pelevator->vortex+j+2*data.elevator.nvert));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"elevatornorm=[");
    for (j=0;j<data.elevator.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->normal+j),*(pelevator->normal+j+data.elevator.nface),*(pelevator->normal+j+2*data.elevator.nface));
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"elevatortangx=[");
    for (j=0;j<data.elevator.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->tangx+j),*(pelevator->tangx+j+data.elevator.nface),*(pelevator->tangx+j+2*data.elevator.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"elevatortangy=[");
    for (j=0;j<data.elevator.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->tangy+j),*(pelevator->tangy+j+data.elevator.nface),*(pelevator->tangy+j+2*data.elevator.nface));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"elevatordxy=[");
    for (j=0;j<data.elevator.nface;j++){
        fprintf(ofp,"%f %f\n",*(pelevator->dxy+j),*(pelevator->dxy+j+data.elevator.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"elevatornsurf=[");
    for (j=0;j<data.elevator.nface;j++){
        fprintf(ofp,"%f\n",*(pelevator->nsurf+j));
    }
    fprintf(ofp,"];\n");     

    fprintf(ofp,"elevatorDeltap=[");
    for (j=0;j<data.elevator.nface;j++){
        fprintf(ofp,"%f\n",*(pelevator->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"elevatorDeltad=[");
    for (j=0;j<data.elevator.nface;j++){
        fprintf(ofp,"%f\n",*(pelevator->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"elevatorwind=[");
    for (j=0;j<data.elevator.nface;j++){
        fprintf(ofp,"%f\n",*(pelevator->wind+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"elevatorgamma=[");
    for (i=0;i<data.elevator.nface;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pelevator->gamma+i+j*data.elevator.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(pelevator->gamma+i+j*data.elevator.nface));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"pelevator.xw=[");
    for (i=0;i<(pelevator->nshed+pelevator->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pelevator->xw+i+j*(pelevator->nshed+pelevator->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pelevator->xw+i+j*(pelevator->nshed+pelevator->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"pelevator.yw=[");
    for (i=0;i<(pelevator->nshed+pelevator->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pelevator->yw+i+j*(pelevator->nshed+pelevator->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pelevator->yw+i+j*(pelevator->nshed+pelevator->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"pelevator.zw=[");
    for (i=0;i<(pelevator->nshed+pelevator->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pelevator->zw+i+j*(pelevator->nshed+pelevator->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pelevator->zw+i+j*(pelevator->nshed+pelevator->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"pelevator.uw=[");
    for (i=0;i<(pelevator->nshed+pelevator->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pelevator->uw+i+j*(pelevator->nshed+pelevator->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pelevator->uw+i+j*(pelevator->nshed+pelevator->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"pelevator.vw=[");
    for (i=0;i<(pelevator->nshed+pelevator->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pelevator->vw+i+j*(pelevator->nshed+pelevator->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pelevator->vw+i+j*(pelevator->nshed+pelevator->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"pelevator.ww=[");
    for (i=0;i<(pelevator->nshed+pelevator->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pelevator->ww+i+j*(pelevator->nshed+pelevator->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pelevator->ww+i+j*(pelevator->nshed+pelevator->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"pelevator.gw=[");
    for (i=0;i<pelevator->nshed;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pelevator->gw+i+j*pelevator->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(pelevator->gw+i+j*pelevator->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"elevatorneighbours=[");
    for (i=0;i<data.elevator.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pelevator->neighbours+i+j*data.elevator.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pelevator->neighbours+i+j*data.elevator.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"pelevator.nwakes=%i ;",pelevator->nwakes);
    fprintf(ofp,"pelevator.nshed=%i ;",pelevator->nshed);
    fprintf(ofp,"pelevator.wakeinds=[");
    for (i=0;i<data.elevator.nshed;i++){
        fprintf(ofp,"%i\n",*(pelevator->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"pelevator.wakelengths=[");
    for (i=0;i<data.elevator.nwakes;i++){
        fprintf(ofp,"%i\n",*(pelevator->wakelengths+i));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"vtailcp=[");
    for (j=0;j<data.vtail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->control+j),*(pvtail->control+j+data.vtail.nface),*(pvtail->control+j+2*data.vtail.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"vtailcv=[");
    for (j=0;j<data.vtail.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->vortex+j),*(pvtail->vortex+j+data.vtail.nvert),*(pvtail->vortex+j+2*data.vtail.nvert));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"vtailnorm=[");
    for (j=0;j<data.vtail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->normal+j),*(pvtail->normal+j+data.vtail.nface),*(pvtail->normal+j+2*data.vtail.nface));
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"vtailtangx=[");
    for (j=0;j<data.vtail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->tangx+j),*(pvtail->tangx+j+data.vtail.nface),*(pvtail->tangx+j+2*data.vtail.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"vtailtangy=[");
    for (j=0;j<data.vtail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->tangy+j),*(pvtail->tangy+j+data.vtail.nface),*(pvtail->tangy+j+2*data.vtail.nface));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"vtaildxy=[");
    for (j=0;j<data.vtail.nface;j++){
        fprintf(ofp,"%f %f\n",*(pvtail->dxy+j),*(pvtail->dxy+j+data.vtail.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"vtailnsurf=[");
    for (j=0;j<data.vtail.nface;j++){
        fprintf(ofp,"%f\n",*(pvtail->nsurf+j));
    }
    fprintf(ofp,"];\n");       

    fprintf(ofp,"vtailDeltap=[");
    for (j=0;j<data.vtail.nface;j++){
        fprintf(ofp,"%f\n",*(pvtail->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"vtailDeltad=[");
    for (j=0;j<data.vtail.nface;j++){
        fprintf(ofp,"%f\n",*(pvtail->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"vtailwind=[");
    for (j=0;j<data.vtail.nface;j++){
        fprintf(ofp,"%f\n",*(pvtail->wind+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"vtailgamma=[");
    for (i=0;i<data.vtail.nface;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pvtail->gamma+i+j*data.vtail.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(pvtail->gamma+i+j*data.vtail.nface));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"pvtail.xw=[");
    for (i=0;i<(pvtail->nshed+pvtail->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pvtail->xw+i+j*(pvtail->nshed+pvtail->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pvtail->xw+i+j*(pvtail->nshed+pvtail->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"pvtail.yw=[");
    for (i=0;i<(pvtail->nshed+pvtail->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pvtail->yw+i+j*(pvtail->nshed+pvtail->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pvtail->yw+i+j*(pvtail->nshed+pvtail->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"pvtail.zw=[");
    for (i=0;i<(pvtail->nshed+pvtail->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pvtail->zw+i+j*(pvtail->nshed+pvtail->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pvtail->zw+i+j*(pvtail->nshed+pvtail->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"pvtail.uw=[");
    for (i=0;i<(pvtail->nshed+pvtail->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pvtail->uw+i+j*(pvtail->nshed+pvtail->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pvtail->uw+i+j*(pvtail->nshed+pvtail->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"pvtail.vw=[");
    for (i=0;i<(pvtail->nshed+pvtail->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pvtail->vw+i+j*(pvtail->nshed+pvtail->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pvtail->vw+i+j*(pvtail->nshed+pvtail->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"pvtail.ww=[");
    for (i=0;i<(pvtail->nshed+pvtail->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pvtail->ww+i+j*(pvtail->nshed+pvtail->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(pvtail->ww+i+j*(pvtail->nshed+pvtail->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"pvtail.gw=[");
    for (i=0;i<pvtail->nshed;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(pvtail->gw+i+j*pvtail->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(pvtail->gw+i+j*pvtail->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");       
    
    fprintf(ofp,"vtailneighbours=[");
    for (i=0;i<data.vtail.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pvtail->neighbours+i+j*data.vtail.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pvtail->neighbours+i+j*data.vtail.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"pvtail.nwakes=%i ;",pvtail->nwakes);
    fprintf(ofp,"pvtail.nshed=%i ;",pvtail->nshed);
    fprintf(ofp,"pvtail.wakeinds=[");
    for (i=0;i<data.vtail.nshed;i++){
        fprintf(ofp,"%i\n",*(pvtail->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"pvtail.wakelengths=[");
    for (i=0;i<data.vtail.nwakes;i++){
        fprintf(ofp,"%i\n",*(pvtail->wakelengths+i));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"ruddercp=[");
    for (j=0;j<data.rudder.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->control+j),*(prudder->control+j+data.rudder.nface),*(prudder->control+j+2*data.rudder.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"ruddercv=[");
    for (j=0;j<data.rudder.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->vortex+j),*(prudder->vortex+j+data.rudder.nvert),*(prudder->vortex+j+2*data.rudder.nvert));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"ruddernorm=[");
    for (j=0;j<data.rudder.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->normal+j),*(prudder->normal+j+data.rudder.nface),*(prudder->normal+j+2*data.rudder.nface));
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"ruddertangx=[");
    for (j=0;j<data.rudder.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->tangx+j),*(prudder->tangx+j+data.rudder.nface),*(prudder->tangx+j+2*data.rudder.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"ruddertangy=[");
    for (j=0;j<data.rudder.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->tangy+j),*(prudder->tangy+j+data.rudder.nface),*(prudder->tangy+j+2*data.rudder.nface));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"rudderdxy=[");
    for (j=0;j<data.rudder.nface;j++){
        fprintf(ofp,"%f %f\n",*(prudder->dxy+j),*(prudder->dxy+j+data.rudder.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"ruddernsurf=[");
    for (j=0;j<data.rudder.nface;j++){
        fprintf(ofp,"%f\n",*(prudder->nsurf+j));
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"rudderDeltap=[");
    for (j=0;j<data.rudder.nface;j++){
        fprintf(ofp,"%f\n",*(prudder->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"rudderDeltad=[");
    for (j=0;j<data.rudder.nface;j++){
        fprintf(ofp,"%f\n",*(prudder->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"rudderwind=[");
    for (j=0;j<data.rudder.nface;j++){
        fprintf(ofp,"%f\n",*(prudder->wind+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"ruddergamma=[");
    for (i=0;i<data.rudder.nface;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(prudder->gamma+i+j*data.rudder.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(prudder->gamma+i+j*data.rudder.nface));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"prudder.xw=[");
    for (i=0;i<(prudder->nshed+prudder->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(prudder->xw+i+j*(prudder->nshed+prudder->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(prudder->xw+i+j*(prudder->nshed+prudder->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"prudder.yw=[");
    for (i=0;i<(prudder->nshed+prudder->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(prudder->yw+i+j*(prudder->nshed+prudder->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(prudder->yw+i+j*(prudder->nshed+prudder->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"prudder.zw=[");
    for (i=0;i<(prudder->nshed+prudder->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(prudder->zw+i+j*(prudder->nshed+prudder->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(prudder->zw+i+j*(prudder->nshed+prudder->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"prudder.uw=[");
    for (i=0;i<(prudder->nshed+prudder->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(prudder->uw+i+j*(prudder->nshed+prudder->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(prudder->uw+i+j*(prudder->nshed+prudder->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"prudder.vw=[");
    for (i=0;i<(prudder->nshed+prudder->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(prudder->vw+i+j*(prudder->nshed+prudder->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(prudder->vw+i+j*(prudder->nshed+prudder->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"prudder.ww=[");
    for (i=0;i<(prudder->nshed+prudder->nwakes);i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(prudder->ww+i+j*(prudder->nshed+prudder->nwakes)));
            }
            else{
                fprintf(ofp,"%f   ",*(prudder->ww+i+j*(prudder->nshed+prudder->nwakes)));
            }
        }
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"prudder.gw=[");
    for (i=0;i<prudder->nshed;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f\n",*(prudder->gw+i+j*prudder->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(prudder->gw+i+j*prudder->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");      
    
    fprintf(ofp,"rudderneighbours=[");
    for (i=0;i<data.rudder.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(prudder->neighbours+i+j*data.rudder.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(prudder->neighbours+i+j*data.rudder.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"prudder.nwakes=%i ;",prudder->nwakes);
    fprintf(ofp,"prudder.nshed=%i ;",prudder->nshed);
    fprintf(ofp,"prudder.wakeinds=[");
    for (i=0;i<data.rudder.nshed;i++){
        fprintf(ofp,"%i\n",*(prudder->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"prudder.wakelengths=[");
    for (i=0;i<data.rudder.nwakes;i++){
        fprintf(ofp,"%i\n",*(prudder->wakelengths+i));
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"AN=[");
    for (i=0;i<data.mtn;i++){
        for (j=0;j<data.mtn;j++){
            if (j == data.mtn-1){
                fprintf(ofp,"%f\n",*(data.AN+i+j*data.mtn));
            }
            else{
                fprintf(ofp,"%f   ",*(data.AN+i+j*data.mtn));
            }
        }
    }    
    fprintf(ofp,"];\n");

    fprintf(ofp,"invAN=[");
    for (i=0;i<data.mtn;i++){
        for (j=0;j<data.mtn;j++){
            if (j == data.mtn-1){
                fprintf(ofp,"%f\n",*(data.invAN+i+j*data.mtn));
            }
            else{
                fprintf(ofp,"%f   ",*(data.invAN+i+j*data.mtn));
            }
        }
    }    
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"BN=[");
    for (i=0;i<data.mtn;i++){
        for (j=0;j<data.mtn;j++){
            if (j == data.mtn-1){
                fprintf(ofp,"%f\n",*(data.BN+i+j*data.mtn));
            }
            else{
                fprintf(ofp,"%f   ",*(data.BN+i+j*data.mtn));
            }
        }
    }    
    fprintf(ofp,"];\n");

    fprintf(ofp,"RHS=[");
    for (j=0;j<data.mtn;j++){
        fprintf(ofp,"%f\n",*(data.RHS+j));
    }
    fprintf(ofp,"];\n");  
    
    fprintf(ofp,"Gammas=[");
    for (j=0;j<data.mtn;j++){
        fprintf(ofp,"%f\n",*(data.Gammas+j));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"wind=[");
    for (j=0;j<data.mtn;j++){
        fprintf(ofp,"%f\n",*(data.wind+j));
    }
    fprintf(ofp,"];\n");     
    /* Close output file */
     fclose(ofp);
     ofp = fopen("outfile.py", "w");
     if (ofp == NULL) {
        fprintf(stderr, "Can't open output file %s!\n","outfile.py");
        exit(1);
    }

    fprintf(ofp,"it=%i\n",it);
    fprintf(ofp,"dt=%f\n",data.dt);
    
    fprintf(ofp,"totalforce=[");
    for (i=0;i<4;i++){
        for (j=0;j<data.ntimes;j++){
            if (j == data.ntimes-1){
                fprintf(ofp,"%f,\n",*(data.totalforce+i+j*4));
            }
            else{
                fprintf(ofp,"%f,",*(data.totalforce+i+j*4));
            }
        }
    }
    fprintf(ofp,"]\n");
    fclose(ofp);
}
