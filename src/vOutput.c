//
//  vOutput.c
//
//  Data output
// mtn: total number of panels
//

#include "vOutput.h"
#include "vLiftsurf.h"
#include <stdio.h>
#include <stdlib.h>

void exportTextOutput(char *Outfile, int ntimes, int it, double dt, double *totalforce, struct liftsurf flap, struct liftsurf *pflap,
                    struct liftsurf aileron, struct liftsurf *paileron, struct liftsurf wing, struct liftsurf *pwing,
                    struct liftsurf htail, struct liftsurf *phtail, struct liftsurf vtail, struct liftsurf *pvtail,
                    struct liftsurf elevator, struct liftsurf *pelevator, struct liftsurf rudder, struct liftsurf *prudder,
                    int mtn, double *AN, double *invAN, double *BN, double *RHS, double *Gammas, double *wind)
{
    int i, j;
    FILE *ofp = fopen(Outfile, "w");
    if (ofp == NULL) {
        fprintf(stderr, "Can't open output file %s!\n",Outfile);
        exit(1);
    }

    fprintf(ofp,"it=%i;\n",it);
    fprintf(ofp,"dt=%f;\n",dt);
    
    fprintf(ofp,"totalforce=[");
    for (i=0;i<4;i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(totalforce+i+j*4));
            }
            else{
                fprintf(ofp,"%f   ",*(totalforce+i+j*4));
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"flap.faces=[");
    for (i=0;i<flap.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pflap->faces+i+j*flap.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pflap->faces+i+j*flap.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"aileron.faces=[");
    for (i=0;i<aileron.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(paileron->faces+i+j*aileron.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(paileron->faces+i+j*aileron.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"wing.faces=[");
    for (i=0;i<wing.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pwing->faces+i+j*wing.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pwing->faces+i+j*wing.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"htail.faces=[");
    for (i=0;i<htail.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(phtail->faces+i+j*htail.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(phtail->faces+i+j*htail.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"elevator.faces=[");
    for (i=0;i<elevator.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pelevator->faces+i+j*elevator.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pelevator->faces+i+j*elevator.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"vtail.faces=[");
    for (i=0;i<vtail.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pvtail->faces+i+j*vtail.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pvtail->faces+i+j*vtail.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"rudder.faces=[");
    for (i=0;i<rudder.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(prudder->faces+i+j*rudder.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(prudder->faces+i+j*rudder.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");


    fprintf(ofp,"flap.vertices=[");
    for (j=0;j<flap.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->vertices+j),*(pflap->vertices+j+flap.nvert),*(pflap->vertices+j+2*flap.nvert));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"aileron.vertices=[");
    for (j=0;j<aileron.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->vertices+j),*(paileron->vertices+j+aileron.nvert),*(paileron->vertices+j+2*aileron.nvert));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"wing.vertices=[");
    for (j=0;j<wing.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->vertices+j),*(pwing->vertices+j+wing.nvert),*(pwing->vertices+j+2*wing.nvert));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"htail.vertices=[");
    for (j=0;j<htail.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->vertices+j),*(phtail->vertices+j+htail.nvert),*(phtail->vertices+j+2*htail.nvert));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"elevator.vertices=[");
    for (j=0;j<elevator.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->vertices+j),*(pelevator->vertices+j+elevator.nvert),*(pelevator->vertices+j+2*elevator.nvert));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"vtail.vertices=[");
    for (j=0;j<vtail.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->vertices+j),*(pvtail->vertices+j+vtail.nvert),*(pvtail->vertices+j+2*vtail.nvert));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"rudder.vertices=[");
    for (j=0;j<rudder.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->vertices+j),*(prudder->vertices+j+rudder.nvert),*(prudder->vertices+j+2*rudder.nvert));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"pwing.uvw=[");
    for (j=0;j<wing.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->uvw+j),*(pwing->uvw+j+wing.nface),*(pwing->uvw+j+2*wing.nface));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"pflap.uvw=[");
    for (j=0;j<flap.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->uvw+j),*(pflap->uvw+j+flap.nface),*(pflap->uvw+j+2*flap.nface));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"paileron.uvw=[");
    for (j=0;j<aileron.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->uvw+j),*(paileron->uvw+j+aileron.nface),*(paileron->uvw+j+2*aileron.nface));
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"phtail.uvw=[");
    for (j=0;j<htail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->uvw+j),*(phtail->uvw+j+htail.nface),*(phtail->uvw+j+2*htail.nface));
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"pelevator.uvw=[");
    for (j=0;j<elevator.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->uvw+j),*(pelevator->uvw+j+elevator.nface),*(pelevator->uvw+j+2*elevator.nface));
    }
    fprintf(ofp,"];\n");     

    fprintf(ofp,"pvtail.uvw=[");
    for (j=0;j<vtail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->uvw+j),*(pvtail->uvw+j+vtail.nface),*(pvtail->uvw+j+2*vtail.nface));
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"prudder.uvw=[");
    for (j=0;j<rudder.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->uvw+j),*(prudder->uvw+j+rudder.nface),*(prudder->uvw+j+2*rudder.nface));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"flap_shedding=[");
    for (j=0;j<flap.nface;j++){
        fprintf(ofp,"%i\n",*(pflap->shedding+j));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"aileron_shedding=[");
    for (j=0;j<aileron.nface;j++){
        fprintf(ofp,"%i\n",*(paileron->shedding+j));
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"wing_shedding=[");
    for (j=0;j<wing.nface;j++){
        fprintf(ofp,"%i\n",*(pwing->shedding+j));
    }    
    fprintf(ofp,"];\n");

    fprintf(ofp,"flapcp=[");
    for (j=0;j<flap.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->control+j),*(pflap->control+j+flap.nface),*(pflap->control+j+2*flap.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"flapcv=[");
    for (j=0;j<flap.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->vortex+j),*(pflap->vortex+j+flap.nvert),*(pflap->vortex+j+2*flap.nvert));
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"flapnorm=[");
    for (j=0;j<flap.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->normal+j),*(pflap->normal+j+flap.nface),*(pflap->normal+j+2*flap.nface));
    }
    fprintf(ofp,"];\n");     

    fprintf(ofp,"flaptangx=[");
    for (j=0;j<flap.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->tangx+j),*(pflap->tangx+j+flap.nface),*(pflap->tangx+j+2*flap.nface));
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"flaptangy=[");
    for (j=0;j<flap.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pflap->tangy+j),*(pflap->tangy+j+flap.nface),*(pflap->tangy+j+2*flap.nface));
    }
    fprintf(ofp,"];\n");      
    
    fprintf(ofp,"flapdxy=[");
    for (j=0;j<flap.nface;j++){
        fprintf(ofp,"%f %f\n",*(pflap->dxy+j),*(pflap->dxy+j+flap.nface));
    }
    fprintf(ofp,"];\n");   
    
    fprintf(ofp,"flapnsurf=[");
    for (j=0;j<flap.nface;j++){
        fprintf(ofp,"%f\n",*(pflap->nsurf+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"flapDeltap=[");
    for (j=0;j<flap.nface;j++){
        fprintf(ofp,"%f\n",*(pflap->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"flapDeltad=[");
    for (j=0;j<flap.nface;j++){
        fprintf(ofp,"%f\n",*(pflap->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"flapwind=[");
    for (j=0;j<flap.nface;j++){
        fprintf(ofp,"%f\n",*(pflap->wind+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"flapgamma=[");
    for (i=0;i<flap.nface;i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(pflap->gamma+i+j*flap.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(pflap->gamma+i+j*flap.nface));
            }
        }
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"pflap.xw=[");
    for (i=0;i<(pflap->nshed+pflap->nwakes);i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(pflap->gw+i+j*pflap->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(pflap->gw+i+j*pflap->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"flapneighbours=[");
    for (i=0;i<flap.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pflap->neighbours+i+j*flap.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pflap->neighbours+i+j*flap.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"pflap.nwakes=%i ;",pflap->nwakes);
    fprintf(ofp,"pflap.nshed=%i ;",pflap->nshed);
    fprintf(ofp,"pflap.wakeinds=[");
    for (i=0;i<flap.nshed;i++){
        fprintf(ofp,"%i\n",*(pflap->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"pflap.wakelengths=[");
    for (i=0;i<flap.nwakes;i++){
        fprintf(ofp,"%i\n",*(pflap->wakelengths+i));
    }
    fprintf(ofp,"];\n");       
    
    fprintf(ofp,"aileroncp=[");
    for (j=0;j<aileron.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->control+j),*(paileron->control+j+aileron.nface),*(paileron->control+j+2*aileron.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"aileroncv=[");
    for (j=0;j<aileron.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->vortex+j),*(paileron->vortex+j+aileron.nvert),*(paileron->vortex+j+2*aileron.nvert));
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"aileronnorm=[");
    for (j=0;j<aileron.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->normal+j),*(paileron->normal+j+aileron.nface),*(paileron->normal+j+2*aileron.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"ailerontangx=[");
    for (j=0;j<aileron.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->tangx+j),*(paileron->tangx+j+aileron.nface),*(paileron->tangx+j+2*aileron.nface));
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"ailerontangy=[");
    for (j=0;j<aileron.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(paileron->tangy+j),*(paileron->tangy+j+aileron.nface),*(paileron->tangy+j+2*aileron.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"ailerondxy=[");
    for (j=0;j<aileron.nface;j++){
        fprintf(ofp,"%f %f\n",*(paileron->dxy+j),*(paileron->dxy+j+aileron.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"aileronnsurf=[");
    for (j=0;j<aileron.nface;j++){
        fprintf(ofp,"%f\n",*(paileron->nsurf+j));
    }
    fprintf(ofp,"];\n");  

    fprintf(ofp,"aileronDeltap=[");
    for (j=0;j<aileron.nface;j++){
        fprintf(ofp,"%f\n",*(paileron->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"aileronDeltad=[");
    for (j=0;j<aileron.nface;j++){
        fprintf(ofp,"%f\n",*(paileron->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"aileronwind=[");
    for (j=0;j<aileron.nface;j++){
        fprintf(ofp,"%f\n",*(paileron->wind+j));
    }
    fprintf(ofp,"];\n");   
    
    fprintf(ofp,"ailerongamma=[");
    for (i=0;i<aileron.nface;i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(paileron->gamma+i+j*aileron.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(paileron->gamma+i+j*aileron.nface));
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"paileron.xw=[");
    for (i=0;i<(paileron->nshed+paileron->nwakes);i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(paileron->gw+i+j*paileron->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(paileron->gw+i+j*paileron->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"aileronneighbours=[");
    for (i=0;i<aileron.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(paileron->neighbours+i+j*aileron.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(paileron->neighbours+i+j*aileron.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"paileron.nwakes=%i ;",paileron->nwakes);
    fprintf(ofp,"paileron.nshed=%i ;",paileron->nshed);
    fprintf(ofp,"paileron.wakeinds=[");
    for (i=0;i<aileron.nshed;i++){
        fprintf(ofp,"%i\n",*(paileron->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"paileron.wakelengths=[");
    for (i=0;i<aileron.nwakes;i++){
        fprintf(ofp,"%i\n",*(paileron->wakelengths+i));
    }
    fprintf(ofp,"];\n");    
    
    
    fprintf(ofp,"wingcp=[");
    for (j=0;j<wing.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->control+j),*(pwing->control+j+wing.nface),*(pwing->control+j+2*wing.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"wingcv=[");
    for (j=0;j<wing.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->vortex+j),*(pwing->vortex+j+wing.nvert),*(pwing->vortex+j+2*wing.nvert));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"wingnorm=[");
    for (j=0;j<wing.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->normal+j),*(pwing->normal+j+wing.nface),*(pwing->normal+j+2*wing.nface));
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"wingtangx=[");
    for (j=0;j<wing.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->tangx+j),*(pwing->tangx+j+wing.nface),*(pwing->tangx+j+2*wing.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"wingtangy=[");
    for (j=0;j<wing.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pwing->tangy+j),*(pwing->tangy+j+wing.nface),*(pwing->tangy+j+2*wing.nface));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"wingdxy=[");
    for (j=0;j<wing.nface;j++){
        fprintf(ofp,"%f %f\n",*(pwing->dxy+j),*(pwing->dxy+j+wing.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"wingnsurf=[");
    for (j=0;j<wing.nface;j++){
        fprintf(ofp,"%f\n",*(pwing->nsurf+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"wingDeltap=[");
    for (j=0;j<wing.nface;j++){
        fprintf(ofp,"%f\n",*(pwing->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"wingDeltad=[");
    for (j=0;j<wing.nface;j++){
        fprintf(ofp,"%f\n",*(pwing->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"wingwind=[");
    for (j=0;j<wing.nface;j++){
        fprintf(ofp,"%f\n",*(pwing->wind+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"winggamma=[");
    for (i=0;i<wing.nface;i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(pwing->gamma+i+j*wing.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(pwing->gamma+i+j*wing.nface));
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"pwing.xw=[");
    for (i=0;i<(pwing->nshed+pwing->nwakes);i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(pwing->gw+i+j*pwing->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(pwing->gw+i+j*pwing->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"wingneighbours=[");
    for (i=0;i<wing.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pwing->neighbours+i+j*wing.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pwing->neighbours+i+j*wing.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"pwing.nwakes=%i ;",pwing->nwakes);
    fprintf(ofp,"pwing.nshed=%i ;",pwing->nshed);
    fprintf(ofp,"pwing.wakeinds=[");
    for (i=0;i<wing.nshed;i++){
        fprintf(ofp,"%i\n",*(pwing->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"pwing.wakelengths=[");
    for (i=0;i<wing.nwakes;i++){
        fprintf(ofp,"%i\n",*(pwing->wakelengths+i));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"htailcp=[");
    for (j=0;j<htail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->control+j),*(phtail->control+j+htail.nface),*(phtail->control+j+2*htail.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"htailcv=[");
    for (j=0;j<htail.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->vortex+j),*(phtail->vortex+j+htail.nvert),*(phtail->vortex+j+2*htail.nvert));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"htailnorm=[");
    for (j=0;j<htail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->normal+j),*(phtail->normal+j+htail.nface),*(phtail->normal+j+2*htail.nface));
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"htailtangx=[");
    for (j=0;j<htail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->tangx+j),*(phtail->tangx+j+htail.nface),*(phtail->tangx+j+2*htail.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"htailtangy=[");
    for (j=0;j<htail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(phtail->tangy+j),*(phtail->tangy+j+htail.nface),*(phtail->tangy+j+2*htail.nface));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"htaildxy=[");
    for (j=0;j<htail.nface;j++){
        fprintf(ofp,"%f %f\n",*(phtail->dxy+j),*(phtail->dxy+j+htail.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"htailnsurf=[");
    for (j=0;j<htail.nface;j++){
        fprintf(ofp,"%f\n",*(phtail->nsurf+j));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"htailDeltap=[");
    for (j=0;j<htail.nface;j++){
        fprintf(ofp,"%f\n",*(phtail->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"htailDeltad=[");
    for (j=0;j<htail.nface;j++){
        fprintf(ofp,"%f\n",*(phtail->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"htailwind=[");
    for (j=0;j<htail.nface;j++){
        fprintf(ofp,"%f\n",*(phtail->wind+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"htailgamma=[");
    for (i=0;i<htail.nface;i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(phtail->gamma+i+j*htail.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(phtail->gamma+i+j*htail.nface));
            }
        }
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"phtail.xw=[");
    for (i=0;i<(phtail->nshed+phtail->nwakes);i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(phtail->gw+i+j*phtail->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(phtail->gw+i+j*phtail->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");       
    
    fprintf(ofp,"htailneighbours=[");
    for (i=0;i<htail.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(phtail->neighbours+i+j*htail.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(phtail->neighbours+i+j*htail.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"phtail.nwakes=%i ;",phtail->nwakes);
    fprintf(ofp,"phtail.nshed=%i ;",phtail->nshed);
    fprintf(ofp,"phtail.wakeinds=[");
    for (i=0;i<htail.nshed;i++){
        fprintf(ofp,"%i\n",*(phtail->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"phtail.wakelengths=[");
    for (i=0;i<htail.nwakes;i++){
        fprintf(ofp,"%i\n",*(phtail->wakelengths+i));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"elevatorcp=[");
    for (j=0;j<elevator.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->control+j),*(pelevator->control+j+elevator.nface),*(pelevator->control+j+2*elevator.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"elevatorcv=[");
    for (j=0;j<elevator.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->vortex+j),*(pelevator->vortex+j+elevator.nvert),*(pelevator->vortex+j+2*elevator.nvert));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"elevatornorm=[");
    for (j=0;j<elevator.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->normal+j),*(pelevator->normal+j+elevator.nface),*(pelevator->normal+j+2*elevator.nface));
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"elevatortangx=[");
    for (j=0;j<elevator.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->tangx+j),*(pelevator->tangx+j+elevator.nface),*(pelevator->tangx+j+2*elevator.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"elevatortangy=[");
    for (j=0;j<elevator.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pelevator->tangy+j),*(pelevator->tangy+j+elevator.nface),*(pelevator->tangy+j+2*elevator.nface));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"elevatordxy=[");
    for (j=0;j<elevator.nface;j++){
        fprintf(ofp,"%f %f\n",*(pelevator->dxy+j),*(pelevator->dxy+j+elevator.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"elevatornsurf=[");
    for (j=0;j<elevator.nface;j++){
        fprintf(ofp,"%f\n",*(pelevator->nsurf+j));
    }
    fprintf(ofp,"];\n");     

    fprintf(ofp,"elevatorDeltap=[");
    for (j=0;j<elevator.nface;j++){
        fprintf(ofp,"%f\n",*(pelevator->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"elevatorDeltad=[");
    for (j=0;j<elevator.nface;j++){
        fprintf(ofp,"%f\n",*(pelevator->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"elevatorwind=[");
    for (j=0;j<elevator.nface;j++){
        fprintf(ofp,"%f\n",*(pelevator->wind+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"elevatorgamma=[");
    for (i=0;i<elevator.nface;i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(pelevator->gamma+i+j*elevator.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(pelevator->gamma+i+j*elevator.nface));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"pelevator.xw=[");
    for (i=0;i<(pelevator->nshed+pelevator->nwakes);i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(pelevator->gw+i+j*pelevator->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(pelevator->gw+i+j*pelevator->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"elevatorneighbours=[");
    for (i=0;i<elevator.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pelevator->neighbours+i+j*elevator.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pelevator->neighbours+i+j*elevator.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"pelevator.nwakes=%i ;",pelevator->nwakes);
    fprintf(ofp,"pelevator.nshed=%i ;",pelevator->nshed);
    fprintf(ofp,"pelevator.wakeinds=[");
    for (i=0;i<elevator.nshed;i++){
        fprintf(ofp,"%i\n",*(pelevator->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"pelevator.wakelengths=[");
    for (i=0;i<elevator.nwakes;i++){
        fprintf(ofp,"%i\n",*(pelevator->wakelengths+i));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"vtailcp=[");
    for (j=0;j<vtail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->control+j),*(pvtail->control+j+vtail.nface),*(pvtail->control+j+2*vtail.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"vtailcv=[");
    for (j=0;j<vtail.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->vortex+j),*(pvtail->vortex+j+vtail.nvert),*(pvtail->vortex+j+2*vtail.nvert));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"vtailnorm=[");
    for (j=0;j<vtail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->normal+j),*(pvtail->normal+j+vtail.nface),*(pvtail->normal+j+2*vtail.nface));
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"vtailtangx=[");
    for (j=0;j<vtail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->tangx+j),*(pvtail->tangx+j+vtail.nface),*(pvtail->tangx+j+2*vtail.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"vtailtangy=[");
    for (j=0;j<vtail.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(pvtail->tangy+j),*(pvtail->tangy+j+vtail.nface),*(pvtail->tangy+j+2*vtail.nface));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"vtaildxy=[");
    for (j=0;j<vtail.nface;j++){
        fprintf(ofp,"%f %f\n",*(pvtail->dxy+j),*(pvtail->dxy+j+vtail.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"vtailnsurf=[");
    for (j=0;j<vtail.nface;j++){
        fprintf(ofp,"%f\n",*(pvtail->nsurf+j));
    }
    fprintf(ofp,"];\n");       

    fprintf(ofp,"vtailDeltap=[");
    for (j=0;j<vtail.nface;j++){
        fprintf(ofp,"%f\n",*(pvtail->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"vtailDeltad=[");
    for (j=0;j<vtail.nface;j++){
        fprintf(ofp,"%f\n",*(pvtail->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"vtailwind=[");
    for (j=0;j<vtail.nface;j++){
        fprintf(ofp,"%f\n",*(pvtail->wind+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"vtailgamma=[");
    for (i=0;i<vtail.nface;i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(pvtail->gamma+i+j*vtail.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(pvtail->gamma+i+j*vtail.nface));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"pvtail.xw=[");
    for (i=0;i<(pvtail->nshed+pvtail->nwakes);i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(pvtail->gw+i+j*pvtail->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(pvtail->gw+i+j*pvtail->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");       
    
    fprintf(ofp,"vtailneighbours=[");
    for (i=0;i<vtail.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(pvtail->neighbours+i+j*vtail.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(pvtail->neighbours+i+j*vtail.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"pvtail.nwakes=%i ;",pvtail->nwakes);
    fprintf(ofp,"pvtail.nshed=%i ;",pvtail->nshed);
    fprintf(ofp,"pvtail.wakeinds=[");
    for (i=0;i<vtail.nshed;i++){
        fprintf(ofp,"%i\n",*(pvtail->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"pvtail.wakelengths=[");
    for (i=0;i<vtail.nwakes;i++){
        fprintf(ofp,"%i\n",*(pvtail->wakelengths+i));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"ruddercp=[");
    for (j=0;j<rudder.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->control+j),*(prudder->control+j+rudder.nface),*(prudder->control+j+2*rudder.nface));
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"ruddercv=[");
    for (j=0;j<rudder.nvert;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->vortex+j),*(prudder->vortex+j+rudder.nvert),*(prudder->vortex+j+2*rudder.nvert));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"ruddernorm=[");
    for (j=0;j<rudder.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->normal+j),*(prudder->normal+j+rudder.nface),*(prudder->normal+j+2*rudder.nface));
    }
    fprintf(ofp,"];\n"); 
    
    fprintf(ofp,"ruddertangx=[");
    for (j=0;j<rudder.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->tangx+j),*(prudder->tangx+j+rudder.nface),*(prudder->tangx+j+2*rudder.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"ruddertangy=[");
    for (j=0;j<rudder.nface;j++){
        fprintf(ofp,"%f %f %f\n",*(prudder->tangy+j),*(prudder->tangy+j+rudder.nface),*(prudder->tangy+j+2*rudder.nface));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"rudderdxy=[");
    for (j=0;j<rudder.nface;j++){
        fprintf(ofp,"%f %f\n",*(prudder->dxy+j),*(prudder->dxy+j+rudder.nface));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"ruddernsurf=[");
    for (j=0;j<rudder.nface;j++){
        fprintf(ofp,"%f\n",*(prudder->nsurf+j));
    }
    fprintf(ofp,"];\n");      

    fprintf(ofp,"rudderDeltap=[");
    for (j=0;j<rudder.nface;j++){
        fprintf(ofp,"%f\n",*(prudder->Deltap+j));
    }
    fprintf(ofp,"];\n");

    fprintf(ofp,"rudderDeltad=[");
    for (j=0;j<rudder.nface;j++){
        fprintf(ofp,"%f\n",*(prudder->Deltad+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"rudderwind=[");
    for (j=0;j<rudder.nface;j++){
        fprintf(ofp,"%f\n",*(prudder->wind+j));
    }
    fprintf(ofp,"];\n");    
    
    fprintf(ofp,"ruddergamma=[");
    for (i=0;i<rudder.nface;i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(prudder->gamma+i+j*rudder.nface));
            }
            else{
                fprintf(ofp,"%f   ",*(prudder->gamma+i+j*rudder.nface));
            }
        }
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"prudder.xw=[");
    for (i=0;i<(prudder->nshed+prudder->nwakes);i++){
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
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
        for (j=0;j<ntimes;j++){
            if (j == ntimes-1){
                fprintf(ofp,"%f\n",*(prudder->gw+i+j*prudder->nshed));
            }
            else{
                fprintf(ofp,"%f   ",*(prudder->gw+i+j*prudder->nshed));
            }
        }
    }
    fprintf(ofp,"];\n");      
    
    fprintf(ofp,"rudderneighbours=[");
    for (i=0;i<rudder.nface;i++){
        for (j=0;j<4;j++){
            if (j == 3){
                fprintf(ofp,"%i\n",*(prudder->neighbours+i+j*rudder.nface)+1);
            }
            else{
                fprintf(ofp,"%i   ",*(prudder->neighbours+i+j*rudder.nface)+1);
            }
        }
    }
    fprintf(ofp,"];\n");    

    fprintf(ofp,"prudder.nwakes=%i ;",prudder->nwakes);
    fprintf(ofp,"prudder.nshed=%i ;",prudder->nshed);
    fprintf(ofp,"prudder.wakeinds=[");
    for (i=0;i<rudder.nshed;i++){
        fprintf(ofp,"%i\n",*(prudder->wakeinds+i)+1);
    }
    fprintf(ofp,"];\n");  
    fprintf(ofp,"prudder.wakelengths=[");
    for (i=0;i<rudder.nwakes;i++){
        fprintf(ofp,"%i\n",*(prudder->wakelengths+i));
    }
    fprintf(ofp,"];\n");     
    
    fprintf(ofp,"AN=[");
    for (i=0;i<mtn;i++){
        for (j=0;j<mtn;j++){
            if (j == mtn-1){
                fprintf(ofp,"%f\n",*(AN+i+j*mtn));
            }
            else{
                fprintf(ofp,"%f   ",*(AN+i+j*mtn));
            }
        }
    }    
    fprintf(ofp,"];\n");

    fprintf(ofp,"invAN=[");
    for (i=0;i<mtn;i++){
        for (j=0;j<mtn;j++){
            if (j == mtn-1){
                fprintf(ofp,"%f\n",*(invAN+i+j*mtn));
            }
            else{
                fprintf(ofp,"%f   ",*(invAN+i+j*mtn));
            }
        }
    }    
    fprintf(ofp,"];\n");
    
    fprintf(ofp,"BN=[");
    for (i=0;i<mtn;i++){
        for (j=0;j<mtn;j++){
            if (j == mtn-1){
                fprintf(ofp,"%f\n",*(BN+i+j*mtn));
            }
            else{
                fprintf(ofp,"%f   ",*(BN+i+j*mtn));
            }
        }
    }    
    fprintf(ofp,"];\n");

    fprintf(ofp,"RHS=[");
    for (j=0;j<mtn;j++){
        fprintf(ofp,"%f\n",*(RHS+j));
    }
    fprintf(ofp,"];\n");  
    
    fprintf(ofp,"Gammas=[");
    for (j=0;j<mtn;j++){
        fprintf(ofp,"%f\n",*(Gammas+j));
    }
    fprintf(ofp,"];\n"); 

    fprintf(ofp,"wind=[");
    for (j=0;j<mtn;j++){
        fprintf(ofp,"%f\n",*(wind+j));
    }
    fprintf(ofp,"];\n");     
    /* Close output file */
     fclose(ofp);    
}
