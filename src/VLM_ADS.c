#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vAssignValues.h"
#include "vControl.h"
#include "vInput.h"
#include "vOutput.h"
#include "vForce.h"
#include "vGeometry.h"
#include "vInfluence.h"
#include "vLiftsurf.h"
#include "vRHS.h"
#include "vVortex.h"
#include "vWake.h"
#include "vHTail.h"
#include "vVTail.h"
#include "vWing.h"


void adaptycoord(double *ypline, double *ypos, int np1, int nypos)
{
    /* Adapt values of vector ypline so that they include the values of vector ypos */
    int i,j,nmin;
    double mindist,dist;
    
    nmin=-1;
    for (i=0;i<nypos;i++){
        mindist=1000.0;
        for (j=nmin+1;j<np1;j++){
            dist=fabs(*(ypos+i)- *(ypline+j)); /* Look for the value of ypline nearest to this element of ypos */
            if (dist <= mindist){
                mindist=dist;
                nmin=j;
            }
        }
        *(ypline+nmin)=*(ypos+i); /* Change the nearest value of ypos */
    }
}

void findneighbours(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, int cwing, int cflap, int caileron)
{
    /* Find all the neighbours of all the wing panels in all wing, flap and aileron panels */
    int i,j,k,l,repeat[4];
    double vertx[4],verty[4],vertz[4],dist[3],totaldist;
    
    pwing->neighbours=(int *)malloc(sizeof(int)*pwing->nface*4);
    for (i=0;i<pwing->nface;i++){
        for (k=0;k<4;k++){
            vertx[k]=*(pwing->vertices+ *(pwing->faces+i+k*pwing->nface));
            verty[k]=*(pwing->vertices+pwing->nvert+ *(pwing->faces+i+k*pwing->nface));
            vertz[k]=*(pwing->vertices+2*pwing->nvert+ *(pwing->faces+i+k*pwing->nface));
            *(pwing->neighbours+i+k*pwing->nface)=-1;
        }
        /* Check if wing panels have neighbours on the wing */
        for (j=0;j<pwing->nface;j++){
            if (i != j){
                for (k=0;k<4;k++){
                    repeat[k]=0;
                    for (l=0;l<4;l++){
                        dist[0]=vertx[k]-*(pwing->vertices+ *(pwing->faces+j+l*pwing->nface));
                        dist[1]=verty[k]-*(pwing->vertices+pwing->nvert+ *(pwing->faces+j+l*pwing->nface));
                        dist[2]=vertz[k]-*(pwing->vertices+2*pwing->nvert+ *(pwing->faces+j+l*pwing->nface));
                        totaldist=sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
                        if (totaldist < 0.001)
                            repeat[k]=j+1; /* +1 so that j=0 does not results in repeat[k]=0 */
                    }
                }
                if (repeat[0] != 0 && repeat[1] !=0){
                    if (i<pwing->nface/2){
                        *(pwing->neighbours+i+pwing->nface)=cwing+j; /* right neighbour looking from leading edge to trailing edge*/
                    }else
                        *(pwing->neighbours+i)=cwing+j; /* left neighbour looking from leading edge to trailing edge*/
                }
                if (repeat[1] != 0 && repeat[2] !=0){
                    *(pwing->neighbours+i+3*pwing->nface)=cwing+j; /* downstream neighbour looking from leading edge to trailing edge*/
                }
                if (repeat[2] != 0 && repeat[3] !=0){
                    if (i<pwing->nface/2){
                        *(pwing->neighbours+i)=cwing+j; /* left neighbour looking from leading edge to trailing edge*/
                    }else
                        *(pwing->neighbours+i+pwing->nface)=cwing+j; /* right neighbour looking from leading edge to trailing edge*/
                }
                if (repeat[3] != 0 && repeat[0] !=0){
                    *(pwing->neighbours+i+2*pwing->nface)=cwing+j; /* upstream neighbour looking from leading edge to trailing edge*/
                }
            }
        }
        /* Now check if wing panels have neighbours on the flap */
        for (j=0;j<pflap->nface;j++){
            for (k=0;k<4;k++){
                repeat[k]=0;
                for (l=0;l<4;l++){
                    dist[0]=vertx[k]-*(pflap->vertices+ *(pflap->faces+j+l*pflap->nface));
                    dist[1]=verty[k]-*(pflap->vertices+pflap->nvert+ *(pflap->faces+j+l*pflap->nface));
                    dist[2]=vertz[k]-*(pflap->vertices+2*pflap->nvert+ *(pflap->faces+j+l*pflap->nface));
                    totaldist=sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
                    if (totaldist < 0.001)
                        repeat[k]=j+1; /* +1 so that j=0 does not result in repeat[k]=0 */
                }
            }
            if (repeat[0] != 0 && repeat[1] !=0){
                if (i<pwing->nface/2){
                    *(pwing->neighbours+i+pwing->nface)=cflap+j; /* right neighbour looking from leading edge to trailing edge*/
                }else
                    *(pwing->neighbours+i)=cflap+j; /* left neighbour looking from leading edge to trailing edge*/
            }
            if (repeat[1] != 0 && repeat[2] !=0){
                *(pwing->neighbours+i+3*pwing->nface)=cflap+j; /* downstream neighbour looking from leading edge to trailing edge*/
            }
            if (repeat[2] != 0 && repeat[3] !=0){
                if (i<pwing->nface/2){
                    *(pwing->neighbours+i)=cflap+j; /* left neighbour looking from leading edge to trailing edge*/
                }else
                    *(pwing->neighbours+i+pwing->nface)=cflap+j; /* right neighbour looking from leading edge to trailing edge*/
            }
            if (repeat[3] != 0 && repeat[0] !=0){
                *(pwing->neighbours+i+2*pwing->nface)=cflap+j; /* upstream neighbour looking from leading edge to trailing edge*/
            }
        }
        /* Now check if wing panels have neighbours on the aileron */
        for (j=0;j<paileron->nface;j++){
            for (k=0;k<4;k++){
                repeat[k]=0;
                for (l=0;l<4;l++){
                    dist[0]=vertx[k]-*(paileron->vertices+ *(paileron->faces+j+l*paileron->nface));
                    dist[1]=verty[k]-*(paileron->vertices+paileron->nvert+ *(paileron->faces+j+l*paileron->nface));
                    dist[2]=vertz[k]-*(paileron->vertices+2*paileron->nvert+ *(paileron->faces+j+l*paileron->nface));
                    totaldist=sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
                    if (totaldist < 0.001)
                        repeat[k]=j+1; /* +1 so that j=0 does not results in repeat[k]=0 */
                }
            }
            if (repeat[0] != 0 && repeat[1] !=0){
                if (i<pwing->nface/2){
                    *(pwing->neighbours+i+pwing->nface)=caileron+j; /* right neighbour looking from leading edge to trailing edge*/
                }else
                    *(pwing->neighbours+i)=caileron+j; /* left neighbour looking from leading edge to trailing edge*/
            }
            if (repeat[1] != 0 && repeat[2] !=0){
                *(pwing->neighbours+i+3*pwing->nface)=caileron+j; /* downstream neighbour looking from leading edge to trailing edge*/
            }
            if (repeat[2] != 0 && repeat[3] !=0){
                if (i<pwing->nface/2){
                    *(pwing->neighbours+i)=caileron+j; /* left neighbour looking from leading edge to trailing edge*/
                }else
                    *(pwing->neighbours+i+pwing->nface)=caileron+j; /* right neighbour looking from leading edge to trailing edge*/
            }
            if (repeat[3] != 0 && repeat[0] !=0){
                *(pwing->neighbours+i+2*pwing->nface)=caileron+j; /* upstream neighbour looking from leading edge to trailing edge*/
            }
        }
    }
}

void gauss(double *a, double *inva, int nrows)
{
    /* Carry out Gaussian elimination on matrix a to calculate its inverse, inva */
    /* Inspired by the algorithm in Numerical Recipes in C */
    int i,j,k;
    double *MI, factor,dummy;
    
    MI = (double *)malloc(sizeof(double)*(2*nrows)*(nrows));
    
    /* Create augmented matrix */
    for (j=0; j<nrows; j++) {
        for (i=0;i<nrows; i++){
            *(MI+i+j*nrows)=*(a+j*nrows+i);
            if (i==j){
                *(MI+i+(j+nrows)*nrows)=1.0;
            } else {
                *(MI+i+(j+nrows)*nrows)=0.0;
            }
        }
    }
    
    /* Create inverse matrix */
    for (i=0; i<nrows; i++) {
        dummy=*(MI+i+i*nrows);
        for (j=0; j<2*nrows; j++) {
            *(MI+i+j*nrows)=*(MI+i+j*nrows)/dummy;
        }
        for (k=0; k<nrows; k++) {
            if (k!=i){
                factor=- *(MI+k+i*nrows)/ *(MI+i+i*nrows);
                for (j=0; j<2*nrows; j++) {
                    *(MI+k+j*nrows) = *(MI+k+j*nrows) + factor* *(MI+i+j*nrows);
                }
            }
        }
    }
    
    /* Extract inverse matrix */
    for (j=0; j<nrows; j++) {
        for(i=0; i<nrows; i++) {
            *(inva+i+nrows*j)=*(MI+i+(j+nrows)*nrows);
        }
    }
    free(MI);
}

int main(int argc, char *argv[])
{
    char *Infile, Wngfile[60], HTailfile[60], VTailfile[60], *Outfile;
    struct liftsurf flap,aileron,wing,htail,elevator,vtail,rudder;
    struct liftsurf *pflap = &flap;
    struct liftsurf *paileron = &aileron;
    struct liftsurf *pwing = &wing;
    struct liftsurf *phtail = &htail;
    struct liftsurf *pelevator = &elevator;
    struct liftsurf *pvtail = &vtail;
    struct liftsurf *prudder = &rudder;
    double *invAN,*AN,*BN,*RHS, *Gammas, *wind;
    double UVW[3],aoa,yaw,dt,timestep_denom,MAC,rho,*totalforce;
    double delta[2],beta[2],eta[2],zeta[2];
    int i,j,mtn,it,ntimes,m,n,freewake,mht,nht,mvt,nvt;
    
    if (argc>1) {
        Infile=argv[1];
    } else {
        Infile="infile.arp";
    }
    
    Outfile="outfile.m";
    
    importInputFile(Infile, UVW, &rho, &aoa, &yaw, &m, &mht, &mvt, &n, &nht, &nvt, &ntimes, &timestep_denom, &freewake,
                    delta, beta, eta, zeta, Wngfile, HTailfile, VTailfile);
        
    totalforce=(double *)calloc(ntimes*4, sizeof(double));
    
    /* Create the lifting surfaces that make up the wing */
    MAC=wingsetup(pflap,paileron,pwing,Wngfile,m,n);
    /* Create the lifting surfaces that make up the horizontal tail */
    htailsetup(phtail,pelevator,HTailfile,mht,nht);
    /* Create the lifting surfaces that make up the vertical tail */
    vtailsetup(pvtail,prudder,VTailfile,mht,nht);
    dt=MAC/UVW[0]/timestep_denom; /* dt is based on the length of the wing's Mean Aerodynamic Chord */
    printf("MAC=%f, dt=%f\n",MAC,dt);
    
    /* Rotate ailerons */
    rotateail(paileron,delta);
    /* Rotate flaps */
    rotateail(pflap,beta);
    /* Rotate elevator */
    rotateail(pelevator,eta);
    /* Rotate rudder */
    rotateail(prudder,zeta);
    
    findneighbours(pwing,pflap,paileron,0,10000,20000);
    findneighbours(pflap,pwing,paileron,10000,0,20000);
    findneighbours(paileron,pwing,pflap,20000,0,10000);
    findneighbours(phtail,pelevator,paileron,30000,40000,20000); /* We don't care about paileron being called here */
    findneighbours(pelevator,phtail,paileron,40000,30000,20000); /* We don't care about paileron being called here */
    findneighbours(pvtail,prudder,paileron,50000,60000,20000); /* We don't care about paileron being called here */
    findneighbours(prudder,pvtail,paileron,60000,50000,20000); /* We don't care about paileron being called here */

    /* Create the wake */
    createwake(pwing,0,ntimes);
    createwake(pflap,10000,ntimes);
    createwake(paileron,20000,ntimes); 
    createwake(phtail,30000,ntimes); 
    createwake(pelevator,40000,ntimes); 
    createwake(pvtail,50000,ntimes); 
    createwake(prudder,60000,ntimes); 
    
    /* Find which bound vortex panels correspond to which wake panel vertices */
    correspshedwake(pwing);
    correspshedwake(pflap);
    correspshedwake(paileron);
    correspshedwake(phtail);
    correspshedwake(pelevator);
    correspshedwake(pvtail);
    correspshedwake(prudder);
    
    /* Shed first wake element */
    shedwake(pwing);
    shedwake(pflap);
    shedwake(paileron);
    shedwake(phtail);
    shedwake(pelevator);
    shedwake(pvtail);
    shedwake(prudder);
    findallwakeneighbours(pwing,pflap,paileron,0,10000,20000);
    findallwakeneighbours(phtail,pelevator,paileron,30000,40000,20000); /* We don't care about paileron being called here */
    findallwakeneighbours(pvtail,prudder,paileron,50000,60000,20000); /* We don't care about paileron being called here */
    
    /* Calculate collocation points and vortex segment lengths */
    colvec(pflap);
    colvec(paileron);
    colvec(pwing);
    colvec(phtail);
    colvec(pelevator);
    colvec(pvtail);
    colvec(prudder);
    /* Calculate normal vectors and surfaces */
    normals(pflap);
    normals(paileron);
    normals(pwing);
    normals(phtail);
    normals(pelevator);
    normals(pvtail);
    normals(prudder);
    /* Calculate tangential vectors */
    tangentials(pflap);
    tangentials(paileron);
    tangentials(pwing);
    tangentials(phtail);
    tangentials(pelevator);
    tangentials(pvtail);
    tangentials(prudder);
    /* Calculate total number of panels */
    mtn=flap.nface+wing.nface+aileron.nface+htail.nface+elevator.nface+vtail.nface+rudder.nface;
    
    /* Assign memory for global matrices */
    AN=(double *)calloc(mtn*mtn, sizeof(double)); /* Normal flow coefficient matrix */
    invAN=(double *)malloc(sizeof(double)*mtn*mtn); /* Inverse of normal flow coefficient matrix */
    BN=(double *)calloc(mtn*mtn, sizeof(double)); /* Downwash coefficient matrix */
    RHS=(double *)malloc(sizeof(double)*mtn); /* Right Hand Side vector */
    Gammas = (double *)malloc(sizeof(double)*mtn); /* Vorticity vector history */
    wind=(double *)malloc(sizeof(double)*mtn); /* Induced windspeeds */

    /* Calculate influence coefficient matrices */
    cycleliftsurf(pwing,pflap,paileron,phtail,pelevator,pvtail,prudder,AN,BN,mtn);    

    /* Invert the influence coefficient matrix */
    gauss(AN,invAN,mtn);

    /* Assign memory for gamma and wind values in every liftsurf */
    pwing->gamma=(double *)calloc(pwing->nface*ntimes, sizeof(double));
    pwing->wind=(double *)calloc(pwing->nface, sizeof(double));
    pflap->gamma=(double *)calloc(pflap->nface*ntimes, sizeof(double));
    pflap->wind=(double *)calloc(pflap->nface, sizeof(double));
    paileron->gamma=(double *)calloc(paileron->nface*ntimes, sizeof(double));
    paileron->wind=(double *)calloc(paileron->nface, sizeof(double));
    phtail->gamma=(double *)calloc(phtail->nface*ntimes, sizeof(double));
    phtail->wind=(double *)calloc(phtail->nface, sizeof(double));
    pelevator->gamma=(double *)calloc(pelevator->nface*ntimes, sizeof(double));
    pelevator->wind=(double *)calloc(pelevator->nface, sizeof(double));
    pvtail->gamma=(double *)calloc(pvtail->nface*ntimes, sizeof(double));
    pvtail->wind=(double *)calloc(pvtail->nface, sizeof(double));
    prudder->gamma=(double *)calloc(prudder->nface*ntimes, sizeof(double));
    prudder->wind=(double *)calloc(prudder->nface, sizeof(double));
    
    /* Assign memory AND INITIALISE uvw values in every liftsurf */
    pwing->uvw=(double *)calloc(pwing->nface * 3, sizeof(double));
    pflap->uvw=(double *)calloc(pflap->nface * 3, sizeof(double));
    paileron->uvw=(double *)calloc(paileron->nface * 3, sizeof(double));
    phtail->uvw = (double *)calloc(phtail->nface * 3, sizeof(double));
    pelevator->uvw=(double *)calloc(pelevator->nface * 3, sizeof(double));
    pvtail->uvw = (double *)calloc(pvtail->nface * 3, sizeof(double));
    prudder->uvw=(double *)calloc(prudder->nface * 3, sizeof(double));
    
    /* Assign memory for Deltap, Deltad and aeroforce on every liftsurf */
    pwing->Deltap=(double *)malloc(sizeof(double)*pwing->nface);
    pflap->Deltap=(double *)malloc(sizeof(double)*pflap->nface);
    paileron->Deltap=(double *)malloc(sizeof(double)*paileron->nface);
    phtail->Deltap=(double *)malloc(sizeof(double)*phtail->nface);
    pelevator->Deltap=(double *)malloc(sizeof(double)*pelevator->nface);
    pvtail->Deltap=(double *)malloc(sizeof(double)*pvtail->nface);
    prudder->Deltap=(double *)malloc(sizeof(double)*prudder->nface);
    pwing->Deltad=(double *)malloc(sizeof(double)*pwing->nface);
    pflap->Deltad=(double *)malloc(sizeof(double)*pflap->nface);
    paileron->Deltad=(double *)malloc(sizeof(double)*paileron->nface);
    phtail->Deltad=(double *)malloc(sizeof(double)*phtail->nface);
    pelevator->Deltad=(double *)malloc(sizeof(double)*pelevator->nface);
    pvtail->Deltad=(double *)malloc(sizeof(double)*pvtail->nface);
    prudder->Deltad=(double *)malloc(sizeof(double)*prudder->nface);
    pwing->aeroforce=(double *)malloc(sizeof(double)*4);
    pflap->aeroforce=(double *)malloc(sizeof(double)*4);
    paileron->aeroforce=(double *)malloc(sizeof(double)*4);
    phtail->aeroforce=(double *)malloc(sizeof(double)*4);
    pelevator->aeroforce=(double *)malloc(sizeof(double)*4);
    pvtail->aeroforce=(double *)malloc(sizeof(double)*4);
    prudder->aeroforce=(double *)malloc(sizeof(double)*4);
    
    for (it=0;it<ntimes-1;it++){
        printf("it=%i\n",it);
        if (it > 0){
            /* Calculate influences of all the wings on all the liftsurfs */
            calcwakeinf(pwing,pflap,paileron,phtail,pelevator,pvtail,prudder,it);
        }
        /* Calculate Right Hand Side vector */
        calcRHS(pwing,pflap,paileron,phtail,pelevator,pvtail,prudder,RHS,UVW);
        
        /* Solve for panel vorticity */
        for (i=0; i < mtn; i++) {
            *(Gammas+i)=0.0;
            for (j=0; j < mtn; j++) {
                *(Gammas+i) += *(invAN+i+j*mtn)* *(RHS+j);
            }
        }

        /* Calculate induced airspeeds */
        for (i=0; i < mtn; i++) {
            *(wind+i)=0.0;
            for (j=0; j < mtn; j++) {
                *(wind+i) += *(BN+i+j*mtn)* *(Gammas+j); /* Downwash by vortices on surface */
            }
        }
        
        /* Assign Gammas and wind onto respective lifting surfaces */
        assignGammas(Gammas,pwing,pflap,paileron,phtail,pelevator,pvtail,prudder,it);
        assignwind(wind,pwing,pflap,paileron,phtail,pelevator,pvtail,prudder);       
        
        /* Propagate wake vortex strength */
        propwakevort(pwing,dt,it,UVW);
        propwakevort(pflap,dt,it,UVW);
        propwakevort(paileron,dt,it,UVW);
        propwakevort(phtail,dt,it,UVW);
        propwakevort(pelevator,dt,it,UVW);
        propwakevort(pvtail,dt,it,UVW);
        propwakevort(prudder,dt,it,UVW);
        
        /* Assign vortex strength of trailing edge panels to the first wake panels */
        wakegamma(pwing,it);
        wakegamma(pflap,it);
        wakegamma(paileron,it);         
        wakegamma(phtail,it);         
        wakegamma(pelevator,it);         
        wakegamma(pvtail,it);         
        wakegamma(prudder,it);         
        
        /* Add influence of latest wake element on wind */
        addlastwind(pwing,pflap,paileron,it,3);
        addlastwind(phtail,pelevator,paileron,it,2); /* paileron is irrelevant */
        addlastwind(pvtail,prudder,paileron,it,2); /* paileron is irrelevant */
        
        /* Calculate the pressure difference across the panels in every liftsurf */
        calcforces(pwing,pwing,pflap,paileron,it,dt,UVW,rho,0,10000,20000);
        calcforces(pflap,pwing,pflap,paileron,it,dt,UVW,rho,0,10000,20000);
        calcforces(paileron,pwing,pflap,paileron,it,dt,UVW,rho,0,10000,20000);
        calcforceshtail(phtail,phtail,pelevator,it,dt,UVW,rho,30000,40000);
        calcforceshtail(pelevator,phtail,pelevator,it,dt,UVW,rho,30000,40000);
        calcforcesvtail(pvtail,pvtail,prudder,it,dt,UVW,rho,50000,60000);
        calcforcesvtail(prudder,pvtail,prudder,it,dt,UVW,rho,50000,60000);
       
        /* Calculate local airspeeds on wake */
        if (it>0 && freewake!=0){
            infonwake(pwing,pwing,it,0);
            infonwake(pwing,pflap,it,1);
            infonwake(pwing,paileron,it,1);
            infonwake(pwing,phtail,it,1);
            infonwake(pwing,pelevator,it,1);
            infonwake(pflap,pwing,it,0);
            infonwake(pflap,pflap,it,1);
            infonwake(pflap,paileron,it,1);
            infonwake(pflap,phtail,it,1);
            infonwake(pflap,pelevator,it,1);
            infonwake(paileron,pwing,it,0);
            infonwake(paileron,pflap,it,1);
            infonwake(paileron,paileron,it,1);
            infonwake(paileron,phtail,it,1);
            infonwake(paileron,pelevator,it,1);
            infonwake(phtail,pwing,it,0);
            infonwake(phtail,pflap,it,1);
            infonwake(phtail,paileron,it,1);
            infonwake(phtail,phtail,it,1);
            infonwake(phtail,pelevator,it,1);
            infonwake(pelevator,pwing,it,0);
            infonwake(pelevator,pflap,it,1);
            infonwake(pelevator,paileron,it,1);
            infonwake(pelevator,phtail,it,1);
            infonwake(pelevator,pelevator,it,1);
            /* vtail and rudder are treated as uncoupled */
            infonwake(pvtail,pvtail,it,0);
            infonwake(pvtail,prudder,it,1);
            infonwake(prudder,pvtail,it,0);
            infonwake(prudder,prudder,it,1);
       }

        /* Propagate wake position */
        propwakexyz(pwing,dt,it,UVW);
        propwakexyz(pflap,dt,it,UVW);
        propwakexyz(paileron,dt,it,UVW);
        propwakexyz(phtail,dt,it,UVW);
        propwakexyz(pelevator,dt,it,UVW);
        propwakexyz(pvtail,dt,it,UVW);
        propwakexyz(prudder,dt,it,UVW);
        
        /* Calculate the total forces on all the liftsurfs */
        for (i=0;i<4;i++){
            *(totalforce+i+it*4)=*(pwing->aeroforce+i)+*(pflap->aeroforce+i)+*(paileron->aeroforce+i)+*(phtail->aeroforce+i)+*(pelevator->aeroforce+i)+*(pvtail->aeroforce+i)+*(prudder->aeroforce+i);
        }
    }   
    
    it--;
    printf("forcex=%f, forcey=%f, forcez=%f, draginduced=%f\n",*(totalforce+0+it*4),*(totalforce+1+it*4),*(totalforce+2+it*4),*(totalforce+3+it*4));
    
    /* Print output */
    exportTextOutput(Outfile, ntimes, it, dt, totalforce, flap, pflap, aileron, paileron, wing, pwing, htail, phtail, vtail, pvtail,
                    elevator, pelevator, rudder, prudder, mtn, AN, invAN, BN, RHS, Gammas, wind);
}