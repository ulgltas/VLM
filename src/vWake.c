//
//  vWake.c
//
//  Wake related functions
//
//

#include "vWake.h"
#include "vLiftsurf.h"
#include "vVLMData.h"
#include "vVortex.h"
#include <stdlib.h>

void storewakeinds(struct liftsurf *pwing)
{
    /* Find out how many panels shed a wake in this liftsurf */
    int i,wakeshed;
    
    wakeshed=0;
    for (i=0;i<pwing->nface;i++){
        if (*(pwing->neighbours+i+pwing->nface*3) == -1){
            wakeshed++;
            *(pwing->wakeinds-1+wakeshed)=i;
            //printf("%i\n",*(pwing->wakeinds-1+wakeshed));
        }
    }
}

int findwakeneighbours(struct liftsurf *pwing, int ipanel, int lrneighs[],int cwing)
{
    /* Find the neighbours of ipanel in liftsurf pwing */
    int j, hasneigh,neigh;
    
    hasneigh=0;
    lrneighs[0]=-1; /* left neighbour */
    lrneighs[1]=-1; /* right neighbour */
    /* Check if left neighbour is shedding */
    neigh=*(pwing->neighbours+*(pwing->wakeinds+ipanel))-cwing;
    for (j=0;j<pwing->nshed;j++){
        if (*(pwing->wakeinds+j)==neigh){
            hasneigh++;
            lrneighs[0]=j;
        }
    }
    /* Check if right neighbour is shedding */
    neigh=*(pwing->neighbours+*(pwing->wakeinds+ipanel)+pwing->nface)-cwing;
    for (j=0;j<pwing->nshed;j++){
        if (*(pwing->wakeinds+j)==neigh){
            hasneigh++;
            lrneighs[1]=j;
        }
    }
    return hasneigh;
}

void setupwakes(struct liftsurf *pwing, int cwing)
{
    /* Find number of distinct wake surfaces shed from a liftsurf */
    int i, j, k,nwakes,lrneighs[2],nneighs,wakelength,beginwake,*dummywl,*dummywi,*dummywio,total,ln,rn;
    
    dummywl=(int *)malloc(sizeof(int)*pwing->nshed); /* Here we store the wake lengths */
    dummywi=(int *)malloc(sizeof(int)*pwing->nshed); /* Copy of wake indices */
    dummywio=(int *)malloc(sizeof(int)*pwing->nshed); /* Ordered wake indices */
    for (i=0;i<pwing->nshed;i++){
        *(dummywi+i)=*(pwing->wakeinds+i);
    }
    
    /* Work out how many contiguous series of wake panels are shed */
    wakelength=1;
    beginwake=0;
    total=0;
    nwakes=1;
    for (i=0;i<pwing->nshed;i++){
        if (*(dummywi+i) != -1){
            *(dummywio+wakelength-1+total)=*(pwing->wakeinds+i); /* Store in the ordered vector this panel */
            *(dummywi+i)=-1; /* Don't consider this panel again */
            nneighs=findwakeneighbours(pwing,i,lrneighs,cwing);
            ln=lrneighs[0]; /* Left neighbour */
            rn=lrneighs[1]; /* Right neighbour */           
            /* Follow the left neighbours first */
            while (lrneighs[0] > -1){
               wakelength++;  
                *(dummywio+wakelength-1+total)=*(pwing->wakeinds+lrneighs[0]); /* Store in the ordered vector the left neighbours */
                k=lrneighs[0];
                if (*(dummywi+k)==-1){
                    lrneighs[0]=-1;
                }else{
                    *(dummywi+lrneighs[0])=-1; /* Do not consider this panel for future wakes */
                    nneighs=findwakeneighbours(pwing,k,lrneighs,cwing); /* Find the neighbours of the neighbour */
                }
            }
            /* Now follow the right neighbours */
            lrneighs[1]=rn;
            while (lrneighs[1] > -1){
                wakelength++;  
                *(dummywio+wakelength-1+total)=*(pwing->wakeinds+lrneighs[1]); /* Store in the ordered vector the left neighbours */
                k=lrneighs[1];
                if (*(dummywi+k)==-1){
                    lrneighs[1]=-1;
                }else{
                    *(dummywi+lrneighs[1])=-1; /* Do not consider this panel for future wakes */
                    nneighs=findwakeneighbours(pwing,k,lrneighs,cwing); /* Find the neighbours of the neighbour */
                }
            }
            /* If we get here, it means that we have found all the neighbours of the neighbours */
            *(dummywl+nwakes-1)=wakelength; /* Store the length of this wake */
            nwakes++; /* Increment number of wakes */
            total+=wakelength;
            wakelength=1;
        }
    }
    /* Set up number of wakes and wake lengths in liftsurf */
    pwing->nwakes=nwakes-1;
    pwing->wakelengths=(int *)malloc(sizeof(int)*nwakes);
    /* Store wake lengths */
    for (i=0;i<nwakes;i++){
        *(pwing->wakelengths+i)=*(dummywl+i);
    }
    for (i=0;i<pwing->nshed;i++){
        *(pwing->wakeinds+i)=*(dummywio+i);
    }
    free(dummywl);
    free(dummywi);
    free(dummywio); 
}

void createwake(struct liftsurf *pwing, double cwing, int ntimes)
{
    /* Create wakes shed from a liftsurf */
    int i,nwakes;
    
    if (pwing->nshed !=0 ){
        pwing->wakeinds=(int *)malloc(sizeof(int)*pwing->nshed);
        storewakeinds(pwing);
    }
    setupwakes(pwing,cwing);
    pwing->xw=(double *)calloc((pwing->nshed + pwing->nwakes)*ntimes, sizeof(double));
    pwing->yw=(double *)calloc((pwing->nshed + pwing->nwakes)*ntimes, sizeof(double));
    pwing->zw=(double *)calloc((pwing->nshed + pwing->nwakes)*ntimes, sizeof(double));
    pwing->uw=(double *)calloc((pwing->nshed + pwing->nwakes)*ntimes, sizeof(double));
    pwing->vw=(double *)calloc((pwing->nshed + pwing->nwakes)*ntimes, sizeof(double));
    pwing->ww=(double *)calloc((pwing->nshed + pwing->nwakes)*ntimes, sizeof(double));
    pwing->gw=(double *)calloc(pwing->nshed * ntimes, sizeof(double));
    /* Set contents of pwing->xw,yw,zw to zero */
    for (i=0;i<(pwing->nshed+pwing->nwakes)*ntimes;i++){
        *(pwing->xw+i)=0.0;
        *(pwing->yw+i)=0.0;
        *(pwing->zw+i)=0.0;
    }
    /* Set contents of pwing->gw to zero */
    for (i=0;i<pwing->nshed*ntimes;i++){
        *(pwing->gw+i)=0.0;
    }
}

void correspshedwake(struct liftsurf *pwing)
{
    /* Shed the first wake element */
    int i,j,total,nother,nother2,total2,panel;
    
    pwing->wakecorrespx=(int *)malloc(sizeof(int)*(pwing->nshed+pwing->nwakes));
    pwing->wakecorrespy=(int *)malloc(sizeof(int)*(pwing->nshed+pwing->nwakes));
    pwing->wakecorrespz=(int *)malloc(sizeof(int)*(pwing->nshed+pwing->nwakes));
    total=0;
    total2=0;
    /* Go through all the wakes in this liftsurf */
    for (i=0;i<pwing->nwakes;i++){
        panel=*(pwing->wakeinds+total); /* first panel in this wake */
        /* Check if the first panel lies on the left half-wing */
        if (panel < pwing->nface/2){
            nother=0;
            nother2=0;
            /* Find how many panels of this wake lie on the right half-wing */
            for (j=1;j<*(pwing->wakelengths+i);j++){
                if (*(pwing->wakeinds+total+j) >= pwing->nface/2){
                    nother++;
                }
            }
            /* Shed both trailing vertices of this panel */
            *(pwing->wakecorrespx+total2+nother)=*(pwing->faces+panel+pwing->nface);
            *(pwing->wakecorrespy+total2+nother)=*(pwing->faces+panel+pwing->nface)+pwing->nvert;
            *(pwing->wakecorrespz+total2+nother)=*(pwing->faces+panel+pwing->nface)+2*pwing->nvert;
            *(pwing->wakecorrespx+total2+1+nother)=*(pwing->faces+panel+2*pwing->nface);
            *(pwing->wakecorrespy+total2+1+nother)=*(pwing->faces+panel+2*pwing->nface)+pwing->nvert;
            *(pwing->wakecorrespz+total2+1+nother)=*(pwing->faces+panel+2*pwing->nface)+2*pwing->nvert;
            /* Continue shedding the other panels */
            for (j=1;j<*(pwing->wakelengths+i);j++){
                panel=*(pwing->wakeinds+total+j);
                if (panel < pwing->nface/2){
                    *(pwing->wakecorrespx+total2+1+nother+j)=*(pwing->faces+panel+2*pwing->nface);
                    *(pwing->wakecorrespy+total2+1+nother+j)=*(pwing->faces+panel+2*pwing->nface)+pwing->nvert;
                    *(pwing->wakecorrespz+total2+1+nother+j)=*(pwing->faces+panel+2*pwing->nface)+2*pwing->nvert;
                    nother2++; /* Increase the number of left side panels we've shed (not counting the first one) */
                }else{
                    *(pwing->wakecorrespx+total2+nother+nother2-j)=*(pwing->faces+panel+2*pwing->nface);
                    *(pwing->wakecorrespy+total2+nother+nother2-j)=*(pwing->faces+panel+2*pwing->nface)+pwing->nvert;
                    *(pwing->wakecorrespz+total2+nother+nother2-j)=*(pwing->faces+panel+2*pwing->nface)+2*pwing->nvert;
                }
            }
        }else{
            /* In this case all the panels will lie on the right, no need to check if there are any on the left */
            /* Order panels so that vertices go from right to left */
            /* Shed both trailing vertices of this panel */
            *(pwing->wakecorrespx+total2+*(pwing->wakelengths+i))=*(pwing->faces+panel+pwing->nface);
            *(pwing->wakecorrespy+total2+*(pwing->wakelengths+i))=*(pwing->faces+panel+pwing->nface)+pwing->nvert;
            *(pwing->wakecorrespz+total2+*(pwing->wakelengths+i))=*(pwing->faces+panel+pwing->nface)+2*pwing->nvert;
            *(pwing->wakecorrespx+total2+*(pwing->wakelengths+i)-1)=*(pwing->faces+panel+2*pwing->nface);
            *(pwing->wakecorrespy+total2+*(pwing->wakelengths+i)-1)=*(pwing->faces+panel+2*pwing->nface)+pwing->nvert;
            *(pwing->wakecorrespz+total2+*(pwing->wakelengths+i)-1)=*(pwing->faces+panel+2*pwing->nface)+2*pwing->nvert;
            /* Continue shedding the other panels */
            for (j=1;j<*(pwing->wakelengths+i);j++){
                panel=*(pwing->wakeinds+total+j);
                *(pwing->wakecorrespx+total2+*(pwing->wakelengths+i)-1-j)=*(pwing->faces+panel+2*pwing->nface);
                *(pwing->wakecorrespy+total2+*(pwing->wakelengths+i)-1-j)=*(pwing->faces+panel+2*pwing->nface)+pwing->nvert;
                *(pwing->wakecorrespz+total2+*(pwing->wakelengths+i)-1-j)=*(pwing->faces+panel+2*pwing->nface)+2*pwing->nvert;
            }
        }
        total+=*(pwing->wakelengths+i);
        total2+=*(pwing->wakelengths+i)+1;
    }
}

void shedwake(struct liftsurf *pwing)
{
    /* Shed the first wake element */
    int i;
    
    for (i=0;i<pwing->nshed+pwing->nwakes;i++){
        *(pwing->xw+i)=*(pwing->vortex+ *(pwing->wakecorrespx+i));
        *(pwing->yw+i)=*(pwing->vortex+ *(pwing->wakecorrespy+i));
        *(pwing->zw+i)=*(pwing->vortex+ *(pwing->wakecorrespz+i));
    }
}

void wakegamma(struct liftsurf *pwing, int it)
{
    /* Assign vortex strength to first row of wake panels */
    int i,j,total,nother,nother2,panel;

    total=0;
    for (i=0;i<pwing->nwakes;i++){
        panel=*(pwing->wakeinds+total); /* first panel in this wake */
        /* Check if the first panel lies on the left half-wing */
        if (panel < pwing->nface/2){
            nother=0;
            nother2=0;
            /* Find how many panels of this wake lie on the right half-wing */
            for (j=1;j<*(pwing->wakelengths+i);j++){
                if (*(pwing->wakeinds+total+j) >= pwing->nface/2){
                    nother++;
                }
            }
            /* Shed vorticity of this panel into the wake */
            *(pwing->gw+total+nother)=*(pwing->gamma+panel+it*pwing->nface);
            /* Continue shedding the other panels */
            for (j=1;j<*(pwing->wakelengths+i);j++){
                panel=*(pwing->wakeinds+total+j);
                if (panel < pwing->nface/2){
                    *(pwing->gw+total+nother+j)=*(pwing->gamma+panel+it*pwing->nface);
                    nother2++; /* Increase the number of left side panels we've shed (not counting the first one) */
                }else{
                    *(pwing->gw+total+nother+nother2-j)= *(pwing->gamma+panel+it*pwing->nface); /* vortex strength must be of the same sign as on the other side */
                }
            }
        }else{
            /* In this case all the panels will lie on the right, no need to check if there are any on the left */
            *(pwing->gw+total+*(pwing->wakelengths+i)-1)= *(pwing->gamma+panel+it*pwing->nface); /* vortex strength must be of the same sign as on the other side */
            /* Continue shedding the other panels */
            for (j=1;j<*(pwing->wakelengths+i);j++){
                panel=*(pwing->wakeinds+total+j);
                *(pwing->gw+total+*(pwing->wakelengths+i)-1-j)= *(pwing->gamma+panel+it*pwing->nface); /* vortex strength must be of the same sign as on the other side */
            }
        }
        total+=*(pwing->wakelengths+i);
    }    
}

void propwakexyz(struct liftsurf *pwing, double dt, int it, double UVW[])
{
    /* Propagate the x,y and z coordinates of the wake in this liftsurf */
    int i,j,npoints;
    
    /* Propagate vertices */
    npoints=pwing->nshed+pwing->nwakes;
    for (j=it;j>-1;j--){
        for (i=0;i<npoints;i++){
//             *(pwing->xw+i+(j+1)*npoints)=*(pwing->xw+i+j*npoints)+UVW[0]*dt;
//             *(pwing->yw+i+(j+1)*npoints)=*(pwing->yw+i+j*npoints)+UVW[1]*dt;
//             *(pwing->zw+i+(j+1)*npoints)=*(pwing->zw+i+j*npoints)+UVW[2]*dt;
            *(pwing->xw+i+(j+1)*npoints)=*(pwing->xw+i+j*npoints)+(UVW[0]+ *(pwing->uw+i+j*npoints))*dt;
            *(pwing->yw+i+(j+1)*npoints)=*(pwing->yw+i+j*npoints)+(UVW[1]+ *(pwing->vw+i+j*npoints))*dt;
            *(pwing->zw+i+(j+1)*npoints)=*(pwing->zw+i+j*npoints)+(UVW[2]+ *(pwing->ww+i+j*npoints))*dt;
        }
    }
}

void propwakevort(struct liftsurf *pwing, double dt, int it, double UVW[])
{
    /* Propagate the vortex strength of the wake in this liftsurf */
    int i,j;
    
    if (it > 0){
        /* Propagate vortex strengths */
        for (j=it-1;j>-1;j--){
            for (i=0;i<pwing->nshed;i++){
                *(pwing->gw+i+(j+1)*pwing->nshed)=*(pwing->gw+i+j*pwing->nshed);
            }
        }
    }
}

void infonwake(struct liftsurf *pwing, struct liftsurf *pflap, int it, int summode)
{
    /* Calculate the influence of pflap (bound and wake vorticity) on the wake of pwing */
    int i, j, i2, j2, m, mp1, start, start2, nw, ppp;  
    double uvw1[3],uvw2[3],uvw3[3],uvw4[3],uvw[3],xc,yc,zc;
    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
        
    m=pwing->nshed+pwing->nwakes;
    mp1=pflap->nshed+pflap->nwakes;
    for (i=0; i < m; i++) {
        for (j=0;j<it+1;j++){
            uvw[0]=0.0;
            uvw[1]=0.0;
            uvw[2]=0.0;
            /* Wake grid point */
            xc=*(pwing->xw+i+j*m);
            yc=*(pwing->yw+i+j*m);
            zc=*(pwing->zw+i+j*m); 
//            printf("%i %i %f  %f %f\n",i,j,xc,yc,zc);
            /* Influence of bound vorticity */
            for (i2=0;i2<pflap->nface;i2++){
                x1=*(pflap->vortex+*(pflap->faces+i2));
                x2=*(pflap->vortex+*(pflap->faces+i2+pflap->nface));
                x3=*(pflap->vortex+*(pflap->faces+i2+2*pflap->nface));
                x4=*(pflap->vortex+*(pflap->faces+i2+3*pflap->nface));
                y1=*(pflap->vortex+pflap->nvert+*(pflap->faces+i2));
                y2=*(pflap->vortex+pflap->nvert+*(pflap->faces+i2+pflap->nface));
                y3=*(pflap->vortex+pflap->nvert+*(pflap->faces+i2+2*pflap->nface));
                y4=*(pflap->vortex+pflap->nvert+*(pflap->faces+i2+3*pflap->nface));
                z1=*(pflap->vortex+2*pflap->nvert+*(pflap->faces+i2));
                z2=*(pflap->vortex+2*pflap->nvert+*(pflap->faces+i2+pflap->nface));
                z3=*(pflap->vortex+2*pflap->nvert+*(pflap->faces+i2+2*pflap->nface));
                z4=*(pflap->vortex+2*pflap->nvert+*(pflap->faces+i2+3*pflap->nface));
                if (i2 < pflap->nface/2){
                    /* Influence of vortex segment 1 */
                    vortex(uvw1,xc,yc,zc,x1,y1,z1,x4,y4,z4,*(pflap->gamma+i2+it*pflap->nface));
                    /* Influence of vortex segment 2 */
                    vortex(uvw2,xc,yc,zc,x4,y4,z4,x3,y3,z3,*(pflap->gamma+i2+it*pflap->nface));
                    /* Influence of vortex segment 3 */
                    vortex(uvw3,xc,yc,zc,x3,y3,z3,x2,y2,z2,*(pflap->gamma+i2+it*pflap->nface));
                    /* Influence of vortex segment 4 */
                    vortex(uvw4,xc,yc,zc,x2,y2,z2,x1,y1,z1,*(pflap->gamma+i2+it*pflap->nface));
                }else{
                    /* Influence of vortex segment 1 */
                    vortex(uvw1,xc,yc,zc,x4,y4,z4,x1,y1,z1,*(pflap->gamma+i2+it*pflap->nface));
                    /* Influence of vortex segment 2 */
                    vortex(uvw2,xc,yc,zc,x1,y1,z1,x2,y2,z2,*(pflap->gamma+i2+it*pflap->nface));
                    /* Influence of vortex segment 3 */
                    vortex(uvw3,xc,yc,zc,x2,y2,z2,x3,y3,z3,*(pflap->gamma+i2+it*pflap->nface));
                    /* Influence of vortex segment 4 */
                    vortex(uvw4,xc,yc,zc,x3,y3,z3,x4,y4,z4,*(pflap->gamma+i2+it*pflap->nface));
                }
                /* Addd the four influences */
                uvw[0]+=uvw1[0]+uvw2[0]+uvw3[0]+uvw4[0];
                uvw[1]+=uvw1[1]+uvw2[1]+uvw3[1]+uvw4[1];
                uvw[2]+=uvw1[2]+uvw2[2]+uvw3[2]+uvw4[2];
            }
            /* Influence of wake vorticity */
            for (nw=0;nw<pflap->nwakes;nw++){
                /* Find start point in pflap wake*/
                start=0;
                start2=0;
                if (nw > 0){
                    for (i2=0;i2<nw;i2++){
                        start+=*(pflap->wakelengths+i2);
                        start2+=*(pflap->wakelengths+i2)+1;
                    }
                }
                for (i2=0; i2 < *(pflap->wakelengths+nw); i2++) {
                    for (j2=0; j2 < it; j2++) {
                        ppp=start+i2+j2*pflap->nshed;
                        /* Influence of vortex segment 1 */
                        vortexblob(uvw1,xc,yc,zc,
                                *(pflap->xw+start2+i2+1+j2*mp1),*(pflap->yw+start2+i2+1+j2*mp1),*(pflap->zw+start2+i2+1+j2*mp1),
                                *(pflap->xw+start2+i2+1+(j2+1)*mp1),*(pflap->yw+start2+i2+1+(j2+1)*mp1),*(pflap->zw+start2+i2+1+(j2+1)*mp1),
                                *(pflap->gw+ppp));
                        /* Influence of vortex segment 2 */
                        vortexblob(uvw2,xc,yc,zc,
                                *(pflap->xw+start2+i2+1+(j2+1)*mp1),*(pflap->yw+start2+i2+1+(j2+1)*mp1),*(pflap->zw+start2+i2+1+(j2+1)*mp1),
                                *(pflap->xw+start2+i2+(j2+1)*mp1),*(pflap->yw+start2+i2+(j2+1)*mp1),*(pflap->zw+start2+i2+(j2+1)*mp1),
                                *(pflap->gw+ppp));
                        /* Influence of vortex segment 3 */
                        vortexblob(uvw3,xc,yc,zc,
                                *(pflap->xw+start2+i2+(j2+1)*mp1),*(pflap->yw+start2+i2+(j2+1)*mp1),*(pflap->zw+start2+i2+(j2+1)*mp1),
                                *(pflap->xw+start2+i2+j2*mp1),*(pflap->yw+start2+i2+j2*mp1),*(pflap->zw+start2+i2+j2*mp1),
                                *(pflap->gw+ppp));
                        /* Influence of vortex segment 4 */
                        vortexblob(uvw4,xc,yc,zc,
                                *(pflap->xw+start2+i2+j2*mp1),*(pflap->yw+start2+i2+j2*mp1),*(pflap->zw+start2+i2+j2*mp1),
                                *(pflap->xw+start2+i2+1+j2*mp1),*(pflap->yw+start2+i2+1+j2*mp1),*(pflap->zw+start2+i2+1+j2*mp1),
                                *(pflap->gw+ppp));
                        /* Addd the four influences */
                        uvw[0]+=uvw1[0]+uvw2[0]+uvw3[0]+uvw4[0];
                        uvw[1]+=uvw1[1]+uvw2[1]+uvw3[1]+uvw4[1];
                        uvw[2]+=uvw1[2]+uvw2[2]+uvw3[2]+uvw4[2];
                        //printf("start=%i %i %i %i %f %f %f %f\n",start,nw,i2,j2,uvw1[0],uvw2[1],uvw3[2],uvw4[2]);
                    }
                }
            }            
            if (summode == 0){
                *(pwing->uw+i+j*m)=uvw[0];
                *(pwing->vw+i+j*m)=uvw[1];
                *(pwing->ww+i+j*m)=uvw[2];
            }else{
                *(pwing->uw+i+j*m)+=uvw[0];
                *(pwing->vw+i+j*m)+=uvw[1];
                *(pwing->ww+i+j*m)+=uvw[2];
            }
        }
    }
}

void wakeinf(struct liftsurf *pwing, struct liftsurf *pflap, int nw, int it, int summode)
{
    /* Calculate flow induced by the vortex rings of wake nw in pflap on pwing */
    double uvw1[3],uvw2[3],uvw3[3],uvw4[3],uvw[3],xc,yc,zc;
    int i, i2, j2, mp1, start, start2, ppp;
    
    mp1=pflap->nshed+pflap->nwakes;
    /* Find start point in pflap wake*/
    start=0;
    start2=0;
    if (nw > 0){
        for (i=0;i<nw;i++){
            start+=*(pflap->wakelengths+i);
            start2+=*(pflap->wakelengths+i)+1;
        }
    }

    /* Cycle through all the control points of pwing and all the vortex segments of pflap */
    for (i=0; i <pwing->nface; i++) {
        uvw[0]= 0.0;
        uvw[1]= 0.0;
        uvw[2]= 0.0;
        /* Control point */
        xc=*(pwing->control+i);
        yc=*(pwing->control+i+pwing->nface);
        zc=*(pwing->control+i+2*pwing->nface);
        //printf("%i %f %f %f\n",i,xc,yc,zc);
        for (i2=0; i2 < *(pflap->wakelengths+nw); i2++) {
            for (j2=0; j2 < it; j2++) {
                ppp=start+i2+j2*pflap->nshed;
                /* Influence of vortex segment 1 */
                vortex(uvw1,xc,yc,zc,
                        *(pflap->xw+start2+i2+1+j2*mp1),*(pflap->yw+start2+i2+1+j2*mp1),*(pflap->zw+start2+i2+1+j2*mp1),
                        *(pflap->xw+start2+i2+1+(j2+1)*mp1),*(pflap->yw+start2+i2+1+(j2+1)*mp1),*(pflap->zw+start2+i2+1+(j2+1)*mp1),
                        *(pflap->gw+ppp));
                /* Influence of vortex segment 2 */
                vortex(uvw2,xc,yc,zc,
                        *(pflap->xw+start2+i2+1+(j2+1)*mp1),*(pflap->yw+start2+i2+1+(j2+1)*mp1),*(pflap->zw+start2+i2+1+(j2+1)*mp1),
                        *(pflap->xw+start2+i2+(j2+1)*mp1),*(pflap->yw+start2+i2+(j2+1)*mp1),*(pflap->zw+start2+i2+(j2+1)*mp1),
                        *(pflap->gw+ppp));
                /* Influence of vortex segment 3 */
                vortex(uvw3,xc,yc,zc,
                         *(pflap->xw+start2+i2+(j2+1)*mp1),*(pflap->yw+start2+i2+(j2+1)*mp1),*(pflap->zw+start2+i2+(j2+1)*mp1),
                         *(pflap->xw+start2+i2+j2*mp1),*(pflap->yw+start2+i2+j2*mp1),*(pflap->zw+start2+i2+j2*mp1),
                         *(pflap->gw+ppp));
                /* Influence of vortex segment 4 */
                vortex(uvw4,xc,yc,zc,
                        *(pflap->xw+start2+i2+j2*mp1),*(pflap->yw+start2+i2+j2*mp1),*(pflap->zw+start2+i2+j2*mp1),
                        *(pflap->xw+start2+i2+1+j2*mp1),*(pflap->yw+start2+i2+1+j2*mp1),*(pflap->zw+start2+i2+1+j2*mp1),
                        *(pflap->gw+ppp));
                /* Addd the four influcences */
                uvw[0]=uvw[0]+uvw1[0]+uvw2[0]+uvw3[0]+uvw4[0];
                uvw[1]=uvw[1]+uvw1[1]+uvw2[1]+uvw3[1]+uvw4[1];
                uvw[2]=uvw[2]+uvw1[2]+uvw2[2]+uvw3[2]+uvw4[2];
            }
        }
        /* Choose whether to sum effects of different wake surfaces or not */
        /* This assumes that the effect of one wake (e.g. wing) has already been calculated */
        if (summode == 0){
            *(pwing->uvw+i)=uvw[0];
            *(pwing->uvw+i+pwing->nface)=uvw[1];
            *(pwing->uvw+i+2*pwing->nface)=uvw[2];
        }else{
            *(pwing->uvw+i)+=uvw[0];
            *(pwing->uvw+i+pwing->nface)+=uvw[1];
            *(pwing->uvw+i+2*pwing->nface)+=uvw[2];
        }
    }
}

void calcwakeinf(struct VLMData *data, int it)
{
    /* Calculate the influence of every wake on every liftsurf */
    /* Treat the vtail and elevator as uncoupled */
    int i;
    
    /* Influence on the wing */
    /* Influence of the wing wake */
    for (i=0;i<(data->wing).nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(&(data->wing),&(data->wing),i,it,i);
    }
    /* Influence of the flap wake */
    for (i=0;i<(data->flap).nwakes;i++){
        wakeinf(&(data->wing),&(data->flap),i,it,1);
    }
    /* Influence of the aileron wake */
    for (i=0;i<(data->aileron).nwakes;i++){
        wakeinf(&(data->wing),&(data->aileron),i,it,1);
    }
    /* Influence of the horizontal tail wake */
    for (i=0;i<(data->htail).nwakes;i++){
        wakeinf(&(data->wing),&(data->htail),i,it,1);
    }
    /* Influence of the elevator wake */
    for (i=0;i<(data->elevator).nwakes;i++){
        wakeinf(&(data->wing),&(data->elevator),i,it,1);
    } 
//     /* Influence of the vertical tail wake */
//     for (i=0;i<(data->vtail).nwakes;i++){
//         wakeinf(&(data->wing),&(data->vtail),i,it,1);
//     }
//     /* Influence of the rudder wake */
//     for (i=0;i<(data->rudder).nwakes;i++){
//         wakeinf(&(data->wing),&(data->rudder),i,it,1);
//     }     
    /* Influence on the flap */
    /* Influence of the wing wake */
    for (i=0;i<(data->wing).nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(&(data->flap),&(data->wing),i,it,i);
    }
    /* Influence of the flap wake */
    for (i=0;i<(data->flap).nwakes;i++){
        wakeinf(&(data->flap),&(data->flap),i,it,1);
    }
    /* Influence of the aileron wake */
    for (i=0;i<(data->aileron).nwakes;i++){
        wakeinf(&(data->flap),&(data->aileron),i,it,1);
    } 
    /* Influence of the horizontal tail wake */
    for (i=0;i<(data->htail).nwakes;i++){
        wakeinf(&(data->flap),&(data->htail),i,it,1);
    } 
    /* Influence of the elevator wake */
    for (i=0;i<(data->elevator).nwakes;i++){
        wakeinf(&(data->flap),&(data->elevator),i,it,1);
    } 
//     /* Influence of the vertical tail wake */
//     for (i=0;i<(data->vtail).nwakes;i++){
//         wakeinf(&(data->flap),&(data->vtail),i,it,1);
//     }
//     /* Influence of the rudder wake */
//     for (i=0;i<(data->rudder).nwakes;i++){
//         wakeinf(&(data->flap),&(data->rudder),i,it,1);
//     }     
    /* Influence on the aileron */
    /* Influence of the wing wake */
    for (i=0;i<(data->wing).nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(&(data->aileron),&(data->wing),i,it,i);
    }
    /* Influence of the flap wake */
    for (i=0;i<(data->flap).nwakes;i++){
        wakeinf(&(data->aileron),&(data->flap),i,it,1);
    }
    /* Influence of the aileron wake */
    for (i=0;i<(data->aileron).nwakes;i++){
        wakeinf(&(data->aileron),&(data->aileron),i,it,1);
    }  
    /* Influence of the horizontal tail wake */
    for (i=0;i<(data->htail).nwakes;i++){
        wakeinf(&(data->aileron),&(data->htail),i,it,1);
    }
    /* Influence of the elevator wake */
    for (i=0;i<(data->elevator).nwakes;i++){
        wakeinf(&(data->aileron),&(data->elevator),i,it,1);
    }  
//     /* Influence of the vertical tail wake */
//     for (i=0;i<(data->vtail).nwakes;i++){
//         wakeinf(&(data->aileron),&(data->vtail),i,it,1);
//     }
//     /* Influence of the rudder wake */
//     for (i=0;i<(data->rudder).nwakes;i++){
//         wakeinf(&(data->aileron),&(data->rudder),i,it,1);
//     } 
    /* Influence on the horizontal tail */
    /* Influence of the wing wake */
    for (i=0;i<(data->wing).nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(&(data->htail),&(data->wing),i,it,i);
    }
    /* Influence of the flap wake */
    for (i=0;i<(data->flap).nwakes;i++){
        wakeinf(&(data->htail),&(data->flap),i,it,1);
    }
    /* Influence of the aileron wake */
    for (i=0;i<(data->aileron).nwakes;i++){
        wakeinf(&(data->htail),&(data->aileron),i,it,1);
    }  
    /* Influence of the horizontal tail wake */
    for (i=0;i<(data->htail).nwakes;i++){
        wakeinf(&(data->htail),&(data->htail),i,it,1);
    }
    /* Influence of the elevator wake */
    for (i=0;i<(data->elevator).nwakes;i++){
        wakeinf(&(data->htail),&(data->elevator),i,it,1);
    } 
//     /* Influence of the vertical tail wake */
//     for (i=0;i<(data->vtail).nwakes;i++){
//         wakeinf(&(data->htail),&(data->vtail),i,it,1);
//     }
//     /* Influence of the rudder wake */
//     for (i=0;i<(data->rudder).nwakes;i++){
//         wakeinf(&(data->htail),&(data->rudder),i,it,1);
//     }    
    /* Influence on the elevator tail */
    /* Influence of the wing wake */
    for (i=0;i<(data->wing).nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(&(data->elevator),&(data->wing),i,it,i);
    }
    /* Influence of the flap wake */
    for (i=0;i<(data->flap).nwakes;i++){
        wakeinf(&(data->elevator),&(data->flap),i,it,1);
    }
    /* Influence of the aileron wake */
    for (i=0;i<(data->aileron).nwakes;i++){
        wakeinf(&(data->elevator),&(data->aileron),i,it,1);
    }  
    /* Influence of the horizontal tail wake */
    for (i=0;i<(data->htail).nwakes;i++){
        wakeinf(&(data->elevator),&(data->htail),i,it,1);
    }
    /* Influence of the elevator wake */
    for (i=0;i<(data->elevator).nwakes;i++){
        wakeinf(&(data->elevator),&(data->elevator),i,it,1);
    }     
//     /* Influence of the vertical tail wake */
//     for (i=0;i<(data->vtail).nwakes;i++){
//         wakeinf(&(data->elevator),&(data->vtail),i,it,1);
//     }
//     /* Influence of the rudder wake */
//     for (i=0;i<(data->rudder).nwakes;i++){
//         wakeinf(&(data->elevator),&(data->rudder),i,it,1);
//     } 
    /* Influence on the vertical tail */
//     /* Influence of the wing wake */
//     for (i=0;i<(data->wing).nwakes;i++){
//         /* We set summode to i so summode=0 the first time we calculate an influence */
//         wakeinf(&(data->vtail),&(data->wing),i,it,i);
//     }
//     /* Influence of the flap wake */
//     for (i=0;i<(data->flap).nwakes;i++){
//         wakeinf(&(data->vtail),&(data->flap),i,it,1);
//     }
//     /* Influence of the aileron wake */
//     for (i=0;i<(data->aileron).nwakes;i++){
//         wakeinf(&(data->vtail),&(data->aileron),i,it,1);
//     }  
//     /* Influence of the horizontal tail wake */
//     for (i=0;i<(data->htail).nwakes;i++){
//         wakeinf(&(data->vtail),&(data->htail),i,it,1);
//     }
//     /* Influence of the elevator wake */
//     for (i=0;i<(data->elevator).nwakes;i++){
//         wakeinf(&(data->vtail),&(data->elevator),i,it,1);
//     } 
    /* Influence of the vertical tail wake */
    for (i=0;i<(data->vtail).nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(&(data->vtail),&(data->vtail),i,it,i);
    }
    /* Influence of the rudder wake */
    for (i=0;i<(data->rudder).nwakes;i++){
        wakeinf(&(data->vtail),&(data->rudder),i,it,1);
    }    
    /* Influence on the rudder */
//     /* Influence of the wing wake */
//     for (i=0;i<(data->wing).nwakes;i++){
//         /* We set summode to i so summode=0 the first time we calculate an influence */
//         wakeinf(&(data->rudder),&(data->wing),i,it,i);
//     }
//     /* Influence of the flap wake */
//     for (i=0;i<(data->flap).nwakes;i++){
//         wakeinf(&(data->rudder),&(data->flap),i,it,1);
//     }
//     /* Influence of the aileron wake */
//     for (i=0;i<(data->aileron).nwakes;i++){
//         wakeinf(&(data->rudder),&(data->aileron),i,it,1);
//     }  
//     /* Influence of the horizontal tail wake */
//     for (i=0;i<(data->htail).nwakes;i++){
//         wakeinf(&(data->rudder),&(data->htail),i,it,1);
//     }
//     /* Influence of the elevator wake */
//     for (i=0;i<(data->elevator).nwakes;i++){
//         wakeinf(&(data->rudder),&(data->elevator),i,it,1);
//     } 
    /* Influence of the vertical tail wake */
    for (i=0;i<(data->vtail).nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(&(data->rudder),&(data->vtail),i,it,i);
    }
    /* Influence of the rudder wake */
    for (i=0;i<(data->rudder).nwakes;i++){
        wakeinf(&(data->rudder),&(data->rudder),i,it,1);
    }
}

void latestwakeinf(struct liftsurf *pwing, struct liftsurf *pflap, int nw)
{
    /* Calculate flow induced by the latest vortex segments of wake nw in pflap on pwing */
    double uvw1[3],xc,yc,zc,uvw;
    int i, i2, m, mp1, start, start2;
    
    m=pflap->nshed;
    mp1=pflap->nshed+pflap->nwakes;
    /* Find start point in pflap wake*/
    start=0;
    start2=0;
    if (nw > 0){
        for (i=0;i<nw;i++){
            start+=*(pflap->wakelengths+i);
            start2+=*(pflap->wakelengths+i)+1;
        }
    }

    /* Cycle through all the control points of pwing and all the vortex segments of pflap */
    for (i=0; i <pwing->nface; i++) {
        uvw = 0.0;
        /* Control point */
        xc=*(pwing->control+i);
        yc=*(pwing->control+i+pwing->nface);
        zc=*(pwing->control+i+2*pwing->nface);
        //printf("%i %f %f %f\n",i,xc,yc,zc);
        for (i2=0; i2 < *(pflap->wakelengths+nw); i2++) {
            /* Influence of vortex segment 1 */
            vortex(uvw1,xc,yc,zc,
                    *(pflap->xw+start2+i2),*(pflap->yw+start2+i2),*(pflap->zw+start2+i2),
                    *(pflap->xw+start2+i2+1),*(pflap->yw+start2+i2+1),*(pflap->zw+start2+i2+1),
                    *(pflap->gw+start+i2));
            uvw+=uvw1[2];
        }
        *(pwing->wind+i)-=uvw;
    }
}

void addlastwind(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, int it, int narg)
{
    /* Add influence of latest unsteady wake element */
    int i;
    
    /* Influence on the wing */
    /* Influence of the wing wake */
    for (i=0;i<pwing->nwakes;i++){
        latestwakeinf(pwing,pwing,i);
    }
    /* Influence of the flap wake */
    for (i=0;i<pflap->nwakes;i++){
        latestwakeinf(pwing,pflap,i);
    }
    /* Influence on the flap */
    /* Influence of the wing wake */
    for (i=0;i<pwing->nwakes;i++){
        latestwakeinf(pflap,pwing,i);
    }
    /* Influence of the flap wake */
    for (i=0;i<pflap->nwakes;i++){
        latestwakeinf(pflap,pflap,i);
    }
    if (narg > 2){
        /* Influence on the wing of the aileron wake */
        for (i=0;i<paileron->nwakes;i++){
            latestwakeinf(pflap,paileron,i);
        }
        /* Influence on the flap of the aileron wake */
        for (i=0;i<paileron->nwakes;i++){
            latestwakeinf(pwing,paileron,i);
        }
        /* Influence on the aileron */
        /* Influence of the wing wake */
        for (i=0;i<pwing->nwakes;i++){
            latestwakeinf(paileron,pwing,i);
        }
        /* Influence of the flap wake */
        for (i=0;i<pflap->nwakes;i++){
            latestwakeinf(paileron,pflap,i);
        }
        /* Influence of the aileron wake */
        for (i=0;i<paileron->nwakes;i++){
            latestwakeinf(paileron,paileron,i);
        }
    }
}

