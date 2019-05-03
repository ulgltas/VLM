#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <string.h>

struct liftsurf{
   int    nface; /* number of panels */
   int    nvert; /* number of vertices */
   int *faces; /* panel definitions */
   int *neighbours; /* neighbouring panels to each panel arranged in the order left, right, upstream, downstream */
   double *vertices; /* panel vertices */
   int nshed; /* Number of trailing edge panels that shed wake */
   int *shedding; /* Trailing edge panel that sheds wake (0 or 1) */
   double *control; /* control point positions */
   double *vortex; /* Vortex panel vertices */
   double *tangx; /* Tangential vectors in x direction */
   double *tangy; /* Tangential vectors in y direction */
   double *normal; /* Normal vectors */
   double *nsurf; /* Normal surfaces of panels */
   double *dxy; /* Vortex segment lengths in x and y directions */
   double *Deltap; /* Pressure difference across each panel */
   double *Deltad; /* Induced drag contribution of each panel */
   double *gamma; /* Vortex strength on each panel */
   double *wind; /* Induced airspeed on each panel */
   double *uvw; /* Local airspeed induced by wakes */
   double *aeroforce; /* Aerodynamic force components in x, y and z directions */
   int nwakes; /* Number of contiguous wakes shed from the lifting surface */
   int *wakeinds; /* Indices of panels that shed wake */
   int *wakelengths; /* Number of wake panels shed at each time step in each wake */
   int *wakecorrespx; /* Correspondence between bound vortex panels and wake panel vertices */
   int *wakecorrespy; /* Correspondence between bound vortex panels and wake panel vertices */
   int *wakecorrespz; /* Correspondence between bound vortex panels and wake panel vertices */
   int *wakecorrespwake; /* Correspondence between wake panel vertices in different wakes */
   double *xw; /* x-coordinates of wake vertices */
   double *yw; /* y-coordinates of wake vertices */
   double *zw; /* z-coordinates of wake vertices */
   double *uw; /* velocity in x-direction of wake vertices */
   double *vw; /* velocity in y-direction of wake vertices */
   double *ww; /* velocity in z-direction of wake vertices */
   double *gw; /* Vortex strength on each wake panel */
};

int compare_function(const void *a,const void *b) {
    /* Comparison function for use with qsort */
    double *x = (double *) a;
    double *y = (double *) b;
    if (*x < *y) return -1;
    else if (*x > *y) return 1; return 0;
}

int checkdoubles(double *ypos, int m)
{
    /* Check if double entries exist in vectors */
    int n,i,j,nj;
    
    n=0;
    for (i=0;i<m;i++){
        for (j=i+1;j<m;j++){
            if (fabs(*(ypos+i)- *(ypos+j)) < 0.01){
                n++;
            }
        }
    }
    return n;
}

int finddoubles(double *ypos, int m)
{
    /* Find double entries in vectors */
    /* If there are double entries, they are pushed towards the end of the vector and overwritten */
    int n,i,j,nj;
    
    n=0;
    for (i=0;i<m;i++){
        nj=0;
        for (j=i+1;j<m;j++){
            if (fabs(*(ypos+i)- *(ypos+j)) < 0.01){
                nj=j;
                n++;
                for (j=nj+1;j<m;j++){
                    *(ypos+j-1)=*(ypos+j); /* Move all values from double to end of vector one element earlier */
                }
                *(ypos+m-1)=-1.0*n; /* Overwrite last element of vector with a non-double value */
                j--;
            }
        }
    }
    return m-n;
}

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

void  createypline(double *ypos,int nypos,double *ypline,int np1)
{
    /* Create regularly spaced spanwise panels on wing that 
     * take into account trapezoidal sections, flaps and ailerons */
    
    int n,nyposm1,i,j,*npanels,donepanels,diffpanels,nadj,iko;
    double *ydiff, *npanelsold, dy, maxdiffold,ddd;
    
    n=np1-1;
    nyposm1=nypos-1;

    ydiff = (double *)malloc(sizeof(double)*nyposm1);
    npanelsold = (double *)malloc(sizeof(double)*nyposm1);
    npanels = (int *)malloc(sizeof(int)*nyposm1);
    
    /* Calculate size of wing sections */
    for (i=0;i<nyposm1;i++){
        *(ydiff+i)=*(ypos+i+1)-*(ypos+i);
    }
    /* Calculate representative panel length */
    dy=*(ypos+nyposm1)/n;
    /* Find number of panels in each section */
    for (i=0;i<nyposm1;i++){
        *(npanelsold+i)=*(ydiff+i)/dy;
        if (*(npanelsold+i) < 1.0)
            *(npanels+i)=1;
        else
            *(npanels+i)=round(*(npanelsold+i));
    }
    /* Calculate if the total number of panels is equal to the desired number of panels */
    diffpanels=-n;
    for (i=0;i<nyposm1;i++){
        diffpanels+=*(npanels+i);
    }
    /* If the total number of panels is greater than the desired number of panels */
    if (diffpanels < 0){
        nadj=0;
        while (nadj < -diffpanels){
            /* Find maximum difference between rounded and raw number of panels */
            maxdiffold=-10;
            iko=-10;
            for (i=0;i<nyposm1;i++){
                ddd=*(npanelsold+i)-*(npanels+i);
                if (ddd > maxdiffold){
                    maxdiffold=ddd;
                    iko=i;
                }
            }
            /* Use ceiling instead of round on the maximum difference */
            if (*(npanelsold+iko) > 1.0){
                *(npanels+iko)=ceil(*(npanelsold+iko));
                nadj++;
            }else{
                *(npanelsold+iko)=1.000000001;
            }        
        }
    }
    /* If the total number of panels is less than the desired number of panels */
    if (diffpanels > 0){
        nadj=0;
        while (nadj < diffpanels){
            /* Find maximum difference between rounded and raw number of panels */
            maxdiffold=-10;
            iko=-10;
            for (i=0;i<nyposm1;i++){
                ddd=*(npanels+i)- *(npanelsold+i);
                if (ddd > maxdiffold){
                    maxdiffold=ddd;
                    iko=i;
                }
            }
            /* Use floor instead of round on the maximum difference */
            if (*(npanelsold+iko) > 1.0){
                *(npanels+iko)=floor(*(npanelsold+iko));
                nadj++;
            }else{
                *(npanelsold+iko)=1.000000001;
            }
        }
    }
    /* Calculate ypanel vector with all the panel endpoints */
    donepanels=0;
    for (i=0;i<nyposm1;i++){
        for (j=0;j<*(npanels+i);j++){
            *(ypline+donepanels+j)=*(ydiff+i)/(*(npanels+i))*j+*(ypos+i);
        }
        donepanels+=*(npanels+i);
    }
    *(ypline+n)=*(ypos+nyposm1);
    
    free(ydiff);
    free(npanelsold);
    free(npanels);
}


int countfaces(int *ijflap, int nflap, int mp1, int np1)
{
    /* Count the number of panels in lifting surface
     * ijflap contains the indices from the complete wing grid that 
     * lie in the current lifting surface. 
     * nflap is the number of vertices in the current lifting surface. */
    int nfaces, i, i1, j1, k, condij[3];
    
    nfaces=0; /* number of panels */
    for (i=0;i<nflap;i++){
        i1=*(ijflap+i); /* chordwise index */
        j1=*(ijflap+i+mp1*np1); /* Spanwise index */
        condij[0]=0;condij[1]=0;condij[2]=0;
        for (k=i+1;k<nflap;k++){
            /* Check if we can create a complete rectangle from point i1,j1 */
            if (*(ijflap+k)==i1+1 && *(ijflap+k+mp1*np1)==j1){
                condij[0]=k;
            }
            if (*(ijflap+k)==i1+1 && *(ijflap+k+mp1*np1)==j1+1){
                condij[1]=k;
            }
            if (*(ijflap+k)==i1 && *(ijflap+k+mp1*np1)==j1+1){
                condij[2]=k;
            }
        }
        /* If we can create the rectangle, increment nfaces */
        if (condij[0] != 0 && condij[1] != 0 && condij[2] != 0){
            nfaces++;
        }
    }
    return nfaces;
}

void arrangefaces(int *ijflap, int mp1, int np1, struct liftsurf *pflap)
{
    /* Create the connections between the vertices of the current lifting surface
     * that form the panels
     * pflap is the lifting surface structure
     * ijflap contains the indices from the complete wing grid that 
     * lie in the current lifting surface.*/
    
    int nfaces, nfaced2, nvertd2, i, i1, j1, k, condij[3],nte;
    
    nfaced2=pflap->nface/2;
    nvertd2=pflap->nvert/2;
    nfaces=0; /* number of panels */
    pflap->nshed=0; /* Initialize the number of shedding panels to 0 */
    for (i=0;i<nvertd2;i++){
        i1=*(ijflap+i); /* chordwise index */
        j1=*(ijflap+i+mp1*np1); /* Spanwise index */
        condij[0]=0;condij[1]=0;condij[2]=0;
        for (k=i+1;k<nvertd2;k++){
           /* Check if we can create a complete rectangle from point i1,j1
            * and store in condij the indices of each panel vertex */
           if (*(ijflap+k)==i1+1 && *(ijflap+k+mp1*np1)==j1){
                condij[0]=k;
            }
            if (*(ijflap+k)==i1+1 && *(ijflap+k+mp1*np1)==j1+1){
                condij[1]=k;
            }
            if (*(ijflap+k)==i1 && *(ijflap+k+mp1*np1)==j1+1){
                condij[2]=k;
            }
        }
        /* If we can create the rectangle, increment nfaces 
         * and store the elements of condij in pflap->faces */
        if (condij[0] != 0 && condij[1] != 0 && condij[2] != 0){
            nfaces++;  
            /* Right lifting surface */
            *(pflap->faces+nfaces-1)=i;
            *(pflap->faces+nfaces-1+pflap->nface)=condij[0];
            *(pflap->faces+nfaces-1+2*pflap->nface)=condij[1];
            *(pflap->faces+nfaces-1+3*pflap->nface)=condij[2];
            *(pflap->shedding+nfaces-1)=0;
            *(pflap->shedding+nfaces-1)=0;
            /* Left lifting surface */
            *(pflap->faces+nfaces-1+nfaced2)=i+nvertd2;
            *(pflap->faces+nfaces-1+nfaced2+pflap->nface)=condij[0]+nvertd2;
            *(pflap->faces+nfaces-1+nfaced2+2*pflap->nface)=condij[1]+nvertd2;
            *(pflap->faces+nfaces-1+nfaced2+3*pflap->nface)=condij[2]+nvertd2;
            *(pflap->shedding+nfaced2+nfaces-1)=0;
            if (i1 == mp1-2){ /* Check if this a trailing edge panel that sheds wake */
                /* Right lifting surface */
                *(pflap->shedding+nfaces-1)=1;
                pflap->nshed++;
                /* Left lifting surface */
                *(pflap->shedding+nfaces-1+nfaced2)=1;
                pflap->nshed++;
            }
        }
    }
}

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

void readarfsize(char *WngArf,int ArfXUNumber[],int ArfXLNumber[])
{
    /* Read size of data in airfoil file */
    FILE *fp2;
    int cond,npoint;
    char line[110],code[7];
    char* arfname=malloc(strlen(WngArf) + strlen(".arf") + 2);
    
    /* Append .arf at the end of the airfoil name to obtain the filename */
    sscanf(WngArf,"%s",arfname);
    strcat(arfname, ".arf");
    
    /* Open airfoil file */
    fp2 = fopen(arfname,"r");
    if(fp2 == NULL) {
        fprintf(stderr,"Error:  Could not open file %s\n",arfname);
        exit(1);
    }
    
    /* Scan the file until lines GMD401 and GMD402 to read the number of points on the two surfaces */
    cond=0;
    while (cond == 0){
        fgets(line, 110, fp2);
        /* Number of points on upper surface */
        if ( strncmp("GMD401",line,6) == 0 ){
            sscanf(line,"%s %i",code,ArfXUNumber);
        }
        /* Number of points on lower surface */
        if ( strncmp("GMD402",line,6) == 0 ){
            sscanf(line,"%s %i",code,ArfXLNumber);
            cond=1;
        }
    }
    /* Close airfoil file */
    fclose(fp2);
}

void readarf(char *WngArf,double *ArfXU, double *ArfYU, double *ArfXL,double *ArfYL,int ArfXUNumber[],int ArfXLNumber[])
{
    /* Read airfoil file data */
    FILE *fp2;
    int cond,npoint;
    char line[110],code[7];
    char* arfname=malloc(strlen(WngArf) + strlen(".arf") + 2);
    
    /* Append .arf at the end of the airfoil name to obtain the filename */
    sscanf(WngArf,"%s",arfname);
    strcat(arfname, ".arf");
    
    /* Open airfoil file */
    fp2 = fopen(arfname,"r");
    if(fp2 == NULL) {
        fprintf(stderr,"Error:  Could not open file %s\n",arfname);
        exit(1);
    }
    
    /* Scan the file until lines GMD15 to read the upper and lower surface data */
    fgets(line, 110, fp2);
    cond=0;
    npoint=-1; /* Set up counter to count how many datapoints have been read */
    while (cond == 0){
        fgets(line, 110, fp2);
        if ( strncmp("GM15",line,4) == 0 ){
            npoint++;
            /* If there are more points on the upper surface than the lower surface */
            if (ArfXUNumber[0] > ArfXLNumber[0]){
                if (npoint < ArfXLNumber[0]){
                    sscanf(line,"%s %lf %lf %lf %lf",code,(ArfXU+npoint),(ArfYU+npoint),(ArfXL+npoint),(ArfYL+npoint));
                }
                if (npoint >= ArfXLNumber[0] && npoint < ArfXUNumber[0]){
                    sscanf(line,"%s %lf %lf",code,(ArfXU+npoint),(ArfYU+npoint));
                }
           /* If there are more points on the lower surface than the upper surface */
           }else if (ArfXUNumber[0] < ArfXLNumber[0]){
                if (npoint < ArfXUNumber[0]){
                    sscanf(line,"%s %lf %lf %lf %lf",code,(ArfXU+npoint),(ArfYU+npoint),(ArfXL+npoint),(ArfYL+npoint));
                }
                if (npoint >= ArfXUNumber[0] && npoint < ArfXLNumber[0]){
                    sscanf(line,"%s %lf %lf",code,(ArfXL+npoint),(ArfYL+npoint));
                }
           /* If the number of points on the upper and lower surface is the same */
            }else{
                if (npoint < ArfXUNumber[0]){
                    sscanf(line,"%s %lf %lf %lf %lf",code,(ArfXU+npoint),(ArfYU+npoint),(ArfXL+npoint),(ArfYL+npoint));
                }
            }
        }
        /* Stop condition for the while loop */
        if ( strncmp("GMD6",line,4) == 0 ){
            cond=1;
        }
    }
    /* Close airfoil file */
    fclose(fp2);
}

void interparf(double *ycamber, double *ArfXU, double *ArfYU, double *ArfXL, double *ArfYL, int ArfXUNumber[], int ArfXLNumber[], double *xpline, int mp1)
{
    /* Use interpolation to obtain airfoil camber data at the desired x-coordinates */
    int i,j;
    double slope;
    
    /* Carry out linear interpolation to find ycamber at every element of xpline */
    for (i=0;i<mp1;i++){
        j=0;
        /* Upper surface interpolation */
        while (*(ArfXU+j)/100.0 <= *(xpline+i) && j <= ArfXUNumber[0]-1){
            j++;
        }
        slope=(*(ArfYU+j)- *(ArfYU+j-1))/(*(ArfXU+j)- *(ArfXU+j-1));
        /* Store the upper surface y-coordinate in ycamber */
        *(ycamber+i)=(slope*(*(xpline+i)- *(ArfXU+j-1)/100)+ *(ArfYU+j-1)/100.0);
        
        /* Lower surface interpolation */
        j=0;
        while (*(ArfXL+j)/100.0 <= *(xpline+i) && j <= ArfXLNumber[0]-1){
            j++;
        }
        slope=(*(ArfYL+j)- *(ArfYL+j-1))/(*(ArfXL+j)- *(ArfXL+j-1));
        /* Add the lower surface y-coordinate in ycamber and divide by two */
        *(ycamber+i)=(*(ycamber+i)+(slope*(*(xpline+i)- *(ArfXL+j-1)/100)+ *(ArfYL+j-1)/100.0))/2.0;
    }
}

void nacafourfivedigit(double *xpline, double *ycamber, int mp1, char WngArf[])
{
    /* Calculate the camber line of a 4- or 5-digit NACA series airfoil */
    int i;
    double xhere,m,p,k1,meanline;
    char *end;
    long number = strtol(WngArf, &end, 10); /* Get airfoil name */
    
    /* Check if it's a 4- or 5-digit series */
    if (number > 9999){
        /* NACA 5-digit series */
        meanline=floor(((double)number)/100.0);
        /* Choose a mean line */
        if (meanline == 210.0){
            p=0.05;
            m=0.0580;
            k1=361.4;
        }else if (meanline == 220.0){
            p=0.1;
            m=0.1260;
            k1=51.64;
        }else if (meanline == 230.0){
            p=0.15;
            m=0.2025;
            k1=15.957;
        }else if (meanline == 240.0){
            p=0.2;
            m=0.29;
            k1=6.643;
        }else if (meanline == 250.0){
            p=0.25;
            m=0.391;
            k1=3.23;
        }
        /* Calculate camber line */
        for (i=0;i<mp1;i++){
            xhere=*(xpline+i);
            if (xhere <= m){
                *(ycamber+i)=1/6.0*k1*(pow(xhere,3)-3*m*pow(xhere,2)+m*m*(3-m)*xhere);
            }else{
                *(ycamber+i)=1/6.0*k1*m*m*m*(1-xhere);
            }
        }
    }else{
        /* NACA 4-digit series */
        /* Choose maximum camber m and maximum camber position p */
        m=floor(((double)number)/1000.0)/100;
        p=floor((((double)number)-m*100000.0)/100.0)/10.0;
        /* Calculate camber line */
        for (i=0;i<mp1;i++){
            if (m == 0){
                *(ycamber+i)=0;
            }else{
                xhere=*(xpline+i);
                if (xhere <= p){
                    *(ycamber+i)=m/(p*p)*(2.0*p*xhere-xhere*xhere);
                }else{
                    *(ycamber+i)=m/(pow((1.0-p),2))*((1-2.0*p)+2.0*p*xhere-xhere*xhere);
                }
            }
        }
    }
}

void treatarfwgl(char line[], double *xpline, double *ycamber, int mp1)
{
    /* Carries out the treatment of the airfoil on the winglet */
    int ArfXUNumber[1],ArfXLNumber[1],j;
    double *ArfXU, *ArfYU, *ArfXL, *ArfYL;
    char code[7], *WglArf, *nacanumber;
    char *token = NULL;
               
    token = strtok(line, "\t,\n");
    WglArf = strtok(NULL, "\t,\n");
    
    /* Inboard airfoil */
    if ( strncmp("NACA",WglArf,4) == 0 ){
        token = strtok(WglArf, " ,\n");
        nacanumber=strtok(NULL, " ,\n");
        nacafourfivedigit(xpline,ycamber,mp1,nacanumber);
    }else{
         readarfsize(WglArf,ArfXUNumber,ArfXLNumber); /* Read airfoil data size */
        /* Set up vectors to store airfoil data */
        ArfXU = (double *)malloc(sizeof(double)*ArfXUNumber[0]);
        ArfYU = (double *)malloc(sizeof(double)*ArfXUNumber[0]);
        ArfXL = (double *)malloc(sizeof(double)*ArfXLNumber[0]);
        ArfYL = (double *)malloc(sizeof(double)*ArfXLNumber[0]);
        /* Read airfoil data */
        readarf(WglArf,ArfXU,ArfYU,ArfXL,ArfYL,ArfXUNumber,ArfXLNumber);
        /* Interpolate airfoil at the desired coordinates xpline */
        interparf(ycamber,ArfXU,ArfYU,ArfXL,ArfYL,ArfXUNumber,ArfXLNumber,xpline,mp1);
    }  
}

void treatarf(char line[], double *xpline, double *ycamber, double *ycamberall, int mp1, int iTS)
{
    /* Carries out the treatment of the airfoils on each end of each trapezoidal section */
    int ArfXUNumber[1],ArfXLNumber[1],j;
    double *ArfXU, *ArfYU, *ArfXL, *ArfYL;
    char code[7], *WngTSRtArf, *WngTSTpArf, *nacanumber;
    char *token = NULL;
               
    token = strtok(line, "\t,\n");
    WngTSRtArf = strtok(NULL, "\t,\n");
    WngTSTpArf = strtok(NULL, "\t,\n");
    
    /* Inboard airfoil */
    if ( strncmp("NACA",WngTSRtArf,4) == 0 ){
        token = strtok(WngTSRtArf, " ,\n");
        nacanumber=strtok(NULL, " ,\n");
        nacafourfivedigit(xpline,ycamber,mp1,nacanumber);
        for (j=0;j<mp1;j++){
            *(ycamberall+j+2*iTS*mp1)=*(ycamber+j);
        }
    }else{
         readarfsize(WngTSRtArf,ArfXUNumber,ArfXLNumber); /* Read airfoil data size */
        /* Set up vectors to store airfoil data */
        ArfXU = (double *)malloc(sizeof(double)*ArfXUNumber[0]);
        ArfYU = (double *)malloc(sizeof(double)*ArfXUNumber[0]);
        ArfXL = (double *)malloc(sizeof(double)*ArfXLNumber[0]);
        ArfYL = (double *)malloc(sizeof(double)*ArfXLNumber[0]);
        /* Read airfoil data */
        readarf(WngTSRtArf,ArfXU,ArfYU,ArfXL,ArfYL,ArfXUNumber,ArfXLNumber);
        /* Interpolate airfoil at the desired coordinates xpline */
        interparf(ycamber,ArfXU,ArfYU,ArfXL,ArfYL,ArfXUNumber,ArfXLNumber,xpline,mp1);
        for (j=0;j<mp1;j++){
            *(ycamberall+j+2*iTS*mp1)=*(ycamber+j);
        }
    }
    
     /* Outboard airfoil */
    if ( strncmp("NACA",WngTSTpArf,4) == 0 ){
        token = strtok(WngTSTpArf, " ,\n");
        nacanumber=strtok(NULL, " ,\n");
        nacafourfivedigit(xpline,ycamber,mp1,nacanumber);
        for (j=0;j<mp1;j++){
            *(ycamberall+j+(2*iTS+1)*mp1)=*(ycamber+j);
        }
    }else{
         readarfsize(WngTSTpArf,ArfXUNumber,ArfXLNumber); /* Read airfoil data size */
        /* Set up vectors to store airfoil data */
        ArfXU = (double *)malloc(sizeof(double)*ArfXUNumber[0]);
        ArfYU = (double *)malloc(sizeof(double)*ArfXUNumber[0]);
        ArfXL = (double *)malloc(sizeof(double)*ArfXLNumber[0]);
        ArfYL = (double *)malloc(sizeof(double)*ArfXLNumber[0]);
        /* Read airfoil data */
        readarf(WngTSTpArf,ArfXU,ArfYU,ArfXL,ArfYL,ArfXUNumber,ArfXLNumber);
        /* Interpolate airfoil at the desired coordinates xpline */
        interparf(ycamber,ArfXU,ArfYU,ArfXL,ArfYL,ArfXUNumber,ArfXLNumber,xpline,mp1);
        for (j=0;j<mp1;j++){
            *(ycamberall+j+(2*iTS+1)*mp1)=*(ycamber+j);
        }
    }   
}

void vortexpanel(double *xv, double *yv, double *zv, double *xp, double *yp, double *zp, double dxw, int m, int n)
{
    /* Calculate the collocation points and vector rings on a panel surface xp,yp,zp */
    int i,j,mp1,np1;
    double angp1;
    
    mp1=m+1;
    np1=n+1;
    /* Calculate vortex positions */
    for (j=0;j<np1;j++){
        for (i=0;i<m;i++){
            *(xv+i+j*mp1)=3.0* *(xp+i+j*mp1)/4.0+ *(xp+(i+1)+j*mp1)/4.0;
            *(yv+i+j*mp1)= *(yp+i+j*mp1);
            *(zv+i+j*mp1)=3.0* *(zp+i+j*mp1)/4.0+ *(zp+(i+1)+j*mp1)/4.0;
        }
        /* Trailing edge vortex panels */
        angp1=atan2(*(zp+m+j*mp1)- *(zp+(m-1)+j*mp1),*(xp+m+j*mp1)- *(xp+(m-1)+j*mp1));
        *(xv+m+j*mp1)= *(xp+m+j*mp1)+dxw* cos(angp1);
        *(yv+m+j*mp1)= *(yp+m+j*mp1);
        *(zv+m+j*mp1)= *(zp+m+j*mp1)+dxw* sin(angp1);
    }
}

int findindex(double *ypline, int np1, double value)
{
    /* Find the index of the element of a vector that is closest to a number */ 
    int i,index;
    double dist,mindist;
    
    mindist=10000.0;
    for (i=0;i<np1;i++){
        dist=fabs(value-*(ypline+i));
        if (dist < mindist){
            mindist=dist;
            index=i;
        }
    }
    return index;
}

void vtailsetup(struct liftsurf *pvtail, struct liftsurf *prudder, char *VTailfile, int m, int n)
{
    /* Set up vertical tail */
    
    double pi;
    int i,j,k,nylim,mxlim,nxpos,nypos,nxyzTS;
    int mp1,np1;
    int vtailhere,rudderhere;
    int vtail_nvert,rudder_nvert,vtail_nface,rudder_nface;
    
    double *ypos,*xpos,RDRypos[2],*ylim,*xlim;
    double *ypline,*xpline,*xpgrid,*ypgrid,*zpgrid,yhere,*xTS,*yTS,*yTS2,*zTS,twistcentre,*ycamber,*ycamberall,*ycamberwgl;
    double *xv,*yv,*zv,dxw,minchord;
    double *xvtail,*yvtail,*zvtail,*xrudder,*yrudder,*zrudder,*chordvec,*levec;
    double *xvvtail,*yvvtail,*zvvtail,*xvrudder,*yvrudder,*zvrudder;
    double zplineRt,zplineTp,twistangle,dihedral,xpTS,ypTS,zpTS;
    int *ijvtail, *ijrudder;

    FILE *fp1;
    int iTS, VTTSNumber, cond, ChkRdr, npTS, npTSp1, ndouble;
    int RDRinds[2][2],dummyint,OptVTFusMounted,OptVTTailMounted,OptVTWingMounted;
    char line[110], code[7], VTType[12], VTArf[12], VTSurfFinish[12], nindex[12];
    double VTTSLength[2], VTTSRtChord[2], VTTSTpChord[2];
    double VTTSSwpLE[2], VTTSDhdr[2], VTTSTR[2];
    double VTSpan, VTLPosFus, VTVPosFus, VTRtChord, VTTpChord, VTSPosFus, VTArea, VTTankVol, VTSwpLE, VTVAngle;
    double VTAR, VTTR, VTVolCoeff,VTRlTpChord,VTRlPosFus,VTRlVPosHT;
    double VTArea_Vs_WngArea,VTRlVPosWng;
    double RdrSpan, RdrPosSpan, RdrGearRatio, RdrArea, RdrHingeLoc,RdrRlChord,RdrRlSpan;
    double RdrMxDDflct, RdrMxUDflct, RdrRtChord, RdrTpChord, RdrSMC;

    /* Read VTail.arp and extract vertical tail description */
    printf("Reading from %s\n",VTailfile);
    fp1 = fopen(VTailfile,"r");
    if(fp1 == NULL) {
        fprintf(stderr,"Error:  Could not open %s\n",VTailfile);
        exit(1);
    }
    
    pi=atan(1.0)*4.0;
    
    mp1=m+1;
    np1=n+1;
    ypline= (double *)malloc(sizeof(double)*np1);
    xpline= (double *)malloc(sizeof(double)*mp1);
    chordvec= (double *)malloc(sizeof(double)*np1);
    levec= (double *)malloc(sizeof(double)*np1);
    ycamber = (double *)malloc(sizeof(double)*mp1);
    xpgrid= (double *)malloc(sizeof(double)*mp1*np1);
    ypgrid= (double *)malloc(sizeof(double)*mp1*np1);
    zpgrid= (double *)malloc(sizeof(double)*mp1*np1);
    xv= (double *)malloc(sizeof(double)*mp1*np1);
    yv= (double *)malloc(sizeof(double)*mp1*np1);
    zv= (double *)malloc(sizeof(double)*mp1*np1);    

    cond=0;
    while (cond == 0){
        fgets(line, 110, fp1);
        if ( strncmp("VTL101",line,6) == 0 ){
            sscanf(line,"%s\t%s",code,VTType);
        }
        if ( strncmp("VTL102",line,6) == 0 ){
            sscanf(line,"%s\t%s\t%s",code,VTArf,VTSurfFinish);
        }
        if ( strncmp("VTL301",line,6) == 0 ){
            sscanf(line,"%s\t%i %i %i",code,&OptVTFusMounted,&OptVTTailMounted,&OptVTWingMounted);
        }
        if ( strncmp("VTL401",line,6) == 0 ){
            sscanf(line,"%s\t%i %i",code,&dummyint,&VTTSNumber);
            cond=1;
        }
    }    
    char *token = NULL;
    cond=0;
    while(cond == 0){
        fgets(line, 110, fp1);
        if (  VTTSNumber >= 1){
            i=0;
            if ( strncmp("VT1502",line,6) == 0 ){
                sscanf(line,"%s\t%lf\t%lf\t%lf",code,&VTTSLength[i],&VTTSRtChord[i],&VTTSTpChord[i]);
            }
            if ( strncmp("VT1503",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&VTTSSwpLE[i],&VTTSDhdr[i]);
            }
            if ( strncmp("VT1504",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&VTTSTR[i]);
            }
        }
        if (  VTTSNumber >= 2){
            i=1;
            if ( strncmp("VT2502",line,6) == 0 ){
                sscanf(line,"%s\t%lf\t%lf\t%lf",code,&VTTSLength[i],&VTTSRtChord[i],&VTTSTpChord[i]);
            }
            if ( strncmp("VT2503",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&VTTSSwpLE[i],&VTTSDhdr[i]);
            }
            if ( strncmp("VT2504",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&VTTSTR[i]);
            }
        }
        if ( strncmp("VTL601",line,6) == 0 ){
            sscanf(line,"%s %lf %lf %lf %lf %lf %lf",code,&VTSpan,&VTLPosFus,&VTVPosFus,&VTRtChord,&VTTpChord,&VTSPosFus);
        }
        if ( strncmp("VTL602",line,6) == 0 ){
            sscanf(line,"%s %lf",code,&VTArea);
        }
        if ( strncmp("VTL604",line,6) == 0 ){
            sscanf(line,"%s %lf %lf",code,&VTSwpLE,&VTVAngle);
        }
        if ( strncmp("VTL605",line,6) == 0 ){
            sscanf(line,"%s %lf %lf %lf %lf %lf %lf",code,&VTAR,&VTTR,&VTVolCoeff,&VTRlTpChord,&VTRlPosFus,&VTRlVPosHT);
        }
        if ( strncmp("RDR201",line,6) == 0 ){
            sscanf(line,"%s %i",code,&ChkRdr);
            cond=1;
        }
    }
    cond=0;
    while(cond == 0){
        fgets(line, 110, fp1);
        if ( ChkRdr == 1 ){
            if ( strncmp("RDR601",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&RdrSpan,&RdrPosSpan);
            }
            if ( strncmp("RDR602",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&RdrArea);
            }
            if ( strncmp("RDR603",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&RdrHingeLoc,&RdrRlChord,&RdrRlSpan);
            }
            if ( strncmp("RDR604",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&RdrMxDDflct,&RdrMxUDflct);
            }
            if ( strncmp("RDR607",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&RdrRtChord,&RdrTpChord,&RdrSMC);
            }
        }
        if ( strncmp("RDT201",line,6) == 0 ){
            cond=1;
        }
    }
    
    /* Convert to rad */
    minchord=100; /* Initialize the calculation of the minimum chord of the wing */
    for (i=0;i<VTTSNumber;i++){
        VTTSSwpLE[i]=VTTSSwpLE[i]*pi/180.0;
        VTTSDhdr[i]=VTTSDhdr[i]*pi/180.0;
        if (VTTSTpChord[i] < minchord){
            minchord=VTTSTpChord[i];
        }
    } 
    
    nypos=VTTSNumber+3;
    nxpos=3;
    nxyzTS=VTTSNumber+1;
    /* Create vector containing y-coords of TS and Rdr */
    ypos= (double *)malloc(sizeof(double)*nypos); 
    xTS= (double *)malloc(sizeof(double)*nxyzTS); 
    yTS= (double *)malloc(sizeof(double)*nxyzTS); 
    yTS2= (double *)malloc(sizeof(double)*nxyzTS); 
    zTS= (double *)malloc(sizeof(double)*nxyzTS); 
    *ypos=0.0;*yTS=0.0;*yTS2=0.0;*xTS=0.0;*zTS=0.0;
    for (i=1;i<VTTSNumber+1;i++){
        *(ypos+i)=VTTSLength[i-1]+*(ypos+i-1);
        *(yTS+i)=*(ypos+i);
        *(yTS2+i)=cos(pi/2.0-VTTSDhdr[i-1])*VTTSLength[i-1]+*(yTS2+i-1);
        *(xTS+i)=VTTSLength[i-1]*tan(VTTSSwpLE[i-1])+*(xTS+i-1);
        *(zTS+i)=sin(pi/2.0-VTTSDhdr[i-1])*VTTSLength[i-1]+*(zTS+i-1);
    }
    *(ypos+i)=RdrPosSpan;i++;
    *(ypos+i)=RdrSpan+RdrPosSpan;
    
    /* Sort ypos vector */
    qsort(ypos, nypos, sizeof(double), compare_function);
    /* Remove repeated values from ypos vector */   
    ndouble=checkdoubles(ypos,nypos);
    while (ndouble > 0){
        nypos=finddoubles(ypos,nypos);
        ndouble=checkdoubles(ypos,nypos);
    }

    /* Check if number of spanwise panels is sufficient */
    if (n <= nypos){
        fprintf(stderr,"Error:  The number of spanwise panels on the vertical tail should be larger than %i\n",nypos);
        exit(1);
    }    
    /* Create vector containing the full spanwise grid */
    createypline(ypos,nypos,ypline,np1);
    /* Check between which elements of ypline lies the elevator */    
    RDRinds[0][1]=findindex(ypline,np1,RdrPosSpan);
    RDRinds[1][1]=findindex(ypline,np1,RdrSpan+RdrPosSpan); 
    
    /* Create vector containing x-coords of TS and Rdr */
    xpos= (double *)malloc(sizeof(double)*nxpos); 
    *(xpos+0)=0.0;
    *(xpos+1)=(100.0-RdrRlChord)/100.0;
    *(xpos+2)=1.0;
    /* Sort xpos vector */
    qsort(xpos, nxpos, sizeof(double), compare_function);
    /* Remove repeated values from xpos vector */   
    ndouble=checkdoubles(xpos,nxpos);
    while (ndouble > 0){
        nxpos=finddoubles(xpos,nxpos);
        ndouble=checkdoubles(xpos,nxpos);
    }
    /* Make sure the trailing edge lies at 1 (it might be 0.994 or something, which can be inaccurate for large chord values) */
    *(xpos+nxpos-1)=1.0;
    /* Check if number of chordwise panels is sufficient */
    if (m <= nxpos){
        fprintf(stderr,"Error:  The number of chordwise panels on the vertical tail should be larger than %i\n",nxpos);
        exit(1);
    }
    /* Create vector containing the full spanwise grid */
    createypline(xpos,nxpos,xpline,mp1);
    /* Check between which elements of xpline lies the elevator */
    RDRinds[0][0]=findindex(xpline,mp1,(100.0-RdrRlChord)/100.0);
    RDRinds[1][0]=m; /* Will always lie on the trailing edge */ 
    
    /* Create matrix of non-dimensional camber lines */
    rewind(fp1);
    ycamberall = (double *)malloc(sizeof(double)*mp1*VTTSNumber*2);

    /* Re-read VTail.arp file to find root and tip airfoil names */
    cond=0;
    while(cond == 0){
        fgets(line, 110, fp1);
        if (  VTTSNumber >= 1){
            i=0;
            if ( strncmp("VT1501",line,6) == 0 ){
                treatarf(line,xpline,ycamber,ycamberall,mp1,i); /* Airfoils */
           }
        }
        if (  VTTSNumber >= 2){
            i=1;
            if ( strncmp("VT2501",line,6) == 0 ){
                treatarf(line,xpline,ycamber,ycamberall,mp1,i); /* Airfoils */
            } 
        }       
        if ( strncmp("RDT201",line,6) == 0 ){
            cond=1;
        }
    }              
    fclose(fp1);

    /* Create complete vertical tail grid, split it into fin and rudder later */
    iTS=0;
    for (j=0;j<np1;j++){
        /* Find out which trapezoidal section we're on */
        yhere=*(ypline+j);
        for (i=1;i<nxyzTS-1;i++){
            if (yhere > *(yTS+i))
                iTS=i;
        }
        for (i=0;i<mp1;i++){
            *(ypgrid+i+j*mp1)=yhere;
            *(chordvec+j)=(VTTSTpChord[iTS]-VTTSRtChord[iTS])/VTTSLength[iTS]*(yhere-*(yTS+iTS))+VTTSRtChord[iTS]; /* Local chord length */
            *(levec+j)=*(xTS+iTS)+(yhere-*(yTS+iTS))*tan(VTTSSwpLE[iTS]); /* Local leading edge position */
            *(xpgrid+i+j*mp1)=*(chordvec+j)* *(xpline+i)+ *(levec+j);
            zplineRt=*(ycamberall+i+2*iTS*mp1)*VTTSRtChord[iTS]; /* Root camber line of trapezoidal section */
            zplineTp=*(ycamberall+i+(2*iTS+1)*mp1)*VTTSTpChord[iTS]; /* Tip camber line of trapezoidal section */
            *(zpgrid+i+j*mp1)=(zplineTp-zplineRt)/VTTSLength[iTS]*(yhere-*(yTS+iTS))+zplineRt;
        }
    } 

    /* wake shedding distance */
    dxw=0.3*minchord/m;
    vortexpanel(xv,yv,zv,xpgrid,ypgrid,zpgrid,dxw,m,n);  
         
    /* Assign wing grid cells to stabilizer and elevator */
    xvtail= (double *)malloc(sizeof(double)*mp1*np1); 
    yvtail= (double *)malloc(sizeof(double)*mp1*np1); 
    zvtail= (double *)malloc(sizeof(double)*mp1*np1); 
    xvvtail= (double *)malloc(sizeof(double)*mp1*np1); 
    yvvtail= (double *)malloc(sizeof(double)*mp1*np1); 
    zvvtail= (double *)malloc(sizeof(double)*mp1*np1); 
    ijvtail= (int *)malloc(sizeof(int)*mp1*np1*2);
    xrudder= (double *)malloc(sizeof(double)*mp1*np1); 
    yrudder= (double *)malloc(sizeof(double)*mp1*np1); 
    zrudder= (double *)malloc(sizeof(double)*mp1*np1); 
    xvrudder= (double *)malloc(sizeof(double)*mp1*np1); 
    yvrudder= (double *)malloc(sizeof(double)*mp1*np1); 
    zvrudder= (double *)malloc(sizeof(double)*mp1*np1); 
    ijrudder= (int *)malloc(sizeof(int)*mp1*np1*2);
    vtail_nvert=0;rudder_nvert=0; /* Initialize the number of vertices in all the lifting surfaces */  
       
    iTS=0;
    for (j=0;j<np1;j++){
        /* Find out which trapezoidal section we're on */
        yhere=*(ypline+j);
        for (i=1;i<nxyzTS-1;i++){
            if (yhere > *(yTS+i))
                iTS=i;
        }
        dihedral=pi/2.0-VTTSDhdr[iTS];
        twistangle=0.0;
        for (i=0;i<mp1;i++){
            vtailhere=0;
            rudderhere=0;
            /* Check if this point lies on the rudder */
            if (i >= RDRinds[0][0] && i <= RDRinds[1][0] && j >= RDRinds[0][1] && j <= RDRinds[1][1]){ /* Elevator */
                rudderhere=1;
                if (i == RDRinds[0][0]){ /* Rudder leading edge */
                    vtailhere=1;
                }
                if (j == RDRinds[0][1]){ /* Rudder inboard edge */
                    vtailhere=1;
                }
                if (j == RDRinds[1][1]){ /* Rudder outboard edge */
                    vtailhere=1;
                }
            }
            if (vtailhere == 1){
                if (rudderhere == 1 && RDRinds[0][1] == 0 && i > RDRinds[0][0]  && j == 0) /* if rudder starts at tail root */
                    vtailhere=0;
                if (rudderhere == 1 && RDRinds[1][1] == np1-1 && i > RDRinds[0][0] && j > RDRinds[0][1]) /* if rudder extends to wingtip */
                    vtailhere=0;
            }else{
                if (rudderhere == 0) /* If this point does not lie on the rudder */
                    vtailhere=1;
            }
            /* Impose dihedral on geometric panels*/
            xpTS=*(xpgrid+i+j*mp1);
            ypTS=*(ypgrid+i+j*mp1);
            zpTS=*(zpgrid+i+j*mp1);
            *(ypgrid+i+j*mp1)=cos(dihedral)*(ypTS-*(yTS+iTS))+*(yTS2+iTS);
            *(zpgrid+i+j*mp1)=sin(dihedral)*(ypTS-*(yTS+iTS))+zpTS+*(zTS+iTS);
            /* Impose dihedral on vortex panels*/
            xpTS=*(xv+i+j*mp1);
            ypTS=*(yv+i+j*mp1);
            zpTS=*(zv+i+j*mp1);
            *(yv+i+j*mp1)=cos(dihedral)*(ypTS-*(yTS+iTS))+*(yTS2+iTS);
            *(zv+i+j*mp1)=sin(dihedral)*(ypTS-*(yTS+iTS))+zpTS+*(zTS+iTS);
            
            /* Store xpgrid, ypgrid, zpgrid, xv, yv, zv, i and j in the corresponding arrays */
            if (rudderhere == 1){
                rudder_nvert++;
                *(xrudder+rudder_nvert-1)=*(xpgrid+i+j*mp1)+VTLPosFus;
                *(yrudder+rudder_nvert-1)=*(ypgrid+i+j*mp1);
                *(zrudder+rudder_nvert-1)=*(zpgrid+i+j*mp1)+VTVPosFus;
                *(xvrudder+rudder_nvert-1)=*(xv+i+j*mp1)+VTLPosFus;
                *(yvrudder+rudder_nvert-1)=*(yv+i+j*mp1);
                *(zvrudder+rudder_nvert-1)=*(zv+i+j*mp1)+VTVPosFus;
                *(ijrudder+rudder_nvert-1)=i;
                *(ijrudder+mp1*np1+rudder_nvert-1)=j;
            }
            if (vtailhere == 1){
                vtail_nvert++;
                *(xvtail+vtail_nvert-1)=*(xpgrid+i+j*mp1)+VTLPosFus;
                *(yvtail+vtail_nvert-1)=*(ypgrid+i+j*mp1);
                *(zvtail+vtail_nvert-1)=*(zpgrid+i+j*mp1)+VTVPosFus;
                *(xvvtail+vtail_nvert-1)=*(xv+i+j*mp1)+VTLPosFus;
                *(yvvtail+vtail_nvert-1)=*(yv+i+j*mp1);
                *(zvvtail+vtail_nvert-1)=*(zv+i+j*mp1)+VTVPosFus;
                *(ijvtail+vtail_nvert-1)=i;
                *(ijvtail+mp1*np1+vtail_nvert-1)=j;
            }
        }
    }     
    
    /* Count number of panels on all halves of lifting surfaces */
    rudder_nface=countfaces(ijrudder,rudder_nvert,mp1,np1);
    vtail_nface=countfaces(ijvtail,vtail_nvert,mp1,np1);   
    
    /* Create all the panel information in the rudder liftsurf structure */
    prudder->nface=2*rudder_nface;
    prudder->nvert=2*rudder_nvert;
    prudder->faces=(int *)malloc(sizeof(int)*prudder->nface*4); 
    prudder->shedding=(int *)malloc(sizeof(int)*prudder->nface); 
    arrangefaces(ijrudder,mp1,np1,prudder);  
    
    /* Create all the panel information in the vtail liftsurf structure */
    pvtail->nface=2*vtail_nface;
    pvtail->nvert=2*vtail_nvert;
    pvtail->faces=(int *)malloc(sizeof(int)*pvtail->nface*4); 
    pvtail->shedding=(int *)malloc(sizeof(int)*pvtail->nface); 
    arrangefaces(ijvtail,mp1,np1,pvtail);

    /* Copy all panel vertex information to the relevant liftsurf structures */
    assignvertices_fin(prudder,xrudder,yrudder,zrudder,xvrudder,yvrudder,zvrudder,OptVTFusMounted,OptVTTailMounted,OptVTWingMounted);
    assignvertices_fin(pvtail,xvtail,yvtail,zvtail,xvvtail,yvvtail,zvvtail,OptVTFusMounted,OptVTTailMounted,OptVTWingMounted);
}


void htailsetup(struct liftsurf *phtail, struct liftsurf *pelevator, char *HTailfile, int m, int n)
{
    /* Set up horizontal tail */
    
    double pi;
    int i,j,k,nylim,mxlim,nxpos,nypos,nxyzTS;
    int mp1,np1;
    int htailhere,elevatorhere;
    int htail_nvert,elevator_nvert,htail_nface,elevator_nface;
    
    double *ypos,*xpos,ELVypos[2],*ylim,*xlim;
    double *ypline,*xpline,*xpgrid,*ypgrid,*zpgrid,yhere,*xTS,*yTS,*yTS2,*zTS,twistcentre,*ycamber,*ycamberall,*ycamberwgl;
    double *xv,*yv,*zv,dxw,minchord;
    double *xhtail,*yhtail,*zhtail,*xelevator,*yelevator,*zelevator,*chordvec,*levec;
    double *xvhtail,*yvhtail,*zvhtail,*xvelevator,*yvelevator,*zvelevator;
    double zplineRt,zplineTp,twistangle,dihedral,xpTS,ypTS,zpTS;
    int *ijhtail, *ijelevator;
       
    FILE *fp1;
    int iTS, HTTSNumber, cond, ChkElev, npTS, npTSp1, ndouble;
    int ELVinds[2][2],dummyint;
    char line[110], code[7], HTType[12], HTArf[12], HTSurfFinish[12], nindex[12];
    double HTTSLength[3], HTTSRtChord[3], HTTSTpChord[3];
    double HTTSLPosFus[3], HTTSSPosFus[3], HTTSVPosFus[3], HTTSSwpMxRlThick[3];
    double HTTSRtInc[3], HTTSTpInc[3], HTTSSwpLE[3], HTTSDhdr[3], HTTSTwist[3], HTTSTR[3];
    double HTSpan, HTLPosFus, HTVPosFus, HTRtChord, HTTpChord, HTArea, HTTankVol, HTSwpLE, HTTwist, HTRlInc, HTInc, HTDhdrl;
    double HTArea_Vs_WngArea, HTAR, HTTR, HTVolCoeff;
    double ElevSpan, ElevPosSpan, ElevGearRatio, ElevArea, ElevHingeLoc,ElevRlChord,ElevRlSpan;
    double ElevMxDDflct, ElevMxUDflct, ElevRtChord, ElevTpChord, ElevSMC;
    
    /* Read Htail.arp and extract horizontal tail description */
    printf("Reading from %s\n",HTailfile);
    fp1 = fopen(HTailfile,"r");
    if(fp1 == NULL) {
        fprintf(stderr,"Error:  Could not open %s\n",HTailfile);
        exit(1);
    }
    
    pi=atan(1.0)*4.0;
    
    mp1=m+1;
    np1=n+1;
    ypline= (double *)malloc(sizeof(double)*np1);
    xpline= (double *)malloc(sizeof(double)*mp1);
    chordvec= (double *)malloc(sizeof(double)*np1);
    levec= (double *)malloc(sizeof(double)*np1);
    ycamber = (double *)malloc(sizeof(double)*mp1);
    xpgrid= (double *)malloc(sizeof(double)*mp1*np1);
    ypgrid= (double *)malloc(sizeof(double)*mp1*np1);
    zpgrid= (double *)malloc(sizeof(double)*mp1*np1);
    xv= (double *)malloc(sizeof(double)*mp1*np1);
    yv= (double *)malloc(sizeof(double)*mp1*np1);
    zv= (double *)malloc(sizeof(double)*mp1*np1);
    
    cond=0;
    while (cond == 0){
        fgets(line, 110, fp1);
        if ( strncmp("HTL101",line,6) == 0 ){
            sscanf(line,"%s\t%s",code,HTType);
        }
        if ( strncmp("HTL102",line,6) == 0 ){
            sscanf(line,"%s\t%s\t%s",code,HTArf,HTSurfFinish);
        }
        if ( strncmp("HTL409",line,6) == 0 ){
            sscanf(line,"%s\t%i %i",code,&dummyint,&HTTSNumber);
            cond=1;
        }
    }
    
    char *token = NULL;
    cond=0;
    while(cond == 0){
        fgets(line, 110, fp1);
        if (  HTTSNumber >= 1){
            i=0;
            if ( strncmp("HT1502",line,6) == 0 ){
                sscanf(line,"%s\t%lf\t%lf\t%lf",code,&HTTSLength[i],&HTTSRtChord[i],&HTTSTpChord[i]);
            }
            if ( strncmp("HT1503",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf %lf %lf",code,&HTTSRtInc[i],&HTTSTpInc[i],&HTTSSwpLE[i],&HTTSDhdr[i],&HTTSTwist[i]);
            }
            if ( strncmp("HT1504",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&HTTSTR[i]);
            }
        }
        if (  HTTSNumber >= 2){
            i=1;
            if ( strncmp("HT2502",line,6) == 0 ){
                sscanf(line,"%s\t%lf\t%lf\t%lf",code,&HTTSLength[i],&HTTSRtChord[i],&HTTSTpChord[i]);
            }
            if ( strncmp("HT2503",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf %lf %lf",code,&HTTSRtInc[i],&HTTSTpInc[i],&HTTSSwpLE[i],&HTTSDhdr[i],&HTTSTwist[i]);
            }
            if ( strncmp("HT2504",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&HTTSTR[i]);
            }
        }
        if (  HTTSNumber >= 3){
            i=2;
            if ( strncmp("HT3502",line,6) == 0 ){
                sscanf(line,"%s\t%lf\t%lf\t%lf",code,&HTTSLength[i],&HTTSRtChord[i],&HTTSTpChord[i]);
            }
            if ( strncmp("HT3503",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf %lf %lf",code,&HTTSRtInc[i],&HTTSTpInc[i],&HTTSSwpLE[i],&HTTSDhdr[i],&HTTSTwist[i]);
           }
            if ( strncmp("HT3504",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&HTTSTR[i]);
            }
        }
        if ( strncmp("HTL601",line,6) == 0 ){
            sscanf(line,"%s %lf %lf %lf %lf %lf",code,&HTSpan,&HTLPosFus,&HTVPosFus,&HTRtChord,&HTTpChord);
        }
        if ( strncmp("HTL602",line,6) == 0 ){
            sscanf(line,"%s %lf",code,&HTArea);
        }
        if ( strncmp("HTL604",line,6) == 0 ){
            sscanf(line,"%s %lf %lf %lf %lf %lf",code,&HTSwpLE,&HTTwist,&HTRlInc,&HTInc,&HTDhdrl);
        }
        if ( strncmp("HTL605",line,6) == 0 ){
            sscanf(line,"%s %lf %lf %lf",code,&HTAR,&HTTR,&HTVolCoeff);
        }
        if ( strncmp("ELV201",line,6) == 0 ){
            sscanf(line,"%s %i",code,&ChkElev);
            cond=1;
        }
    }
    cond=0;
    while(cond == 0){
        fgets(line, 110, fp1);
        if ( ChkElev == 1 ){
            if ( strncmp("ELV601",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&ElevSpan,&ElevPosSpan);
            }
            if ( strncmp("ELV602",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&ElevArea);
            }
            if ( strncmp("ELV603",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&ElevHingeLoc,&ElevRlChord,&ElevRlSpan);
            }
            if ( strncmp("ELV604",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&ElevMxDDflct,&ElevMxUDflct);
            }
            if ( strncmp("ELV607",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&ElevRtChord,&ElevTpChord,&ElevSMC);
            }
        }
        if ( strncmp("ELT201",line,6) == 0 ){
            cond=1;
        }
    }
    
    /* Convert to rad */
    minchord=100; /* Initialize the calculation of the minimum chord of the wing */
    for (i=0;i<HTTSNumber;i++){
        HTTSSwpLE[i]=HTTSSwpLE[i]*pi/180.0;
        HTTSDhdr[i]=HTTSDhdr[i]*pi/180.0;
        HTTSTwist[i]=HTTSTwist[i]*pi/180.0;
        HTTSRtInc[i]=HTTSRtInc[i]*pi/180.0;
        HTTSTpInc[i]=HTTSTpInc[i]*pi/180.0;
        if (HTTSTpChord[i] < minchord){
            minchord=HTTSTpChord[i];
        }
    }
    
    nypos=HTTSNumber+3;
    nxpos=3;
    nxyzTS=HTTSNumber+1;
    /* Create vector containing y-coords of TS and Elev */
    ypos= (double *)malloc(sizeof(double)*nypos); 
    xTS= (double *)malloc(sizeof(double)*nxyzTS); 
    yTS= (double *)malloc(sizeof(double)*nxyzTS); 
    yTS2= (double *)malloc(sizeof(double)*nxyzTS); 
    zTS= (double *)malloc(sizeof(double)*nxyzTS); 
    *ypos=0.0;*yTS=0.0;*yTS2=0.0;*xTS=0.0;*zTS=0.0;
    for (i=1;i<HTTSNumber+1;i++){
        *(ypos+i)=HTTSLength[i-1]+*(ypos+i-1);
        *(yTS+i)=*(ypos+i);
        *(yTS2+i)=cos(HTTSDhdr[i-1])*HTTSLength[i-1]+*(yTS2+i-1);
        *(xTS+i)=HTTSLength[i-1]*tan(HTTSSwpLE[i-1])+*(xTS+i-1);
        *(zTS+i)=sin(HTTSDhdr[i-1])*HTTSLength[i-1]+*(zTS+i-1);
    }
    *(ypos+i)=ElevPosSpan;i++;
    *(ypos+i)=ElevSpan/2.0+ElevPosSpan; /* ElevSpan is the total elevator span */
    
    /* Sort ypos vector */
    qsort(ypos, nypos, sizeof(double), compare_function);
    /* Remove repeated values from ypos vector */   
    ndouble=checkdoubles(ypos,nypos);
    while (ndouble > 0){
        nypos=finddoubles(ypos,nypos);
        ndouble=checkdoubles(ypos,nypos);
    }
    
    /* Check if number of spanwise panels is sufficient */
    if (n <= nypos){
        fprintf(stderr,"Error:  The number of spanwise panels on the horizontal tail should be larger than %i\n",nypos);
        exit(1);
    }    
    /* Create vector containing the full spanwise grid */
    createypline(ypos,nypos,ypline,np1);
    /* Check between which elements of ypline lies the elevator */    
    ELVinds[0][1]=findindex(ypline,np1,ElevPosSpan);
    ELVinds[1][1]=findindex(ypline,np1,ElevSpan/2+ElevPosSpan);

    /* Create vector containing x-coords of TS and Elev */
    xpos= (double *)malloc(sizeof(double)*nxpos); 
    *(xpos+0)=0.0;
    *(xpos+1)=(100.0-ElevRlChord)/100.0;
    *(xpos+2)=1.0;
    /* Sort xpos vector */
    qsort(xpos, nxpos, sizeof(double), compare_function);
    /* Remove repeated values from ypos vector */   
    ndouble=checkdoubles(xpos,nxpos);
    while (ndouble > 0){
        nxpos=finddoubles(xpos,nxpos);
        ndouble=checkdoubles(xpos,nxpos);
    }
    /* Make sure the trailing edge lies at 1 (it might be 0.994 or something, which can be inaccurate for large chord values) */
    *(xpos+nxpos-1)=1.0;
    /* Check if number of chordwise panels is sufficient */
    if (m <= nxpos){
        fprintf(stderr,"Error:  The number of chordwise panels on the horizontal tail should be larger than %i\n",nxpos);
        exit(1);
    }
    /* Create vector containing the full spanwise grid */
    createypline(xpos,nxpos,xpline,mp1);
    /* Check between which elements of xpline lies the elevator */
    ELVinds[0][0]=findindex(xpline,mp1,(100.0-ElevRlChord)/100.0);
    ELVinds[1][0]=m; /* Will always lie on the trailing edge */    

    /* Create matrix of non-dimensional camber lines */
    rewind(fp1);
    ycamberall = (double *)malloc(sizeof(double)*mp1*HTTSNumber*2);

    /* Re-read HTail.arp file to find root and tip airfoil names */
    cond=0;
    while(cond == 0){
        fgets(line, 110, fp1);
        if (  HTTSNumber >= 1){
            i=0;
            if ( strncmp("HT1501",line,6) == 0 ){
                treatarf(line,xpline,ycamber,ycamberall,mp1,i); /* Airfoils */
           }
        }
        if (  HTTSNumber >= 2){
            i=1;
            if ( strncmp("HT2501",line,6) == 0 ){
                treatarf(line,xpline,ycamber,ycamberall,mp1,i); /* Airfoils */
            } 
        }
        if (  HTTSNumber >= 3){
            i=2;
            if ( strncmp("HT3501",line,6) == 0 ){
                treatarf(line,xpline,ycamber,ycamberall,mp1,i); /* Airfoils */
            } 
        }       
        if ( strncmp("ELT201",line,6) == 0 ){
            cond=1;
        }
    }              
    fclose(fp1);

    /* Create complete horizontal tail grid, split it into stabilizer and elevator later */
    iTS=0;
    for (j=0;j<np1;j++){
        /* Find out which trapezoidal section we're on */
        yhere=*(ypline+j);
        for (i=1;i<nxyzTS-1;i++){
            if (yhere > *(yTS+i))
                iTS=i;
        }
        for (i=0;i<mp1;i++){
            *(ypgrid+i+j*mp1)=yhere;
            *(chordvec+j)=(HTTSTpChord[iTS]-HTTSRtChord[iTS])/HTTSLength[iTS]*(yhere-*(yTS+iTS))+HTTSRtChord[iTS]; /* Local chord length */
            *(levec+j)=*(xTS+iTS)+(yhere-*(yTS+iTS))*tan(HTTSSwpLE[iTS]); /* Local leading edge position */
            *(xpgrid+i+j*mp1)=*(chordvec+j)* *(xpline+i)+ *(levec+j);
            zplineRt=*(ycamberall+i+2*iTS*mp1)*HTTSRtChord[iTS]; /* Root camber line of trapezoidal section */
            zplineTp=*(ycamberall+i+(2*iTS+1)*mp1)*HTTSTpChord[iTS]; /* Tip camber line of trapezoidal section */
            *(zpgrid+i+j*mp1)=(zplineTp-zplineRt)/HTTSLength[iTS]*(yhere-*(yTS+iTS))+zplineRt;
        }
    }    

    /* wake shedding distance */
    dxw=0.3*minchord/m;
    vortexpanel(xv,yv,zv,xpgrid,ypgrid,zpgrid,dxw,m,n);  
         
    /* Assign wing grid cells to stabilizer and elevator */
    xhtail= (double *)malloc(sizeof(double)*mp1*np1); 
    yhtail= (double *)malloc(sizeof(double)*mp1*np1); 
    zhtail= (double *)malloc(sizeof(double)*mp1*np1); 
    xvhtail= (double *)malloc(sizeof(double)*mp1*np1); 
    yvhtail= (double *)malloc(sizeof(double)*mp1*np1); 
    zvhtail= (double *)malloc(sizeof(double)*mp1*np1); 
    ijhtail= (int *)malloc(sizeof(int)*mp1*np1*2);
    xelevator= (double *)malloc(sizeof(double)*mp1*np1); 
    yelevator= (double *)malloc(sizeof(double)*mp1*np1); 
    zelevator= (double *)malloc(sizeof(double)*mp1*np1); 
    xvelevator= (double *)malloc(sizeof(double)*mp1*np1); 
    yvelevator= (double *)malloc(sizeof(double)*mp1*np1); 
    zvelevator= (double *)malloc(sizeof(double)*mp1*np1); 
    ijelevator= (int *)malloc(sizeof(int)*mp1*np1*2);
    htail_nvert=0;elevator_nvert=0; /* Initialize the number of vertices in all the lifting surfaces */

    iTS=0;
    for (j=0;j<np1;j++){
        /* Find out which trapezoidal section we're on */
        yhere=*(ypline+j);
        for (i=1;i<nxyzTS-1;i++){
            if (yhere > *(yTS+i))
                iTS=i;
        }
        if (iTS < HTTSNumber){
            twistangle=(HTTSTpInc[iTS]-HTTSRtInc[iTS])/HTTSLength[iTS]* (*(ypgrid+j*mp1)-*(yTS+iTS))+HTTSRtInc[iTS];
            dihedral=HTTSDhdr[iTS];
        }else{
            /* Winglet is untwisted and has the same twist angle as wingtip */
            dihedral=HTTSDhdr[iTS];
        }
        for (i=0;i<mp1;i++){
            htailhere=0;
            elevatorhere=0;
            /* Check if this point lies on the elevator */
            if (i >= ELVinds[0][0] && i <= ELVinds[1][0] && j >= ELVinds[0][1] && j <= ELVinds[1][1]){ /* Elevator */
                elevatorhere=1;
                if (i == ELVinds[0][0]){ /* Elevator leading edge */
                    htailhere=1;
                }
                if (j == ELVinds[0][1]){ /* Elevator inboard edge */
                    htailhere=1;
                }
                if (j == ELVinds[1][1]){ /* Elevator outboard edge */
                    htailhere=1;
                }
            }
            if (htailhere == 1){
                if (elevatorhere == 1 && ELVinds[0][1] == 0 && i > ELVinds[0][0]  && j == 0) /* if elevator starts at tail root */
                    htailhere=0;
                if (elevatorhere == 1 && ELVinds[1][1] == np1-1 && i > ELVinds[0][0] && j > ELVinds[0][1]) /* if elevator extends to wingtip */
                    htailhere=0;
            }else{
                if (elevatorhere == 0) /* If this point does not lie on the elevator */
                    htailhere=1;
            }
            /* Impose twist on geometric panels*/
            xpTS=*(xpgrid+i+j*mp1);
            ypTS=*(ypgrid+i+j*mp1);
            zpTS=*(zpgrid+i+j*mp1);
//            twistcentre=*(levec+j)+*(chordvec+j)*0.5; /* The wing twists around the half-chord */
            twistcentre=0.0;
            *(xpgrid+i+j*mp1)=cos(twistangle)*(xpTS-twistcentre)+sin(twistangle)*zpTS+twistcentre;
            *(zpgrid+i+j*mp1)=-sin(twistangle)*(xpTS-twistcentre)+cos(twistangle)*zpTS;
            zpTS=*(zpgrid+i+j*mp1);
            /* Impose dihedral on geometric panels*/
            *(ypgrid+i+j*mp1)=cos(dihedral)*(ypTS-*(yTS+iTS))+*(yTS2+iTS);
            *(zpgrid+i+j*mp1)=sin(dihedral)*(ypTS-*(yTS+iTS))+zpTS+*(zTS+iTS);
            /* Impose twist on vortex panels*/
            xpTS=*(xv+i+j*mp1);
            ypTS=*(yv+i+j*mp1);
            zpTS=*(zv+i+j*mp1);
            *(xv+i+j*mp1)=cos(twistangle)*(xpTS-twistcentre)+sin(twistangle)*zpTS+twistcentre;
            *(zv+i+j*mp1)=-sin(twistangle)*(xpTS-twistcentre)+cos(twistangle)*zpTS;
            zpTS=*(zv+i+j*mp1);
            /* Impose dihedral on vortex panels*/
            *(yv+i+j*mp1)=cos(dihedral)*(ypTS-*(yTS+iTS))+*(yTS2+iTS);
            *(zv+i+j*mp1)=sin(dihedral)*(ypTS-*(yTS+iTS))+zpTS+*(zTS+iTS);
            
            /* Store xpgrid, ypgrid, zpgrid, xv, yv, zv, i and j in the corresponding arrays */
            if (elevatorhere == 1){
                elevator_nvert++;
                *(xelevator+elevator_nvert-1)=*(xpgrid+i+j*mp1)+HTLPosFus;
                *(yelevator+elevator_nvert-1)=*(ypgrid+i+j*mp1);
                *(zelevator+elevator_nvert-1)=*(zpgrid+i+j*mp1)+HTVPosFus;
                *(xvelevator+elevator_nvert-1)=*(xv+i+j*mp1)+HTLPosFus;
                *(yvelevator+elevator_nvert-1)=*(yv+i+j*mp1);
                *(zvelevator+elevator_nvert-1)=*(zv+i+j*mp1)+HTVPosFus;
                *(ijelevator+elevator_nvert-1)=i;
                *(ijelevator+mp1*np1+elevator_nvert-1)=j;
            }
            if (htailhere == 1){
                htail_nvert++;
                *(xhtail+htail_nvert-1)=*(xpgrid+i+j*mp1)+HTLPosFus;
                *(yhtail+htail_nvert-1)=*(ypgrid+i+j*mp1);
                *(zhtail+htail_nvert-1)=*(zpgrid+i+j*mp1)+HTVPosFus;
                *(xvhtail+htail_nvert-1)=*(xv+i+j*mp1)+HTLPosFus;
                *(yvhtail+htail_nvert-1)=*(yv+i+j*mp1);
                *(zvhtail+htail_nvert-1)=*(zv+i+j*mp1)+HTVPosFus;
                *(ijhtail+htail_nvert-1)=i;
                *(ijhtail+mp1*np1+htail_nvert-1)=j;
            }
        }
    }     

    /* Count number of panels on all halves of lifting surfaces */
    elevator_nface=countfaces(ijelevator,elevator_nvert,mp1,np1);
    htail_nface=countfaces(ijhtail,htail_nvert,mp1,np1);   
    
    /* Create all the panel information in the elevator liftsurf structure */
    pelevator->nface=2*elevator_nface;
    pelevator->nvert=2*elevator_nvert;
    pelevator->faces=(int *)malloc(sizeof(int)*pelevator->nface*4); 
    pelevator->shedding=(int *)malloc(sizeof(int)*pelevator->nface); 
    arrangefaces(ijelevator,mp1,np1,pelevator);  
    
    /* Create all the panel information in the htail liftsurf structure */
    phtail->nface=2*htail_nface;
    phtail->nvert=2*htail_nvert;
    phtail->faces=(int *)malloc(sizeof(int)*phtail->nface*4); 
    phtail->shedding=(int *)malloc(sizeof(int)*phtail->nface); 
    arrangefaces(ijhtail,mp1,np1,phtail);

    /* Copy all panel vertex information to the relevant liftsurf structures */
    assignvertices(pelevator,xelevator,yelevator,zelevator,xvelevator,yvelevator,zvelevator);
    assignvertices(phtail,xhtail,yhtail,zhtail,xvhtail,yvhtail,zvhtail);     
}

double wingsetup(struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *pwing, char *Wngfile, int m, int n)
{
    double pi;
    int i,j,k,nylim,mxlim,nxpos,nypos,nxyzTS;
    int mp1,np1;
    int flaphere,winghere,aileronhere;
    int flap_nvert,wing_nvert,aileron_nvert,flap_nface,wing_nface,aileron_nface;
    
    double *ypos,*xpos,Ailypos[2],WTEDypos[2],*ylim,*xlim;
    double *ypline,*xpline,*xpgrid,*ypgrid,*zpgrid,yhere,*xTS,*yTS,*yTS2,*zTS,twistcentre,*ycamber,*ycamberall,*ycamberwgl;
    double *xv,*yv,*zv,dxw,minchord;
    double *xwing,*ywing,*zwing,*xflap,*yflap,*zflap,*xaileron,*yaileron,*zaileron,*chordvec,*levec;
    double *xvwing,*yvwing,*zvwing,*xvflap,*yvflap,*zvflap,*xvaileron,*yvaileron,*zvaileron;
    double zplineRt,zplineTp,twistangle,dihedral,xpTS,ypTS,zpTS;
    int *ijwing, *ijflap, *ijaileron;
       
    FILE *fp1;
    int iTS, lstWngTSNumber, cond, ChkAil, ChkWLED, ChkWTED, ChkWGL, WTEDStopPointsNumber, npTS, npTSp1, ndouble;
    int Ailinds[2][2],WTEDinds[2][2],Wglinds[2][2];
    char line[110], code[7], WngArf[12], nindex[12];
    double WngTSLength[4], WngTSRtChord[4], WngTSTpChord[4], WngTSLPosFus[4], WngTSSPosFus[4], WngTSVPosFus[4];
    double WngTSRtInc[4], WngTSTpInc[4], WngTSSwpLE[4], WngTSDhdr[4], WngTSTwist[4], WngTSTR[4], WngTSRtIncZLA[4], WngTSTpIncZLA[4];
    double WngTSSwp25Prct[4], WngTSSwp50Prct[4];
    double WngSpan, WngLPosFus, WngVPosFus, WngRtChord, WngTpChord, WngArea, WngSwpLE, WngTwist, WngDhdrl, WngInc;
    double WngAR, WngTR, WngAReff, WngMAC, WngMACPosX, WngMACPosY, WngMACPosZ, WngMACLPosFus;
    double AilSpan, AilPosSpan, AilArea, AilHingeRlPos,AilRlChord,AilRlSpan;
    double AilMxDDflct, AilMxUDflct, AilDDflct, AilLocSpan, AilArea_Vs_WngArea, AilSpan_Vs_WngSpan;
    double WLEDSpan, WLEDPosSpan, WLEDMxExtdChord, WLEDLocSpan, WLEDRlSpan, WLEDRlChord, WLEDSpan_Vs_WngSpan;
    double WTEDSpan, WTEDPosSpan, WTEDHingeRlPos, WTEDRlChord, WTEDRlSpan, WTEDMxDDflct, WTEDMxUDflct, WTEDLocSpan, WTEDMxExtdChord, WTEDEfficiency;
    double WTEDSpan_Vs_WngSpan, *WTEDStopPoint;
    double WglSpan,WglRtChord,WglTpChord,WglArea,WglSwpLE,WglDhdrl,WglLeOffset,WglTpr;
    
    /* Read Wing.arp and extract wing description */
    printf("Reading from %s\n",Wngfile);
    fp1 = fopen(Wngfile,"r");
    if(fp1 == NULL) {
        fprintf(stderr,"Error:  Could not open %s\n",Wngfile);
        exit(1);
    }
    
    pi=atan(1.0)*4.0;
    
    mp1=m+1;
    np1=n+1;
    ypline= (double *)malloc(sizeof(double)*np1);
    xpline= (double *)malloc(sizeof(double)*mp1);
    chordvec= (double *)malloc(sizeof(double)*np1);
    levec= (double *)malloc(sizeof(double)*np1);
    ycamber = (double *)malloc(sizeof(double)*mp1);
    ycamberwgl = (double *)malloc(sizeof(double)*mp1);
    xpgrid= (double *)malloc(sizeof(double)*mp1*np1);
    ypgrid= (double *)malloc(sizeof(double)*mp1*np1);
    zpgrid= (double *)malloc(sizeof(double)*mp1*np1);
    xv= (double *)malloc(sizeof(double)*mp1*np1);
    yv= (double *)malloc(sizeof(double)*mp1*np1);
    zv= (double *)malloc(sizeof(double)*mp1*np1);
    
    cond=0;
    while (cond == 0){
        fgets(line, 110, fp1);
        if ( strncmp("WNG101",line,6) == 0 ){
            sscanf(line,"%s\t%s",code,WngArf);
        }
        if ( strncmp("WNG401",line,6) == 0 ){
            sscanf(line,"%s\t%s\t%i",code,nindex,&lstWngTSNumber);
            cond=1;
        }
    }
    
    char *token = NULL;
    cond=0;
    while(cond == 0){
        fgets(line, 110, fp1);
        if (  lstWngTSNumber >= 1){
            i=0;
            if ( strncmp("WN1502",line,6) == 0 ){
                sscanf(line,"%s\t%lf\t%lf\t%lf",code,&WngTSLength[i],&WngTSRtChord[i],&WngTSTpChord[i]);
            }
            if ( strncmp("WN1503",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf %lf %lf",code,&WngTSRtInc[i],&WngTSTpInc[i],&WngTSSwpLE[i],&WngTSDhdr[i],&WngTSTwist[i]);
            }
            if ( strncmp("WN1504",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&WngTSTR[i]);
            }
            if ( strncmp("WN1505",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&WngTSRtIncZLA[i],&WngTSTpIncZLA[i]);
            }
            if ( strncmp("WN1603",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&WngTSSwp25Prct[i],&WngTSSwp50Prct[i]);
            }
            if ( strncmp("WN1604",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&WngTSLPosFus[i],&WngTSSPosFus[i],&WngTSVPosFus[i]);
            }
        }
        if (  lstWngTSNumber >= 2){
            i=1;
            if ( strncmp("WN2502",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&WngTSLength[i],&WngTSRtChord[i],&WngTSTpChord[i]);
            }
            if ( strncmp("WN2503",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf %lf %lf",code,&WngTSRtInc[i],&WngTSTpInc[i],&WngTSSwpLE[i],&WngTSDhdr[i],&WngTSTwist[i]);
            }
            if ( strncmp("WN2504",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&WngTSTR[i]);
            }
            if ( strncmp("WN2505",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&WngTSRtIncZLA[i],&WngTSTpIncZLA[i]);
            }
            if ( strncmp("WN2603",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&WngTSSwp25Prct[i],&WngTSSwp50Prct[i]);
            }
            if ( strncmp("WN2604",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&WngTSLPosFus[i],&WngTSSPosFus[i],&WngTSVPosFus[i]);
            }
        }
        if (  lstWngTSNumber >= 3){
            i=2;
            if ( strncmp("WN3502",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&WngTSLength[i],&WngTSRtChord[i],&WngTSTpChord[i]);
            }
            if ( strncmp("WN3503",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf %lf %lf",code,&WngTSRtInc[i],&WngTSTpInc[i],&WngTSSwpLE[i],&WngTSDhdr[i],&WngTSTwist[i]);
            }
            if ( strncmp("WN3504",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&WngTSTR[i]);
            }
            if ( strncmp("WN3505",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&WngTSRtIncZLA[i],&WngTSTpIncZLA[i]);
            }
            if ( strncmp("WN3603",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&WngTSSwp25Prct[i],&WngTSSwp50Prct[i]);
            }
            if ( strncmp("WN3604",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&WngTSLPosFus[i],&WngTSSPosFus[i],&WngTSVPosFus[i]);
            }
        }
        if (  lstWngTSNumber >= 4){
            i=3;          
            if ( strncmp("WN4502",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&WngTSLength[i],&WngTSRtChord[i],&WngTSTpChord[i]);
            }
            if ( strncmp("WN4503",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf %lf %lf",code,&WngTSRtInc[i],&WngTSTpInc[i],&WngTSSwpLE[i],&WngTSDhdr[i],&WngTSTwist[i]);
            }
            if ( strncmp("WN4504",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&WngTSTR[i]);
            }
            if ( strncmp("WN4505",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&WngTSRtIncZLA[i],&WngTSTpIncZLA[i]);
            }
            if ( strncmp("WN4603",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&WngTSSwp25Prct[i],&WngTSSwp50Prct[i]);
            }
            if ( strncmp("WN4604",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&WngTSLPosFus[i],&WngTSSPosFus[i],&WngTSVPosFus[i]);
            }
        }
        if ( strncmp("WNG601",line,6) == 0 ){
            sscanf(line,"%s %lf %lf %lf %lf %lf",code,&WngSpan,&WngLPosFus,&WngVPosFus,&WngRtChord,&WngTpChord);
        }
        if ( strncmp("WNG602",line,6) == 0 ){
            sscanf(line,"%s %lf",code,&WngArea);
        }
        if ( strncmp("WNG604",line,6) == 0 ){
            sscanf(line,"%s %lf %lf %lf %lf",code,&WngSwpLE,&WngTwist,&WngDhdrl,&WngInc);
        }
        if ( strncmp("WNG605",line,6) == 0 ){
            sscanf(line,"%s %lf %lf %lf",code,&WngAR,&WngTR,&WngAReff);
        }
        if ( strncmp("WNG607",line,6) == 0 ){
            sscanf(line,"%s %lf %lf %lf %lf %lf",code,&WngMAC,&WngMACPosX,&WngMACPosY,&WngMACPosZ,&WngMACLPosFus);
        }
        if ( strncmp("AIL201",line,6) == 0 ){
            sscanf(line,"%s %i",code,&ChkAil);
            cond=1;
        }
    }
    cond=0;
    while(cond == 0){
        fgets(line, 110, fp1);
        if ( ChkAil == 1 ){
            if ( strncmp("AIL601",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&AilSpan,&AilPosSpan); /* AilPosSpan is the distance of the of inboard part of the aileron from the wing root?  */
            }
            if ( strncmp("AIL602",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&AilArea); /* AilArea is the area of the two ailerons together? */
            }
            if ( strncmp("AIL603",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&AilHingeRlPos,&AilRlChord,&AilRlSpan); /*AilRlChord is the aileron chord as a percentage of the local chord? */
            }
            if ( strncmp("AIL604",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&AilMxDDflct,&AilMxUDflct,&AilDDflct);
            }
            if ( strncmp("AIL608",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&AilLocSpan); /* No idea what this is */
            }
            if ( strncmp("AIL610",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&AilArea_Vs_WngArea,&AilSpan_Vs_WngSpan); /* AilSpan_Vs_WngSpan is the percentage of the span of one aileron to the total span of the wing? */ 
            }
        }
        if ( strncmp("LED201",line,6) == 0 ){
            sscanf(line,"%s %i",code,&ChkWLED);
            cond=1;
        }
    }
    cond=0;
    while( cond == 0 ){
        fgets(line, 110, fp1);
        if ( ChkWLED == 1 ){
            if ( strncmp("LED601",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&WLEDSpan,&WLEDPosSpan);
            }
            if ( strncmp("LED602",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&WLEDMxExtdChord);
            }
            if ( strncmp("LED604",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&WLEDLocSpan,&WLEDRlSpan,&WLEDRlChord);
            }
            if ( strncmp("LED605",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&WLEDSpan_Vs_WngSpan);
            }
        }
        if ( strncmp("TED201",line,6) == 0 ){
            sscanf(line,"%s %i",code,&ChkWTED);
            cond=1;
        }
    }
    cond=0;
    while( cond == 0 ){
        fgets(line, 110, fp1);
        if ( ChkWTED == 1 ){
            if ( strncmp("TED601",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&WTEDSpan,&WTEDPosSpan);
            }
            if ( strncmp("TED603",line,6) == 0 ){
                token = strtok(line, ":,\t,\n"); /* This puts the line code in token */
                token = strtok(NULL, "\t,\n");
                sscanf(token,"%lf",&WTEDHingeRlPos);
                token = strtok(NULL, "\t,\n");
                sscanf(token,"%lf",&WTEDRlChord);
                token = strtok(NULL, "\t,\n");
                sscanf(token,"%lf",&WTEDRlSpan);
                token = strtok(NULL, "\t,\n");
                sscanf(token,"%lf",&WTEDMxExtdChord);
            }
            if ( strncmp("TED604",line,6) == 0 ){
                sscanf(line,"%s %lf %lf",code,&WTEDMxDDflct,&WTEDMxUDflct);
            }
            if ( strncmp("TED607",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&WTEDLocSpan);
            }
            if ( strncmp("TED608",line,6) == 0 ){
                sscanf(line,"%s %lf",code,&WTEDSpan_Vs_WngSpan);
            }
            if ( strncmp("TED611",line,6) == 0 ){
                sscanf(line,"%s %i",code,&WTEDStopPointsNumber);
                WTEDStopPoint = (double *)malloc(sizeof(double)*WTEDStopPointsNumber);
                fscanf(fp1,"%s",code);
                for (i=0;i<WTEDStopPointsNumber;i++){
                    fscanf(fp1,"%lf",(WTEDStopPoint+i));
                }
            }
        }
        if ( strncmp("WGL201",line,6) == 0 ){
            sscanf(line,"%s %i",code,&ChkWGL);
            cond=1;
        }        
    }
    while(fgets(line, 110, fp1) != NULL){
        if ( ChkWGL == 1 ){
            if ( strncmp("WGL601",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&WglSpan,&WglRtChord,&WglTpChord);
            }
            if ( strncmp("WGL602",line,6) == 0 ){
                sscanf(line,"%s %lf %lf %lf",code,&WglArea,&WglSwpLE,&WglDhdrl);
            }
        }
    }
    
    /* Convert to rad */
    minchord=100; /* Initialize the calculation of the minimum chord of the wing */
    for (i=0;i<lstWngTSNumber;i++){
        WngTSSwpLE[i]=WngTSSwpLE[i]*pi/180.0;
        WngTSDhdr[i]=WngTSDhdr[i]*pi/180.0;
        WngTSTwist[i]=WngTSTwist[i]*pi/180.0;
        WngTSRtInc[i]=WngTSRtInc[i]*pi/180.0;
        WngTSTpInc[i]=WngTSTpInc[i]*pi/180.0;
        if (WngTSTpChord[i] < minchord){
            minchord=WngTSTpChord[i];
        }
    }
    nypos=lstWngTSNumber+5;
    nxpos=4;
    nxyzTS=lstWngTSNumber+1;
    if (  ChkWGL == 1){
        WglLeOffset=0.01;
        WglSwpLE=WglSwpLE*pi/180.0;
        WglDhdrl=WglDhdrl*pi/180.0;
        nypos=lstWngTSNumber+6;
        nxpos=6;
        nxyzTS=lstWngTSNumber+2;
        WglTpr=WglTpChord/WglRtChord; /* Winglet taper ratio */
    }
    
    /* Create vector containing y-coords of TS, Ail, WTED and Wgl */
    ypos= (double *)malloc(sizeof(double)*nypos); 
    xTS= (double *)malloc(sizeof(double)*nxyzTS); 
    yTS= (double *)malloc(sizeof(double)*nxyzTS); 
    yTS2= (double *)malloc(sizeof(double)*nxyzTS); 
    zTS= (double *)malloc(sizeof(double)*nxyzTS); 
    *ypos=0.0;*yTS=0.0;*yTS2=0.0;*xTS=0.0;*zTS=0.0;
    for (i=1;i<lstWngTSNumber+1;i++){
        *(ypos+i)=WngTSLength[i-1]+*(ypos+i-1);
        *(yTS+i)=*(ypos+i);
        *(yTS2+i)=cos(WngTSDhdr[i-1])*WngTSLength[i-1]+*(yTS2+i-1);
        *(xTS+i)=WngTSLength[i-1]*tan(WngTSSwpLE[i-1])+*(xTS+i-1);
        *(zTS+i)=sin(WngTSDhdr[i-1])*WngTSLength[i-1]+*(zTS+i-1);
    }
    *(ypos+i)=WTEDPosSpan;i++;
    *(ypos+i)=WTEDSpan+WTEDPosSpan;i++;
    *(ypos+i)=AilPosSpan;i++;
    *(ypos+i)=AilSpan+AilPosSpan;i++;
    if (  ChkWGL == 1){
        *(ypos+i)=*(yTS+lstWngTSNumber)+WglSpan;
        *(yTS+lstWngTSNumber+1)=*(ypos+i);
        *(yTS2+lstWngTSNumber+1)=WglSpan*cos(WngTSDhdr[i-1])+*(yTS2+lstWngTSNumber);
        *(xTS+lstWngTSNumber+1)=WglSpan*tan(WglSwpLE)+*(xTS+lstWngTSNumber);
        *(zTS+lstWngTSNumber+1)=WglSpan*sin(WglDhdrl)+*(zTS+lstWngTSNumber);
    }                
    /* Sort ypos vector */
    qsort(ypos, nypos, sizeof(double), compare_function);
    /* Remove repeated values from ypos vector */   
    ndouble=checkdoubles(ypos,nypos);
    while (ndouble > 0){
        nypos=finddoubles(ypos,nypos);
        ndouble=checkdoubles(ypos,nypos);
    }
    /* Check if number of spanwise panels is sufficient */
    if (n <= nypos){
        fprintf(stderr,"Error:  The number of spanwise panels on the wing should be larger than %i\n",nypos);
        exit(1);
    }    
    /* Create vector containing the full spanwise grid */
    createypline(ypos,nypos,ypline,np1);
    /* Check between which elements of ypline lie the aileron, flap and winglet */
    WTEDinds[0][1]=findindex(ypline,np1,WTEDPosSpan);
    WTEDinds[1][1]=findindex(ypline,np1,WTEDSpan+WTEDPosSpan);
    Ailinds[0][1]=findindex(ypline,np1,AilPosSpan);
    Ailinds[1][1]=findindex(ypline,np1,AilSpan+AilPosSpan);
    Wglinds[0][1]=findindex(ypline,np1,*(yTS+lstWngTSNumber));
    Wglinds[1][1]=findindex(ypline,np1,*(yTS+lstWngTSNumber)+WglSpan);
    
    /* Create vector containing y-coords of TS, Ail, WTED and Wgl */
    xpos= (double *)malloc(sizeof(double)*nxpos); 
    *(xpos+0)=0.0;
    *(xpos+1)=(100.0-WTEDRlChord)/100.0;
    *(xpos+2)=(100.0-AilRlChord)/100.0;
    *(xpos+3)=1.0;
    if (  ChkWGL == 1){
        *(xpos+4)=WglLeOffset/ WngTSTpChord[lstWngTSNumber-1];
        *(xpos+5)=(WglLeOffset+WglRtChord)/ WngTSTpChord[lstWngTSNumber-1];
    }
    /* Sort xpos vector */
    qsort(xpos, nxpos, sizeof(double), compare_function);
    /* Remove repeated values from ypos vector */   
    ndouble=checkdoubles(xpos,nxpos);
    while (ndouble > 0){
        nxpos=finddoubles(xpos,nxpos);
        ndouble=checkdoubles(xpos,nxpos);
    }
    /* Make sure the trailing edge lies at 1 (it might be 0.994 or something, which can be inaccurate for large chord values) */
    *(xpos+nxpos-1)=1.0;
    /* Check if number of chordwise panels is sufficient */
    if (m <= nxpos){
        fprintf(stderr,"Error:  The number of chordwise panels on the wing should be larger than %i\n",nxpos);
        exit(1);
    }
    /* Create vector containing the full spanwise grid */
    createypline(xpos,nxpos,xpline,mp1);
    /* Check between which elements of ypline lie the aileron, flap and winglet */
    WTEDinds[0][0]=findindex(xpline,mp1,(100.0-WTEDRlChord)/100.0);
    WTEDinds[1][0]=m; /* Will always lie on the trailing edge */
    Ailinds[0][0]=findindex(xpline,mp1,(100.0-AilRlChord)/100.0);
    Ailinds[1][0]=m; /* Will always lie on the trailing edge */
    Wglinds[0][0]=findindex(xpline,mp1,WglLeOffset/ WngTSTpChord[lstWngTSNumber-1]);
    Wglinds[1][0]=findindex(xpline,mp1,(WglLeOffset+WglRtChord)/ WngTSTpChord[lstWngTSNumber-1]);

    /* Create matrix of non-dimensional camber lines */
    rewind(fp1);
    ycamberall = (double *)malloc(sizeof(double)*mp1*lstWngTSNumber*2);
    
    /* Re-read Wing.arp file to find root and tip airfoil names */
    cond=0;
    while(cond == 0){
        fgets(line, 110, fp1);
        if (  lstWngTSNumber >= 1){
            i=0;
            if ( strncmp("WN1501",line,6) == 0 ){
                treatarf(line,xpline,ycamber,ycamberall,mp1,i); /* Airfoils */
           }
        }
        if (  lstWngTSNumber >= 2){
            i=1;
            if ( strncmp("WN2501",line,6) == 0 ){
                treatarf(line,xpline,ycamber,ycamberall,mp1,i); /* Airfoils */
            } 
        }
        if (  lstWngTSNumber >= 3){
            i=2;
            if ( strncmp("WN3501",line,6) == 0 ){
                treatarf(line,xpline,ycamber,ycamberall,mp1,i); /* Airfoils */
            } 
        }
        if (  lstWngTSNumber >= 4){
            i=3;
            if ( strncmp("WN4501",line,6) == 0 ){
                treatarf(line,xpline,ycamber,ycamberall,mp1,i); /* Airfoils */
            } 
        }
        if (  ChkWGL == 1){
            if ( strncmp("WGL102",line,6) == 0 ){
                treatarfwgl(line,xpline,ycamberwgl,mp1); /* Winglet airfoil */
            } 
        }        
        if ( strncmp("WGL601",line,6) == 0 ){
            cond=1;
        }
    }              
    fclose(fp1);
    
    /* Create complete wing grid, split it into wing, flap and aileron later */
    iTS=0;
    for (j=0;j<np1;j++){
        /* Find out which trapezoidal section we're on */
        yhere=*(ypline+j);
        for (i=1;i<nxyzTS-1;i++){
            if (yhere > *(yTS+i))
                iTS=i;
        }
        for (i=0;i<mp1;i++){
            *(ypgrid+i+j*mp1)=yhere;
            if (iTS < lstWngTSNumber){
                *(chordvec+j)=(WngTSTpChord[iTS]-WngTSRtChord[iTS])/WngTSLength[iTS]*(yhere-*(yTS+iTS))+WngTSRtChord[iTS]; /* Local chord length */
                *(levec+j)=*(xTS+iTS)+(yhere-*(yTS+iTS))*tan(WngTSSwpLE[iTS]); /* Local leading edge position */
                *(xpgrid+i+j*mp1)=*(chordvec+j)* *(xpline+i)+ *(levec+j);
                zplineRt=*(ycamberall+i+2*iTS*mp1)*WngTSRtChord[iTS]; /* Root camber line of trapezoidal section */
                zplineTp=*(ycamberall+i+(2*iTS+1)*mp1)*WngTSTpChord[iTS]; /* Tip camber line of trapezoidal section */
                *(zpgrid+i+j*mp1)=(zplineTp-zplineRt)/WngTSLength[iTS]*(yhere-*(yTS+iTS))+zplineRt;
            }else{
                *(chordvec+j)=WngTSTpChord[lstWngTSNumber-1]*(WglTpr-1.0)/WglSpan*(yhere-*(yTS+iTS))+WngTSTpChord[lstWngTSNumber-1]; /* Local chord length */
                *(levec+j)=*(xTS+iTS)+(yhere-*(yTS+iTS))*tan(WglSwpLE); /* Local leading edge position */
                *(xpgrid+i+j*mp1)=*(chordvec+j)* *(xpline+i)+ *(levec+j);
                zplineRt=*(ycamberwgl+i)*WngTSTpChord[lstWngTSNumber-1]; /* Root camber line of trapezoidal section */
                zplineTp=*(ycamberwgl+i)*WngTSTpChord[lstWngTSNumber-1]*WglTpr; /* Tip camber line of trapezoidal section */
                *(zpgrid+i+j*mp1)=(zplineTp-zplineRt)/WglSpan*(yhere-*(yTS+iTS))+zplineRt;
            }
        }
    }    
    
    /* wake shedding distance */
    dxw=0.3*minchord/m;
    vortexpanel(xv,yv,zv,xpgrid,ypgrid,zpgrid,dxw,m,n);  
         
    /* Assign wing grid cells to wing, aileron and flap */
    xwing= (double *)malloc(sizeof(double)*mp1*np1); 
    ywing= (double *)malloc(sizeof(double)*mp1*np1); 
    zwing= (double *)malloc(sizeof(double)*mp1*np1); 
    xvwing= (double *)malloc(sizeof(double)*mp1*np1); 
    yvwing= (double *)malloc(sizeof(double)*mp1*np1); 
    zvwing= (double *)malloc(sizeof(double)*mp1*np1); 
    ijwing= (int *)malloc(sizeof(int)*mp1*np1*2);
    xflap= (double *)malloc(sizeof(double)*mp1*np1); 
    yflap= (double *)malloc(sizeof(double)*mp1*np1); 
    zflap= (double *)malloc(sizeof(double)*mp1*np1); 
    xvflap= (double *)malloc(sizeof(double)*mp1*np1); 
    yvflap= (double *)malloc(sizeof(double)*mp1*np1); 
    zvflap= (double *)malloc(sizeof(double)*mp1*np1); 
    ijflap= (int *)malloc(sizeof(int)*mp1*np1*2);
    xaileron= (double *)malloc(sizeof(double)*mp1*np1); 
    yaileron= (double *)malloc(sizeof(double)*mp1*np1); 
    zaileron= (double *)malloc(sizeof(double)*mp1*np1); 
    xvaileron= (double *)malloc(sizeof(double)*mp1*np1); 
    yvaileron= (double *)malloc(sizeof(double)*mp1*np1); 
    zvaileron= (double *)malloc(sizeof(double)*mp1*np1); 
    ijaileron= (int *)malloc(sizeof(int)*mp1*np1*2);
    wing_nvert=0;flap_nvert=0;aileron_nvert=0; /* Initialize the number of vertices in all the lifting surfaces */
    
    iTS=0;
    for (j=0;j<np1;j++){
        /* Find out which trapezoidal section we're on */
        yhere=*(ypline+j);
        for (i=1;i<nxyzTS-1;i++){
            if (yhere > *(yTS+i))
                iTS=i;
        }
        if (iTS < lstWngTSNumber){
            twistangle=(WngTSTpInc[iTS]-WngTSRtInc[iTS])/WngTSLength[iTS]* (*(ypgrid+j*mp1)-*(yTS+iTS))+WngTSRtInc[iTS];
            dihedral=WngTSDhdr[iTS];
        }else{
            /* Winglet is untwisted and has the same twist angle as wingtip */
            dihedral=WglDhdrl;
        }
        for (i=0;i<mp1;i++){
            flaphere=0;
            winghere=0;
            aileronhere=0;
            /* Check if this point lies on the flap */
            if (i >= WTEDinds[0][0] && i <= WTEDinds[1][0] && j >= WTEDinds[0][1] && j <= WTEDinds[1][1]){ /* Flap */
                flaphere=1;
                if (i == WTEDinds[0][0]){ /* Flap leading edge */
                    winghere=1;
                }
                if (j == WTEDinds[0][1]){ /* Flap inboard edge */
                    winghere=1;
                }
                if (j == WTEDinds[1][1]){ /* Flap outboard edge */
                    winghere=1;
                }
            }
            /* Check if this point lies on the aileron */
            if (i >= Ailinds[0][0] && i <= Ailinds[1][0] && j >= Ailinds[0][1] && j <= Ailinds[1][1]){ /* Aileron */
                aileronhere=1;
                if (i == Ailinds[0][0]){ /* Aileron leading edge */
                    winghere=1;
                }
                if (j == Ailinds[0][1]){ /* Aileron inboard edge */
                    winghere=1;
                }
                if (j == Ailinds[1][1]){ /* Aileron outboard edge */
                    winghere=1;
                }
            }
            /* Check if this point lies on the wiglet */
            if (  ChkWGL == 1){
                if (i >= Wglinds[0][0] && i <= Wglinds[1][0] && j >= Wglinds[0][1] && j <= Wglinds[1][1]){ /* Aileron */
                    winghere=1;
                }
            }
            if (winghere == 1){
                if (aileronhere == 1 && flaphere == 1){ /* If this point belongs to both the aileron and the flap */
                    if (WTEDinds[0][0] > Ailinds[0][0] && i > Ailinds[0][0]){
                        winghere=0;
                    }else if (WTEDinds[0][0] <= Ailinds[0][0] && i > WTEDinds[0][0]){
                        winghere=0;
                    }
                }
                if (flaphere == 1 && WTEDinds[0][1] == 0 && i > WTEDinds[0][0] && j == 0) /* if flap starts at wing root */
                    winghere=0;
                if (aileronhere == 1 && Ailinds[1][1] == np1-1 && i > Ailinds[0][0] && j > Ailinds[0][1]) /* if aileron extends to wingtip */
                    winghere=0;
            }else{
                if (  ChkWGL == 1){
                    if (aileronhere == 0 && flaphere == 0 && j <= Wglinds[0][1]) /* If this point does not lie on the flap, aileron or winglet */
                        winghere=1;
                }else{
                    if (aileronhere == 0 && flaphere == 0) /* If this point lies neither on the flap nor on the aileron */
                        winghere=1;
                }
            }
            /* Impose twist on geometric panels*/
            xpTS=*(xpgrid+i+j*mp1);
            ypTS=*(ypgrid+i+j*mp1);
            zpTS=*(zpgrid+i+j*mp1);
//            twistcentre=*(levec+j)+*(chordvec+j)*0.5; /* The wing twists around the half-chord */
            twistcentre=0.0;
            *(xpgrid+i+j*mp1)=cos(twistangle)*(xpTS-twistcentre)+sin(twistangle)*zpTS+twistcentre;
            *(zpgrid+i+j*mp1)=-sin(twistangle)*(xpTS-twistcentre)+cos(twistangle)*zpTS;
            zpTS=*(zpgrid+i+j*mp1);
            /* Impose dihedral on geometric panels*/
            *(ypgrid+i+j*mp1)=cos(dihedral)*(ypTS-*(yTS+iTS))+*(yTS2+iTS);
            *(zpgrid+i+j*mp1)=sin(dihedral)*(ypTS-*(yTS+iTS))+zpTS+*(zTS+iTS);
            /* Impose twist on vortex panels*/
            xpTS=*(xv+i+j*mp1);
            ypTS=*(yv+i+j*mp1);
            zpTS=*(zv+i+j*mp1);
            *(xv+i+j*mp1)=cos(twistangle)*(xpTS-twistcentre)+sin(twistangle)*zpTS+twistcentre;
            *(zv+i+j*mp1)=-sin(twistangle)*(xpTS-twistcentre)+cos(twistangle)*zpTS;
            zpTS=*(zv+i+j*mp1);
            /* Impose dihedral on vortex panels*/
            *(yv+i+j*mp1)=cos(dihedral)*(ypTS-*(yTS+iTS))+*(yTS2+iTS);
            *(zv+i+j*mp1)=sin(dihedral)*(ypTS-*(yTS+iTS))+zpTS+*(zTS+iTS);
            
            /* Store xpgrid, ypgrid, zpgrid, xv, yv, zv, i and j in the corresponding arrays */
            if (flaphere == 1){
                flap_nvert++;
                *(xflap+flap_nvert-1)=*(xpgrid+i+j*mp1)+WngLPosFus;
                *(yflap+flap_nvert-1)=*(ypgrid+i+j*mp1);
                *(zflap+flap_nvert-1)=*(zpgrid+i+j*mp1)+WngVPosFus;
                *(xvflap+flap_nvert-1)=*(xv+i+j*mp1)+WngLPosFus;
                *(yvflap+flap_nvert-1)=*(yv+i+j*mp1);
                *(zvflap+flap_nvert-1)=*(zv+i+j*mp1)+WngVPosFus;
                *(ijflap+flap_nvert-1)=i;
                *(ijflap+mp1*np1+flap_nvert-1)=j;
            }
            if (aileronhere == 1){
                aileron_nvert++;
                *(xaileron+aileron_nvert-1)=*(xpgrid+i+j*mp1)+WngLPosFus;
                *(yaileron+aileron_nvert-1)=*(ypgrid+i+j*mp1);
                *(zaileron+aileron_nvert-1)=*(zpgrid+i+j*mp1)+WngVPosFus;
                *(xvaileron+aileron_nvert-1)=*(xv+i+j*mp1)+WngLPosFus;
                *(yvaileron+aileron_nvert-1)=*(yv+i+j*mp1);
                *(zvaileron+aileron_nvert-1)=*(zv+i+j*mp1)+WngVPosFus;
                *(ijaileron+aileron_nvert-1)=i;
                *(ijaileron+mp1*np1+aileron_nvert-1)=j;
            }
            if (winghere == 1){
                wing_nvert++;
                *(xwing+wing_nvert-1)=*(xpgrid+i+j*mp1)+WngLPosFus;
                *(ywing+wing_nvert-1)=*(ypgrid+i+j*mp1);
                *(zwing+wing_nvert-1)=*(zpgrid+i+j*mp1)+WngVPosFus;
                *(xvwing+wing_nvert-1)=*(xv+i+j*mp1)+WngLPosFus;
                *(yvwing+wing_nvert-1)=*(yv+i+j*mp1);
                *(zvwing+wing_nvert-1)=*(zv+i+j*mp1)+WngVPosFus;
                *(ijwing+wing_nvert-1)=i;
                *(ijwing+mp1*np1+wing_nvert-1)=j;
            }
        }
    }    
    
    /* Count number of panels on all halves of lifting surfaces */
    flap_nface=countfaces(ijflap,flap_nvert,mp1,np1);
    aileron_nface=countfaces(ijaileron,aileron_nvert,mp1,np1); 
    wing_nface=countfaces(ijwing,wing_nvert,mp1,np1);   
    
    /* Create all the panel information in the flap liftsurf structure */
    pflap->nface=2*flap_nface;
    pflap->nvert=2*flap_nvert;
    pflap->faces=(int *)malloc(sizeof(int)*pflap->nface*4); 
    pflap->shedding=(int *)malloc(sizeof(int)*pflap->nface); 
    arrangefaces(ijflap,mp1,np1,pflap);  

    /* Create all the panel information in the aileron liftsurf structure */
    paileron->nface=2*aileron_nface;
    paileron->nvert=2*aileron_nvert;
    paileron->faces=(int *)malloc(sizeof(int)*paileron->nface*4); 
    paileron->shedding=(int *)malloc(sizeof(int)*paileron->nface); 
    arrangefaces(ijaileron,mp1,np1,paileron);
    
    /* Create all the panel information in the wing liftsurf structure */
    pwing->nface=2*wing_nface;
    pwing->nvert=2*wing_nvert;
    pwing->faces=(int *)malloc(sizeof(int)*pwing->nface*4); 
    pwing->shedding=(int *)malloc(sizeof(int)*pwing->nface); 
    arrangefaces(ijwing,mp1,np1,pwing);

    /* Copy all panel vertex information to the relevant liftsurf structures */
    assignvertices(pflap,xflap,yflap,zflap,xvflap,yvflap,zvflap);
    assignvertices(paileron,xaileron,yaileron,zaileron,xvaileron,yvaileron,zvaileron);
    assignvertices(pwing,xwing,ywing,zwing,xvwing,yvwing,zvwing);     
    
    return WngMAC;
}

void rotateail(struct liftsurf *paileron,double delta[])
{
    /* Rotate ailerons by angle delta */
    int i,j,k,nfaced2,nvertd2,te,nchordpanels,maxj;
    double x1[3],x2[3],uvw[3],T[16];
    double L,sindelta[2],cosdelta[2],sqrtL;
    double u,v,w,a,b,c,u2,v2,w2;
    double x,y,z,xnew,ynew,znew,dx,dy,dz;
    
    nfaced2=paileron->nface/2;
    nvertd2=paileron->nvert/2;
    nchordpanels=paileron->nface/paileron->nshed;
    
    sindelta[0]=sin(delta[0]);
    cosdelta[0]=cos(delta[0]);
    sindelta[1]=sin(-delta[1]);
    cosdelta[1]=cos(-delta[1]);         
    
    te=1; /* 1 if the previous panel was a trailing edge panel */
    for (k=0;k<2;k++){ /* k=0 for left aileron, k=1 for right aileron */        
        for (i=0;i<nfaced2;i++){
            if (te == 1 ){
                /* Find axis of rotation */
                /* The axis is given by the leading edge points of the first panel */
                for (j=0;j<3;j++){
                    x1[j]=*(paileron->vertices+j*paileron->nvert+*(paileron->faces+k*nfaced2+i));
                    x2[j]=*(paileron->vertices+j*paileron->nvert+*(paileron->faces+k*nfaced2+i+3*paileron->nface));
                    uvw[j]=x2[j]-x1[j]; /* Calculate direction vector of rotation axis */
                }
                u=uvw[0];
                v=uvw[1];
                w=uvw[2];
                u2=u*u;
                v2=v*v;
                w2=w*w;
                L=u2+v2+w2;
                sqrtL=sqrt(L);
                a=x1[0];
                b=x1[1];
                c=x1[2];
                /* Calculate how much the leading edge would have moved if it had been rotated */
                x=*(paileron->vertices+0*paileron->nvert+*(paileron->faces+k*nfaced2+i));
                y=*(paileron->vertices+1*paileron->nvert+*(paileron->faces+k*nfaced2+i));
                z=*(paileron->vertices+2*paileron->nvert+*(paileron->faces+k*nfaced2+i));
                xnew=(a*(v2+w2)-u*(b*v+c*w-u*x-v*y-w*z)*(1-cosdelta[k])+L*x*cosdelta[k]+sqrtL*(-c*v+b*w-w*y+v*z)*sindelta[k])/L;
                ynew=(b*(u2+w2)-v*(a*u+c*w-u*x-v*y-w*z)*(1-cosdelta[k])+L*y*cosdelta[k]+sqrtL*(c*u-a*w+w*x-u*z)*sindelta[k])/L;
                znew=(c*(u2+v2)-w*(a*u+b*v-u*x-v*y-w*z)*(1-cosdelta[k])+L*z*cosdelta[k]+sqrtL*(-b*u+a*v-v*x+u*y)*sindelta[k])/L;
                dx=xnew-a;
                dy=ynew-b;
                dz=znew-c;                
            }   
            /* For all panels, only rotate first trailing vertex */
            /* For tip panes also rotate second trailing vertex */
            maxj=2;
            if (i>=nfaced2-nchordpanels)
                maxj=3;
            for (j=1;j<maxj;j++){
                /* Rotate geometric panels */
                x=*(paileron->vertices+0*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface));
                y=*(paileron->vertices+1*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface));
                z=*(paileron->vertices+2*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface));
                /* Equations below from http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ */
                xnew=(a*(v2+w2)-u*(b*v+c*w-u*x-v*y-w*z)*(1-cosdelta[k])+L*x*cosdelta[k]+sqrtL*(-c*v+b*w-w*y+v*z)*sindelta[k])/L;
                ynew=(b*(u2+w2)-v*(a*u+c*w-u*x-v*y-w*z)*(1-cosdelta[k])+L*y*cosdelta[k]+sqrtL*(c*u-a*w+w*x-u*z)*sindelta[k])/L;
                znew=(c*(u2+v2)-w*(a*u+b*v-u*x-v*y-w*z)*(1-cosdelta[k])+L*z*cosdelta[k]+sqrtL*(-b*u+a*v-v*x+u*y)*sindelta[k])/L;
                /* Assign rotated values to paileron->vertices */
                *(paileron->vertices+0*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface))=xnew-dx;
                *(paileron->vertices+1*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface))=ynew-dy;
                *(paileron->vertices+2*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface))=znew-dz;
                /* Rotate vortex panels */
                x=*(paileron->vortex+0*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface));
                y=*(paileron->vortex+1*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface));
                z=*(paileron->vortex+2*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface));
                /* Equations below from http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ */
                xnew=(a*(v2+w2)-u*(b*v+c*w-u*x-v*y-w*z)*(1-cosdelta[k])+L*x*cosdelta[k]+sqrtL*(-c*v+b*w-w*y+v*z)*sindelta[k])/L;
                ynew=(b*(u2+w2)-v*(a*u+c*w-u*x-v*y-w*z)*(1-cosdelta[k])+L*y*cosdelta[k]+sqrtL*(c*u-a*w+w*x-u*z)*sindelta[k])/L;
                znew=(c*(u2+v2)-w*(a*u+b*v-u*x-v*y-w*z)*(1-cosdelta[k])+L*z*cosdelta[k]+sqrtL*(-b*u+a*v-v*x+u*y)*sindelta[k])/L;
                /* Assign rotated values to paileron->vertices */
                *(paileron->vortex+0*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface))=xnew-dx;
                *(paileron->vortex+1*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface))=ynew-dy;
                *(paileron->vortex+2*paileron->nvert+*(paileron->faces+k*nfaced2+i+j*paileron->nface))=znew-dz;
            }
            /* Is this a trailing edge panel? */
            te=*(paileron->shedding+k*nfaced2+i); 
        }
    }
}

void colvec(struct liftsurf *pflap)
{
    /* Calculate the collocation points and the vortex segment lengths on a lifting surface */
    int i,j;
    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    
    /* Calculate collocation points */  
    pflap->control=(double *)malloc(sizeof(double)*pflap->nface*3); 
    pflap->dxy=(double *)malloc(sizeof(double)*pflap->nface*2); 
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
        /* Calculate collocation points */ 
        *(pflap->control+i)=(x1+x4)/8.0+3.0*(x2+x3)/8.0; /* x */
        *(pflap->control+i+pflap->nface)=(y1+y2+y3+y4)/4.0; /* y */
        *(pflap->control+i+2*pflap->nface)=(z1+z4)/8.0+3.0*(z2+z3)/8.0; /* z */
        /* Calculate vortex segment lengths */
        *(pflap->dxy+i+pflap->nface)=sqrt((x4-x1)*(x4-x1)+(y4-y1)*(y4-y1)+(z4-z1)*(z4-z1)); /* vortex segment length in y */
        *(pflap->dxy+i)=sqrt((x3-x4)*(x3-x4)+(y3-y4)*(y3-y4)+(z3-z4)*(z3-z4)); /* vortex segment length in x */
    }
}

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
    int i,j;
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

void vortexblob(double uvw[], double x, double y, double z, double x1, double y1, double z1,
        double x2, double y2, double z2,double gama)
{
    /* Calculate the flow induced at point P(x,y,z) by a vortex blob segment with endpoints P1(x1,y1,z1) and P2(x2,y2,z2) */
    /* Inspired by the code in Katz and Plotkin */
    double rsquared,pi,rcut,r0[3],r1[3],r2[3],ir1,ir2,r1r2[3],square,coeff,blob,epsilon;
    
    pi=3.14159265358979;
    epsilon=0.1; /* Blob radius. This value seems to work very well */
    rcut=1.0e-11; /* Cut-off radius; if distance between point and segment smaller than this, flow is zero */
    /* Calculate vectors r0=P1-P2, r1=P-P1 and r2=P-P2 */
    r0[0]=x2-x1;
    r0[1]=y2-y1;
    r0[2]=z2-z1;
    r1[0]=x-x1;
    r1[1]=y-y1;
    r1[2]=z-z1;
    r2[0]=x-x2;
    r2[1]=y-y2;
    r2[2]=z-z2;
    /* Calculate lengths of r1 and r2 */
    ir1=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]);
    ir2=sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]);
    
    /* Calculate cross product between r1 and r2 */
    cross(r1r2,r1,r2);
    /* Calculate square of length of r1 x r2 */
    square=r1r2[0]*r1r2[0]+r1r2[1]*r1r2[1]+r1r2[2]*r1r2[2];
    /* Calculate blob factor */
    blob=(1.0-exp(-square/(epsilon*epsilon)));
    /* Check if P is too close to the vortex segment */
    if ((ir1 < rcut) || (ir2 < rcut) || (sqrt(square) < rcut)){
        uvw[0]=0.0;
        uvw[1]=0.0;
        uvw[2]=0.0;
    }else{
        /* Calculate flow induced by vortex segment on P */
        coeff=blob*gama/4/pi/square*(r0[0]*(r1[0]/ir1-r2[0]/ir2)+r0[1]*(r1[1]/ir1-r2[1]/ir2)+r0[2]*(r1[2]/ir1-r2[2]/ir2));
        uvw[0]=coeff*r1r2[0];
        uvw[1]=coeff*r1r2[1];
        uvw[2]=coeff*r1r2[2];
    }
}

void vortex(double uvw[], double x, double y, double z, double x1, double y1, double z1,
        double x2, double y2, double z2,double gama)
{
    /* Calculate the flow induced at point P(x,y,z) by a vortex segment with endpoints P1(x1,y1,z1) and P2(x2,y2,z2) */
    /* Inspired by the code in Katz and Plotkin */
    double rsquared,pi,rcut,r0[3],r1[3],r2[3],ir1,ir2,r1r2[3],square,coeff;
    
    pi=3.14159265358979;
    rcut=1.0e-11; /* Cut-off radius; if distance between point and segment smaller than this, flow is zero */
    /* Calculate vectors r0=P1-P2, r1=P-P1 and r2=P-P2 */
    r0[0]=x2-x1;
    r0[1]=y2-y1;
    r0[2]=z2-z1;
    r1[0]=x-x1;
    r1[1]=y-y1;
    r1[2]=z-z1;
    r2[0]=x-x2;
    r2[1]=y-y2;
    r2[2]=z-z2;
    /* Calculate lengths of r1 and r2 */
    ir1=sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]);
    ir2=sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]);
    
    /* Calculate cross product between r1 and r2 */
    cross(r1r2,r1,r2);
    /* Calculate square of length of r1 x r2 */
    square=r1r2[0]*r1r2[0]+r1r2[1]*r1r2[1]+r1r2[2]*r1r2[2];
    /* Check if P is too close to the vortex segment */
    if ((ir1 < rcut) || (ir2 < rcut) || (square < rcut)){
        uvw[0]=0.0;
        uvw[1]=0.0;
        uvw[2]=0.0;
    }else{
        /* Calculate flow induced by vortex segment on P */
        coeff=gama/4/pi/square*(r0[0]*(r1[0]/ir1-r2[0]/ir2)+r0[1]*(r1[1]/ir1-r2[1]/ir2)+r0[2]*(r1[2]/ir1-r2[2]/ir2));
        uvw[0]=coeff*r1r2[0];
        uvw[1]=coeff*r1r2[1];
        uvw[2]=coeff*r1r2[2];
    }
}

void infcoeff(struct liftsurf *plift1, struct liftsurf *plift2, double *AN, double *BN, int istart, int jstart,int mtn)
{
    /* Calculate the matrices of flow influence coefficients */
    /* This is the influence of plift2 on plift1 */
    /* AN is the matrix of normal flow influence coefficients */
    /* BN is the matrix of downwash influence coefficients */
    int i,j;
    double xc,yc,zc;
    double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
    double uvw1[3],uvw2[3],uvw3[3],uvw4[3],uvw[3],uvwstar[3];
    
    for (i=0; i < plift1->nface; i++) {
        /* Control point */
        xc=*(plift1->control+i);
        yc=*(plift1->control+i+plift1->nface);
        zc=*(plift1->control+i+2*plift1->nface);
        /* Vortex rings are quadrilaterals with four vortex segments each */
        for (j=0; j < plift2->nface; j++) {
            /* Vortex panel vertices */
            x1=*(plift2->vortex+*(plift2->faces+j));
            x2=*(plift2->vortex+*(plift2->faces+j+plift2->nface));
            x3=*(plift2->vortex+*(plift2->faces+j+2*plift2->nface));
            x4=*(plift2->vortex+*(plift2->faces+j+3*plift2->nface));
            y1=*(plift2->vortex+plift2->nvert+*(plift2->faces+j));
            y2=*(plift2->vortex+plift2->nvert+*(plift2->faces+j+plift2->nface));
            y3=*(plift2->vortex+plift2->nvert+*(plift2->faces+j+2*plift2->nface));
            y4=*(plift2->vortex+plift2->nvert+*(plift2->faces+j+3*plift2->nface));
            z1=*(plift2->vortex+2*plift2->nvert+*(plift2->faces+j));
            z2=*(plift2->vortex+2*plift2->nvert+*(plift2->faces+j+plift2->nface));
            z3=*(plift2->vortex+2*plift2->nvert+*(plift2->faces+j+2*plift2->nface));
            z4=*(plift2->vortex+2*plift2->nvert+*(plift2->faces+j+3*plift2->nface));
            if (j<plift2->nface/2){
                /* Influence of vortex segment 1 */
                vortex(uvw1,xc,yc,zc,x1,y1,z1,x4,y4,z4,1.0);
                /* Influence of vortex segment 2 */
                vortex(uvw2,xc,yc,zc,x4,y4,z4,x3,y3,z3,1.0);
                /* Influence of vortex segment 3 */
                vortex(uvw3,xc,yc,zc,x3,y3,z3,x2,y2,z2,1.0);
                /* Influence of vortex segment 4 */
                vortex(uvw4,xc,yc,zc,x2,y2,z2,x1,y1,z1,1.0);
            }else{
                /* Influence of vortex segment 1 */
                vortex(uvw1,xc,yc,zc,x4,y4,z4,x1,y1,z1,1.0);
                /* Influence of vortex segment 2 */
                vortex(uvw2,xc,yc,zc,x1,y1,z1,x2,y2,z2,1.0);
                /* Influence of vortex segment 3 */
                vortex(uvw3,xc,yc,zc,x2,y2,z2,x3,y3,z3,1.0);
                /* Influence of vortex segment 4 */
                vortex(uvw4,xc,yc,zc,x3,y3,z3,x4,y4,z4,1.0);
            }
            /* Addd the four influcences to the normal flow */
            uvw[0]=uvw1[0]+uvw2[0]+uvw3[0]+uvw4[0];
            uvw[1]=uvw1[1]+uvw2[1]+uvw3[1]+uvw4[1];
            uvw[2]=uvw1[2]+uvw2[2]+uvw3[2]+uvw4[2];
            /* Addd the two influcences to the downwash */
            uvwstar[0]=uvw2[0]+uvw4[0];
            uvwstar[1]=uvw2[1]+uvw4[1];
            uvwstar[2]=uvw2[2]+uvw4[2];
            /* Create the AN and BN matrices */
            *(AN+istart+i+(jstart+j)*mtn)=uvw[0]* *(plift1->normal+i)+uvw[1]* *(plift1->normal+i+plift1->nface) +
                    uvw[2]* *(plift1->normal+i+2*plift1->nface);
            *(BN+istart+i+(jstart+j)*mtn)=uvwstar[0]* *(plift1->normal+i)+uvwstar[1]* *(plift1->normal+i+plift1->nface) +
                    uvwstar[2]* *(plift1->normal+i+2*plift1->nface);
        }
    }
}

void cycleliftsurf(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder, double *AN, double *BN, int mtn)
{
    int m,n;
    
    /* Cycle through the lifting surfaces in the order wing, flap, aileron, htail, elevator, vtail, rudder */
    /* Treat the vtal and rudder problem as uncoupled */
    m=0;
    n=0;
    /* Influence on the wing from the wing */
    infcoeff(pwing,pwing,AN,BN,m,n,mtn);
    /* Influence on the wing from the flap */
    n+=pwing->nface;
    infcoeff(pwing,pflap,AN,BN,m,n,mtn);
    /* Influence on the wing from the aileron */
    n+=pflap->nface;
    infcoeff(pwing,paileron,AN,BN,m,n,mtn);
    /* Influence on the wing from the horizontal tail */
    n+=paileron->nface;
    infcoeff(pwing,phtail,AN,BN,m,n,mtn);
    /* Influence on the wing from the elevator */
    n+=phtail->nface;
    infcoeff(pwing,pelevator,AN,BN,m,n,mtn);
//     /* Influence on the wing from the vertical tail */
//     n+=pelevator->nface;
//     infcoeff(pwing,pvtail,AN,BN,m,n,mtn);
//     /* Influence on the wing from the rudder */
//     n+=pvtail->nface;
//     infcoeff(pwing,prudder,AN,BN,m,n,mtn);

    m+=pwing->nface;
    n=0;
    /* Influence on the flap from the wing */
    infcoeff(pflap,pwing,AN,BN,m,n,mtn);
    /* Influence on the flap from the flap */
    n+=pwing->nface;
    infcoeff(pflap,pflap,AN,BN,m,n,mtn);
    /* Influence on the flap from the aileron */
    n+=pflap->nface;
    infcoeff(pflap,paileron,AN,BN,m,n,mtn);
    /* Influence on the flap from the horizontal tail */
    n+=paileron->nface;
    infcoeff(pflap,phtail,AN,BN,m,n,mtn);
    /* Influence on the flap from the elevator */
    n+=phtail->nface;
    infcoeff(pflap,pelevator,AN,BN,m,n,mtn);
//     /* Influence on the flap from the vertical tail */
//     n+=pelevator->nface;
//     infcoeff(pflap,pvtail,AN,BN,m,n,mtn);
//     /* Influence on the flap from the rudder */
//     n+=pvtail->nface;
//     infcoeff(pflap,prudder,AN,BN,m,n,mtn);
    
    m+=pflap->nface;
    n=0;
    /* Influence on the aileron from the wing */
    infcoeff(paileron,pwing,AN,BN,m,n,mtn);
    /* Influence on the aileron from the flap */
    n+=pwing->nface;
    infcoeff(paileron,pflap,AN,BN,m,n,mtn);
    /* Influence on the aileron from the aileron */
    n+=pflap->nface;
    infcoeff(paileron,paileron,AN,BN,m,n,mtn);
    /* Influence on the aileron from the horizontal tail */
    n+=paileron->nface;
    infcoeff(paileron,phtail,AN,BN,m,n,mtn);
    /* Influence on the aileron from the elevator */
    n+=phtail->nface;
    infcoeff(paileron,pelevator,AN,BN,m,n,mtn);
//     /* Influence on the aileron from the vertical tail */
//     n+=pelevator->nface;
//     infcoeff(paileron,pvtail,AN,BN,m,n,mtn);
//     /* Influence on the aileron from the rudder */
//     n+=pvtail->nface;
//     infcoeff(paileron,prudder,AN,BN,m,n,mtn);
    
    m+=paileron->nface;
    n=0;
    /* Influence on the horizontal tail from the wing */
    infcoeff(phtail,pwing,AN,BN,m,n,mtn);
    /* Influence on the horizontal tail from the flap */
    n+=pwing->nface;
    infcoeff(phtail,pflap,AN,BN,m,n,mtn);
    /* Influence on the horizontal tail from the aileron */
    n+=pflap->nface;
    infcoeff(phtail,paileron,AN,BN,m,n,mtn);
    /* Influence on the horizontal tail from the horizontal tail */
    n+=paileron->nface;
    infcoeff(phtail,phtail,AN,BN,m,n,mtn);
    /* Influence on the horizontal tail from the elevator */
    n+=phtail->nface;
    infcoeff(phtail,pelevator,AN,BN,m,n,mtn);
//     /* Influence on the horizontal tail from the vertical tail */
//     n+=pelevator->nface;
//     infcoeff(phtail,pvtail,AN,BN,m,n,mtn);
//     /* Influence on the horizontal tail from the rudder */
//     n+=pvtail->nface;
//     infcoeff(phtail,prudder,AN,BN,m,n,mtn);
    
    m+=phtail->nface;
    n=0;
    /* Influence on the elevator from the wing */
    infcoeff(pelevator,pwing,AN,BN,m,n,mtn);
    /* Influence on the elevator from the flap */
    n+=pwing->nface;
    infcoeff(pelevator,pflap,AN,BN,m,n,mtn);
    /* Influence on the elevator from the aileron */
    n+=pflap->nface;
    infcoeff(pelevator,paileron,AN,BN,m,n,mtn);
    /* Influence on the elevator from the horizontal tail */
    n+=paileron->nface;
    infcoeff(pelevator,phtail,AN,BN,m,n,mtn);
    /* Influence on the elevator from the elevator */
    n+=phtail->nface;
    infcoeff(pelevator,pelevator,AN,BN,m,n,mtn);
//     /* Influence on the elevator from the vertical tail */
//     n+=pelevator->nface;
//     infcoeff(pelevator,pvtail,AN,BN,m,n,mtn);
//     /* Influence on the elevator from the rudder */
//     n+=pvtail->nface;
//     infcoeff(pelevator,prudder,AN,BN,m,n,mtn);

    m+=pelevator->nface;
    n=0;
    /* Influence on the vertical tail from the wing */
//    infcoeff(pvtail,pwing,AN,BN,m,n,mtn);
    /* Influence on the vertical tail from the flap */
    n+=pwing->nface;
//    infcoeff(pvtail,pflap,AN,BN,m,n,mtn);
    /* Influence on the vertical tail from the aileron */
    n+=pflap->nface;
//    infcoeff(pvtail,paileron,AN,BN,m,n,mtn);
    /* Influence on the vertical tail from the horizontal tail */
    n+=paileron->nface;
//    infcoeff(pvtail,phtail,AN,BN,m,n,mtn);
    /* Influence on the vertical tail from the elevator */
    n+=phtail->nface;
//    infcoeff(pvtail,pelevator,AN,BN,m,n,mtn);
    /* Influence on the vertical tail from the vertical tail */
    n+=pelevator->nface;
    infcoeff(pvtail,pvtail,AN,BN,m,n,mtn);
    /* Influence on the vertical tail from the rudder */
    n+=pvtail->nface;
    infcoeff(pvtail,prudder,AN,BN,m,n,mtn);
    
    m+=pvtail->nface;
    n=0;
    /* Influence on the rudder from the wing */
//    infcoeff(prudder,pwing,AN,BN,m,n,mtn);
    /* Influence on the rudder from the flap */
    n+=pwing->nface;
//    infcoeff(prudder,pflap,AN,BN,m,n,mtn);
    /* Influence on the rudder from the aileron */
    n+=pflap->nface;
//    infcoeff(prudder,paileron,AN,BN,m,n,mtn);
    /* Influence on the rudder from the horizontal tail */
    n+=paileron->nface;
//    infcoeff(prudder,phtail,AN,BN,m,n,mtn);
    /* Influence on the rudder from the elevator */
    n+=phtail->nface;
//    infcoeff(prudder,pelevator,AN,BN,m,n,mtn);
    /* Influence on the rudder from the vertical tail */
    n+=pelevator->nface;
    infcoeff(prudder,pvtail,AN,BN,m,n,mtn);
    /* Influence on the rudder from the rudder */
    n+=pvtail->nface;
    infcoeff(prudder,prudder,AN,BN,m,n,mtn);
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

void calcRHS(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder, double *RHS, double UVW[])
{   
    /* Calculate the Right Hand Side vector */
    /* Cycle through the lifting surfaces in the order wing, flap, aileron, horizontal tail, elevator */
    
    int i,m;

    m=0;
    for (i=0; i < pwing->nface; i++){
        *(RHS+i)=-(*(pwing->uvw+i)+UVW[0])* *(pwing->normal+i)-(*(pwing->uvw+i+pwing->nface)+UVW[1])* *(pwing->normal+i+pwing->nface) -
                (*(pwing->uvw+i+2*pwing->nface)+UVW[2])* *(pwing->normal+i+2*pwing->nface);
    }
    m+=pwing->nface;
    for (i=0; i < pflap->nface; i++){
        *(RHS+i+m)=-(*(pflap->uvw+i)+UVW[0])* *(pflap->normal+i)-(*(pflap->uvw+i+pflap->nface)+UVW[1])* *(pflap->normal+i+pflap->nface) -
                (*(pflap->uvw+i+2*pflap->nface)+UVW[2])* *(pflap->normal+i+2*pflap->nface);
    }
    m+=pflap->nface;
    for (i=0; i < paileron->nface; i++){
        *(RHS+i+m)=-(*(paileron->uvw+i)+UVW[0])* *(paileron->normal+i)-(*(paileron->uvw+i+paileron->nface)+UVW[1])* *(paileron->normal+i+paileron->nface) -
                (*(paileron->uvw+i+2*paileron->nface)+UVW[2])* *(paileron->normal+i+2*paileron->nface);
    }
    m+=paileron->nface;
    for (i=0; i < phtail->nface; i++){
        *(RHS+i+m)=-(*(phtail->uvw+i)+UVW[0])* *(phtail->normal+i)-(*(phtail->uvw+i+phtail->nface)+UVW[1])* *(phtail->normal+i+phtail->nface) -
                (*(phtail->uvw+i+2*phtail->nface)+UVW[2])* *(phtail->normal+i+2*phtail->nface);
    }
    m+=phtail->nface;
    for (i=0; i < pelevator->nface; i++){
        *(RHS+i+m)=-(*(pelevator->uvw+i)+UVW[0])* *(pelevator->normal+i)-(*(pelevator->uvw+i+pelevator->nface)+UVW[1])* *(pelevator->normal+i+pelevator->nface) -
                (*(pelevator->uvw+i+2*pelevator->nface)+UVW[2])* *(pelevator->normal+i+2*pelevator->nface);
    }
    m+=pelevator->nface;
    for (i=0; i < pvtail->nface; i++){
        *(RHS+i+m)=-(*(pvtail->uvw+i)+UVW[0])* *(pvtail->normal+i)-(*(pvtail->uvw+i+pvtail->nface)+UVW[1])* *(pvtail->normal+i+pvtail->nface) -
                (*(pvtail->uvw+i+2*pvtail->nface)+UVW[2])* *(pvtail->normal+i+2*pvtail->nface);
    }
    m+=pvtail->nface;
    for (i=0; i < prudder->nface; i++){
        *(RHS+i+m)=-(*(prudder->uvw+i)+UVW[0])* *(prudder->normal+i)-(*(prudder->uvw+i+prudder->nface)+UVW[1])* *(prudder->normal+i+prudder->nface) -
                (*(prudder->uvw+i+2*prudder->nface)+UVW[2])* *(prudder->normal+i+2*prudder->nface);
    }
}

void assignGammas(double *Gammas, struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder, int it)
{
    /* Assign the vortex strengths in vector Gammas back onto the panels of all the lifting surfaces */
    int i,m;
        
    /* Assign wing vortex strengths */
    m=0;
    for (i=0;i<pwing->nface;i++){
        *(pwing->gamma+i+it*pwing->nface)=*(Gammas+i);
    }
    /* Assign flap vortex strengths */
    m+=pwing->nface;
    for (i=0;i<pflap->nface;i++){
        *(pflap->gamma+i+it*pflap->nface)=*(Gammas+i+m);
    }
    /* Assign aileron vortex strengths */
    m+=pflap->nface;
    for (i=0;i<paileron->nface;i++){
        *(paileron->gamma+i+it*paileron->nface)=*(Gammas+i+m);
    }
    m+=paileron->nface;
    /* Assign horizontal tail vortex strengths */
    for (i=0;i<phtail->nface;i++){
        *(phtail->gamma+i+it*phtail->nface)=*(Gammas+i+m);
    }  
    m+=phtail->nface;
    /* Assign elevator vortex strengths */
    for (i=0;i<pelevator->nface;i++){
        *(pelevator->gamma+i+it*pelevator->nface)=*(Gammas+i+m);
    }
    m+=pelevator->nface;
    /* Assign vertical tail vortex strengths */
    for (i=0;i<pvtail->nface;i++){
        *(pvtail->gamma+i+it*pvtail->nface)=*(Gammas+i+m);
    }
    m+=pvtail->nface;
    /* Assign rudder vortex strengths */
    for (i=0;i<prudder->nface;i++){
        *(prudder->gamma+i+it*prudder->nface)=*(Gammas+i+m);
    }    
}

void assignwind(double *wind, struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder)
{
    /* Assign the vortex strengths in vector Gammas back onto the panels of all the lifting surfaces */
    int i,m;
        
    /* Assign wing vortex strengths */
    m=0;
    for (i=0;i<pwing->nface;i++){
        *(pwing->wind+i)=*(wind+i);
    }
    /* Assign flap vortex strengths */
    m+=pwing->nface;
    for (i=0;i<pflap->nface;i++){
        *(pflap->wind+i)=*(wind+i+m);
    }
    /* Assign aileron vortex strengths */
    m+=pflap->nface;
    for (i=0;i<paileron->nface;i++){
        *(paileron->wind+i)=*(wind+i+m);
    }
    m+=paileron->nface;
    /* Assign horizontal tail vortex strengths */
    for (i=0;i<phtail->nface;i++){
        *(phtail->wind+i)=*(wind+i+m);
    }  
    m+=phtail->nface;
    /* Assign elevator vortex strengths */
    for (i=0;i<pelevator->nface;i++){
        *(pelevator->wind+i)=*(wind+i+m);
    }
    m+=pelevator->nface;
    /* Assign vertical tail vortex strengths */
    for (i=0;i<pvtail->nface;i++){
        *(pvtail->wind+i)=*(wind+i+m);
    }  
    m+=pvtail->nface;
    /* Assign rudder vortex strengths */
    for (i=0;i<prudder->nface;i++){
        *(prudder->wind+i)=*(wind+i+m);
    }      
}

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

void createwake(struct liftsurf *pwing, double cwing,int ntimes)
{
    /* Create wakes shed from a liftsurf */
    int i,nwakes;
    
    if (pwing->nshed !=0 ){
        pwing->wakeinds=(int *)malloc(sizeof(int)*pwing->nshed);
        storewakeinds(pwing);
    }
    setupwakes(pwing,cwing);
    pwing->xw=(double *)malloc(sizeof(double)*(pwing->nshed+pwing->nwakes)*ntimes);
    pwing->yw=(double *)malloc(sizeof(double)*(pwing->nshed+pwing->nwakes)*ntimes);
    pwing->zw=(double *)malloc(sizeof(double)*(pwing->nshed+pwing->nwakes)*ntimes);
    pwing->uw=(double *)malloc(sizeof(double)*(pwing->nshed+pwing->nwakes)*ntimes);
    pwing->vw=(double *)malloc(sizeof(double)*(pwing->nshed+pwing->nwakes)*ntimes);
    pwing->ww=(double *)malloc(sizeof(double)*(pwing->nshed+pwing->nwakes)*ntimes);
    pwing->gw=(double *)malloc(sizeof(double)*pwing->nshed*ntimes);
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

void calcdistwake(struct liftsurf *pwing, struct liftsurf *pflap, int i, int cflap)
{
    /* Calculate the distance between vertex i in the first shed element of the wake in pwing and
     all the vertices in the first shed element of the wake in pflap */
    
    int j;
    double distx,disty,distz,dist;
    
    for (j=0;j<pflap->nshed+pflap->nwakes;j++){
        distx=(*(pwing->xw+i)-*(pflap->xw+j));
        disty=(*(pwing->yw+i)-*(pflap->yw+j));
        distz=(*(pwing->zw+i)-*(pflap->zw+j));
        dist=distx*distx+disty*disty+distz*distz;
        if (dist < 0.0000001){
            *(pwing->wakecorrespwake+i)=j+cflap;
            if (dist > 0){
                *(pflap->xw+j)=*(pwing->xw+i);
                *(pflap->yw+j)=*(pwing->yw+i);
                *(pflap->zw+j)=*(pwing->zw+i);
            }
        }
    }
}

void findallwakeneighbours(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, int cwing, int cflap, int caileron)
{
    /* Check first shed wake element to see if the different lifsurfs have coincident points */
    int i,j;
    double dist;
    
    pwing->wakecorrespwake=(int *)malloc(sizeof(int)*(pwing->nshed+pwing->nwakes));
    pflap->wakecorrespwake=(int *)malloc(sizeof(int)*(pflap->nshed+pflap->nwakes));
    paileron->wakecorrespwake=(int *)malloc(sizeof(int)*(paileron->nshed+paileron->nwakes));
    
    /* Wing */
    for (i=0;i<pwing->nshed+pwing->nwakes;i++){
        *(pwing->wakecorrespwake+i)=-1;
        /* Flap */
        calcdistwake(pwing,pflap,i,cflap);
        /* Aileron */
        calcdistwake(pwing,paileron,i,caileron);
    }
    /* Flap */
    for (i=0;i<pflap->nshed+pflap->nwakes;i++){
        *(pflap->wakecorrespwake+i)=-1;
        /* Wing */
        calcdistwake(pflap,pwing,i,cwing);
        /* Aileron */
        calcdistwake(pflap,paileron,i,caileron);
    }  
    /* Aileron */
    for (i=0;i<paileron->nshed+paileron->nwakes;i++){
        *(paileron->wakecorrespwake+i)=-1;
        /* Wing */
        calcdistwake(paileron,pwing,i,cwing);
        /* Flap */
        calcdistwake(paileron,pflap,i,cflap);
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

void calcwakeinf(struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *phtail, struct liftsurf *pelevator, struct liftsurf *pvtail, struct liftsurf *prudder, int it)
{
    /* Calculate the influence of every wake on every liftsurf */
    /* Treat the vtail and elevator as uncoupled */
    int i;
    
    /* Influence on the wing */
    /* Influence of the wing wake */
    for (i=0;i<pwing->nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(pwing,pwing,i,it,i);
    }
    /* Influence of the flap wake */
    for (i=0;i<pflap->nwakes;i++){
        wakeinf(pwing,pflap,i,it,1);
    }
    /* Influence of the aileron wake */
    for (i=0;i<paileron->nwakes;i++){
        wakeinf(pwing,paileron,i,it,1);
    }
    /* Influence of the horizontal tail wake */
    for (i=0;i<phtail->nwakes;i++){
        wakeinf(pwing,phtail,i,it,1);
    }
    /* Influence of the elevator wake */
    for (i=0;i<pelevator->nwakes;i++){
        wakeinf(pwing,pelevator,i,it,1);
    } 
//     /* Influence of the vertical tail wake */
//     for (i=0;i<pvtail->nwakes;i++){
//         wakeinf(pwing,pvtail,i,it,1);
//     }
//     /* Influence of the rudder wake */
//     for (i=0;i<prudder->nwakes;i++){
//         wakeinf(pwing,prudder,i,it,1);
//     }     
    /* Influence on the flap */
    /* Influence of the wing wake */
    for (i=0;i<pwing->nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(pflap,pwing,i,it,i);
    }
    /* Influence of the flap wake */
    for (i=0;i<pflap->nwakes;i++){
        wakeinf(pflap,pflap,i,it,1);
    }
    /* Influence of the aileron wake */
    for (i=0;i<paileron->nwakes;i++){
        wakeinf(pflap,paileron,i,it,1);
    } 
    /* Influence of the horizontal tail wake */
    for (i=0;i<phtail->nwakes;i++){
        wakeinf(pflap,phtail,i,it,1);
    } 
    /* Influence of the elevator wake */
    for (i=0;i<pelevator->nwakes;i++){
        wakeinf(pflap,pelevator,i,it,1);
    } 
//     /* Influence of the vertical tail wake */
//     for (i=0;i<pvtail->nwakes;i++){
//         wakeinf(pflap,pvtail,i,it,1);
//     }
//     /* Influence of the rudder wake */
//     for (i=0;i<prudder->nwakes;i++){
//         wakeinf(pflap,prudder,i,it,1);
//     }     
    /* Influence on the aileron */
    /* Influence of the wing wake */
    for (i=0;i<pwing->nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(paileron,pwing,i,it,i);
    }
    /* Influence of the flap wake */
    for (i=0;i<pflap->nwakes;i++){
        wakeinf(paileron,pflap,i,it,1);
    }
    /* Influence of the aileron wake */
    for (i=0;i<paileron->nwakes;i++){
        wakeinf(paileron,paileron,i,it,1);
    }  
    /* Influence of the horizontal tail wake */
    for (i=0;i<phtail->nwakes;i++){
        wakeinf(paileron,phtail,i,it,1);
    }
    /* Influence of the elevator wake */
    for (i=0;i<pelevator->nwakes;i++){
        wakeinf(paileron,pelevator,i,it,1);
    }  
//     /* Influence of the vertical tail wake */
//     for (i=0;i<pvtail->nwakes;i++){
//         wakeinf(paileron,pvtail,i,it,1);
//     }
//     /* Influence of the rudder wake */
//     for (i=0;i<prudder->nwakes;i++){
//         wakeinf(paileron,prudder,i,it,1);
//     } 
    /* Influence on the horizontal tail */
    /* Influence of the wing wake */
    for (i=0;i<pwing->nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(phtail,pwing,i,it,i);
    }
    /* Influence of the flap wake */
    for (i=0;i<pflap->nwakes;i++){
        wakeinf(phtail,pflap,i,it,1);
    }
    /* Influence of the aileron wake */
    for (i=0;i<paileron->nwakes;i++){
        wakeinf(phtail,paileron,i,it,1);
    }  
    /* Influence of the horizontal tail wake */
    for (i=0;i<phtail->nwakes;i++){
        wakeinf(phtail,phtail,i,it,1);
    }
    /* Influence of the elevator wake */
    for (i=0;i<pelevator->nwakes;i++){
        wakeinf(phtail,pelevator,i,it,1);
    } 
//     /* Influence of the vertical tail wake */
//     for (i=0;i<pvtail->nwakes;i++){
//         wakeinf(phtail,pvtail,i,it,1);
//     }
//     /* Influence of the rudder wake */
//     for (i=0;i<prudder->nwakes;i++){
//         wakeinf(phtail,prudder,i,it,1);
//     }    
    /* Influence on the elevator tail */
    /* Influence of the wing wake */
    for (i=0;i<pwing->nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(pelevator,pwing,i,it,i);
    }
    /* Influence of the flap wake */
    for (i=0;i<pflap->nwakes;i++){
        wakeinf(pelevator,pflap,i,it,1);
    }
    /* Influence of the aileron wake */
    for (i=0;i<paileron->nwakes;i++){
        wakeinf(pelevator,paileron,i,it,1);
    }  
    /* Influence of the horizontal tail wake */
    for (i=0;i<phtail->nwakes;i++){
        wakeinf(pelevator,phtail,i,it,1);
    }
    /* Influence of the elevator wake */
    for (i=0;i<pelevator->nwakes;i++){
        wakeinf(pelevator,pelevator,i,it,1);
    }     
//     /* Influence of the vertical tail wake */
//     for (i=0;i<pvtail->nwakes;i++){
//         wakeinf(pelevator,pvtail,i,it,1);
//     }
//     /* Influence of the rudder wake */
//     for (i=0;i<prudder->nwakes;i++){
//         wakeinf(pelevator,prudder,i,it,1);
//     } 
    /* Influence on the vertical tail */
//     /* Influence of the wing wake */
//     for (i=0;i<pwing->nwakes;i++){
//         /* We set summode to i so summode=0 the first time we calculate an influence */
//         wakeinf(pvtail,pwing,i,it,i);
//     }
//     /* Influence of the flap wake */
//     for (i=0;i<pflap->nwakes;i++){
//         wakeinf(pvtail,pflap,i,it,1);
//     }
//     /* Influence of the aileron wake */
//     for (i=0;i<paileron->nwakes;i++){
//         wakeinf(pvtail,paileron,i,it,1);
//     }  
//     /* Influence of the horizontal tail wake */
//     for (i=0;i<phtail->nwakes;i++){
//         wakeinf(pvtail,phtail,i,it,1);
//     }
//     /* Influence of the elevator wake */
//     for (i=0;i<pelevator->nwakes;i++){
//         wakeinf(pvtail,pelevator,i,it,1);
//     } 
    /* Influence of the vertical tail wake */
    for (i=0;i<pvtail->nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(pvtail,pvtail,i,it,i);
    }
    /* Influence of the rudder wake */
    for (i=0;i<prudder->nwakes;i++){
        wakeinf(pvtail,prudder,i,it,1);
    }    
    /* Influence on the rudder */
//     /* Influence of the wing wake */
//     for (i=0;i<pwing->nwakes;i++){
//         /* We set summode to i so summode=0 the first time we calculate an influence */
//         wakeinf(prudder,pwing,i,it,i);
//     }
//     /* Influence of the flap wake */
//     for (i=0;i<pflap->nwakes;i++){
//         wakeinf(prudder,pflap,i,it,1);
//     }
//     /* Influence of the aileron wake */
//     for (i=0;i<paileron->nwakes;i++){
//         wakeinf(prudder,paileron,i,it,1);
//     }  
//     /* Influence of the horizontal tail wake */
//     for (i=0;i<phtail->nwakes;i++){
//         wakeinf(prudder,phtail,i,it,1);
//     }
//     /* Influence of the elevator wake */
//     for (i=0;i<pelevator->nwakes;i++){
//         wakeinf(prudder,pelevator,i,it,1);
//     } 
    /* Influence of the vertical tail wake */
    for (i=0;i<pvtail->nwakes;i++){
        /* We set summode to i so summode=0 the first time we calculate an influence */
        wakeinf(prudder,pvtail,i,it,i);
    }
    /* Influence of the rudder wake */
    for (i=0;i<prudder->nwakes;i++){
        wakeinf(prudder,prudder,i,it,1);
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

void calcforces(struct liftsurf *plift, struct liftsurf *pwing, struct liftsurf *pflap, struct liftsurf *paileron, int it, double dt, double UVW[], double rho, int cwing, int cflap, int caileron)
{
    /* Calculate aerodynamic forces acting on all liftsurfs */
    int i,neigh;
    double dGdt, Ux, Uy, Ui, Gx, Gy;
    
    *(plift->aeroforce)=0.0;
    *(plift->aeroforce+1)=0.0;
    *(plift->aeroforce+2)=0.0;
    *(plift->aeroforce+3)=0.0;
    
    for (i=0; i < plift->nface; i++){
        /* Calculate change of vorticity in time */
        if (it == 0){
            dGdt=*(plift->gamma+i+it*plift->nface)/dt;
        }else{
            dGdt=(*(plift->gamma+i+it*plift->nface)-*(plift->gamma+i+(it-1)*plift->nface))/dt;
        }
        
        Ux=(*(plift->uvw+i)+UVW[0])* *(plift->tangx+i)+
                (*(plift->uvw+i+plift->nface)+UVW[1])* *(plift->tangx+i+plift->nface)+
                (*(plift->uvw+i+2*plift->nface)+UVW[2])* *(plift->tangx+i+2*plift->nface);
        Uy=(*(plift->uvw+i)+UVW[0])* *(plift->tangy+i)+
                (*(plift->uvw+i+plift->nface)+UVW[1])* *(plift->tangy+i+plift->nface)+
                (*(plift->uvw+i+2*plift->nface)+UVW[2])* *(plift->tangy+i+2*plift->nface);
        Ui= *(plift->wind+i)+ *(plift->uvw+i+2*plift->nface);

        /* Calculate lift contribution of each panel */
        /* Separate calculations for each half-wing */
        
        neigh=*(plift->neighbours+i+2*plift->nface);
        Gx=0;
        Gy=0;
        if (neigh == -1){
            Gx=*(plift->gamma+i+it*plift->nface);
        }else if (neigh < cflap && neigh > cwing-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(pwing->gamma+neigh-cwing+it*pwing->nface);
        }else if (neigh < caileron && neigh > cflap-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(pflap->gamma+neigh-cflap+it*pflap->nface);
        }else if (neigh > caileron-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(paileron->gamma+neigh-caileron+it*paileron->nface);
        }      
        if (i < plift->nface/2){
            neigh=*(plift->neighbours+i);
            if (neigh == -1){
                Gy=*(plift->gamma+i+it*plift->nface);
            }else if (neigh < cflap && neigh > cwing-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pwing->gamma+neigh-cwing+it*pwing->nface);
            }else if (neigh < caileron && neigh > cflap-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pflap->gamma+neigh-cflap+it*pflap->nface);
            }else if (neigh > caileron-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(paileron->gamma+neigh-caileron+it*paileron->nface);
            }
            *(plift->Deltap+i)=rho*(Ux*Gx/ *(plift->dxy+i)+Uy*Gy/ *(plift->dxy+i+plift->nface)+dGdt);
        }else{
            neigh=*(plift->neighbours+i+plift->nface);
            if (neigh == -1){
                Gy=*(plift->gamma+i+it*plift->nface);
            }else if (neigh < cflap && neigh > cwing-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pwing->gamma+neigh-cwing+it*pwing->nface);
            }else if (neigh < caileron && neigh > cflap-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pflap->gamma+neigh-cflap+it*pflap->nface);
            }else if (neigh > caileron-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(paileron->gamma+neigh-caileron+it*paileron->nface);
            }
            *(plift->Deltap+i)=rho*(Ux*Gx/ *(plift->dxy+i)-Uy*Gy/ *(plift->dxy+i+plift->nface)+dGdt);
        }
        /* Calculate induced drag contribution of each panel */
        *(plift->Deltad+i)=rho*(-Ui*Gx* *(plift->dxy+i+plift->nface)+dGdt* *(plift->nsurf+i));
        
        *(plift->aeroforce) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i);
        *(plift->aeroforce+1) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i+1*plift->nface);
        *(plift->aeroforce+2) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i+2*plift->nface);
        *(plift->aeroforce+3) += *(plift->Deltad+i);
    }
}

void calcforceshtail(struct liftsurf *plift, struct liftsurf *phtail, struct liftsurf *pelevator, int it, double dt, double UVW[], double rho, int chtail, int celevator)
{
    /* Calculate aerodynamic forces acting on all liftsurfs */
    int i,neigh;
    double dGdt, Ux, Uy, Ui, Gx, Gy;

    *(plift->aeroforce)=0.0;
    *(plift->aeroforce+1)=0.0;
    *(plift->aeroforce+2)=0.0;
    *(plift->aeroforce+3)=0.0;

    for (i=0; i < plift->nface; i++){
        /* Calculate change of vorticity in time */
        if (it == 0){
            dGdt=*(plift->gamma+i+it*plift->nface)/dt;
        }else{
            dGdt=(*(plift->gamma+i+it*plift->nface)-*(plift->gamma+i+(it-1)*plift->nface))/dt;
        }

        Ux=(*(plift->uvw+i)+UVW[0])* *(plift->tangx+i)+
                (*(plift->uvw+i+plift->nface)+UVW[1])* *(plift->tangx+i+plift->nface)+
                (*(plift->uvw+i+2*plift->nface)+UVW[2])* *(plift->tangx+i+2*plift->nface);
        Uy=(*(plift->uvw+i)+UVW[0])* *(plift->tangy+i)+
                (*(plift->uvw+i+plift->nface)+UVW[1])* *(plift->tangy+i+plift->nface)+
                (*(plift->uvw+i+2*plift->nface)+UVW[2])* *(plift->tangy+i+2*plift->nface);
        Ui= *(plift->wind+i)+ *(plift->uvw+i+2*plift->nface);

        /* Calculate lift contribution of each panel */
        /* Separate calculations for each half-wing */

        neigh=*(plift->neighbours+i+2*plift->nface);
        Gx=0;
        Gy=0;
        if (neigh < chtail){
            Gx=*(plift->gamma+i+it*plift->nface);
        }else if (neigh < celevator && neigh > chtail-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(phtail->gamma+neigh-chtail+it*phtail->nface);
        }else if (neigh > celevator-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(pelevator->gamma+neigh-celevator+it*pelevator->nface);
        }
        if (i < plift->nface/2){
            neigh=*(plift->neighbours+i);
            if (neigh < chtail){
                Gy=*(plift->gamma+i+it*plift->nface);
            }else if (neigh < celevator && neigh > chtail-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(phtail->gamma+neigh-chtail+it*phtail->nface);
            }else if (neigh > celevator-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pelevator->gamma+neigh-celevator+it*pelevator->nface);
            }
            *(plift->Deltap+i)=rho*(Ux*Gx/ *(plift->dxy+i)+Uy*Gy/ *(plift->dxy+i+plift->nface)+dGdt);
        }else{
            neigh=*(plift->neighbours+i+plift->nface);
            if (neigh < chtail){
                Gy=*(plift->gamma+i+it*plift->nface);
            }else if (neigh < celevator && neigh > chtail-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(phtail->gamma+neigh-chtail+it*phtail->nface);
            }else if (neigh > celevator-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pelevator->gamma+neigh-celevator+it*pelevator->nface);
            }
            *(plift->Deltap+i)=rho*(Ux*Gx/ *(plift->dxy+i)-Uy*Gy/ *(plift->dxy+i+plift->nface)+dGdt);
        }
        /* Calculate induced drag contribution of each panel */
        *(plift->Deltad+i)=rho*(-Ui*Gx* *(plift->dxy+i+plift->nface)+dGdt* *(plift->nsurf+i));

        *(plift->aeroforce) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i);
        *(plift->aeroforce+1) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i+1*plift->nface);
        *(plift->aeroforce+2) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i+2*plift->nface);
        *(plift->aeroforce+3) += *(plift->Deltad+i);
    }
}

void calcforcesvtail(struct liftsurf *plift, struct liftsurf *pvtail, struct liftsurf *prudder, int it, double dt, double UVW[], double rho, int cvtail, int crudder)
{
    /* Calculate aerodynamic forces acting on all liftsurfs */
    int i,neigh;
    double dGdt, Ux, Uy, Ui, Gx, Gy;

    *(plift->aeroforce)=0.0;
    *(plift->aeroforce+1)=0.0;
    *(plift->aeroforce+2)=0.0;
    *(plift->aeroforce+3)=0.0;

    for (i=0; i < plift->nface; i++){
        /* Calculate change of vorticity in time */
        if (it == 0){
            dGdt=*(plift->gamma+i+it*plift->nface)/dt;
        }else{
            dGdt=(*(plift->gamma+i+it*plift->nface)-*(plift->gamma+i+(it-1)*plift->nface))/dt;
        }

        Ux=(*(plift->uvw+i)+UVW[0])* *(plift->tangx+i)+
                (*(plift->uvw+i+plift->nface)+UVW[1])* *(plift->tangx+i+plift->nface)+
                (*(plift->uvw+i+2*plift->nface)+UVW[2])* *(plift->tangx+i+2*plift->nface);
        Uy=(*(plift->uvw+i)+UVW[0])* *(plift->tangy+i)+
                (*(plift->uvw+i+plift->nface)+UVW[1])* *(plift->tangy+i+plift->nface)+
                (*(plift->uvw+i+2*plift->nface)+UVW[2])* *(plift->tangy+i+2*plift->nface);
        Ui= *(plift->wind+i)+ *(plift->uvw+i+2*plift->nface);

        /* Calculate lift contribution of each panel */
        /* Separate calculations for each half-wing */

        neigh=*(plift->neighbours+i+2*plift->nface);
        Gx=0;
        Gy=0;
        if (neigh < cvtail){
            Gx=*(plift->gamma+i+it*plift->nface);
        }else if (neigh < crudder && neigh > cvtail-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(pvtail->gamma+neigh-cvtail+it*pvtail->nface);
        }else if (neigh > crudder-1){
            Gx=*(plift->gamma+i+it*plift->nface)- *(prudder->gamma+neigh-crudder+it*prudder->nface);
        }
        if (i < plift->nface/2){
            neigh=*(plift->neighbours+i);
            if (neigh < cvtail){
                Gy=*(plift->gamma+i+it*plift->nface);
            }else if (neigh < crudder && neigh > cvtail-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pvtail->gamma+neigh-cvtail+it*pvtail->nface);
            }else if (neigh > crudder-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(prudder->gamma+neigh-crudder+it*prudder->nface);
            }
            *(plift->Deltap+i)=rho*(Ux*Gx/ *(plift->dxy+i)+Uy*Gy/ *(plift->dxy+i+plift->nface)+dGdt);
        }else{
            neigh=*(plift->neighbours+i+plift->nface);
            if (neigh < cvtail){
                Gy=*(plift->gamma+i+it*plift->nface);
            }else if (neigh < crudder && neigh > cvtail-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(pvtail->gamma+neigh-cvtail+it*pvtail->nface);
            }else if (neigh > crudder-1){
                Gy=*(plift->gamma+i+it*plift->nface)- *(prudder->gamma+neigh-crudder+it*prudder->nface);
            }
            *(plift->Deltap+i)=rho*(Ux*Gx/ *(plift->dxy+i)-Uy*Gy/ *(plift->dxy+i+plift->nface)+dGdt);
        }
        /* Calculate induced drag contribution of each panel */
        *(plift->Deltad+i)=rho*(-Ui*Gx* *(plift->dxy+i+plift->nface)+dGdt* *(plift->nsurf+i));

        *(plift->aeroforce) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i);
        *(plift->aeroforce+1) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i+1*plift->nface);
        *(plift->aeroforce+2) += *(plift->Deltap+i)* *(plift->nsurf+i)* *(plift->normal+i+2*plift->nface);
        *(plift->aeroforce+3) += *(plift->Deltad+i);
    }
}

int main(void)
{
    FILE *ofp, *ifp;
    char *Infile, Wngfile[20], HTailfile[20], VTailfile[20], *Outfile, line[110], code[7];
    struct liftsurf flap,aileron,wing,htail,elevator,vtail,rudder;
    struct liftsurf *pflap = &flap;
    struct liftsurf *paileron = &aileron;
    struct liftsurf *pwing = &wing;
    struct liftsurf *phtail = &htail;
    struct liftsurf *pelevator = &elevator;
    struct liftsurf *pvtail = &vtail;
    struct liftsurf *prudder = &rudder;
    double *invAN,*AN,*BN,*RHS, *Gammas, *wind;
    double UVW[3],Q,aoa,aoadegrees,yaw,yawdegrees,pi,dt,timestep_denom,MAC,rho,*totalforce;
    double delta[2],beta[2],eta[2],zeta[2],deltadegrees,betadegrees,etadegrees,zetadegrees;
    int i,j,mtn,it,ntimes,m,n,freewake,mht,nht,mvt,nvt;
    
    Infile="infile.arp";
    Outfile="outfile.m";

    printf("Reading parameters from %s\n",Infile);
    ifp = fopen(Infile,"r");
    if(ifp == NULL) {
        fprintf(stderr,"Error:  Could not open %s\n",Infile);
        exit(1);
    }

    while(fgets(line, 110, ifp) != NULL)
    {
        if ( strncmp("INP1",line,4) == 0 ){
            sscanf(line,"%s\t%lf\t%lf\t%lf\t%lf",code,&Q,&rho,&aoadegrees,&yawdegrees);
            printf ("%f %f %f %f\n",Q,rho,aoadegrees,yawdegrees);
        }
        if ( strncmp("INP2",line,4) == 0 ){
            sscanf(line,"%s\t%i\t%i\t%i",code,&m,&mht,&mvt);
            printf ("%i %i %i\n",m,mht,mvt);
        }
        if ( strncmp("INP3",line,4) == 0 ){
            sscanf(line,"%s\t%i\t%i\t%i",code,&n,&nht,&nvt);
            printf ("%i %i %i\n",n,nht,nvt);
        }
        if ( strncmp("INP4",line,4) == 0 ){
            sscanf(line,"%s\t%i\t%lf\t%i",code,&ntimes,&timestep_denom,&freewake);
            printf ("%i %f %i\n",ntimes,timestep_denom,freewake);
        }
        if ( strncmp("INP5",line,4) == 0 ){
            sscanf(line,"%s\t%lf\t%lf\t%lf\t%lf",code,&deltadegrees,&betadegrees,&etadegrees,&zetadegrees);
            printf ("%f %f %f %f\n",deltadegrees,betadegrees,etadegrees,zetadegrees);
        }
        if ( strncmp("INP6",line,4) == 0 ){
            sscanf(line,"%s\t%s",code,Wngfile);
            printf ("%s\n",Wngfile);
        }
        if ( strncmp("INP7",line,4) == 0 ){
            sscanf(line,"%s\t%s",code,HTailfile);
            printf ("%s\n",HTailfile);
        }
        if ( strncmp("INP8",line,4) == 0 ){
            sscanf(line,"%s\t%s",code,VTailfile);
            printf ("%s\n",VTailfile);
        }
    }
    
    fclose(ifp);    
        
    totalforce=(double *)malloc(sizeof(double)*ntimes*4);
    pi=atan(1.0)*4.0;
    aoa=aoadegrees*pi/180.0;
    yaw=yawdegrees*pi/180.0;
    delta[0]=(-deltadegrees)*pi/180.0;    /* left aileron angle; negative means up*/
    delta[1]=deltadegrees*pi/180.0;     /* right aileron angle; negative means up */
    beta[0]=betadegrees*pi/180.0;      /* left flap angle; negative means up*/
    beta[1]=beta[0];            /* right flap angle; equal to left */
    eta[0]=etadegrees*pi/180.0;    /* left elevator angle; negative means up*/
    eta[1]=eta[0];              /* right elevator angle; equal to left */
    zeta[0]=zetadegrees*pi/180.0;      /* Real rudder */
    zeta[1]=zeta[0];            /* Rudder mirror image */
    /* Calculate free stream flow components */
    UVW[0]=Q*cos(aoa)*cos(yaw);
    UVW[1]=-Q*sin(yaw);
    UVW[2]=Q*sin(aoa)*cos(yaw);
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
    AN=(double *)malloc(sizeof(double)*mtn*mtn); /* Normal flow coefficient matrix */
    invAN=(double *)malloc(sizeof(double)*mtn*mtn); /* Inverse of normal flow coefficient matrix */
    BN=(double *)malloc(sizeof(double)*mtn*mtn); /* Downwash coefficient matrix */
    RHS=(double *)malloc(sizeof(double)*mtn); /* Right Hand Side vector */
    Gammas = (double *)malloc(sizeof(double)*mtn); /* Vorticity vector history */
    wind=(double *)malloc(sizeof(double)*mtn); /* Induced windspeeds */

    /* Calculate influence coefficient matrices */
    cycleliftsurf(pwing,pflap,paileron,phtail,pelevator,pvtail,prudder,AN,BN,mtn);    

    /* Invert the influence coefficient matrix */
    gauss(AN,invAN,mtn);

    /* Assign memory for gamma and wind values in every liftsurf */
    pwing->gamma=(double *)malloc(sizeof(double)*pwing->nface*ntimes);
    pwing->wind=(double *)malloc(sizeof(double)*pwing->nface);
    pflap->gamma=(double *)malloc(sizeof(double)*pflap->nface*ntimes);
    pflap->wind=(double *)malloc(sizeof(double)*pflap->nface);
    paileron->gamma=(double *)malloc(sizeof(double)*paileron->nface*ntimes);
    paileron->wind=(double *)malloc(sizeof(double)*paileron->nface);
    phtail->gamma=(double *)malloc(sizeof(double)*phtail->nface*ntimes);
    phtail->wind=(double *)malloc(sizeof(double)*phtail->nface);
    pelevator->gamma=(double *)malloc(sizeof(double)*pelevator->nface*ntimes);
    pelevator->wind=(double *)malloc(sizeof(double)*pelevator->nface);
    pvtail->gamma=(double *)malloc(sizeof(double)*pvtail->nface*ntimes);
    pvtail->wind=(double *)malloc(sizeof(double)*pvtail->nface);
    prudder->gamma=(double *)malloc(sizeof(double)*prudder->nface*ntimes);
    prudder->wind=(double *)malloc(sizeof(double)*prudder->nface);
    
    /* Assign memory for uvw values in every liftsurf */
    pwing->uvw=(double *)malloc(sizeof(double)*pwing->nface*3);
    pflap->uvw=(double *)malloc(sizeof(double)*pflap->nface*3);
    paileron->uvw=(double *)malloc(sizeof(double)*paileron->nface*3);
    phtail->uvw=(double *)malloc(sizeof(double)*phtail->nface*3);
    pelevator->uvw=(double *)malloc(sizeof(double)*pelevator->nface*3);
    pvtail->uvw=(double *)malloc(sizeof(double)*pvtail->nface*3);
    prudder->uvw=(double *)malloc(sizeof(double)*prudder->nface*3);
    
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
    ofp = fopen(Outfile, "w");
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
    
    /* Close outpout file */
     fclose(ofp);    
}