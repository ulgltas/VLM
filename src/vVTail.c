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
//  vVTail.c
//
//  Vertical tail setup
//
//

#include "vVTail.h"
#include "vAirfoil.h"
#include "vLiftsurf.h"
#include "vPanel.h"
#include "vSpanwisePanels.h"
#include "vUtilities.h"
#include "vVertices.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void vtailsetup(struct liftsurf *pvtail, struct liftsurf *prudder, char *VTailfile, int m, int n, int *ChkRdr)
{
    /* Set up vertical tail */
    
    double pi;
    int i,j,nxpos,nypos,nxyzTS;
    int mp1,np1;
    int vtailhere,rudderhere;
    int vtail_nvert,rudder_nvert,vtail_nface,rudder_nface;
    
    double *ypos,*xpos;
    double *ypline,*xpline,*xpgrid,*ypgrid,*zpgrid,yhere,*xTS,*yTS,*yTS2,*zTS,*ycamber,*ycamberall;
    double *xv,*yv,*zv,dxw,minchord;
    double *xvtail,*yvtail,*zvtail,*xrudder,*yrudder,*zrudder,*chordvec,*levec;
    double *xvvtail,*yvvtail,*zvvtail,*xvrudder,*yvrudder,*zvrudder;
    double zplineRt,zplineTp,twistangle,dihedral,xpTS,ypTS,zpTS;
    int *ijvtail, *ijrudder;

    FILE *fp1;
    int iTS, VTTSNumber, cond, ndouble;
    int RDRinds[2][2],dummyint,OptVTFusMounted,OptVTTailMounted,OptVTWingMounted;
    char line[110], code[8], VTType[12], VTArf[12], VTSurfFinish[12];
    double VTTSLength[2], VTTSRtChord[2], VTTSTpChord[2];
    double VTTSSwpLE[2], VTTSDhdr[2], VTTSTR[2];
    double VTSpan, VTLPosFus, VTVPosFus, VTRtChord, VTTpChord, VTSPosFus, VTArea, VTSwpLE, VTVAngle;
    double VTAR, VTTR, VTVolCoeff,VTRlTpChord,VTRlPosFus,VTRlVPosHT;
    double RdrSpan, RdrPosSpan, RdrArea, RdrHingeLoc,RdrRlChord,RdrRlSpan;
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
            sscanf(line,"%s %i",code, ChkRdr);
            cond=1;
        }
    }
    cond=0;
    while(cond == 0){
        fgets(line, 110, fp1);
        if ( *ChkRdr == 1 ){
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
    
    nypos=VTTSNumber+1;
    nxpos=2;
    if (*ChkRdr == 1)
    {
        nypos += 2;
        nxpos++;
    }

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
    if (*ChkRdr == 1)
    {
        *(ypos+i)=RdrPosSpan;
        i++;
        *(ypos+i)=RdrSpan+RdrPosSpan;
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
        fprintf(stderr,"Error:  The number of spanwise panels on the vertical tail should be larger than %i\n",nypos);
        exit(1);
    }    
    /* Create vector containing the full spanwise grid */
    createypline(ypos,nypos,ypline,np1);
    /* Check between which elements of ypline lies the elevator */    
    RDRinds[0][1]=findindex(ypline,np1,RdrPosSpan);
    RDRinds[1][1]=findindex(ypline,np1,RdrSpan+RdrPosSpan); 
    
    /* Create vector containing x-coords of TS and Rdr */
    xpos = (double *)malloc(sizeof(double)*nxpos);
    *(xpos+0)=0.0;
    i = 1;
    if (*ChkRdr == 1)
    {
        *(xpos+i)=(100.0-RdrRlChord)/100.0;
        i++; 
    }
    *(xpos+i)=1.0;
    
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
    free(ypos);
    free(xpos);
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
    free(xpline);
    free(ycamberall);
    free(chordvec);
    free(levec);
    free(ycamber);
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
            if (*ChkRdr==1 && i >= RDRinds[0][0] && i <= RDRinds[1][0] && j >= RDRinds[0][1] && j <= RDRinds[1][1]){ /* Elevator */
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
    free(ypline);
    free(xv);
    free(yv);
    free(zv);
    free(xpgrid);
    free(ypgrid);
    free(zpgrid);
    free(xTS);
    free(yTS);
    free(yTS2);
    free(zTS);
    /* Count number of panels on all halves of lifting surfaces */
    rudder_nface=countfaces(ijrudder,rudder_nvert,mp1,np1);
    vtail_nface=countfaces(ijvtail,vtail_nvert,mp1,np1);   
    
    /* Create all the panel information in the rudder liftsurf structure */
    prudder->nface=2*rudder_nface;
    prudder->nvert=2*rudder_nvert;
    prudder->faces=(int *)malloc(sizeof(int)*prudder->nface*4); 
    prudder->shedding=(int *)malloc(sizeof(int)*prudder->nface); 
    arrangefaces(ijrudder,mp1,np1,prudder);  
    free(ijrudder);
    /* Create all the panel information in the vtail liftsurf structure */
    pvtail->nface=2*vtail_nface;
    pvtail->nvert=2*vtail_nvert;
    pvtail->faces=(int *)malloc(sizeof(int)*pvtail->nface*4); 
    pvtail->shedding=(int *)malloc(sizeof(int)*pvtail->nface); 
    arrangefaces(ijvtail,mp1,np1,pvtail);
    free(ijvtail);
    /* Copy all panel vertex information to the relevant liftsurf structures */
    assignvertices_fin(prudder,xrudder,yrudder,zrudder,xvrudder,yvrudder,zvrudder,OptVTFusMounted,OptVTTailMounted,OptVTWingMounted);
    free(xrudder);
    free(yrudder);
    free(zrudder);
    free(xvrudder);
    free(yvrudder);
    free(zvrudder);
    assignvertices_fin(pvtail,xvtail,yvtail,zvtail,xvvtail,yvvtail,zvvtail,OptVTFusMounted,OptVTTailMounted,OptVTWingMounted);
    free(xvtail);
    free(yvtail);
    free(zvtail);
    free(xvvtail);
    free(yvvtail);
    free(zvvtail);
}