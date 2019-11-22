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
//  vHTail.c
//
//  Horizontal tail setup
//
//

#include "vHTail.h"
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

void htailsetup(struct liftsurf *phtail, struct liftsurf *pelevator, char *HTailfile, int m, int n, int *ChkElev)
{
    /* Set up horizontal tail */
    
    double pi;
    int i,j,nxpos,nypos,nxyzTS;
    int mp1,np1;
    int htailhere,elevatorhere;
    int htail_nvert,elevator_nvert,htail_nface,elevator_nface;
    
    double *ypos,*xpos;
    double *ypline,*xpline,*xpgrid,*ypgrid,*zpgrid,yhere,*xTS,*yTS,*yTS2,*zTS,twistcentre,*ycamber,*ycamberall;
    double *xv,*yv,*zv,dxw,minchord;
    double *xhtail,*yhtail,*zhtail,*xelevator,*yelevator,*zelevator,*chordvec,*levec;
    double *xvhtail,*yvhtail,*zvhtail,*xvelevator,*yvelevator,*zvelevator;
    double zplineRt,zplineTp,twistangle,dihedral,xpTS,ypTS,zpTS;
    int *ijhtail, *ijelevator;
       
    FILE *fp1;
    int iTS, HTTSNumber, cond, npTS, ndouble;
    int ELVinds[2][2],dummyint;
    char line[110], code[8], HTType[12], HTArf[12], HTSurfFinish[12];
    double HTTSLength[3], HTTSRtChord[3], HTTSTpChord[3];
    double HTTSRtInc[3], HTTSTpInc[3], HTTSSwpLE[3], HTTSDhdr[3], HTTSTwist[3], HTTSTR[3];
    double HTSpan, HTLPosFus, HTVPosFus, HTRtChord, HTTpChord, HTArea, HTSwpLE, HTTwist, HTRlInc, HTInc, HTDhdrl;
    double HTAR, HTTR, HTVolCoeff;
    double ElevSpan, ElevPosSpan, ElevArea, ElevHingeLoc,ElevRlChord,ElevRlSpan;
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
            sscanf(line,"%s %i",code, ChkElev);
            cond=1;
        }
    }
    cond=0;
    while(cond == 0){
        fgets(line, 110, fp1);
        if ( *ChkElev == 1 ){
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
    
    nypos=HTTSNumber+1;
    nxpos=2;
    if (*ChkElev == 1)
    {
        nypos += 2;
        nxpos++;
    }
    
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
    if (*ChkElev == 1)
    {
        *(ypos+i)=ElevPosSpan;i++;
        *(ypos+i)=ElevSpan/2.0+ElevPosSpan; /* ElevSpan is the total elevator span */
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
        fprintf(stderr,"Error:  The number of spanwise panels on the horizontal tail should be larger than %i\n",nypos);
        exit(1);
    }    
    /* Create vector containing the full spanwise grid */
    createypline(ypos,nypos,ypline,np1);
    free(ypos);
    /* Check between which elements of ypline lies the elevator */    
    ELVinds[0][1]=findindex(ypline,np1,ElevPosSpan);
    ELVinds[1][1]=findindex(ypline,np1,ElevSpan/2+ElevPosSpan);

    /* Create vector containing x-coords of TS and Elev */
    xpos= (double *)malloc(sizeof(double)*nxpos); 
    *(xpos+0)=0.0;
    i = 1;
    if (*ChkElev == 1)
    {
        *(xpos+1)=(100.0-ElevRlChord)/100.0;
        i++;
    }
    *(xpos+i)=1.0;
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
    free(xpos);
    /* Check between which elements of xpline lies the elevator */
    if (*ChkElev == 1)
    {
        ELVinds[0][0]=findindex(xpline,mp1,(100.0-ElevRlChord)/100.0);
        ELVinds[1][0]=m; /* Will always lie on the trailing edge */
    }

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
    free(xpline);
    free(chordvec);
    free(levec);
    free(ycamber);
    free(ycamberall);
    free(xTS);
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
            if (*ChkElev == 1 && i >= ELVinds[0][0] && i <= ELVinds[1][0] && j >= ELVinds[0][1] && j <= ELVinds[1][1]){ /* Elevator */
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
    free(ypline);
    free(xv);
    free(yv);
    free(zv);
    free(xpgrid);
    free(ypgrid);
    free(zpgrid);
    free(yTS);
    free(yTS2);
    free(zTS);
    /* Count number of panels on all halves of lifting surfaces */
    elevator_nface=countfaces(ijelevator,elevator_nvert,mp1,np1);
    htail_nface=countfaces(ijhtail,htail_nvert,mp1,np1);   
    
    /* Create all the panel information in the elevator liftsurf structure */
    pelevator->nface=2*elevator_nface;
    pelevator->nvert=2*elevator_nvert;
    pelevator->faces=(int *)malloc(sizeof(int)*pelevator->nface*4); 
    pelevator->shedding=(int *)malloc(sizeof(int)*pelevator->nface); 
    arrangefaces(ijelevator,mp1,np1,pelevator);  
    free(ijelevator);
    /* Create all the panel information in the htail liftsurf structure */
    phtail->nface=2*htail_nface;
    phtail->nvert=2*htail_nvert;
    phtail->faces=(int *)malloc(sizeof(int)*phtail->nface*4); 
    phtail->shedding=(int *)malloc(sizeof(int)*phtail->nface); 
    arrangefaces(ijhtail,mp1,np1,phtail);
    free(ijhtail);
    /* Copy all panel vertex information to the relevant liftsurf structures */
    assignvertices(pelevator,xelevator,yelevator,zelevator,xvelevator,yvelevator,zvelevator);
    free(xelevator);
    free(yelevator);
    free(zelevator);
    free(xvelevator);
    free(yvelevator);
    free(zvelevator);
    assignvertices(phtail,xhtail,yhtail,zhtail,xvhtail,yvhtail,zvhtail);
    free(xhtail);
    free(yhtail);
    free(zhtail);
    free(xvhtail);
    free(yvhtail);
    free(zvhtail);
}