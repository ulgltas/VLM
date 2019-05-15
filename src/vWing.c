//
//  vWing.c
//
//  Wing setup
//
//

#include "vWing.h"
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

double wingsetup(struct liftsurf *pflap, struct liftsurf *paileron, struct liftsurf *pwing, char *Wngfile, int m, int n, int *ChkAil, int *ChkWTED)
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
    int iTS, lstWngTSNumber, cond, ChkWLED, ChkWGL, WTEDStopPointsNumber, npTS, npTSp1, ndouble;
    int Ailinds[2][2],WTEDinds[2][2],Wglinds[2][2];
    char line[110], code[8], WngArf[12], nindex[12];
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
            sscanf(line,"%s %i",code,ChkAil);
            cond=1;
        }
    }
    cond=0;
    while(cond == 0){
        fgets(line, 110, fp1);
        if ( *ChkAil == 1 ){
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
    if (*ChkAil == 0)
    {
        AilSpan = 0.0;
        AilPosSpan = 0.0;
        AilRlChord = 0.0;
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
            sscanf(line,"%s %i",code, ChkWTED);
            cond=1;
        }
    }
    cond=0;
    while( cond == 0 ){
        fgets(line, 110, fp1);
        if ( *ChkWTED == 1 ){
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
    if (*ChkWTED == 0)
    {
        WTEDSpan = 0.0;
        WTEDPosSpan = 0.0;
        WTEDRlChord = 0.0;
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
    nypos=lstWngTSNumber+1;
    nxpos=4;
    nxyzTS=lstWngTSNumber+1;
    if (*ChkAil == 1)
    {
        nypos++;
    }
    if (*ChkWTED == 1)
    {
        nypos++;
    }

    if (  ChkWGL == 1){
        WglLeOffset=0.01;
        WglSwpLE=WglSwpLE*pi/180.0;
        WglDhdrl=WglDhdrl*pi/180.0;
        nypos+=3; // Three panels in the winglet
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
    if (*ChkWTED == 1)
    {
        *(ypos + i) = WTEDPosSpan; i++;
        *(ypos + i) = WTEDSpan + WTEDPosSpan; i++;
    }
    if (*ChkAil == 1)
    {
        *(ypos + i) = AilPosSpan; i++;
        *(ypos + i) = AilSpan + AilPosSpan; i++;
    }
    
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

    if (*ChkWTED == 1)
    {
        WTEDinds[0][1] = findindex(ypline, np1, WTEDPosSpan);
        WTEDinds[1][1] = findindex(ypline, np1, WTEDSpan + WTEDPosSpan);
    }
    
    if (*ChkAil == 1)
    {
        Ailinds[0][1] = findindex(ypline, np1, AilPosSpan);
        Ailinds[1][1] = findindex(ypline, np1, AilSpan + AilPosSpan);
    }
    
    if (ChkWGL == 1)
    {
        Wglinds[0][1] = findindex(ypline, np1, *(yTS + lstWngTSNumber));
        Wglinds[1][1] = findindex(ypline, np1, *(yTS + lstWngTSNumber) + WglSpan);
    }
    
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
    if (ChkWGL==1)
    {
        Wglinds[0][0] = findindex(xpline, mp1, WglLeOffset / WngTSTpChord[lstWngTSNumber - 1]);
        Wglinds[1][0] = findindex(xpline, mp1, (WglLeOffset + WglRtChord) / WngTSTpChord[lstWngTSNumber - 1]);
    }
    

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
            if (*ChkWTED == 1)
            {
                if (i >= WTEDinds[0][0] && i <= WTEDinds[1][0] && j >= WTEDinds[0][1] && j <= WTEDinds[1][1])
                { /* Flap */
                    flaphere = 1;
                    if (i == WTEDinds[0][0])
                    { /* Flap leading edge */
                        winghere = 1;
                    }
                    if (j == WTEDinds[0][1])
                    { /* Flap inboard edge */
                        winghere = 1;
                    }
                    if (j == WTEDinds[1][1])
                    { /* Flap outboard edge */
                        winghere = 1;
                    }
                }
            }

            /* Check if this point lies on the aileron */
            if (*ChkAil == 1)
            {
                if (i >= Ailinds[0][0] && i <= Ailinds[1][0] && j >= Ailinds[0][1] && j <= Ailinds[1][1])
                { /* Aileron */
                    aileronhere = 1;
                    if (i == Ailinds[0][0])
                    { /* Aileron leading edge */
                        winghere = 1;
                    }
                    if (j == Ailinds[0][1])
                    { /* Aileron inboard edge */
                        winghere = 1;
                    }
                    if (j == Ailinds[1][1])
                    { /* Aileron outboard edge */
                        winghere = 1;
                    }
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