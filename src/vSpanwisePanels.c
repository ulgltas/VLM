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
//  vSpanwisePanels.c
//  
//
//
//
//

#include "vSpanwisePanels.h"
#include <stdlib.h>
#include <math.h>

void createypline(double *ypos,int nypos,double *ypline,int np1)
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
