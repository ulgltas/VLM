#include "VLM_ADS.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vGeometrySetup.h"
#include "vSetup.h"
#include "vInput.h"
#include "vOutput.h"
#include "vInfluence.h"
#include "vIteration.h"
#include "vLiftsurf.h"
#include "vVLMData.h"

int main(int argc, char *argv[])
{
    char *Infile, *Outfile;
    
    
    if (argc>1) {
        Infile=argv[1];
    } else {
        Infile="infile.arp";
    }

    if (argc>2) {
        Outfile=argv[2];
    } else {
        Outfile="outfile.m";
    }
    run(Infile, Outfile);
    return 0;
}
void run(char *Infile, char *Outfile)
{
    char Wngfile[60], HTailfile[60], VTailfile[60];
    struct VLMData data;
    int i, j, it;
    
    setup(Infile, &data);
    geometry_setup(&data);

    /* Calculate influence coefficient matrices */
    cycleliftsurf(&data);
    memory_setup(&data);
    for (it=0;it<(data.ntimes)-1;it++)
    {
        iteration(&data, it);
    }
    it--;
    printf("forcex=%f, forcey=%f, forcez=%f, draginduced=%f\n",*(data.totalforce+0+it*4),*(data.totalforce+1+it*4),*(data.totalforce+2+it*4),*(data.totalforce+3+it*4));
    
    /* Print output */
    exportTextOutput(Outfile, it, data);
}