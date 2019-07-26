//
//  vAssignValues.h
//
//  Assign values back to the Liftsurf struct
//
//
#ifndef vAssignValues_h
#define vAssignValues_h

#include "vLiftsurf.h"
#include "vVLMData.h"

void assignGammas(struct VLMData *data, int it);

void assignwind(struct VLMData *data);

#endif /* vAssignValues_h */