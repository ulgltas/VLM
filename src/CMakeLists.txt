# Copyright 2019 University of Liège
# CMake file for the VLM source
# Configure the compilation of the C core of the VLM code.
#
# Authors : M. SANCHEZ MARTINEZ.

ADD_LIBRARY(VLMCore STATIC vAirfoil.h vAirfoil.c vSpanwisePanels.h vSpanwisePanels.c vUtilities.h vUtilities.c vLiftsurf.h vPanel.h vPanel.c vVertices.h vVertices.c vVTail.h vVTail.c vHTail.h vHTail.c vWing.h vWing.c vGeometry.h vGeometry.c vVortex.h vVortex.c vWake.h vWake.c)

ADD_EXECUTABLE(VLM VLM_ADS.c)

TARGET_LINK_LIBRARIES(VLM VLMCore)

install (TARGETS VLM DESTINATION "${PROJECT_SOURCE_DIR}/bin")