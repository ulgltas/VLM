# Copyright 2019 University of Liège
# Main CMake file.
# Configure the compilation of the C core of the VLM code.
#
# Authors : M. SANCHEZ MARTINEZ.

ADD_LIBRARY(VLMCore STATIC vAirfoil.h vAirfoil.c)

ADD_EXECUTABLE(VLM VLM_ADS.c)

TARGET_LINK_LIBRARIES(VLM VLMCore)