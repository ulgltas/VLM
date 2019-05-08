# Copyright 2019 University of Liège
# CMake file for the VLM source
# Configure the compilation of the C core of the VLM code.
#
# Authors : M. SANCHEZ MARTINEZ.

FILE(GLOB SRCS *.h v*.c) # Stores these files in variable SRCS

ADD_LIBRARY(VLMCore STATIC ${SRCS})

ADD_EXECUTABLE(VLM VLM_ADS.c)

TARGET_LINK_LIBRARIES(VLM VLMCore)

install (TARGETS VLM DESTINATION "${PROJECT_SOURCE_DIR}/bin")

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    TARGET_LINK_LIBRARIES(VLM m) # Adds -lm option for linking
endif(CMAKE_C_COMPILER_ID STREQUAL "GNU")
