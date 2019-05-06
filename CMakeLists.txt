# Copyright 2019 University of Liège
# Main CMake file.
# Configure the compilation of the C core of the VLM code.
#
# Some parts based on the CUPyDO CMake file by D. THOMAS
# 
# Authors : M. SANCHEZ MARTINEZ.

CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (VLM)

# Type of build
IF( NOT CMAKE_BUILD_TYPE )
  SET( CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: None Debug Release" FORCE )
ENDIF( NOT CMAKE_BUILD_TYPE )
MESSAGE(STATUS "Build type : ${CMAKE_BUILD_TYPE}")

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")
# macros/fcts
INCLUDE(VLMMacros)

# Source subdirectory
ADD_SUBDIRECTORY( src )


# --- Install ---
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    TARGET_LINK_LIBRARIES(VLM m) # Adds -lm option for linking
endif(CMAKE_C_COMPILER_ID STREQUAL "GNU")

if(CMAKE_BUILD_TYPE MATCHES DEBUG )
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
endif(CMAKE_BUILD_TYPE MATCHES DEBUG )

SET(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")
MESSAGE(STATUS "Installation directory : ${CMAKE_INSTALL_PREFIX}")


# --- Enable CTest ---
ENABLE_TESTING()
ADD_SUBDIRECTORY(tests)