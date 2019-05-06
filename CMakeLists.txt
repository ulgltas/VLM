# Copyright 2019 University of Liège
# Main CMake file.
# Configure the compilation of the C core of the VLM code.
#
# Authors : M. SANCHEZ MARTINEZ.

CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (VLM)

LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")
# macros/fcts
INCLUDE(VLMMacros)

# Source subdirectory
ADD_SUBDIRECTORY( src )


# --- Install ---
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    TARGET_LINK_LIBRARIES(VLM m) # Adds -lm option for linking
endif(CMAKE_C_COMPILER_ID STREQUAL "GNU")

SET(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")
MESSAGE(STATUS "Installation directory : ${CMAKE_INSTALL_PREFIX}")
install (TARGETS VLM DESTINATION bin)

# --- Enable CTest ---
ENABLE_TESTING()
ADD_SUBDIRECTORY(tests)