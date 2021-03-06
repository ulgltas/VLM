# Copyright 2019 Université de Liège
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# Main CMake file.
# Configure the compilation of the C core of the VLM code.
#
# Some parts based on the CUPyDO CMake file by D. THOMAS
# 
# Authors : M. SANCHEZ MARTINEZ.

CMAKE_MINIMUM_REQUIRED (VERSION 3.0)
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


# Python wrapper with SWIG
option(PYTHON_WRAPPER "Build Python wrapper" ON)

if(PYTHON_WRAPPER)
  FIND_PACKAGE(PythonInterp 3 REQUIRED) # Interpreter
  FIND_PACKAGE(PythonLibs 3 REQUIRED)
  FIND_PACKAGE(NUMPY REQUIRED)
  FIND_PACKAGE(SWIG REQUIRED)
  MESSAGE(STATUS "PythonLibs=${PYTHON_LIBRARIES}")
  MESSAGE(STATUS "Python=${PYTHON_EXECUTABLE}")
  MESSAGE(STATUS "SWIG_EXECUTABLE=${SWIG_EXECUTABLE}")
  INCLUDE(${SWIG_USE_FILE})
  INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src
                      ${PYTHON_INCLUDE_PATH}
                      ${NUMPY_INCLUDE_DIR})
endif(PYTHON_WRAPPER)

# Source subdirectory
ADD_SUBDIRECTORY( src )


# --- Install ---

if(CMAKE_BUILD_TYPE MATCHES DEBUG )
  if(MSVC)
    ADD_COMPILE_OPTIONS(/W4)
  else(MSVC)
    ADD_COMPILE_OPTIONS(-Wall -Wextra -pedantic)
  endif(MSVC)
endif(CMAKE_BUILD_TYPE MATCHES DEBUG )

SET(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")
MESSAGE(STATUS "Installation directory : ${CMAKE_INSTALL_PREFIX}")


# --- Enable CTest ---
ENABLE_TESTING()
ADD_SUBDIRECTORY(tests)
