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

# CMake file for the VLM source
# Configure the compilation of the C core of the VLM code.
#
# Authors : M. SANCHEZ MARTINEZ.

FILE(GLOB SRCS *.h *.c) # Stores these files in variable SRCS
SET(SWIGSRCS VLM_ADS.i)


ADD_LIBRARY(VLMCore STATIC ${SRCS})

ADD_EXECUTABLE(VLM VLM_ADS.c)

TARGET_LINK_LIBRARIES(VLM VLMCore)

set(TARGS VLM VLMCore)
set(CMAKE_SWIG_FLAGS )

IF(${CMAKE_VERSION} VERSION_LESS "3.8.0" AND PYTHON_WRAPPER)
  SET(CMAKE_SWIG_OUTDIR "${PROJECT_SOURCE_DIR}/bin")
  SWIG_ADD_MODULE(CVLM python ${SWIGSRCS} ${SRCS})
  SWIG_LINK_LIBRARIES(CVLM ${PYTHON_LIBRARIES})
  set(TARGS ${TARGS} _CVLM)
ELSEIF(PYTHON_WRAPPER)
  SET(CMAKE_SWIG_OUTDIR "${PROJECT_SOURCE_DIR}/bin")
  SWIG_ADD_LIBRARY(CVLM
                   LANGUAGE python
                   OUTPUT_DIR "${PROJECT_SOURCE_DIR}/bin"
                   SOURCES ${SWIGSRCS} ${SRCS})
  SWIG_LINK_LIBRARIES(CVLM ${PYTHON_LIBRARIES})
  set(TARGS ${TARGS} _CVLM)
ENDIF()

install (TARGETS ${TARGS} DESTINATION "${PROJECT_SOURCE_DIR}/bin")

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    TARGET_LINK_LIBRARIES(VLM m) # Adds -lm option for linking
endif(CMAKE_C_COMPILER_ID STREQUAL "GNU")

