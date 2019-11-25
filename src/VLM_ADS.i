// Copyright 2019 Université de Liège
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

//
// VLM_ADS.i
//
// SWIG interface file for VLM code
// 

%module CVLM

%typemap(in) (int argc, char *argv[]) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $1 = PyList_Size($input);
  $2 = (char **) malloc(($1+1)*sizeof(char *));
  for (i = 0; i < $1; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (!PyString_Check(s)) {
        free($2);
        PyErr_SetString(PyExc_ValueError, "List items must be strings");
        return NULL;
    }
    $2[i] = PyString_AsString(s);
  }
  $2[i] = 0;
}

%typemap(freearg) (int argc, char *argv[]) {
   if ($2) free($2);
}

%{
    #define SWIG_FILE_WITH_INIT
    #include "VLM_ADS.h"
    #include "vAirfoil.h"
    #include "vLiftsurf.h"
    #include "vVLMData.h"
    #include "vSetup.h"
    #include "vGeometrySetup.h"
    #include "vInfluence.h"
    #include "vIteration.h"
    #include "vOutput.h"
    #include "vControl.h"
    #include "vVortex.h"
%}

%include "vLiftsurf.h"
%include "vVLMData.h"

%include "numpy.i"
%include "typemaps.i"
%init %{
import_array();
%}
%include "carrays.i"
%array_functions(double, Deltap)
%array_functions(double, Deltad)
%array_functions(double, vertices)
%array_functions(double, vortex)
%array_functions(double, aeroforce)
%array_functions(int, neighbours)
%array_functions(double, normal)
%array_functions(int, faces)
%array_functions(double, UVW)
%array_functions(double, nsurf)
%apply (double IN_ARRAY1[ANY]) {(double delta[2])}

%typemap(in, numinputs = 1) char* Infile {
  $1 = PyString_AsString($input);
}
%typemap(in) (struct VLMData *data) (int temp) {
  temp = SWIG_ConvertPtr($input, (void **) &$1, $1_descriptor , 0);
}

%typemap(in) (struct liftsurf *wing) (int temp) {
  temp = SWIG_ConvertPtr($input, (void **) &$1, $1_descriptor , 0);
}

%include "vSetup.h"
%include "vControl.h"
%include "vGeometrySetup.h"
%include "vInfluence.h"
%include "vIteration.h"
%include "vOutput.h"
%include "vControl.h"
%include "vVortex.h"

%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* xpline)}
%apply (int DIM1, double* ARGOUT_ARRAY1) {(int len2, double* ycamber)}
%rename (nacafourfivedigit) naca_four_five;
%exception naca_four_five {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%inline %{
void naca_four_five(int len1, double* xpline, int len2, double* ycamber, char airfoil[]) {
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
    }
    return nacafourfivedigit(xpline, ycamber, len1, "NACA 0001");
}
%}

%include "VLM_ADS.h"
%include "vAirfoil.h"
