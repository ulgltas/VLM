//
// VLM_ADS.i
//
// SWIG interface file for VLM code
// 

%module PyVLM

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
%}
%include "numpy.i"
%init %{
import_array();
%}

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
