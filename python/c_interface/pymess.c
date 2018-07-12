//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
//               2009-2018
//

/**
 * @addtogroup interfaces_python
 * @{
 * @file python/c_interface/pymess.c
 * @brief Functions which should be called from @python.
 *
 * @author @nitin
 * @author @mbehr
 *
 */

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL MESS_VECTOR_MATRIX_PYTHON_C_API
#include <numpy/arrayobject.h>
#include "mess/mess.h"
#include "pymess.h"


#if PY_MAJOR_VERSION < 3
#define ConvStringtoC(X) PyBytes_AsString(X)
#else
#define ConvStringtoC(X) PyBytes_AsString(PyUnicode_AsUTF8String(X))
#endif

#if PY_MAJOR_VERSION >= 3
static int init_numpy(void) {
    import_array1(-1);
    return 0;
}
#else
static void init_numpy(void){
    import_array();
}
#endif



///////////////////////////////////////////////////////////////////////////////////////////////////
// MODULE INITIALISATIONS
// This part of the code is exceptionally messy to make sure that
// it works with both Python 2 and Python 3
///////////////////////////////////////////////////////////////////////////////////////////////////

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
/*
   static PyObject * error_out(PyObject *m) {
   struct module_state *st = GETSTATE(m);
   PyErr_SetString(st->error, "Something went wrong!");
   return NULL;
   }
   */
#endif


    PyObject* pymess_mess_is_debug(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_zlib(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_bzip2(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_umfpack(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_amd(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_colamd(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_cholmod(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_csparse(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_superlu(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_mklpardiso(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_arpack(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_matio(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_openmp(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_mess64(PyObject *self, PyObject* args);
    PyObject* pymess_mess_git_id(PyObject *self, PyObject* args);
    PyObject* pymess_mess_git_branch(PyObject *self, PyObject* args);



/**
 * @brief Function which are callable from @python.
 * Edit below function if you need to add more functions that you wish to be called from @python
 * Avoid editing any other piece of code in the Initialisation section
 *
 * */
static PyMethodDef interface_python_methods[] = {
    {"test_matrix",             (PyCFunction) pymess_test_matrix,           METH_VARARGS, "Test function to check matrix conversions"},
    {"test_vector",             (PyCFunction) pymess_test_vector,           METH_VARARGS, "Test function to check array conversions"},
    {"test_equation",           (PyCFunction) pymess_test_equation,         METH_VARARGS, "Test function to check eqn transfer"},
    {"_eps",                    (PyCFunction) pymess_eps,                   METH_VARARGS, "Returns mess_eps()"},
    {"_mess_version",           (PyCFunction) pymess_mess_version,          METH_VARARGS, "Calls the mess_version function"},
    {"_mess_version_verbose",   (PyCFunction) pymess_mess_version_verbose,  METH_VARARGS, "Calls the mess_version_verbose function"},
    {"_mess_version_major",     (PyCFunction) pymess_mess_version_major,    METH_VARARGS, "Calls the mess_version_major function"},
    {"_mess_version_minor",     (PyCFunction) pymess_mess_version_minor,    METH_VARARGS, "Calls the mess_version_minor function"},
    {"_mess_version_patch",     (PyCFunction) pymess_mess_version_patch,    METH_VARARGS, "Calls the mess_version_patch function"},
    {"_mess_is_debug",          (PyCFunction) pymess_mess_is_debug,         METH_VARARGS, "Calls the mess_is_debug function"},
    {"_mess_have_zlib",         (PyCFunction) pymess_mess_have_zlib,        METH_VARARGS, "Calls the mess_have_zlib function"},
    {"_mess_have_bzip2",        (PyCFunction) pymess_mess_have_bzip2,       METH_VARARGS, "Calls the mess_have_bzip2 function"},
    {"_mess_have_umfpack",      (PyCFunction) pymess_mess_have_umfpack,     METH_VARARGS, "Calls the mess_have_umfpack function"},
    {"_mess_have_amd",          (PyCFunction) pymess_mess_have_amd,         METH_VARARGS, "Calls the mess_have_amd function"},
    {"_mess_have_colamd",       (PyCFunction) pymess_mess_have_colamd,      METH_VARARGS, "Calls the mess_have_colamd function"},
    {"_mess_have_cholmod",      (PyCFunction) pymess_mess_have_cholmod,     METH_VARARGS, "Calls the mess_have_cholmod function"},
    {"_mess_have_csparse",      (PyCFunction) pymess_mess_have_csparse,     METH_VARARGS, "Calls the mess_have_csparse function"},
    {"_mess_have_superlu",      (PyCFunction) pymess_mess_have_superlu,     METH_VARARGS, "Calls the mess_have_superlu function"},
    {"_mess_have_mklpardiso",   (PyCFunction) pymess_mess_have_mklpardiso,  METH_VARARGS, "Calls the mess_have_mklpardiso function"},
    {"_mess_have_arpack",       (PyCFunction) pymess_mess_have_arpack,      METH_VARARGS, "Calls the mess_have_arpack function"},
    {"_mess_have_matio",        (PyCFunction) pymess_mess_have_matio,       METH_VARARGS, "Calls the mess_have_matio function"},
    {"_mess_have_openmp",       (PyCFunction) pymess_mess_have_openmp,      METH_VARARGS, "Calls the mess_have_openmp function"},
    {"_mess_git_id",            (PyCFunction) pymess_mess_git_id,           METH_VARARGS, "Calls the mess_git_id function"},
    {"_mess_git_branch",        (PyCFunction) pymess_mess_git_branch,       METH_VARARGS, "Calls the mess_git_branch function"},
    {"_mess_have_mess64",       (PyCFunction) pymess_mess_have_mess64,      METH_VARARGS, "Calls the mess_have_mess64 function"},
    {"_direct_select",          (PyCFunction) pymess_direct_select,         METH_VARARGS, "Interface to mess_direct_select"},
    {"_direct_chol_select",     (PyCFunction) pymess_direct_chol_select,    METH_VARARGS, "Interface to mess_direct_chol_select"},
    {"_multidirect_select",     (PyCFunction) pymess_multidirect_select,    METH_VARARGS, "Interface to mess_multidirect_select"},
    {"_set_errorlevel",         (PyCFunction) pymess_set_errorlevel,        METH_VARARGS, "Interface to mess_set_errorlevel."},
    {"_dense_nm_gmpare",        (PyCFunction) pymess_dense_nm_gmpare,       METH_VARARGS | METH_KEYWORDS, "Interface to the C-MESS function mess_dense_nm_gmpare"},
    {"_lradi",                  (PyCFunction) pymess_lradi,                 METH_VARARGS | METH_KEYWORDS, "Interface to the C-MESS lrcfadi"},
    {"_lrnm",                   (PyCFunction) pymess_lrnm,                  METH_VARARGS | METH_KEYWORDS, "Interface to the C-MESS lrnm"},
    {"_lyap",                   (PyCFunction) pymess_lyap,                  METH_VARARGS | METH_KEYWORDS, "Interface to the high-level function mess_lyap"},
    {"_sylvester_sparsedense",  (PyCFunction) pymess_sylvester_sparsedense, METH_VARARGS | METH_KEYWORDS, "Returns solution from mess_direct_create_sylvester_sparsedense solver"},
    {"_care",                   (PyCFunction) pymess_care,                  METH_VARARGS | METH_KEYWORDS, "Interface to the high-level function mess_care"},
    {"_glyap",                  (PyCFunction) pymess_glyap,                 METH_VARARGS | METH_KEYWORDS, "Interface to mess_glyap"},
    {"_gstein",                 (PyCFunction) pymess_gstein,                METH_VARARGS | METH_KEYWORDS, "Interface to mess_gstein"},
    {NULL, NULL, 0, NULL}
};




/*-----------------------------------------------------------------------------
 *  Python >= 3
 *-----------------------------------------------------------------------------*/
#if PY_MAJOR_VERSION >= 3

static int interface_python_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int interface_python_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_c_interface",
    NULL,
    sizeof(struct module_state),
    interface_python_methods,
    NULL,
    interface_python_traverse,
    interface_python_clear,
    NULL
};

#define INITERROR return NULL
#else
#define INITERROR return
#endif


#if PY_MAJOR_VERSION >= 3
    PyObject *
PyInit__c_interface(void)
#else
    void
init_c_interface(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("_c_interface", interface_python_methods);
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("pymess.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }
    init_numpy();
    mess_set_errorlevel(1);
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

/** @}  */
