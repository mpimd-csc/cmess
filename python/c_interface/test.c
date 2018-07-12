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

#include <Python.h>
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "pymess.h"


/**
 * @file python/c_interface/test.c
 * @brief Test functions which can be called from @python.
 * @author @nitin
 * @author @mbehr
 */


/**
 * @brief Test functions that can be called from @python to ensure that all the conversions are working.
 *
 */
PyObject* pymess_test_matrix(PyObject *self, PyObject* args){
    int ret = 0;
    PyObject* data          = NULL;
    PyObject* Return_data   = NULL;

    mess_matrix c_matrix    = NULL;
    mess_matrix c_matrix_T  = NULL;

    /*-----------------------------------------------------------------------------
     *  parse input arguments
     *-----------------------------------------------------------------------------*/
    PYCMESS_ERROR(!PyArg_ParseTuple(args, "O", &data),"The call sequence is wrong.");

    /*-----------------------------------------------------------------------------
     *  prepare memlist
     *-----------------------------------------------------------------------------*/
    mess_freelist mem;
    mess_freelist_init(&mem);

    /*-----------------------------------------------------------------------------
     * transfer data from Python to C-M.E.S.S.
     *-----------------------------------------------------------------------------*/
    c_matrix = matrix_to_c(data);
    if (!c_matrix){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem,c_matrix);


    /*-----------------------------------------------------------------------------
     *  compute conjugate transpose
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&c_matrix_T);
    mess_freelist_add_mess_matrix(mem,c_matrix_T);

    ret = mess_matrix_ctranspose(c_matrix, c_matrix_T);
    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"C-M.E.S.S.: mess_matrix_ctranspose returned an error.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  transfer data from C-M.E.S.S. to Python
     *-----------------------------------------------------------------------------*/
    Return_data = matrix_to_python(c_matrix_T);

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    mess_freelist_clear(&mem);
    return Return_data;

error:
    mess_freelist_clear(&mem);
    return NULL;
}

/**
 * @brief Test functions that can be called from @python to ensure that all the conversions are working.
 *
 */
PyObject* pymess_test_vector(PyObject *self, PyObject* args){
    PyObject* data          = NULL;
    PyObject* Return_data   = NULL;

    mess_vector c_vector = NULL;

    /*-----------------------------------------------------------------------------
     *  parse input arguments
     *-----------------------------------------------------------------------------*/
    PYCMESS_ERROR(!PyArg_ParseTuple(args, "O", &data),"The call sequence is wrong.");

    /*-----------------------------------------------------------------------------
     *  prepare memlist
     *-----------------------------------------------------------------------------*/
    mess_freelist mem;
    mess_freelist_init(&mem);

    /*-----------------------------------------------------------------------------
     * transfer data from Python to C-M.E.S.S.
     *-----------------------------------------------------------------------------*/
    c_vector = vector_to_c(data);
    if(!c_vector){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer vector from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_vector(mem,c_vector);

    /*-----------------------------------------------------------------------------
     * transfer data from C-M.E.S.S. to Python
     *-----------------------------------------------------------------------------*/
    Return_data = vector_to_python(c_vector,0);

    mess_freelist_clear(&mem);
    return Return_data;

error:
    mess_freelist_clear(&mem);
    return NULL;
}



PyObject * pymess_test_equation(PyObject *self, PyObject *args) {
    PyObject* return_value  = NULL;
    PyObject* data          = NULL;
    PyObject *pymess_module = NULL;
    PyObject *pymess_eqn_base = NULL;
    int ret = 0;

    if (!PyArg_ParseTuple(args, "O", &data))
        Py_RETURN_NONE;

    pymess_module   = PyImport_ImportModule((char *) "pymess");
    pymess_eqn_base = PyObject_GetAttrString(pymess_module,(char*)"Equation");

    if ( pymess_module == NULL || pymess_eqn_base == NULL) {
        Py_XDECREF(pymess_eqn_base);
        Py_XDECREF(pymess_module);
        PyErr_SetString(PyExc_Exception, "Cannot import pymess module or the equation class.");
        return NULL;
    }

    ret = PyObject_IsInstance(data, pymess_eqn_base);

    if ( ret == 1) {
        MSG_PRINT("Input is instance or subclass of pymess\n");
    } else if ( ret == 0  ){
        MSG_PRINT("Input is not an instance\n");
        return_value = Py_None;
        goto fail;
    } else {
        MSG_PRINT("An error occured.\n");
        return_value  = Py_None;
        goto fail;
    }
    return_value = Py_None;
    return_value =  PyObject_CallMethod(data, (char*) "_test", "ii", 1, 2);
fail:
    Py_XDECREF(pymess_eqn_base);
    Py_XDECREF(pymess_module);

    return return_value;
}

