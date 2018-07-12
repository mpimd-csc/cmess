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
 * @file python/c_interface/conv/vector.c
 * @brief Convert a @python matrix/array to a @c C matrix/array or vice-versa.
 *
 * @author @koehlerm
 * @author @nitin
 * @author @mbehr
 *
 *
 */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL MESS_VECTOR_MATRIX_PYTHON_C_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "../pymess.h"
#include "cscutils/error_message.h"

#if PY_MAJOR_VERSION < 3
#define ConvStringtoC(X) PyBytes_AsString(X)
#else
#define ConvStringtoC(X) PyBytes_AsString(PyUnicode_AsUTF8String(X))
#endif


#if NPY_API_VERSION < 0x00000007
#define NPY_ARRAY_OWNDATA NPY_OWNDATA
#define PyArray_ENABLEFLAGS(arr,flag)  ((PyArrayObject*)(arr))->flags |= (flag)
#endif



#define PyArray_CALL_ISSIGNED(x) PyArray_ISSIGNED((PyArrayObject*)(x))
#define PyArray_CALL_ISUNSIGNED(x) PyArray_ISUNSIGNED((PyArrayObject*)(x))
#define PyArray_CALL_DATA(x) PyArray_DATA((PyArrayObject*)(x))


/**
 * @brief Converts @python vector to @ref mess_vector.
 * @param[in] data input @python vector
 * @return converted vector as @ref mess_vector
 *
 * The @ref vector_to_c function supports the numpy ndarray and the following datatypes
 *  \li int8
 *  \li int16
 *  \li int32
 *  \li int64
 *  \li float32
 *  \li float64
 *  \li complex64
 *  \li complex128
 *
 */
mess_vector vector_to_c(PyObject* data) {
    PyObject *temp = NULL;
    PyObject *temp2 = NULL;
    char *dtype_str = NULL;
    mess_int_t i = 0 ;
    mess_vector c_vector;
    int ret = 0;
    npy_intp stride = 1 ;


    if(!PyArray_CheckExact(data))  {
        PyErr_SetString(PyExc_TypeError,"Argument type should be a 1-d numpy array");
        return NULL;
    }

    PY_GET_LONG(i, data ,"ndim");
    if (i != 1) {
        PyErr_SetString(PyExc_TypeError,"Argument type should be a 1-d numpy array");
        return NULL;
    }

    PY_GET_LONG(i, data, "size");
    temp = (PyObject*) PyObject_GetAttrString(data, "dtype");
    temp2 = (PyObject*) PyObject_GetAttrString(temp, "name");
    dtype_str = ConvStringtoC(temp2);

    if (strcmp(dtype_str,"float64") == 0 || strcmp(dtype_str, "float32") == 0 || PyArray_CALL_ISSIGNED(data) || PyArray_CALL_ISUNSIGNED(data)){
        PyObject *values = NULL ;
        ret = mess_vector_init(&c_vector);
        ret = mess_vector_alloc(c_vector, i, MESS_REAL);
        if ( ret != 0) {
            PyErr_SetString(PyExc_RuntimeError,
                    "Cannot allocate vector.");
            Py_DECREF(temp);
            Py_DECREF(temp2);
            return NULL;
        }
        if ( strcmp(dtype_str, "float64") != 0 ) {
            PyErr_Clear();
            values = PyArray_FROM_OT(data, NPY_DOUBLE);
        } else {
            values = data;
        }
        stride = PyArray_STRIDE((PyArrayObject *) values,0) / sizeof(double);
        if ( stride == 1 ) {
            memcpy(c_vector->values, PyArray_CALL_DATA(values), sizeof(double)*c_vector->dim);
        }  else {
            double *v = PyArray_CALL_DATA(values);
            for (i = 0; i < c_vector->dim; i++) {
                c_vector->values[i] = v[i*stride];
            }
        }
        if ( values != data ) {
            Py_DECREF(values);
        }
    }else if (strcmp(dtype_str,"complex128") == 0 || strcmp(dtype_str, "complex64") == 0){
        PyObject *values = NULL ;
        ret = mess_vector_init(&c_vector);
        ret = mess_vector_alloc(c_vector, i, MESS_COMPLEX);
        if ( ret != 0) {
            PyErr_SetString(PyExc_RuntimeError, "Cannot allocate vector.");
            Py_DECREF(temp);
            Py_DECREF(temp2);
            return NULL;
        }
        if ( strcmp(dtype_str, "complex128") != 0 ) {
            PyErr_Clear();
            values = PyArray_FROM_OT(data, NPY_COMPLEX128);
        } else {
            values = data;
        }
        stride = PyArray_STRIDE((PyArrayObject*) values,0) / sizeof(mess_double_cpx_t);
        if ( stride == 1 ) {
            memcpy(c_vector->values_cpx, PyArray_CALL_DATA(values), sizeof(mess_double_cpx_t)*c_vector->dim);
        }  else {
            mess_double_cpx_t *v = PyArray_CALL_DATA(values);
            for (i = 0; i < c_vector->dim; i++) {
                c_vector->values_cpx[i] = v[i*stride];
            }
        }
        if ( values != data ) {
            Py_DECREF(values);
        }
    }else  {
        Py_DECREF(temp);
        Py_DECREF(temp2);
        PyErr_SetString(PyExc_TypeError,"Argument Matrices should have either float64 or complex128 entries");
        return NULL;
    }
    //Clean up
    Py_DECREF(temp);
    Py_DECREF(temp2);

    return c_vector;
}

/**
 * @brief Converts mess_vector to @python vector.
 * @param[in] c_vector input @ref mess_vector
 * @param[in] copy     input set to true of the vector should be copied to @python
 * @return vector as @python Object
 *
 * The @ref vector_to_python function converts a mess_vector a numpy one. If the copy
 * argument is true the data is copied otherwise move to @python.
 *
 */
PyObject* vector_to_python(mess_vector c_vector, int copy)
{
    PyObject* ReturnArray = NULL;
    npy_intp dim;
    if ( copy ) {
        if(c_vector->data_type == MESS_REAL) {
            double * tmp ;
            dim = (npy_intp ) c_vector -> dim;
            ReturnArray = PyArray_SimpleNew(1, &dim, NPY_FLOAT64);
            tmp = PyArray_CALL_DATA(ReturnArray);
            memcpy(tmp, c_vector->values, sizeof(double) * dim);
        }
        else if(c_vector->data_type == MESS_COMPLEX) {
            mess_double_cpx_t *tmp;
            dim = (npy_intp ) c_vector -> dim;
            ReturnArray = PyArray_SimpleNew(1, &dim, NPY_COMPLEX128);
            tmp = PyArray_CALL_DATA(ReturnArray);
            memcpy(tmp, c_vector->values_cpx, sizeof(mess_double_cpx_t) * dim);
        }
        else {
            PyErr_SetString(PyExc_TypeError, "C-M.E.S.S. Vectors can only be Real or Complex");
            return NULL;
        }
    } else {
        if(c_vector->data_type == MESS_REAL) {
            dim = (npy_intp ) c_vector -> dim;
            ReturnArray = PyArray_SimpleNewFromData(1, &dim, NPY_FLOAT64, c_vector->values);
            PyArray_ENABLEFLAGS((PyArrayObject *)ReturnArray, NPY_ARRAY_OWNDATA);
            c_vector->values = NULL;
            c_vector->dim = 0;
        }
        else if(c_vector->data_type == MESS_COMPLEX) {
            dim = (npy_intp ) c_vector -> dim;
            ReturnArray = PyArray_SimpleNewFromData(1, &dim, NPY_COMPLEX128, c_vector->values_cpx);
            PyArray_ENABLEFLAGS((PyArrayObject *)ReturnArray, NPY_ARRAY_OWNDATA);
            c_vector->values_cpx = NULL;
            c_vector->dim = 0;
        }
        else {
            PyErr_SetString(PyExc_TypeError, "C-M.E.S.S. Vectors can only be Real or Complex");
            return NULL;
        }
    }
    return ReturnArray;
}

