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
 * @file python/c_interface/conv/options.c
 * @brief Convert a @python options class to @ref mess_options and vice-versa.
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
 * @brief Parse the attributes for the lrcfadi algorithm from PyObject options argument.
 * @param[in] options input PyObject with opt for lrcfadi algorithm
 * @return parsed mess_options from PyObject options
 *
 *
 */
mess_options mess_options_from_python(PyObject *options) {
    mess_options opt = NULL;
    PyObject *adi_opt = NULL;
    PyObject *adi_shifts_opt = NULL;
    PyObject *nm_opt = NULL;
    int ret = 0;
    int optype=0;
    int pt = 0;

    /*-----------------------------------------------------------------------------
     *  init options
     *-----------------------------------------------------------------------------*/
    ret = mess_options_init(&opt);
    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"Cannot initialize mess_options structure.");
        return NULL;
    }

    /*-----------------------------------------------------------------------------
     *  Get Substructures
     *-----------------------------------------------------------------------------*/
    adi_opt         = PyObject_GetAttrString(options,"adi");
    adi_shifts_opt  = PyObject_GetAttrString(adi_opt, "shifts");
    nm_opt          = PyObject_GetAttrString(options,"nm");
    if(adi_opt == NULL || adi_shifts_opt == NULL || nm_opt == NULL){
        PyErr_SetString(PyExc_TypeError, "One of the option sub structures adi, adi.shift or nm is missing");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  get operation type
     *-----------------------------------------------------------------------------*/
    opt->type = MESS_OP_NONE;
    PY_GET_LONG(optype,options,"type");
    if(optype==0){
        opt->type=MESS_OP_NONE;
    }else if(optype==1){
        opt->type=MESS_OP_TRANSPOSE;
    }else if(optype==2){
        opt->type=MESS_OP_HERMITIAN;
    }else{
        PyErr_SetString(PyExc_TypeError, "Unknown transpose operations in type");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  get residual method
     *-----------------------------------------------------------------------------*/
    PY_GET_LONG(opt->residual_method,options,"residual_method");

    /*-----------------------------------------------------------------------------
     *  Extract Values from ADI
     *-----------------------------------------------------------------------------*/
    PY_GET_LONG(opt->adi_maxit, adi_opt, "maxit");
    PY_GET_DOUBLE(opt->adi_res2_tol,adi_opt,"res2_tol");
    PY_GET_DOUBLE(opt->adi_res2c_tol,adi_opt,"res2c_tol");
    PY_GET_DOUBLE(opt->adi_rel_change_tol,adi_opt,"rel_change_tol");
    PY_GET_LONG(opt->adi_output,adi_opt, "output");
    PY_GET_LONG(opt->memory_usage,adi_opt, "memory_usage");

    /*-----------------------------------------------------------------------------
     *  Extract Values from ADI.SHIFTS
     *-----------------------------------------------------------------------------*/
    PY_GET_VECTOR(opt->adi_shifts_p, adi_shifts_opt, "p");
    PY_GET_LONG(opt->adi_shifts_arp_p, adi_shifts_opt, "arp_p");
    PY_GET_LONG(opt->adi_shifts_arp_m, adi_shifts_opt, "arp_m");
    PY_GET_LONG(opt->adi_shifts_l0, adi_shifts_opt, "l0");
    PY_GET_LONG(pt, adi_shifts_opt, "paratype");
    switch(pt) {
        case 0:
            opt->adi_shifts_paratype = MESS_LRCFADI_PARA_MINMAX;
            break;
        case 1:
            opt->adi_shifts_paratype = MESS_LRCFADI_PARA_MINMAX_REAL;
            break;
        case 2:
            opt->adi_shifts_paratype = MESS_LRCFADI_PARA_WACHSPRESS;
            break;
        case 3:
            opt->adi_shifts_paratype = MESS_LRCFADI_PARA_ADAPTIVE_V;
            break;
        case 4:
            opt->adi_shifts_paratype = MESS_LRCFADI_PARA_ADAPTIVE_Z;
            break;
        default:
            PyErr_SetString(PyExc_Exception, "The requested parameter type is not supported.");
            goto error;
    }
    PY_GET_VECTOR(opt->adi_shifts_b0, adi_shifts_opt, "b0");

    /*-----------------------------------------------------------------------------
     *  Extract Values from NM
     *-----------------------------------------------------------------------------*/
    PY_GET_LONG(opt->nm_maxit, nm_opt, "maxit");
    PY_GET_DOUBLE(opt->nm_res2_tol,nm_opt,"res2_tol");
    PY_GET_LONG(opt->nm_gpStep,nm_opt, "gpstep");
    PY_GET_LONG(opt->nm_output,nm_opt, "output");
    PY_GET_LONG(opt->nm_singleshifts,nm_opt, "singleshifts");
    PY_GET_LONG(opt->nm_linesearch,nm_opt, "linesearch");
    PY_GET_MATRIX(opt->K0,nm_opt,"k0");

    /*-----------------------------------------------------------------------------
     *  clear and error
     *-----------------------------------------------------------------------------*/
    MESS_Py_XDECREF(adi_opt, adi_shifts_opt, nm_opt);
    return opt;

error:
    MESS_Py_XDECREF(adi_opt, adi_shifts_opt, nm_opt);
    mess_options_clear(&opt);
    return NULL;
}



