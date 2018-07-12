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
 * @file python/c_interface/conv/eqn/so1.c
 * @brief Equation conversion for so1 function handles.
 * @author @mbehr
 *
 */

#include <Python.h>
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "../../pymess.h"

/*-----------------------------------------------------------------------------
 *  Convert a pymess.EquationGLyapSO1 to mess_equation_glyap_so1
 *-----------------------------------------------------------------------------*/
mess_equation eqn_conv_lyap_so1(PyObject *obj, mess_freelist mem) {
    mess_matrix M = NULL, D = NULL, K = NULL, B=NULL;
    mess_options opts = NULL;
    mess_equation eqn = NULL;
    double lowerbound=0.0, upperbound=0.0;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  get attributes
     *-----------------------------------------------------------------------------*/
    PyObject *pyM = PyObject_GetAttrString(obj,"m");
    PyObject *pyD = PyObject_GetAttrString(obj,"d");
    PyObject *pyK = PyObject_GetAttrString(obj,"k");
    PyObject *pyB = PyObject_GetAttrString(obj,"b");
    PyObject *pyOpt = PyObject_GetAttrString(obj,"options");
    PY_GET_DOUBLE(lowerbound, obj, "lowerbound");
    PY_GET_DOUBLE(upperbound, obj, "upperbound");

    if ( pyM == NULL || pyD == NULL || pyK == NULL || pyB==NULL || pyOpt ==NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot import data from the equation class.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  transfer data from Python to C-M.E.S.S.
     *-----------------------------------------------------------------------------*/
    M = matrix_to_c(pyM);
    if(!M){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix m from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, M);

    D = matrix_to_c(pyD);
    if(!D){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix d from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, D);

    K = matrix_to_c(pyK);
    if(!K){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix k from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, K);

    B = matrix_to_c(pyB);
    if(!B){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix b from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, B);

    opts = mess_options_from_python(pyOpt);
    if(!opts){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer options from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_options(mem,opts);

    /*-----------------------------------------------------------------------------
     *  create equation
     *-----------------------------------------------------------------------------*/
    mess_equation_init(&eqn);
    mess_freelist_add_mess_equation(mem,eqn);
    ret = mess_equation_glyap_so1(eqn, opts, M, D, K, B, lowerbound, upperbound);
    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"Cannot setup Lyapunov equation.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  clear and return
     *-----------------------------------------------------------------------------*/
    MESS_Py_XDECREF(pyM, pyD, pyK, pyB, pyOpt);
    return eqn;

error:
    MESS_Py_XDECREF(pyM, pyD, pyK, pyB, pyOpt);
    return NULL;
}

/*-----------------------------------------------------------------------------
 *  Convert a pymess.EquationGRiccatiSO1 to mess_equation_riccati_so1
 *-----------------------------------------------------------------------------*/
mess_equation eqn_conv_riccati_so1(PyObject *obj, mess_freelist mem) {
    mess_matrix M = NULL, D = NULL, K = NULL, B = NULL, C=NULL;
    mess_options opts = NULL;
    mess_equation eqn = NULL;
    double lowerbound=0, upperbound=0;
    int ret = 0;


    /*-----------------------------------------------------------------------------
     *  get attributes
     *-----------------------------------------------------------------------------*/
    PyObject *pyM = PyObject_GetAttrString(obj,"m");
    PyObject *pyD = PyObject_GetAttrString(obj,"d");
    PyObject *pyK = PyObject_GetAttrString(obj,"k");
    PyObject *pyB = PyObject_GetAttrString(obj,"b");
    PyObject *pyC = PyObject_GetAttrString(obj,"c");
    PyObject *pyOpt = PyObject_GetAttrString(obj,"options");
    PY_GET_DOUBLE(upperbound,obj,"upperbound");
    PY_GET_DOUBLE(lowerbound,obj,"lowerbound");

    if ( pyM == NULL || pyD == NULL || pyK == NULL || pyB == NULL || pyC == NULL || pyOpt == NULL ) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot import data from the equation class.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  transfer data from Python to C-M.E.S.S.
     *-----------------------------------------------------------------------------*/
    M = matrix_to_c(pyM);
    if(!M){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix m from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, M);

    D = matrix_to_c(pyD);
    if(!D){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix d from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, D);

    K = matrix_to_c(pyK);
    if(!K){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix k from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, K);

    B = matrix_to_c(pyB);
    if(!K){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix b from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, B);

    C = matrix_to_c(pyC);
    if(!C){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix c from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, C);

    opts = mess_options_from_python(pyOpt);
    if(!opts){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer options from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_options(mem,opts);

    /*-----------------------------------------------------------------------------
     *  create equation
     *-----------------------------------------------------------------------------*/
    mess_equation_init(&eqn);
    mess_freelist_add_mess_equation(mem,eqn);
    ret = mess_equation_griccati_so1(eqn, opts, M, D, K, B, C, lowerbound, upperbound);
    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"Cannot setup Riccati equation.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  clear and return
     *-----------------------------------------------------------------------------*/
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyM, pyD, pyK, pyB, pyC, pyOpt);
    return eqn;

error:
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyM, pyD, pyK, pyB, pyC, pyOpt);
    return NULL;
}

