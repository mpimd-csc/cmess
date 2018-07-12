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
 * @file python/c_interface/functions/dense_nm_gmpare.c
 * @brief Interface to @ref mess_dense_nm_gmpare.
 *
 * @author @koehlerm
 * @author @dykstra
 * @author @mbehr
 *
 */

#include <Python.h>
#include "mess/mess.h"
#include "../pymess.h"


/**
 * @brief Wrapper around the interface @ref mess_dense_nm_gmpare for use with @python.
 * @param[in]  self input
 * @param[in]  args input
 * @return @python object containing X or None in case of an error / exception.
 *
 * The @ref pymess_dense_nm_gmpare function provides a wrapper around the \ref mess_dense_nm_gmpare
 * function in order to use from @python.
 *
 */
PyObject* pymess_dense_nm_gmpare(PyObject * self, PyObject * args){
    int ret = 0;
    PyObject *mA = NULL;
    PyObject *mX0 = NULL;
    PyObject *mE = NULL;
    PyObject *mQ = NULL;
    PyObject *mG = NULL;
    PyObject *mX = NULL;
    PyObject *pyResult = NULL;

    mess_matrix A = NULL;
    mess_matrix X0 = NULL;
    mess_matrix E = NULL;
    mess_matrix Q = NULL;
    mess_matrix G = NULL;
    mess_matrix X = NULL;
    double absres_tol, relres_tol, absres, relres;
    mess_int_t plus, linesearch, maxit;
    mess_operation_t trans;
    mess_norm_t nrm;


    /*-----------------------------------------------------------------------------
     *  Parse the function header
     *-----------------------------------------------------------------------------*/
    PYCMESS_ERROR(!PyArg_ParseTuple(args, "OOOOOiiiiddi", &mX0, &mA, &mE, &mQ, &mG,
                &plus, &linesearch, &trans, &maxit, &absres_tol, &relres_tol, &nrm),"The call sequence is wrong.");

    /*-----------------------------------------------------------------------------
     *  prepare memlist
     *-----------------------------------------------------------------------------*/
    mess_freelist mem;
    mess_freelist_init(&mem);

    mess_matrix_init(&X);
    mess_freelist_add_mess_matrix(mem,X);

    /*-----------------------------------------------------------------------------
     *  tranfer data
     *-----------------------------------------------------------------------------*/
    A  = matrix_to_c(mA);
    if (!A) {
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem,A);

    Q  = matrix_to_c(mQ);
    if (!Q) {
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix q from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem,Q);

    G  = matrix_to_c(mG);
    if (!G) {
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix g from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem,G);


    if (mE != NULL && mE != Py_None) {
        E = matrix_to_c(mE);
        if (!E) {
            PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix e from Python to C-M.E.S.S.");
            goto error;
        }
        mess_freelist_add_mess_matrix(mem,E);
    }


    if (mX0 != NULL && mX0 != Py_None) {
        X0 = matrix_to_c(mX0);
        if (!X0) {
            PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix x0 from Python to C-M.E.S.S.");
            goto error;
        }
        mess_freelist_add_mess_matrix(mem,X0);
    }

    /*-----------------------------------------------------------------------------
     *  solve equation
     *-----------------------------------------------------------------------------*/
    ret = mess_dense_nm_gmpare(X0, A, E, Q, G, plus, linesearch, trans, maxit, nrm, absres_tol, relres_tol, &absres, &relres, NULL, X);
    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"C-M.E.S.S.: mess_dense_nm_gmpare returned an error.");
        goto error;
    }


    /*-----------------------------------------------------------------------------
     * convert result to Python
     *-----------------------------------------------------------------------------*/
    mX = matrix_to_python(X);
    pyResult = Py_BuildValue("Odd", mX, absres, relres);

    if(!pyResult){
        Py_XDECREF(mX);
        PyErr_SetString(PyExc_RuntimeError, "Cannot build result");
        goto error;
    }else{
        Py_XDECREF(mX);
    }


    /*-----------------------------------------------------------------------------
     *  clear data and return
     *-----------------------------------------------------------------------------*/
    mess_freelist_clear(&mem);
    return pyResult;

    /*-----------------------------------------------------------------------------
     *  error handling
     *-----------------------------------------------------------------------------*/
error:
    mess_freelist_clear(&mem);
    return NULL;
}



