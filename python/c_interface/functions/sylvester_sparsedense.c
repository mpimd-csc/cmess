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
 * @file python/c_interface/functions/sylvester_sparsedense.c
 * @brief Interface to @ref mess_direct_create_sylvester_sparsedense.
 *
 * @author @koehlerm
 * @author @dykstra
 *
 */

#include <Python.h>
#include "mess/mess.h"
#include "../pymess.h"


/**
 * @brief Wrapper around the interface \ref mess_direct_create_sylvester_sparsedense for use with @python.
 * @param[in]  self input
 * @param[in]  args input
 * @return @python object containing X or None in case of an error / exception.
 *
 * The @ref pymess_sylvester_sparsedense function provides a wrapper around the \ref mess_direct_create_sylvester_sparsedense
 * function in order to use from @python. The function is called as
 * \code[py]
 * Z = pymess.sylvester_sparsedense(a, None, None, b, m)
 * \endcode
 * to solve the standard Sylvester Equation \f$ AX+XB+M=0 \f$, as
 * \code[py]
 * Z = pymess.sylvester_sparsedense(a, None, e, b, m)
 * \endcode
 * to solve the semi-generalized Sylvester Equation \f$ AX+EXB+M=0 \f$ and as
 * \code[py]
 * Z = pymess.sylvester_sparsedense(a, f, e, b, m)
 * \endcode
 * to solve the generalized Sylvester Equation \f$ AXF+EXB+M=0 \f$.
 *
 */
PyObject* pymess_sylvester_sparsedense(PyObject * self, PyObject * args){
    //PyObject* pymess_lyap(PyObject *self, PyObject *args, PyObject *kwrds){
    int ret = 0;
    PyObject *mA = NULL;
    PyObject *mF = NULL;
    PyObject *mE = NULL;
    PyObject *mB = NULL;
    PyObject *mM = NULL;
    PyObject *mX = NULL;

    mess_matrix A = NULL;
    mess_matrix F = NULL;
    mess_matrix E = NULL;
    mess_matrix B = NULL;
    mess_matrix M = NULL;
    mess_matrix X = NULL;
    
    mess_direct solver = NULL;
    mess_direct_init(&solver);

    /*-----------------------------------------------------------------------------
     *  Parse the function header
     *-----------------------------------------------------------------------------*/
    PYCMESS_ERROR(!PyArg_ParseTuple(args, "OOOOO", &mA, &mF, &mE, &mB, &mM),"The call sequence is wrong.");


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

    B  = matrix_to_c(mB);
    if (!B) {
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix b from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem,B);

    M  = matrix_to_c(mM);
    if (!M) {
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix m from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem,M);


    if (mE != NULL && mE != Py_None) {
        E = matrix_to_c(mE);
        if (!E) {
            PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix e from Python to C-M.E.S.S.");
            goto error;
        }
        mess_freelist_add_mess_matrix(mem,E);
    }
    

    if (mF != NULL && mF != Py_None) {
        F = matrix_to_c(mF);
        if (!F) {
            PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix f from Python to C-M.E.S.S.");
            goto error;
        }
        mess_freelist_add_mess_matrix(mem,F);
    }

    /*-----------------------------------------------------------------------------
     *  solve equation
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_create_sylvester_sparsedense(A,F,E,B,solver);
    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"C-M.E.S.S.: mess_direct_create_sylvester_sparsedense returned an error.");
        goto error;
    }
    ret = mess_direct_solvem(MESS_OP_NONE, solver, M, X);
    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"C-M.E.S.S.: mess_direct_solvem returned an error.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     * convert result to Python
     *-----------------------------------------------------------------------------*/
    mX = matrix_to_python(X);

    /*-----------------------------------------------------------------------------
     *  clear data and return
     *-----------------------------------------------------------------------------*/
    mess_freelist_clear(&mem);
    mess_direct_clear(&solver);
    return mX;

    /*-----------------------------------------------------------------------------
     *  error handling
     *-----------------------------------------------------------------------------*/
error:
    mess_freelist_clear(&mem);
    return NULL;
}



