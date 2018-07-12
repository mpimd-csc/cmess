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
 * @file python/c_interface/functions/care.c
 * @brief Interface to @ref mess_care .
 *
 * @author @koehlerm
 * @author @nitin
 * @author @mbehr
 *
 *
 */

#include <Python.h>
#include "mess/mess.h"
#include "../pymess.h"

/**
 * @brief Wrapper around the easy to use interface \ref mess_care for use with @python.
 * @param[in] self input
 * @param[in] args input
 * @return @python object containing Z or None in case of an error and an exception.
 *
 * The @ref pymess_care function provides and wrapper around the \ref mess_care function  in order
 * to use from @python. The function is called as
 * \code[py]
 * Z = pymess.care(a, None, b, c )
 * \endcode
 * to solve the standard Riccati Equation \f$ A^TX+TXA-XBB^TX+C^TC=0 \f$ and as
 * \code[py]
 * Z = pymess.care(a, e, b, c)
 * \endcode
 * in case of the generalized Riccati Equation \f$ A^TXE+E^TXA-E^TXBB^TXE+C^TC= 0 \f$.
 *
 */
PyObject* pymess_care(PyObject *self, PyObject *args){
    int ret = 0;
    PyObject *mA = NULL;
    PyObject *mE = NULL;
    PyObject *mB = NULL;
    PyObject *mC = NULL;
    PyObject *mZ = NULL;


    mess_matrix A = NULL;
    mess_matrix E = NULL;
    mess_matrix B = NULL;
    mess_matrix C = NULL;
    mess_matrix Z = NULL;

    /*-----------------------------------------------------------------------------
     *  Parse the function header
     *-----------------------------------------------------------------------------*/
    PYCMESS_ERROR(!PyArg_ParseTuple(args, "OOOO", &mA, &mE, &mB, &mC),"The call sequence is wrong");

    /*-----------------------------------------------------------------------------
     *  prepare memlist
     *-----------------------------------------------------------------------------*/
    mess_freelist mem;
    mess_freelist_init(&mem);

    mess_matrix_init(&Z);
    mess_freelist_add_mess_matrix(mem,Z);

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

    C  = matrix_to_c(mC);
    if (!B) {
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix c from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem,C);

    if (mE != NULL && mE != Py_None) {
        E = matrix_to_c(mE);
        if (!E) {
            PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix e from Python to C-M.E.S.S.");
            goto error;
        }
        mess_freelist_add_mess_matrix(mem,E);
    }

    /*-----------------------------------------------------------------------------
     *  solve equation
     *-----------------------------------------------------------------------------*/
    ret = mess_care(A, E, B, C, Z);
    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"C-M.E.S.S.: mess_care returned an error.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     * convert result to Python
     *-----------------------------------------------------------------------------*/
    mZ = matrix_to_python(Z);

    /*-----------------------------------------------------------------------------
     *  clear data and return
     *-----------------------------------------------------------------------------*/
    mess_freelist_clear(&mem);
    return mZ;

    /*-----------------------------------------------------------------------------
     *  error handling
     *-----------------------------------------------------------------------------*/
error:
    mess_freelist_clear(&mem);
    return NULL;
}

