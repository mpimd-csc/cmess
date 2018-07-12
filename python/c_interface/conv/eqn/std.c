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
 * @file python/c_interface/conv/eqn/std.c
 * @brief Equation conversion for standard lyapunov riccati equation.
 * @author @mbehr
 *
 */

#include <Python.h>
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "../../pymess.h"

/*-----------------------------------------------------------------------------
 *  Convert a pymess.EquationGLyap to mess_equation_lyap
 *-----------------------------------------------------------------------------*/
mess_equation eqn_conv_lyap(PyObject *obj, mess_freelist mem) {
    mess_matrix A = NULL, E = NULL, RHS = NULL;
    mess_options opts = NULL;
    mess_equation eqn=NULL;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  get attributes
     *-----------------------------------------------------------------------------*/
    PyObject *pyA = PyObject_GetAttrString(obj,"a");
    PyObject *pyE = PyObject_GetAttrString(obj,"e");
    PyObject *pyRHS = PyObject_GetAttrString(obj,"rhs");
    PyObject *pyOpt = PyObject_GetAttrString(obj,"options");


    if( pyA == NULL || pyE==NULL || pyRHS == NULL || pyOpt == NULL){
        PyErr_SetString(PyExc_RuntimeError, "Cannot import data from the equation class.");
        goto error;
    }


    /*-----------------------------------------------------------------------------
     *  transfer data from Python to C-M.E.S.S.
     *-----------------------------------------------------------------------------*/
    A = matrix_to_c(pyA);
    if (!A) {
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, A);

    if (pyE != Py_None){
        E = matrix_to_c(pyE);
        if(!E){
            PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix e from Python to C-M.E.S.S.");
            goto error;
        }
        mess_freelist_add_mess_matrix(mem, E);
    }


    RHS = matrix_to_c(pyRHS);
    if(!RHS){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix rhs from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, RHS);

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
    ret = mess_equation_lyap(eqn, opts, A, E, RHS);

    if ( ret != 0 ) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot setup Lyapunov equation.");
        goto error;
    }

    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyA, pyE, pyRHS, pyOpt);
    return eqn;
error:
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyA, pyE, pyRHS, pyOpt);
    return NULL;
}

/*-----------------------------------------------------------------------------
 *  Convert a pymess.EquationGRiccati to mess_equation_riccati
 *-----------------------------------------------------------------------------*/
mess_equation eqn_conv_riccati(PyObject *obj, mess_freelist mem) {
    mess_matrix A = NULL, E = NULL, B = NULL, C = NULL;
    mess_options opts = NULL;
    mess_equation eqn = NULL;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  get attributes
     *-----------------------------------------------------------------------------*/
    PyObject *pyA = PyObject_GetAttrString(obj,"a");
    PyObject *pyE = PyObject_GetAttrString(obj,"e");
    PyObject *pyB = PyObject_GetAttrString(obj,"b");
    PyObject *pyC = PyObject_GetAttrString(obj,"c");
    PyObject *pyOpt = PyObject_GetAttrString(obj,"options");


    if ( pyA == NULL || pyB == NULL || pyC == NULL || pyOpt == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot import data from the equation class.");
        goto error;
    }


    /*-----------------------------------------------------------------------------
     *  transfer data from Python to C-M.E.S.S.
     *-----------------------------------------------------------------------------*/
    A = matrix_to_c(pyA);
    if(!A){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, A);

    if (pyE != Py_None){
        E = matrix_to_c(pyE);
        if(!E){
            PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix e from Python to C-M.E.S.S.");
            goto error;
        }
        mess_freelist_add_mess_matrix(mem, E);
    }

    B = matrix_to_c(pyB);
    if(!B){
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
    ret = mess_equation_riccati(eqn, opts, A, E, B, C);

    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"Cannot setup Riccati equation.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  clear and return
     *-----------------------------------------------------------------------------*/
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyA, pyE, pyB, pyC, pyOpt);
    return eqn;

error:
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyA, pyE, pyB, pyC, pyOpt);
    return NULL;

}

