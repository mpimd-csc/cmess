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
 * @file python/c_interface/conv/eqn/dae1.c
 * @brief Equation conversion for dae1 function handles.
 * @author @mbehr
 *
 */


#include <Python.h>
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "../../pymess.h"

/*-----------------------------------------------------------------------------
 *  Convert a pymess.EquationGLyapDAE1 to mess_equation_glyap_dae1
 *-----------------------------------------------------------------------------*/
mess_equation eqn_conv_lyap_dae1(PyObject *obj, mess_freelist mem) {
    mess_matrix E11 = NULL, A11 = NULL, A12 = NULL, A21 = NULL, A22 = NULL, B = NULL;
    mess_options opts = NULL;
    mess_equation eqn = NULL;
    mess_int_t ret = 0;


    /*-----------------------------------------------------------------------------
     *  get attributes
     *-----------------------------------------------------------------------------*/
    PyObject *pyE11 = PyObject_GetAttrString(obj,"e11");
    PyObject *pyA11 = PyObject_GetAttrString(obj,"a11");
    PyObject *pyA12 = PyObject_GetAttrString(obj,"a12");
    PyObject *pyA21 = PyObject_GetAttrString(obj,"a21");
    PyObject *pyA22 = PyObject_GetAttrString(obj,"a22");
    PyObject *pyB   = PyObject_GetAttrString(obj,"b");
    PyObject *pyOpt = PyObject_GetAttrString(obj,"options");

    if ( pyE11 == NULL || pyA11 == NULL || pyA12 == NULL || pyA21 == NULL || pyA22 == NULL || pyB == NULL || pyOpt == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot import data from the equation class.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  transfer data from Python to C-M.E.S.S.
     *-----------------------------------------------------------------------------*/
    E11 = matrix_to_c(pyE11);
    if(!E11){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix e11 from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, E11);

    A11 = matrix_to_c(pyA11);
    if(!A11){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a11 from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, A11);

    A12 = matrix_to_c(pyA12);
    if(!A12){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a12 from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, A12);

    A21 = matrix_to_c(pyA21);
    if(!A21){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a21 from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, A21);

    A22 = matrix_to_c(pyA22);
    if(!A22){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a22 from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, A22);

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
    ret = mess_equation_glyap_dae1(eqn, opts, E11, A11, A12, A21, A22, B);

    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"Cannot setup Lyapunov equation.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  clear and return
     *-----------------------------------------------------------------------------*/
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyE11, pyA11, pyA12, pyA21, pyA22, pyB, pyOpt);
    return eqn;

error:
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyE11, pyA11, pyA12, pyA21, pyA22, pyB, pyOpt);
    return NULL;
}


/*-----------------------------------------------------------------------------
 *  Convert a pymess.EquationGRiccatiDAE1 to mess_equation_riccati_dae1
 *-----------------------------------------------------------------------------*/
mess_equation eqn_conv_riccati_dae1(PyObject *obj, mess_freelist mem) {
    mess_matrix E11 = NULL, A11 = NULL, A12 = NULL, A21 = NULL, A22 = NULL, B = NULL, C = NULL;
    mess_options opts = NULL;
    mess_equation eqn = NULL;
    mess_int_t ret = 0;

    /*-----------------------------------------------------------------------------
     *  get attributes
     *-----------------------------------------------------------------------------*/
    PyObject *pyE11 = PyObject_GetAttrString(obj,"e11");
    PyObject *pyA11 = PyObject_GetAttrString(obj,"a11");
    PyObject *pyA12 = PyObject_GetAttrString(obj,"a12");
    PyObject *pyA21 = PyObject_GetAttrString(obj,"a21");
    PyObject *pyA22 = PyObject_GetAttrString(obj,"a22");
    PyObject *pyB   = PyObject_GetAttrString(obj,"b");
    PyObject *pyC   = PyObject_GetAttrString(obj,"c");
    PyObject *pyOpt = PyObject_GetAttrString(obj,"options");

    if (pyE11==NULL || pyA11==NULL || pyA12==NULL || pyA21==NULL || pyA22==NULL || pyB==NULL || pyC==NULL ||  pyOpt == NULL ){
        PyErr_SetString(PyExc_RuntimeError, "Cannot import data from the equation class.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  transfer data from Python to C-M.E.S.S.
     *-----------------------------------------------------------------------------*/
    E11 = matrix_to_c(pyE11);
    if(!E11){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix e11 from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, E11);

    A11 = matrix_to_c(pyA11);
    if(!A11){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a11 from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, A11);

    A12 = matrix_to_c(pyA12);
    if(!A12){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a12 from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, A12);

    A21 = matrix_to_c(pyA21);
    if(!A21){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a21 from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, A21);

    A22 = matrix_to_c(pyA22);
    if(!A22){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a22 from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, A22);

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
    ret = mess_equation_griccati_dae1(eqn, opts, E11, A11, A12, A21, A22, B, C);
    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"Cannot setup Riccati equation.");
        goto error;
    }



    /*-----------------------------------------------------------------------------
     *  clear and return
     *-----------------------------------------------------------------------------*/
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyE11, pyA11, pyA12, pyA21, pyA22, pyB, pyC, pyOpt);
    return eqn;

error:
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyE11, pyA11, pyA12, pyA21, pyA22, pyB, pyC, pyOpt);
    return NULL;
}

