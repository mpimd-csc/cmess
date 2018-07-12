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
 * @file python/c_interface/conv/eqn/dae2.c
 * @brief Equation conversion for dae2 function handles.
 * @author @mbehr
 *
 */


#include <Python.h>
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "../../pymess.h"

/*-----------------------------------------------------------------------------
 *  Convert a pymess.EquationGLyapDAE2 to mess_equation_glyap_dae2
 *-----------------------------------------------------------------------------*/
mess_equation eqn_conv_lyap_dae2(PyObject *obj, mess_freelist mem) {
    MSG_FNAME(__func__);
    mess_matrix M = NULL, A = NULL, G = NULL, B=NULL, K0=NULL;
    mess_options opts = NULL;
    mess_equation eqn = NULL, eqn_stable = NULL;
    double delta =0;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  get attributes
     *-----------------------------------------------------------------------------*/
    PyObject *pyM = PyObject_GetAttrString(obj,"m");
    PyObject *pyA = PyObject_GetAttrString(obj,"a");
    PyObject *pyG = PyObject_GetAttrString(obj,"g");
    PyObject *pyB = PyObject_GetAttrString(obj,"b");
    PyObject *pyK0 = PyObject_GetAttrString(obj,"k0");
    PyObject *pyOpt = PyObject_GetAttrString(obj,"options");
    PY_GET_DOUBLE(delta,obj,"delta");

    if ( pyM == NULL || pyA == NULL || pyG == NULL || pyB == NULL || pyK0 == NULL || pyOpt == NULL) {
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

    A = matrix_to_c(pyA);
    if(!A){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, A);

    G = matrix_to_c(pyG);
    if(!G){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix g from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, G);

    B = matrix_to_c(pyB);
    if(!B){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix b from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, B);

    if (pyK0 != Py_None){
        K0 = matrix_to_c(pyK0);
        if(!K0){
            PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix k0 from Python to C-M.E.S.S.");
            goto error;
        }
        mess_freelist_add_mess_matrix(mem, K0);
    }

    opts = mess_options_from_python(pyOpt);
    if(!opts){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer options from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_options(mem,opts);

    /*-----------------------------------------------------------------------------
     *  create equation
     *-----------------------------------------------------------------------------*/
    mess_equation_init(&eqn_stable);
    mess_equation_init(&eqn);
    ret = mess_equation_glyap_dae2(eqn, opts, M, A, G, B, delta);
    //it is important that eqn_stable is inserted before eqn, because of the cleaning method afterwards
    mess_freelist_add_mess_equation(mem, eqn_stable);
    mess_freelist_add_mess_equation(mem, eqn);

    if ( ret != 0 ) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot setup Lyapunov equation.");
        goto error;
    }


    /*-----------------------------------------------------------------------------
     *  stabilize equation if K0 was given
     *-----------------------------------------------------------------------------*/
    if(K0){
        MSG_INFO("Stabilize system with given initial feedback.");
        if(opts->type == MESS_OP_NONE){
            ret = mess_equation_stable(eqn_stable,opts,eqn,eqn->B,K0);
        }else{
            ret = mess_equation_stable(eqn_stable,opts,eqn,K0,eqn->B);
        }
        if (ret){
            PyErr_SetString(PyExc_RuntimeError,"Cannot stabilize Lyapunov Equation with given initial feedback.");
            goto error;
        }
    }

    /*-----------------------------------------------------------------------------
     *  clear and return
     *-----------------------------------------------------------------------------*/
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyM, pyA, pyG, pyB, pyK0, pyOpt);
    return K0?eqn_stable:eqn;

error:
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyM, pyA, pyG, pyB, pyK0, pyOpt);
    return NULL;
}


/*-----------------------------------------------------------------------------
 *  Convert a pymess.EquationGRiccatiDAE2 to mess_equation_riccati_dae2
 *-----------------------------------------------------------------------------*/
mess_equation eqn_conv_riccati_dae2(PyObject *obj, mess_freelist mem) {
    mess_matrix M = NULL, A = NULL, G = NULL, B = NULL, C=NULL ;
    mess_options opts = NULL;
    mess_equation eqn = NULL;
    double delta=0;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  get attributes
     *-----------------------------------------------------------------------------*/
    PyObject *pyM = PyObject_GetAttrString(obj,"m");
    PyObject *pyA = PyObject_GetAttrString(obj,"a");
    PyObject *pyG = PyObject_GetAttrString(obj,"g");
    PyObject *pyB = PyObject_GetAttrString(obj,"b");
    PyObject *pyC = PyObject_GetAttrString(obj,"c");
    PyObject *pyOpt = PyObject_GetAttrString(obj,"options");
    PY_GET_DOUBLE(delta,obj,"delta");

    if ( pyM == NULL || pyA == NULL || pyG == NULL || pyB == NULL || pyC == NULL ||  pyOpt == NULL ) {
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

    A = matrix_to_c(pyA);
    if(!A){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix a from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, A);

    G = matrix_to_c(pyG);
    if(!G){
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix g from Python to C-M.E.S.S.");
        goto error;
    }
    mess_freelist_add_mess_matrix(mem, G);

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
    ret = mess_equation_griccati_dae2(eqn, opts, M,A,G,B,C,delta);

    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"Cannot setup Riccati equation.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  clear and return
     *-----------------------------------------------------------------------------*/
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyM, pyA, pyG, pyB, pyC, pyOpt);
    return eqn;

error:
    //we do not clear mem, because it was not allocated by us
    MESS_Py_XDECREF(pyM, pyA, pyG, pyB, pyC, pyOpt);
    return NULL;
}


