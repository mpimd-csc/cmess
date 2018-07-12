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
 * @file python/c_interface/functions/glyap.c
 * @brief Interface to @ref mess_glyap.
 *
 * @author @mbehr
 *
 */

#include <Python.h>
#include "mess/mess.h"
#include "../pymess.h"


/**
 * @brief Interface to \ref mess_glyap for use with @python.
 * @param[in]  self input
 * @param[in]  args input
 * @return @python object containing X and optional the schur decomposition or None in case of an error / exception.
 *
 * The @ref pymess_glyap function provides an interface to \ref mess_glyap function  in order
 * to use from @python. The function is called as
 *
 */
PyObject* pymess_glyap(PyObject * self, PyObject * args){
    unsigned int optype = 0;
    unsigned int schur_given = 0;
    mess_operation_t op = MESS_OP_NONE;
    PyObject *mA = NULL, *mE = NULL, *mY = NULL, *mAhat = NULL, *mEhat = NULL, *mQA = NULL, *mQE = NULL;
    PyObject *resX = NULL, *resAhat = NULL, *resEhat = NULL, *resQA = NULL, *resQE = NULL, *result = NULL;
    mess_matrix A = NULL, E = NULL, Y = NULL, Ahat = NULL, Ehat = NULL, QA = NULL, QE = NULL, X = NULL;

    /*-----------------------------------------------------------------------------
     *  Parse the function header
     *-----------------------------------------------------------------------------*/
    PYCMESS_ERROR(!PyArg_ParseTuple(args, "OOOIOOOO", &mA, &mE, &mY, &optype, &mAhat, &mEhat, &mQA, &mQE),"The call sequence is wrong.");

    /*-----------------------------------------------------------------------------
     *  prepare memlist
     *-----------------------------------------------------------------------------*/
    mess_freelist mem;
    mess_freelist_init(&mem);

    mess_matrix_init(&X);
    mess_freelist_add_mess_matrix(mem,X);

    /*-----------------------------------------------------------------------------
     *  set operation type
     *-----------------------------------------------------------------------------*/
    op = optype?MESS_OP_TRANSPOSE:MESS_OP_NONE;
    if(optype >= 3){
        PyErr_SetString(PyExc_TypeError, "Unknown operation  in op is given.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  tranfer data
     *-----------------------------------------------------------------------------*/
    PYMESS_READ_MATRIX(mem,A,mA,a,error);
    PYMESS_READ_MATRIX(mem,Y,mY,y,error);
    if (mE != NULL && mE != Py_None) {
        PYMESS_READ_MATRIX(mem,E,mE,e,error);
    }

    /*-----------------------------------------------------------------------------
     *  check for given schur decomposition is fully given or not
     *-----------------------------------------------------------------------------*/
    if(mAhat != Py_None && mQA != Py_None && ( !E ||  (mEhat != Py_None && mQE != Py_None))){
        /*-----------------------------------------------------------------------------
         *  read given schur decomposition
         *-----------------------------------------------------------------------------*/
        PYMESS_READ_MATRIX(mem,Ahat,mAhat,ahat,error);
        PYMESS_READ_MATRIX(mem,QA,mQA,qa,error);
        if(E){
            PYMESS_READ_MATRIX(mem,Ehat,mEhat,ehat,error);
            PYMESS_READ_MATRIX(mem,QE,mQE,qe,error);
        }
        schur_given = 1;
    }else if(mAhat == Py_None && mQA == Py_None && ( !E ||  (mEhat == Py_None && mQE == Py_None))){
        /*-----------------------------------------------------------------------------
         *  schur decomposition is computed and returned to PY-M.E.S.S.
         *-----------------------------------------------------------------------------*/
        if(E){
            MESS_INIT_MATRICES(&Ahat,&Ehat,&QA,&QE);
        }else{
            MESS_INIT_MATRICES(&Ahat,&QA);
        }
        schur_given = 0;
    }else{
        PyErr_SetString(PyExc_RuntimeError,"Either all Schur transformation matrices must be given or none.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  solve equation
     *-----------------------------------------------------------------------------*/
    if(schur_given){
        if(mess_tglyap(op, Ahat, QA, Ehat, QE, Y, X)){
            PyErr_SetString(PyExc_RuntimeError,"C-M.E.S.S.: mess_tglyap returned an error.");
            goto error;
        }
    }else{
        if(mess_glyap(op, A, E, Y, Ahat, QA, Ehat, QE, X)){
            PyErr_SetString(PyExc_RuntimeError,"C-M.E.S.S.: mess_glyap returned an error.");
            goto error;
        }
    }

    /*-----------------------------------------------------------------------------
     * convert result to Python
     *-----------------------------------------------------------------------------*/
    if(schur_given){
        result = matrix_to_python(X);
    }else{
        resX    = matrix_to_python(X);
        resAhat = matrix_to_python(Ahat);
        resQA   = matrix_to_python(QA);
        if(E){
            resEhat = matrix_to_python(Ehat);
            resQE   = matrix_to_python(QE);
            result  = Py_BuildValue("OOOOO", resX, resAhat, resEhat, resQA, resQE);
            Py_DECREF(resEhat);
            Py_DECREF(resQE);
        }else{
            result  = Py_BuildValue("OOO", resX, resAhat, resQA);
        }
        Py_DECREF(resX);
        Py_DECREF(resAhat);
        Py_DECREF(resQA);
    }

    /*-----------------------------------------------------------------------------
     *  clear data and return
     *-----------------------------------------------------------------------------*/
    mess_freelist_clear(&mem);
    return result;

    /*-----------------------------------------------------------------------------
     *  error handling
     *-----------------------------------------------------------------------------*/
error:
    mess_freelist_clear(&mem);
    return NULL;
}



