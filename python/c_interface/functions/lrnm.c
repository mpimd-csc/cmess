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
 * @file python/c_interface/functions/lrnm.c
 * @brief Interface to @ref mess_lrnm.
 *
 * @author @koehlerm
 * @author @nitin
 * @author @mbehr
 *
 */


#include <Python.h>
#include "mess/mess.h"
#include "../pymess.h"


PyObject* pymess_lrnm(PyObject * self, PyObject * args) {
    int need_callback = 1;
    int ret = 0;
    PyObject *pyEqn     = NULL;
    PyObject *pyOpts    = NULL;
    PyObject *pyZ       = NULL;
    PyObject *pyStatus  = NULL;
    PyObject *pyResult  = NULL;
    mess_equation eqn       = NULL;
    mess_options  opt       = NULL;
    mess_status   status    = NULL;
    mess_matrix   Z         = NULL;

    mess_freelist mem   = NULL;;
    mess_freelist_init(&mem);


    /*-----------------------------------------------------------------------------
     *  parse input arguments
     *-----------------------------------------------------------------------------*/
    if (!PyArg_ParseTuple(args, "OO", &pyEqn, &pyOpts)){
        PyErr_SetString(PyExc_RuntimeError, "Cannot parse input arguments.");
        goto error;
    }


    /*-----------------------------------------------------------------------------
     *  get equation from python
     *-----------------------------------------------------------------------------*/
    eqn = mess_equation_from_python(pyEqn, mem, MESS_EQN_GRICCATI, &need_callback);
    if ( eqn == NULL ) {
        PyErr_SetString(PyExc_Exception, "Cannot setup  Equation.");
        goto error;
    }
    //equation was added by mess_equation_from_python
    //mess_freelist_add_mess_equation(mem,eqn);


    /*-----------------------------------------------------------------------------
     * get options from python
     *-----------------------------------------------------------------------------*/
    opt = mess_options_from_python (pyOpts);
    if(!opt){
        PyErr_SetString(PyExc_Exception, "Cannot copy the options from Python.");
        goto error;
    }
    mess_freelist_add_mess_options(mem,opt);

    /*-----------------------------------------------------------------------------
     *  call lrnm
     *-----------------------------------------------------------------------------*/
    mess_matrix_init(&Z);
    mess_freelist_add_mess_matrix(mem,Z);
    mess_status_init(&status);
    mess_freelist_add_mess_status(mem,status);
    if ( need_callback ) {
        ret =  mess_lrnm(eqn, opt, status, Z);
    } else {
        Py_BEGIN_ALLOW_THREADS;
        ret = mess_lrnm(eqn, opt, status, Z);
        Py_END_ALLOW_THREADS;
    }

    if(ret){
        PyErr_SetString(PyExc_RuntimeError,"C-M.E.S.S.: mess_lrnm returned an error.");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  transfer results to Python
     *-----------------------------------------------------------------------------*/
    pyZ = matrix_to_python(Z);
    if (!pyZ) {
        PyErr_SetString(PyExc_RuntimeError, "Cannot transfer matrix a from C to Python.");
        goto error;
    }

    pyStatus = mess_status_to_python(status);
    pyResult = Py_BuildValue("OO", pyZ, pyStatus);
    if(!pyResult){
        Py_XDECREF(pyZ);
        Py_XDECREF(pyStatus);
        PyErr_SetString(PyExc_RuntimeError, "Cannot build result");
        goto error;
    }

    /*-----------------------------------------------------------------------------
     *  clear data and return
     *-----------------------------------------------------------------------------*/
    Py_DECREF(pyZ);
    Py_DECREF(pyStatus);

    mess_freelist_clear(&mem);

    return pyResult;

    /*-----------------------------------------------------------------------------
     *  error handling
     *-----------------------------------------------------------------------------*/
error:
    mess_freelist_clear(&mem);
    return NULL;
}


