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
 * @file python/c_interface/conv/status.c
 * @brief Convert a @ref mess_status structure to a python status object.
 *
 * @author @koehlerm
 * @author @nitin
 * @author @mbehr
 *
 */

#include <Python.h>
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "../pymess.h"
#include "cscutils/error_message.h"

/**
 * @brief Convert a @ref mess_status object from C to @python.
 * @param[in] status    input  @ref mess_status object
 * @return The @ref mess_status object convertet to an instance of pymess.status or Py_None.
 *
 * The @ref mess_status_to_python function converts a @ref mess_status object to an instance
 * of pymess.status.
 */
PyObject *mess_status_to_python(mess_status status) {
    if ( status == NULL ){
        return Py_None;
    }

    /*-----------------------------------------------------------------------------
     *  Create new Status object
     *-----------------------------------------------------------------------------*/
    PyObject *PyCMESS = PyImport_ImportModule((char*) "pymess");
    if ( PyCMESS == NULL ) {
        MESS_Py_XDECREF(PyCMESS);
        PYCMESS_ERROR(1, "Failed to import pymess\n");
    }
    PyObject *PyCMESS_status = PyObject_GetAttrString(PyCMESS,(char*)"Status");
    if ( PyCMESS_status == NULL ) {
        MESS_Py_XDECREF(PyCMESS, PyCMESS_status);
        PYCMESS_ERROR(1, "Failed to find the status object");
    }
    PyObject *return_status  = PyObject_CallObject(PyCMESS_status, NULL);
    if ( return_status == NULL ) {
        MESS_Py_XDECREF(PyCMESS, PyCMESS_status, return_status);
        PYCMESS_ERROR(1, "Failed to create an instance of pymess.Status\n");
    }

    /*-----------------------------------------------------------------------------
     *  fill values
     *-----------------------------------------------------------------------------*/
    PY_SET_VECTOR(return_status, "res2_norms", status->res2_norms);
    PY_SET_VECTOR(return_status, "rel_changes", status->rel_changes);
    PY_SET_DOUBLE(return_status, "res2_norm", status->res2_norm);
    PY_SET_LONG(  return_status, "it", status->it);
    PY_SET_DOUBLE(return_status, "res2_change", status->res2_change);
    PY_SET_DOUBLE(return_status, "res2_0", status->res2_0);
    PY_SET_DOUBLE(return_status, "rel_change", status->rel_change);
    PY_SET_BOOL(  return_status, "stop_res2", status->stop_res2);
    PY_SET_BOOL(  return_status, "stop_res2c", status->stop_res2c);
    PY_SET_BOOL(  return_status, "stop_rel", status->stop_rel);
    PY_SET_BOOL(  return_status, "stop_user", status->stop_user);
    PY_SET_DOUBLE(return_status, "time_all", status->time_all);
    PY_SET_DOUBLE(return_status, "time_adi", status->time_adi);
    PY_SET_BOOL(  return_status, "unstable", status->unstable);

    if (status->n_internal_status> 0 ) {
        int i;
        PY_SET_LONG( return_status, "n_internal_status", status->n_internal_status);
        PyObject *status_list = PyList_New(status->n_internal_status);
        for (i = 0; i < status->n_internal_status; i++ ) {
            PyObject * sub_status = mess_status_to_python(status->internal_status[i]);
            PyList_SetItem(status_list, i, sub_status);
        }
        PyObject_SetAttrString(return_status, "internal_status", status_list);
        MESS_Py_XDECREF(status_list);
    } else {
        PY_SET_LONG( return_status, "n_internal_status", 0);
        PyObject_SetAttrString(return_status, "internal_status", Py_None);
    }

    MESS_Py_XDECREF(PyCMESS, PyCMESS_status);
    return return_status;
}


