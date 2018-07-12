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
 * @file python/c_interface/reimplemented.c
 * @brief Check if instance is reimplemented.
 *
 * @author @koehlerm
 *
 */
#include <Python.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "pymess.h"


/**
 * @brief Check if a function is reimplemented by a derived  object or equal to the one contained in the base class.
 * @param[in]   base        input PyObject containing the base class
 * @param[in]   instance    input Instance of a object
 * @param[in]   attr        input Name of the attribute to check
 * @return -1 in case of an error, 0 if the function only exists in the base class, or one if it is reimplemented.
 *
 * The @ref pymess_is_reimplemented function checks if a given instance reimplements the function named by "attr" or if
 * the one implemented in the base class is used.
 */
int pymess_is_reimplemented(PyObject * base, PyObject *instance, const char *attr ) {
    PyObject *py_base_func  = NULL;
    PyObject *py_inst_class = NULL;
    PyObject *py_inst_func  = NULL;
    int ret = 0;

    if ( base == NULL || instance == NULL || attr == NULL ) {
        return -1;
    }
    py_base_func = PyObject_GetAttrString(base, attr);
    if ( py_base_func == NULL ) {
        return -1;
    }

    py_inst_class = PyObject_GetAttrString(instance, "__class__");
    if ( py_inst_class == NULL ) {
        Py_DECREF(py_base_func);
        return -1;
    }

    py_inst_func = PyObject_GetAttrString(py_inst_class, attr);
    if ( py_inst_func == NULL ) {
        Py_DECREF(py_base_func);
        Py_DECREF(py_inst_class);
        return -1;
    }

    ret = PyObject_RichCompareBool(py_base_func, py_inst_func, Py_EQ);
    Py_DECREF(py_base_func);
    Py_DECREF(py_inst_class);
    Py_DECREF(py_inst_func);

    switch (ret) {
        case -1:
            return -1;
            break;
        case 0:
            return 1;
            break;
        case 1:
            return 0;
            break;
        default:
            ;
    }
    return -1;
}


