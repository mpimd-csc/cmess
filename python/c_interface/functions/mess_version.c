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
 * @file python/c_interface/functions/mess_version.c
 * @brief Interface to @ref mess_version.
 *
 * @author @mbehr
 *
 *
 */


#include <Python.h>
#include "mess/mess.h"
#include "../pymess.h"


/**
 * @brief Call the @ref mess_version function.
 *  The @ref pymess_mess_version is callable as @c mess_version and calls internally the @ref mess_version function.
 */
PyObject* pymess_mess_version(PyObject *self, PyObject* args){
    mess_version();
    Py_RETURN_NONE;
}


/**
 * @brief Call the @ref mess_version_verbose function.
 *  The @ref pymess_mess_version_verbose is callable as @c mess_version_verbose and calls internally the @ref mess_version_verbose function.
 */
PyObject* pymess_mess_version_verbose(PyObject *self, PyObject* args){
    mess_version_verbose();
    Py_RETURN_NONE;
}


/**
 * @brief Call the @ref mess_version_major function and returns the major version number.
 */
PyObject* pymess_mess_version_major(PyObject *self, PyObject* args){
    return Py_BuildValue("i",mess_version_major());
}


/**
 * @brief Call the @ref mess_version_minor function and returns the minor version number.
 */
PyObject* pymess_mess_version_minor(PyObject *self, PyObject* args){
    return Py_BuildValue("i",mess_version_minor());
}


/**
 * @brief Call the @ref mess_version_patch function and returns the patch version number.
 */
PyObject* pymess_mess_version_patch(PyObject *self, PyObject* args){
    return Py_BuildValue("i",mess_version_patch());
}

/**
 * @brief Call the @ref mess_git_id function and returns git id, which was used to configure the C-M.E.S.S. build.
 */
PyObject* pymess_mess_git_id(PyObject *self, PyObject* args){
    return PyUnicode_FromString(mess_git_id());
}

/**
 * @brief Call the @ref mess_git_branch function and returns git branch, which was used to configure the C-M.E.S.S. build.
 */
PyObject* pymess_mess_git_branch(PyObject *self, PyObject* args){
    return PyUnicode_FromString(mess_git_branch());
}


#undef PYMESS_MESS_HAVE_GENERATOR
#define PYMESS_MESS_HAVE_GENERATOR(NAME)                                \
    PyObject* pymess_## NAME (PyObject *self, PyObject* args){          \
        if( NAME ()){Py_RETURN_TRUE;}                                   \
        Py_RETURN_FALSE;                                                \
    }


PYMESS_MESS_HAVE_GENERATOR(mess_is_debug);
PYMESS_MESS_HAVE_GENERATOR(mess_have_zlib);
PYMESS_MESS_HAVE_GENERATOR(mess_have_bzip2);
PYMESS_MESS_HAVE_GENERATOR(mess_have_umfpack);
PYMESS_MESS_HAVE_GENERATOR(mess_have_amd);
PYMESS_MESS_HAVE_GENERATOR(mess_have_colamd);
PYMESS_MESS_HAVE_GENERATOR(mess_have_cholmod);
PYMESS_MESS_HAVE_GENERATOR(mess_have_csparse);
PYMESS_MESS_HAVE_GENERATOR(mess_have_superlu);
PYMESS_MESS_HAVE_GENERATOR(mess_have_mklpardiso);
PYMESS_MESS_HAVE_GENERATOR(mess_have_arpack);
PYMESS_MESS_HAVE_GENERATOR(mess_have_matio);
PYMESS_MESS_HAVE_GENERATOR(mess_have_openmp);
PYMESS_MESS_HAVE_GENERATOR(mess_have_mess64);


