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
 * @file python/c_interface/functions/mess_set_errorlevel.c
 * @brief Interface to @ref mess_set_errorlevel.
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
 * @brief Wrapper around \ref mess_set_errorlevel for use with @python.
 * @param[in]  self input
 * @param[in]  args input
 * @return @python None.
 *
 * The @ref pymess_set_errorlevel function provides and wrapper arount the \ref mess_set_errorlevel function  in order
 * to use from @python. The function is called as
 * \code[py]
 * pymess.set_errorlevel(0)
 * pymess.set_errorlevel(1)
 * pymess.set_errorlevel(2)
 * pymess.set_errorlevel(3)
 * \endcode
 *
 *
 *
 */
PyObject* pymess_set_errorlevel(PyObject*self, PyObject *args){
    int lvl;

    /*-----------------------------------------------------------------------------
     *  parse input
     *-----------------------------------------------------------------------------*/
    if(!PyArg_ParseTuple(args, "i", &lvl)){
        PYCMESS_ERROR(1,"Cannot get argument.");
    }

    /*-----------------------------------------------------------------------------
     *  call mess_set_error_level
     *-----------------------------------------------------------------------------*/
    mess_set_errorlevel(lvl);

    /*-----------------------------------------------------------------------------
     *  return None
     *-----------------------------------------------------------------------------*/
    Py_RETURN_NONE;
}

