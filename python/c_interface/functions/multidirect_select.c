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
 * @file python/c_interface/functions/multidirect_select.c
 * @brief Interface to @ref mess_multidirect_select.
 *
 * @author @mbehr
 *
 */

#include <Python.h>
#include "mess/mess.h"
#include "../pymess.h"


/**
 * @brief Wrapper around \ref mess_multidirect_select for use with @python.
 * @param[in]  self input
 * @param[in]  args input
 * @return @python None.
 *
 * The @ref pymess_multidirect_select function provides and wrapper arount the \ref mess_multidirect_select function  in order
 * to use from @python. The function is called as
 * \code[py]
 * Z = pymess.multidirect_select( MESS_MULTIDIRECT_UMFPACK_LU )
 * \endcode
 * to use @umfpack as a multi direct solver.
 *
 *
 */
PyObject* pymess_multidirect_select(PyObject*self, PyObject *args){
    int ret = 0;
    mess_multidirect_t mdirect;

    /*-----------------------------------------------------------------------------
     *  parse input
     *-----------------------------------------------------------------------------*/
    if(!PyArg_ParseTuple(args, "i", &mdirect)){
        PYCMESS_ERROR(1,"Cannot convert to mess_multidirect_t.");
    }

    /*-----------------------------------------------------------------------------
     *  call lu select
     *-----------------------------------------------------------------------------*/
    ret = mess_multidirect_select(mdirect); PYCMESS_ERROR(ret, "mess_multidirect_select returned an error.");

    /*-----------------------------------------------------------------------------
     *  return None
     *-----------------------------------------------------------------------------*/
    Py_RETURN_NONE;

}



