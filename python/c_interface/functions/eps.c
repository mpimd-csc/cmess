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
 * @file python/c_interface/functions/eps.c
 * @brief Interface to @ref mess_eps.
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
 * @brief  Return machine epsilon.
 *  The @ref pymess_eps is callable as @c eps and calls internally the @ref mess_eps function and returns the machine epsilon.
 * */
PyObject* pymess_eps(PyObject *self, PyObject* args){
    return Py_BuildValue("d",mess_eps());
}

