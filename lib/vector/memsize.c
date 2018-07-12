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
 * @file lib/vector/memsize.c
 * @brief Memory size function of @ref mess_vector structures.
 * @author @mbehr
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>
#include "blas_defs.h"
#ifdef _OPENMP_H
#include <omp.h>
#endif

/**
 * On MacOS X 10.6 the zdotc implementation
 * is defect, this is the work around
 */
#ifdef MESS_USE_APPLE_BLAS
#include <Accelerate/Accelerate.h>
#endif





/**
 * @brief Return the size of a @ref mess_vector in bytes.
 * @param[in] v input vector
 * @return approximate size of the @p v in bytes.
 *
 * The @ref mess_vector_memsize function determines the size of @p v in bytes. \n
 *
 * @sa mess_matrix_memsize
 *
 */
mess_int_t mess_vector_memsize(mess_vector v){
    MSG_FNAME(__func__);

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(v);


    if ( MESS_IS_COMPLEX(v)) {
        return v->dim*sizeof(mess_double_cpx_t);
    } else if ( MESS_IS_REAL(v)){
        return v->dim* sizeof(double);
    } else {
        MSG_ERROR("Unknown datatype.");
        return MESS_ERROR_DATATYPE;
    }
}

