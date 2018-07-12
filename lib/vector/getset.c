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
 * @file lib/vector/getset.c
 * @brief Get/set of @ref mess_vector structures.
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
 * @brief Get entry @c i of a @ref mess_vector instance.
 * @param[in] v         input @ref mess_vector instance
 * @param[in] i         input index of entry
 * @param[out] val      output the @p i -th entry of @p v.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_get function gets the @p i -th entry of @p v.
 * If the index @p i is out of range or @p v points to @c NULL an error is returned.
 * We assume zero based indexing.
 *
 *
 */
int mess_vector_get(mess_vector v, mess_int_t i, mess_double_cpx_t* val) {
    MSG_FNAME(__func__);

    mess_check_nullpointer(v);
    mess_check_nullpointer(val);
    mess_check_real_or_complex(v);

    /*-----------------------------------------------------------------------------
     *  check index i
     *-----------------------------------------------------------------------------*/
    if( i < 0  || v->dim <= i){
        return MESS_ERROR_ARGUMENTS;
    }


    if(MESS_IS_REAL(v)){
        *val = v->values[i];
    }else{
        *val = v->values_cpx[i];
    }

    return 0;
} /* -----  end of function mess_vector_get  ----- */


