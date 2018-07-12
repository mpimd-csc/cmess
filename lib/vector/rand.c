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
 * @file lib/vector/rand.c
 * @brief Rand functions of @ref mess_vector structures.
 * @author @koehlerm
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
 * @brief Generate a positive random floating point number between -1 and 1.
 * @return A positive random number between -1 and 1
 *
 * The @ref __drand function generates an uniformly distributed random number between
 * -1 and 1.
 *
 */
static double __drand()
{
    double ret = -1.0 + 2.0*( (double) rand() / ( (double) RAND_MAX ) );
    return ret;
}

/**
 * @brief Init random number generator with given seed.
 * @param[in] seed   pointer to @ref mess_int_t for given seed, if @c NULL current calender time is used.
 * @return zero in every case
 *
 * The @ref mess_vector_rand_init initializes the @c rand  function from @c libc for generating random numbers. \n
 * If @p seed points to @c NULL the current calender time is used for initialization of @c rand function.
 * This function is usefull in combination with:
 * This is a wrapper around @ref mess_matrix_rand_init.
 *
 * @sa mess_vector_rand
 *
 */
int mess_vector_rand_init(mess_int_t *seed){ return mess_matrix_rand_init(seed);}

/**
 * @brief Fill a vector with random numbers.
 * @param[in,out] v input/ouput vector \f$ v\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_rand function fills a vector \f$v\f$ with random entries.
 */
int mess_vector_rand(mess_vector v) {
    MSG_FNAME(__func__);
    mess_int_t i;
    // int ret = 0;
    mess_check_nullpointer(v);
    if (MESS_IS_REAL(v)) {
        for (i = 0; i < v->dim; i++)
            v->values[i] = __drand();
    } else if (MESS_IS_COMPLEX(v)) {
        for (i = 0; i < v->dim; i++)
            v->values_cpx[i] = __drand() + I * __drand();
    } else {
        MSG_ERROR("unknown datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}



