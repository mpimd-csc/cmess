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
 * @file lib/vector/ones.c
 * @brief Fill functions of @ref mess_vector structures.
 * @author @koehlerm
 *
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
 * @brief Fill a vector with ones.
 * @param[in,out] v input/output vector \f$v\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_ones function fills a vector \f$v\f$ with one entries.
 */
int mess_vector_ones(mess_vector v) {
    MSG_FNAME(__func__);
    mess_int_t i;
    // int ret = 0;

    mess_check_nullpointer(v);
    if (MESS_IS_REAL(v)) {
        for (i = 0; i < v->dim; i++)
            v->values[i] = 1.0;
    } else if (MESS_IS_COMPLEX(v)) {
        for (i = 0; i < v->dim; i++)
            v->values_cpx[i] = 1.0;
    } else {
        MSG_ERROR("unknown datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}




/**
 * @brief Fill a vector with zeros.
 * @param[in,out] v input/output vector \f$v\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_zeros function fills a vector \f$v\f$ with zeros.
 */
int mess_vector_zeros(mess_vector v) {
    MSG_FNAME(__func__);
    mess_int_t i;
    // int ret = 0;

    mess_check_nullpointer(v);
    if (MESS_IS_REAL(v)) {
        for (i = 0; i < v->dim; i++)
            v->values[i] = 0.0;
    } else if (MESS_IS_COMPLEX(v)) {
        for (i = 0; i < v->dim; i++)
            v->values_cpx[i] = 0.0;
    } else {
        MSG_ERROR("unknown datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}

