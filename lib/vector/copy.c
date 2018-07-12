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
 * @file lib/vector/copy.c
 * @brief Copy functions of @ref mess_vector structures.
 * @author @koehlerm
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


/*< one a variable to use with BLAS/Fortran code as a pointer */
static mess_int_t __ONE = 1;



/**
 * @brief Copy @ref mess_vector structure \f$ x \f$ to a @ref mess_vector \f$ y \f$ structure.
 * @param[in] x     input vector
 * @param[in,out] y input/output vector
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_vector_copy function copies a @ref mess_vector structure \f$ x \f$ to a @ref mess_vector structure \f$ y \f$. If
 * the size of \f$y\f$ is not correct it will be resized. The data type of \f$ y \f$ is changed to
 * the data type of \f$ x \f$. In case of errors the return
 * values will be:
 * \li @ref MESS_ERROR_NULLPOINTER -\f$ x\f$ or \f$y\f$ points to @c NULL
 * \li @ref MESS_ERROR_DATATYPE - \f$ x\f$ is not real or complex
 *
 */
int mess_vector_copy(mess_vector x, mess_vector y) {
    MSG_FNAME(__func__);
    mess_int_t dim;
    int ret = 0;

    mess_check_nullpointer(x);
    mess_check_nullpointer(y);

    if (x->dim != y->dim) {
        // MSG_WARN("dimension mismatch. resize\n");
        ret = mess_vector_resize(y, x->dim);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
    }

    if (x->data_type != y->data_type) {
        // MSG_WARN("data type mismatch in=%u, out=%u\n", x->data_type, y->data_type);
        // return MESS_ERROR_NOCOMPLEX;
        if (x->data_type == MESS_COMPLEX) {
            ret = mess_vector_tocomplex(y);         FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_tocomplex);
        } else {
            ret = mess_vector_toreal_nowarn(y);     FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_toreal_nowarn);
        }
    }
    dim = x->dim;
    if (MESS_IS_REAL(x)) {
        F77_GLOBAL(dcopy,DCOPY)(&dim, x->values, &__ONE, y->values, &__ONE);
    } else if (MESS_IS_COMPLEX(x)) {
        F77_GLOBAL(zcopy,ZCOPY)(&dim, x->values_cpx, &__ONE, y->values_cpx, &__ONE);
    } else {
        MSG_ERROR("unknown/unsupported data type\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
} /* -----  end of function mess_vector_copy  ----- */

/**
 * @brief Copy a @ref mess_vector strucutre into a complex @ref mess_vector structure.
 * @param[in] in        input vector
 * @param[in,out] out   input/output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_copy_tocomplex function copies a @ref mess_vector structure into a
 * complex @ref mess_vector structure, e.g. it does not change the data type if you copy a
 * real vector into a complex one like @ref mess_vector_copy it do.
 *
 */
int mess_vector_copy_tocomplex(mess_vector in, mess_vector out) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i;

    mess_check_nullpointer(in);
    mess_check_nullpointer(out);

    if ( MESS_IS_COMPLEX(in) && !MESS_IS_COMPLEX(out)) {
        ret = mess_vector_tocomplex(out);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_tocomplex);
    }
    if (in->dim != out->dim) {
        ret = mess_vector_resize(out, in->dim);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
    }
    if (MESS_IS_REAL(in)) {
        for (i = 0; i < in->dim; i++) {
            out->values_cpx[i] = in->values[i];
        }
    } else {
        for (i = 0; i < in->dim; i++) {
            out->values_cpx[i] = in->values_cpx[i];
        }
    }
    return 0;
} /* -----  end of function mess_vector_copy_tocomplex  ----- */


