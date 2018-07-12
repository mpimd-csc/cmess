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
 * @file lib/vector/axpy.c
 * @brief axpy function of @ref mess_vector structures.
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
 * @internal
 * @brief Internal axpy for real/complex mixing.
 * @param[in] a     input mess_double_cpx_t calar
 * @param[in] x     input   first vector
 * @param[in,out] y input/output    second vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __axpy_rc function computes
 * \f[ y \leftarrow ax + y \f]
 * with \f[ x \f] real and  \f[ y \f] complex.
 *
 * @attention Internal use only.
 */
static int __axpy_rc ( mess_double_cpx_t a, mess_vector x, mess_vector y )
{
    mess_int_t i;
#ifdef _OPENMP
#pragma omp parallel for private(i) default(shared)
#endif
    for (i = 0; i < x->dim; i++) {
        y->values_cpx[i]+=a*x->values[i];
    }

    return 0;
} /* -----  end of function __axpy_rc  ----- */

/**
 * @brief Compute  \f$ y \leftarrow ax + y \f$ (real).
 * @param[in] a     input scaling scalar \f$a\f$
 * @param[in] x     input vector \f$x\f$
 * @param[in,out] y input/output vector \f$y\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_axpy function adds \f$a\f$ times a vector  \f$ x  \f$  to  \f$ y  \f$ and store it in  \f$ y  \f$
 * \f[y \leftarrow a\cdot x+y . \f]
 *
 * @sa mess_vector_axpyc
 */
int mess_vector_axpy(double a, mess_vector x, mess_vector y) {
    MSG_FNAME(__func__);
    mess_int_t dim;
    // float alpha  =0;

    mess_check_nullpointer(x);
    mess_check_nullpointer(y);

    if (x->dim != y->dim) {
        MSG_ERROR(
                "dimension mismatch x->dim=" MESS_PRINTF_INT " , y->dim=" MESS_PRINTF_INT "\n",
                x->dim, y->dim);
        return MESS_ERROR_DIMENSION;
    }
    /* if (x->data_type != y->data_type){
       MSG_ERROR("data type mismatch, x->dt=%s , y->dt=%s\n", mess_datatype_t_str(x->data_type), mess_datatype_t_str(y->data_type));
       return MESS_ERROR_DATATYPE;
       }
       */
    if ( MESS_IS_REAL(x) && MESS_IS_COMPLEX(y)) {
        return __axpy_rc((mess_double_cpx_t)a, x, y);
    }
    if ( MESS_IS_REAL(y) && MESS_IS_COMPLEX(x)) {
        mess_vector_tocomplex(y);
        return mess_vector_axpyc((mess_double_cpx_t) a, x, y);
    }
    if (MESS_IS_COMPLEX(x)) {
        return mess_vector_axpyc((mess_double_cpx_t) a, x, y);
    }

    dim = x->dim;
    if (!MESS_IS_REAL(x)) {
        MSG_ERROR("unknown/unsupported data type\n");
        return MESS_ERROR_DATATYPE;

    }
    F77_GLOBAL(daxpy,DAXPY)(&dim, &a, x->values, &__ONE, y->values, &__ONE);
    return 0;
} /* -----  end of function mess_vector_axpy  ----- */

/**
 * @brief Compute   \f$ y \leftarrow ax + y  \f$ (complex).
 * @param[in]  a    input scalar (double or mess_double_cpx_t) \f$a\f$
 * @param[in] x     input vector
 * @param[in,out] y input/output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_axpyc function adds \f$a\f$ times a vector  \f$ x  \f$ to  \f$ y  \f$ and store it in  \f$ y \f$:
 * \f[y \leftarrow a\cdot x+y . \f]
 *
 * @sa mess_vector_axpy
 *
 */
int mess_vector_axpyc(mess_double_cpx_t a, mess_vector x, mess_vector y) {
    MSG_FNAME(__func__);
    mess_int_t dim;

    mess_check_nullpointer(x);mess_check_nullpointer(y);

    if (x->dim != y->dim){
        MSG_ERROR("dimension mismatch\n");return MESS_ERROR_DIMENSION;
    }
    if ( MESS_IS_REAL(x)&& MESS_IS_REAL(y)&& cimag(a) == 0.0) {
        return mess_vector_axpy(creal(a), x,y);
    }
    if ( MESS_IS_REAL(y)){
        mess_vector_tocomplex(y);
    }

    if ( MESS_IS_REAL(x)) {
        return __axpy_rc(a, x, y);
    }

    dim = x->dim;
    F77_GLOBAL(zaxpy,ZAXPY)(&dim, &a, x->values_cpx, &__ONE, y->values_cpx, &__ONE);
    return 0;
} /* -----  end of function mess_vector_axpyc  ----- */

