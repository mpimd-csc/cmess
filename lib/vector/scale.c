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
 * @file lib/vector/scale.c
 * @brief Scale functions of @ref mess_vector structures.
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
 * @brief Scale a vector with a real scalar.
 * @param[in]  a        input scalar \f$a\f$
 * @param[in,out]  x    input/output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_scale function scales a vector \f$ x \f$ with a scalar \f$ a \f$ such that
 * \f[ x \leftarrow a\cdot x. \f]
 *
 *
 */
int mess_vector_scale(double a, mess_vector x) {
    MSG_FNAME(__func__);
    mess_int_t dim;
    // float alpha = 0;
    mess_check_nullpointer(x);

    if (MESS_IS_COMPLEX(x)) {
        return mess_vector_scalec((mess_double_cpx_t)a, x);
    }
    dim = x->dim;
    if (MESS_IS_REAL(x)) {
        F77_GLOBAL(dscal,DSCAL)(&dim, &a, x->values, &__ONE);
    } else {
        MSG_ERROR("unknown/unsupported data type\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}

/**
 * @brief Scale a vector with a complex scalar.
 * @param[in] a input scalar \f$a\f$
 * @param[in] x input/output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_scalec function scales a vector \f$ x \f$ with a complex scalar \f$ a \f$ such that
 * \f[ x\leftarrow a\cdot x. \f]
 *
 */
int mess_vector_scalec(mess_double_cpx_t a, mess_vector x) {
    MSG_FNAME(__func__);
    mess_int_t dim;
    int ret = 0;

    mess_check_nullpointer(x);
    mess_check_real_or_complex(x);

    if(MESS_IS_REAL(x)){
        if(!cimag(a)){
            return mess_vector_scale(creal(a), x);
        }else{
            ret = mess_vector_tocomplex(x);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        }
    }

    dim = x->dim;
    F77_GLOBAL(zscal,ZSCAL)(&dim, &a, x->values_cpx, &__ONE);
    return 0;
}

