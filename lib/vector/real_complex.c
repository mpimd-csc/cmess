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
 * @file lib/vector/real_complex.c
 * @brief Realpart Imagpart functions of @ref mess_vector structures.
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





/**
 * @brief Get a real value from a vector.
 * @param[in] vec input vector
 * @param[in] index input index of the value
 * @return pointer value at the index or a non-zero error value
 *
 * The @ref mess_vector_valuereal function gets a pointer to value of the given index.
 * If the vector is not long enough it is resized.
 *
 */
double *mess_vector_valuereal(mess_vector vec, mess_int_t index) {
    MSG_FNAME(__func__);
    static double null = 0.0;
    int err = 0;
    if (vec == NULL) {
        MSG_ERROR("input vector points to NULL\n");
        return &(null);
    }
    if (!MESS_IS_REAL(vec)) {
        MSG_ERROR("input vector isn't real\n");
        return &(null);
    }
    if (index < 0) {
        MSG_ERROR("index < 0 \n");
        return &(null);
    }
    if (index >= vec->dim) {
        if ((err = mess_vector_resize(vec, index)) != 0) {
            MSG_ERROR("mess_vector_resize returned with %d - %s\n", err, mess_get_error(err));
            return &(null);
        }
    }
    return &(vec->values[index]);
}

/**
 * @brief Get a complex value from a vector.
 * @param[in] vec input vector
 * @param[in] index input index of the value
 * @return pointer value at the index or a non-zero error value
 *
 * The @ref mess_vector_valuecomplex function gets a pointer to value of the given index.
 * If the vector is not long enough it is resized.
 *
 */
mess_double_cpx_t *mess_vector_valuecomplex(mess_vector vec, mess_int_t index){
    MSG_FNAME(__func__);static mess_double_cpx_t null = 0.0;
    int err = 0;
    if ( vec == NULL ){
        MSG_ERROR("input vector points to NULL\n");return &(null);
    }
    if ( !MESS_IS_COMPLEX(
                vec)){
        MSG_ERROR("input vector isn't real\n");return &(null);
    }
    if ( index < 0 ) {
        MSG_ERROR("index < 0 \n");return &(null);
    }
    if ( index >= vec->dim) {
        if ( (err=mess_vector_resize(vec, index)) != 0) {
            MSG_ERROR(
                    "mess_vector_resize returned with %d - %s\n", err, mess_get_error(err)); return &(null);
        }
    }
    return &(vec->values_cpx[index]);
}




/**
 * @brief Get the real part of a vector.
 * @param[in] in        input vector
 * @param[in,out] out   output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_realpart function gets the real part of a vector.
 *
 */
int mess_vector_realpart(mess_vector in, mess_vector out) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);

    ret = mess_vector_toreal(out);
    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_toreal);

    if (MESS_IS_REAL(in)) {
        ret = mess_vector_copy(in, out);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_copy);
        return 0;
    } else if (MESS_IS_COMPLEX(in)) {
        ret = mess_vector_resize(out, in->dim);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
        for (i = 0; i < in->dim; i++) {
            out->values[i] = creal(in->values_cpx[i]);
        }
        return 0;
    } else {
        MSG_ERROR("unknown datatype.\n");
        return MESS_ERROR_DATATYPE;
    }

    return 0;
} /* -----  end of function mess_vector_realpart  ----- */

/**
 * @brief Get the imaginary part of a vector.
 * @param[in] in        input vector
 * @param[in,out] out   output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_imagpart function gets the imaginary part of a vector.
 *
 */
int mess_vector_imagpart(mess_vector in, mess_vector out) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);

    ret = mess_vector_toreal(out);
    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_toreal);
    ret = mess_vector_resize(out, in->dim);
    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);

    if (MESS_IS_REAL(in)) {
        ret = mess_vector_zeros(out);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_zeros);
        return 0;
    } else if (MESS_IS_COMPLEX(in)) {
        for (i = 0; i < in->dim; i++) {
            out->values[i] = cimag(in->values_cpx[i]);
        }
        return 0;
    } else {
        MSG_ERROR("unknown datatype.\n");
        return MESS_ERROR_DATATYPE;
    }

    return 0;
} /* -----  end of function mess_vector_realpart  ----- */

/**
 * @brief Generate a complex vector from real and imaginary part.
 * @param[in] xr        input real part
 * @param[in] xc        input imaginary part
 * @param[out] x        output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_complex_from_parts function creates a complex vector from real and imaginary part.
 *
 */
int mess_vector_complex_from_parts(mess_vector xr, mess_vector xc, mess_vector x) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(xr);
    mess_check_nullpointer(xc);
    mess_check_nullpointer(x);
    mess_check_real(xr);
    mess_check_real(xc);
    if (xr->dim != xc->dim) {
        MSG_ERROR("xr and xc must have the same dimension. dim(xr) = " MESS_PRINTF_INT " \t dim(xc) = " MESS_PRINTF_INT "\n",
                    xr->dim, xc->dim);
        return MESS_ERROR_DIMENSION;
    }
    ret = mess_vector_tocomplex(x);
    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_tocomplex);
    ret = mess_vector_resize(x, xr->dim);
    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
    for (i = 0; i < xr->dim; i++) {
        x->values_cpx[i] = xr->values[i] + I * xc->values[i];
    }
    return 0;
} /* -----  end of function mess_vector_complex_from_parts  ----- */


/**
 * @brief Check if a vector is real or if the imaginary part is  \f$ 0  \f$.
 * @param[in] v input vector
 * @return one if v is real, otherwise zero
 *
 * The @ref mess_vector_isreal function checks if a vector is real or if the imaginary part is \f$ 0 \f$.
 *
 */
int mess_vector_isreal(mess_vector v) {
    // MSG_FNAME(__func__);
    mess_int_t i;
    mess_int_t c = 0;
    double eps = mess_eps();
    if (v == NULL)
        return 0;
    if (MESS_IS_REAL(v))
        return 1;
    for (i = 0; i < v->dim; i++) {
        if (fabs(cimag(v->values_cpx[i])) > eps) {
            c++;
        }
    }
    if (c == 0)
        return 1;
    else
        return 0;
    return 0;
} /* -----  end of function mess_vector_isreal  ----- */



/**
 * @brief Complex conjugatation of a vector.
 * @param[in,out] vector input/output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_conj function performs a complex conjugation on
 * a vector.
 *
 */
int mess_vector_conj(mess_vector vector) {
    MSG_FNAME(__func__);
    mess_int_t i;

    mess_check_nullpointer(vector);
    mess_check_real_or_complex(vector);
    if (MESS_IS_REAL(vector))
        return 0;

    for (i = 0; i < vector->dim; i++) {
        vector->values_cpx[i] = conj(vector->values_cpx[i]);
    }

    return 0;
} /* -----  end of function mess_vector_conj  ----- */




