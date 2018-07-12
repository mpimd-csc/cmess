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
 * @file lib/vector/linspace.c
 * @brief Linspace functions of @ref mess_vector structures.
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
 * @brief Create a linspace vector.
 * @param[out] vect     output vector
 * @param[in] a         input   lower bound
 * @param[in] b         input   upper bound
 * @param[in] nsample   input number of samples
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_linspace function creates a linearly spaced vector.
 *
 */
int mess_vector_linspace(mess_vector vect, double a, double b,
        mess_int_t nsample) {
    MSG_FNAME(__func__);
    double h;
    mess_int_t i;
    int ret = 0;

    if(nsample <= 1 ){
        MSG_ERROR("nsample <= 1.");
        return MESS_ERROR_ARGUMENTS;
    }
    mess_check_nullpointer(vect);
    mess_check_positive(nsample);
    ret = mess_vector_toreal(vect);             FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_toreal);
    ret = mess_vector_resize(vect, nsample);    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);

    h = (b - a) / (double) (nsample - 1);
    for (i = 0; i < nsample; i++) {
        vect->values[i] = a + (double) i * h;
    }

    return 0;
} /* -----  end of function mess_vector_linspace  ----- */

/**
 * @brief Create a  \f$ \log_{10}  \f$ space vector.
 * @param[out] vect     output vector
 * @param[in] a         input lower bound for the exponential
 * @param[in] b         input upper bound for the exponential
 * @param[in] nsample   input number of samples
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_logspace10 function creates a logarithmically spaced vector to base \f$ 10  \f$.
 *
 */
int mess_vector_logspace10(mess_vector vect, double a, double b, mess_int_t nsample) {
    MSG_FNAME(__func__);
    double h;
    mess_int_t i;
    int ret = 0;

    if(nsample <= 1 ){
        MSG_ERROR("nsample <= 1.");
        return MESS_ERROR_ARGUMENTS;
    }

    mess_check_nullpointer(vect);
    mess_check_positive(nsample);
    ret = mess_vector_toreal(vect);             FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_toreal);
    ret = mess_vector_resize(vect, nsample);    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);

    h = (b - a) / (double) (nsample - 1);
    for (i = 0; i < nsample; i++) {
        vect->values[i] = pow(10, a + (double) i * h);
    }

    return 0;

} /* -----  end of function mess_vector_logspace10  ----- */

/**
 * @brief Create a  \f$ \log_{\mathrm{e}} \f$ space.
 * @param[out] vect     output vector
 * @param[in] a         input lower bound for the exponential
 * @param[in] b         input upper bound for the exponential
 * @param[in] nsample   input number of samples
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_logspacee function creates a logarithmically spaced vector to base \f$ \mathrm{e} \f$.
 *
 */
int mess_vector_logspacee(mess_vector vect, double a, double b,
        mess_int_t nsample) {
    MSG_FNAME(__func__);
    double h;
    mess_int_t i;
    int ret = 0;

    if(nsample <= 1 ){
        MSG_ERROR("nsample <= 1.");
        return MESS_ERROR_ARGUMENTS;
    }

    mess_check_nullpointer(vect);
    mess_check_positive(nsample);
    ret = mess_vector_toreal(vect);             FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_toreal);
    ret = mess_vector_resize(vect, nsample);    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);

    h = (b - a) / (double) (nsample - 1);
    for (i = 0; i < nsample; i++) {
        vect->values[i] = exp(a + (double) i * h);
    }

    return 0;

} /* -----  end of function mess_vector_logspacee  ----- */

/**
 * @brief Create a  \f$ \log_2 \f$ space.
 * @param[out] vect     output vector
 * @param[in] a         input lower bound for the exponential
 * @param[in] b         input upper bound for the exponential
 * @param[in] nsample   input number of samples
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_logspacee function creates a logarithmically spaced vector to base \f$ 2  \f$.
 *
 */
int mess_vector_logspace2(mess_vector vect, double a, double b,
        mess_int_t nsample) {
    MSG_FNAME(__func__);
    double h;
    mess_int_t i;
    int ret = 0;

    if(nsample <= 1 ){
        MSG_ERROR("nsample <= 1.");
        return MESS_ERROR_ARGUMENTS;
    }

    mess_check_nullpointer(vect);
    mess_check_positive(nsample);
    ret = mess_vector_toreal(vect);             FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_toreal);
    ret = mess_vector_resize(vect, nsample);    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);

    h = (b - a) / (double) (nsample - 1);
    for (i = 0; i < nsample; i++) {
        vect->values[i] = pow(2, a + (double) i * h);
    }

    return 0;

} /* -----  end of function mess_vector_logspace2  ----- */

