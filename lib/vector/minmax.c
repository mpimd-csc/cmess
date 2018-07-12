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
 * @file lib/vector/minmax.c
 * @brief MinMax functions of @ref mess_vector structures.
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
 * @brief Return the largest entry of a vector.
 * @param[in] v         input vector
 * @param[out] maxval   output the largest absolute value
 * @param[out] maxind   output index of the maxval element
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_max function returns the largest absolute value and its
 * position in the vector.
 *
 */
int mess_vector_max(mess_vector v, double * maxval, mess_int_t *maxind) {
    MSG_FNAME(__func__);
    double mv;
    mess_int_t mi;
    mess_int_t i;

    mess_check_nullpointer(v);
    mess_check_nullpointer(maxval);
    mess_check_nullpointer(maxind);
    if (v->dim < 1) {
        MSG_ERROR("The dimension of the vector must at least be one.");
        return MESS_ERROR_DIMENSION;
    }
    if (MESS_IS_REAL(v)) {
        mv = fabs(v->values[0]);
        mi = 0;
        for (i = 1; i < v->dim; i++) {
            if (mv < fabs(v->values[i])) {
                mv = fabs(v->values[i]);
                mi = i;
            }
        }
    } else if (MESS_IS_COMPLEX(v)) {
        mv = cabs(v->values_cpx[0]);
        mi = 0;
        for (i = 1; i < v->dim; i++) {
            if (mv < cabs(v->values_cpx[i])) {
                mv = cabs(v->values_cpx[i]);
                mi = i;
            }
        }
    } else {
        MSG_ERROR("unknown datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    *maxval = mv;
    *maxind = mi;
    return 0;
} /* -----  end of function mess_vector_max  ----- */

/**
 * @brief Get the minimal value in a vector.
 * @param[in] v         input vector
 * @param[out] min      output minimum value
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_minvalue function gets the minimum value of a vector. \n
 * In case of a complex vector we only use the real part.
 *
 */
int mess_vector_minvalue(mess_vector v, double *min) {
    MSG_FNAME(__func__);
    mess_int_t i;
    double m;

    mess_check_nullpointer(v);
    mess_check_real_or_complex(v);
    if (v->dim == 0) {
        *min = 0;
        return 0;
    }
    if (v->dim == 1) {
        if (MESS_IS_REAL(v))
            *min = v->values[0];
        else
            *min = creal(v->values_cpx[0]);
        return 0;
    }

    if (MESS_IS_REAL(v)) {
        m = v->values[0];
        for (i = 0; i < v->dim; i++) {
            m = (m < v->values[i]) ? m : v->values[i];
        }
    } else {
        m = v->values_cpx[0];
        for (i = 0; i < v->dim; i++) {
            m = (m < creal(v->values_cpx[i])) ? m : creal(v->values_cpx[i]);
        }
    }
    *min = m;
    return 0;
} /* -----  end of function mess_vector_min  ----- */

/**
 * @brief Get the maximum value in a vector.
 * @param[in] v     input vector
 * @param[out] max  output maximum value
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_max function gets the maximum value of a vector. \n
 * In case of a complex vector we only use the real part.
 *
 */
int mess_vector_maxvalue(mess_vector v, double *max) {
    MSG_FNAME(__func__);
    mess_int_t i;
    double m;

    mess_check_nullpointer(v);
    mess_check_real_or_complex(v);
    if (v->dim == 0) {
        *max = 0;
        return 0;
    }
    if (v->dim == 1) {
        if (MESS_IS_REAL(v))
            *max = v->values[0];
        else
            *max = creal(v->values_cpx[0]);
        return 0;
    }

    if (MESS_IS_REAL(v)) {
        m = v->values[0];
        for (i = 0; i < v->dim; i++) {
            m = (m > v->values[i]) ? m : v->values[i];
        }
    } else {
        m = v->values_cpx[0];
        for (i = 0; i < v->dim; i++) {
            m = (m > creal(v->values_cpx[i])) ? m : creal(v->values_cpx[i]);
        }
    }
    *max = m;
    return 0;
} /* -----  end of function mess_vector_max  ----- */

