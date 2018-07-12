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
 * @file lib/vector/alloc.c
 * @brief Allocation function of @ref mess_vector structures.
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
 * @brief Allocate a @ref MESS_REAL or @ref MESS_COMPLEX the @ref mess_vector structure.
 * @param[in,out] vect  input/output allocated @ref mess_vector instance
 * @param[in] dim       input dimension of the new vector
 * @param[in] dtype     input data type of the new vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_alloc function allocates  a vector @p vect  with length @p dim and
 * @ref mess_datatype_t @p dtype. \n
 * If no memory is left it returns @ref MESS_ERROR_MEMORY.
 *
 */
int mess_vector_alloc(mess_vector vect, mess_int_t dim, mess_datatype_t dtype) {
    MSG_FNAME(__func__);

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(vect);

    /*-----------------------------------------------------------------------------
     *  alloc fields
     *-----------------------------------------------------------------------------*/
    vect->dim = dim;
    if (dtype == MESS_REAL) {
        mess_try_alloc(vect->values, double *, sizeof(double) * dim);
        memset(vect->values, 0, sizeof(double)*dim);
        vect->values_cpx = NULL;
        vect->data_type = MESS_REAL;
    } else if (dtype == MESS_COMPLEX) {
        mess_try_alloc(vect->values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * dim);
        memset(vect->values_cpx, 0, sizeof(mess_double_cpx_t)*dim);
        vect->values = NULL;
        vect->data_type = MESS_COMPLEX;
    } else {
        MSG_ERROR("unknown data type\n");
        return MESS_ERROR_DATATYPE;
    }

    return 0;
} /* -----  end of function mess_vector_alloc  ----- */



/**
 * @brief Clean up the @ref mess_vector structure to use it again.
 * @param[in,out] vect input/outout pointer to the @ref mess_vector structure
 * @return always zero
 *
 * The @ref mess_vector_clear function removes a vector from memory and cleans up the @ref mess_vector structure.
 */
int mess_vector_clear(mess_vector *vect) {
    // MSG_FNAME(__func__);
    if (vect == NULL)
        return 0;
    if (*vect == NULL) {
        // MSG_ERROR("vect points to NULL\n");
        return 0;
    }
    if ((*vect)->values != NULL)
        mess_free((*vect)->values);
    if ((*vect)->values_cpx != NULL)
        mess_free((*vect)->values_cpx);

    mess_free(*vect);
    *vect = NULL;
    return 0;
} /* -----  end of function mess_vector_clear  ----- */




/**
 *
 * @brief Resets a @ref mess_vector object.
 * @param[in,out] v  input/output vector object to reset
 * @return zeros on success or a non zero error value.
 *
 * The @ref mess_vector_reset function resets a \ref mess_vector object that it behaves like a newly initialized one.
 * In contrast to @ref mess_vector_clear it does only free the internal data and not the surrounding structure.
 */
int mess_vector_reset(mess_vector v){
    MSG_FNAME(__func__);
    mess_check_nullpointer(v);

    mess_free(v->values);
    mess_free(v->values_cpx);
    v->data_type = MESS_REAL;
    v->dim = 0;
    return 0;
}/* -----  end of function mess_vector_reset  ----- */





