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
 * @file lib/vector/cat.c
 * @brief Cat and Resize functions of @ref mess_vector structures.
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
 * @brief Resize a vector.
 * @param[in,out] v input/outout vector to resize
 * @param[in] dim   input new dimension of the vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_resize function resizes a given vector  \f$ v  \f$ to dimension dim. \n
 * If the new dimension is smaller the  \f$ v \to dim  \f$ elements after \f$ dim \f$ are lost. \n
 * The following errors can be returned:
 * <center>
 *  | Error                         |   Description                                                 |
 *  |:-----------------------------:|:-------------------------------------------------------------:|
 *  | @ref MESS_ERROR_DIMENSION     | the new dimension is less or equal to zero                    |
 *  | @ref MESS_ERROR_NULLPOINTER   | the input points to @c NULL                                   |
 *  | @ref MESS_ERROR_MEMORY        | reallocating memory caused an error                           |
 *  | @ref MESS_ERROR_NOCOMPLEX     | no support for complex numbers available                      |
 *  | @ref MESS_ERROR_DATATYPE      | the data type in \f$v\f$ is wrong                             |
 *  </center>
 */
int mess_vector_resize(mess_vector v, mess_int_t dim) {
    MSG_FNAME(__func__);
    mess_int_t i;

    if (dim < 0) {
        MSG_ERROR("wrong dimension. dim = " MESS_PRINTF_INT " \n", dim);
        return MESS_ERROR_DIMENSION;
    }
    mess_check_nullpointer(v);
    if (v->dim == dim)
        return 0;

    if (MESS_IS_REAL(v)) {
        mess_try_realloc(v->values, double *, sizeof(double) * dim);
        if (v->dim < dim) {
            for (i = v->dim; i < dim; i++)
                v->values[i] = 0;
        }
    } else if (MESS_IS_COMPLEX(v)) {
        mess_try_realloc(v->values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t )* dim);
        for (i = v->dim; i < dim; i++)
            v->values_cpx[i] = 0;
    } else {
        MSG_ERROR("unknown/unsupported datatype = %u \n", v->data_type);
        return MESS_ERROR_DATATYPE;
    }
    v->dim = dim;
    return 0;
}



/**
 * @brief Concatenate two vectors.
 * @param[in] x1    input vector \f$x_1\f$
 * @param[in] x2    input vector \f$x_2\f$
 * @param[out] x    output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_cat function concatenates two vectors, i.e.
 * \f[ x= \left[
 *  \begin{array}{cc}
 * x_1 & x_2
 * \end{array}
 *  \right].
 * \f]
 *
 */
int mess_vector_cat(mess_vector x1, mess_vector x2, mess_vector x) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(x1);
    mess_check_nullpointer(x2);
    mess_check_nullpointer(x);

    if ( MESS_IS_COMPLEX(x1) || MESS_IS_COMPLEX(x2)) {
        ret = mess_vector_tocomplex(x);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_tocomplex);
    }
    ret = mess_vector_resize(x, x1->dim + x2->dim);
    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);

    if (MESS_IS_COMPLEX(x)) {
        if (MESS_IS_REAL(x1)) {
            for (i = 0; i < x1->dim; i++) {
                x->values_cpx[i] = x1->values[i];
            }
        } else {
            for (i = 0; i < x1->dim; i++) {
                x->values_cpx[i] = x1->values_cpx[i];
            }
        }
        if (MESS_IS_REAL(x2)) {
            for (i = 0; i < x2->dim; i++) {
                x->values_cpx[i + x1->dim] = x2->values[i];
            }
        } else {
            for (i = 0; i < x2->dim; i++) {
                x->values_cpx[i + x1->dim] = x2->values_cpx[i];
            }
        }
    } else {
        for (i = 0; i < x1->dim; i++) {
            x->values[i] = x1->values[i];
        }
        for (i = 0; i < x2->dim; i++) {
            x->values[x1->dim + i] = x2->values[i];
        }
    }
    return 0;
}

/**
 * @brief Split a vector after the n-th entry.
 * @param[in] input     input vector
 * @param[in] n         input size of the first output \f$n\f$
 * @param[out] x1       first output \f$x_1\f$
 * @param[out] x2       second output \f$x_2\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_split function splits a vector into two vectors. The split is done
 * after the \f$n\f$-th element. In this way \f$ dim(x_1)=n \f$ and
 * \f[ \begin{array}{ccc}
 *  x_1 &=& input(0:n-1), \\
 *  x_2 &=& input(n:end).
 * \end{array} \f]
 * If  \f$ x_1  \f$ points to @c NULL,  \f$ x_1  \f$ is neglected and  \f$ x_2 = input(n:end) \f$. \n
 * If  \f$ x_2  \f$ points to @c NULL,  \f$ x_2  \f$ is neglected and  \f$ x_1 = input(0:n-1)  \f$.
 *
 */
int mess_vector_split(mess_vector input, mess_int_t n, mess_vector x1,mess_vector x2) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t j;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(input);
    if(!x1 && !x2){
        MSG_ERROR("x1 and x2 point to null!\n.");
        return MESS_ERROR_NULLPOINTER;
    }

    n = (n > input->dim) ? (input->dim) : n;

    if(x1){
        if ( MESS_IS_COMPLEX(input) && !MESS_IS_COMPLEX(x1)) {
            ret = mess_vector_tocomplex(x1);
            FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_tocomplex);
        }
        ret = mess_vector_resize(x1, n);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
        if (MESS_IS_REAL(input) && MESS_IS_REAL(x1)) {
            for (j = 0; j < n; j++) {
                x1->values[j] = input->values[j];
            }
        } else if (MESS_IS_REAL(input) && MESS_IS_COMPLEX(x1)) {
            for (j = 0; j < n; j++) {
                x1->values_cpx[j] = input->values[j];
            }
        } else {
            for (j = 0; j < n; j++) {
                x1->values_cpx[j] = input->values_cpx[j];
            }
        }
    }

    if(x2){
        if ( MESS_IS_COMPLEX(input) && !MESS_IS_COMPLEX(x2)) {
            ret = mess_vector_tocomplex(x2);
            FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_tocomplex);
        }
        ret = mess_vector_resize(x2, input->dim - n);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
        if (MESS_IS_REAL(input) && MESS_IS_REAL(x2)) {
            for (j = 0; j < input->dim - n; j++) {
                x2->values[j] = input->values[n + j];
            }
        } else if (MESS_IS_REAL(input) && MESS_IS_COMPLEX(x2)) {
            for (j = 0; j < input->dim - n; j++) {
                x2->values_cpx[j] = input->values[n + j];
            }
        } else {
            for (j = 0; j < input->dim - n; j++) {
                x2->values_cpx[j] = input->values_cpx[n + j];
            }
        }
    }
    return 0;
} /* -----  end of function mess_vector_split  ----- */




/**
 * @brief Add \f$ n \f$ zeros to vector.
 * @param[in]       in input vector
 * @param[in]        n input positive number of zeros to add
 * @param[out]     out output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_lift function adds \f$ n  \f$ zeros at the end of the vector \f$ in  \f$. \n
 * The result is in vector \f$ out \f$.
 *
 */
int mess_vector_lift(mess_vector in, mess_int_t n,  mess_vector out ){
    MSG_FNAME(__func__);
    int ret=0;

    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_positive(n);
    mess_check_nullpointer(in);
    mess_check_real_or_complex(in);
    mess_check_nullpointer(out);

    /*-----------------------------------------------------------------------------
     *  perform [in;zeros(n,1)]->out
     *-----------------------------------------------------------------------------*/
    //ret = mess_vector_zeros(out);                         FUNCTION_FAILURE_HANDLE(ret,(ret=0),mess_vector_zeros);

    if(out->dim!=n+in->dim){
        ret = mess_vector_resize(out,n+in->dim);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);
    }

    if(MESS_IS_REAL(in)){
        if(MESS_IS_COMPLEX(out)){
            ret = mess_vector_toreal_nowarn(out);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal_nowarn);
        }
        F77_GLOBAL(dcopy,DCOPY)(&(in->dim), in->values, &__ONE, out->values, &__ONE);
    }else{
        if(MESS_IS_REAL(out)){
            ret = mess_vector_tocomplex(out);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        }
        F77_GLOBAL(zcopy,ZCOPY)(&(in->dim), in->values_cpx, &__ONE, out->values_cpx, &__ONE);
    }
    return ret;
}

