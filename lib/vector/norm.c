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
 * @file lib/vector/norm.c
 * @brief Norm functions of @ref mess_vector structures.
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
 * @brief Compute the \f$ 2 \f$-norm of a vector.
 * @param[in] x         input vector
 * @param[out] nrm      output computed \f$2\f$-norm of \f$x\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_norm2 function computes the \f$2\f$-norm of a vector  \f$x\f$
 * \f[nrm = \Vert x \Vert_2 . \f]
 *
 * @sa mess_vector_norm1
 * @sa mess_vector_norminf
 */
int mess_vector_norm2(mess_vector x, double* nrm) {
    MSG_FNAME(__func__);
    mess_int_t dim;

    mess_check_nullpointer(x);

    dim = x->dim;
    if (MESS_IS_REAL(x)) {
        *nrm = F77_GLOBAL(dnrm2,DNRM2)(&dim, x->values, &__ONE);
    } else if (MESS_IS_COMPLEX(x)) {
        *nrm = F77_GLOBAL(dznrm2,DZNRM2)(&dim, x->values_cpx, &__ONE);
    } else {
        MSG_ERROR("unknown/unsupported data type\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
} /* -----  end of function mess_vector_norm2  ----- */

/**
 * @brief Compute the \f$ 1 \f$-norm of a vector.
 * @param[in] x         input vector
 * @param[out] nrm   output computed  \f$1\f$-norm of \f$x\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_norm1 function computes the \f$1\f$-norm of a vector \f$x\f$
 * \f[ nrm =  \Vert x \Vert_1 . \f]
 */
int mess_vector_norm1(mess_vector x, double* nrm) {
    MSG_FNAME(__func__);
    mess_int_t dim;
    mess_int_t j = 0;

    mess_check_nullpointer(x);
    mess_check_nullpointer(nrm);

    dim = x->dim;
    *nrm = 0;
    if (MESS_IS_REAL(x)) {
        for (j = 0; j < dim; ++j) {
            *nrm += fabs(x->values[j]);
        }
    } else if (MESS_IS_COMPLEX(x)) {
        for (j = 0; j < dim; ++j) {
            *nrm += cabs(x->values_cpx[j]);
        }
    } else {
        MSG_ERROR("unknown/unsupported data type\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}

/**
 * @brief Compute the \f$\infty\f$-norm of a vector.
 * @param[in] x         input vector
 * @param[out] nrm      output computed \f$ \infty \f$-norm
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_norminf function computes the \f$\infty\f$-norm of a vector \f$x\f$
 * \f[nrm =  \Vert x \Vert_{\infty} . \f]
 */
int mess_vector_norminf(mess_vector x, double *nrm) {
    MSG_FNAME(__func__);
    mess_int_t dim, i;
    double mx = 0;

    mess_check_nullpointer(x);
    mess_check_nullpointer(nrm);

    dim = x->dim;
    if (MESS_IS_REAL(x)) {
        mx = fabs(x->values[0]);
        for (i = 0; i < dim; i++) {
            if (mx < fabs(x->values[i]))
                mx = fabs(x->values[i]);
        }
    } else if (MESS_IS_COMPLEX(x)) {
        mx = cabs(x->values_cpx[0]);
        for (i = 0; i < dim; i++) {
            if (mx < cabs(x->values_cpx[i]))
                mx = cabs(x->values_cpx[i]);
        }
    } else {
        MSG_ERROR("unknown/unsupported data type\n");
        return MESS_ERROR_DATATYPE;
    }

    *nrm = mx;
    return 0;
}


/**
 * @brief Compute the norm of a vector.
 * @param[in] x         input vector \f$x\f$
 * @param[in] nrm_t     input @ref mess_norm_t the desired type of norm
 * @param[out] nrm      output  \f$ \Vert x \Vert \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_vector_norm function calculates the norm of a @ref mess_vector
 * Supported norms are:
 *
 * @ref MESS_2_NORM, @ref MESS_FROBENIUS_NORM, @ref MESS_1_NORM, @ref MESS_INF_NORM.
 *
 * Here @ref MESS_FROBENIUS_NORM and @ref MESS_2_NORM is treated as the same.
 *
 * @sa mess_vector_norm2
 * @sa mess_vector_norminf
 * @sa mess_vector_norm1
 */
int mess_vector_norm(mess_vector x, mess_norm_t nrm_t, double *nrm){
    MSG_FNAME(__func__);
    int ret;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(x);
    mess_check_real_or_complex(x);
    mess_check_nullpointer(nrm);

    /*-----------------------------------------------------------------------------
     *  compute norm
     *-----------------------------------------------------------------------------*/
    switch (nrm_t) {
        case MESS_2_NORM:
        case MESS_FROBENIUS_NORM:
            ret = mess_vector_norm2(x, nrm);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_norm2);
            break;
        case MESS_1_NORM:
            ret = mess_vector_norm1(x, nrm);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_norm1);
            break;
        case MESS_INF_NORM:
            ret = mess_vector_norminf(x, nrm);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_norminf);
            break;
        default:
            MSG_ERROR("unkown/unsupported norm type: %s\n", mess_norm_t_str(nrm_t));
            return MESS_ERROR_NOSUPPORT;
    }
    return(ret);
}



/**
 * @brief Compute the dot product of two vectors.
 * @param[in] x input vector \f$x\f$
 * @param[in] y input vector \f$y\f$
 * @param[out] dot output dot product
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_dot function computes the dot product of a real vector \f$ x \f$ and a vector \f$ y \f$, i.e.
 * \f[dot=x^T y. \f]
 *
 *
 */
int mess_vector_dot(mess_vector x, mess_vector y, double*  dot) {
    MSG_FNAME(__func__);
    mess_int_t dim;

    mess_check_nullpointer(x);
    mess_check_nullpointer(y);

    if (x->dim != y->dim) {
        MSG_ERROR("dimension mismatch\n");
        return MESS_ERROR_DIMENSION;
    }
    if (x->data_type != y->data_type) {
        MSG_ERROR("data type mismatch\n");
        return MESS_ERROR_DATATYPE;
    }
    if (MESS_IS_COMPLEX(x)) {
        MSG_ERROR("to compute the complex dot product use mess_vector_dotc\n");
        return 0;
    }

    dim = x->dim;
    if (MESS_IS_REAL(x)) {
        *dot = F77_GLOBAL(ddot,DDOT)(&dim, x->values, &__ONE, y->values, &__ONE);
        return 0;
    } else {
        MSG_ERROR("unknown/unsupported data type\n");
        return MESS_ERROR_DATATYPE;
    }
}

/**
 * @brief Compute the dot product of two vectors (complex).
 * @param[in] x input vector \f$x\f$
 * @param[in] y input vector \f$y\f$
 * @param[out] dot output dot product
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_dotc function computes the dot product of two vectors, i.e.
 * \f[dot = x^H y. \f]
 */
int mess_vector_dotc(mess_vector x, mess_vector y, mess_double_cpx_t * dot) {
    MSG_FNAME(__func__);mess_int_t dim;

    mess_check_nullpointer(x);mess_check_nullpointer(y);

    if (x->dim != y->dim){
        MSG_ERROR("dimension mismatch\n");return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  handle different cases
     *-----------------------------------------------------------------------------*/

    dim = x->dim;
    if( MESS_IS_REAL(x) && MESS_IS_REAL(y) ){
        *dot = F77_GLOBAL(ddot,DDOT)(&dim, x->values, &__ONE, y->values, &__ONE);
    }else if ( MESS_IS_REAL(x) && MESS_IS_COMPLEX(y) ) {
        *dot = F77_GLOBAL(dzdotu,DZDOTU)(&dim, x->values, &__ONE, y->values_cpx, &__ONE);
    }else if ( MESS_IS_COMPLEX(x) && MESS_IS_REAL(y) ) {
        *dot = F77_GLOBAL(zddotc,ZDDOTC)(&dim, x->values_cpx, &__ONE, y->values, &__ONE);
    }else if ( MESS_IS_COMPLEX(x) && MESS_IS_COMPLEX(y) ) {
#ifdef MESS_USE_APPLE_BLAS
        cblas_zdotc_sub(dim, x->values_cpx, 1, y->values_cpx,1, dot);
#else
#ifdef ZDOTC_MKL
        F77_GLOBAL(zdotc,ZDOTC)(dot, &dim, x->values_cpx, &__ONE, y->values_cpx, &__ONE);
#else
        *dot = F77_GLOBAL(zdotc,ZDOTC)(&dim, x->values_cpx, &__ONE, y->values_cpx, &__ONE);
#endif
#endif
    }else{
        MSG_ERROR("unknown/unsupported data type\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}


/**
 * @brief Compute the dot product of two vectors (complex).
 * @param[in] x input vector \f$x\f$
 * @param[in] y input vector \f$y\f$
 * @param[out] dot output dot product
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_dotu function computes the dot product of two vectors, i.e.
 * \f[dot = x^T y. \f]
 */
int mess_vector_dotu(mess_vector x, mess_vector y, mess_double_cpx_t * dot) {
    MSG_FNAME(__func__);
    mess_int_t dim;

    mess_check_nullpointer(x);
    mess_check_nullpointer(y);

    if (x->dim != y->dim){
        MSG_ERROR("dimension mismatch\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  handle different cases
     *-----------------------------------------------------------------------------*/

    dim = x->dim;
    if( MESS_IS_REAL(x) && MESS_IS_REAL(y) ){
        *dot = F77_GLOBAL(ddot,DDOT)(&dim, x->values, &__ONE, y->values, &__ONE);
    }else if ( MESS_IS_REAL(x) && MESS_IS_COMPLEX(y) ) {
        *dot = F77_GLOBAL(dzdotu,DZDOTU)(&dim, x->values, &__ONE, y->values_cpx, &__ONE);
    }else if ( MESS_IS_COMPLEX(x) && MESS_IS_REAL(y) ) {
        *dot = F77_GLOBAL(zddotu,ZDDOTU)(&dim, x->values_cpx, &__ONE, y->values, &__ONE);
    }else if ( MESS_IS_COMPLEX(x) && MESS_IS_COMPLEX(y) ) {
#ifdef MESS_USE_APPLE_BLAS
        cblas_zdotu_sub(dim, x->values_cpx, 1, y->values_cpx,1, dot);
#else
#ifdef ZDOTC_MKL
        F77_GLOBAL(zdotu,ZDOTU)(dot, &dim, x->values_cpx, &__ONE, y->values_cpx, &__ONE);
#else
        *dot = F77_GLOBAL(zdotu,ZDOTU)(&dim, x->values_cpx, &__ONE, y->values_cpx, &__ONE);
#endif
#endif
    }else{
        MSG_ERROR("unknown/unsupported data type\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}




/**
 * @brief Compute \f$ \Vert x_1-x_2 \Vert_2 \f$.
 * @param[in] x1        input vector \f$x_1\f$
 * @param[in] x2        input vector \f$x_2\f$
 * @param[out] diff     output difference
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_diffnorm function computes
 * \f[ diff = \Vert x_1-x_2 \Vert_2 . \f]
 *
 */
int mess_vector_diffnorm(mess_vector x1, mess_vector x2, double *diff) {
    MSG_FNAME(__func__);
    mess_int_t i;
    double sum = 0;
    double t = 0;

    *diff = 0;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(x1);
    mess_check_nullpointer(x2);
    mess_check_nullpointer(diff);

    if (x1->dim != x2->dim) {
        MSG_ERROR("inputs must have the same dimension.\n");
        return MESS_ERROR_DIMENSION;
    }

    if ( MESS_IS_COMPLEX(x1) && MESS_IS_COMPLEX(x2)) {
        for (i = 0; i < x1->dim; i++) {
            t = cabs(x1->values_cpx[i] - x2->values_cpx[i]);
            sum += (t * t);
        }
    } else if ( MESS_IS_COMPLEX(x1) && MESS_IS_REAL(x2)) {
        for (i = 0; i < x1->dim; i++) {
            t = cabs(x1->values_cpx[i] - x2->values[i]);
            sum += (t * t);
        }
    } else if ( MESS_IS_REAL(x1) && MESS_IS_COMPLEX(x2)) {
        for (i = 0; i < x1->dim; i++) {
            t = cabs((mess_double_cpx_t)x1->values[i] - x2->values_cpx[i]);
            sum += (t * t);
        }
    } else if ( MESS_IS_REAL(x1) && MESS_IS_REAL(x2)) {
        for (i = 0; i < x1->dim; i++) {
            t = fabs(x1->values[i] - x2->values[i]);
            sum += (t * t);
        }
    } else {
        MSG_ERROR("unknown datatype\n.");
        return MESS_ERROR_DATATYPE;
    }

    *diff = sqrt(sum);

    return 0;
} /* -----  end of function mess_vector_diffnorm  ----- */

/**
 * @brief Compute \f$ \Vert x_1-x_2 \Vert_{\infty} \f$.
 * @param[in] x1        input vector \f$x_1\f$
 * @param[in] x2        input vector two
 * @param[out] diff     output difference
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_diffnorminf function computes
 * \f[ diff = \Vert x_1-x_2 \Vert_{\infty} . \f]
 *
 */
int mess_vector_diffnorminf(mess_vector x1, mess_vector x2, double *diff) {
    MSG_FNAME(__func__);
    mess_int_t i;
    double sum = 0;
    double t = 0;

    *diff = 0;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(x1);
    mess_check_nullpointer(x2);
    mess_check_nullpointer(diff);

    if (x1->dim != x2->dim) {
        MSG_ERROR("inputs must have the same dimension.\n");
        return MESS_ERROR_DIMENSION;
    }

    if ( MESS_IS_COMPLEX(x1) && MESS_IS_COMPLEX(x2)) {
        for (i = 0; i < x1->dim; i++) {
            t = cabs(x1->values_cpx[i] - x2->values_cpx[i]);
            sum = MESS_MAX(sum, t);
        }
    } else if ( MESS_IS_COMPLEX(x1) && MESS_IS_REAL(x2)) {
        for (i = 0; i < x1->dim; i++) {
            t = cabs(x1->values_cpx[i] - x2->values[i]);
            sum = MESS_MAX(sum, t);
        }
    } else if ( MESS_IS_REAL(x1) && MESS_IS_COMPLEX(x2)) {
        for (i = 0; i < x1->dim; i++) {
            t = cabs(x1->values[i] - x2->values_cpx[i]);
            sum = MESS_MAX(sum, t);
        }
    } else if ( MESS_IS_REAL(x1) && MESS_IS_REAL(x2)) {
        for (i = 0; i < x1->dim; i++) {
            t = fabs(x1->values[i] - x2->values[i]);
            sum = MESS_MAX(sum, t);
        }
    } else {
        MSG_ERROR("unknown datatype\n.");
        return MESS_ERROR_DATATYPE;
    }

    *diff = sum;

    return 0;
} /* -----  end of function mess_vector_diffnorm  ----- */


