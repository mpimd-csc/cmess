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
 * @file lib/matrix/orth.c
 * @brief Orthogonalize a matrix.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/eigenvalue.h"
#include "mess/error_macro.h"
#include <complex.h>

/**
 * @brief Compute an orthogonal basis for the range of a matrix.
 * @param[in] A  input matrix \f$A\f$
 * @param[out] Q output \f$Q\f$  orthogonal basis for the range of \f$A\f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_orth function computes an orthogonal basis for the range
 * of a matrix using the SVD.
 *
 * @sa mess_matrix_rank
 * @sa mess_matrix_null
 *
 */
int mess_matrix_orth ( mess_matrix A, mess_matrix Q )
{
    MSG_FNAME(__func__);
    int conv = 0;
    mess_matrix work;
    mess_vector sigma;
    int ret = 0;
    mess_int_t r=0, i;
    double s,tol;
    double eps = mess_eps();


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(Q);
    mess_check_real_or_complex(A);
    MESS_MATRIX_RESET(Q);
    MESS_MATRIX_CHECKFORMAT(A, work, conv, MESS_DENSE);

    ret = mess_vector_init(&sigma);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sigma, MESS_MAX(work->cols, work->rows), MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_eigen_svd_econ(work, sigma, Q, NULL);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_svd_econ);

    s=sigma->values[0];
    tol = MESS_MAX(work->rows, work->cols) * s * eps;
    for (i=0 ; i < sigma->dim; i++){
        if ( sigma->values[i] > tol ) r++;
        else break;
    }
    r=MESS_MAX(1, r) ;
    ret = mess_matrix_resize(Q, Q->rows, r);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_resize);
    mess_vector_clear(&sigma);

    if ( conv == 0) {
        mess_matrix_clear(&work);
    }

    return(0);
}       /* -----  end of function mess_matrix_orth  ----- */

/**
 * @brief Compute an orthogonal basis for the nullspace of a matrix.
 * @param[in] A  input matrix \f$A\f$
 * @param[out] Z output \f$Z\f$  orthogonal basis for the nullspace of \f$A \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_null function computes an orthogonal basis for the nullspace
 * of a matrix using the SVD.
 *
 * @sa mess_matrix_rank
 * @sa mess_matrix_orth
 */
int mess_matrix_null ( mess_matrix A, mess_matrix Z )
{
    MSG_FNAME(__func__);
    int conv = 0;
    mess_matrix work;
    mess_matrix Q;
    mess_vector sigma;
    int ret = 0;
    mess_int_t r=0, i;
    double s,tol;
    double eps = mess_eps();


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(Z);
    mess_check_real_or_complex(A);
    MESS_MATRIX_RESET(Z);
    MESS_MATRIX_CHECKFORMAT(A, work, conv, MESS_DENSE);

    ret = mess_matrix_init(&Q);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_vector_init(&sigma);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sigma, MESS_MAX(work->cols, work->rows), MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_eigen_svd_econ(work, sigma, NULL, Q);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_svd_econ);

    s=sigma->values[0];
    tol = MESS_MAX(work->rows, work->cols) * s * eps;
    for (i=0 ; i < sigma->dim; i++){
        if ( sigma->values[i] > tol ) r++;
        else break;
    }
    r=MESS_MAX(1, r) ;


    ret = mess_matrix_colsub(Q,r,Q->cols-1, Z); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_null);

    if ( conv == 0) {
        mess_matrix_clear(&work);
    }
    mess_matrix_clear(&Q);
    mess_vector_clear(&sigma);

    return(0);
}

