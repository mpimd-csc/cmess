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
 * @file lib/matrix/colrowsums.c
 * @brief Compute the column or row sum of a matrix.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


/**
 * @brief Compute the row sums of a matrix.
 * @param[in] A     input matrix \f$A\f$
 * @param[out] x    output vector \f$x\f$ with the row sums of \f$A\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_rowsums function computes the sum of each row of a matrix by evaluating
 * \f[ x= A e \f]
 * where \f$ e = [1,...,1]^T \f$.
 * @sa mess_matrix_colsums
 * */
int mess_matrix_rowsums(mess_matrix A, mess_vector x){
    MSG_FNAME(__func__);
    mess_vector tmp;
    int ret = 0;

    mess_check_nullpointer(A);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(A);


    if (x->dim != A->rows){
        MSG_WARN("dimension mismatch x = " MESS_PRINTF_INT " , matrix-rows = " MESS_PRINTF_INT "\n", x->dim, A->rows);
        ret = mess_vector_resize(x, A->rows); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }
    if ( MESS_IS_REAL (A) ){
        ret = mess_vector_toreal_nowarn(x);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        ret = mess_vector_init(&tmp);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(tmp, A->cols, A->data_type );   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_ones(tmp);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);
        ret = mess_matrix_mvp(MESS_OP_NONE,A, tmp, x);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    }else {
        ret = mess_vector_tocomplex(x);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_init(&tmp);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(tmp, A->cols, A->data_type );   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_ones(tmp);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);
        ret = mess_matrix_mvp(MESS_OP_NONE,A, tmp, x);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    }
    mess_vector_clear(&tmp);
    return(0);
}

/**
 * @brief Compute the column sums of a matrix.
 * @param[in] A     input matrix \f$A\f$
 * @param[out] x    output vector \f$x\f$ with the column sums of \f$A\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_colsums function computes the sum of each column of a matrix by evaluating
 * \f[ x= A^T e \f]
 * where \f$ e = [1,...,1]^T \f$.
 * @sa mess_matrix_rowsums
 */
int mess_matrix_colsums(mess_matrix A, mess_vector x){
    MSG_FNAME(__func__);
    mess_vector tmp;
    int ret = 0;

    mess_check_nullpointer(A);
    mess_check_nullpointer(x);

    if (x->dim != A->cols){
        MSG_WARN("dimension mismatch x= " MESS_PRINTF_INT " matrix->cols = " MESS_PRINTF_INT "  \n", x->dim, A->cols );
        ret = mess_vector_resize(x, A->rows); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }

    if ( MESS_IS_REAL (A) ){
        ret = mess_vector_toreal_nowarn(x);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        ret = mess_vector_init(&tmp);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(tmp, A->rows, A->data_type );   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_ones(tmp);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);
        ret = mess_matrix_mvp(MESS_OP_HERMITIAN,A, tmp, x);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    }else {
        ret = mess_vector_tocomplex ( x ) ;                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_init(&tmp);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mesS_vector_init);
        ret = mess_vector_alloc(tmp, A->rows, A->data_type );   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_ones(tmp);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);
        ret = mess_matrix_mvp(MESS_OP_HERMITIAN,A, tmp, x);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    }
    mess_vector_clear(&tmp);
    return(0);
}

