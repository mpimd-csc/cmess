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
 * @file lib/direct/res2/res2.c
 * @brief \f$ 2 \f$-norm residual of \f$ r=Ax-b \f$ and similar equations.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "complex.h"


/**
 * @brief Compute the (relative) residual of \f$ Ax=b \f$.
 * @param[in] op  input operation on matrix A
 * @param[in] A  input system matrix A
 * @param[in] x  input x vector
 * @param[in] b  input b vector
 * @param[out] resid residual
 * @param[out] relresid relative residual
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_res2 function computes the  absolute and the relative residual of
 * \f[ op(A)x=b \f]
 * where \f$ op \f$ can be
 * <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 *  |@ref MESS_OP_HERMITIAN      |   \f$ op(A)=A^H\f$         |
 * </center>
 *
 * \sa mess_direct_res2_shifted
 * \sa mess_direct_res2m
 */
int  mess_direct_res2( mess_operation_t op,  mess_matrix A, mess_vector x, mess_vector b, double *resid, double * relresid )
{
    MSG_FNAME(__func__);
    mess_vector t1, t2;
    int ret = 0;
    double nrmb = 0;

    mess_check_nullpointer(A);
    mess_check_nullpointer(x);
    mess_check_nullpointer(b);
    mess_check_nullpointer(resid);
    mess_check_nullpointer(relresid);

    if ( x->dim != A->rows) {
        MSG_ERROR("x has the wrong dimension. \n");
        return MESS_ERROR_DIMENSION;
    }
    if ( b->dim != A->rows) {
        MSG_ERROR("b has the wrong dimension. \n");
        return MESS_ERROR_DIMENSION;
    }

    ret = mess_vector_norm2(b, &nrmb);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    MESS_INIT_VECTORS(&t1,&t2);
    ret = mess_vector_alloc(t1, x->dim, x->data_type);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(t2, x->dim, x->data_type);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);

    ret = mess_matrix_mvp(op,A,x,t1);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_copy(b,t2);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
    ret = mess_vector_axpy(-1,t1,t2);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
    ret = mess_vector_norm2(t2, resid);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_norm2);
    *relresid = *resid / nrmb;

    mess_vector_clear(&t1);
    mess_vector_clear(&t2);
    return 0;
}       /* -----  end of function mess_direct_res2  ----- */


/**
 * @brief Compute the (relative) residual of \f$ Ax=b \f$.
 * @param[in] op  input operation on matrix A
 * @param[in] A  input system matrix A
 * @param[in] x  input x matrix
 * @param[in] b  input b matrix
 * @param[out] resid residual
 * @param[out] relresid relative residual
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_res2m function computes the  absolute and the relative residual of
 * \f[ op(A)x=b \f]
 * where \f$ op \f$ can be
 * <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 *  |@ref MESS_OP_HERMITIAN      |   \f$ op(A)=A^H\f$         |
 * </center>
 *
 * \sa mess_direct_res2_shifted
 * \sa mess_direct_res2
 *
 */
int mess_direct_res2m(mess_operation_t op, mess_matrix A, mess_matrix x, mess_matrix b, double *resid, double * relresid)
{
    MSG_FNAME(__func__);
    mess_matrix t1, t2;
    int ret = 0;
    double nrmb = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(x);
    mess_check_nullpointer(b);
    mess_check_nullpointer(resid);
    mess_check_nullpointer(relresid);

    if ( x->rows != A->rows) {
        MSG_ERROR("x has the wrong nuber of rows.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( b->rows != A->rows) {
        MSG_ERROR("b has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  compute residuals
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_norm2(b, &nrmb);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm2);
    MESS_INIT_MATRICES(&t1, &t2);
    ret = mess_matrix_copy(x, t1);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_zeros(t1);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);

    ret = mess_matrix_copy(x, t2);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_zeros(t2);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);

    ret = mess_matrix_multiply(op, A, MESS_OP_NONE, x, t1);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_copy(b,t2);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_add(-1,t1,1.0,t2);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    ret = mess_matrix_norm2(t2, resid);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm2);
    *relresid = *resid / nrmb;

    MESS_CLEAR_MATRICES(&t1, &t2);
    return 0;
}       /* -----  end of function mess_direct_res2m  ----- */


