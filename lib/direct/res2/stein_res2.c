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
 * @file lib/direct/res2/stein_res2.c
 * @brief Compute the \f$ 2 \f$-norm residual of a Stein equation.
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Compute the \f$ 2 \f$-norm residual of a Stein Equation (dense).
 * @param[in] op    input operation on matrix A
 * @param[in] A     input system matrix of the Stein equation
 * @param[in] M     input right hand side of the Stein equation
 * @param[in] X     input solution of the Stein equation
 * @param[out] res two-norm of the residual
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_stein_res2 function computes the \f$ 2 \f$-norm residual of a Stein Equation:
 * \f[ op(A)X op(A)^T - X  + M =0 \f]
 * that means it computes:
 * \f[ res= \Vert  op(A)X op(A)^T - X  + M  \Vert_2,\f]
 * where \f$ op \f$ can be
 * <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 * </center>
 * All computations are done in dense arithmetics. Please do not use this function for
 * large-scale equations.
 *
 */
int mess_direct_stein_res2 ( mess_operation_t op, mess_matrix A, mess_matrix M, mess_matrix X, double *res )
{
    MSG_FNAME(__func__);
    mess_matrix T;
    mess_matrix T2;
    int ret = 0;
    double nrm = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(M);
    mess_check_nullpointer(X);
    mess_check_nullpointer(res);
    mess_check_real(A);
    mess_check_real(X);
    mess_check_real(M);
    mess_check_dense(A);
    mess_check_dense(X);
    mess_check_dense(M);
    mess_check_square(A);
    mess_check_square(M);
    mess_check_square(X);

    mess_check_same_size(A,M);
    mess_check_same_size(A,X);

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&T,&T2);
    ret = mess_matrix_copy(M, T);                                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);

    if (op == MESS_OP_NONE) {
        ret = mess_matrix_multiply(MESS_OP_NONE, A, MESS_OP_NONE, X, T);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, T, MESS_OP_TRANSPOSE, A, T2);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_add(-1.0,X,1.0,T2);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_copy(M,T);                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        ret = mess_matrix_add(1.0,T2,1.0,T);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    } else {
        ret = mess_matrix_multiply(MESS_OP_TRANSPOSE, A, MESS_OP_NONE, X, T);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, T, MESS_OP_NONE, A, T2);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_add(-1.0,X,1.0,T2);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_copy(M,T);                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        ret = mess_matrix_add(1.0,T2,1.0,T);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    }

    /*-----------------------------------------------------------------------------
     *  compute the norm
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_norm2(T, &nrm);                                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm2);
    *res = nrm;

    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&T,&T2);

    return 0;
}

/**
 * @brief Compute the \f$ 2 \f$-norm residual of a generalized Stein Equation.
 * @param[in] op  input operation on matrices \f$ A \f$ and \f$ E \f$
 * @param[in] A  input matrix A
 * @param[in] E  input matrix E
 * @param[in] M  input right hand side
 * @param[in] X  input solution
 * @param[out] res2 output residual
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_generalized_stein_res2 function computes the \f$ 2 \f$-norm resiudal
 * of a generalized Stein Equation:
 * \f[ op(A)Xop(A)^T - op(E)Xop(E)^T +M=0 \f]
 * that means it computes:
 * \f[ res2= \Vert op(A)Xop(A)^T - op(E)Xop(E)^T + M \Vert_2, \f]
 * where \f$ op \f$ can be
 * <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 * </center>
 *
 */
int  mess_direct_generalized_stein_res2 ( mess_operation_t op, mess_matrix A, mess_matrix E, mess_matrix M , mess_matrix X, double *res2 )
{
    MSG_FNAME(__func__);
    int ret;
    double r;
    mess_matrix Y, T, T1, T2;

    /*-----------------------------------------------------------------------------
     *  call standard residual routine if E points to NULL
     *-----------------------------------------------------------------------------*/
    if(!E){
        return mess_direct_stein_res2(op, A, M, X, res2);
    }

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(E);
    mess_check_nullpointer(M);
    mess_check_nullpointer(X);
    mess_check_nullpointer(res2);
    mess_check_real(A);
    mess_check_real(E);
    mess_check_real(M);
    mess_check_real(X);
    mess_check_dense(A);
    mess_check_dense(E);
    mess_check_dense(M);
    mess_check_dense(X);
    mess_check_square(A);
    mess_check_square(E);
    mess_check_square(M);
    mess_check_square(X);
    mess_check_same_size(A,E);
    mess_check_same_size(A,M);
    mess_check_same_size(A,X);

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&Y,&T,&T1,&T2);

    if (op == MESS_OP_NONE) {
        ret = mess_matrix_multiply(MESS_OP_NONE, A, MESS_OP_NONE, X, T);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, T, MESS_OP_TRANSPOSE, A, T2);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, E, MESS_OP_NONE, X, T);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, T, MESS_OP_TRANSPOSE, E, Y);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_add(-1.0,Y,1.0,T2);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_copy(M,T);                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        ret = mess_matrix_add(1.0,T2,1.0,T);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    } else {
        ret = mess_matrix_multiply(MESS_OP_TRANSPOSE, A, MESS_OP_NONE, X, T);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, T, MESS_OP_NONE, A, T2);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_TRANSPOSE, E, MESS_OP_NONE, X, T);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, T, MESS_OP_NONE, E, Y);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_add(-1.0,Y,1.0,T2);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_copy(M,T);                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        ret = mess_matrix_add(1.0,T2,1.0,T);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    }


    /*-----------------------------------------------------------------------------
     *  compute the norm
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_norm2(T, &r);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm2);
    *res2 = r;

    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&Y,&T,&T1,&T2);


    return 0;
}       /* -----  end of function mess_direct_generalized_stein_res2  ----- */


