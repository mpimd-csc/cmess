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
 * @file lib/direct/res2/lyap_res2.c
 * @brief Compute the \f$ 2 \f$-norm residual of a Lyapunov equation.
 * @author @koehlerm
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
 * @brief Compute the \f$ 2 \f$-norm residual of a Lyapunov Equation (DENSE).
 * @param[in] op  input operation on matrix A
 * @param[in] A  input system matrix of the Lyapunov equation
 * @param[in] M  input right hand side of the Lyapunov Equation
 * @param[in] X  input solution of the Lyapunov Equation
 * @param[out] res two-norm of the residual
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_lyapunov_res2 function computes the \f$ 2 \f$-norm residual of a Lyapunov Equation:
 * \f[ op(A)X + X op(A)^T +M =0 \f]
 * that means it computes:
 * \f[ res= \Vert op(A)X + Xop(A)^T + M \Vert_2,\f]
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
int mess_direct_lyapunov_res2 ( mess_operation_t op, mess_matrix A, mess_matrix M, mess_matrix X, double *res )
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


    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&T);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0) ,mess_matrix_init);
    ret = mess_matrix_init(&T2);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0) ,mess_matrix_init);
    ret = mess_matrix_copy(M, T);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);
    if (op == MESS_OP_NONE) {
        ret = mess_matrix_multiply(MESS_OP_NONE, A, MESS_OP_NONE, X, T2);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
        ret = mess_matrix_add(1.0,T2,1.0,T);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_multiply(MESS_OP_NONE, X, MESS_OP_TRANSPOSE, A, T2);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_add(1.0,T2,1.0,T);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    } else {
        ret = mess_matrix_multiply(MESS_OP_TRANSPOSE, A, MESS_OP_NONE, X, T2);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
        ret = mess_matrix_add(1.0,T2,1.0,T);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_multiply(MESS_OP_NONE, X, MESS_OP_NONE, A, T2);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_add(1.0,T2,1.0,T);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    }
    ret = mess_matrix_norm2(T, &nrm);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm2);
    *res = nrm;

    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
    mess_matrix_clear(&T);
    mess_matrix_clear(&T2);

    return 0;
}

/**
 * @brief Compute the \f$ 2 \f$-norm residual of a generalized Lyapunov Equation.
 * @param[in] op  input operation on matrices \f$ A \f$ and \f$ E \f$
 * @param[in] A  input matrix A
 * @param[in] E  input matrix E
 * @param[in] M  input right hand side
 * @param[in] X  input solution
 * @param[out] res2 output residual
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_generalized_lyapunov_res2 function computes the \f$ 2 \f$-norm resiudal
 * of a generalized Lyapunov Equation:
 * \f[ op(A)Xop(E)^T +op(E)Xop(A)^T +M=0 \f]
 * that means it computes:
 * \f[ res2= \Vert op(A)Xop(E)^T + op(E)Xop(A)^T + M \Vert_2, \f]
 * where \f$ op \f$ can be
 * <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 * </center>
 *
 */
int  mess_direct_generalized_lyapunov_res2( mess_operation_t op, mess_matrix A, mess_matrix E, mess_matrix M , mess_matrix X, double *res2 )
{
    MSG_FNAME(__func__);
    int ret;
    double r;
    mess_matrix Y, T, T1, T2;

    /*-----------------------------------------------------------------------------
     *  call standard residual routine if E points to NULL
     *-----------------------------------------------------------------------------*/
    if(!E){
        return mess_direct_lyapunov_res2(op, A, M,  X, res2);
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
    mess_check_square(A);
    mess_check_square(E);
    mess_check_square(M);
    mess_check_square(X);

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&Y);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&T);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&T1);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&T2);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_convert(M,Y,MESS_DENSE);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    if ( op == MESS_OP_NONE ) {
        ret = mess_matrix_multiply(MESS_OP_NONE, A, MESS_OP_NONE, X, T);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, T, MESS_OP_TRANSPOSE, E, T1);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, E, MESS_OP_TRANSPOSE, T, T2);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_add(1.0,T1,1.0,Y);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_add);
        ret = mess_matrix_add(1.0,T2,1.0,Y);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_add);
    } else {
        ret = mess_matrix_multiply(MESS_OP_TRANSPOSE, A, MESS_OP_NONE, X, T);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, T, MESS_OP_NONE, E, T1);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_TRANSPOSE, E, MESS_OP_TRANSPOSE, T, T2);FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_add(1.0,T1,1.0,Y);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_add);
        ret = mess_matrix_add(1.0,T2,1.0,Y);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_add);
    }
    /*-----------------------------------------------------------------------------
     *  compute the norm
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_norm2(Y, &r); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm2);
    *res2 = r;

    mess_matrix_clear(&T);
    mess_matrix_clear(&T1);
    mess_matrix_clear(&T2);
    mess_matrix_clear(&Y);

    return 0;
}       /* -----  end of function mess_direct_generalized_lyapunov_res2  ----- */


