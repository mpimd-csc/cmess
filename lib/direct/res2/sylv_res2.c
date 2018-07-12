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
 * @file lib/direct/res2/sylv_res2.c
 * @brief Residual computations for dense Sylvester equations.
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
 * @brief Compute the \f$ 2 \f$-norm residual of a standard Sylvester equation.
 * @param[in] op  input operation on matrices \f$ A \f$ and \f$ H \f$
 * @param[in] A  input first matrix
 * @param[in] H  input second matrix
 * @param[in] M  input right hand side matrix
 * @param[in] X  input solution matrix
 * @param[out] res  norm
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_sylvester_res2 function computes the 2-Norm residual of the Sylvester equation
 * \f[ R:= op(A)X+Xop(H)+M =0 \f]
 * that means it computes
 * \f[ res= \Vert op(A)X+Xop(H)+M \Vert_2, \f]
 * where \f$ op \f$ can be
 * <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 *  |@ref MESS_OP_HERMITIAN      |   \f$ op(A)=A^H\f$         |
 * </center>
 * It is implemented as dense computation. The matrix \f$ R \f$ is computed
 * inside the function. If the input matrices are too large, this function
 * will be time consuming.
 *
 */
int mess_direct_sylvester_res2 ( mess_operation_t op, mess_matrix A, mess_matrix H, mess_matrix M, mess_matrix X, double *res )
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
    mess_check_nullpointer(H);
    mess_check_nullpointer(M);
    mess_check_nullpointer(X);
    mess_check_nullpointer(res);
    mess_check_square(A);
    mess_check_square(H);
    mess_check_same_size(M,X);
    mess_check_same_rows(A,M);
    mess_check_same_cols(H,M);

    /*-----------------------------------------------------------------------------
     *  perpare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&T);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0) ,mess_matrix_init);
    ret = mess_matrix_init(&T2);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0) ,mess_matrix_init);

    /*-----------------------------------------------------------------------------
     *  compute residual
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_copy(M, T);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_multiply(op, A, MESS_OP_NONE , X, T2);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_add(1.0,T2,1.0,T);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    ret = mess_matrix_multiply(MESS_OP_NONE, X, op, H, T2);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_add(1.0,T2,1.0,T);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    ret = mess_matrix_norm2(T, &nrm);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm2);
    *res = nrm;

    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
    mess_matrix_clear(&T);
    mess_matrix_clear(&T2);

    return 0;
}


/**
 * @brief Compute the \f$ 2 \f$-norm residual of a semi-generalized Sylvester equation.
 * @param[in] op  input operation on the matrices \f$ A, E \f$ and \f$ H \f$
 * @param[in] A  input first matrix
 * @param[in] E  input second matrix
 * @param[in] H  input third matrix
 * @param[in] M  input right hand side matrix
 * @param[in] X  input solution matrix
 * @param[out] res  norm
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_sylvestersg_res2 function computes the 2-Norm residual of the semi-generalized Sylvester equation
 * \f[  R:= op(A)X+op(E)Xop(H)+M = 0 \f]
 * that means it computes
 * \f[ \Vert op(A)X + op(E)Xop(H) + M \Vert_2 , \f]
 * where \f$ op \f$ can be
 * <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 *  |@ref MESS_OP_HERMITIAN      |   \f$ op(A)=A^H\f$         |
 * </center>
 * It is implemented as dense computation. The matrix \f$ R \f$ is computed
 * inside the function. If the input matrices are too large, this function
 * will be time consuming.
 *
 */
int mess_direct_sylvestersg_res2 (mess_operation_t op,  mess_matrix A, mess_matrix E, mess_matrix H, mess_matrix M, mess_matrix X, double *res )
{
    MSG_FNAME(__func__);
    mess_matrix T;
    mess_matrix T2;
    mess_matrix T3;
    int ret = 0;
    double nrm = 0;

    /*-----------------------------------------------------------------------------
     *  call standard residual routine if E points to NULL
     *-----------------------------------------------------------------------------*/
    if(!E){
        return mess_direct_sylvester_res2( op, A, H, M, X, res);
    }


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(E);
    mess_check_nullpointer(H);
    mess_check_nullpointer(M);
    mess_check_nullpointer(X);
    mess_check_nullpointer(res);
    mess_check_same_size(A,E);
    mess_check_square(A);
    mess_check_square(H);
    mess_check_same_size(M,X);
    mess_check_same_rows(A,M);
    mess_check_same_cols(H,M);


    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&T);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0) ,mess_matrix_init);
    ret = mess_matrix_init(&T2);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0) ,mess_matrix_init);
    ret = mess_matrix_init(&T3);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0) ,mess_matrix_init);
    ret = mess_matrix_copy(M, T);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);

    ret = mess_matrix_multiply(op, A, MESS_OP_NONE, X, T2); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
    ret = mess_matrix_add(1.0,T2,1.0,T);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);

    ret = mess_matrix_multiply(op, E, MESS_OP_NONE, X, T3); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, T3, op, H, T2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_add(1.0,T2,1.0,T);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);

    ret = mess_matrix_norm2(T, &nrm); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm2);
    *res = nrm;

    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
    mess_matrix_clear(&T);
    mess_matrix_clear(&T2);
    mess_matrix_clear(&T3);

    return 0;
}

/**
 * @brief Compute the \f$ 2 \f$-norm residual of a generalized Sylvester equation.
 * @param[in] op  input operation on the matrices \f$ A, F, E \f$ and \f$ H \f$
 * @param[in] A  input first matrix
 * @param[in] F  input second matrix
 * @param[in] E  input third matrix
 * @param[in] H  input fourth matrix
 * @param[in] M  input right hand side matrix
 * @param[in] X  input solution matrix
 * @param[out] res  norm
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_generalized_sylvester_res2 function computes the 2-Norm residual of the generalized Sylvester equation
 * \f[ R:=op(A) X op(F)+op(E) X op(H)+M = 0 \f]
 * that means it computes
 * \f[ \Vert op(A) Xop(F)+op(E)Xop(H)+M \Vert_2, \f]
 * where \f$ op \f$ can be
 * <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 *  |@ref MESS_OP_HERMITIAN      |   \f$ op(A)=A^H\f$         |
 * </center>
 * It is implemented as dense computation. The matrix \f$ R \f$ is computed
 * inside the function. If the input matrices are too large, this function
 * will be time consuming.
 *
 */
int mess_direct_generalized_sylvester_res2 ( mess_operation_t op,  mess_matrix A, mess_matrix F, mess_matrix E, mess_matrix H, mess_matrix M, mess_matrix X, double *res )
{
    MSG_FNAME(__func__);
    mess_matrix T;
    mess_matrix T2;
    mess_matrix T3;
    int ret = 0;
    double nrm = 0;

    /*-----------------------------------------------------------------------------
     *  call standard residual routine if E points to NULL
     *-----------------------------------------------------------------------------*/
    if(!E && !F){
        return mess_direct_sylvester_res2( op, A, H, M, X, res);
    }

    if(E && !F){
        return mess_direct_sylvestersg_res2(op, A, E, H, M, X, res);
    }


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(E);
    mess_check_nullpointer(H);
    mess_check_nullpointer(F);
    mess_check_nullpointer(M);
    mess_check_nullpointer(X);
    mess_check_nullpointer(res);
    mess_check_same_size(A,E);
    mess_check_same_size(H,F);
    mess_check_square(A);
    mess_check_square(H);
    mess_check_same_size(M,X);
    mess_check_same_rows(A,M);
    mess_check_same_cols(H,M);

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&T);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0) ,mess_matrix_init);
    ret = mess_matrix_init(&T2);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0) ,mess_matrix_init);
    ret = mess_matrix_init(&T3);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0) ,mess_matrix_init);
    ret = mess_matrix_copy(M, T);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);

    ret = mess_matrix_multiply(op, A, MESS_OP_NONE, X, T2); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, T2, op, F, T3); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_add(1.0,T3,1.0,T);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);


    ret = mess_matrix_multiply(op, E, MESS_OP_NONE, X, T3); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, T3, op, H, T2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_add(1.0,T2,1.0,T);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);

    ret = mess_matrix_norm2(T, &nrm); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm2);
    *res = nrm;

    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
    mess_matrix_clear(&T);
    mess_matrix_clear(&T2);
    mess_matrix_clear(&T3);

    return 0;
}
