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
 * @file lib/direct/care.c
 * @brief Solve the algebraic Riccati Equation.
 * @author @koehlerm
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>
#include <math.h>


/**
 * @brief Compute the solution of a (generalized) Riccati Equation using a dense newton method.
 * @param[in] A     input \f$ A \f$ matrix of the Riccati Equation
 * @param[in] E     input \f$ E \f$ matrix of the Riccati Equation (NULL if \f$ E \f$ does not exist)
 * @param[in] Q     input \f$ Q \f$ matrix of the Riccati Equation
 * @param[in] G     input \f$ G \f$ matrix of the Riccati Equation
 * @param[out] X    output \f$ X \f$ solution of the Riccati Equation
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_direct_care function solves either
 * \f[
 *  A^T X + XA - XGX + Q = 0
 * \f]
 * or
 * \f[
 *  A^T X E+E^T X A - E^T XGX E + Q = 0
 * \f]
 * depending on the existence of the matrix \f$ E \f$. \n
 * If the matrix object for \f$ E  \f$ is set to @c NULL the standard Riccati Equation is solved, otherwise the
 * generalized one is solved.
 *
 * The dense Newton method from \ref mess_dense_nm_gmare  is used.
 *
 * \attention
 * \li If the inputs are sparse matrices, they are converted to dense ones. This might cause an out of memory.
 */
int mess_direct_care (mess_matrix A, mess_matrix E, mess_matrix Q, mess_matrix G, mess_matrix X ){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_matrix sA, sE, sG, sQ;
    mess_int_t convA=-1, convE=-1, convG=-1, convQ=-1;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(Q);
    mess_check_nullpointer(G);
    mess_check_nullpointer(X);
    mess_check_square(A);
    mess_check_square(Q);
    mess_check_square(G);
    mess_check_real(A);
    mess_check_real(G);
    mess_check_real(Q);
    mess_check_same_size(A,Q);
    mess_check_same_size(A,G);

    if(E){
        mess_check_same_size(A,E);
    }


    /*-----------------------------------------------------------------------------
     *  convert matrices
     *-----------------------------------------------------------------------------*/
    MESS_MATRIX_CHECKFORMAT(A,sA,convA,MESS_DENSE);         FUNCTION_FAILURE_HANDLE(convA,(convA>0),mess_mess_matrix_convert);
    MESS_MATRIX_CHECKFORMAT(G,sG,convG,MESS_DENSE);         FUNCTION_FAILURE_HANDLE(convG,(convG>0),mess_mess_matrix_convert);
    MESS_MATRIX_CHECKFORMAT(Q,sQ,convQ,MESS_DENSE);         FUNCTION_FAILURE_HANDLE(convQ,(convQ>0),mess_mess_matrix_convert);

    if(E){
        MESS_MATRIX_CHECKFORMAT(E,sE,convE,MESS_DENSE);     FUNCTION_FAILURE_HANDLE(convE,(convE>0),mess_matrix_convert);
    }

    /*-----------------------------------------------------------------------------
     *  Solve using Netwon Method
     *-----------------------------------------------------------------------------*/
    ret = mess_dense_nm_gmare(sA, E? sE : NULL, sQ, sG, X);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_dense_nm_gmare );

    /*-----------------------------------------------------------------------------
     *  clear memory if matrix was converted
     *-----------------------------------------------------------------------------*/
    if (convA==0) mess_matrix_clear(&sA);
    if (convG==0) mess_matrix_clear(&sG);
    if (convQ==0) mess_matrix_clear(&sG);
    if (convE==0) mess_matrix_clear(&sE);

    return ret;
}       /* -----  end of function mess_direct_care  ----- */



/**
 * @brief Compute the \f$ 2 \f$-norm residual of the algebraic (generalized) Riccati equation.
 * @param[in] A         input matrix
 * @param[in] E         input matrix of the Riccati Equation (NULL if \f$ E \f$ does not exist)
 * @param[in] Q         input matrix
 * @param[in] G         input matrix
 * @param[in] X         input solution of the Riccati Equation
 * @param[out] res2     output \f$ 2 \f$-norm residual
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_care_res2 function computes the \f$ 2 \f$-norm residual of the dense Riccati Equation:
 * \f[ res2= \Vert Q + A^T X+XA -XGX \Vert_2 . \f]
 *
 * If @p E is not @c NULL then we compute the \f$ 2 \f$-norm residual of the generalized Riccati Equation
 *
 * \f[ res2=\Vert Q + A^T XE+E^T XA-E^T XGXE \Vert_2. \f]
 *
 *
 *
 */
int mess_direct_care_res2 ( mess_matrix A, mess_matrix E,  mess_matrix Q, mess_matrix G, mess_matrix X, double *res2 )
{
    MSG_FNAME(__func__);
    int ret =0;

    ret = mess_dense_res_gmpare( A, E, Q, G, X, 0, MESS_OP_TRANSPOSE, MESS_2_NORM, res2, NULL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_pcare_res);

    return ret;
}       /* -----  end of function mess_direct_care_res2  ----- */





