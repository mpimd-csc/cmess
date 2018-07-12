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
 * @file lib/direct/backslash.c
 * @brief Implementation of a @matlab  backslash like function call.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#ifdef _OPENMP_H
#include <omp.h>
#endif

/**
 * @brief Solve \f$ op(A)x=b \f$  with \f$ b \f$ and \f$ x \f$ vectors.
 * @param[in] op    input operation on \f$ A \f$
 * @param[in] A     input system matrix  \f$A\f$
 * @param[in] b     input right hand side \f$ b\f$
 * @param[out] x    output solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_backslash function solves the system \f$ op(A)x=b \f$. \n
 * The solver (sparse, dense) is detected automatically.\n
 * If you plan to use the matrix \f$ A \f$ more than one time to solve
 * system please use a @ref mess_direct instance. \n
 * The operation \f$ op \f$ can be one of
 * \li \ref MESS_OP_NONE, that means \f$ op(A)=A \f$
 * \li \ref MESS_OP_TRANSPOSE, that means \f$ op(A)=A^T \f$
 * \li \ref MESS_OP_HERMITIAN, that means \f$ op(A)=A^H. \f$
 *
 * \attention
 * This function only works if \f$ b \f$ and \f$ x \f$ are vectors. If the right hand side
 * is a matrix please use \ref mess_matrix_backslashm.
 *
 */
int mess_matrix_backslash(mess_operation_t op, mess_matrix A, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    mess_direct solver;
    int ret = 0;

    mess_check_nullpointer(A);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    ret = mess_direct_init(&solver);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_direct_init);
    if ( A->rows != A->cols ) {
        if ( MESS_IS_DENSE(A)) {
            ret = mess_direct_create_lapack_qr(A,solver);   FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_direct_create_lapack_qr);
        } else {
            MSG_ERROR("solver not available for sparse matrix\n");
            return MESS_ERROR_NOSUPPORT;
        }
    } else {
        if ( MESS_IS_DENSE(A)) {
            ret = mess_direct_create_lapack_lu(A, solver);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_lapack_lu);
        } else {
            ret = mess_direct_lu(A, solver);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);
        }
    }
    ret = mess_direct_solve(op,solver,b,x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solve);
    mess_direct_clear(&solver);
    return 0;
}


/**
 * @brief Solve \f$ op(A)x=b \f$ with \f$ b \f$ and \f$ x \f$ matrices.
 * @param[in] op    input operation on \f$ A \f$
 * @param[in] A     input system matrix  \f$ A \f$
 * @param[in] b     input right hand side \f$ b\f$
 * @param[out] x    output solution \f$x\f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_backslashm function solves the system \f$ op(A)x=b \f$. \n
 * The solver (sparse, dense) is detected automatically. \n
 * If you plan to use the matrix \f$ A \f$ more than one time to solve
 * system please use a @ref mess_direct instance. \n
 * The operation \f$ op \f$ can be one of
 * \li \ref MESS_OP_NONE, that means \f$ op(A)=A \f$
 * \li \ref MESS_OP_TRANSPOSE, that means \f$ op(A)=A^T \f$
 * \li \ref MESS_OP_HERMITIAN, that means \f$ op(A)=A^H. \f$
 *
 * \attention
 * This function only works if \f$ b \f$ and \f$ x \f$ are matrices. If the right hand side
 * is a vector please use \ref mess_matrix_backslash.
 *
 */
int mess_matrix_backslashm(mess_operation_t op, mess_matrix A, mess_matrix b, mess_matrix x) {
    MSG_FNAME(__func__);
    mess_direct solver;
    int ret = 0;

    mess_check_nullpointer(A);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    ret = mess_direct_init(&solver);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
    if ( A->rows != A->cols ) {
        if ( MESS_IS_DENSE(A)) {
            ret = mess_direct_create_lapack_qr(A,solver);   FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_direct_create_lapack_qr);
        } else {
            MSG_ERROR("solver not available for sparse\n");
            return MESS_ERROR_NOSUPPORT;
        }
    } else {
        if ( MESS_IS_DENSE(A)) {
            ret = mess_direct_create_lapack_lu(A, solver);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_lapack_lu);
        } else {
            ret = mess_direct_lu(A, solver);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);
        }
    }
    ret = mess_direct_solvem(op,solver,b,x);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
    mess_direct_clear(&solver);
    return 0;
}
