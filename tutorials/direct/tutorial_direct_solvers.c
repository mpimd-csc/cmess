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
 * @file tutorials/direct/tutorial_direct_solvers.c
 * @brief Tutorial about linear equation systems.
 * @author  @mbehr
 *
 * # Tutorial About Linear Equation Systems
 * In this Tutorial we present functions for solving linear equation systems.
 *
 * ## 1. Include Header Files
 * We include standard @c C header files and in addtion the @c mess.h header file
 * to make the @mess library available.
 * @snippet "tutorials/direct/tutorial_direct_solvers.c" HEADER
 * \n
 *
 * ## 2. Generate Sparse Random Matrix
 * At first we want to generate a sparse diagonal dominant matrix \f$A\f$, to make sure that
 * \f$ A \f$ has full rank. We perform the following steps
 * \f[A\leftarrow I\in \mathbb{R}^{n\times n },\f]
 * \f[A\leftarrow (n+1)\cdot A + rand \in \mathbb{R}^{n\times n}. \f]
 * @snippet "tutorials/direct/tutorial_direct_solvers.c" SPARSE
 * \n
 *
 * ## 3. Solve Sparse Linear Systems
 * In this section we show different functions to solve a sparse linear system.
 *
 * ### 3.1 Solve Sparse Linear Systems with mess_matrix_backslash
 * The simplest way to solve a linear equation system is to use @ref mess_matrix_backslash.
 * We generate a random right hand side \f$b\in \mathbb{R}^{n}\f$ and solve \f$ Ax_0=b\f$ and \f$ A^Tx_0=b\f$.
 * After that we want to check our solution \f$x\f$ and compute the residual
 * \f[y\leftarrow Ax,\f]
 * \f[y\leftarrow A^Tx,\f]
 * \f[res \leftarrow \Vert y-b \Vert_2.\f]
 * We also measure the wall-clock-time.
 * @snippet "tutorials/direct/tutorial_direct_solvers.c" SPARSESOLVE
 * \n
 *
 * ### 3.2 Solve Sparse Linear Systems with mess_direct structure
 * @ref mess_matrix_backslash and @ref mess_matrix_backslashm recompute the decomposition of \f$A\f$ in every
 * call. If you want to solve a lot of systems with different right hand sides you should use @ref mess_direct structure.
 *
 * @snippet "tutorials/direct/tutorial_direct_solvers.c" SPARSEDIRECT
 * At this point please compare the wall-clock-time between @ref mess_direct_lu and and @ref mess_matrix_backslash.
 * \n
 *
 * ### 3.3 Change LU solver
 * You can change the direct solver with @ref mess_direct_lu_select. The residual and relative residual
 * is easily computed by @ref mess_direct_res2.
 * @snippet "tutorials/direct/tutorial_direct_solvers.c" SPARSESELECT
 * \n
 *
 * ## 4. Symmetric Positive Definite System
 * Remember that \f$ A \f$ is diagonal dominant and all diagonal entries are real and positive.
 * We compute
 * \f[ A\leftarrow A+A^{T} \f]
 * and obtain a symmetric positive definite matrix.
 * Generating a cholesky decomposition based direct solver is similiar to the previous part.
 * @snippet "tutorials/direct/tutorial_direct_solvers.c" SPARSEPOS
 * \n
 *
 * ## 5. Dense Systems
 * In the case of a dense matrix, solving linear system works almost the same in @mess.
 *
 * @snippet "tutorials/direct/tutorial_direct_solvers.c" DENSE
 * \n
 * If you want to perform a QR decomposition feel free to take a look at @ref mess_direct_create_lapack_qr.
 *
 * ## 6. Clear Memory
 *
 * @snippet "tutorials/direct/tutorial_direct_solvers.c" CLEAR
 * \n
 *
 * @sa @ref tutorials
 */

///@cond
///[HEADER]
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "mess/mess.h"
///[HEADER]

int main ( int argc, char **argv){

    ///[SPARSE]
    mess_int_t n = 1000;
    double fill_rate = 0.05;
    mess_matrix A,rand;
    mess_matrix_init(&A);
    mess_matrix_init(&rand);
    mess_matrix_eye(A,n,n,MESS_CSR);
    mess_matrix_rand(rand,n,n,MESS_CSR,MESS_REAL,fill_rate);
    mess_matrix_add(1.0,rand,n+1,A);
    ///[SPARSE]

    ///[SPARSESOLVE]
    double res, walltime;
    mess_vector b, x, y;
    MESS_INIT_VECTORS(&b,&x,&y);
    mess_vector_alloc(b,n,MESS_REAL);
    mess_vector_alloc(x,n,MESS_REAL);
    mess_vector_alloc(y,n,MESS_REAL);
    mess_vector_rand(b);

    walltime = mess_wtime();
    mess_matrix_backslash(MESS_OP_NONE,A,b,x);
    mess_matrix_mvp(MESS_OP_NONE,A,x,y);
    mess_vector_diffnorm(y,b,&res);
    printf("||A*x-b||_2 \t= %10.5e\n",res);

    mess_matrix_backslash(MESS_OP_TRANSPOSE,A,b,x);
    mess_matrix_mvp(MESS_OP_TRANSPOSE,A,x,y);
    mess_vector_diffnorm(y,b,&res);
    printf("||A^T*x-b||_2 \t= %10.5e\n",res);

    printf("Wall Time (mess_matrix_backslash) = %10.5f\n\n",mess_wtime()-walltime);
    ///[SPARSESOLVE]

    ///[SPARSEDIRECT]
    walltime = mess_wtime();

    mess_direct solver;
    mess_direct_init(&solver);
    mess_direct_lu(A,solver);

    mess_direct_solve(MESS_OP_NONE,solver,b,x);
    mess_matrix_mvp(MESS_OP_NONE,A,x,y);
    mess_vector_diffnorm(y,b,&res);
    printf("||A*x-b||_2 \t= %10.5e\n",res);

    mess_direct_solve(MESS_OP_TRANSPOSE,solver,b,x);
    mess_matrix_mvp(MESS_OP_TRANSPOSE,A,x,y);
    mess_vector_diffnorm(y,b,&res);
    printf("||A^T*x-b||_2 \t= %10.5e\n",res);

    mess_direct_clear(&solver);
    printf("Wall Time (mess_matrix_backslash) = %10.5f\n\n",mess_wtime()-walltime);
    ///[SPARSEDIRECT]

    ///[SPARSESELECT]
    double rel_res;

    mess_direct_lu_select(MESS_DIRECT_SPARSE_LU);
    mess_direct_init(&solver);
    mess_direct_lu(A,solver);
    mess_direct_solve(MESS_OP_NONE,solver,b,x);
    mess_direct_res2(MESS_OP_NONE,A,x,b,&res,&rel_res);
    printf("MESS_DIRECT_SPARSE_LU\n");
    printf("||A*x-b||_2 \t= %10.5e\n",res);
    printf("||A*x-b||_2/||b||_2 \t= %10.5e\n",rel_res);
    printf("\n");
    mess_direct_clear(&solver);


    mess_direct_lu_select(MESS_DIRECT_UMFPACK_LU);
    mess_direct_init(&solver);
    mess_direct_lu(A,solver);
    mess_direct_solve(MESS_OP_NONE,solver,b,x);
    mess_direct_res2(MESS_OP_NONE,A,x,b,&res,&rel_res);
    printf("MESS_DIRECT_UMFPACK_LU\n");
    printf("||A*x-b||_2 \t= %10.5e\n",res);
    printf("||A*x-b||_2/||b||_2 \t= %10.5e\n",rel_res);
    printf("\n");
    mess_direct_clear(&solver);


    mess_direct_lu_select(MESS_DIRECT_CSPARSE_LU);
    mess_direct_init(&solver);
    mess_direct_lu(A,solver);
    mess_direct_solve(MESS_OP_NONE,solver,b,x);
    mess_direct_res2(MESS_OP_NONE,A,x,b,&res,&rel_res);
    printf("MESS_DIRECT_CSPARSE_LU\n");
    printf("||A*x-b||_2 \t= %10.5e\n",res);
    printf("||A*x-b||_2/||b||_2 \t= %10.5e\n",rel_res);
    printf("\n");
    mess_direct_clear(&solver);
    ///[SPARSESELECT]

    ///[SPARSEPOS]
    mess_matrix_ctranspose(A,rand);
    mess_matrix_add(1.0,rand,1.0,A);

    mess_direct_chol_select(MESS_DIRECT_CSPARSE_CHOLESKY);
    mess_direct_init(&solver);
    mess_direct_chol(A,solver);
    mess_direct_solve(MESS_OP_NONE,solver,b,x);
    mess_direct_res2(MESS_OP_NONE,A,x,b,&res,&rel_res);
    printf("MESS_DIRECT_CSPARSE_CHOLESKY\n");
    printf("||A*x-b||_2 \t= %10.5e\n",res);
    printf("||A*x-b||_2/||b||_2 \t= %10.5e\n",rel_res);
    printf("\n");
    mess_direct_clear(&solver);
    ///[SPARSEPOS]


    ///[DENSE]
    mess_int_t size2 = 50;
    mess_matrix_rand(A,size2,size2,MESS_DENSE,MESS_REAL,0.0);

    mess_direct_lu_select(MESS_DIRECT_LAPACK_LU);
    mess_direct_init(&solver);
    mess_direct_lu(A,solver);

    mess_vector_resize(b,size2);
    mess_vector_resize(x,size2);
    mess_direct_solve(MESS_OP_NONE,solver,b,x);
    mess_direct_res2(MESS_OP_NONE,A,x,b,&res,&rel_res);
    printf("MESS_DIRECT_LAPACK_LU\n");
    printf("||A*x-b||_2 \t= %10.5e\n",res);
    printf("||A*x-b||_2/||b||_2 \t= %10.5e\n",rel_res);
    printf("\n");
    mess_direct_clear(&solver);
    ///[DENSE]

    ///[CLEAR]
    mess_matrix_clear(&A);
    mess_matrix_clear(&rand);
    mess_vector_clear(&b);
    mess_vector_clear(&x);
    mess_vector_clear(&y);
    ///[CLEAR]

    return 0;
}
///@endcond

