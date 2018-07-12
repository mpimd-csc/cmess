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
 *
 * @file tutorials/direct/tutorial_sylvester.c
 * @brief Demonstrate how to solve a sparse-dense Sylvester equation.
 * @author @dykstra
 *
 *
 * # Tutorial Solving a (generalized) sparse-dense Sylvester Equation
 * This tutorial demonstrates how to solve sparse-dense Sylvester equations
 * \f[ op(A)X[op(F)] + [op(E)]Xop(H) + M = 0 , \f]
 * where \f$ A \f$ and \f$ E \f$ are large and sparse matrices, \f$ F \f$ and \f$ H \f$ are small and dense matrices and \f$ M \f$ is a matrix,
 * by generating a solver and using the \ref mess_direct_solvem  function. \n
 * \f$ E \f$ and \f$ F \f$ may or may not be used, there are solvers for three cases (further info below). \n
 * Operations \f$ op(.) \f$ on matrices \f$ A \f$ and \f$ H \f$ can be
 * <ul>
 * <li> \ref MESS_OP_NONE (\f$ op(C) = C \f$),
 * <li> \ref MESS_OP_TRANSPOSE (\f$ op(C) = C^T \f$),
 * <li> \ref MESS_OP_HERMITIAN (\f$ op(C)= C^H \f$).
 * </ul>
 * In this program we have to differentiate between the three cases of sparse-dense Sylvester equations the \ref mess_direct_create_sylvester_sparsedense
 * function can generate solvers for :\n
 * <ul>
 * <li> \f$ op(A)X+Xop(H)+M=0 \f$
 * <li> \f$ op(A)X+op(E)Xop(H)+M=0 \f$
 * <li> \f$ op(A)Xop(F)+op(E)Xop(H)+M=0 \f$
 * </ul>
 * After solving for X we compute the relative residual
 * \f$ \|op(A)X[op(F)] + [op(E)]Xop(H) + M\|_2/ \|\| M \|_2 \f$
 * and clear the memory.
 *
 * ## 1. Include Header Files
 * We include standard @c C header files and in addtion the @c mess.h header file
 * to make the @mess library available and the @c complex.h header file for handling of complex matrices.
 * @snippet "tutorials/direct/tutorial_sylvester.c" HEADER
 * \n
 *
 * ## 2. Print Title (optional)
 * Print the M.E.S.S. version and the title of the demo.
 * @snippet "tutorials/direct/tutorial_sylvester.c" TITLE
 * \n
 *
 * ## 3. Declare Matrices and Check the Number of Input Arguments.
 * We use the \ref mess_matrix structure to declare our matrices \f$ A, F, E, H, M\f$ and \f$ X1, X2\f$ for operations \ref MESS_OP_NONE
 * and \ref MESS_OP_HERMITIAN. We also declare a \ref mess_direct structure to store the generated Sylvester solver.
 * @snippet "tutorials/direct/tutorial_sylvester.c" DECLARE
 * By using the \c argc variable we can check the number of input arguments.
 * @snippet "tutorials/direct/tutorial_sylvester.c" INPUT
 * \n
 *
 * ## 4. Init Matrices and Read Data
 * Before we are able to use we have to init the fields of the \ref mess_matrix structure variables using \ref mess_matrix_init.
 * Be careful to catch the different input cases.\n
 * @snippet "tutorials/direct/tutorial_sylvester.c" INIT
 * \n
 *
 * Now we can use \ref mess_matrix_read or \ref mess_matrix_read_formated to read our matrices from a Matrix Market File Format and convert \f$ A \f$ and \f$ E\f$ to
 * the storage type \ref MESS_CSR if needed. Also differentiate the input cases here.\n
 * @snippet "tutorials/direct/tutorial_sylvester.c" READ
 * \n
 *
 * ## 5. Generate the Sylvester Solver and Solve
 * As already announced \ref mess_direct_create_sylvester_sparsedense creates a solver for the Sylvester equation.
 * We can use this solver with \ref mess_direct_solvem to solve the (generalized) Sylvester Equation. \n
 * We measure the time using \ref mess_wtime. \n
 * The solver is initialized using \ref mess_direct_init and created using \ref mess_direct_create_sylvester_sparsedense.
 * We can indicate the case of the standard Sylvester Equation by providing  @c NULL as  the second  and third argument to
 * \ref mess_direct_create_sylvester_sparsedense, as well as the semi-generalized case \f$ op(A)X+op(E)Xop(H)+M=0 \f$ by
 * providing NULL as second argument.\n
 * @snippet "tutorials/direct/tutorial_sylvester.c" SOLVE
 * \n
*
* ## 6. Compute the Relative Residual
* \ref mess_matrix_norm2 computes the norm of M, which we divide from \ref mess_direct_sylvester_res2, \ref mess_direct_sylvestersg_res2
* or \ref mess_direct_generalized_sylvester_res2 to get the absolute residual. The relative residual residual can easily obtained by computing the norm of the right hand side.
* @snippet "tutorials/direct/tutorial_sylvester.c" RES
* \n
*
* ## 7. Don't Forget to Clear Memory
* A rule of thumb is that every structure which needs an @c init call before usage must be cleared at the
* at the end of your code. In our cases we have used \ref mess_matrix_init which means we have to use
* \ref mess_matrix_clear to clear dynamical allocated memory (we already cleared the \ref mess_direct instance.
*
* @snippet "tutorials/direct/tutorial_sylvester.c" CLEAR
* \n
*
* For a special dense Sylvester solve take a look at  @ref mess_direct_create_sylvester_dense.
*
* @sa @ref mess_direct_create_sylvester_sparsedense
* @sa @ref mess_direct_create_sylvester_sparsedense
* @sa @ref mess_direct_solvem
*/

///@cond
///[HEADER]
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/mess.h"
#include <complex.h>
///[HEADER]


int main ( int argc, char **argv){
    ///[TITLE]
    mess_version();
    printf("mess sylvester solver demo\n");
    printf("==========================\n\n");
    ///[TITLE]

    ///[DECLARE]
    mess_matrix A,F=NULL,E=NULL,H,M,X1,X2;
    mess_direct sylv;
    ///[DECLARE]

    ///[INPUT]
    if ( !( 3 <= argc && argc <= 5)){
        printf("usage: %s A.mtx [F.mtx] [E.mtx] H.mtx\n", argv[0]);
        return -1;
    }
    ///[INPUT]

    ///[INIT]
    mess_matrix_init(&A);
    mess_matrix_init(&F);
    mess_matrix_init(&E);
    mess_matrix_init(&H);
    mess_matrix_init(&M);
    mess_matrix_init(&X1);
    mess_matrix_init(&X2);
    ///[INIT]

    ///[READ]
    mess_matrix_read_formated(argv[1], A, MESS_CSR);
    if( argc == 3 ){
        mess_matrix_read_formated(argv[2], H, MESS_DENSE);
    } else if( argc == 4 ){
        mess_matrix_read_formated(argv[2], E, MESS_CSR);
        mess_matrix_read(argv[3], H);
    } else {
        mess_matrix_read_formated(argv[2], F, MESS_DENSE);
        mess_matrix_read_formated(argv[3], E, MESS_CSR);
        mess_matrix_read_formated(argv[4], H, MESS_DENSE);
    }
    mess_int_t seed = 1234;
    mess_matrix_rand_init(&seed);
    mess_matrix_rand_dense(M, A->rows, H->rows, MESS_REAL);
    ///[READ]

    ///[SOLVE]
    double ts, te;
    ts = mess_wtime();
    mess_direct_init(&sylv);
    mess_direct_create_sylvester_sparsedense(A,NULL,NULL,H,sylv);
    mess_direct_solvem(MESS_OP_NONE,sylv,M,X1);
    mess_direct_solvem(MESS_OP_HERMITIAN,sylv,M,X2);
    mess_direct_clear(&sylv);
    te = mess_wtime();
    printf("time: %.3f\n\n", te-ts);
    ///[SOLVE]

    ///[RES]
    double res = 0, res0 = 0;
    mess_matrix_norm2(M, &res0);
    if(argc == 3){
        mess_direct_sylvester_res2(MESS_OP_NONE,A,H,M,X1,&res);
        printf("||AX+XH+M||_2 / ||M||_2 = %.10e\n", res/res0);
        mess_direct_sylvester_res2(MESS_OP_HERMITIAN,A,H,M,X2,&res);
        printf("||A^HX+XH^H+M||_2 / ||M||_2 = %.10e\n", res/res0);
    } else if (argc == 4){
        mess_direct_sylvestersg_res2(MESS_OP_NONE,A,E,H,M,X1,&res);
        printf("||AX+EXH+M||_2 / ||M||_2 = %.10e\n", res/res0);
        mess_direct_sylvestersg_res2(MESS_OP_HERMITIAN,A,E,H,M,X2,&res);
        printf("||A^HX+E^HXH^H+M||_2 / ||M||_2 = %.10e\n", res/res0);
    } else {
        mess_direct_generalized_sylvester_res2(MESS_OP_NONE,A,F,E,H,M,X1,&res);
        printf("||AXF+EXH+M||_2 || / ||M||_2 = %.10e\n", res/res0);
        mess_direct_generalized_sylvester_res2(MESS_OP_HERMITIAN,A,F,E,H,M,X2,&res);
        printf("||A^HXF^H+E^HXH^H+M||_2 / ||M||_2 = %.10e\n", res/res0);
    }
    ///[RES]

    ///[CLEAR]
    mess_matrix_init(&A);
    if(F)mess_matrix_init(&F);
    if(E)mess_matrix_init(&E);
    mess_matrix_init(&H);
    mess_matrix_init(&M);
    mess_matrix_init(&X1);
    mess_matrix_init(&X2);
    ///[CLEAR]

    return 0;
}
///@endcond
