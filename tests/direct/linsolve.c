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
 * @addtogroup test_direct
 * @{
 * @file tests/direct/linsolve.c
 * @brief Check the solution of a linear system of equations.
 * @author @koehlerm
 * @author @mbehr
 * @test
 *
 * This function checks the @ref  mess_direct_solve function defined in direct.c for a linear system of equations that means
 * it checks if the solution of
 * \f[op(A) x= b ,\f]
 * is computed correctly for a given matrix \f$ A \f$ and a given vector \f$ b \f$. \n
 * Operations \f$ op (.)\f$ on matrix \f$ A \f$ can be
 * <ul>
 * <li> \ref MESS_OP_NONE (\f$ op(A)= A \f$),
 * <li> \ref MESS_OP_TRANSPOSE (\f$ op(A) = A^T \f$),
 * <li> \ref MESS_OP_HERMITIAN (\f$ op(A) = A^H \f$).
 * </ul>
 * The solver can be selected by the third argument for the second input:
 *
 *
 * |   Input    |           Description                                                 |
 * |:----------:|:---------------------------------------------------------------------:|
 * |    0       | mess_direct_create_sparse_lu                                          |
 * |    1       |  mess_direct_create_csparse_lu                                        |
 * |    2       |  mess_direct_create_umfpack                                           |
 * |    3       |  mess_direct_create_bicgstab                                          |
 * |    4       |  mess_direct_create_lapack_lu                                         |
 * |    5       |  mess_direct_create_lapack_qr                                         |
 * |    6       |  mess_direct_create_cholesky                                          |
 * |    7       |  mess_direct_create_csparse_cholesky                                  |
 * |    8       |  mess_direct_create_cholmod_cholesky                                  |
 * |    9       |  mess_direct_create_banded                                            |
 * |   10       |  mess_direct_create_superlu                                           |
 * |   11       |  mess_direct_create_mklpardiso                                        |
 *
 *
 *
 * @}
 *
 */
#include "../call_macro.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"



/**
 * @brief Select a solver.
 * @param[in] i      input solver to choose
 * @param[in] A          input matrix
 * @param[out] sol  generated solver
 * @return zero on success or a non-zero error value otherwise
 *
 * The  select_solver function generates a solver for a linear system of equations
 * \f[A x= b \f]
 * Solvers are choosen with \f$ i \f$. Possible values are:
 *
 *
 * |   Input    |           Description                                                 |
 * |:----------:|:---------------------------------------------------------------------:|
 * |    0       |  mess_direct_create_sparse_lu                                         |
 * |    1       |  mess_direct_create_csparse_lu                                        |
 * |    2       |  mess_direct_create_umfpack                                           |
 * |    3       |  mess_direct_create_bicgstab                                          |
 * |    4       |  mess_direct_create_lapack_lu                                         |
 * |    5       |  mess_direct_create_lapack_qr                                         |
 * |    6       |  mess_direct_create_cholesky                                          |
 * |    7       |  mess_direct_create_csparse_cholesky                                  |
 * |    8       |  mess_direct_create_cholmod_cholesky                                  |
 * |    9       |  mess_direct_create_banded                                            |
 * |   10       |  mess_direct_create_superlu                                           |
 * |   11       |  mess_direct_create_mklpardiso                                        |
 *
 *
 */
int select_solver(int i, mess_matrix A, mess_direct sol) {
    int ret =0;
    switch(i){
        case 0:
            CALL(mess_direct_create_sparse_lu(A,sol));
            break;
#ifdef MESS_HAVE_CSPARSE
        case 1:
            CALL(mess_direct_create_csparse_lu(A, sol));
            break;
#endif
#ifdef MESS_HAVE_UMFPACK
        case 2:
            CALL(mess_direct_create_umfpack(A,sol));
            break;
#endif
        case 3:
            CALL(mess_direct_create_bicgstab(A,sol));
            break;
        case 4:
            CALL(mess_direct_create_lapack_lu(A,sol));
            break;
        case 5:
            CALL(mess_direct_create_lapack_qr(A,sol));
            break;
        case 6:
            CALL(mess_direct_create_cholesky(A,sol));
            break;
#ifdef MESS_HAVE_CSPARSE
        case 7:
            CALL(mess_direct_create_csparse_cholesky(A,sol));
            break;
#endif
#ifdef MESS_HAVE_CHOLMOD
        case 8:
            CALL(mess_direct_create_cholmod_cholesky(A,sol));
            break;
#endif
        case 9:
            CALL(mess_direct_create_banded(A,sol));
            break;
#ifdef MESS_HAVE_SUPERLU
        case 10:
            CALL(mess_direct_create_superlu(A,sol));
            break;
#endif

#ifdef MESS_HAVE_MKLPARDISO
        case 11:
            CALL(mess_direct_create_mklpardiso(A,sol));
            break;

#endif
        default:
            return 1;
    }
    return 0;
}
/**
 * @brief Test solving a linear system of equations.
 * @param[in] i      input solver to choose
 * @param[in] A          input matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The  testmat function solves a linear system of equations
 * \f[op(A) x= b \f]
 * and checks if the solution is computed correctly.\n
 * Operations \f$ op (.)\f$ on matrix \f$ A \f$ can be
 * <ul>
 * <li> \ref MESS_OP_NONE (\f$ op(A)= A \f$),
 * <li> \ref MESS_OP_TRANSPOSE (\f$ op(A) = A^T \f$),
 * <li> \ref MESS_OP_HERMITIAN (\f$ op(A) = A^H \f$).
 * </ul>
 * Parameter \f$ i \f$ is only used to choose the solver. Possible values are:
 *
 * |   Input    |           Description                                                 |
 * |:----------:|:---------------------------------------------------------------------:|
 * |    0       |  mess_direct_create_sparse_lu                                         |
 * |    1       |  mess_direct_create_csparse_lu                                        |
 * |    2       |  mess_direct_create_umfpack                                           |
 * |    3       |  mess_direct_create_bicgstab                                          |
 * |    4       |  mess_direct_create_lapack_lu                                         |
 * |    5       |  mess_direct_create_lapack_qr                                         |
 * |    6       |  mess_direct_create_cholesky                                          |
 * |    7       |  mess_direct_create_csparse_cholesky                                  |
 * |    8       |  mess_direct_create_cholmod_cholesky                                  |
 * |    9       |  mess_direct_create_banded                                            |
 * |   10       |  mess_direct_create_superlu                                           |
 * |   11       |  mess_direct_pardiso                                                  |
 *
 *
 */
int testmat(int i, mess_matrix A) {
    int ret = 0;
    double res=1, relres=1, tol = sqrt(mess_eps());
    mess_vector x, b;
    mess_direct sol;

    mess_operation_t ops [] = {MESS_OP_NONE, MESS_OP_TRANSPOSE, MESS_OP_HERMITIAN};
    const int length_ops = sizeof(ops)/sizeof(mess_operation_t);
    int i_op;

    mess_datatype_t dts [] = {MESS_REAL, MESS_COMPLEX};
    const int length_dts = sizeof(dts)/sizeof(mess_datatype_t);
    int i_dtb,i_dtA;


    /*-----------------------------------------------------------------------------
     *  create solver
     *-----------------------------------------------------------------------------*/
    CALL(mess_direct_init(&sol));

    /*-----------------------------------------------------------------------------
     *  init vectors
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&x,&b);
    CALL(mess_vector_alloc(b, A->rows, A->data_type));
    CALL(mess_vector_rand(b));

    /*-----------------------------------------------------------------------------
     *  solve and compute residual
     *-----------------------------------------------------------------------------*/
    for(i_dtA = 0; i_dtA < length_dts; ++i_dtA){

        //cholesky based solvers need hermitian
        if(i!=6 && i != 7){
            if(dts[i_dtA]==MESS_COMPLEX){
                CALL(mess_matrix_scalec(1+2*I,A));
            }else{
                CALL(mess_matrix_toreal(A));
            }
        }

        CALL(select_solver(i,A,sol));

        for(i_op = 0; i_op < length_ops;++i_op){

            for(i_dtb = 0; i_dtb < length_dts;++i_dtb){

                /*-----------------------------------------------------------------------------
                 *  create right hand side
                 *-----------------------------------------------------------------------------*/
                if(dts[i_dtb]==MESS_COMPLEX){
                    CALL(mess_vector_tocomplex(b));
                    CALL(mess_vector_scalec(1+1*I,b));
                }else{
                    CALL(mess_vector_toreal(b));
                }

                /*-----------------------------------------------------------------------------
                 *  solve system and compute residual
                 *-----------------------------------------------------------------------------*/
                switch(i_op){
                    case 0:
                        if(MESS_HAVE_SOLVE(sol)){
                            CALL(mess_direct_solve(MESS_OP_NONE,sol,b,x));
                            CALL(mess_direct_res2(MESS_OP_NONE,A,x,b,&res,&relres));
                        }else{
                            printf("solver->solve not available.\n");
                        }
                        break;
                    case 1:
                        if(MESS_HAVE_SOLVET(sol)){
                            CALL(mess_direct_solve(MESS_OP_TRANSPOSE,sol,b,x));
                            CALL(mess_direct_res2(MESS_OP_TRANSPOSE,A,x,b,&res,&relres));
                        }else{
                            printf("solver->solvet not available.\n");

                        }
                        break;
                    case 2:
                        if(MESS_HAVE_SOLVEH(sol)){
                            CALL(mess_direct_solve(MESS_OP_HERMITIAN,sol,b,x));
                            CALL(mess_direct_res2(MESS_OP_HERMITIAN,A,x,b,&res,&relres));
                        }else{
                            printf("solver->solveh not available.\n");
                        }
                        break;
                    default:
                        printf("unknown operation type\n");
                        return 1;

                }

                /*-----------------------------------------------------------------------------
                 *  check residual
                 *-----------------------------------------------------------------------------*/
                if ( relres > tol || res > tol){
                    printf("Failed: ");
                    printf("op      = %s\n",mess_operation_t_str(ops[i_op]));
                    printf("relres  = %e\n",relres);
                    printf("res     = %e\n",res);
                    printf("tol     = %e\n",tol);
                    printf("Matrix A:\n"); mess_matrix_printinfo(A);
                    printf("Vector b:\n"); mess_vector_printinfo(b);
                    printf("Vector x:\n"); mess_vector_printinfo(x);
                    return 1;
                }
                printf("\n");

                //clear and init to check that solver handle correctly a possible resize and type conversion
                CALL(mess_vector_clear(&x));
                CALL(mess_vector_init(&x));
            }
        }

        mess_direct_clear(&sol);
        mess_direct_init(&sol);
    }

    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    mess_direct_clear(&sol);
    MESS_CLEAR_VECTORS(&x,&b);
    return 0;
}


int main ( int argc, char ** argv) {
    int ret = 0, solver=0;
    mess_error_level = 2;
    mess_matrix mat_coord, mat_csr, mat_csc, mat_dense;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if ( argc != 3) {
        printf("check the direct linear solvers\n");
        printf("usage: %s matrix.mtx op\n", argv[0]);
        printf(" op: 0 -> internal, 1 -> csparse, 2->umfpack, 3->bicgstab, 4->lapacklu, 5->lapackqr\n");
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  init and read matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&mat_coord,&mat_csr,&mat_dense,&mat_csc);
    CALL(mess_matrix_read_formated(argv[1], mat_coord, MESS_COORD));
    CALL(mess_matrix_convert(mat_coord, mat_csr, MESS_CSR));
    CALL(mess_matrix_convert(mat_coord, mat_csc, MESS_CSC));
    CALL(mess_matrix_convert(mat_coord, mat_dense, MESS_DENSE));
    solver = atoi(argv[2]);

    /*-----------------------------------------------------------------------------
     *  test with CSR
     *-----------------------------------------------------------------------------*/
    CALL(testmat(solver,mat_csr));
    CALL(testmat(solver,mat_csc));
    CALL(testmat(solver,mat_coord));
    CALL(testmat(solver,mat_dense));

    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&mat_coord,&mat_csr,&mat_dense,&mat_csc);

    return 0;
}

