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
 * @file tests/direct/linsolvem.c
 * @brief Check the solution of a linear system of equations (matrix version).
 * @author @koehlerm
 * @test
 * This function checks the  mess_direct_solvem function defined in direct.c for a linear system of equations that means
 * it checks if the solution of
 * \f[op(A) X= B ,\f]
 * is computed correctly for a given matrices \f$ A \f$ and \f$ B \f$. \n
 * Operations \f$ op (.)\f$ on matrix \f$ A \f$ can be
 * <ul>
 * <li> \ref MESS_OP_NONE (\f$ op(A)= A \f$),
 * <li> \ref MESS_OP_TRANSPOSE (\f$ op(A) = A^T \f$),
 * <li> \ref MESS_OP_HERMITIAN (\f$ op(A) = A^H \f$).
 * </ul>
 * The solver can be selected by the third argument for the second input:
 *
 * |   Input    |           Description                                                 |
 * |:----------:|:---------------------------------------------------------------------:|
 * |    0       |  mess_direct_create_sparse_lu         defined in newlu.c              |
 * |    1       |  mess_direct_create_csparse_lu        defined in csparse.c            |
 * |    2       |  mess_direct_create_umfpack       defined in umfpack.c                |
 * |    3       |  mess_direct_create_bicgstab      defined in bicgstab.c               |
 * |    4       |  mess_direct_create_lapack_lu         defined in lapack.c             |
 * |    5       |  mess_direct_create_lapack_qr         defined in lapackqr.c           |
 * |    6       |  mess_direct_create_cholesky      defined in cholesky.c               |
 * |    7       |  mess_direct_create_csparse_cholesky  defined in csparse_chol.c       |
 * |    8       |  mess_direct_create_cholmod_cholesky  defined in cholmod_chol.c       |
 * |    9       |  mess_direct_create_banded        defined in banded.c                 |
 * |   10       |  mess_direct_create_superlu       defined in superlu.c                |
 * |   11       |  mess_direct_create_mklpardiso        defined in mklpardiso.c         |
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
#include "mess/error_macro.h"

/**
 * @brief Select a solver.
 * @param[in] i      input solver to choose
 * @param[in] A          input matrix
 * @param[out] sol  generated solver
 * @return zero on success or a non-zero error value otherwise
 *
 * The  select_solver function generates a solver for a linear system of equations
 * \f[A X= B \f]
 * Solvers are choosen with \f$ i \f$. Possible values are:
 *
 * |   Input    |           Description                                                 |
 * |:----------:|:---------------------------------------------------------------------:|
 * |    0       |  mess_direct_create_sparse_lu         defined in newlu.c              |
 * |    1       |  mess_direct_create_csparse_lu        defined in csparse.c            |
 * |    2       |  mess_direct_create_umfpack       defined in umfpack.c                |
 * |    3       |  mess_direct_create_bicgstab      defined in bicgstab.c               |
 * |    4       |  mess_direct_create_lapack_lu         defined in lapack.c             |
 * |    5       |  mess_direct_create_lapack_qr         defined in lapackqr.c           |
 * |    6       |  mess_direct_create_cholesky      defined in cholesky.c               |
 * |    7       |  mess_direct_create_csparse_cholesky  defined in csparse_chol.c       |
 * |    8       |  mess_direct_create_cholmod_cholesky  defined in cholmod_chol.c       |
 * |    9       |  mess_direct_create_banded        defined in banded.c                 |
 * |   10       |  mess_direct_create_superlu       defined in superlu.c                |
 * |   11       |  mess_direct_create_mklpardiso        defined in mklpardiso.c         |
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
 * @brief Test solving a linear system of equations (matrix version).
 * @param[in] i      input solver to choose
 * @param[in] A          input matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The  testmat function solves a linear system of equations
 * \f[op(A) X= B \f]
 * and checks if the solution is computed correctly.\n
 * Operations \f$ op (.)\f$ on matrix \f$ A \f$ can be
 * <ul>
 * <li> \ref MESS_OP_NONE (\f$ op(A)= A \f$),
 * <li> \ref MESS_OP_TRANSPOSE (\f$ op(A) = A^T \f$),
 * <li> \ref MESS_OP_HERMITIAN (\f$ op(A) = A^H \f$).
 * </ul>
 * Parameter \f$ i \f$ is only used to choose the solver. Possible values are:
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
int testmat(int i, mess_matrix A) {
    mess_direct sol;
    mess_matrix x1, x2, x3,  b1, b2, b3;
    double res1,res2, relres1, relres2;
    double res3, relres3;
    double eps = mess_eps();
    double nrm;
    int ret = 0 ;
    int err = 0;
    mess_int_t cols=10;
    mess_int_t j  ;

    CALL(mess_direct_init(&sol));
    CALL(mess_matrix_init(&x1));
    CALL(mess_matrix_init(&x2));
    CALL(mess_matrix_init(&x3));
    CALL(mess_matrix_init(&b1));
    CALL(mess_matrix_init(&b2));
    CALL(mess_matrix_init(&b3));

    CALL(mess_matrix_alloc(x1, A->rows, cols, cols*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(x2, A->rows, cols, cols*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(x3, A->rows, cols, cols*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(b1, A->rows, cols, cols*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(b2, A->rows, cols, cols*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(b3, A->rows, cols, cols*A->rows, MESS_DENSE, MESS_REAL));



    CALL(mess_matrix_norm2(A,&nrm));
    nrm = nrm *10;
    if ( MESS_IS_REAL(A) && !(i<0) ){
        for (j=0; j<A->rows*cols; j++){ x1->values[j]=j+1;}
    } else {
        CALL(mess_matrix_tocomplex(x1));
        for (j=0; j<A->rows*cols; j++){ x1->values_cpx[j]=j+1+(j*I);}
    }

    // CALL(mess_matrix_rowsums(A,b1));
    // CALL(mess_matrix_colsums(A,b2));
    CALL( mess_matrix_multiply(MESS_OP_NONE, A, MESS_OP_NONE, x1, b1));
    CALL( mess_matrix_multiply(MESS_OP_HERMITIAN, A, MESS_OP_NONE, x1, b2));
    CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE, A, MESS_OP_NONE, x1, b3));

    // mess_matrix_printshort(b2);
    // mess_matrix_printshort(b3);

    nrm *= 10;

    printf("Relative Error Bound=%lg\n",nrm*eps);
    if (i<0) {
        CALL(mess_matrix_tocomplex(b1));
        CALL(mess_matrix_tocomplex(b2));
        CALL(mess_matrix_tocomplex(b3));
        CALL(select_solver(abs(i)-1,A,sol));
    } else {
        CALL(select_solver(i,A,sol));
    }

     if ( MESS_HAVE_SOLVEM (sol )) {
        CALL(mess_direct_solvem(MESS_OP_NONE,sol,b1,x1));
        CALL(mess_direct_res2m(MESS_OP_NONE, A,x1,b1,&res1,&relres1));
        if ( relres1 > nrm*eps){printf("Failed: "); err++;}
        printf("A:   res1 = %lg \t\t relres1 = %lg\n", res1, relres1);
     } else {
        printf("solver->solvem not available.\n");
     }
    if ( MESS_HAVE_SOLVEMT(sol)){
        CALL(mess_direct_solvem(MESS_OP_TRANSPOSE,sol,b2,x2));
        CALL(mess_direct_res2m(MESS_OP_TRANSPOSE, A,x2,b2,&res2,&relres2));
        if ( relres2 > nrm*eps){printf("Failed: "); err++;}
        printf("A^T:   res2 = %lg \t\t relres2 = %lg\n", res2, relres2);
    } else {
        printf("solver->solvemt not available.\n");
    }
    if ( MESS_HAVE_SOLVEMH(sol)){
        CALL(mess_direct_solvem(MESS_OP_HERMITIAN,sol,b3,x3));
        CALL(mess_direct_res2m(MESS_OP_HERMITIAN,A,x3,b3,&res3,&relres3));
        if ( relres3 > nrm*eps){printf("Failed: "); err++;}
        printf("A^H:   res3 = %lg \t\t relres3 = %lg\n", res3, relres3);
    } else {
        printf("solver->solvemh not available.\n");
    }

    mess_direct_clear(&sol);
    mess_matrix_clear(&x1);
    mess_matrix_clear(&x2);
    mess_matrix_clear(&x3);
    mess_matrix_clear(&b1);
    mess_matrix_clear(&b2);
    mess_matrix_clear(&b3);

    return err ;
}


int main ( int argc, char ** argv) {
    mess_matrix mat_coord;
    mess_matrix mat_csr;
    mess_matrix mat_csc;
    mess_matrix mat_dense;
    int ret = 0;

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
     *  init matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_init(&mat_coord));
    CALL(mess_matrix_init(&mat_csr));
    CALL(mess_matrix_init(&mat_dense));
    CALL(mess_matrix_init(&mat_csc));


    /*-----------------------------------------------------------------------------
     *  read matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_read (argv[1], mat_coord));
    CALL(mess_matrix_convert(mat_coord, mat_csr, MESS_CSR));
    CALL(mess_matrix_convert(mat_coord, mat_csc, MESS_CSC));
    CALL(mess_matrix_convert(mat_coord, mat_dense, MESS_DENSE));


    /*-----------------------------------------------------------------------------
     *  test with CSR
     *-----------------------------------------------------------------------------*/
    CALL(testmat(atoi(argv[2]),mat_csr));
    CALL(testmat(atoi(argv[2]),mat_csc));
    CALL(testmat(atoi(argv[2]),mat_coord));
    CALL(testmat(atoi(argv[2]),mat_dense));


    mess_matrix_clear(&mat_coord);
    mess_matrix_clear(&mat_csr);
    mess_matrix_clear(&mat_dense);
    mess_matrix_clear(&mat_csc);

    return 0;
}

