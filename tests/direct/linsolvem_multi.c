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
 * @file tests/direct/linsolvem_multi.c
 * @brief Check the solution of a shifted linear system of equations (matrix version).
 * @author @koehlerm
 * @test
 * This function checks the @ref mess_multidirect_solvem function defined in multidirect.c for a linear system of equations
 * that means it checks if the solution of
 * \f[op(A + p I) X= B ,\f]
 * is computed correctly for a given matrices \f$ A \f$ and \f$ B \f$, identity matrix \f$ I \f$ and
 * shifts \f$ p= 1 \f$ (real case), \f$ p = 4-0.2 i \f$ (complex case). \n
 * Operations \f$ op (.)\f$ on matrix \f$ A + p I \f$ can be
 * * \ref MESS_OP_NONE (\f$ op(A + p I)= A + p I \f$),
 * * \ref MESS_OP_TRANSPOSE (\f$ op(A + p I) = A^T + p I \f$),
 * * \ref MESS_OP_HERMITIAN (\f$ op(A) = A^H + p I \f$).
 * The solver can be selected by the third argument for the second input:
 * *
 * * 0 - \ref mess_multidirect_create_sparse_lu solver defined in multilu_advanced.c,
 * * 1 - \ref mess_multidirect_create_umfpack solver defined in multilu_umfpack.c,
 * *
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
 * @param[in] shifts     input shifts
 * @param[out] sol  generated solver
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref select_solver function generates a solver for a shifted linear system of equations
 * \f[(A + shifts I) X= B \f]
 * Solvers are choosen with \f$ i \f$. Possible values are:
 * <ul>
 * <li> 0 - \ref mess_multidirect_create_sparse_lu,
 * <li> 1 - \ref mess_multidirect_create_umfpack,
 * </ul>
 *
 */
int select_solver(int i, mess_matrix A, mess_vector shifts, mess_multidirect sol) {
    int ret =0;
    switch(i){
        case 0:
            CALL(mess_multidirect_create_sparse_lu(A, NULL, shifts, sol,NULL,NULL));
            break;
#ifdef MESS_HAVE_UMFPACK
        case 1:
            CALL(mess_multidirect_create_umfpack(A, NULL, shifts, sol,NULL,NULL));
            break;
#else
            return 0;
#endif
        default:
            return 1;

    }
    return 0;
}


/**
 * @brief Test solving a shifted linear system of equations (matrix version).
 * @param[in] i      input solver to choose
 * @param[in] cpx    input true if complex values are used
 * @param[in] A          input matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref linsolvem_multi_testmat function solves a shifted linear system of equations
 * \f[op(A + p I) X= B \f]
 * and checks if the solution is computed correctly.\n
 * Operations \f$ op (.)\f$ on matrix \f$ A + p I \f$ can be
 * <ul>
 * <li> \ref MESS_OP_NONE (\f$ op(A + p I)= A + p I \f$),
 * <li> \ref MESS_OP_TRANSPOSE (\f$ op(A + p I) = A^T + p I \f$),
 * <li> \ref MESS_OP_HERMITIAN (\f$ op(A) = A^H + p I \f$).
 * </ul>
 * Parameter \f$ i \f$ is only used to choose the solver. Possible values are:
 * <ul>
 * <li> 0 - \ref mess_multidirect_create_sparse_lu,
 * <li> 1 - \ref mess_multidirect_create_umfpack,
 * </ul>
 *
 */
int linsolvem_multi_testmat(int i, int cpx, mess_matrix A) {
    mess_multidirect sol;
    mess_matrix x11, x12, x13,  b11, b12, b13;
    mess_matrix x21, x22, x23,  b21, b22, b23;
    double res11,res12, relres11, relres12 ,res13, relres13;
    double res21,res22, relres21, relres22 ,res23, relres23;

    double eps = mess_eps();
    double nrm1, nrm2;
    int ret = 0 ;
    int err = 0;
    mess_int_t j  ;
    mess_matrix A1, A2, Identity;
    mess_vector shifts;

    CALL(mess_matrix_init(&A1));
    CALL(mess_matrix_init(&A2));
    CALL(mess_matrix_init(&Identity));
    MESS_INIT_VECTORS(&shifts);
    CALL(mess_vector_alloc(shifts, 2, MESS_COMPLEX));

    CALL(mess_multidirect_init(&sol));
    CALL(mess_matrix_init(&x11));
    CALL(mess_matrix_init(&x12));
    CALL(mess_matrix_init(&x13));
    CALL(mess_matrix_init(&b11));
    CALL(mess_matrix_init(&b12));
    CALL(mess_matrix_init(&b13));
    CALL(mess_matrix_init(&x21));
    CALL(mess_matrix_init(&x22));
    CALL(mess_matrix_init(&x23));
    CALL(mess_matrix_init(&b21));
    CALL(mess_matrix_init(&b22));
    CALL(mess_matrix_init(&b23));

    CALL(mess_matrix_alloc(x11, A->rows, 4, 4*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(x12, A->rows, 4, 4*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(x13, A->rows, 4, 4*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(b11, A->rows, 4, 4*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(b12, A->rows, 4, 4*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(b13, A->rows, 4, 4*A->rows, MESS_DENSE, MESS_REAL));

    CALL(mess_matrix_alloc(x21, A->rows, 4, 4*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(x22, A->rows, 4, 4*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(x23, A->rows, 4, 4*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(b21, A->rows, 4, 4*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(b22, A->rows, 4, 4*A->rows, MESS_DENSE, MESS_REAL));
    CALL(mess_matrix_alloc(b23, A->rows, 4, 4*A->rows, MESS_DENSE, MESS_REAL));

    shifts->values_cpx[0] = 1;
    shifts->values_cpx[1] = 4-0.2*I;

    CALL(select_solver(i,A,shifts,sol));

    /*-----------------------------------------------------------------------------
     *  Assemble matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_copy(A,A1));
    CALL(mess_matrix_copy(A,A2));
    CALL(mess_matrix_eye(Identity,A->rows,A->rows,A->store_type));
    CALL(mess_matrix_add(1,Identity,1,A1));
    CALL(mess_matrix_eyec(Identity,A->rows,A->rows,A->store_type));
    CALL(mess_matrix_addc(4-0.2*I,Identity,1,A2));
    CALL(mess_matrix_norm2(A1,&nrm1));
    CALL(mess_matrix_norm2(A2,&nrm2));

    nrm1 = nrm1 *1000;
    nrm2 = nrm2 *1000;

    if ( MESS_IS_REAL(A) && !cpx ){
        for (j=0; j<A->rows*4; j++){ x11->values[j]=j+1;}
    } else {
        CALL(mess_matrix_tocomplex(x11));
        for (j=0; j<A->rows*4; j++){ x11->values_cpx[j]=j+1+(j*I);}
    }

    CALL( mess_matrix_multiply(MESS_OP_NONE, A1, MESS_OP_NONE, x11, b11));
    CALL( mess_matrix_multiply(MESS_OP_HERMITIAN, A1, MESS_OP_NONE, x11, b12));
    CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE, A1, MESS_OP_NONE, x11, b13));
    CALL( mess_matrix_multiply(MESS_OP_NONE, A2, MESS_OP_NONE, x11, b21));
    CALL( mess_matrix_multiply(MESS_OP_HERMITIAN, A2, MESS_OP_NONE, x11, b22));
    CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE, A2, MESS_OP_NONE, x11, b23));

    if (cpx) {
        CALL(mess_matrix_tocomplex(b11));
        CALL(mess_matrix_tocomplex(b12));
        CALL(mess_matrix_tocomplex(b13));
        CALL(mess_matrix_tocomplex(b21));
        CALL(mess_matrix_tocomplex(b22));
        CALL(mess_matrix_tocomplex(b23));
    }

    if ( MESS_HAVE_SOLVEM (sol )) {
        CALL(mess_multidirect_solvem(MESS_OP_NONE,sol,0,b11,x11));
        CALL(mess_direct_res2m(MESS_OP_NONE, A1,x11,b11,&res11,&relres11));
        CALL(mess_multidirect_solvem(MESS_OP_NONE,sol,1,b21,x21));
        CALL(mess_direct_res2m(MESS_OP_NONE, A2,x21,b21,&res21,&relres21));
    } else {
        res11=0;
        relres11=0;
        res21=0;
        relres21=0;
    }

    if ( MESS_HAVE_SOLVEMH (sol )) {
        CALL(mess_multidirect_solvem(MESS_OP_HERMITIAN,sol,0,b12,x12));
        CALL(mess_direct_res2m(MESS_OP_HERMITIAN, A1,x12,b12,&res12,&relres12));
        CALL(mess_multidirect_solvem(MESS_OP_HERMITIAN,sol,1,b22,x22));
        CALL(mess_direct_res2m(MESS_OP_HERMITIAN, A2,x22,b22,&res22,&relres22));
    } else {
        res12=0;
        relres12=0;
        res22=0;
        relres22=0;
    }

    if ( MESS_HAVE_SOLVEMT (sol )) {
        CALL(mess_multidirect_solvem(MESS_OP_TRANSPOSE,sol,0,b13,x13));
        CALL(mess_direct_res2m(MESS_OP_TRANSPOSE, A1,x13,b13,&res13,&relres13));
        CALL(mess_multidirect_solvem(MESS_OP_TRANSPOSE,sol,1,b23,x23));
        CALL(mess_direct_res2m(MESS_OP_TRANSPOSE, A2,x23,b23,&res23,&relres23));
    } else {
        res13=0;
        relres13=0;
        res23=0;
        relres23=0;
    }

    printf("A(0):   res1 = %lg \t\t relres1 = %lg\n", res11, relres11);
    printf("A(0)^H: res2 = %lg \t\t relres2 = %lg\n", res12, relres12);
    printf("A(0)^T: res3 = %lg \t\t relres3 = %lg\n", res13, relres13);
    printf("A(1):   res1 = %lg \t\t relres1 = %lg\n", res21, relres21);
    printf("A(1)^H: res2 = %lg \t\t relres2 = %lg\n", res22, relres22);
    printf("A(1)^T: res3 = %lg \t\t relres3 = %lg\n", res23, relres23);

    if ( relres11 > nrm1*eps) err++;
    if ( relres12 > nrm1*eps) err++;
    if ( relres13 > nrm1*eps) err++;
    if ( relres21 > nrm2*eps) err++;
    if ( relres22 > nrm2*eps) err++;
    if ( relres23 > nrm2*eps) err++;

    mess_multidirect_clear(&sol);
    mess_matrix_clear(&x11);
    mess_matrix_clear(&x12);
    mess_matrix_clear(&x13);
    mess_matrix_clear(&b11);
    mess_matrix_clear(&b12);
    mess_matrix_clear(&b13);
    mess_matrix_clear(&x21);
    mess_matrix_clear(&x22);
    mess_matrix_clear(&x23);
    mess_matrix_clear(&b21);
    mess_matrix_clear(&b22);
    mess_matrix_clear(&b23);

    mess_matrix_clear(&Identity);
    mess_matrix_clear(&A1);
    mess_matrix_clear(&A2);
    mess_vector_clear(&shifts);


    return err ;
}


int main ( int argc, char ** argv) {
    mess_init();
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
    CALL(linsolvem_multi_testmat(atoi(argv[2]),0,mat_csr));
    // CALL(linsolvem_multi_testmat(atoi(argv[2]),0,mat_csc));
    // CALL(linsolvem_multi_testmat(atoi(argv[2]),0,mat_coord));
    // CALL(linsolvem_multi_testmat(atoi(argv[2]),0,mat_dense));


    mess_matrix_clear(&mat_coord);
    mess_matrix_clear(&mat_csr);
    mess_matrix_clear(&mat_dense);
    mess_matrix_clear(&mat_csc);

    return 0;
}

