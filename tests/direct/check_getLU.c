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
 * @file tests/direct/check_getLU.c
 * @brief Check the getL and getU functions of the solver.
 * @author  @mbehr
 * @test
 * This function checks the correctness of the decomposition and permutation
 * of a @ref mess_matrix structure defined by a @ref mess_direct solver.
 *
 * @}
 */
#include "../call_macro.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"


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
    default:
        return 1;

    }
    return 0;
}


int testmat(int i, mess_matrix A) {
    mess_direct sol;
    mess_matrix L,U,LU,Acopy;
    double diff;
    double eps = 1e-10;
    int ret = 0;
    int err = 0;

    CALL(mess_direct_init(&sol));
    MESS_INIT_MATRICES(&L,&U,&LU,&Acopy);
    CALL(mess_matrix_copy(A,Acopy));
    CALL(select_solver(i,A,sol));


    if ( MESS_HAVE_GETL(sol) && MESS_HAVE_GETU(sol)) {
        CALL(mess_direct_getU(sol,U));
        CALL(mess_direct_getL(sol,L));
        CALL(mess_matrix_multiply(MESS_OP_NONE,L,MESS_OP_NONE,U,LU));
        mess_int_t *p=NULL;
        mess_int_t *q=NULL;
        if(MESS_HAVE_GETPERMP(sol)){
            p = (mess_int_t *) malloc(sizeof(mess_int_t)*(A->rows));
            CALL(mess_direct_getpermp(sol,p));
        }
        if(MESS_HAVE_GETPERMQ(sol)){
            q = (mess_int_t *) malloc(sizeof(mess_int_t)*(A->cols));
            CALL(mess_direct_getpermq(sol,q));
        }

        mess_vector r=NULL,c=NULL;

        if(MESS_HAVE_GETSCALEROW(sol)){
            ret = mess_vector_init(&r);
            mess_vector_alloc(r,A->rows, MESS_REAL);
            CALL(mess_direct_getscalerow(sol,r));
            CALL(mess_matrix_rowscalem(r,Acopy));
        }

        if(MESS_HAVE_GETSCALECOL(sol)){
            ret = mess_vector_init(&c);
            mess_vector_alloc(c,A->cols, MESS_REAL);
            CALL(mess_direct_getscalecol(sol,c));
            CALL(mess_matrix_colscalem(c,Acopy));
        }

        CALL(mess_matrix_perm(Acopy,p,q));
        CALL(mess_matrix_diffnorm(LU,Acopy,&diff));
        printf("|A-LU|:   diff = %lg \n", diff);
        if ( diff > eps){
            err++;
            printf("Failed with:\n");
            CALL(mess_matrix_print(Acopy));
            CALL(mess_matrix_print(L));
            CALL(mess_matrix_print(U));
        }

        if(p)free(p);
        if(q)free(q);
        if(r)MESS_CLEAR_VECTORS(&r);
        if(c)MESS_CLEAR_VECTORS(&c);

    }else{
        printf("GETL or GETU not available.\n");
        return MESS_ERROR_ARGUMENTS;
    }


    CALL(mess_direct_clear(&sol));
    CALL(mess_matrix_clear(&L));
    CALL(mess_matrix_clear(&U));
    CALL(mess_matrix_clear(&LU));
    CALL(mess_matrix_clear(&Acopy));


    return err ;
}



int main ( int argc, char ** argv) {
    mess_matrix mat_csr;
    mess_matrix mat_csc;
    mess_matrix mat_dense;
    int ret = 0;
    int i =0;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if ( argc != 3) {
        printf("check the direct linear solvers\n");
        printf("usage: %s matrix.mtx op\n", argv[0]);
        printf(" op: 0->internal, 1->csparse, 2->umfpack, 3->bicgstab, 4->lapacklu, 5->lapackqr  6->lapack cholesky 7->csparse cholesky 8->cholmod cholesky\n");
        return 1;
    }
    /*-----------------------------------------------------------------------------
     *  init matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_init(&mat_csr));
    CALL(mess_matrix_init(&mat_dense));
    CALL(mess_matrix_init(&mat_csc));


    /*-----------------------------------------------------------------------------
     *  read matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_read_formated(argv[1],mat_csr,MESS_CSR));
    CALL(mess_matrix_convert(mat_csr, mat_csc, MESS_CSC));
    CALL(mess_matrix_convert(mat_csr, mat_dense, MESS_DENSE));

    /*-----------------------------------------------------------------------------
     *  test with CSR/CSC/DENSE
     *-----------------------------------------------------------------------------*/
    i = atoi(argv[2]);
    printf("SOLVER i = %d\n",i);
    CALL(testmat(i,mat_csr));
    CALL(testmat(i,mat_csc));
    CALL(testmat(i,mat_dense));

    mess_matrix_clear(&mat_csr);
    mess_matrix_clear(&mat_dense);
    mess_matrix_clear(&mat_csc);

    return 0;
}

