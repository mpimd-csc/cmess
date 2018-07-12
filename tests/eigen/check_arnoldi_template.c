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
 * @addtogroup test_eigen
 * @{
 * @file tests/eigen/check_arnoldi_template.c
 * @brief Check the @ref mess_eigen_arnoldi_template function.
 * @author @mbehr
 * @test
 *
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mess/mess.h"

#include "../call_macro.h"



#define CHECK_ARNOLDI(K, A, MVPA, SV, H, V){                                                \
                                                                                            \
        /*check if k is to large*/                                                          \
        if(K>=A->rows){                                                                     \
            printf("k="MESS_PRINTF_INT" is larger than the order of A=" MESS_PRINTF_INT     \
                    ", therefore skip this test.\n", K, A->rows);                           \
        }else{                                                                              \
        CALL(mess_eigen_arnoldi_template(MVPA, K, SV, H, V));                               \
                                                                                            \
        /* V'*V = I*/                                                                       \
        CALL(mess_matrix_eye(eye, V->cols, V->cols, MESS_DENSE));                           \
        CALL(mess_matrix_multiply(MESS_OP_HERMITIAN, V, MESS_OP_NONE, V, Tmp1));            \
        CALL(mess_matrix_diffnorm(Tmp1, eye, &diff));                                       \
        printf("k="MESS_PRINTF_INT",    || V' * V - eye ||_2        = %e\n", K, diff);      \
                                                                                            \
        if(diff>tol){                                                                       \
            printf("Failed, tol =%e\n",tol);                                                \
            ret = 1;                                                                        \
            goto clear;                                                                     \
        }                                                                                   \
                                                                                            \
        /* A*V(:,1:k) = V*H */                                                              \
        if(H->rows != H->cols){                                                             \
            CALL(mess_matrix_colsub(V,0,K-1,Vtilde));                                       \
        }else{                                                                              \
            /* early breakdown */                                                           \
            CALL(mess_matrix_copy(V,Vtilde));                                               \
        }                                                                                   \
        CALL(mess_matrix_multiply(MESS_OP_NONE, A, MESS_OP_NONE, Vtilde, Tmp1));            \
        CALL(mess_matrix_multiply(MESS_OP_NONE, V, MESS_OP_NONE, H, Tmp2));                 \
        CALL(mess_matrix_diffnorm(Tmp1, Tmp2, &diff));                                      \
                                                                                            \
        printf("k="MESS_PRINTF_INT",    || A*V(:,1:k) - H*V ||_2    = %e\n", K, diff);      \
                                                                                            \
        if(diff>tol){                                                                       \
            printf("Failed, tol =%e\n",tol);                                                \
            ret = 1;                                                                        \
            goto clear;                                                                     \
        }                                                                                   \
            if(A->data_type != H->data_type || A->data_type != V->data_type){               \
                printf("H or V have wrong datatype\n.");                                    \
                mess_matrix_printinfo(A);                                                   \
                mess_matrix_printinfo(H);                                                   \
                mess_matrix_printinfo(V);                                                   \
                ret = 1;                                                                    \
                goto clear;                                                                 \
            }                                                                               \
                                                                                            \
        /*check size of matrices*/                                                          \
        if( !(H->rows == H->cols && V->cols == H->rows) &&                                  \
            !(H->rows -1 == H->cols && H->cols == K && V->cols == H->rows)                  \
          ){                                                                                \
            printf("Failed:H or V has wrong size\n");                                       \
            printf("k="MESS_PRINTF_INT"\n",K);                                              \
            printf("H:\n");                                                                 \
            mess_matrix_printinfo(H);                                                       \
            printf("V:\n");                                                                 \
            mess_matrix_printinfo(V);                                                       \
            ret = 1;                                                                        \
            goto clear;                                                                     \
        }                                                                                   \
        }                                                                                   \
        printf("\n");                                                                       \
}


int main ( int argc, char ** argv) {

    mess_init();
    mess_error_level=3;
    mess_int_t ret = 0;
    double diff = 0, tol = sqrt(mess_eps());

    mess_vector sv, svc;
    mess_matrix A, Ac, H, V, Vtilde, Tmp1, Tmp2, eye;
    mess_mvpcall mvpA, mvpAc;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 2) {
        printf("usage: %s A.mtx \n", argv[0]);
        return 1;
    }


    /*-----------------------------------------------------------------------------
     *  init matrices and vectors
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&Ac,&H,&V,&Vtilde,&Tmp1,&Tmp2,&eye);
    MESS_INIT_VECTORS(&sv,&svc);

    CALL(mess_matrix_read_formated(argv[1], A, MESS_CSC));
    CALL(mess_mvpcall_matrix(&mvpA, MESS_OP_NONE, A));

    CALL(mess_matrix_copy(A,Ac));
    CALL(mess_matrix_scalec(1+2*I,Ac));
    CALL(mess_mvpcall_matrix(&mvpAc, MESS_OP_NONE, Ac));

    CALL(mess_vector_alloc(sv, A->rows, MESS_REAL));
    CALL(mess_vector_rand(sv));

    CALL(mess_vector_alloc(svc, A->rows, MESS_COMPLEX));
    CALL(mess_vector_rand(svc));


    /*-----------------------------------------------------------------------------
     *  call function and check results
     *-----------------------------------------------------------------------------*/
    printf("Real Matrix\n");
    CHECK_ARNOLDI(1,    A, mvpA, sv, H, V)
    CHECK_ARNOLDI(2,    A, mvpA, sv, H, V)
    CHECK_ARNOLDI(5,    A, mvpA, sv, H, V)
    CHECK_ARNOLDI(13,   A, mvpA, sv, H, V)
    CHECK_ARNOLDI(300,  A, mvpA, sv, H, V)

    printf("Complex Matrix\n");
    CHECK_ARNOLDI(1,    Ac, mvpAc, svc, H, V)
    CHECK_ARNOLDI(2,    Ac, mvpAc, svc, H, V)
    CHECK_ARNOLDI(5,    Ac, mvpAc, svc, H, V)
    CHECK_ARNOLDI(13,   Ac, mvpAc, svc, H, V)
    CHECK_ARNOLDI(300,  Ac, mvpAc, svc, H, V)


    /*-----------------------------------------------------------------------------
     * add block matrix and force early breakdown
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_rand(Tmp1,10,10,MESS_DENSE, MESS_REAL,1.0));
    CALL(mess_matrix_cat(Tmp1, NULL, NULL, A, MESS_CSC, Tmp2));
    CALL(mess_matrix_copy(Tmp2, A));
    CALL(mess_mvpcall_clear(&mvpA));
    CALL(mess_mvpcall_matrix(&mvpA, MESS_OP_NONE, A));
    CALL(mess_vector_resize(sv,A->rows));
    CALL(mess_vector_zeros(sv));
    sv->values[0]=1;

    CALL(mess_matrix_rand(Tmp1,10,10,MESS_DENSE, MESS_COMPLEX,1.0));
    CALL(mess_matrix_cat(Tmp1, NULL, NULL, Ac, MESS_CSC, Tmp2));
    CALL(mess_matrix_copy(Tmp2, Ac));
    CALL(mess_mvpcall_clear(&mvpAc));
    CALL(mess_mvpcall_matrix(&mvpAc, MESS_OP_NONE, Ac));
    CALL(mess_vector_resize(svc,Ac->rows));
    CALL(mess_vector_zeros(svc));
    svc->values_cpx[0]=1;

    printf("Check for early breakdown\n");
    printf("Real Matrix\n");
    CHECK_ARNOLDI(11,   A, mvpA, sv, H, V)
    CHECK_ARNOLDI(50,   A, mvpA, sv, H, V)
    CHECK_ARNOLDI(100,  A, mvpA, sv, H, V)

    printf("Complex Matrix\n");
    CHECK_ARNOLDI(11,   Ac, mvpAc, svc, H, V)
    CHECK_ARNOLDI(50,   Ac, mvpAc, svc, H, V)
    CHECK_ARNOLDI(100,  Ac, mvpAc, svc, H, V)


    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
clear:
    MESS_CLEAR_MATRICES(&A,&Ac,&H,&V,&Vtilde,&Tmp1,&Tmp2,&eye);
    MESS_CLEAR_VECTORS(&sv,&svc);
    CALL(mess_mvpcall_clear(&mvpA));
    CALL(mess_mvpcall_clear(&mvpAc));

    return ret;
}

