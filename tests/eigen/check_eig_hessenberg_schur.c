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
 * @file tests/eigen/check_eig_hessenberg_schur.c
 * @brief Check the @ref mess_eigen_hessenberg_schur, @ref mess_eigen_schur_to_evd.
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



int main(int argc, char ** argv){
    mess_init();
    mess_error_level = 3;
    mess_int_t ret = 0, seed = 1, dim = 50, k_triu=-1;
    double diff = 0, tol = sqrt(mess_eps());
    mess_datatype_t dt = MESS_REAL;
    mess_vector ev, ev_test;
    mess_matrix hessA, hessQ, hessT, hessV, Tmp1, Tmp2, Tmp3, eye;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    if ( argc != 4 ) {
        fprintf (stderr, "Usage: %s dim cpx k_triu\n", argv[0]) ;
        exit(1);
    }


    dim     = atoi(argv[1]);
    dt      = atoi(argv[2])?MESS_COMPLEX:MESS_REAL;
    k_triu  = atoi(argv[3]);

    if(k_triu<-1){
        printf("ktriu must be larger or equal to -1, otherwise no Hessenberg matrix is generated\n");
        ret = 1;
        goto clear;
    }

    /*-----------------------------------------------------------------------------
     *  init matrices and vectors
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&hessA, &hessQ, &hessT, &hessV, &Tmp1, &Tmp2, &Tmp3, &eye);
    MESS_INIT_VECTORS(&ev, &ev_test);

    /*-----------------------------------------------------------------------------
     *  create upper hessenberg matrix
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_rand_init(&seed));
    CALL(mess_matrix_rand(hessA, dim, dim, MESS_DENSE, dt, 1.0));
    CALL(mess_matrix_triu(hessA, k_triu));
    CALL(mess_matrix_eye(eye, hessA->rows, hessA->cols, MESS_DENSE));

    /*-----------------------------------------------------------------------------
     *  call and check mess_eigen hessenberg schur function
     *-----------------------------------------------------------------------------*/
    CALL(mess_eigen_hessenberg_schur(hessA, ev, hessQ, hessT));

    //call again without schur transformation as a result
    CALL(mess_eigen_hessenberg_schur(hessA, ev_test, NULL, NULL));

    //check if Q is unitary
    CALL(mess_matrix_multiply(MESS_OP_HERMITIAN, hessQ, MESS_OP_NONE, hessQ, Tmp1));
    CALL(mess_matrix_diffnorm(Tmp1, eye, &diff));

    printf("||Q^H Q  - I ||_2 =%e\n",diff);
    if(diff>=tol){
        printf("Failed;\n");
        ret = 1;
        goto clear;
    }

    //check if decomposition is correct
    CALL(mess_matrix_multiply(MESS_OP_NONE, hessQ, MESS_OP_NONE,  hessT, Tmp1));
    CALL(mess_matrix_multiply(MESS_OP_NONE, Tmp1, MESS_OP_HERMITIAN, hessQ, Tmp2));
    CALL(mess_matrix_diffnorm(Tmp2, hessA, &diff));

    printf("||hessA - Q T Q^H || _2 = %e\n", diff);
    if(diff>=tol){
        printf("Failed:\n");
        ret = 1;
        goto clear;
    }

    //compare eigenvalues of both calls
    CALL(mess_vector_diffnorm(ev, ev_test, &diff));
    printf("||ev-ev_test||_2 = %e\n", diff);
    if(diff>=tol){
        printf("Failed:\n");
        printf("ev\n");         mess_vector_print(ev);
        printf("ev_test\n");    mess_vector_print(ev_test);
        ret = 1;
        goto clear;
    }

    /*-----------------------------------------------------------------------------
     *  call mess_eigen_schur_to_evd and check and check function
     *-----------------------------------------------------------------------------*/
    CALL(mess_eigen_schur_to_evd(hessT, hessQ, hessV));
    CALL(mess_matrix_multiply(MESS_OP_NONE, hessA, MESS_OP_NONE, hessV, Tmp1));
    CALL(mess_matrix_diag_from_vector(ev, Tmp2));
    CALL(mess_matrix_multiply(MESS_OP_NONE, hessV, MESS_OP_NONE, Tmp2, Tmp3));
    CALL(mess_matrix_diffnorm(Tmp3, Tmp1, &diff));

    //check eigenvalue/eigenvector property
    printf("|| A V  -  V D ||_2 = %e\n", diff);
    if(diff>=tol){
        printf("Failed:\n");
        ret = 1;
        goto clear;
    }

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
clear:
    MESS_CLEAR_MATRICES(&hessA, &hessQ, &hessT, &hessV, &Tmp1, &Tmp2, &Tmp3, &eye);
    MESS_CLEAR_VECTORS(&ev, &ev_test);

    return ret;
}

