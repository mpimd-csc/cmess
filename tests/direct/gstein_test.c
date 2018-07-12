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
 * @file tests/direct/gstein_test.c
 * @brief Check the solution of generalized Stein Equations.
 * @test
 *
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "../call_macro.h"


int main ( int argc, char **argv){

    mess_init();
    mess_error_level = 3;

    int ret = 0;
    double res2, res20, tol = sqrt(mess_eps());

    mess_matrix A=NULL, E=NULL, B=NULL, C=NULL, BB=NULL, CC=NULL, X=NULL;
    mess_direct gstein, stein;

    /*-----------------------------------------------------------------------------
     *  read matrices
     *-----------------------------------------------------------------------------*/
    if ( argc != 5 ) {
        fprintf(stderr, "usage: %s A.mtx E.mtx B.mtx C.mtx\n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  read matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A, &E, &B, &C, &BB, &CC, &X);

    // read data
    CALL(mess_matrix_read_formated(argv[1], A, MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[2], E, MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[3], B, MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[4], C, MESS_DENSE));

    // Build BB and CC
    CALL( mess_matrix_multiply(MESS_OP_NONE, B, MESS_OP_TRANSPOSE, B, BB));
    CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE, C, MESS_OP_NONE, C, CC));

    /*-----------------------------------------------------------------------------
     *  create solver and solve generalized stein equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_direct_init(&gstein));
    CALL(mess_direct_create_generalized_stein(A, E, gstein));

    CALL(mess_direct_solvem(MESS_OP_NONE, gstein, BB, X));
    CALL(mess_direct_generalized_stein_res2(MESS_OP_NONE,A,E,BB,X,&res2));
    CALL(mess_matrix_norm2(B, &res20));

    printf("res2            = %e\n",res2);
    printf("res2 / res20    = %e\n",res2/res20);
    if((res2 > tol) || ( res2/res20 > tol)){
        printf("FAILED (1):\n");
        printf("tol             = %e\n",tol);
       return 1;
    }

    CALL(mess_direct_solvem(MESS_OP_TRANSPOSE, gstein, CC, X));
    CALL(mess_direct_generalized_stein_res2(MESS_OP_TRANSPOSE,A,E,CC,X,&res2));
    CALL(mess_matrix_norm2(C, &res20));

    printf("res2            = %e\n",res2);
    printf("res2 / res20    = %e\n",res2/res20);
    if((res2 > tol) || ( res2/res20 > tol)){
        printf("FAILED: (2):\n");
        printf("tol             = %e\n",tol);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  create solver and solve standard stein equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_direct_init(&stein));
    CALL(mess_direct_create_generalized_stein(A, NULL, stein));

    CALL(mess_direct_solvem(MESS_OP_NONE, stein, BB, X));
    CALL(mess_direct_generalized_stein_res2(MESS_OP_NONE,A,NULL,BB,X,&res2));
    CALL(mess_matrix_norm2(B, &res20));

    printf("res2            = %e\n",res2);
    printf("res2 / res20    = %e\n",res2/res20);
    if((res2 > tol) || ( res2/res20 > tol)){
        printf("FAILED (3):\n");
        printf("tol             = %e\n",tol);
        return 1;
    }

    CALL(mess_direct_solvem(MESS_OP_TRANSPOSE, stein, CC, X));
    CALL(mess_direct_generalized_stein_res2(MESS_OP_TRANSPOSE,A,NULL,CC,X,&res2));
    CALL(mess_matrix_norm2(C, &res20));

    printf("res2            = %e\n",res2);
    printf("res2 / res20    = %e\n",res2/res20);
    if((res2 > tol) || ( res2/res20 > tol)){
        printf("FAILED (4):\n");
        printf("tol             = %e\n",tol);
        return 1;
    }


    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    mess_direct_clear(&stein);
    mess_direct_clear(&gstein);
    MESS_CLEAR_MATRICES(&A, &E, &B, &C, &BB, &CC, &X);

    return 0;
}

