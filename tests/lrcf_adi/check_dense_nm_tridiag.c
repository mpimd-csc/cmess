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
 * @addtogroup test_lrcfadi
 * @{
 *
 * @file tests/lrcf_adi/check_dense_nm_tridiag.c
 * @brief Check the solution of a generalized Riccati Equation in transposed and nontransposed version using @ref mess_dense_nm_gmpare.
 * @test
 * This function checks the @ref mess_dense_nm_gmpare function for small input, that means it checks if the Riccati Equation
 * \f[ A^T X + X A +/- X G X + Q  = 0 \f]
 * using a dense Newton method.
 *
 * @attention The Newton method does not converge for this instance using the positive Riccati equation.
 *
 * @}
 */



#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "../call_macro.h"


int main ( int argc, char **argv){

    char plus_char='+';
    mess_int_t dim, plus, linesearch, ret;
    mess_operation_t op;
    double alpha=5;
    double res, rel;
    mess_matrix A, X, Q, G;

    mess_version();
    mess_init();
    mess_error_level = 3;


    const double absres_tol = sqrt(mess_eps()), relres_tol = sqrt(mess_eps());
    const mess_int_t maxit = 50;
    const mess_norm_t nrm = MESS_2_NORM;

    /*-----------------------------------------------------------------------------
     * check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 5) {
        fprintf(stderr, "Usage: %s dim plus linesearch op\nplus, linesearch and op element of {0,1}\n", argv[0]);
        return 1;
    }


    /*-----------------------------------------------------------------------------
     * read data input arguments
     *-----------------------------------------------------------------------------*/
    dim  = atoi(argv[1]);
    plus  = atoi(argv[2]);
    linesearch  = atoi(argv[3]);
    op = atoi(argv[4]) ? MESS_OP_TRANSPOSE:MESS_OP_NONE;
    plus_char = plus ? '+' : '-';

    /*-----------------------------------------------------------------------------
     * create matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A, &X, &Q, &G);
    CALL(mess_matrix_tridiag(A, dim, dim, MESS_DENSE, MESS_REAL, alpha, -1, -alpha));
    CALL(mess_matrix_alloc(Q,dim,dim,dim*dim,MESS_DENSE,MESS_REAL));
    CALL(mess_matrix_alloc(G,dim,dim,dim*dim,MESS_DENSE,MESS_REAL));
    CALL(mess_matrix_ones(Q));
    CALL(mess_matrix_ones(G));


    /*-----------------------------------------------------------------------------
     * solve Riccati Equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_dense_nm_gmpare(NULL, A, NULL, Q ,G, plus, linesearch, op, maxit, nrm, absres_tol, relres_tol, NULL, NULL, NULL, X));


    /*-----------------------------------------------------------------------------
     * compute residual
     *-----------------------------------------------------------------------------*/
    CALL(mess_dense_res_gmpare(A, NULL, Q, G, X, plus, op, nrm, &res, &rel));


    /*-----------------------------------------------------------------------------
     *  print results
     *-----------------------------------------------------------------------------*/
    printf("A: "MESS_PRINTF_INT"-x-"MESS_PRINTF_INT"\n", A->rows, A->cols);

    printf("abs./rel. residual tolerance = %e / %e\n", absres_tol, relres_tol);
    if(op == MESS_OP_TRANSPOSE){
        printf("A^T*X + X*A %c X*G*X + Q = 0, linesearch=%d, abs./rel. %s residual = %e / %e\n", plus_char, linesearch, mess_norm_t_str(nrm), res, rel);
    } else {
        printf("A*X + X*A^T %c X*G*X + Q = 0, linesearch=%d, abs./rel. %s residual = %e / %e\n", plus_char, linesearch, mess_norm_t_str(nrm), res, rel);
    }


    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&X,&Q,&G);

    return rel >= relres_tol || res >= absres_tol;

}

