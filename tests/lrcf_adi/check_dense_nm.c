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
 * @file tests/lrcf_adi/check_dense_nm.c
 * @brief Check the dense Newton Method for dense algebraic Riccati equations.
 * @test
 * This function checks the @ref mess_dense_nm_gmpare function that means it checks if the Riccati Equation
 *  \f[
 *      A^T X E + E^T X A -/+ E^T X BB^T X E + C^TC = 0
 *  \f]
 * is correctly solved by the dense Newton Method.
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

    mess_error_level = 3;
    mess_operation_t op = MESS_OP_NONE;
    char plus_char='+';
    int plus, linesearch, ret;
    double res, rel;
    const double absres_tol = sqrt(mess_eps()), relres_tol = sqrt(mess_eps());
    const mess_int_t maxit = 50;
    const mess_norm_t nrm = MESS_2_NORM;
    mess_matrix X0 = NULL, A, E = NULL, B, C, X, Q, G;

    mess_version();
    mess_init();


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 8  && argc !=9){
        fprintf(stderr, "Usage: %s A.mtx [E.mtx] B.mtx C.mtx plus linesearch op", argv[0]);
        return 1;
    }


    /*-----------------------------------------------------------------------------
     *  init, read and prepare matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A, &X, &B, &C, &Q, &G);

    if(argc==8){
        CALL(mess_matrix_read_formated(argv[1], A, MESS_DENSE));
        CALL(mess_matrix_read_formated(argv[2], B, MESS_DENSE));
        CALL(mess_matrix_read_formated(argv[3], C, MESS_DENSE));
        plus = atoi(argv[4]);
        linesearch  = atoi(argv[5]);
        op  = atoi(argv[6]) ? MESS_OP_TRANSPOSE:MESS_OP_NONE;
    }else{
        CALL(mess_matrix_init(&E));
        CALL(mess_matrix_read_formated(argv[1], A, MESS_DENSE));
        CALL(mess_matrix_read_formated(argv[2], E, MESS_DENSE));
        CALL(mess_matrix_read_formated(argv[3], B, MESS_DENSE));
        CALL(mess_matrix_read_formated(argv[4], C, MESS_DENSE));
        plus = atoi(argv[5]);
        linesearch  = atoi(argv[6]);
        op  = atoi(argv[7]) ? MESS_OP_TRANSPOSE:MESS_OP_NONE;
    }

    plus_char = plus ? '+':'-';

    CALL(mess_matrix_multiply(MESS_OP_NONE, B, MESS_OP_HERMITIAN, B, G));
    CALL(mess_matrix_multiply(MESS_OP_HERMITIAN, C, MESS_OP_NONE, C, Q));



    /*-----------------------------------------------------------------------------
     *  Setup and Solve Riccati Equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_dense_nm_gmpare(X0, A, E, Q, G, plus, linesearch, op, maxit, nrm, absres_tol, relres_tol, NULL, NULL, NULL, X));

    /*-----------------------------------------------------------------------------
     *  calculate and print residual
     *-----------------------------------------------------------------------------*/
    CALL(mess_dense_res_gmpare(A, E, Q, G, X, plus, op, nrm, &res, &rel));

    if(E){
        if(op == MESS_OP_NONE){
            printf("A*X*E^T + E*X*A^T %c E*X*G*X*E^T + Q = 0, X0=%p, linesearch=%d, abs./rel. %s residual: %e / %e\n", plus_char, X0, linesearch, mess_norm_t_str(nrm), res, rel);
        }else{
            printf("A^T*X*E + E^T*X*A %c E^T*X*G*X*E + Q = 0, X0=%p, linesearch=%d, abs./rel. %s residual: %e / %e\n", plus_char, X0, linesearch, mess_norm_t_str(nrm), res, rel);
        }
    }else{
        if(op == MESS_OP_NONE){
            printf("A*X + X*A^T %c X*G*X + Q = 0, X0=%p, linesearch=%d, abs./rel. %s residual: %e / %e\n", plus_char, X0, linesearch, mess_norm_t_str(nrm), res, rel);
        }else{
            printf("A^T*X + X*A %c X*G*X + Q = 0, X0=%p, linesearch=%d, abs./rel. %s residual: %e / %e\n", plus_char, X0, linesearch, mess_norm_t_str(nrm), res, rel);
        }
    }

    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&B,&C,&X,&Q,&G);
    if(E) mess_matrix_clear(&E);
    if(X0) mess_matrix_clear(&X0);

    return rel >= relres_tol || res >= absres_tol;

}

