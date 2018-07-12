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
 * @file tests/direct/check_hessenberg.c
 * @brief Check the @ref mess_direct_create_hessenberg_lu solver.
 * @author  @mbehr
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


int main ( int argc, char ** argv) {
    mess_init();
    int ret = 0;
    mess_int_t dim = 0, seed = 1234, i_dts = 0, i_ops = 0, cols = 3;
    double tol = sqrt(mess_eps()), res = 0, relres = 0;
    mess_datatype_t dt = MESS_REAL;
    mess_datatype_t dts [] = {MESS_REAL, MESS_COMPLEX};
    mess_operation_t ops [] = {MESS_OP_NONE, MESS_OP_TRANSPOSE, MESS_OP_HERMITIAN};

    mess_matrix H, B, X;
    mess_vector x, b;
    mess_direct hessenberg;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if ( argc != 3) {
        printf("usage: %s dim cpx\n", argv[0]);
        return 1;
    }

    dim = atoi(argv[1]);
    dt = atoi(argv[2])?MESS_COMPLEX:MESS_REAL;

    /*-----------------------------------------------------------------------------
     * create matrix and rhs
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&H, &B, &X);

    mess_matrix_rand_init(&seed);
    mess_matrix_alloc(H, dim, dim, dim*dim, MESS_DENSE, dt);
    mess_matrix_rand(H, dim, dim, MESS_DENSE, dt, 1.0);
    mess_matrix_triu(H,-1);

    MESS_INIT_VECTORS(&x,&b);

    /*-----------------------------------------------------------------------------
     *  create hessenberg solver and solve
     *-----------------------------------------------------------------------------*/
    MESS_INIT_DIRECTS(&hessenberg);
    mess_direct_create_hessenberg_lu(H, hessenberg);

    for(i_ops=0;i_ops<3;++i_ops){

        for(i_dts=0;i_dts<2;++i_dts){

            /*-----------------------------------------------------------------------------
             *  check solve* functions
             *-----------------------------------------------------------------------------*/
            //reset right hand sides
            mess_vector_reset(b);
            mess_vector_alloc(b, dim, dts[i_dts]);
            mess_vector_rand(b);

            //solve and compute residual
            mess_direct_solve(ops[i_ops], hessenberg, b, x);

            mess_direct_res2(ops[i_ops], H, x, b, &res, &relres);
            printf("Operation type          = %s\n", mess_operation_t_str(ops[i_ops]));
            printf("tol                     = %e\n", tol);
            printf("abs. / rel. residual    = %e / %e\n\n", res, relres);

            if(res > tol || relres >tol){
                printf("Failed solve*:\n");
                printf("Operation type          = %s\n", mess_operation_t_str(ops[i_ops]));
                printf("right hand side;\n");
                mess_vector_printinfo(b);
                printf("Hessenberg Matrix:\n");
                mess_matrix_printinfo(H);
                printf("tol                     = %e\n", tol);
                printf("abs. / rel. residual    = %e / %e\n", res, relres);
                ret = 1;
                goto clear;
            }

            /*-----------------------------------------------------------------------------
             *  check solvem* functions
             *-----------------------------------------------------------------------------*/
            //reset right hand sides
            mess_matrix_reset(B);
            mess_matrix_rand(B, dim, cols, MESS_DENSE, dts[i_dts], 1.0);

            //solve and compute residual
            mess_direct_solvem(ops[i_ops], hessenberg, B, X);
            mess_direct_res2m(ops[i_ops], H, X, B, &res, &relres);

            printf("Operation type          = %s\n", mess_operation_t_str(ops[i_ops]));
            printf("tol                     = %e\n", tol);
            printf("abs. / rel. residual    = %e / %e\n\n", res, relres);

            if(res > tol || relres >tol){
                printf("Failed solvem*:\n");
                printf("Operation type          = %s\n", mess_operation_t_str(ops[i_ops]));
                printf("right hand side;\n");
                mess_matrix_printinfo(B);
                printf("Hessenberg Matrix:\n");
                mess_matrix_printinfo(H);
                printf("tol                     = %e\n", tol);
                printf("abs. / rel. residual    = %e / %e\n", res, relres);
                ret = 1;
                goto clear;
            }

        }

    }

    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/
clear:
    MESS_CLEAR_MATRICES(&H, &B, &X);
    MESS_CLEAR_VECTORS(&x,&b);
    MESS_CLEAR_DIRECTS(&hessenberg);

    return ret;
}
