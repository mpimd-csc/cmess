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
 * @file tests/direct/check_dgeqp3.c
 * @brief Check the @ref mess_direct_dgeqp3 function.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_direct_dgeqp3 function defined in dgeqp3.c that means it checks
 * if the decomposition results in an orthonormal matrix and checks the correctness of the decomposition.
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

    double eps = mess_eps(), diff=0;
    int ret = 0, err = 0;
    mess_matrix A, Q, R, temp, eye;
    mess_int_t * perm;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if (argc !=2) {
        printf("usage: %s [mat]\n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  init and read matrix matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A, &Q, &R, &temp, &eye);
    CALL(mess_matrix_read_formated(argv[1], A, MESS_DENSE));
    perm = (mess_int_t * )  malloc (sizeof(mess_int_t)*A->cols);

    /*-----------------------------------------------------------------------------
     *  compute decomposition
     *-----------------------------------------------------------------------------*/
    CALL(mess_direct_dgeqp3(A, Q, R, perm));

    /*-----------------------------------------------------------------------------
     *  check decomposition
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_multiply(MESS_OP_NONE, Q, MESS_OP_NONE, R, temp));
    CALL(mess_matrix_perm(A, NULL, perm));
    CALL(mess_matrix_diffnormf(A,temp,&diff));
    if(diff>eps){
        printf("Failed:\n");
        mess_matrix_printinfo(A);
        printf("||AP-QR||_F=%e\n",diff);
    }

    /*-----------------------------------------------------------------------------
     *  check orthogonal
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_eye(eye,A->rows, A->cols, MESS_DENSE));
    CALL(mess_matrix_multiply(MESS_OP_TRANSPOSE, Q, MESS_OP_NONE, Q, temp));
    CALL(mess_matrix_diffnormf(temp,eye,&diff));
    if(diff>eps){
        printf("Failed:\n");
        printf("||QQ^T - eye(n,n)||_F = %e \n",diff);
    }

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&Q, &R, &temp, &eye);
    free(perm);

    return (err>0)?(1):(0);
}

