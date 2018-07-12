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
 * @file tests/direct/check_direct_care.c
 * @brief Check the @ref mess_direct_care and @ref mess_direct_care_res2 function.
 * @author @mbehr
 * @test
 *
 * Checks the @ref mess_direct_care and @ref mess_direct_care_res2 function.
 *
 * @}
 *
 */
#include "../call_macro.h"
#include <stdio.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"



int main ( int argc, char ** argv) {
    mess_init();
    mess_error_level=2;
    int ret, err=0;
    double tol = 1e-4, res2_1, res2_2;
    mess_matrix A, E, B, C, Q, G, X1, X2, TMP;

    /*-----------------------------------------------------------------------------
     * check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 5) {
        printf("usage: %s A.mtx E.mtx B.mtx C.mtx\n", argv[0]);
        return 1;
    }

    /*----------------------------------------------------------------------------
     *  init and read matriceas
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&E,&B,&C,&Q,&G,&X1,&X2,&TMP);
    CALL(mess_matrix_read_formated(argv[1], A, MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[2], E, MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[3], B, MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[4], C, MESS_DENSE));

    /*-----------------------------------------------------------------------------
     *  compute Q and G
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_multiply(MESS_OP_NONE, B, MESS_OP_TRANSPOSE, B, Q));
    CALL(mess_matrix_multiply(MESS_OP_TRANSPOSE, C, MESS_OP_NONE, C, G));


    /*-----------------------------------------------------------------------------
     *  solve riccati equation and compute residuals
     *-----------------------------------------------------------------------------*/
    CALL(mess_direct_care(A, NULL, Q, G, X1));
    CALL(mess_direct_care_res2(A, NULL, Q, G, X1, &res2_1));
    CALL(mess_direct_care(A, E, Q, G, X2));
    CALL(mess_direct_care_res2(A, E, Q, G, X2, &res2_2));

    printf("res2_1 = %e,    tolerance = %e\n", res2_1, tol);
    if(res2_1>tol){
        err++;
    }

    printf("res2_2 = %e,    tolerance = %e\n", res2_2, tol);
    if(res2_2>tol){
        err++;
    }


    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&E,&B,&C,&Q,&G,&X1,&X2,&TMP);

    return err;
}

