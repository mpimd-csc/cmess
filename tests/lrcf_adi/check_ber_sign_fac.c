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
 *
 * @addtogroup test_lrcfadi
 * @{
 * @file tests/lrcf_adi/check_ber_sign_fac.c
 * @brief Check the computation of the solution of Bernoulli equations using matrix sign function and its residual.
 * @author @mbehr
 *
 * This function checks the solution of the generalized Bernoulli equation
 * \f[ AX E^T + EXA^T - E XBB^TXE^T = 0 \f]
 * and the Bernoulli equations
 * \f[ A^TX E + E^TXA - E^T XBB^TXE = 0 \f]
 * by computing the corresponding resiudal. <br>
 * Time is measured for computing the solution of one Bernoulli equation and its residual.
 *
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "../call_macro.h"

int  main(int argc , char **argv) {
    mess_version();
    mess_init();
    mess_error_level = 3;

    double tol = sqrt(mess_eps()),err_tol = 1e-5, nrm, nrmB;
    mess_matrix A, B, E, Z, X, temp1, temp2, temp3, temp4;
    int ret = 0;
    mess_int_t err = 0, maxit=50,j;


    if ( argc != 4 ) {
        printf("Usage: %s A.mtx B.mtx E.mtx \n", argv[0]);
        return MESS_ERROR_ARGUMENTS;
    }
    /*-----------------------------------------------------------------------------
     *  init/read matrice and add instable blocks
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&B,&E,&Z,&X,&temp1,&temp2,&temp3,&temp4);
    CALL(mess_matrix_read_formated(argv[1],A,MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[2],B,MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[3],E,MESS_DENSE));
    CALL(mess_matrix_dynormf(B,&nrmB));
    printf("nrmB = %lg\n", nrmB);

    /*-----------------------------------------------------------------------------
     *  compute solution factor Z of A'X + XA - XBB'X = 0
     *-----------------------------------------------------------------------------*/

    for ( j=0;j<3;++j ) {
        tol = sqrt(mess_eps());
        switch ( j ) {
            case 0:
                CALL(mess_lrcf_signber(A,B,&maxit,&tol,Z,MESS_SIGN_SCALE_NONE));
                printf("No Scaling.\n");
                break;
            case 1:
                printf ( "Frobenius Scaling\n" );
                CALL(mess_lrcf_signber(A,B,&maxit,&tol,Z,MESS_SIGN_SCALE_FRO));
                break;
            case 2:
                printf ( "Determinant Scaling\n" );
                CALL(mess_lrcf_signber(A,B,&maxit,&tol,Z,MESS_SIGN_SCALE_DET));
                break;
            default:
                break;
        }
        CALL( mess_matrix_multiply(MESS_OP_NONE,Z,MESS_OP_TRANSPOSE,Z,X));
        CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE,A,MESS_OP_NONE,X,temp1));
        CALL( mess_matrix_multiply(MESS_OP_NONE,X,MESS_OP_NONE,A,temp2));
        CALL( mess_matrix_multiply(MESS_OP_NONE,X,MESS_OP_NONE,B,temp3));
        CALL( mess_matrix_multiply(MESS_OP_NONE,temp3,MESS_OP_TRANSPOSE,B,temp4));
        CALL( mess_matrix_multiply(MESS_OP_NONE,temp4,MESS_OP_NONE,X,temp3));
        CALL( mess_matrix_add(1.0,temp1,1.0,temp2));
        CALL( mess_matrix_add(1.0,temp2,-1.0,temp3));
        CALL( mess_matrix_normf(temp3,&nrm));

        printf("||A'X+XA-XBB'X||_F             =%e\n",nrm);
        printf("||A'X+XA-XBB'X||_F/||BB'||_F   =%e\n",nrm/nrmB);
        maxit = 50;
        tol = sqrt(mess_eps());
        if (nrm>err_tol){
            printf("Test failed\n");
            err++;
            return err;
        }
    }

    /*-----------------------------------------------------------------------------
     *  compute solution factor Z of A'XE + E'XA - E'XBB'XE = 0
     *-----------------------------------------------------------------------------*/

    for ( j=0;j<3;++j ) {
        tol = sqrt(mess_eps());
        switch ( j ) {
            case 1:
                printf("No Scaling.\n");
                CALL(mess_lrcf_gsignber(A,E,B,&maxit,&tol,Z,MESS_SIGN_SCALE_NONE));
                break;
            case 2:
                printf ( "Frobenius Scaling\n" );
                CALL(mess_lrcf_gsignber(A,E,B,&maxit,&tol,Z,MESS_SIGN_SCALE_FRO));
                break;
            case 0:
                printf ( "Determinant Scaling\n" );
                CALL(mess_lrcf_gsignber(A,E,B,&maxit,&tol,Z,MESS_SIGN_SCALE_DET));
                break;
            default:
                break;
        }
        CALL( mess_matrix_multiply(MESS_OP_NONE,Z,MESS_OP_TRANSPOSE,Z,X));
        CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE,A,MESS_OP_NONE,X,temp4));
        CALL( mess_matrix_multiply(MESS_OP_NONE,temp4,MESS_OP_NONE,E,temp1));

        CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE,E,MESS_OP_NONE,X,temp4));
        CALL( mess_matrix_multiply(MESS_OP_NONE,temp4,MESS_OP_NONE,A,temp2));

        CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE,E,MESS_OP_NONE,X,temp4));
        CALL( mess_matrix_multiply(MESS_OP_NONE,temp4,MESS_OP_NONE,B,temp3));
        CALL( mess_matrix_multiply(MESS_OP_NONE,temp3,MESS_OP_TRANSPOSE,B,temp4));
        CALL( mess_matrix_multiply(MESS_OP_NONE,temp4,MESS_OP_NONE,X,temp3));
        CALL( mess_matrix_multiply(MESS_OP_NONE,temp3,MESS_OP_NONE,E,temp4));
        CALL( mess_matrix_copy(temp4,temp3));

        CALL( mess_matrix_add(1.0,temp1,1.0,temp2));
        CALL( mess_matrix_add(1.0,temp2,-1.0,temp3));
        CALL( mess_matrix_normf(temp3,&nrm));

        printf("||A'XE+E'XA-E'XBB'XE||_F             =%e\n",nrm);
        printf("||A'XE+E'XA-E'XBB'XE||_F/||BB'||_F   =%e\n",nrm/nrmB);

        maxit = 50;
        tol = sqrt(mess_eps());
        if (nrm>err_tol){
            printf("Test failed\n");
            err++;
            return err;
        }
    }


    /*-----------------------------------------------------------------------------
     * clear
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&B,&E,&Z,&X,&temp1,&temp2,&temp3,&temp4);
    mess_exit();
    return (err);
}




