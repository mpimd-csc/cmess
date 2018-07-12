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
 * @file tests/eigen/check_gsign.c
 * @brief Check the @ref mess_eigen_gsign function.
 * @author @mbehr
 * @test
 *
 * This function checks the @ref mess_eigen_gsign function defined in sign.c that means it checks if the sign of
 * a matrix is computed correctly, by the identity
 * \f[ (E^{-1}sign(A,B))^2 =  I_n\f]
 * is computed correctly.
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
    mess_version();
    mess_error_level = 3;
    double tol = sqrt(mess_eps()), abs_err;
    mess_matrix A,E,Z,tmp,tmp2,tmp3,eye;
    mess_int_t j, k;
    int ret = 0, err= 0;

    if ( argc != 3 ) {
        printf("Usage: %s A.mtx E.mtx\n", argv[0]);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  init/read matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&E,&Z,&tmp,&tmp2,&tmp3,&eye);
    CALL(mess_matrix_read_formated(argv[1],A,MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[2],E,MESS_DENSE));
    CALL(mess_matrix_eye(eye,A->rows,A->rows,MESS_DENSE));

    /*-----------------------------------------------------------------------------
     *  compute sign of matrix
     *-----------------------------------------------------------------------------*/
    for( k=0;k<2;++k){
        if(k==1){
            CALL(mess_matrix_tocomplex(A));
            CALL(mess_matrix_tocomplex(E));
            CALL(mess_matrix_tocomplex(eye));
            CALL(mess_matrix_scalec(1.0/sqrt(2)*(1+1*I),A));
            CALL(mess_matrix_scalec(1.0/sqrt(2)*(1+1*I),E));
        }
        for ( j=0;j<3;++j ) {
            tol = sqrt(mess_eps());
            switch ( j ) {
                case 0:
                    printf("No Scaling.\n");
                    CALL(mess_eigen_gsign(A,E,Z,MESS_SIGN_SCALE_NONE));
                    break;
                case 1:
                    printf ( "Frobenius Scaling\n" );
                    CALL(mess_eigen_gsign(A,E,Z,MESS_SIGN_SCALE_FRO));
                    break;
                case 2:
                    printf ( "Determinant Scaling\n" );
                    CALL(mess_eigen_gsign(A,E,Z,MESS_SIGN_SCALE_DET));
                    break;
                default:
                    break;
            }

            //compute error
            CALL(mess_matrix_copy(Z,tmp));
            CALL(mess_matrix_backslashm(MESS_OP_NONE,E,Z,tmp2));
            CALL(mess_matrix_copy(tmp2,tmp));
            CALL(mess_matrix_multiply(MESS_OP_NONE,tmp,MESS_OP_NONE,tmp2,tmp3));
            CALL(mess_matrix_diffnorm(tmp3,eye,&abs_err));
            printf("abs_err=%e\n",abs_err);
            if (abs_err>tol){
                CALL(mess_matrix_printinfo(A));
                CALL(mess_matrix_printinfo(E));
                printf("Test failed\n");
                printf("abs_err=%e\n",abs_err);
                return ++err;
            }
        }
    }

    /*-----------------------------------------------------------------------------
     * clear
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&Z,&tmp,&tmp2,&E,&eye,&tmp3);
    mess_exit();
    return err;
}

