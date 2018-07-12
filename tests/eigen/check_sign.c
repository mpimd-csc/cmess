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
 * @file tests/eigen/check_sign.c
 * @brief Check the @ref mess_eigen_sign function.
 * @author @mbehr
 * @test
 *
 * This function checks the @ref mess_eigen_sign function defined in sign.c
 * that means it checks if the sign of
 * a matrix is computed correctly, by the identity
 * \f[ sign(A)^2 =  I\f].
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
    mess_matrix A,Z,eye,tmp,tmp2;
    int ret = 0, err = 0;
    mess_int_t j, k;


    if ( argc != 2 ) {
        printf("Usage: %s A.mtx\n", argv[0]);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  init/read matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&Z,&eye,&tmp,&tmp2);
    CALL(mess_matrix_read_formated(argv[1],A,MESS_DENSE));
    CALL(mess_matrix_eye(eye,A->rows,A->rows,MESS_DENSE));

    /*-----------------------------------------------------------------------------
     *  compute sign of matrix
     *-----------------------------------------------------------------------------*/
    for( k=0;k<2;++k){
        if(k==1){CALL(mess_matrix_tocomplex(A)); CALL(mess_matrix_scalec(1.0/sqrt(2)*(1+1*I),A));}
        for ( j=0;j<3;++j ) {
            tol = sqrt(mess_eps());
            switch ( j ) {
                case 0:
                    printf("No Scaling.\n");
                    CALL(mess_eigen_sign(A,Z,MESS_SIGN_SCALE_NONE));
                    break;
                case 1:
                    printf ( "Frobenius Scaling\n" );
                    CALL(mess_eigen_sign(A,Z,MESS_SIGN_SCALE_FRO));
                    break;
                case 2:
                    printf ( "Determinant Scaling\n" );
                    CALL(mess_eigen_sign(A,Z,MESS_SIGN_SCALE_DET));
                    break;
                default:
                    break;
            }

            //compute error
            CALL(mess_matrix_copy(Z,tmp));
            CALL(mess_matrix_multiply(MESS_OP_NONE,Z,MESS_OP_NONE,tmp,tmp2));
            CALL(mess_matrix_diffnorm(tmp2,eye,&abs_err));
            printf("abs_err=%e\n",abs_err);
            if (abs_err>tol){
                CALL(mess_matrix_printinfo(A));
                printf("Test failed\n");
                printf("abs_err=%e\n",abs_err);
                return ++err;
            }
        }
    }
    /*-----------------------------------------------------------------------------
     * clear
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&Z,&tmp,&tmp2,&eye);
    mess_exit();
    return err;
}

