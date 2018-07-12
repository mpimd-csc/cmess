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
 * @file tests/eigen/check_schur.c
 * @brief Check the @ref mess_eigen_schur and @ref mess_eigen_schur_complex functions.
 * @author @mbehr
 * @test
 *
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


int main (int argc, char ** argv){
    int ret = 0;
    mess_int_t seed = 1234;
    mess_int_t dim;
    double tol = sqrt(mess_eps()), err;
    mess_vector EV1,EV2;
    mess_matrix A,T1,U1,T2,T3,U3,T4,Temp1,Temp2,Eye;
    mess_error_level = 2;
    mess_int_t trial=0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    if ( argc != 2 ) {
        printf("Usage: %s dim \n", argv[0]);
        return 1;
    }

    dim = atoi(argv[1]);
    /*-----------------------------------------------------------------------------
     *  init/create matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&T1,&U1,&T2,&T3,&U3,&T4,&Temp1,&Temp2,&Eye);
    MESS_INIT_VECTORS(&EV1,&EV2);

    CALL(mess_matrix_rand_init(&seed));
    CALL(mess_matrix_rand(A,dim,dim,MESS_DENSE,MESS_REAL,1.0));

    CALL(mess_matrix_eye(Eye,dim,dim,MESS_DENSE));

    /*-----------------------------------------------------------------------------
     * call mess_eigen_schur / mess_eigen_schur_complex with different possible calling sequences
     *-----------------------------------------------------------------------------*/
    for(trial=0;trial<5;++trial){
        if(trial==0){
            CALL(mess_eigen_schur(A,T1,U1,EV1));
            CALL(mess_eigen_schur(A,T2,NULL,EV2));
            CALL(mess_eigen_schur(A,T3,U3,NULL));
            CALL(mess_eigen_schur(A,T4,NULL,NULL));
        }else if (trial==1){
            CALL(mess_eigen_schur_complex_post(A,T1,U1,EV1));
            CALL(mess_eigen_schur_complex_post(A,T2,NULL,EV2));
            CALL(mess_eigen_schur_complex_post(A,T3,U3,NULL));
            CALL(mess_eigen_schur_complex_post(A,T4,NULL,NULL));
        }else{
            CALL(mess_matrix_toreal(A));
            if(trial==4){
                CALL(mess_matrix_tocomplex(A));
                CALL(mess_matrix_scalec(1+1*I,A));
            }
            CALL(mess_eigen_schur_complex(A,T1,U1,EV1));
            CALL(mess_eigen_schur_complex(A,T2,NULL,EV2));
            CALL(mess_eigen_schur_complex(A,T3,U3,NULL));
            CALL(mess_eigen_schur_complex(A,T4,NULL,NULL));
        }
        /*-----------------------------------------------------------------------------
         *  check first call of mess_eigen_schur
         *-----------------------------------------------------------------------------*/
        CALL(mess_matrix_multiply(MESS_OP_HERMITIAN,U1,MESS_OP_NONE,U1,Temp1));
        CALL(mess_matrix_diffnorm(Temp1,Eye,&err));
        if(err > tol){printf("Failed:\n"); printf("|| U1'*U1 - I ||_2    = %e\n",err); return 1;}

        CALL(mess_matrix_multiply(MESS_OP_HERMITIAN,U1,MESS_OP_NONE,A,Temp1));
        CALL(mess_matrix_multiply(MESS_OP_NONE,Temp1,MESS_OP_NONE,U1,Temp2));
        CALL(mess_matrix_diffnorm(Temp2,T1,&err));
        if(err > tol){printf("Failed:\n"); printf("|| U1'*A*U1 - T1 ||_2  = %e\n",err); return 1;}

        /*-----------------------------------------------------------------------------
         *  check second call of eigenschur
         *-----------------------------------------------------------------------------*/
        CALL(mess_vector_diffnorm(EV1,EV2,&err));
        if(err > tol){printf("Failed:\n"); printf("|| EV1 - EV2 ||_2  = %e\n",err); return 1;}

        CALL(mess_matrix_diffnorm(T1,T2,&err));
        if(err > tol){printf("Failed:\n"); printf("|| T1 - T2 ||_2   = %e\n",err); return 1;}

        /*-----------------------------------------------------------------------------
         *  check third call of eigenschur
         *-----------------------------------------------------------------------------*/
        CALL(mess_matrix_diffnorm(T1,T3,&err));
        if(err > tol){printf("Failed:\n"); printf("|| T1 - T3 ||_2   = %e\n",err); return 1;}

        CALL(mess_matrix_multiply(MESS_OP_HERMITIAN,U3,MESS_OP_NONE,A,Temp1));
        CALL(mess_matrix_multiply(MESS_OP_NONE,Temp1,MESS_OP_NONE,U3,Temp2));
        CALL(mess_matrix_diffnorm(Temp2,T3,&err));
        if(err > tol){printf("Failed:\n"); printf("|| U1'*A*U1 - T1 ||_2  = %e\n",err); return 1;}

        /*-----------------------------------------------------------------------------
         *  check fourth call of eigenschur
         *-----------------------------------------------------------------------------*/
        CALL(mess_matrix_diffnorm(T1,T4,&err));
        if(err > tol){printf("Failed:\n"); printf("|| T1 - T3 ||_2   = %e\n",err); return 1;}
    }

    /*-----------------------------------------------------------------------------
     * clear
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&T1,&U1,&T2,&T3,&U3,&T4,&Temp1,&Temp2,&Eye);
    MESS_CLEAR_VECTORS(&EV1,&EV2);
    return 0;
}


