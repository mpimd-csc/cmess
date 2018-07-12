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
 * @addtogroup test_matrix
 * @{
 * @file tests/matrix/check_dynorm_indefinitef.c
 * @brief Check the @ref mess_matrix_indefinite_dynormf and @ref mess_matrix_indefinite_dynormf2 functions.
 * @test
 *
 * This function checks the @ref mess_matrix_indefinite_dynormf and @ref mess_matrix_indefinite_dynormf2 functions defined in dynorm.c.
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

static int test_cases(mess_matrix W, mess_matrix K, double tol){
    double err1, err2, nrm,nrm2, sol;
    int ret =0;
    mess_matrix Temp1, Temp2;
    MESS_INIT_MATRICES(&Temp1,&Temp2);
    CALL(mess_matrix_multiply(MESS_OP_NONE, W, MESS_OP_HERMITIAN, W, Temp1));
    CALL(mess_matrix_multiply(MESS_OP_NONE, K, MESS_OP_HERMITIAN, K, Temp2));
    CALL(mess_matrix_add(1.0,Temp1,-1.0,Temp2));
    CALL(mess_matrix_normf2(Temp2,&sol));
    MESS_CLEAR_MATRICES(&Temp1,&Temp2);

    CALL(mess_matrix_indefinite_dynormf2(W,K,&nrm));
    CALL(mess_matrix_indefinite_dynormf(W,K,&nrm2));
    nrm2*=nrm2;

    err1= fabs(sol-nrm)/fabs(sol);
    err2= fabs(sol-nrm2)/fabs(sol);

    if(err1>tol || err2>tol){
        printf("nrm=%e\tsol=%e\terr=%e\n",nrm,sol,err1);
        printf("nrm2=%e\tsol=%e\terr=%e\n",nrm,sol,err2);
        printf("failed with:\n");
        CALL(mess_matrix_printinfo(W));
        CALL(mess_matrix_printinfo(K));
        CALL(mess_matrix_printshort(W));
        CALL(mess_matrix_printshort(K));
        return 1;
    }else{
        printf("nrm=%e\tsol=%e\terr1=%e\n",nrm,sol,err1);
        printf("nrm=%e\tsol=%e\terr2=%e\n",nrm,sol,err2);
    }
    return 0;

}

int main ( int argc, char **argv){
    mess_init();
    mess_error_level = 3;
    mess_version();
    int ret = 0;
    mess_matrix K, W, Kc, Wc;
    mess_int_t Kcols=10, Wcols=10, Kcol=1, Wcol=1;
    mess_int_t rows=100;
    double tol= sqrt(mess_eps());

    MESS_INIT_MATRICES(&K,&W,&Kc,&Wc);

    for(Kcol=1;Kcol<Kcols;++Kcol){
        for(Wcol=1;Wcol<Wcols;++Wcol){

            /*-----------------------------------------------------------------------------
             *  generate matrices
             *-----------------------------------------------------------------------------*/

            CALL(mess_matrix_rand(K, rows, Kcols, MESS_DENSE, MESS_REAL, 1.0));
            CALL(mess_matrix_rand(W, rows, Wcols, MESS_DENSE, MESS_REAL, 1.0));
            CALL(mess_matrix_scale(13, W));
            CALL(mess_matrix_rand(Kc,rows, Kcols, MESS_DENSE, MESS_COMPLEX, 1.0));
            CALL(mess_matrix_rand(Wc,rows, Wcols, MESS_DENSE, MESS_COMPLEX, 1.0));
            CALL(mess_matrix_scalec(13+12*I,Wc));

            /*-----------------------------------------------------------------------------
             *  check cases
             *-----------------------------------------------------------------------------*/
            if(test_cases(W , K , tol)) return 1;
            if(test_cases(Wc, K , tol)) return 1;
            if(test_cases(W , Kc, tol)) return 1;
            if(test_cases(Wc, Kc, tol)) return 1;

        }
    }

    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/

    MESS_CLEAR_MATRICES(&K,&W,&Kc,&Wc);

    return 0;
}

