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
 * @file tests/matrix/check_decomp.c
 * @brief Check the @ref mess_matrix_decomp function.
 * @author @mbehr
 * @test
 * This function checks the @ref mess_matrix_decomp function defined in basic/decomp.c that means it checks
 * if the identity \f$ A = A_{sym}+ A_{skewsym}\f$ holds.
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
    mess_error_level = 2;
    mess_version();
    mess_init();
    int err = 0, n = 100,i=0,j=0,ret;
    double p = 0.5, eps = sqrt(mess_eps()), abs_err=0;
    mess_matrix A, Asym, Askewsym;
    mess_datatype_t dt;
    mess_storage_t st;

    /*-----------------------------------------------------------------------------
     *  init, create and test
     *-----------------------------------------------------------------------------*/
    for(i=0;i<3;++i){
        for(j=0;j<2;++j){
            switch(j){
                case 0:
                    dt = MESS_REAL;
                    break;
                default:
                    dt = MESS_COMPLEX;
                    break;
            }

            switch(i){
                case 0:
                    st = MESS_DENSE;
                    break;
                case 1:
                    st = MESS_CSC;
                    break;
                default:
                    st = MESS_CSR;
                    break;
            }
            MESS_INIT_MATRICES(&A,&Asym,&Askewsym);
            CALL(mess_matrix_rand(A,n,n,st,dt,p));
            CALL(mess_matrix_decomp(A,Asym,Askewsym));
            CALL(mess_matrix_add(1.0,Asym,1.0,Askewsym));
            CALL(mess_matrix_diffnorm(Askewsym,A,&abs_err));
            MESS_CLEAR_MATRICES(&A,&Asym,&Askewsym);

            printf("absolute error=%e\n",abs_err);
            if(abs_err > eps){
                printf("FAILED with:%s\t%s\n",mess_datatype_t_str(dt),mess_storage_t_str(st));
                return ++err;
            }
        }
    }


    return err;
}

