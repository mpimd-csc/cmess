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
 * @file tests/matrix/check_scalem.c
 * @brief Check the @ref mess_matrix_rowscalem and @ref mess_matrix_colscalem.
 *
 * @test
 *
 * This function checks the @ref mess_matrix_rowscalem and @ref mess_matrix_colscalem function defined in scal.c that means  it
 * checks the results of the function against each other.
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


int main ( int argc, char **argv){
    mess_init();
    int ret =0;
    mess_int_t i,j,k, rows, cols;
    double tol = sqrt(mess_eps()), diff=1.0, p = 0.5;
    mess_datatype_t dts [] = {MESS_REAL, MESS_COMPLEX};
    mess_datatype_t dts2 [] = {MESS_REAL, MESS_COMPLEX};
    mess_storage_t sts [] = {MESS_DENSE, MESS_CSR, MESS_CSC, MESS_COORD};
    mess_matrix A, temp1, temp2, temp3;
    mess_vector scale;

    mess_init();
    mess_error_level = 3;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if (argc !=3) {
        printf("usage: %s rows cols\n", argv[0]);
        return 1;
    }

    rows = atoi(argv[1]);
    cols = atoi(argv[2]);

    /*-----------------------------------------------------------------------------
     *  init
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&temp1,&temp2,&temp3);

    CALL(mess_matrix_rand_init(NULL));
    CALL(mess_vector_rand_init(NULL));

    /*-----------------------------------------------------------------------------
     *  test different cases
     *-----------------------------------------------------------------------------*/
    for(i=0;i<2;++i){
        for(j=0;j<2;++j){
            for(k=0;k<4;++k){
                //print data and storage types
                printf("A:%s\t%s\t\tscale:%s\n",mess_datatype_t_str(dts[j]),mess_storage_t_str(sts[k]),mess_datatype_t_str(dts[i]));

                //generate random matrix
                CALL(mess_matrix_rand(A, rows, cols, sts[k], dts[j], p));

                //generate random vector
                MESS_INIT_VECTORS(&scale);
                CALL(mess_vector_alloc(scale,rows,dts2[i]));
                CALL(mess_vector_rand(scale));

                //row scale
                CALL(mess_matrix_copy(A,temp1));
                CALL(mess_matrix_rowscalem(scale,A));

                //transpose, compute col scale, transpose back
                CALL(mess_matrix_ctranspose(temp1,temp2));
                CALL(mess_matrix_conj(temp2));
                CALL(mess_matrix_colscalem(scale,temp2));
                CALL(mess_matrix_ctranspose(temp2,temp3));
                CALL(mess_matrix_conj(temp3));

                //compute result
                CALL(mess_matrix_diffnormf(A,temp3,&diff));

                if(diff>tol){
                    printf("FAILED:\n");
                    printf("diff=%e\n",diff);
                    printf("tol =%e\n",tol);
                    mess_matrix_printinfo(A);
                    mess_vector_printinfo(scale);
                    mess_matrix_print(A);
                    mess_matrix_print(temp3);
                    return 1;
                }

                //reset clear
                MESS_RESET_MATRICES(A);
                //MESS_CLEAR_MATRICES(&A);
                //MESS_INIT_MATRICES(&A);
                MESS_CLEAR_VECTORS(&scale);
            }
        }
    }


    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&temp1,&temp2,&temp3);

    return 0;

}
