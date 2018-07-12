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
 * @addtogroup test_io
 * @{
 * @file tests/io/check_matrixio.c
 * @brief Check IO operations on @ref mess_matrix instances.
 * @author @koehlerm
 * @test
 *
 * This function checks the @ref mess_matrix_read and @ref mess_matrix_write functions defined in write.c and read.c that means it
 * checks if the reading and writing functions are working for vectors.
 *
 * @}
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "../call_macro.h"
#include "mess/mess.h"


int main (int argc, char ** argv) {

    mess_init();
    int ret=0;
    mess_int_t i=0, j=0, k=0, l=0, rows=20, cols=33;
    double eps = sqrt(mess_eps()), diff=0, p=1;
    mess_matrix A, B;
    mess_datatype_t  dts [] = {MESS_REAL,MESS_COMPLEX};
    mess_storage_t sts []  = {MESS_DENSE,MESS_CSR,MESS_CSC,MESS_COORD};


    /*-----------------------------------------------------------------------------
     *  init matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&B);

    for(i=1;i<rows;++i){
        for(j=1;j<cols;++j){
            for(k=0;k<2;++k){
                for(l=0;l<4;++l){

                    /*-----------------------------------------------------------------------------
                     *  generate random matrix and test
                     *-----------------------------------------------------------------------------*/
                    CALL(mess_matrix_rand(A,i,j,sts[l],dts[k],p));
                    CALL(mess_matrix_copy(A,B));
                    CALL(mess_matrix_write("check_matrix_io.mtx",B));
                    MESS_MATRIX_RESET(B);
                    CALL(mess_matrix_read("check_matrix_io.mtx",B));

                    if(A->rows != B->rows || A->cols != B->cols || A->data_type != B->data_type ||
                            (MESS_IS_DENSE(A)  && !MESS_IS_DENSE(B)) ||
                            (!MESS_IS_DENSE(A) && MESS_IS_DENSE(B))
                      ){
                        printf("FAILED\n");
                        mess_matrix_printinfo(A);
                        mess_matrix_printinfo(B);
                        return 1;
                    }



                    CALL(mess_matrix_diffnormf(A,B,&diff));
                    if(diff>eps){
                        printf("FAILED\n");
                        printf("diff=%e\n",diff);
                        mess_matrix_printinfo(A);
                        mess_matrix_printinfo(B);
                        mess_matrix_print(A);
                        mess_matrix_print(B);
                        return 1;
                    }

                    MESS_RESET_MATRICES(A,B);

                }
            }
        }
    }
    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&B);

    return 0;
}

