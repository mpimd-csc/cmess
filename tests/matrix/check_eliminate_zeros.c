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
 * @file tests/matrix/check_eliminate_zeros.c
 * @brief Check the @ref mess_matrix_eliminate_zeros function.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_eliminate_zeros  function defined in dupl.c that means it checks
 * if all zero entries from the matrix are removed and the matrices are still equal after removing the zero
 * entries from the data structure.
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


static double f_real(double val){ return fabs(val)<0.4? 0.0:val ; }
static mess_double_cpx_t f_cpx(mess_double_cpx_t val){ return cabs(val)<0.4? 0.0:val; }
static mess_int_t iszero_real(double val){ return val==0.0; }
static mess_int_t iszero_cpx(mess_double_cpx_t val){ return val==0.0;}



int main ( int argc, char ** argv) {
    mess_init();

    double eps = mess_eps(), diff=0;
    int ret = 0, err = 0, rows = 1, cols = 4, i=0, has_zero=0;
    mess_int_t seed = 1;
    mess_matrix mat, temp;
    mess_vector test;


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
     *  init matrices and random number generator
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&test);
    mess_vector_alloc(test,0,MESS_REAL);
    CALL(mess_matrix_rand_init(&seed)); //init with calender time

    MESS_INIT_MATRICES(&mat, &temp);

    /*-----------------------------------------------------------------------------
     * iterate over data and storage types
     *-----------------------------------------------------------------------------*/
    mess_storage_t  storage_types [] = {MESS_CSR, MESS_CSC, MESS_COORD};
    mess_datatype_t data_types [] = {MESS_REAL, MESS_COMPLEX};
    mess_int_t st=0, dt=0;


    for(st=0;st<3;++st){
        for(dt=0;dt<2;++dt){

            /*-----------------------------------------------------------------------------
             *  generate random matrix, insert zeros
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_rand(mat, rows, cols, storage_types[st], data_types[dt], 1.0));
            CALL(mess_matrix_map(mat,f_real,f_cpx,0));

            /*-----------------------------------------------------------------------------
             *  eliminate zeros
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_copy(mat,temp));
            CALL(mess_matrix_eliminate_zeros(temp));

            /*-----------------------------------------------------------------------------
             *  check if matrices are still equal
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_diffnormf(mat,temp,&diff));

            if(diff>eps){
                printf("FAILED: diff=%e\n",diff);
                CALL(mess_matrix_printinfo(mat));
                return 1;
            }

            /*-----------------------------------------------------------------------------
             *  check if there is still a zero value
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_any(temp,iszero_real, iszero_cpx,0,test));
            for(i=0;i<test->dim;++i){ if(test->values[i]!=0.0) has_zero=1;};
            if(has_zero){
                printf("Failed: Matrix has still zero entries.\n");
                printf("Before\n");
                CALL(mess_matrix_print(mat));
                printf("After:\n");
                CALL(mess_matrix_print(temp));
                return 1;
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&mat,&temp);
    MESS_CLEAR_VECTORS(&test);

    return (err>0)?(1):(0);
}

