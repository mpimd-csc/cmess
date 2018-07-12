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
 * @file tests/matrix/check_one_value.c
 * @brief Check the @ref mess_matrix_one_value and @ref mess_matrix_one_valuec functions.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_one_value and @ref mess_matrix_one_valuec functions defined in ones.c that means it uses
 * @ref mess_matrix_map  to compare the results.
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


#define VAL 2.1

static double f_one_value(double val){return  VAL;}
static mess_double_cpx_t f_one_valuec(mess_double_cpx_t valc){return VAL;}

int main ( int argc, char ** argv) {

    mess_init();
    double eps=mess_eps(), diff, p = 0.4;
    int ret=0, err=0;
    mess_matrix mat, test, test_dense, temp;
    mess_storage_t  storage_types [] = {MESS_CSR, MESS_CSC, MESS_COORD,MESS_DENSE};
    mess_datatype_t data_types [] = {MESS_REAL, MESS_COMPLEX};
    mess_int_t st=0, dt=0;
    mess_int_t rows, cols;

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
     *  generate matrix with some sparsity pattern
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&mat,&test, &test_dense, &temp);
    CALL(mess_matrix_rand(mat, rows, cols,  MESS_CSR, MESS_REAL,p));

    //make sure that test has the same sparsitiy pattern
    CALL(mess_matrix_copy(mat,test));

    //generate a dense matrix because one value works only on the sparsity pattern
    CALL(mess_matrix_alloc(test_dense, rows, cols, rows*cols, MESS_DENSE, MESS_REAL));

    /*-----------------------------------------------------------------------------
     *  generate sparse/dense  matrix with correct value for comparision
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_map(test,f_one_value,f_one_valuec,0));
    CALL(mess_matrix_map(test_dense,f_one_value,f_one_valuec,0));

    /*-----------------------------------------------------------------------------
     *  check different storage and datatypes
     *-----------------------------------------------------------------------------*/

    for(st=0;st<3;++st){
        for(dt=0;dt<2;++dt){

            /*-----------------------------------------------------------------------------
             *  conver to datatype
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_totype(mat,data_types[dt]));

            /*-----------------------------------------------------------------------------
             *  convert to storage type
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_convert(mat, temp, storage_types[st]));
            CALL(mess_matrix_copy(temp,mat));

            /*-----------------------------------------------------------------------------
             *  apply function mess_matrix_one_value
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_one_value(mat, VAL));

            /*-----------------------------------------------------------------------------
             *  check correctness
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_diffnormf(mat, (storage_types[st]==MESS_DENSE)?test_dense:test, &diff));
            if(diff > eps){
                printf("FAILED mess_matrix_one_value:\n");
                mess_matrix_printinfo(mat); mess_matrix_print(mat);
                mess_matrix_printinfo((storage_types[st]==MESS_DENSE)?test_dense:test);
                mess_matrix_print((storage_types[st]==MESS_DENSE)?test_dense:test);
                printf("diff=%e\n",diff);
                return 1;
            }

            /*-----------------------------------------------------------------------------
             *  apply function mess_matrix_one_value
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_one_valuec(mat, VAL));

            /*-----------------------------------------------------------------------------
             *  check correctness
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_diffnormf(mat, (storage_types[st]==MESS_DENSE)?test_dense:test, &diff));
            if(diff > eps){
                printf("FAILED mess_matrix_one_valuec:\n");
                mess_matrix_printinfo(mat); mess_matrix_print(mat);
                mess_matrix_printinfo((storage_types[st]==MESS_DENSE)?test_dense:test); mess_matrix_print((storage_types[st]==MESS_DENSE)?test_dense:test);
                printf("diff=%e\n",diff);
                return 1;
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&test,&test_dense, &temp, &mat);

    return (err>0)?(1):(0);
}

