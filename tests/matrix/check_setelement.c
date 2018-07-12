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
 * @file tests/matrix/check_setelement.c
 * @brief Check the overwriting of a matrix element.
 * @test
 *
 * This function checks the @ref mess_matrix_setelement and @ref mess_matrix_setelement_complex functions defined in setelement.c
 * that means it checks if an element of a matrix is overwritten correctly with a given element by use of @ref mess_matrix_getelement.
 *
 * @}
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../call_macro.h"
#include "mess/mess.h"

int main (int argc, char *argv[])
{
    mess_init();
    int rows, cols, err=0,i,j,ret;
    double p = 0.5, val_set = 2, val_set2 = 0, val_get=0;
    mess_double_cpx_t val_setc = 2+1*I, val_set2c=0, val_getc = 0;
    mess_matrix mat, matc,temp;
    mess_int_t round, rounds=10;
    mess_storage_t storage_types[] = {MESS_CSR, MESS_COORD, MESS_CSC, MESS_DENSE};
    mess_int_t st;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    if(argc!=3){
        printf("Usage: %s rows cols\n", argv[0]);
    }
    rows = atoi(argv[1]);
    cols = atoi(argv[2]);

    /*-----------------------------------------------------------------------------
     *  generate random matrix
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&mat, &matc, &temp);

    CALL(mess_matrix_rand(mat, rows, cols, MESS_CSR, MESS_REAL, p));
    CALL(mess_matrix_rand(matc, rows, cols, MESS_CSR, MESS_COMPLEX, p));

    /*-----------------------------------------------------------------------------
     *  check different cases
     *-----------------------------------------------------------------------------*/
    for(st=0;st<4; ++st){
        for(round=0;round<rounds;++round){

            /*-----------------------------------------------------------------------------
             *  convert to datatype
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_convert(mat, temp, storage_types[st]));
            CALL(mess_matrix_copy(temp, mat));
            CALL(mess_matrix_convert(matc, temp, storage_types[st]));
            CALL(mess_matrix_copy(temp, matc));

            /*-----------------------------------------------------------------------------
             *  choose random row and column indices
             *-----------------------------------------------------------------------------*/
            i = rand() % rows;
            j = rand() % cols;

            /*-----------------------------------------------------------------------------
             *  check setelement (MESS_REAL)
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_setelement(mat, i, j, val_set));
            CALL(mess_matrix_getelement(mat, i, j, &val_get, NULL));
            if(val_set!=val_get){
                printf("Failed (val_set):\n");
                CALL(mess_matrix_printinfo(mat));
                printf("i=%d\t j=%d\n", i, j);
                return 1;
            }

            CALL(mess_matrix_setelement(mat, i, j, val_set2));
            CALL(mess_matrix_getelement(mat, i, j, &val_get, NULL));
            if(val_set2!=val_get){
                printf("Failed (val_set2):\n");
                CALL(mess_matrix_printinfo(mat));
                printf("i=%d\t j=%d\n", i, j);
                return 1;
            }


            /*-----------------------------------------------------------------------------
             *  check setelement complex (MESS_COMPLEX)
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_setelement_complex(matc, i, j, val_setc));
            CALL(mess_matrix_getelement(matc, i, j, NULL, &val_getc));
            if(val_setc!=val_getc){
                printf("Failed (val_setc):\n");
                CALL(mess_matrix_printinfo(matc));
                printf("i=%d\t j=%d\n", i, j);
                return 1;
            }

            CALL(mess_matrix_setelement_complex(matc, i, j, val_set2c));
            CALL(mess_matrix_getelement(matc, i, j, NULL, &val_getc));
            if(val_set2c!=val_getc){
                printf("Failed (val_set2c):\n");
                CALL(mess_matrix_printinfo(matc));
                printf("i=%d\t j=%d\n", i, j);
                return 1;
            }


        }
    }

    /*-----------------------------------------------------------------------------
     *  clear memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&mat,&matc,&temp);

    return (err>0);
}

