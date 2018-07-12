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
 * @file tests/matrix/check_bandwidth.c
 * @brief Check @ref mess_matrix_bandwidth function.
 * @author @koehlerm
 * @test
 *
 * This function checks if @ref mess_matrix_bandwidth works correctly.
 *
 * @}
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mess/mess.h"
#include "../call_macro.h"

int main ( int argc, char **argv){
    mess_init();
    mess_matrix A, B ;
    mess_int_t rowptr[6] = {0,2,5,7,11,13} ;
    mess_int_t colptr[13] = {0,1,0,1,2,0,2,1,2,3,4,3,4};
    double values[13] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    mess_int_t kl, ku;
    int ret ;
    int err = 0 ;

    /*-----------------------------------------------------------------------------
     *  read data
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_csr(A, 5, 5, rowptr, colptr, values, NULL));
    CALL(mess_matrix_sort(A));
    CALL(mess_matrix_print(A));
    CALL(mess_matrix_bandwidth(A, &kl, &ku));

    printf("kl = "MESS_PRINTF_INT" \t ku = "MESS_PRINTF_INT" \n", kl, ku);
    if ( kl != 2 && ku != 1) {
        printf("Bandwidth wrong\n");
        err ++;
    }
    CALL(mess_matrix_init(&B));
    CALL(mess_matrix_convert(A, B, MESS_CSC));

    CALL(mess_matrix_bandwidth(B, &kl, &ku));

    printf("kl = "MESS_PRINTF_INT" \t ku = "MESS_PRINTF_INT" \n", kl, ku);
    if ( kl != 2 && ku != 1) {
        printf("Bandwidth wrong\n");
        err ++;
    }


    CALL(mess_matrix_clear(&B));
    CALL(mess_matrix_clear(&A));
    return 0;
}




