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
 * @file tests/matrix/check_rank.c
 * @brief Check rank estimation.
 * @author  @mbehr
 * @test
 *
 * This function checks the @ref mess_matrix_rank function defined in rank.c that means it checks if the rank of a matrix
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
    int ret=0, err=0;
    int  sol;
    mess_int_t rank;

    mess_matrix A;

    if(argc!=3){
        printf("Usage: check_rank MATRIX.mtx rank\n");
        ++err;
        return (err>0)?(1):(0);
    }

    //init matrix
    CALL(mess_matrix_init(&A));

    //read matrix
    CALL(mess_matrix_read(argv[1],A));

    //read rank
    sol = atoi(argv[2]);

    CALL(mess_matrix_rank(A,&rank));

    if(rank!=sol){
        printf("Failed with %s correct rank=%d \t computed rank="MESS_PRINTF_INT"\n",argv[1],sol,rank);
        ++err;
    }

    CALL(mess_matrix_clear(&A));

    return (err>0)?(1):(0);
}
