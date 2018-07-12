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
 * @file tests/matrix/check_getelement.c
 * @brief Check if a matrix element is returned correctly.
 * @test
 * This function checks the @ref mess_matrix_getelement element function defined in getelement.c that means it checks if
 * a wanted matrix element is returned correctly.
 *
 * @}
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/mess.h"

int main (int argc, char *argv[])
{
    mess_init();
    mess_matrix mat_coord, mat_csr;
    mess_int_t i, j;
    double val1, val2;
    int err = 0 ;

    mess_matrix_init(&mat_coord);
    mess_matrix_init(&mat_csr);

    if ( argc != 2) {
        printf("compute the row sums of a matrix\n");
        printf("usage: %s matrix.mtx\n", argv[0]);
    }
    mess_matrix_read(argv[1], mat_coord);
    mess_matrix_convert(mat_coord, mat_csr, MESS_CSR);
    for ( i = 0 ; i < mat_csr->rows; i++) {
        for (j= 0; j < mat_csr->cols; j++) {
            mess_matrix_getelement(mat_csr, i, j , &val1, NULL);
            mess_matrix_getelement(mat_coord, i, j, &val2, NULL);
            if ( val1 != val2) err++;
        }
    }
    mess_matrix_clear(&mat_coord);
    mess_matrix_clear(&mat_csr);

    return (err==0)?0:1;
}

