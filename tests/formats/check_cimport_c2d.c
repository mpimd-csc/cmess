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
 * @addtogroup test_format
 * @{
 * @file tests/formats/check_cimport_c2d.c
 * @brief Check the creation of a dense mess_matrix from a C \f$ 2 \f$ dimensional array.
 * @author @koehlerm
 * @test
 * This function checks the @ref mess_matrix_dense_from_carray function defined in cimport.c
 * that means it checks if a @ref MESS_DENSE @ref mess_matrix structure is created correctly.
 *
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mess/mess.h"
#include "../call_macro.h"

/**
 * @brief Check the creation of a dense matrix from a real C \f$ 2 \f$ dimensional array.
 * @return zero on success or a non zero error code
 *
 * The @ref check_real function creates a dense mess_matrix from a given real C \f$ 2 \f$ dimensional array and checks
 * if this matrix is in the desired storage format.
 */
int check_real() {
    int ret,err;
    mess_matrix A;
    double **data;
    mess_int_t i,j;
    data = (double**) malloc(sizeof(double*)*100);
    for( i= 0 ; i < 100; i++){
        data[i] = malloc(sizeof(double)*100);
    }
    for ( i = 0; i < 10000 ; i++){
        data[i/100][i%100] = rand();
    }
    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_dense_from_carray(A, 100,100,data, NULL));
    err = 0;

    for ( j = 0; j < A->cols; j++ ) {
        for ( i = 0; i < A->rows; i++ ) {
            if ( data[i][j] != A->values[i+j*A->ld] ) err++;
        }
    }
    if ( err ) {
        fprintf(stderr, "Real Matrix not converted correctly err=%d\n",err);
    }
    mess_matrix_clear(&A);
    for ( i = 0; i < 100;i++){
        free(data[i]);
    }
    free(data);
    return err;
}


/**
 * @brief Check the creation of a dense matrix from a complex C \f$ 2 \f$ dimensional array.
 * @return zero on success or a non zero error code
 *
 * The @ref check_complex function creates a dense mess_matrix from a given complex C \f$ 2 \f$ dimensional array and checks
 * if this matrix is in the desired storage format.
 */
int check_complex() {
    int ret,err;
    mess_matrix A;
    mess_int_t i,j;
    mess_double_cpx_t **data;
    data = (mess_double_cpx_t**) malloc(sizeof(mess_double_cpx_t *)*100);
    for( i= 0 ; i < 100; i++){
        data[i] = malloc(sizeof(mess_double_cpx_t)*100);
    }

    for ( i = 0; i < 10000 ; i++){
        data[i/100][i%100] = rand() + rand()*I;
    }
    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_dense_from_carray(A, 100,100,NULL, data));
    err = 0;
    for ( j = 0; j < A->cols; j++ ) {
        for ( i = 0; i < A->rows; i++ ) {
            if ( data[i][j] != A->values_cpx[i+j*A->ld] ) err++;
        }
    }

    if ( err ) {
        fprintf(stderr, "Complex Matrix not converted correctly err=%d\n",err);
    }
    mess_matrix_clear(&A);
    for ( i = 0; i < 100;i++){
        free(data[i]);
    }
    free(data);

    return err;
}


int main( int argc, char **argv) {
    mess_init();
    int ret;
    srand(time(NULL));

    CALL(check_real());
    CALL(check_complex());


    return ret;
}

