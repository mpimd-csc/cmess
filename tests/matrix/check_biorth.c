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
 * @file tests/matrix/check_biorth.c
 * @brief Check the construction of a biorthonormal basis for two matrices.
 * @author @koehlerm
 * @test
 * This function checks if the @ref mess_matrix_biorth function defined in biorth.c works, that means
 * it checks if it constructs for given matrices \f$ V_{in} \f$ and \f$ W_{in} \f$ a biorthonormal basis such that
 * \f[ W_{out}^T V_{out}= I, \f]
 * where denotes the identity matrix of appropriate dimension.
 * @}
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/mess.h"

#include "../call_macro.h"

int main (int argc, char *argv[])
{
    mess_init();
    mess_matrix V, W, Vin, Win, eye,t;
    double diff, tol;

    int ret = 0;
    mess_int_t rows, cols;
    int err = 0;
    mess_int_t seed = 1234;
    double nrmV, nrmW;
    double eps = mess_eps();

    if ( argc != 3 ) {
        printf("usage %s rows cols\n", argv[0]);
        return -1;
    }

    CALL(mess_matrix_init(&V));
    CALL(mess_matrix_init(&W));
    CALL(mess_matrix_init(&Vin));
    CALL(mess_matrix_init(&Win));
    CALL(mess_matrix_init(&eye));
    CALL(mess_matrix_init(&t));

    rows = atoi(argv[1]);
    cols = atoi(argv[2]);

    CALL(mess_matrix_rand_init(&seed));
    CALL(mess_matrix_rand(Vin, rows, cols, MESS_DENSE,MESS_REAL, 1));
    CALL(mess_matrix_rand(Win, rows, cols, MESS_DENSE,MESS_REAL, 1));
    CALL(mess_matrix_norm2(Vin, &nrmV));
    CALL(mess_matrix_norm2(Win, &nrmW));


    CALL(mess_matrix_biorth(Vin, Win, V, W));
    CALL(mess_matrix_eye(eye, cols, cols, MESS_DENSE));
    CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE, W, MESS_OP_NONE, V, t));
    CALL(mess_matrix_diffnorm(t,eye,&diff));
    tol = sqrt(eps);

    err = diff > tol;
    printf("nrmV = %lg\nnrmW = %lg\n",nrmV, nrmW);
    printf("diff = %lg \n", diff);
    printf("tol  = %e\n",tol);

    mess_matrix_clear(&Vin);
    mess_matrix_clear(&Win);
    mess_matrix_clear(&eye);
    mess_matrix_clear(&V);
    mess_matrix_clear(&W);
    mess_matrix_clear(&t);

    return err;
}

