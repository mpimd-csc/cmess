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
 * @file tests/io/check_vectorio.c
 * @brief Check IO operations on @ref mess_vector instances.
 * @author @koehlerm
 * @test
 *
 * This function checks the @ref mess_vector_read and @ref mess_vector_write functions defined in vector_read.c that means it
 * checks if the reading and writing functions are working for vectors.
 *
 * @}
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "../call_macro.h"
#include "mess/mess.h"


int main (int argc, char ** argv) {
    mess_init();
    int ret =0;
    mess_vector a, b;
    double diffnorn, nrma;
    char *fn ;
    if ( argc != 2 )
        fn = "test.mtx";
    else
        fn = argv[1];


    MESS_INIT_VECTORS(&a,&b);
    CALL(mess_vector_alloc(a, 2000, MESS_REAL));
    CALL(mess_vector_alloc(b, 2000, MESS_REAL));
    CALL(mess_vector_rand(a));

    printf("Random vector:\n");
    CALL(mess_vector_print(a));
    CALL(mess_vector_write(fn, a ));
    CALL(mess_vector_read(fn, b));

    printf("Read Vector:\n");
    CALL(mess_vector_print(b));
    CALL(mess_vector_norm2(a,&nrma));
    CALL(mess_vector_diffnorm(a,b,&diffnorn));
    printf("a= %lg, diffnrm=%lg \n", nrma, diffnorn);

    CALL(mess_vector_clear(&a));
    CALL(mess_vector_clear(&b));

    if ( diffnorn < mess_eps() * nrma) return 0; else return 1;
}

