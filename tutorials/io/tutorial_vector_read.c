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
 *
 * @file tutorials/io/tutorial_vector_read.c
 * @brief Demonstrate a simple input/ouput example for vectors (generate random vectors, read, write, print).
 *
 * @sa mess_vector_rand
 * @sa mess_vector_write
 * @sa mess_vector_read
 * @sa mess_vector_print
 *
 * # Tutorial: Read, Write and Generate a Random Vector
 *
 * This function demonstrates how to generate a random vector, to write it to a file, to read this vector from a file
 * and to print it.
 *
 * @snippet "tutorials/io/tutorial_vector_read.c" CODE
 *
 */

///@cond
///[CODE]

#include <stdio.h>
#include <stdlib.h>
#include "mess/mess.h"
#include <complex.h>

int main (int argc, char ** argv) {
    mess_vector a, b;
    mess_vector_init(&a);
    mess_vector_init(&b);
    mess_vector_alloc(a, 10, MESS_REAL);
    mess_vector_alloc(b, 10, MESS_REAL);
    mess_vector_rand(a);
    printf("Random vector:\n");
    mess_vector_print(a);
    mess_vector_write("test.mtx", a );
    mess_vector_read("test.mtx", b);
    printf("Read Vector:\n");
    mess_vector_print(b);
    mess_vector_clear(&a);
    mess_vector_clear(&b);
    return 0;
}
///[CODE]
///@endcond
