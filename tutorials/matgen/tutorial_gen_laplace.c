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
 * @file tutorials/matgen/tutorial_gen_laplace.c
 * @brief Demonstrate how to discretize a convection-diffusion operator on the unit square.
 *
 * @author @koehlerm
 * @sa mess_matgen_fdmmatrix
 *
 * # Tutorial: Matrix Generators
 *
 * This function demonstrates how a convection-diffusion operator can be discretized on the unit square.
 *
 * @snippet "tutorials/matgen/tutorial_gen_laplace.c" CODE
 *
 */

///@cond
///[CODE]

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mess/mess.h"

int main ( int argc, char ** argv) {
    mess_matrix A;
    mess_int_t n;
    double ts, te;

    if (argc !=3 ) {
        printf("FDM Generator: LAPLACE on Unit Square\n");
        printf("Usage: %s n0 output.mtx\n", argv[0]);
        return -1;
    }
    n = atoi ( argv[1] );
    mess_matrix_init(&A);
    printf("Generate LAPLACE on (0,1) x (0,1)...");
    ts = mess_wtime();
    mess_matgen_fdmmatrix(A, n, NULL, NULL, NULL);
    te = mess_wtime();
    printf("%lgs\n", te-ts);
    printf("Write it to %s ...", argv[2]);
    mess_matrix_write(argv[2], A);
    printf("done.\n");
    mess_matrix_clear(&A);
    return 0;
}

///[CODE]
///@endcond
