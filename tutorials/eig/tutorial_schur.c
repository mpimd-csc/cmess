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
 * @file tutorials/eig/tutorial_schur.c
 * @brief Demonstrate how to compute the (complex) Schur decomposition of a matrix.
 *
 * @sa mess_eigen_schur
 * @sa mess_eigen_schur_complex
 *
 * # Tutorial: Compute Schur Decomposition
 *
 * This function demonstrates computing the Schur decomposition of a matrices \f$ A \f$
 * \f[ T= U^T A U \f]
 * where \f$ T \f$ is the Schur form of \f$ A \f$ and \f$ U \f$ is an orthogonal matrix
 * and the eigenvalues of the eigenvalue problem:
 * \f[ A x = \lambda x. \f]
 * using @lapack's DGEES subroutine.
 *
 * @snippet "tutorials/eig/tutorial_schur.c" CODE
 *
 *
 */

///@cond
///[CODE]
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/mess.h"
#include <complex.h>

int main ( int argc, char **argv){
    mess_matrix A, T, U;
    mess_vector EV;
    int ret = 0;
    int type = 1;
    mess_version();

    printf("mess schur decomposition demo\n");
    printf("==============================\n");

    if ( argc!=3){
        printf("Usage: %s matrix.mtx type\n", argv[0]);
        printf(" type = 1 real schur form\n");
        printf(" type = 2 complex schur form\n");
        exit(1);
    }
    type = atoi(argv[2]);
    mess_matrix_init(&A);
    mess_matrix_init(&T);
    mess_matrix_init(&U);
    mess_vector_init(&EV);
    mess_vector_alloc(EV, 1, MESS_COMPLEX);

    if ( mess_matrix_read(argv[1], A) != 0){
        printf("error reading matrix\n");
    }
    if ( type == 1) {
        ret = mess_eigen_schur(A,T,U, EV);
    } else if ( type ==2){
        ret = mess_eigen_schur_complex(A,T,U, EV);
    }else{
        ret =1;
    }
    if ( ret !=0) {
        printf("an error occured. ret = %d\n", ret);
        return ret;
    } else{
        printf("T=\n");
        mess_matrix_print(T);
        printf("\nU=\n");
        mess_matrix_print(U);
        printf("\nEigenvalues:\n");
        mess_vector_print(EV);
    }

    mess_matrix_clear(&A);
    mess_matrix_clear(&T);
    mess_matrix_clear(&U);
    mess_vector_clear(&EV);
    return 0;
}
///[CODE]
///@endcond

