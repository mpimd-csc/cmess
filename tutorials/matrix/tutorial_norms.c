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
 * @file tutorials/matrix/tutorial_norms.c
 * @brief Demonstrate how to compute different matrix norms with different storage formats.
 * @author @mbehr
 * @author @koehlerm
 *
 * # Tutorial for Computing Different Norms of a Matrix
 * This function demonstrates computing different matrix norms for matrix \f$ A \f$
 * <ul>
 * <li> \f$ 1 \f$-norm (column-sum-norm)
 * <li> \f$ 2 \f$-norm (largest eigenvalue of \f$ A^HA \f$)
 * <li> \f$ \infty \f$-norm (row-sum-norm)
 * <li> Frobenius-norm
 * </ul>
 * with different storage formats (CSC, CSC, COORD).\n
 * Additionally the rank of a matrix is computed for CSR storage format.
 *
 * @sa mess_matrix_norm1
 * @sa mess_matrix_norm2
 * @sa mess_matrix_norminf
 * @sa mess_matrix_normf
 * @sa mess_matrix_rank
 *
 * @snippet "tutorials/matrix/tutorial_norms.c" CODE
 *
 */

///@cond
///[CODE]
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/mess.h"

int main (int argc, char *argv[])
{
    mess_matrix mat_coord;
    mess_matrix mat_csr;
    mess_matrix mat_csc;
    double nrm1,nrm2,nrmf, nrminf;
    mess_int_t rank =0;

    MESS_INIT_MATRICES(&mat_coord,&mat_csr,&mat_csc);

    mess_error_level = 1;


    if ( argc != 2) {
        printf("compute the norms of a matrix\n");
        printf("usage: %s matrix.mtx\n", argv[0]);
    }

    mess_matrix_read(argv[1], mat_coord);
    mess_matrix_convert(mat_coord, mat_csr, MESS_CSR);
    mess_matrix_convert(mat_csr, mat_csc, MESS_CSC);


    printf ("CSR: \n");
    mess_matrix_norm1(mat_csr, &nrm1);
    mess_matrix_norminf(mat_csr, &nrminf);
    mess_matrix_rank(mat_csr, &rank);
    mess_matrix_norm2(mat_csr, &nrm2);
    mess_matrix_normf(mat_csr, &nrmf);
    printf (" 1-Norm:  \t %g\n", nrm1);
    printf (" Inf-Norm:\t %g\n", nrminf);
    printf (" F-Norm:  \t %g\n", nrmf);
    printf (" 2-Norm:  \t %g\n", nrm2);
    printf (" Rank:    \t " MESS_PRINTF_INT "\n", rank);

    printf ("CSC: \n");
    mess_matrix_norm1(mat_csc, &nrm1);
    mess_matrix_norminf(mat_csc, &nrminf);
    mess_matrix_norm2(mat_csc, &nrm2);
    mess_matrix_normf(mat_csc, &nrmf);

    printf (" 1-Norm:  \t %g\n", nrm1);
    printf (" Inf-Norm:\t %g\n", nrminf);
    printf (" F-Norm:  \t %g\n", nrmf);
    printf (" 2-Norm:  \t %g\n", nrm2);


    printf ("COORD: \n");
    mess_matrix_norm1(mat_coord, &nrm1);
    mess_matrix_norminf(mat_coord, &nrminf);
    mess_matrix_norm2(mat_coord, &nrm2);
    mess_matrix_normf(mat_coord, &nrmf);

    printf (" 1-Norm:  \t %g\n", nrm1);
    printf (" Inf-Norm:\t %g\n", nrminf);
    printf (" F-Norm:  \t %g\n", nrmf);
    printf (" 2-Norm:  \t %g\n", nrm2);

    printf("\n2-Norm inv(A)\n");
    mess_matrix_norm2inv(mat_csr, &nrm2);
    printf(" = %lg\n", nrm2);

    MESS_CLEAR_MATRICES(&mat_coord,&mat_csr,&mat_csc);

    return 0;
}
///[CODE]
///@endcond
