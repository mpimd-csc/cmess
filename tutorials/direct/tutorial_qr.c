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
 * @file tutorials/direct/tutorial_qr.c
 * @brief Demonstrate how to use a QR solver.
 *
 * @sa mess_direct_create_lapack_qr
 * @sa mess_direct_solve
 * @sa mess_direct_solvem
 *
 * # Tutorial: QR solver
 *
 * This function demonstrates how to generate a QR solver to solve equation
 * \f[ Ax = b \f]
 * or
 * \f[ AX = B, \f]
 * where \f$ A \f$ and \f$ B \f$ are matrices and \f$ b \f$ is a vector.
 *
 *  @snippet "tutorials/direct/tutorial_qr.c" CODE
 */


///@cond
///[CODE]

#include <stdio.h>
#include <stdlib.h>
#include "mess/mess.h"

int main ( int argc, char **argv) {
    mess_matrix A,RHS,X;
    mess_vector b,x;
    mess_direct qr;

    mess_set_errorlevel(3);
    mess_version();

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 2 ) {
        printf("Usage: %s file\n", argv[0]);
        return -1;
    }

    /*-----------------------------------------------------------------------------
     *  init matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&RHS,&X);

    mess_matrix_read_formated(argv[1],A,MESS_DENSE);

    /*-----------------------------------------------------------------------------
     *  init vectors
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&x,&b);
    mess_vector_alloc(x, A->rows, A->data_type);
    mess_vector_alloc(b, A->rows, A->data_type);

    mess_matrix_rowsums(A,b);

    /*-----------------------------------------------------------------------------
     *  create qr solver
     *-----------------------------------------------------------------------------*/
    mess_direct_init(&qr);
    mess_direct_create_lapack_qr(A,qr);

    /*-----------------------------------------------------------------------------
     *  solve system and plot solution
     *-----------------------------------------------------------------------------*/
    mess_direct_solve(MESS_OP_NONE,qr, b, x);
    printf("Ax=b:\n");
    mess_vector_printshort(x);

    /*-----------------------------------------------------------------------------
     *  create rhs and plot solution
     *-----------------------------------------------------------------------------*/
    mess_matrix_alloc(RHS, A->rows, 2, A->rows*2, MESS_DENSE, A->data_type);
    mess_matrix_setcol(RHS, 0, b);
    mess_vector_scalee(2.0, b);
    mess_matrix_setcol(RHS, 1, b);
    mess_direct_solvem(MESS_OP_NONE,qr,RHS, X);

    printf("AX=B:\n");
    mess_matrix_print(X);

    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&X,&RHS);
    MESS_CLEAR_VECTORS(&x,&b);
    mess_direct_clear(&qr);

    return 0;
}
///[CODE]
///@endcond

