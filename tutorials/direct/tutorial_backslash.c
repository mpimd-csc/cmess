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
 * @file tutorials/direct/tutorial_backslash.c
 * @brief Demostrate how to solve a linear system using a \f$ \backslash \f$-operator.
 *
 * @sa mess_matrix_backslash
 *
 * # Tutorial: Use mess_matrix_backslash to solve a linear system.
 *
 * This function demonstrates how to solve a linear system
 * \f[ A   x = b \f]
 * or
 * \f[ A^T x = b \f]
 * using a backslash operator.
 *
 * @snippet "tutorials/direct/tutorial_backslash.c" CODE
 *
 */

///@cond
///[CODE]
#include <stdio.h>
#include <stdlib.h>
#include "mess/mess.h"

int main ( int argc, char **argv) {
    mess_matrix A;
    mess_vector x,x2,b,b2;

    mess_version();

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 2 ) {
        printf("Usage: %s file\n", argv[0]);
        return 1;
    }
    mess_set_errorlevel(3);

    /*-----------------------------------------------------------------------------
     *  init and read matrix
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A);
    mess_matrix_read_formated(argv[1],A,MESS_DENSE);

    /*-----------------------------------------------------------------------------
     *  init vectors
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&x,&x2,&b,&b2);
    mess_vector_alloc(x,  A->cols, A->data_type);
    mess_vector_alloc(x2, A->rows, A->data_type);
    mess_vector_alloc(b,  A->rows, A->data_type);
    mess_vector_alloc(b2, A->cols, A->data_type);

    /*-----------------------------------------------------------------------------
     *  create right hand side
     *-----------------------------------------------------------------------------*/
    mess_matrix_rowsums(A,b);
    mess_matrix_colsums(A,b2);

    /*-----------------------------------------------------------------------------
     *  solve systems
     *-----------------------------------------------------------------------------*/
    mess_matrix_backslash(MESS_OP_NONE, A, b, x);
    mess_matrix_backslash(MESS_OP_TRANSPOSE,A, b2, x2);

    /*-----------------------------------------------------------------------------
     *  print solution should be [1,..,1]
     *-----------------------------------------------------------------------------*/
    printf("Ax=b:\n");
    mess_vector_printshort(x);
    printf("A^Tx=b:\n");
    mess_vector_printshort(x2);

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A);
    MESS_CLEAR_VECTORS(&x,&x2,&b,&b2);

    return 0;
}
///[CODE]
///@endcond

