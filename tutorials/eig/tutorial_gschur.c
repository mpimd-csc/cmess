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
 * @file tutorials/eig/tutorial_gschur.c
 * @brief Demonstrate how to compute the generalized Schur decomposition of a matrix pair \f$ (A,B) \f$.
 * @author @mbehr
 * @author @koehlerm
 *
 * @sa mess_eigen_gschur
 * @sa mess_eigen_gschur_complex
 *
 *  # Tutorial: Computing generalized Schur Decomposition
 *
 * This function demonstrates computing the generalized Schur decomposition of two matrices \f$ A \f$ and \f$ B \f$
 * \f[
 * \begin{array}{ccc}
 * S &=& U^T A V, \\
 * T &=& U^T B V, \\
 * \end{array}
 * \f]
 * where \f$ S \f$ is an upper quasi-triangular matrix and \f$ T \f$ is an upper triangular matrix
 * and the eigenvalues of the generalized eigenvalue problem:
 * \f[ A x = \lambda B x. \f]
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
    mess_matrix A,B, T,S;
    mess_vector EVA, EVB;
    int type = 1;

    mess_version();

    printf("mess generalized schur decomposition demo\n");
    printf("==============================\n");

    /*-----------------------------------------------------------------------------
     *  check input args
     *-----------------------------------------------------------------------------*/

    if ( argc!=4){
        printf("Usage: %s A.mtx B.mtx type\n", argv[0]);
        printf(" type = 1 real schur form\n");
        printf(" type = 2 complex schur form\n");
        return -1;
    }
    type = atoi(argv[3]);

    /*-----------------------------------------------------------------------------
     *  init and read matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&B,&T,&S);
    mess_vector_init(&EVA);
    mess_vector_init(&EVB);
    mess_vector_alloc(EVA, 1, MESS_COMPLEX);
    mess_vector_alloc(EVB, 1, MESS_REAL);

    mess_matrix_read_formated(argv[1],A,MESS_DENSE);
    mess_matrix_read_formated(argv[2],B,MESS_DENSE);

    /*-----------------------------------------------------------------------------
     *  compute schur decomposition
     *-----------------------------------------------------------------------------*/
    if ( type == 1) {
        mess_eigen_gschur(A,B,S,T,NULL,NULL, EVA, EVB,NULL);
    } else if ( type ==2){
        mess_eigen_gschur_complex(A,B,S,T,NULL, NULL,EVA,EVB,NULL);
    }else{
        printf("FAILED: wrong input argument.");
        return 1 ;
    }

    printf("\nEigenvalues A:\n");
    mess_vector_print(EVA);
    printf("\nEigenvalues B:\n");
    mess_vector_print(EVB);

    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/

    MESS_CLEAR_MATRICES(&A,&B,&S,&T);
    MESS_CLEAR_VECTORS(&EVA,&EVB);

    return 0;
}
///[COND]
///@endcond

