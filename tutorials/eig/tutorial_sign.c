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
 * @file tutorials/eig/tutorial_sign.c
 * @brief Tutorial for demonstrating how to compute the sign of a matrix using @ref mess_eigen_sign.
 * @author @mbehr
 * @author @koehlerm
 *
 * # Tutorial for computing the sign of a matrix.
 * In this tutorial we show how to use the @ref mess_eigen_sign function to
 * compute the sign of a function. We compute the eigenvalues of
 * \f$ Z=sign(A)\f$. We also check if \f$ Z^2=I_n\f$.

 * ## 1. Include Header Files
 * We include standard @c C header files and in addtion the @c mess.h header file
 * to make the @mess library available.
 * @snippet "tutorials/eig/tutorial_sign.c" HEADER
 * \n
 *
 *
 * ## 2. Check Input arguments, Init Matrix and Read Data
 * We declare four @ref mess_matrix instances, one @ref mess_vector instance and
 * three double value variables. \n
 * We check the number of input arguments by using the @c argc variable and use @ref mess_matrix_read_formated
 * to read from a @mm file directly into a @ref mess_matrix with storage type @ref MESS_DENSE.
 * @snippet "tutorials/eig/tutorial_sign.c" DECLAREINIT
 * \n
 *
 * ## 3. Compute the sign
 * The method @ref mess_eigen_sign takes the matrix \f$ A \f$ as argument, stores the result \f$Z=sign(A)\f$ and
 * the last argument is the @ref mess_sign_scale_t enum type to speedup the the computation.  \n
 * Using @ref mess_wtime we can measure the elapsed wall clock time.
 * @snippet "tutorials/eig/tutorial_sign.c" SIGN
 * \n
 *
 *
 * ## 4. CHECK
 * The eigenvalues of \f$ Z\f$ should be either \f$ 1 \f$ or \f$ -1\f$ and \f$ Z \f$ should fullfill
 * \f$ Z^2-I_n=0_n\f$. \n
 * We want to check these properties. Therefore we allocate a @ref mess_vector instance for @ref mess_eigen_eig to
 * compute the eigenvalues. Since we are not interested in the eigenvectors the last argument of @ref mess_eigen_eig
 * is @c NULL. \n
 * In order to check if \f$ Z^2-I_n = 0_n\f$, we use @ref mess_matrix_eye to create a unit matrix @c eye.
 * The second step is to compute \f$ Z^2\f$ and store the result in @c Z2. We check the identity
 * \f$ Z^2-I_n\f$ by computing \f$||Z^2-I_n||_2\f$. The spectral norm of the difference of two matrices can be computed
 *  with @ref mess_matrix_diffnorm and the frobenius norm of the difference with @ref mess_matrix_diffnormf.
 * @snippet "tutorials/eig/tutorial_sign.c" CHECK
 * \n
 *
 * ## 5. Clear the memory
 * Allocated memory has to be cleared. @ref MESS_CLEAR_MATRICES and @ref MESS_CLEAR_VECTORS do this job.
 *
 * @snippet "tutorials/eig/tutorial_sign.c" CLEAR
 *
 *
 * @sa @ref tutorials
 *
 */

//cond -> code is ignored by doxygen documentation, but snippet references are working

///@cond
///[HEADER]
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/mess.h"
///[HEADER]

int main ( int argc,char ** argv ) {

    /// [DECLAREINIT]
    mess_matrix A, Z, Z2, eye;
    mess_vector ev;
    double ts, te, diff;

    if ( argc != 2 ) {
        printf("usage: %s <FILE>\n", argv[0]);
        exit(1);
    }

    MESS_INIT_MATRICES(&A,&Z,&Z2,&eye);
    mess_matrix_read_formated(argv[1],A,MESS_DENSE);
    /// [DECLAREINIT]

    /// [SIGN]
    ts = mess_wtime();
    mess_eigen_sign(A,Z,MESS_SIGN_SCALE_NONE);
    te = mess_wtime();
    printf("Sign function tooks %lg seconds\n", te-ts);
    /// [SIGN]

    /// [CHECK]
    mess_vector_init(&ev);
    mess_vector_alloc(ev,A->rows,MESS_COMPLEX);
    mess_eigen_eig(Z,ev,NULL);
    printf("Eigenvalues of A\n");
    mess_vector_printshort(ev);

    mess_matrix_eye(eye,Z->rows,Z->cols,MESS_DENSE);
    mess_matrix_multiply(MESS_OP_NONE,Z,MESS_OP_NONE,Z,Z2);
    mess_matrix_diffnorm(Z2,eye,&diff);
    printf("||sign(A)^2-I||_2=%e\n",diff);
    /// [CHECK]

    /// [CLEAR]
    MESS_CLEAR_MATRICES(&A,&Z,&Z2,&eye);
    MESS_CLEAR_VECTORS(&ev);
    /// [CLEAR]

    return 0;
}
///@endcond




