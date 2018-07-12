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
 * @file tutorials/matrix/tutorial_matrix_operations.c
 * @brief Tutorial on matrix operations.
 * @author @mbehr
 *
 *  # Tutorial on Matrix Operations
 * In this tutorial we show some algebraic matrix operations. We will perform the following steps
 * \f[ B \leftarrow A^{T}, \f]
 * \f[ B \leftarrow 2A + A^{T}, \f]
 * \f[ B \leftarrow (1+i)B, \f]
 * \f[ B \leftarrow \overline{B}, \f]
 * \f[ C \leftarrow A*B^{H}, \f]
 * \f[ D \leftarrow \operatorname{Re}(C), \f]
 * \f[ D \leftarrow \operatorname{Im}(C). \f]
 * \f[ || A ||_F, \; ||A||_\infty, \ldots \f]
 *
 * ## 1. Include Header Files
 * We include standard @c C header files and in addtion the @c mess.h header file
 * to make the @mess library available.
 * @snippet "tutorials/matrix/tutorial_matrix_operations.c" HEADER
 * \n
 *
 * ## 2. (Optional) Print @mess Status Informations
 * You can use @ref mess_version to print status informations about your @mess version.
 * @snippet "tutorials/matrix/tutorial_matrix_operations.c" MESSVERSION
 * \n
 *
 * ## 3. Declare, Init and Read Matrices
 * At first we declare some @ref mess_matrix structures.
 * @snippet "tutorials/matrix/tutorial_matrix_operations.c" DECLARE
 * Now we check the number of input arguments and init the @ref mess_matrix structures.
 * @snippet "tutorials/matrix/tutorial_matrix_operations.c" INITMATRICES
 * Reading from @mm file format works the following way.
 * @snippet "tutorials/matrix/tutorial_matrix_operations.c" READMATRICES
 * Note that the required memory is allocated by @ref mess_matrix_read.
 * @ref mess_matrix_convert converts to @ref MESS_CSR storage format and the required memory is allocated by @ref mess_matrix_convert.
 * \n
 *
 * ## 4. Transpose, Add and Print
 * Transposing, Adding, Scaling and Printing matrices is performed by @ref mess_matrix_ctranspose, @ref mess_matrix_add, @ref mess_matrix_scale / @ref mess_matrix_scalec and @ref mess_matrix_print.
 * Most functions of @mess follow the convention that the @a "input" or @a "source" argument(s) is (are) the first one(s) and the
 * @a "output" or @a "target" argument(s) is (are) the last one(s). The prefix @a mess_matrix  indicates that the function is related to  @ref mess_matrix structures.
 * Now let us do the following operations:
 * \f[ B \leftarrow A^{T}, \f]
 * \f[ B \leftarrow 2A + A^{T}. \f]
 *
 * @snippet "tutorials/matrix/tutorial_matrix_operations.c" ADD
 * \n
 *
 * ## 5. Scale, Complex Conjugate, Multiply and Real/Imaginary Part
 * Next, we perform
 * \f[ B \leftarrow (1+i)B, \f]
 * \f[ B \leftarrow \overline{B}, \f]
 * \f[ C \leftarrow A*B^{H}, \f]
 * \f[ A \leftarrow \operatorname{Re}(C), \f]
 * \f[ A \leftarrow \operatorname{Im}(C). \f]
 * @snippet "tutorials/matrix/tutorial_matrix_operations.c" CONJANDMUL
 * \n
 *
 * ## 6. Compute Martix Norms
 * @mess contains the following functions to compute matrix norms:
 * \li @ref mess_matrix_norm1, to compute the 1-norm,
 * \li @ref mess_matrix_norm2, to compute the 2-norm,
 * \li @ref mess_matrix_norminf, to compute the \f$\infty\f$-norm and
 * \li @ref mess_matrix_normf, to compute the Frobenius norm.
 * Depending on the storage scheme of the matrix the 2-norm computation may get expensive.
 * The example computes all four norms of a previously generated matrix \f$C\f$.
 *
 * @snippet "tutorials/matrix/tutorial_matrix_operations.c" NORMS
 *
* ## 7. Clear Memory
* @snippet "tutorials/matrix/tutorial_matrix_operations.c" CLEAR
*
* \n
*
*
* @sa @ref tutorials
*/

///@cond
///[HEADER]
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/mess.h"
///[HEADER]

int main (int argc, char **argv) {

    ///[MESSVERSION]
    mess_version();
    ///[MESSVERSION]

    ///[DECLARE]
    mess_matrix input, A, B, C, D;
    ///[DECLARE]

    ///[INITMATRICES]
    if ( argc != 2){
        printf("usage: %s <file.mtx>\n", argv[0]);
        return 1;
    }
    mess_matrix_init(&A);
    mess_matrix_init(&B);
    mess_matrix_init(&C);
    mess_matrix_init(&D);
    mess_matrix_init(&input);
    ///[INITMATRICES]

    ///[READMATRICES]
    mess_matrix_read(argv[1], input);
    mess_matrix_convert(input, A, MESS_CSR);
    ///[READMATRICES]

    ///[ADD]
    mess_matrix_ctranspose(A, B);

    mess_matrix_add(2.0, A, 1.0, B);
    ///[ADD]

    ///[CONJANDMUL]
    mess_matrix_scalec(1+1*I,B);

    mess_matrix_conj(B);

    mess_matrix_multiply(MESS_OP_NONE,A,MESS_OP_HERMITIAN, B, C);

    mess_matrix_realpart(C,D);

    mess_matrix_imagpart(C,D);
    ///[CONJANDMUL]

    ///[NORMS]
    double nrm1,nrm2,nrmf, nrminf;
    mess_matrix_norm1(A, &nrm1);
    mess_matrix_norminf(A, &nrminf);
    mess_matrix_norm2(A, &nrm2);
    mess_matrix_normf(A, &nrmf);

    printf (" 1-Norm:  \t %g\n", nrm1);
    printf (" Inf-Norm:\t %g\n", nrminf);
    printf (" F-Norm:  \t %g\n", nrmf);
    printf (" 2-Norm:  \t %g\n", nrm2);
    ///[NORMS]

    ///[CLEAR]
    mess_matrix_clear(&A);
    mess_matrix_clear(&B);
    mess_matrix_clear(&C);
    mess_matrix_clear(&D);
    mess_matrix_clear(&input);
    return 0;
    ///[CLEAR]
}
///@endcond

