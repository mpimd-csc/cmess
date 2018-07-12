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
 * @file tutorials/h2/tutorial_h2norm.c
 * @brief Demonstrate how to compute the \f$ \mathcal{H}_2 \f$-norm of a LTI system.
 *
 * # Tutorial: Compute \f$ \mathcal{H}_2 \f$ error.
 *
 * This function demonstrates computing the \f$ \mathcal{H}_2 \f$-norm of a linear time invariant system
 * \f[
 * \begin{array}{ccc}
 * \dot{x} &= A x & + B u \\
 *       y &= C x &
 *  \end{array} .
 * \f]
 *
 * @sa mess_h2_norm
 *
 *
 * @snippet "tutorials/h2/tutorial_h2norm.c" CODE
 */


///@cond
///[CODE]
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


int main ( int argc, char **argv){
    mess_matrix A,B,C;
    mess_dynsys lti;
    double norm = 0;
    int ret = 0;
    double ts, te;
    mess_version();
    mess_set_errorlevel(2);

    printf("mess - h2norm\n");
    printf("==============\n");


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc!=4) {
        printf("usage: %s A.mtx B.mtx C.mtx \n", argv[0]);
        return 0;
    }

    /*-----------------------------------------------------------------------------
     *  init matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&B,&C);

    mess_matrix_read_formated(argv[1], A, MESS_CSR);
    mess_matrix_read_formated(argv[2], B, MESS_DENSE);
    mess_matrix_read_formated(argv[3], C, MESS_DENSE);


    printf("Dimension: " MESS_PRINTF_INT "\n", A->rows);
    printf("Inputs:    " MESS_PRINTF_INT "\n", B->cols);
    printf("Outputs:   " MESS_PRINTF_INT "\n", C->rows);
    mess_dynsys_init(&lti);
    mess_dynsys_lti_copy(lti, A, B, C);

    ts = mess_wtime();
    ret = mess_h2_norm(lti, &norm);
    te = mess_wtime();
    if ( ret !=0) {
        printf("error computing h2 norm: %s\n", mess_get_error(ret));
        return 1;
    } else {
        printf("h2 norm: %lg\n", norm);
    }
    printf("time: %lg\n", te-ts);


    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&B,&C);
    mess_dynsys_clear(&lti);

    return 0;
}
///[CODE]
///@endcond

