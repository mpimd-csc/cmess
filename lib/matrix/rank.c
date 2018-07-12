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
 * @file lib/matrix/rank.c
 * @brief Estimate the rank of a matrix.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/eigenvalue.h"
#include "mess/error_macro.h"
#include <complex.h>


/**
 * @brief Estimate the rank of a matrix.
 * @param[in] A         input matrix \f$A\f$
 * @param[out] rank     output estimated rank of \f$A\f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_rank function estimates the rank of a matrix.
 * It uses a similar algorithm like it is done by @matlab.
 \verbatim
 S=SVD(A)
 tol = max(size(A))*eps*max(S)
 rank = sum ( s>tol)
 \endverbatim
 *
 */
int  mess_matrix_rank ( mess_matrix A, mess_int_t *rank )
{
    MSG_FNAME(__func__);
    mess_matrix work;
    mess_vector sigma;
    int conv = 0;
    int ret = 0;
    double eps = mess_eps();
    double tol = 0.0;
    mess_int_t i = 0;
    mess_int_t r = 0;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(rank);
    *rank=0;
    MESS_MATRIX_CHECKFORMAT(A, work, conv, MESS_DENSE);

    ret = mess_vector_init(&sigma);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sigma,MESS_MAX(work->rows, work->cols), MESS_REAL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_eigen_svd(work,sigma,NULL,NULL);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_svd);

    tol = MESS_MAX(work->cols, work->rows) * eps * sigma->values[0];
    // MSG_ERROR("tol: %lg\n",tol);
    for ( i = 0; i<sigma->dim; i++){
        if ( sigma->values[i] > tol) {
            r++;
        } else {
            break;
        }
    }
    *rank = r;

    mess_vector_clear(&sigma);
    if ( conv == 0) {
        mess_matrix_clear(&work);
    }

    return(0);
}       /* -----  end of function mess_matrix_rank  ----- */

