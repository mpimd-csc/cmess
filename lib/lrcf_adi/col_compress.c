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
 * @file lib/lrcf_adi/col_compress.c
 * @brief Column compression for the ADI process.
 * @author @koehlerm
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Perform a column compression on a low rank factor.
 * @param[in,out] Z  current low rank factor
 * @param[in] ccTol input      cut-off tolerance in the rank decission
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref mess_lrcfadi_colcompress function performs a column compression on a low rank factor
 * iterate. The detailed description can be found in \cite Saa09 . This function
 * is only a wrapper around \ref mess_lrcfadi_ccsvd in order to change its implementation
 * in future times.
 *
 * \sa mess_lrcfadi_ccsvd
 *
 */
int mess_lrcfadi_colcompress(mess_matrix Z, double ccTol) {
    return mess_lrcfadi_ccsvd(Z, ccTol);
}



/**
 * @brief Perform a SVD based column compression on low rank factor.
 * @param[in,out] Z     input/output current ADI factor
 * @param[in] ccTol     input      cut-off tolerance for the singular values
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_ccsvd function performs a column compression on a low rank factor.
 * It is the current backend behind \ref mess_lrcfadi_colcompress.
 *
 * It works using a SVD on the input matrix \f$Z\f$ and decide uppon the rank of the
 * matrix using the cut-off tolerance @c ccTol.
 * All singular values that are smaller than
 * \f[ MAX(rows(Z),cols(Z)) \cdot ccTol \cdot \sigma_1\f]
 * are neglected, where \f$ \sigma_1 \f$ is the largest singular value.
 *
 * \sa mess_lrcfadi_colcompress
 *
 */
int mess_lrcfadi_ccsvd ( mess_matrix Z, double ccTol )
{
    MSG_FNAME(__func__);
    mess_matrix U, S;
    mess_vector sigma;
    int ret =0;
    mess_int_t r;
    mess_int_t M,N,i;
    double tol;

    mess_check_nullpointer(Z);
    if (ccTol < 0) ccTol = mess_eps();

    M= Z->rows;
    N= Z->cols;
    ret = mess_matrix_init(&U);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&S);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_vector_init(&sigma);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);

    ret = mess_vector_alloc(sigma, MESS_MIN(M,N), MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_init);

    ret = mess_eigen_svd_econ(Z,sigma,U,NULL);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_svd_econ);
    tol = MESS_MAX(M,N) * ccTol * sigma->values[0];
    r = 0;
    for (i=0; i < sigma->dim; i++){ if ( sigma->values[i] > tol ) r++; else break; }
    MSG_INFO("est. rank: " MESS_PRINTF_INT "\n", r);

    ret = mess_matrix_alloc(S,r,r,r*r, MESS_DENSE, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    for (i=0; i < r; i++){
        S->values[i+i*r] = sigma->values[i];
    }
    ret= mess_matrix_resize(U,U->rows, r);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_resize);
    ret= mess_matrix_multiply(MESS_OP_NONE, U, MESS_OP_NONE, S, Z); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    mess_matrix_clear(&U);
    mess_matrix_clear(&S);
    mess_vector_clear(&sigma);

    return 0;
}       /* -----  end of function mess_lrcfadi_ccsvd  ----- */

