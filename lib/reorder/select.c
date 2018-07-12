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
 * @file lib/reorder/select.c
 * @brief Select reordering by name.
 * @author @koehlerm
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Select and compute a matrix reordering.
 * @param[in] reorder   input @ref mess_reorder_t type of the reordering
 * @param[in] A         input matrix \f$A\f$ to reorder
 * @param[out] p        output row permutation vector of length  rows of \f$A\f$
 * @param[out] q        output column permutation vector of length columns of \f$A\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_reorder function selects and computes a reordering given by the value of the reorder argument.
 * The possible reorderings are listed in the \ref mess_reorder_t enumeration. If the input matrix \f$A\f$ is dense,
 * the identity permutation is returned.
 *
 * Both output vectors \f$ p \f$ and \f$ q \f$ need to be preallocated to the right size before.
 *
 * @sa mess_matrix_reorder_amd
 * @sa mess_matrix_reorder_colamd
 * @sa mess_matrix_reorder_rcm
 *
 */
int mess_matrix_reorder ( mess_reorder_t reorder, mess_matrix A, mess_int_t *p, mess_int_t *q ){
    MSG_FNAME(__func__);
    mess_int_t i;
    int ret = 0;

    mess_check_nullpointer(A);
    mess_check_nullpointer(p);
    mess_check_nullpointer(q);

    if ( reorder == MESS_REORDER_NONE || MESS_IS_DENSE(A)) {
        for (i=0; i < A->rows; i++) p[i]=i;
        for (i=0; i < A->cols; i++) q[i]=i;
        return 0;
    }

#ifdef MESS_HAVE_AMD
    if ( reorder == MESS_REORDER_AMD) {
        mess_check_square(A);
        ret = mess_matrix_reorder_amd(A, p);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_reorder_amd);
        memcpy(q, p, sizeof(mess_int_t)*A->rows);
        return 0;
    }
#endif
#ifdef MESS_HAVE_COLAMD
    if ( reorder == MESS_REORDER_COLAMD) {
        mess_check_square(A);
        ret = mess_matrix_reorder_colamd(A, q);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_reorder_colamd);
        memcpy(p, q, sizeof(mess_int_t)*A->rows);
        return 0;
    }
#endif

    if ( reorder == MESS_REORDER_RCM) {
        mess_check_square(A);
        ret = mess_matrix_reorder_rcm(A, q);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_reorder_rcm);
        memcpy(p, q, sizeof(mess_int_t)*A->rows);
        return 0;
    }
    MSG_WARN("Unknown or unsupported reordering. Use identity instead.\n");
    for (i=0; i < A->rows; i++) p[i]=i;
    for (i=0; i < A->cols; i++) q[i]=i;

    return 0;
}       /* -----  end of function mess_matrix_reorder  ----- */


