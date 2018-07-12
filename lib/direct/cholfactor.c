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
 * @file lib/direct/cholfactor.c
 * @brief Compute a Cholesky factor \f$ X=ZZ^T \f$ of a matrix.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Compute a Cholesky factor of a given matrix.
 * @param[in] ZZ    input symmetric matrix
 * @param[out] Z    output Cholesky factor such that \f$ Z_{out} Z_{out}^T  = ZZ \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_direct_cholfactor function computes a Cholesky factor of the
 * matrix \f$ Z \f$ using the SVD. \n
 * It works for singular matrices, too.
 */
int mess_direct_cholfactor(mess_matrix ZZ, mess_matrix Z)
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_vector S;
    mess_int_t i;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(ZZ);
    mess_check_nullpointer(Z);
    mess_check_real_or_complex(ZZ);

    ret = mess_vector_init(&S);                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_init);
    ret = mess_vector_alloc(S, ZZ->rows, MESS_REAL);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_init);
    ret = mess_eigen_svd(ZZ, S, Z, NULL);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_svd);
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for (i = 0; i < Z->cols; i++) {
        mess_matrix_colscale(Z, i, sqrt(S->values[i]));
    }
    mess_vector_clear(&S);
    return 0;
}
