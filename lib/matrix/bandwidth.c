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
 * @file lib/matrix/bandwidth.c
 * @brief Determine the bandwidth of a matrix.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Determine the bandwidth of a matrix.
 * @param[in] matrix    input matrix
 * @param[out] kl       output bandwidth below the diagonal
 * @param[out] ku       output bandwidth above the diagonal
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_bandwidth function determines the bandwidth of a matrix. The bandwidth
 * is returned as bandwidth below and above the diagonal like it is necessary for the banded
 * matrix subroutines in @lapack. If the input matrix is dense both bandwidths are set to the number of rows -1 or
 * respectively the number of columns -1.
 * 
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 */
int mess_matrix_bandwidth(mess_matrix matrix, mess_int_t *kl, mess_int_t *ku)
{
    MSG_FNAME(__func__);
    mess_int_t i;
    mess_int_t klm, kum;
    mess_int_t pos;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  Check Input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(kl);
    mess_check_nullpointer(ku);

    if (MESS_IS_DENSE(matrix)) {
        *kl = matrix->rows -1;
        *ku = matrix->cols -1;
        return 0;
    }

    *kl = 0;
    *ku = 0;
    klm = kum = 0;

    if (MESS_IS_CSR(matrix)) {
        ret = mess_matrix_sort(matrix);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_sort );
        for (i = 0; i < matrix->rows; i++) {
            /* Left most element   */
            pos = matrix->colptr[matrix->rowptr[i]];
            if ( pos < i && ( (i - pos) > klm )) {
                klm = i - pos;
            }
            /* Right most element  */
            pos = matrix->colptr[matrix->rowptr[i+1]-1];
            if ( pos > i && ( (pos-i) > kum )) {
                kum =  pos - i;
            }
        }
        *kl = klm;
        *ku = kum;
        return 0;
    } else if  (MESS_IS_CSC(matrix)) {
        ret = mess_matrix_sort(matrix);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_sort );
        for (i = 0; i < matrix->cols; i++) {
            /* Left most element   */
            pos = matrix->rowptr[matrix->colptr[i]];
            if ( pos < i && ( (i - pos) > kum )) {
                kum = i - pos;
            }
            /* Right most element  */
            pos = matrix->rowptr[matrix->colptr[i+1]-1];
            if ( pos > i && ( (pos-i) > klm )) {
                klm =  pos - i;
            }
        }
        *kl = klm;
        *ku = kum;
        return 0;
    }

    MSG_ERROR("Matrix storage type %s is not supported.\n", mess_storage_t_str(matrix->store_type));
    return MESS_ERROR_STORAGETYPE;
}

