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
 * @file lib/matrix/eye.c
 * @brief Generate identity matrices.
 * @author @koehlerm
 */


#include <stdlib.h>
#include <stdio.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#define MIN(A,B)    ((A)<(B))?(A):(B)

/**
 * @brief Generate a real identity matrix.
 * @param[out] matrix   output matrix
 * @param[in] rows      input number of rows
 * @param[in] cols      input number of columns
 * @param[in] store     input storage type of the matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_eye function generates a real identity matrix of size \f$ (rows \times cols) \f$.
 * @see mess_matrix_eyec
 */
int mess_matrix_eye(mess_matrix matrix, mess_int_t rows, mess_int_t cols, mess_storage_t store){
    MSG_FNAME(__func__);
    mess_int_t mn = 0;
    int ret = 0;
    mess_int_t i = 0;

    mess_check_nullpointer(matrix);
    MESS_MATRIX_RESET(matrix);

    if (  rows <= 0 || cols <=0 ){
        MSG_ERROR("rows or cols has an invalid value: rows = " MESS_PRINTF_INT "  cols = " MESS_PRINTF_INT "\n", rows, cols);
        return (MESS_ERROR_DIMENSION);
    }
    mn = MIN(rows, cols);

    ret = mess_matrix_alloc(matrix, rows, cols, mn, store, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);


    switch ( store ) {
        case MESS_DENSE:
            for (i = 0;i<mn;i++){
                matrix->values[i*matrix->ld+i] = 1;
            }
            break;
        case MESS_CSR:
        case MESS_CSC:
        case MESS_COORD:
            for (i = 0;i<mn;i++){
                matrix->values[i] = 1;
                matrix->colptr[i] = i;
                matrix->rowptr[i] = i;
            }
            break;
        default:
            MSG_ERROR("unknown storage type: %s \n", mess_storage_t_str(store));
            return(MESS_ERROR_STORAGETYPE);
            break;
    }               /* -----  end switch  ----- */

    return(0);
}

/**
 * @brief Generate a complex identity matrix.
 * @param[out] matrix       output matrix
 * @param[in] rows      input number of rows
 * @param[in] cols      input number of columns
 * @param[in] store     input storage type of the matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_eyec function generates a complex identity matrix of size \f$ (rows \times cols ) \f$.
 * @see mess_matrix_eye
 */
int mess_matrix_eyec(mess_matrix matrix, mess_int_t rows, mess_int_t cols, mess_storage_t store){
    MSG_FNAME(__func__);
    mess_int_t mn = 0;
    int ret = 0;
    mess_int_t i = 0;

    mess_check_nullpointer(matrix);
    MESS_MATRIX_RESET(matrix);
    if (  rows <= 0 || cols <=0 ){
        MSG_ERROR("rows or cols has an invalid value: rows = " MESS_PRINTF_INT "  cols = " MESS_PRINTF_INT "\n", rows, cols);
        return(MESS_ERROR_DIMENSION);
    }
    mn = MIN(rows, cols);

    ret = mess_matrix_alloc(matrix, rows, cols, mn, store, MESS_COMPLEX);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);


    switch ( store ) {
        case MESS_DENSE:
            for (i = 0;i<mn;i++){
                matrix->values_cpx[i*matrix->ld+i] = 1;
            }
            break;
        case MESS_CSR:
        case MESS_CSC:
        case MESS_COORD:
            for (i = 0;i<mn;i++){
                matrix->values_cpx[i] = 1;
                matrix->colptr[i] = i;
                matrix->rowptr[i] = i;
            }
            break;
        default:
            MSG_ERROR("unknown storage type: %s\n", mess_storage_t_str(store));
            return (MESS_ERROR_STORAGETYPE);
            break;
    }               /* -----  end switch  ----- */

    return(0);
}

