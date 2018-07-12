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
 * @file lib/matrix/ones.c
 * @brief Generate a matrix containing one number over and over again.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Set all entries of a matrix to one value.
 * @param[in,out] matrix    input/output matrix to be overwritten
 * @param[in]     value     input value to be written to all matrix entries
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_one_value function sets all entries in a given matrix @p matrix  to a given
 * value @p value. \n
 * The matrix needs to be preallocated. If you want to set a complex value
 * use \ref mess_matrix_one_valuec instead.
 *
 * @sa mess_matrix_one_valuec
 *
 */
int mess_matrix_one_value ( mess_matrix matrix, double value )
{
    MSG_FNAME(__func__);
    mess_int_t i,j;

    mess_check_nullpointer(matrix);
    mess_check_real_or_complex(matrix);

    if ( MESS_IS_DENSE(matrix)) {
        if ( MESS_IS_REAL(matrix)){
            for ( j = 0 ; j < matrix->cols; j++ ) {
                for ( i = 0 ; i < matrix->rows; i++ ) {
                    matrix->values[i+j*matrix->ld] = value;
                }

            }
        } else {
            for ( j = 0 ; j < matrix->cols; j++ ) {
                for ( i = 0 ; i < matrix->rows; i++ ) {
                    matrix->values_cpx[i+j*matrix->ld] = value;
                }
            }
        }
    } else if ( MESS_IS_SPARSE(matrix)){
        if ( MESS_IS_REAL(matrix)) {
            for ( i = 0 ; i < matrix->nnz; i++) {
                matrix->values[i] = value;
            }
        } else {
            for ( i = 0 ; i < matrix->nnz; i++) {
                matrix->values_cpx[i] = value;
            }
        }
    } else {
        MSG_ERROR("Unknown storage type (%d - %s)\n", (int) matrix->store_type, mess_storage_t_str(matrix->store_type));
        return MESS_ERROR_STORAGETYPE;
    }

    return 0;
}       /* -----  end of function mess_matrix_one_value  ----- */


/**
 * @brief Set all entries of a matrix to one value (complex version).
 * @param[in,out] matrix    matrix to be overwritten
 * @param[in]     value  input value to be written to all matrix entries
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_one_valuec function sets all entries in a given matrix \f$ matrix \f$ to the given
 * value \f$ value \f$. \n
 * The matrix \f$ matrix \f$ needs to be preallocated. The value is complex in this function.
 * If the matrix is real only the real part of the value is set. \n
 * If you want to avoid this, convert the matrix to a complex one before. \n
 * For real matrices and a real value use \ref mess_matrix_one_value instead.
 *
 * @sa mess_matrix_one_value
 *
 */
int mess_matrix_one_valuec ( mess_matrix matrix, mess_double_cpx_t value )
{
    MSG_FNAME(__func__);
    mess_int_t i,j;

    mess_check_nullpointer(matrix);
    mess_check_real_or_complex(matrix);
    if ( MESS_IS_REAL(matrix) ) {
        return mess_matrix_one_value(matrix, creal(value));
    }

    if ( MESS_IS_DENSE(matrix)) {
        for ( j = 0 ; j < matrix->cols; j++ ) {
            for ( i = 0 ; i < matrix->rows; i++ ) {
                matrix->values_cpx[i+j*matrix->ld] = value;
            }
        }
    } else if ( MESS_IS_SPARSE(matrix)){
        for ( i = 0 ; i < matrix->nnz; i++) {
            matrix->values_cpx[i] = value;
        }
    } else {
        MSG_ERROR("Unknown storage type (%d - %s)\n", (int) matrix->store_type, mess_storage_t_str(matrix->store_type));
        return MESS_ERROR_STORAGETYPE;
    }

    return 0;
}       /* -----  end of function mess_matrix_one_value  ----- */



/**
 * @brief Set all entries in a matrix to one.
 * @param[in,out] matrix    matrix to be set to one
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_ones function sets all entries in a matrix to one. \n
 * The matrix needs to be preallocated. \n
 * It relies on \ref mess_matrix_one_value.
 *
 * @sa mess_matrix_one_value
 */
int mess_matrix_ones ( mess_matrix matrix )
{
    MSG_FNAME(__func__);
    mess_check_nullpointer(matrix);
    mess_check_real_or_complex(matrix);
    return mess_matrix_one_value(matrix, 1.0);
}       /* -----  end of function mess_matrix_ones  ----- */

/**
 * @brief Set all entries in a matrix to zero.
 * @param[in,out] matrix    matrix to be set to zero
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_zeros function sets all entries in a matrix to zero.\n
 * The matrix needs to be preallocated. \n
 * It relies on \ref mess_matrix_one_value.
 *
 * @sa mess_matrix_one_value
 *
 */
int mess_matrix_zeros ( mess_matrix matrix )
{
    MSG_FNAME(__func__);
    mess_int_t ret = 0;

    /*-----------------------------------------------------------------------------
     *   check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_real_or_complex(matrix);

    /*-----------------------------------------------------------------------------
     *  perfrom zero a matrix handle sparse case different from dense case
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_DENSE(matrix)){
        return mess_matrix_one_value(matrix, 0.0);
    }else{
        mess_storage_t st = matrix->store_type;
        mess_datatype_t dt = matrix->data_type;
        mess_int_t rows = matrix->rows;
        mess_int_t cols = matrix->cols;
        MESS_MATRIX_RESET(matrix);
        ret = mess_matrix_alloc(matrix, rows, cols, 0, st, dt);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
        return ret;
    }
}       /* -----  end of function mess_matrix_ones  ----- */

