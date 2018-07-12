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
 * @file lib/matrix/getelement.c
 * @brief Read access to matrix elements.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


/**
 * @brief Get an element from a matrix.
 * @param[in] matrix    input matrix
 * @param[in] row       input the row index
 * @param[in] col       input the column index
 * @param[out] val      output the real value, @c NULL if not need
 * @param[out] valc     output the complex value if need, @c NULL if not wanted
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_getelement function gets an element from a matrix.\n
 * Please do not use it for traversing the whole matrix because this getting slow especially on sparse matrices. \n
 * If the matrix is real and both @p val  and @p valc  pointer are given, the  (@p row, @p col)  entry of the
 * matrix is stored to both. \n
 * If the matrix is complex and only the @p val  argument is given, the real part of the entry is returned.
 * If both @p val  and @p valc  are given and the matrix is complex the @p valc  argument returns the
 * complex entry and the @p val  argument only the real part.
 *
 * @sa mess_matrix_setelement_complex
 * @sa mess_matrix_setelement
 *
 */
int mess_matrix_getelement ( mess_matrix matrix, mess_int_t row, mess_int_t col, double *val, mess_double_cpx_t *valc )
{
    MSG_FNAME(__func__);
    int ret_val = 0;
    int ret_valc = 0;

    mess_check_nullpointer(matrix);
    mess_check_real_or_complex(matrix);

    if ( val != NULL ) ret_val = 1;
    if ( valc != NULL ) ret_valc = 1;
    if ( row >= matrix->rows || row < 0 ) {
        MSG_ERROR("row index out of range row = " MESS_PRINTF_INT " , matrix->rows=" MESS_PRINTF_INT "\n", row, matrix->rows);
        return(MESS_ERROR_DIMENSION);
    }
    if ( col >= matrix->cols || col < 0 ) {
        MSG_ERROR("column index out of range col = " MESS_PRINTF_INT " , matrix->cols=" MESS_PRINTF_INT "\n", col, matrix->cols);
        return(MESS_ERROR_DIMENSION);
    }

    if ( MESS_IS_DENSE(matrix)) {
        mess_double_cpx_t v;
        if ( MESS_IS_COMPLEX(matrix)) {
            v = matrix->values_cpx[col*matrix->ld+row];
        } else  {
            v = matrix->values[col*matrix->ld+row];
        }
        if ( ret_val ) *val = creal ( v) ;
        if ( ret_valc) *valc = v;
    } else if ( MESS_IS_CSR ( matrix) ){
        mess_int_t i;
        mess_double_cpx_t v = 0 ;
        for ( i = matrix->rowptr[row]; i < matrix->rowptr[row+1]; i++) {
            if ( matrix->colptr[i]==col) {
                v = (matrix->data_type == MESS_COMPLEX) ?  matrix->values_cpx[i] : matrix->values[i];
            }
        }
        if ( ret_val ) *val = creal ( v) ;
        if ( ret_valc) *valc = v;
    } else if (MESS_IS_CSC(matrix)) {
        mess_int_t i;
        mess_double_cpx_t v = 0 ;
        for ( i = matrix->colptr[col]; i < matrix->colptr[col+1]; i++) {
            if ( matrix->rowptr[i]==row) {
                v = (matrix->data_type == MESS_COMPLEX) ?  matrix->values_cpx[i] : matrix->values[i];
            }
        }
        if ( ret_val ) *val = creal ( v) ;
        if ( ret_valc) *valc = v;
    } else if (MESS_IS_COORD(matrix)){
        mess_double_cpx_t v = 0 ;
        mess_int_t i = 0 ;

        for ( i = 0; i < matrix->nnz; i++) {
            if ( matrix->rowptr[i] == row && matrix->colptr[i]==col){
                v = (matrix->data_type == MESS_COMPLEX) ?  matrix->values_cpx[i] : matrix->values[i];
            }
        }
        if ( ret_val ) *val = creal ( v) ;
        if ( ret_valc) *valc = v;

    } else {
        MSG_ERROR("Storage type not supported/not known. storage type=%s\n", mess_storage_t_str(matrix->store_type));
        return(MESS_ERROR_STORAGETYPE);
    }
    return(0);
}       /* -----  end of function mess_matrix_getelement  ----- */

