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
 * @file lib/matrix/diag.c
 * @brief Operations on the diagonal elements of a matrix.
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
 * @brief Get the main diagonal elements of from matrix.
 * @param[in] matrix    input matrix
 * @param[out] d        output vector with the main diagonal elements
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_diag function extracts the main diagonal elements of a matrix and copy them to \f$ d \f$.
 * If \f$ d \f$ has not the right size it is resized.
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 */
int mess_matrix_diag( mess_matrix matrix, mess_vector d){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i, j;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(d);
    mess_check_real_or_complex(matrix);

    if ( MESS_IS_REAL(matrix)) {
        ret = mess_vector_toreal_nowarn(d);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    } else {
        ret = mess_vector_tocomplex(d);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
    }
    ret = mess_vector_resize(d,MESS_MIN(matrix->rows,matrix->cols));    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);

    if (MESS_IS_DENSE(matrix)){
        if (MESS_IS_REAL(matrix)){
            for ( i =0 ; i < MESS_MIN(matrix->rows,matrix->cols); i++ ) {
                d->values[i] = matrix->values[i+i*matrix->ld];
            }
        } else {
            for ( i =0 ; i < MESS_MIN(matrix->rows,matrix->cols); i++ ) {
                d->values_cpx[i] = matrix->values_cpx[i+i*matrix->ld];
            }
        }
    } else if (MESS_IS_CSR(matrix)){
        if ( MESS_IS_REAL(matrix)){
            for ( i = 0; i < matrix->rows; i++){
                d->values[i] = 0.0;
                for ( j = matrix->rowptr[i]; j<matrix->rowptr[i+1]; j++){
                    if ( matrix->colptr[j] == i){
                        d->values[i] = matrix->values[j];
                    }
                }
            }
        } else {
            for ( i = 0; i < matrix->rows; i++){
                d->values_cpx[i] = 0.0;
                for ( j = matrix->rowptr[i]; j<matrix->rowptr[i+1]; j++){
                    if ( matrix->colptr[j] == i){
                        d->values_cpx[i] = matrix->values_cpx[j];
                    }
                }
            }
        }
    } else if (MESS_IS_CSC(matrix)){
        if ( MESS_IS_REAL(matrix)){
            for ( i = 0; i < matrix->cols; i++){
                d->values[i] = 0.0;
                for ( j = matrix->colptr[i]; j<matrix->colptr[i+1]; j++){
                    if ( matrix->rowptr[j] == i){
                        d->values[i] = matrix->values[j];
                    }
                }
            }
        } else {
            for ( i = 0; i < matrix->cols; i++){
                d->values_cpx[i] = 0.0;
                for ( j = matrix->colptr[i]; j<matrix->colptr[i+1]; j++){
                    if ( matrix->rowptr[j] == i){
                        d->values_cpx[i] = matrix->values_cpx[j];
                    }
                }
            }
        }
    } else {
        MSG_ERROR("Unknown/Unsupported Storage Type: %s\n", mess_storage_t_str(matrix->store_type));
        return MESS_ERROR_STORAGETYPE;
    }

    return (0) ;
}

/**
 * @brief Find the position of the main diagonal elements.
 * @param[in] matrix input matrix
 * @param[out] pos   output array containing the postion of the diagonal elements, at least @c MIN(rows,cols) elements long
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_diagpos function determines the position of the main diagonal elements in the
 * values array. \n
 * If an element does not exist \f$ pos[e] \f$ is set to \f$ -1 \f$.
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 */
int mess_matrix_diagpos(mess_matrix matrix, mess_int_t *pos) {
    MSG_FNAME(__func__);
    mess_int_t i, j;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(pos);

    if (MESS_IS_DENSE(matrix)){
        for ( i = 0; i < MESS_MIN(matrix->rows,matrix->cols); i++) {
            pos[i] = i * matrix->ld + i;
        }
    } else if ( MESS_IS_CSR(matrix)){
        for ( i = 0; i < matrix->rows; i++){
            pos[i] = -1;
            for ( j = matrix->rowptr[i]; j < matrix->rowptr[i+1]; j++){
                if ( matrix->colptr[j] == i ) pos[i] = j;
            }
        }
    } else if ( MESS_IS_CSC(matrix)){
        for ( i = 0; i < matrix->cols; i++){
            for ( j = matrix->colptr[i]; j < matrix->colptr[i+1]; j++){
                if ( matrix->rowptr[j] == i ) pos[i] = j;
            }
        }
    } else {
        MSG_ERROR("Unknown/Unsupported Storage Type: %s\n", mess_storage_t_str(matrix->store_type));
        return MESS_ERROR_STORAGETYPE;
    }
    return(0);
}


/**
 * @brief Create a diagonal @ref mess_matrix from @ref mess_vector.
 * @param[in]   v       input vector with the main diagonal elements
 * @param[in]   diag    output matrix diagonal matrix with entries from v
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_diag_from_vector function creates a diagonal matrix with enties from @p v.
 * A @ref MESS_DENSE matrix is created and the datatype is determined by the datatype of @p v.
 *
 */
int mess_matrix_diag_from_vector(mess_vector v, mess_matrix diag){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i;

    /*-----------------------------------------------------------------------------
     *  check input args
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(v);
    mess_check_nullpointer(diag);
    mess_check_real_or_complex(diag);

    /*-----------------------------------------------------------------------------
     *  reset and set entries
     *-----------------------------------------------------------------------------*/
    MESS_MATRIX_RESET(diag);
    ret = mess_matrix_alloc(diag, v->dim, v->dim, (v->dim)*(v->dim), MESS_DENSE, v->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    if(MESS_IS_REAL(diag)){
        for(i=0;i<v->dim;++i){
            diag->values[i+diag->ld*i] = v->values[i];
        }
    }else{
        for(i=0;i<v->dim;++i){
            diag->values_cpx[i+diag->ld*i] = v->values_cpx[i];
        }
    }

    return (0) ;
}


