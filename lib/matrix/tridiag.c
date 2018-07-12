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
 * @file lib/matrix/tridiag.c
 * @brief Generate tridiagonal matrices.
 * @author @mbehr
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
 * @brief Create a tridiagonal matrix.
 * @param[in,out]   matrix      input/output matrix to be overwritten
 * @param[in]       rows        input number of rows
 * @param[in]       cols        input number of columns
 * @param[in]       store_t     input desired storage type
 * @param[in]       data_t      input desired data type
 * @param[in]       lower       input value for lower diagonal
 * @param[in]       diag        input value for main diagonal
 * @param[in]       upper       input value for upper diagonal
 *
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_tridiag function creates a tridigonal @p matrix using the values @p lower, @p diag, @p upper. \n
 * If a real matrix is wanted then the imaginary part of @p lower, @p diag, @p upper is ignored.
 *
 */

int mess_matrix_tridiag(mess_matrix matrix, mess_int_t rows, mess_int_t cols, mess_storage_t store_t, mess_datatype_t data_t, mess_double_cpx_t lower, mess_double_cpx_t diag, mess_double_cpx_t upper){
    MSG_FNAME(__func__);
    mess_int_t i=0, j=0, nnz=0, ret = 0, has_lower = 0, has_upper = 0, has_diag = 0;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_storage_type(store_t);
    mess_check_datatype(data_t);
    mess_check_positive(rows);
    mess_check_positive(cols);

    if( (cimag(lower) || cimag(diag) || cimag(upper)) && (data_t == MESS_REAL)){
        MSG_INFO("lower, diag or upper have a non zero imaginary part, but data_t is MESS_REAL. The imaginary part will be ignored and a real tridiagonal matrix is generated\n");
    }

    /*-----------------------------------------------------------------------------
     *  check which values are nonzero and compute number of nonzero values
     *-----------------------------------------------------------------------------*/
    if(data_t == MESS_REAL){
        has_lower   = creal(lower)  ? 1:0;
        has_diag    = creal(diag)   ? 1:0;
        has_upper   = creal(upper)  ? 1:0;
    }else{
        has_lower   = lower ? 1:0;
        has_diag    = diag  ? 1:0;
        has_upper   = upper ? 1:0;
    }

    nnz = (store_t == MESS_DENSE) ? rows*cols : has_lower*MESS_MIN(rows-1,cols) + has_diag*MESS_MIN(rows,cols)  + has_upper*MESS_MIN(rows,cols-1);

    /*-----------------------------------------------------------------------------
     * reset and alloc matrix
     *-----------------------------------------------------------------------------*/
    MESS_MATRIX_RESET(matrix);
    ret = mess_matrix_alloc(matrix, rows, cols, nnz, store_t, data_t);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    /*-----------------------------------------------------------------------------
     *  fill matrix
     *-----------------------------------------------------------------------------*/
    switch ( store_t ) {
        case MESS_DENSE:
            // fill lower diagonal
            if(has_lower){
                if(MESS_IS_REAL(matrix)){
                    for (i = 0; i<MESS_MIN(rows-1,cols); i++){
                        matrix->values[(i+1)+matrix->ld*i] = creal(lower);
                    }
                }else{
                    for (i = 0; i<MESS_MIN(rows-1,cols); i++){
                        matrix->values_cpx[(i+1)+matrix->ld*i] = lower;
                    }
                }
            }

            // fill diagonal
            if(has_diag){
                if(MESS_IS_REAL(matrix)){
                    for (i = 0; i<MESS_MIN(rows,cols); i++){
                        matrix->values[i+matrix->ld*i] = creal(diag);
                    }
                }else{
                    for (i = 0; i<MESS_MIN(rows,cols); i++){
                        matrix->values_cpx[i+matrix->ld*i] = diag;
                    }
                }
            }

            // fill upper diagonal
            if(has_upper){
                if(MESS_IS_REAL(matrix)){
                    for (i = 0; i<MESS_MIN(rows,cols-1); i++){
                        matrix->values[i+matrix->ld*(i+1)] = creal(upper);
                    }
                }else{
                    for (i = 0; i<MESS_MIN(rows,cols-1); i++){
                        matrix->values_cpx[i+matrix->ld*(i+1)] = upper;
                    }
                }
            }

            break;
        case MESS_CSR:

            matrix->rowptr[0]=0;
            for(j=0; j<matrix->rows;++j){

                matrix->rowptr[j+1] = matrix->rowptr[j];

                //fill lower
                if(has_lower && 1<=j && j-1<matrix->cols){
                    if(MESS_IS_REAL(matrix)){
                        matrix->values[i] = creal(lower);
                    }else{
                        matrix->values_cpx[i] = lower;
                    }
                    matrix->colptr[i] = j-1;
                    i++;
                    matrix->rowptr[j+1]++;
                }

                //fill diag
                if(has_diag && j<matrix->cols){
                    if(MESS_IS_REAL(matrix)){
                        matrix->values[i] = creal(diag);
                    }else{
                        matrix->values_cpx[i] = diag;
                    }
                    matrix->colptr[i] = j;
                    i++;
                    matrix->rowptr[j+1]++;
                }

                //fill upper
                if(has_upper && j < matrix->cols -1 ){
                    if(MESS_IS_REAL(matrix)){
                        matrix->values[i] = creal(upper);
                    }else{
                        matrix->values_cpx[i] = upper;
                    }
                    matrix->colptr[i] = j+1;
                    i++;
                    matrix->rowptr[j+1]++;
                }
            }

            break;
        case MESS_CSC:

            matrix->colptr[0]=0;
            for(j=0;j<matrix->cols;++j){

                matrix->colptr[j+1] = matrix->colptr[j];

                //fill lower
                if(has_lower && j+1<matrix->rows){
                    if(MESS_IS_REAL(matrix)){
                        matrix->values[i] = creal(lower);
                    }else{
                        matrix->values_cpx[i] = lower;
                    }
                    matrix->rowptr[i] = j+1;
                    i++;
                    matrix->colptr[j+1]++;
                }

                //fill diag
                if(has_diag && j<matrix->rows){
                    if(MESS_IS_REAL(matrix)){
                        matrix->values[i] = creal(diag);
                    }else{
                        matrix->values_cpx[i] = diag;
                    }
                    matrix->rowptr[i] = j;
                    i++;
                    matrix->colptr[j+1]++;
                }

                //fill upper
                if(has_upper && 1<=j && j-1<matrix->rows){
                    if(MESS_IS_REAL(matrix)){
                        matrix->values[i] = creal(upper);
                    }else{
                        matrix->values_cpx[i] = upper;
                    }
                    matrix->rowptr[i] = j-1;
                    i++;
                    matrix->colptr[j+1]++;
                }
            }
            break;
        case MESS_COORD:
            j = 0;
            for (i = 0; i<MESS_MIN(rows,cols); i++){

                //fill diagonal
                if(has_diag){
                    if(MESS_IS_REAL(matrix)){
                        matrix->values[j]=creal(diag);
                    }else{
                        matrix->values_cpx[j]=diag;
                    }
                    matrix->rowptr[j] = i;
                    matrix->colptr[j] = i;
                    j++;
                }

                //fill upper
                if(has_upper && i+1<matrix->cols){
                    if(MESS_IS_REAL(matrix)){
                        matrix->values[j]=creal(upper);
                    }else{
                        matrix->values_cpx[j]=upper;
                    }
                    matrix->rowptr[j] = i;
                    matrix->colptr[j] = i+1;
                    j++;
                }

                //fill lower
                if(has_lower && i+1<matrix->rows){
                    if(MESS_IS_REAL(matrix)){
                        matrix->values[j]=creal(lower);
                    }else{
                        matrix->values_cpx[j]=lower;
                    }
                    matrix->rowptr[j] = i+1;
                    matrix->colptr[j] = i;
                    j++;
                }
            }

            break;
        default:
            MSG_ERROR("unknown storage type: %s\n", mess_storage_t_str(store_t));
            return (MESS_ERROR_STORAGETYPE);
            break;
    }   /* -----  end switch  ----- */

    return ret;
}

