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
 * @file lib/matrix/sub.c
 * @brief Get a submatrix from a matrix.
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
#ifdef _OPENMP
#include <omp.h>
#endif


/**
 * @brief Copy a row-block out of a matrix.
 * @param[in] input     input matrix
 * @param[in] srow       input start row of the block
 * @param[in] erow       input end row of the block
 * @param[out] out      output matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_rowsub function gets a row submatrix out of
 * a given matrix. \n
 * It is equal to @matlab @verbatim out=input(srow:erow,:)@endverbatim.
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 */
int  mess_matrix_rowsub ( mess_matrix input, mess_int_t srow, mess_int_t erow, mess_matrix out )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i,j;

    mess_check_nullpointer(input);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(input);

    if ( srow < 0 || srow >= input->rows) {
        MSG_ERROR("srow is out of range\n");
        return(MESS_ERROR_DIMENSION);
    }
    if ( erow < 0 || erow >= input->rows) {
        MSG_ERROR("erow is out of range\n");
        return(MESS_ERROR_DIMENSION);
    }
    if ( erow < srow) {
        MSG_ERROR("erow must be larger or equal to srow\n");
        return(MESS_ERROR_DIMENSION);
    }

    MESS_MATRIX_RESET(out);
    ret = mess_matrix_alloc(out, erow-srow+1, input->cols, 0, MESS_DENSE, input->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    if(MESS_IS_DENSE(input)){
        if ( MESS_IS_REAL(input)) {
            for ( j=0;j<input->cols;j++){
                for ( i = srow; i <=erow; i++){
                    out->values[j*out->ld+(i-srow)] = input->values[j*input->ld+i];
                }
            }
        } else if ( MESS_IS_COMPLEX(input)) {
            for ( j=0;j<input->cols;j++){
                for ( i = srow; i <=erow; i++){
                    out->values_cpx[j*out->ld+(i-srow)] = input->values_cpx[j*input->ld+i];
                }
            }
        }
    }else if(MESS_IS_CSR(input)){
        if(MESS_IS_REAL(input)){
            for(i=srow;i<=erow;++i){
                for(j=input->rowptr[i];j<input->rowptr[i+1];++j){
                    out->values[input->colptr[j]*out->ld+(i-srow)]=input->values[j];
                }
            }
        }else if(MESS_IS_COMPLEX(input)){
            for(i=srow;i<=erow;++i){
                for(j=input->rowptr[i];j<input->rowptr[i+1];++j){
                    out->values_cpx[input->colptr[j]*out->ld+(i-srow)]=input->values_cpx[j];
                }
            }
        }
    }else if(MESS_IS_CSC(input)){
        if(MESS_IS_REAL(input)){
            for(i=0;i<input->cols;++i){
                for(j=input->colptr[i];j<input->colptr[i+1];++j){
                    if(srow<=input->rowptr[j] && input->rowptr[j]<=erow){
                        out->values[i*out->ld+input->rowptr[j]-srow]=input->values[j];
                    }
                }
            }
        }else if(MESS_IS_COMPLEX(input)){
            for(i=0;i<input->cols;++i){
                for(j=input->colptr[i];j<input->colptr[i+1];++j){
                    if(srow<=input->rowptr[j] && input->rowptr[j]<=erow){
                        out->values_cpx[i*out->ld+input->rowptr[j]-srow]=input->values_cpx[j];
                    }
                }
            }
        }

    }else{
        MSG_ERROR("Storagetype not supported!\n");
        return MESS_ERROR_STORAGETYPE;
    }

    return(0);
}       /* -----  end of function mess_matrix_rowsub  ----- */

/**
 * @brief Copy a column block out of a matrix.
 * @param[in] input     input matrix
 * @param[in] scol       input start column of the block
 * @param[in] ecol       input end column of the block
 * @param[out] out      output matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_colsub function gets a column submatrix out of
 * a given matrix. \n
 * It is equal to @matlab  @verbatim out=input(:,scol:ecol) @endverbatim.
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 */
int  mess_matrix_colsub ( mess_matrix input, mess_int_t scol, mess_int_t ecol, mess_matrix out )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i,j;

    mess_check_nullpointer(input);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(input);

    if ( scol == input->cols) {
        MESS_MATRIX_RESET(out);
        return(0);
    }
    if ( scol < 0 || scol > input->cols) {
        MSG_ERROR("scol is out of range ( scol = " MESS_PRINTF_INT ", cols = " MESS_PRINTF_INT " ) \n", scol, input->cols);
        return(MESS_ERROR_DIMENSION);
    }
    if ( ecol < 0 || ecol >= input->cols) {
        MSG_ERROR("ecol is out of range ( ecol = " MESS_PRINTF_INT ", cols = " MESS_PRINTF_INT " ) \n", ecol, input->cols);

        return(MESS_ERROR_DIMENSION);
    }
    if ( ecol < scol) {
        MSG_ERROR("erow must be larger or equal to srow\n");
        return(MESS_ERROR_DIMENSION);
    }

    ret = mess_matrix_alloc(out, input->rows, ecol-scol+1, 0, MESS_DENSE, input->data_type);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    if(MESS_IS_DENSE(input)){
        if ( MESS_IS_REAL(input)) {
            for ( i = scol; i <=ecol; i++){
                for (j=0; j < input->rows; j++) {
                    out->values[(i-scol)*out->ld+j] = input->values[i*input->ld+j];
                }
            }
        } else if (MESS_IS_COMPLEX(input)) {
            for ( i = scol; i <=ecol; i++){
                for (j=0; j < input->rows; j++) {
                    out->values_cpx[(i-scol)*out->ld+j] = input->values_cpx[i*input->ld+j];
                }
            }
        }
    }else if(MESS_IS_CSC(input)){
        if(MESS_IS_REAL(input)){
            for(i=scol;i<=ecol;++i){
                for(j=input->colptr[i];j<input->colptr[i+1];++j){
                    out->values[(i-scol)*out->ld+input->rowptr[j]]=input->values[j];
                }
            }
        }else if(MESS_IS_COMPLEX(input)){
            for(i=scol;i<=ecol;++i){
                for(j=input->rowptr[i];j<input->rowptr[i+1];++j){
                    out->values_cpx[(i-scol)*out->ld+input->rowptr[j]]=input->values_cpx[j];
                }
            }
        }
    }else if(MESS_IS_CSR(input)){
        if(MESS_IS_REAL(input)){
            for(i=0;i<input->rows;++i){
                for(j=input->rowptr[i];j<input->rowptr[i+1];++j){
                    if(scol<=input->colptr[j] && input->colptr[j]<=ecol){
                        out->values[(input->colptr[j]-scol)*out->ld+i]=input->values[j];
                    }
                }
            }
        }else if(MESS_IS_COMPLEX(input)){
            for(i=0;i<input->rows;++i){
                for(j=input->rowptr[i];j<input->rowptr[i+1];++j){
                    if(scol<=input->colptr[j] && input->colptr[j]<=ecol){
                        out->values_cpx[(input->colptr[j]-scol)*out->ld+i]=input->values_cpx[j];
                    }
                }
            }
        }

    }else{
        MSG_ERROR("Storagetype not supported!\n");
        return MESS_ERROR_STORAGETYPE;
    }
    return( 0 ) ;
}       /* ----  end of function mess_matrix_rowsub  ----- */


/**
 * @brief Copy a sub matrix out of a matrix.
 * @param[in] in input matrix
 * @param[in] rowS  input starting row
 * @param[in] rowE  input end row
 * @param[in] colS  input stating column
 * @param[in] colE  input end column
 * @param[out] out ouput block of the matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_sub function gets a sub matrix out of the matrix in. \n
 * It works like  @verbatim out = in (rowS:rowE,colS:colE) @endverbatim in @matlab .
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 */
int  mess_matrix_sub ( mess_matrix in, mess_int_t rowS, mess_int_t rowE, mess_int_t colS, mess_int_t colE, mess_matrix out )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t m,n ;
    mess_int_t mb,nb;
    mess_int_t i, j;


    /*-----------------------------------------------------------------------------
     *  check inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    m = in->rows;
    n = in->cols;

    if ( rowS < 0 || rowE < 0 || colE < 0 || colS< 0) {
        MSG_ERROR("One of the postion arguments is smaller than 0 (rowS=" MESS_PRINTF_INT ", rowE=" MESS_PRINTF_INT ", colS=" MESS_PRINTF_INT ", colE=" MESS_PRINTF_INT ")\n", rowS, rowE, colS, colE);
        return(MESS_ERROR_ARGUMENTS);
    }
    if ( rowS >= m || rowE >= m || colE >= n || colS >= n) {
        MSG_ERROR("One of the postion arguments is out of range (rowS=" MESS_PRINTF_INT ", rowE=" MESS_PRINTF_INT ", colS=" MESS_PRINTF_INT ", colE=" MESS_PRINTF_INT ", m=" MESS_PRINTF_INT ", n=" MESS_PRINTF_INT ")\n", rowS, rowE, colS, colE,m,n);
        return(MESS_ERROR_ARGUMENTS);
    }
    if ( rowS > rowE ) {
        MSG_ERROR("rowS is larger than rowE: ( rowS = " MESS_PRINTF_INT ", rowE = " MESS_PRINTF_INT " ) \n", rowS, rowE);
        return(MESS_ERROR_ARGUMENTS);
    }
    if ( colS > colE ) {
        MSG_ERROR("colS is larger than colE: ( colS = " MESS_PRINTF_INT ", colE = " MESS_PRINTF_INT " ) \n", colS, colE);
        return(MESS_ERROR_ARGUMENTS);
    }

    mb = rowE-rowS+1;
    nb = colE-colS+1;

    ret = mess_matrix_alloc(out, mb, nb, mb*nb, MESS_DENSE, in->data_type);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    if(MESS_IS_DENSE(in)){

        if ( MESS_IS_REAL(in)) {
            for ( j = colS; j <= colE; j++){
                for (i = rowS; i<=rowE; i++) {
                    out->values[(j-colS)*out->ld+(i-rowS)]=in->values[j*in->ld+i];
                }
            }
        } else {
            for ( j = colS; j <= colE; j++){
                for (i = rowS; i<=rowE; i++) {
                    out->values_cpx[(j-colS)*out->ld+(i-rowS)]=in->values_cpx[j*in->ld+i];
                }
            }
        }
    }else if(MESS_IS_CSR(in)){
        if(MESS_IS_REAL(in)){
            for(i=rowS;i<=rowE;++i){
                for(j=in->rowptr[i];j<in->rowptr[i+1];++j){
                    if(colS<=in->colptr[j] && in->colptr[j]<=colE){
                        out->values[(in->colptr[j]-colS)*out->ld+i-rowS]=in->values[j];
                    }
                }
            }
        }else if(MESS_IS_COMPLEX(in)){
            for(i=rowS;i<=rowE;++i){
                for(j=in->rowptr[i];j<in->rowptr[i+1];++j){
                    if(colS<=in->colptr[j] && in->colptr[j]<=colE){
                        out->values_cpx[(in->colptr[j]-colS)*out->ld+i-rowS]=in->values_cpx[j];
                    }
                }
            }
        }

    }else if(MESS_IS_CSC(in)){
        if(MESS_IS_REAL(in)){
            for(i=colS;i<=colE;++i){
                for(j=in->colptr[i];j<in->colptr[i+1];++j){
                    if(rowS<=in->rowptr[j] && in->rowptr[j]<=rowE){
                        out->values[(i-colS)*out->ld+in->rowptr[j]-rowS]=in->values[j];
                    }
                }
            }
        }else if(MESS_IS_COMPLEX(in)){
            for(i=colS;i<=colE;++i){
                for(j=in->colptr[i];j<in->colptr[i+1];++j){
                    if(rowS<=in->rowptr[j] && in->rowptr[j]<=rowE){
                        out->values_cpx[(i-colS)*out->ld+in->rowptr[j]-rowS]=in->values_cpx[j];
                    }
                }
            }
        }
    }else{
        MSG_ERROR("Storagetype not supported!\n");
        return MESS_ERROR_STORAGETYPE;
    }
    return(0);
}       /* -----  end of function mess_matrix_sub  ----- */

