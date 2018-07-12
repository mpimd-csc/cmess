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
 * @file lib/matrix/transpose.c
 * @brief Transpose a matrix.
 * @author @koehlerm
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/eigenvalue.h"
#include "mess/error_macro.h"
#include <complex.h>

/**
 * @brief Transpose a matrix.
 * @param[in] A             input matrix \f$A\f$
 * @param[out] AT           output transposed matrix \f$A^{T}\ (A^{H})\f$
 * @param[in] hermitian     input if true then hermitian transpose is computed
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_xtranspose function creates the transpose/hermitian transpose of a given matrix. \n
 * It supports all internal matrix storage formats.
 */
int mess_matrix_xtranspose (mess_matrix A, mess_matrix AT, int hermitian){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(AT);
    MESS_MATRIX_RESET(AT);


    if ( MESS_IS_CSR(A)) {
        mess_int_t tmp;
        mess_int_t *ptr;
        ret = mess_matrix_convert_csr_csc(A, AT);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_convert_csr_csc);
        // swap row and cols
        tmp = AT->rows;
        AT->rows = AT->cols;
        AT->cols = tmp;
        // swap row ptr;
        ptr = AT->colptr;
        AT->colptr = AT->rowptr;
        AT->rowptr = ptr;
        AT->store_type = MESS_CSR;
        if ( MESS_IS_COMPLEX(AT) && hermitian ) {
            for ( i = 0 ; i < AT->nnz; i++){
                AT->values_cpx[i]=conj(AT->values_cpx[i]);
            }
        }
    } else if ( MESS_IS_CSC(A)) {
        mess_int_t tmp;
        mess_int_t *ptr;
        ret = mess_matrix_convert_csc_csr(A, AT); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert_csc_csr);
        // swap row and cols
        tmp = AT->rows;
        AT->rows = AT->cols;
        AT->cols = tmp;
        // swap row ptr;
        ptr = AT->colptr;
        AT->colptr = AT->rowptr;
        AT->rowptr = ptr;
        AT->store_type = MESS_CSC;
        if ( MESS_IS_COMPLEX(AT) && hermitian ) {
            for ( i = 0 ; i < AT->nnz; i++){
                AT->values_cpx[i]=conj(AT->values_cpx[i]);
            }
        }

    } else if ( MESS_IS_DENSE(A)) {
        // MSG_INFO("-> transpose dense matrix\n");
        mess_check_real_or_complex(A);
        mess_int_t  j;
        ret = mess_matrix_alloc (AT, A->cols, A->rows, A->rows*A->cols, MESS_DENSE, A->data_type);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
        if ( MESS_IS_REAL(A)){
            for ( i = 0; i < A->cols; i++) {
                for (j=0; j < A->rows; j++){
                    AT->values[j*AT->ld+i]=A->values[i*A->ld+j];

                }
            }
        }else {
            if(hermitian)
            for ( i = 0; i < A->cols; i++) {
                for (j=0; j < A->rows; j++){
                    AT->values_cpx[j*AT->ld+i]=conj(A->values_cpx[i*A->ld+j]);
                }
            }
            else{
            for ( i = 0; i < A->cols; i++) {
                for (j=0; j < A->rows; j++){
                    AT->values_cpx[j*AT->ld+i]=A->values_cpx[i*A->ld+j];
                }
            }
            }
        }

    } else if ( MESS_IS_COORD(A)) {
        mess_int_t *ptr, tmp ;
        ret = mess_matrix_copy(A, AT); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);
        ptr = AT->rowptr;
        AT->rowptr = AT->colptr;
        AT->colptr = ptr;
        if ( MESS_IS_COMPLEX(AT) && hermitian ) {
            for ( i = 0 ; i < AT->nnz; i++){
                AT->values_cpx[i]=conj(AT->values_cpx[i]);
            }
        }
        tmp = AT->cols;
        AT->cols = AT->rows;
        AT->rows = tmp;
    } else {
        mess_int_t tmp;
        mess_int_t *ptr;
        mess_matrix work=NULL;
        ret = mess_matrix_init(&work);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_convert(A, work, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
        ret = mess_matrix_convert_csr_csc(A, AT); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert_csr_csc);
        // swap row and cols
        tmp = AT->rows;
        AT->rows = AT->cols;
        AT->cols = tmp;
        // swap row ptr;
        ptr = AT->colptr;
        AT->colptr = AT->rowptr;
        AT->rowptr = ptr;
        AT->store_type = MESS_CSR;
        if ( MESS_IS_COMPLEX(AT) && hermitian ) {
            for ( i = 0 ; i < AT->nnz; i++){
                AT->values_cpx[i]=conj(AT->values_cpx[i]);
            }
        }

        mess_matrix_clear(&work);
    }
    return(0);
}



/**
 * @brief Transpose a matrix.
 * @param[in] A     input matrix \f$A\f$
 * @param[out] AH   output hermitian transposed matrix \f$A^{H}\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_ctranspose function creates the hermitian transpose of a given matrix. \n
 * It supports all internal matrix storage formats
 */
int mess_matrix_ctranspose (mess_matrix A, mess_matrix AH){
   return mess_matrix_xtranspose(A,AH,1);
}


/**
 * @brief Transpose a matrix.
 * @param[in] A     input matrix \f$A\f$
 * @param[out] AT   output transposed matrix \f$A^{T}\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_transpose function creates the transpose of a given matrix. \n
 * It supports all internal matrix storage formats.
 */
int mess_matrix_transpose (mess_matrix A, mess_matrix AT){
   return mess_matrix_xtranspose(A,AT,0);
}


