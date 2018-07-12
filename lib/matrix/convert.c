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
 * @file lib/matrix/convert.c
 * @brief Uniform matrix converting interface.
 * @author @koehlerm
 *
 * This file provides the basic matrix converting routines and calls the user defined functions
 * of a matrix if they are needed.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>

// Local forward definitions
static int __conv_coord_csr ( mess_matrix inmatrix, mess_matrix outmatrix );
static int __conv_coord_csc ( mess_matrix  inmatrix, mess_matrix outmatrix );
static int  __conv_coord_dense ( mess_matrix  inmatrix, mess_matrix outmatrix );
static int __conv_dense_csr(mess_matrix dense, mess_matrix csr);
static int __conv_dense_csc(mess_matrix dense, mess_matrix csc);
static int __conv_csr_dense(mess_matrix csr, mess_matrix dense);
static int __conv_csc_dense(mess_matrix csr, mess_matrix dense);
static int __conv_csr_coord(mess_matrix csr, mess_matrix coord);
static int __conv_csc_coord(mess_matrix csc, mess_matrix coord);
static int __conv_dense_coord(mess_matrix dense, mess_matrix coord);

#include "shellsort.c"


/**
 * @brief Convert the storage type of a matrix (universal converter interface).
 * @param[in] input input matrix
 * @param[out] output output matrix
 * @param[in] outtype   input outputtype
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_convert functions converts a given matrix
 * to an other storage scheme. \n
 * It supports dense and sparse matrices in various formats.
 * If the input matrix has the same storage
 * format than the desired one, it is copied instead of being converted.
 * If a matrix is converted to dense it will have the default leading dimension of
 * @ref mess_matrix_alloc. The support for user supplied matrices is integrated but not used
 * in @mess at the moment.
 *
 * Supported converters:
 * \li CSR \f$ \to \f$ CSC
 * \li CSR \f$ \to \f$ DENSE
 * \li CSR \f$ \to \f$ COORD
 * \li CSC \f$ \to \f$ CSR
 * \li CSC \f$ \to \f$ DENSE
 * \li CSC \f$ \to \f$ COORD
 * \li COORD \f$ \to \f$ CSR
 * \li COORD \f$ \to \f$ CSC
 * \li COORD \f$ \to \f$ DENSE
 * \li DENSE \f$ \to \f$ CSR
 * \li DENSE \f$ \to \f$ CSC
 * \li DENSE \f$ \to \f$ COORD
 *
 * @sa mess_matrix_alloc
 * @sa mess_matrix_read
 * @sa mess_matrix_read_formated
 *
 */
int mess_matrix_convert(mess_matrix input, mess_matrix output, mess_storage_t outtype) {
    MSG_FNAME(__func__);
    unsigned short intype;
    int ret = 0;

    mess_check_nullpointer(input);
    mess_check_nullpointer(output);

    MESS_MATRIX_RESET(output);
    intype = input->store_type;

    if (intype == outtype) {
        // In case of same types copy the matrix
        ret = mess_matrix_copy(input, output);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        return 0;
    }

    switch ( intype ){

        /*-----------------------------------------------------------------------------
         *  convert from COORD
         *-----------------------------------------------------------------------------*/
        case MESS_COORD:{
                            switch (outtype){
                                case MESS_CSR:
                                    ret = __conv_coord_csr(input, output);
                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __conv_coord_csr);
                                    break;
                                case MESS_CSC:
                                    ret =__conv_coord_csc(input, output);
                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __conv_coord_csc);
                                    break;
                                case MESS_DENSE:
                                    ret= __conv_coord_dense(input, output);
                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __conv_coord_dense);
                                    break;
                                default:
                                    MSG_ERROR("convert COORD in format %d(%s) is not supported at the moment\n", outtype, mess_storage_t_str(outtype));
                                    return MESS_ERROR_STORAGETYPE;
                            }
                            return 0;
                        }

                        /*-----------------------------------------------------------------------------
                         *  convert from CSR
                         *-----------------------------------------------------------------------------*/
        case MESS_CSR:{
                          switch (outtype){
                              case MESS_CSC:
                                  ret = mess_matrix_convert_csr_csc(input, output);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert_csr_csc);
                                  break;
                              case MESS_COORD:
                                  ret = __conv_csr_coord(input, output);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __conv_csr_coord);
                                  break;
                              case MESS_DENSE:
                                  ret =__conv_csr_dense(input, output);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __conv_csr_dense);
                                  break;
                              default:
                                  MSG_ERROR("convert CSR in format %d(%s) is not supported at the moment\n", outtype, mess_storage_t_str(outtype));
                                  return MESS_ERROR_STORAGETYPE;
                          }
                          return 0;
                      }

                      /*-----------------------------------------------------------------------------
                       *  convert from CSC
                       *-----------------------------------------------------------------------------*/
        case MESS_CSC:{
                          switch (outtype) {
                              case MESS_CSR:
                                  ret = mess_matrix_convert_csc_csr(input, output);
                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert_csc_csr);
                                  break;
                              case MESS_COORD:
                                  ret = __conv_csc_coord(input, output); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __conv_csc_coord);
                                  break;
                              case MESS_DENSE:
                                  ret = __conv_csc_dense(input, output);
                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __conv_csc_dense);
                                  break;

                              default:
                                  MSG_ERROR("convert CSC in format %d(%s) is not supported at the moment\n", outtype, mess_storage_t_str(outtype));
                                  return MESS_ERROR_STORAGETYPE;
                          }
                          return 0;
                      }

                      /*-----------------------------------------------------------------------------
                       *  convert from DENSE
                       *-----------------------------------------------------------------------------*/
        case MESS_DENSE:{
                            switch ( outtype) {
                                case MESS_CSR:
                                    ret= __conv_dense_csr(input, output);
                                    FUNCTION_FAILURE_HANDLE( ret, (ret!=0), __conv_dense_csr);
                                    break;
                                case MESS_CSC:
                                    ret = __conv_dense_csc(input, output);
                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __conv_dense_csc);
                                    break;
                                case MESS_COORD:
                                    ret = __conv_dense_coord(input, output);
                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __conv_dense_coord);
                                    break;
                                default:
                                    MSG_ERROR("convert DENSE in format %d(%s) is not supported at the moment\n", outtype, mess_storage_t_str(outtype));
                                    return MESS_ERROR_STORAGETYPE;
                            }
                            return 0;
                        }
        default:
                        MSG_ERROR("unknown input type: %s\n", mess_storage_t_str(intype));
                        return MESS_ERROR_STORAGETYPE;
    }

    return 0;

}

/**
 * @internal
 * @brief Convert a coordinate matrix to CSR.
 * @param[in] inmatrix      input coordinate matrix
 * @param[out] outmatrix    output Compressed Sparse Row matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __conv_coord_csr function converts a matrix from coordinate to
 * Compressed Sparse Row storage.
 *
 * @attention Internal use only.
 */
static int  __conv_coord_csr ( mess_matrix inmatrix, mess_matrix outmatrix )
{
    MSG_FNAME(__func__);
    mess_int_t rows, cols, nnz, i, pos;
    mess_int_t *rowcount;
    int ret = 0;

    mess_check_nullpointer(inmatrix);
    mess_check_nullpointer(outmatrix);
    mess_check_coord(inmatrix);


    rows = inmatrix->rows;
    cols = inmatrix->cols;
    nnz  = inmatrix->nnz;
    if ( MESS_IS_INTEGER(inmatrix)) inmatrix->data_type = MESS_REAL;

    ret = mess_matrix_alloc(outmatrix, rows, cols, nnz, MESS_CSR, inmatrix->data_type);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    mess_try_alloc(rowcount, mess_int_t *, rows *sizeof( mess_int_t));
    for (i = 0; i < rows; i++) rowcount[i] = 0;
    // count entries in rows
    for ( i = 0 ; i < nnz; i++) rowcount[inmatrix->rowptr[i]]++;

    // fill rowptr
    outmatrix->rowptr[0] = 0;
    for ( i = 0; i < rows; i++){
        outmatrix->rowptr[i+1] = outmatrix->rowptr[i] + rowcount[i];
        rowcount[i] = outmatrix->rowptr[i];
    }
    // fill values
    for ( i = 0; i < nnz; i++){
        pos =  rowcount[inmatrix->rowptr[i]];
        outmatrix->colptr[pos] = inmatrix->colptr[i];


        if ( MESS_IS_REAL (inmatrix)){
            outmatrix->values[pos] = inmatrix->values[i];
        }

        if ( MESS_IS_COMPLEX(inmatrix)){
            outmatrix->values_cpx[pos] = inmatrix->values_cpx[i];
        }
        rowcount[inmatrix->rowptr[i]]++;
    }
    outmatrix->rows  = rows;
    outmatrix->cols  = cols;
    outmatrix->nnz = nnz;
    outmatrix->symmetry = inmatrix->symmetry;
    outmatrix->data_type = inmatrix->data_type;
    outmatrix->store_type = MESS_CSR;

    // Sort rows
    for ( i = 0; i < rows; i++){
        __shellsort(outmatrix->colptr, outmatrix->values, outmatrix->values_cpx, outmatrix->rowptr[i], outmatrix->rowptr[i+1]-1, outmatrix->data_type);
    }
    mess_free(rowcount);
    return 0;
}


/**
 * @internal
 * @brief Convert a coordinate matrix to CSC.
 * @param[in] inmatrix      input coordinate matrix
 * @param[out] outmatrix    output Compressed Sparse Column matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __conv_coord_csc function converts a matrix from coordinate to
 * Compressed Sparse Column storage. It uses the property that the CSC storage
 * is the transpose of the CSR storage. In this way the matrix is transposed
 * in the coordinate format and then converted to Compressed Sparse Row storage.
 * Afterwards the data structure is modified to fit in a Compressed Sparse
 * Column storage.
 *
 * @attention Internal use only.
 */
static int  __conv_coord_csc ( mess_matrix  inmatrix, mess_matrix outmatrix )
{
    MSG_FNAME(__func__);
    mess_matrix in_tmp;
    mess_matrix csr_tmp;
    int ret = 0;

    mess_check_nullpointer(inmatrix);
    mess_check_nullpointer(outmatrix);

    ret =   mess_matrix_init(&in_tmp);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret =   mess_matrix_init(&csr_tmp);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    // transpose
    in_tmp->store_type      = inmatrix->store_type;
    in_tmp->data_type   = inmatrix->data_type;
    in_tmp->symmetry    = inmatrix->symmetry;
    in_tmp->cols    = inmatrix->rows;
    in_tmp->rows    = inmatrix->cols;
    in_tmp->nnz     = inmatrix->nnz;
    in_tmp->rowptr  = inmatrix->colptr;
    in_tmp->colptr  = inmatrix->rowptr;
    if (MESS_IS_REAL(inmatrix))
        in_tmp->values  = inmatrix->values;
    if (MESS_IS_COMPLEX(inmatrix))
        in_tmp->values_cpx = inmatrix->values_cpx;

    // convert
    ret = __conv_coord_csr(in_tmp, csr_tmp);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __conv_coord_csr);

    // transpose
    outmatrix->store_type = MESS_CSC;
    outmatrix->cols = csr_tmp->rows;
    outmatrix->rows = csr_tmp->cols;
    outmatrix->nnz = csr_tmp->nnz;
    outmatrix->data_type = csr_tmp->data_type;
    outmatrix->symmetry =  csr_tmp->symmetry;
    outmatrix->rowptr = csr_tmp->colptr;
    outmatrix->colptr = csr_tmp->rowptr;
    if ( MESS_IS_REAL (csr_tmp))
        outmatrix->values = csr_tmp->values;
    if ( MESS_IS_COMPLEX(csr_tmp)){
        outmatrix->values_cpx = csr_tmp->values_cpx;
    }

    // only free not clear.
    mess_free(in_tmp);
    mess_free(csr_tmp);
    return 0;
}


/**
 *
 * @brief Convert a Compressed Sparse Row matrix to a Compressed Sparse Column matrix.
 * @param[in] inmatrix input matrix
 * @param[out] outmatrix output matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_convert_csr_csc converts a Compressed Sparse Row
 * matrix to a Compressed Sparse Column matrix.
 *
 * @sa mess_matrix_convert_csc_csr
 * @sa mess_matrix_convert
 */
int mess_matrix_convert_csr_csc(mess_matrix inmatrix, mess_matrix outmatrix){
    MSG_FNAME(__func__);
    mess_int_t rows, cols, nnz, i, pos, j;
    mess_int_t *colcount;
    int ret = 0;

    // MSG_WARN("CSR->CSC\n");

    mess_check_nullpointer(inmatrix);
    mess_check_nullpointer(outmatrix);

    if ( !MESS_IS_CSR(inmatrix) ) return MESS_ERROR_STORAGETYPE;
    MESS_MATRIX_RESET(outmatrix);

    rows = inmatrix->rows;
    cols = inmatrix->cols;
    nnz  = inmatrix->nnz;

    ret = mess_matrix_alloc(outmatrix, rows, cols, nnz, MESS_CSC, inmatrix->data_type);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    mess_try_alloc(colcount,mess_int_t *,cols*sizeof(mess_int_t));
    for ( i = 0; i < cols ; i++ ) { colcount[i] = 0; }
    // count entries in rows
    for ( i = 0 ; i < nnz; i++) colcount[inmatrix->colptr[i]]++;

    // fill rowptr
    outmatrix->colptr[0] = 0;
    for ( i = 0; i < cols; i++){
        outmatrix->colptr[i+1] = outmatrix->colptr[i] + colcount[i];
        colcount[i] = outmatrix->colptr[i];
    }
    // fill values
    for ( i = 0; i < inmatrix->rows; i++){
        for ( j = inmatrix->rowptr[i]; j < inmatrix->rowptr[i+1]; j++){
            pos = colcount[inmatrix->colptr[j]]++;
            outmatrix->rowptr[pos] = i;
            if ( MESS_IS_REAL (inmatrix)){
                outmatrix->values[pos] = inmatrix->values[j];
            }
            if ( MESS_IS_COMPLEX(inmatrix)){
                outmatrix->values_cpx[pos] = inmatrix->values_cpx[j];
            }
        }
    }
    outmatrix->rows  = rows;
    outmatrix->cols  = cols;
    outmatrix->nnz = nnz;
    outmatrix->symmetry = inmatrix->symmetry;
    outmatrix->data_type = inmatrix->data_type;
    outmatrix->store_type = MESS_CSC;

    mess_free(colcount);
    return 0;
}

/**
 * @brief Convert a Compressed Sparse Column matrix to a Compressed Sparse Row matrix.
 * @param[in]  inmatrix input matrix
 * @param[out] outmatrix output matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_convert_csc_csr converts a Compressed Sparse Column
 * matrix to a Compressed Sparse Row matrix.
 *
 * @sa mess_matrix_convert_csr_csc
 * @sa mess_matrix_convert
 */
int mess_matrix_convert_csc_csr(mess_matrix inmatrix, mess_matrix outmatrix){
    MSG_FNAME(__func__);
    mess_int_t rows, cols, nnz, i, pos, j;
    mess_int_t *rowcount;
    int ret = 0;

    mess_check_nullpointer(inmatrix);
    mess_check_nullpointer(outmatrix);

    if ( !MESS_IS_CSC(inmatrix)) return MESS_ERROR_STORAGETYPE;
    MESS_MATRIX_RESET(outmatrix);

    rows = inmatrix->rows;
    cols = inmatrix->cols;
    nnz  = inmatrix->nnz;


    ret = mess_matrix_alloc(outmatrix, rows, cols, nnz, MESS_CSR, inmatrix->data_type);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    mess_try_alloc(rowcount, mess_int_t *,rows*sizeof(mess_int_t));
    for ( i= 0; i < rows; i++) rowcount[i] = 0;
    // count entries in rows
    for ( i = 0 ; i < nnz; i++) rowcount[inmatrix->rowptr[i]]++;

    // fill rowptr
    outmatrix->rowptr[0] = 0;
    for ( i = 0; i < rows; i++){
        outmatrix->rowptr[i+1] = outmatrix->rowptr[i] + rowcount[i];
        rowcount[i] = outmatrix->rowptr[i];
    }
    // fill values
    for ( i = 0; i < inmatrix->cols; i++){
        for ( j = inmatrix->colptr[i]; j < inmatrix->colptr[i+1]; j++){
            // MSG_PRINT ("rowptr[%lu] = %lu\n", j, inmatrix->rowptr[j]);
            pos = rowcount[inmatrix->rowptr[j]];
            outmatrix->colptr[pos] = i;
            if ( MESS_IS_REAL (inmatrix) ){
                outmatrix->values[pos] = inmatrix->values[j];
            }
            if ( MESS_IS_COMPLEX(inmatrix)){
                outmatrix->values_cpx[pos] = inmatrix->values_cpx[j];
            }
            rowcount[inmatrix->rowptr[j]]++;
        }
    }
    outmatrix->rows  = rows;
    outmatrix->cols  = cols;
    outmatrix->nnz = nnz;
    outmatrix->symmetry = inmatrix->symmetry;
    outmatrix->data_type = inmatrix->data_type;
    outmatrix->store_type = MESS_CSR;

    mess_free(rowcount);
    return 0;
}




/**
 * @internal
 * @brief Convert a dense matrix to a Compressed Sparse Row matrix.
 * @param[in] dense input matrix
 * @param[out]  csr   output matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __conv_dense_csr converts a Dense matrix
 * matrix to a Compressed Sparse Row matrix.
 *
 * @attention Internal use only.
 */
static int __conv_dense_csr(mess_matrix dense, mess_matrix csr){
    MSG_FNAME(__func__);
    mess_int_t nnz = 0;
    mess_int_t i=0, j=0;
    mess_int_t pos = 0;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(dense);
    mess_check_nullpointer(csr);
    mess_check_real_or_complex(dense);
    mess_check_dense(dense);

    // count nonnz
    if ( MESS_IS_REAL(dense) ){
        for (j = 0; j < dense->cols; j++){
            for ( i = 0; i < dense->rows ; i++) {
                if ( dense->values[i+j*dense->ld] != 0.0) nnz++;
            }
        }
    } else  if ( MESS_IS_COMPLEX(dense)){
        for (j = 0; j < dense->cols; j++){
            for ( i = 0; i < dense->rows ; i++) {
                if ( dense->values_cpx[i+j*dense->ld] != 0.0) nnz++;
            }
        }
    }

    if ( MESS_IS_REAL( dense )) {
        ret = mess_matrix_alloc(csr, dense->rows, dense->cols, nnz, MESS_CSR, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0) , mess_matrix_alloc);
        csr->rowptr[0]= 0;
        for ( i = 0; i < dense->rows; i++ ){
            for ( j = 0; j < dense->cols; j++){
                if ( dense->values[j*dense->ld+i] != 0) {
                    csr->values[pos] = dense->values[j*dense->ld+i];
                    csr->colptr[pos] = j;
                    pos++;
                }
            }
            csr->rowptr[i+1] = pos;
        }
    } else if ( MESS_IS_COMPLEX(dense)) {
        ret = mess_matrix_alloc(csr, dense->rows, dense->cols, nnz, MESS_CSR, MESS_COMPLEX);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        csr->rowptr[0]= 0;
        for ( i = 0; i < dense->rows; i++ ){
            for ( j = 0; j < dense->cols; j++){
                if ( dense->values_cpx[j*dense->ld+i] != 0) {
                    csr->values_cpx[pos] = dense->values_cpx[j*dense->ld+i];
                    csr->colptr[pos] = j;
                    pos++;
                }
            }
            csr->rowptr[i+1] = pos;
        }
    }

    return 0;
}

/**
 * @internal
 * @brief Convert a Dense matrix to a Compressed Sparse Column  matrix.
 * @param[in] dense input matrix
 * @param[out] csc   output matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __conv_dense_csc converts a Dense matrix
 * matrix to a Compressed Sparse Column one.
 *
 * @attention Internal use only.
 */
static int __conv_dense_csc(mess_matrix dense, mess_matrix csc){
    MSG_FNAME(__func__);
    mess_int_t nnz = 0;
    mess_int_t i=0, j=0;
    mess_int_t pos = 0;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(dense);
    mess_check_nullpointer(csc);
    mess_check_real_or_complex(dense);
    mess_check_dense(dense);

    // count nonnz
    if ( MESS_IS_REAL(dense) ){
        for (j = 0; j < dense->cols; j++){
            for ( i = 0; i < dense->rows ; i++) {
                if ( dense->values[i+j*dense->ld] != 0.0) nnz++;
            }
        }
    } else  if ( MESS_IS_COMPLEX(dense)){
        for (j = 0; j < dense->cols; j++){
            for ( i = 0; i < dense->rows ; i++) {
                if ( dense->values_cpx[i+j*dense->ld] != 0.0) nnz++;
            }
        }
    }

    if ( MESS_IS_REAL( dense )) {
        ret = mess_matrix_alloc(csc, dense->rows, dense->cols, nnz, MESS_CSC, MESS_REAL);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        csc->colptr[0]= 0;
        for ( j = 0 ; j < dense->cols; j++){
            for ( i = 0; i < dense->rows; i++){
                if ( dense->values[j*dense->ld+i]!=0){
                    csc->values[pos] = dense->values[j*dense->ld+i];
                    csc->rowptr[pos] = i;
                    pos++;
                }
            }
            csc->colptr[j+1] = pos;
        }
    } else if ( MESS_IS_COMPLEX(dense)) {
        ret = mess_matrix_alloc(csc, dense->rows, dense->cols, nnz, MESS_CSC, MESS_COMPLEX);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        csc->colptr[0]= 0;
        for ( j = 0 ; j < dense->cols; j++){
            for ( i = 0; i < dense->rows; i++){
                if ( dense->values_cpx[j*dense->ld+i]!=0){
                    csc->values_cpx[pos] = dense->values_cpx[j*dense->ld+i];
                    csc->rowptr[pos] = i;
                    pos++;
                }
            }
            csc->colptr[j+1] = pos;
        }
    }

    return 0;
}


/**
 * @internal
 * @brief Convert a Compressed Sparse Row matrix to a Dense matrix.
 * @param[in] csr   input matrix
 * @param[out] dense output matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __conv_csr_dense converts a Compressed Sparse Row
 * matrix to a Dense one.
 *
 * @attention Internal use only.
 */
static int __conv_csr_dense(mess_matrix csr, mess_matrix dense){
    mess_int_t i, j;
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(csr);
    mess_check_nullpointer(dense);
    if (!MESS_IS_CSR(csr)) {
        MSG_ERROR("input is not CSR\n");
        return MESS_ERROR_STORAGETYPE;
    }
    if ( MESS_IS_REAL(csr)){
        ret = mess_matrix_alloc(dense, csr->rows, csr->cols, 0, MESS_DENSE, MESS_REAL);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_zeros(dense); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);
        for ( i = 0; i < csr->rows; i++){
            for ( j = csr->rowptr[i]; j < csr->rowptr[i+1]; j++){
                dense->values[dense->ld*csr->colptr[j]+i] = csr->values[j];
            }
        }
    } else if (MESS_IS_COMPLEX(csr)) {
        ret = mess_matrix_alloc(dense, csr->rows, csr->cols, 0, MESS_DENSE, MESS_COMPLEX);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_zeros(dense); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);
        for ( i = 0; i < csr->rows; i++){
            for ( j = csr->rowptr[i]; j < csr->rowptr[i+1]; j++){
                dense->values_cpx[dense->ld*csr->colptr[j]+i] = csr->values_cpx[j];
            }
        }
    } else {
        MSG_ERROR("unknown/unsupported datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}

/**
 * @internal
 * @brief Convert a Compressed Sparse Column matrix to a Dense matrix.
 * @param[in] csc   input matrix
 * @param[out] dense output matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __conv_csc_dense converts a Compressed Sparse Column
 * matrix to a Dense one.
 *
 * @attention Internal use only.
 */
static int __conv_csc_dense(mess_matrix csc, mess_matrix dense){
    mess_int_t i, j;
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(csc);
    mess_check_nullpointer(dense);

    if (!MESS_IS_CSC(csc)) {
        MSG_ERROR("input is not CSR\n");
        return MESS_ERROR_STORAGETYPE;
    }
    if ( MESS_IS_REAL(csc)){
        ret = mess_matrix_alloc(dense, csc->rows, csc->cols, 0, MESS_DENSE, MESS_REAL);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        for (j=0; j < csc->cols;j++){
            for (i=0; i < csc->rows;i++){
                dense->values[dense->ld*j+i] = 0;
            }
        }
        for ( i = 0; i < csc->cols; i++){
            for ( j = csc->colptr[i]; j < csc->colptr[i+1]; j++){
                dense->values[dense->ld*i+csc->rowptr[j]] = csc->values[j];
            }
        }
    } else if (MESS_IS_COMPLEX(csc)) {
        ret = mess_matrix_alloc(dense, csc->rows, csc->cols, 0, MESS_DENSE, MESS_COMPLEX);
        FUNCTION_FAILURE_HANDLE( ret, (ret!=0), mess_matrix_alloc);
        for (j=0; j < csc->cols;j++){
            for (i=0; i < csc->rows;i++){
                dense->values_cpx[dense->ld*j+i] = 0;
            }
        }
        for ( i = 0; i < csc->cols; i++){
            for ( j = csc->colptr[i]; j < csc->colptr[i+1]; j++){
                dense->values_cpx[dense->ld*i+csc->rowptr[j]] = csc->values_cpx[j];
            }
        }
    } else {
        MSG_ERROR("unknown/unsupported datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}


/**
 * @internal
 * @brief Convert a Coordinate Matrix to a Dense matrix.
 * @param[in]  inmatrix   input matrix
 * @param[out] outmatrix     output matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __conv_coord_dense converts a Coordinate
 * matrix to a Dense one.
 *
 * @attention Internal use only.
 */
static int  __conv_coord_dense ( mess_matrix  inmatrix, mess_matrix outmatrix )
{
    MSG_FNAME(__func__);
    mess_int_t i;
    int ret;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(inmatrix);
    mess_check_nullpointer(outmatrix);
    mess_check_real_or_complex(inmatrix);
    mess_check_coord(inmatrix);


    ret = mess_matrix_alloc (outmatrix, inmatrix->rows, inmatrix->cols, inmatrix->rows * inmatrix->cols, MESS_DENSE, inmatrix->data_type);
    FUNCTION_FAILURE_HANDLE( ret, (ret!=0), mess_matrix_alloc);

    if ( MESS_IS_REAL(inmatrix)) {
        for (i=0;i < inmatrix->nnz;i++){
            outmatrix->values[inmatrix->colptr[i]*outmatrix->ld+inmatrix->rowptr[i]]=inmatrix->values[i];
        }
    } else if (MESS_IS_COMPLEX(inmatrix)){
        for (i=0;i < inmatrix->nnz;i++){
            outmatrix->values_cpx[inmatrix->colptr[i]*outmatrix->ld+inmatrix->rowptr[i]]=inmatrix->values_cpx[i];
        }
    } else {
        MSG_ERROR(" only for real and complex matrices.\n");
        MESS_MATRIX_RESET(outmatrix);
    }
    return 0;
}


/**
 * @internal
 * @brief Convert a Compressed Sparse Row matrix to a Coordinate matrix.
 * @param[in] csr   input matrix
 * @param[out] coord  output matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __conv_csr_coord converts a Compressed Sparse Row
 * matrix to a Coordinate  one.
 *
 * @attention Internal use only.
 */
static int __conv_csr_coord(mess_matrix csr, mess_matrix coord){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i, j;

    mess_check_nullpointer(csr);
    mess_check_nullpointer(coord);
    mess_check_real_or_complex(csr);

    MESS_MATRIX_RESET(coord);
    ret = mess_matrix_alloc(coord, csr->rows, csr->cols, csr->nnz, MESS_COORD, csr->data_type);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    if ( MESS_IS_REAL(csr) ) {
        for ( i = 0 ; i < csr->rows; i++) {
            for ( j = csr ->rowptr[i]; j < csr->rowptr[i+1]; j++){
                coord->rowptr[j] = i;
                coord->colptr[j] = csr->colptr[j];
                coord->values[j] = csr->values[j];
            }
        }
    } else if ( MESS_IS_COMPLEX(csr)) {
        for ( i = 0 ; i < csr->rows; i++) {
            for ( j = csr ->rowptr[i]; j < csr->rowptr[i+1]; j++){
                coord->rowptr[j] = i;
                coord->colptr[j] = csr->colptr[j];
                coord->values_cpx[j] = csr->values_cpx[j];
            }
        }

    }
    return 0;
}

/**
 * @internal
 * @brief Convert a Compressed Sparse Column  matrix to a Coordinate matrix.
 * @param[in] csc   input matrix
 * @param[out] coord  output matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __conv_csr_coord converts a Compressed Sparse Column
 * matrix to a Coordinate  one.
 *
 * @attention Internal use only.
 */
static int __conv_csc_coord(mess_matrix csc, mess_matrix coord){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i, j;

    mess_check_nullpointer(csc);
    mess_check_nullpointer(coord);
    mess_check_real_or_complex(csc);

    MESS_MATRIX_RESET(coord);
    ret = mess_matrix_alloc(coord, csc->rows, csc->cols, csc->nnz, MESS_COORD, csc->data_type);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    if ( MESS_IS_REAL(csc) ) {
        for ( i = 0 ; i < csc->cols; i++) {
            for ( j = csc ->colptr[i]; j < csc->colptr[i+1]; j++){
                coord->rowptr[j] = csc->rowptr[j];
                coord->colptr[j] = i;
                coord->values[j] = csc->values[j];
            }
        }
    } else if ( MESS_IS_COMPLEX(csc)) {
        for ( i = 0 ; i < csc->cols; i++) {
            for ( j = csc ->colptr[i]; j < csc->colptr[i+1]; j++){
                coord->rowptr[j] = csc->rowptr[j];
                coord->colptr[j] = i;
                coord->values_cpx[j] = csc->values_cpx[j];
            }
        }
    }
    return 0;
}


/**
 * @internal
 * @brief Convert a Dense matrix to a Coordinate matrix.
 * @param[in] dense   input matrix
 * @param[out] coord  output matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __conv_dense_coord converts a Dense
 * matrix to a Coordinate  one.
 *
 * @attention Internal use only.
 */
static int __conv_dense_coord(mess_matrix dense, mess_matrix coord){
    MSG_FNAME(__func__);
    mess_int_t nnz = 0;
    mess_int_t i=0, j=0;
    mess_int_t pos = 0;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(dense);
    mess_check_nullpointer(coord);
    mess_check_real_or_complex(dense);
    mess_check_dense(dense);

    MESS_MATRIX_RESET(coord);
    // count nonnz
    if ( MESS_IS_REAL(dense) ){
        for (j = 0; j < dense->cols; j++){
            for ( i = 0; i < dense->rows ; i++) {
                if ( dense->values[i+j*dense->ld] != 0.0) nnz++;
            }
        }
    } else  if ( MESS_IS_COMPLEX(dense)){
        for (j = 0; j < dense->cols; j++){
            for ( i = 0; i < dense->rows ; i++) {
                if ( dense->values_cpx[i+j*dense->ld] != 0.0) nnz++;
            }
        }
    }

    ret = mess_matrix_alloc(coord, dense->rows, dense->cols, nnz, MESS_COORD, dense->data_type);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    if ( MESS_IS_REAL( dense )) {
        for ( j = 0 ; j < dense->cols; j++){
            for ( i = 0; i < dense->rows; i++){
                if ( dense->values[j*dense->ld+i]!=0){
                    coord->values[pos] = dense->values[j*dense->ld+i];
                    coord->rowptr[pos] = i;
                    coord->colptr[pos] = j;
                    pos++;
                }
            }
        }
    } else if ( MESS_IS_COMPLEX(dense)) {
        for ( j = 0 ; j < dense->cols; j++){
            for ( i = 0; i < dense->rows; i++){
                if ( dense->values_cpx[j*dense->ld+i]!=0){
                    coord->values_cpx[pos] = dense->values_cpx[j*dense->ld+i];
                    coord->rowptr[pos] = i;
                    coord->colptr[pos] = j;
                    pos++;
                }
            }
        }
    }

    return 0;
}

