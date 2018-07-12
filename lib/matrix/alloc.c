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
 * @file lib/matrix/alloc.c
 * @brief Allocate and resize matrices.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


#define NEW_LD

#ifdef LD_MULTPLICITY
static mess_int_t ld_multiplicity = LD_MULTIPLICITY;
#else
static mess_int_t ld_multiplicity = 1;
#endif

/**
 * @internal
 * @brief Sets the multiplicity of the leading dimension for dense matrices.
 * @param[in] ld_m  Multiplicity of the leading dimension.
 *
 * The @ref __mess_matrix_alloc_setld function sets the multiplicity of the leading dimension of
 * a matrix. Everytime a dense matrix gets allocated by \ref mess_matrix_alloc the leading
 * dimension is selected to be a multiple of @c ld_m. If ld_m <= 0 the multiplicity is set 1.
 * Changing the multiplicity of the leading dimension may influence the behaviour of level-3
 * @blas routines.
 * @attention Internal use only.
 */
void __mess_matrix_alloc_setld( mess_int_t ld_m)
{
    if ( ld_m > 0 )
        ld_multiplicity = ld_m;
    else
        ld_multiplicity = 1;
}

static mess_int_t __select_ld(mess_int_t rows, size_t es){
    mess_int_t ld;
#ifdef NEW_LD
    if ( rows % ld_multiplicity == 0 ) {
        ld = rows;
    } else {
        ld = (rows/ld_multiplicity+1)*ld_multiplicity;
    }
#else
    ld = rows;
#endif
    return ld;

}

/**
 * @brief Check if a new allocation of a matrix is necessary.
 * @param[in] matrix     input matrix to check
 * @param[in] rows         input number of rows
 * @param[in] cols        input number of cols
 * @param[in] nnz        input number of nonzero elements (irrelevant for dense matrices)
 * @param[in] storetype   input type of the storage format
 * @param[in] datatype      input data type of the matrix
 * @return 1 if the allocation is necessary, 0 otherwise.
 *
 * The @ref mess_matrix_need_alloc function checks if a new allocation of a matrix is necessary
 * or not. That means the function checks if the matrix has already the internal structure that
 * the corresponding \ref mess_matrix_alloc call will return. \n
 * It can be used to skip an allocation if the matrix was already allocated previously. \n
 * It is used inside \ref mess_matrix_alloc to reduce the number of allocations.
 *
 * @sa mess_matrix_alloc
 */
int mess_matrix_need_alloc(  mess_matrix matrix,
        mess_int_t rows,
        mess_int_t cols,
        mess_int_t nnz,
        mess_storage_t storetype,
        mess_datatype_t datatype)
{
    if (matrix == NULL) return 1;
    if (matrix->rows != rows) return 1;
    if (matrix->cols != cols) return 1;
    if (matrix->nnz != nnz) return 1;
    if (matrix->store_type != storetype) return 1;
    if (matrix->data_type  != datatype) return 1;
    return 0;
}


#ifdef MESS_DEBUG
static long count_alloc = 0;
static void print_alloc(){
    fprintf(stderr, "+++++++++ DEBUG +++++++++\n");
    fprintf(stderr, "Number of mess_matrix_alloc calls: %ld\n",count_alloc);
}
#endif
/**
 * @brief Allocate a sparse or dense  \f$ (m \times n) \f$ matrix.
 * @param[out] matrix       output matrix to allocate
 * @param[in] rows           input number of rows
 * @param[in] cols          input number of cols
 * @param[in] nnz          input number of nonzero elements (irrelevant for dense matrices)
 * @param[in] storetype   input type of the storage format
 * @param[in] datatype      input data type of the matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_alloc function allocates a matrix. \n
 * It is possible to allocate dense, compressed sparse row, compressed sparse column and coordinate matrices. \n
 * This function checks if the matrix is already allocated and has the right size. If this is the case the matrix
 * is only set to zero. Otherwise all data structures are freed and allocated correctly. \n
 * If @mess is configured in debug mode the allocations are counted and printed to stderr.  \n
 * In case you allocate a dense matrix the leading dimension is set to the number
 * of rows. This can change in the future.
 *
 */
int mess_matrix_alloc(  mess_matrix matrix,
        mess_int_t rows,
        mess_int_t cols,
        mess_int_t nnz,
        mess_storage_t storetype,
        mess_datatype_t  datatype)
{
    MSG_FNAME ( __func__ ) ;
    mess_int_t i ;
    mess_int_t ld = 0;
    int ret = 0;
#ifdef MESS_DEBUG
    static int first = 1;
#endif
    /*-----------------------------------------------------------------------------
     *  check inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    if ( rows < 0 || cols < 0 || nnz < 0 ){
        MSG_ERROR("wrong arguments: rows = " MESS_PRINTF_INT ", cols = " MESS_PRINTF_INT ", nnz = " MESS_PRINTF_INT "\n", rows, cols, nnz);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  Check if the matrix is already correctly allocated. If it is
     *  return imediately
     *-----------------------------------------------------------------------------*/
#ifdef MESS_DEBUG
    if (first) {
        atexit(print_alloc);
        first = 0;
    }
#endif
    if ( !mess_matrix_need_alloc(matrix, rows, cols, nnz, storetype, datatype) && (matrix->values!= NULL || matrix->values_cpx!=NULL)) {
        ret = mess_matrix_zeros(matrix);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);
        return 0;
    }
#ifdef MESS_DEBUG
    count_alloc++;
#endif
    if ( matrix->values != NULL) {    mess_free( matrix->values); matrix->values = NULL; }
    if ( matrix->values_cpx != NULL)  {    mess_free( matrix->values_cpx); matrix->values_cpx = NULL; }
    if ( matrix->colptr != NULL)  {    mess_free( matrix->colptr); matrix->colptr = NULL; }
    if ( matrix->rowptr != NULL)  {    mess_free( matrix->rowptr); matrix->rowptr = NULL; }
    matrix->data_type = 0;
    matrix->nnz = 0;
    matrix->cols = 0;
    matrix->rows = 0;
    matrix->store_type =0;
    matrix->symmetry = 0;

    switch (storetype) {
        case MESS_CSC:
            mess_try_calloc(matrix->rowptr, mess_int_t*, nnz+1, sizeof(mess_int_t));
            mess_try_calloc(matrix->colptr, mess_int_t*, cols+1, sizeof(mess_int_t));
            for ( i = 0 ; i < matrix->cols+1; i++) {
                matrix->colptr[i] = 0;
            }

            if ( datatype == MESS_REAL ){
                mess_try_alloc(matrix->values, double*, sizeof(double)*nnz+1);
                for (i=0; i< nnz;i++) matrix->values[i]= 0.0;
            }
            if ( datatype == MESS_COMPLEX ){
                mess_try_alloc(matrix->values_cpx, mess_double_cpx_t*, sizeof (mess_double_cpx_t) * (nnz)+1);
                for (i=0; i< nnz;i++) matrix->values_cpx[i]= 0.0;
            }
            matrix->colptr[cols] = nnz;
            break;

        case MESS_CSR:
            matrix->rowptr = NULL;
            matrix->colptr = NULL;
            mess_try_calloc(matrix->rowptr, mess_int_t*,(rows+1),  sizeof(mess_int_t));
            mess_try_calloc(matrix->colptr, mess_int_t*, nnz+1,  sizeof(mess_int_t));
            for ( i = 0 ; i < matrix->rows+1; i++) {
                matrix->rowptr[i] = 0;
            }
            if ( datatype == MESS_REAL ) {
                mess_try_alloc(matrix->values, double*, sizeof(double)*(nnz));
                for (i=0; i< nnz;i++) matrix->values[i]= 0.0;
            }
            if ( datatype == MESS_COMPLEX ){
                mess_try_alloc(matrix->values_cpx, mess_double_cpx_t*, sizeof (mess_double_cpx_t) * (nnz));
                for (i=0; i< nnz;i++) matrix->values_cpx[i]= 0.0;
            }
            matrix->rowptr[rows]=nnz;
            break;

        case MESS_DENSE:
            if (datatype == MESS_REAL) {
                ld = __select_ld(rows, sizeof(double));
            } else {
                ld = __select_ld(rows, sizeof(mess_double_cpx_t));
            }
            nnz = rows*cols;
            if ( datatype == MESS_REAL ) {
                if ( ld * cols == 0 ) {
                    matrix->values = NULL;
                } else {
                    mess_try_alloc(matrix->values, double *, sizeof (double ) * (ld*cols));
                }
                for (i =0; i < ld * cols; i++) {
                    matrix->values[i] = 0;
                }
            }
            if ( datatype == MESS_COMPLEX ){
                if ( ld * cols == 0 ) {
                    matrix->values_cpx = NULL;
                } else {
                    mess_try_alloc(matrix->values_cpx, mess_double_cpx_t*, sizeof (mess_double_cpx_t) * (ld*cols));
                }
                for (i =0; i < ld*cols; i++) {
                    matrix->values_cpx[i] = 0;
                }
            }

            break;

        case MESS_COORD:
            mess_try_alloc(matrix->rowptr, mess_int_t*, sizeof(mess_int_t)*nnz+1);
            mess_try_alloc(matrix->colptr, mess_int_t*, sizeof(mess_int_t)*nnz+1);
            if ( datatype == MESS_REAL ){
                mess_try_alloc(matrix->values, double*, sizeof(double)*nnz+1);
                for (i=0; i< nnz;i++) {
                    matrix->values[i]= 0.0;
                    matrix->rowptr[i]= 0;
                    matrix->colptr[i]= 0;
                }
            }
            if ( datatype == MESS_COMPLEX ){
                mess_try_alloc(matrix->values_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*nnz+1);
                for (i=0; i< nnz;i++) {
                    matrix->values_cpx[i]= 0.0;
                    matrix->rowptr[i] =0 ;
                    matrix->colptr[i] =0 ;
                }
            }
            break;
        default:
            MSG_ERROR("storage type not supported\n");
            return MESS_ERROR_STORAGETYPE;
    }

    matrix->rows = rows;
    matrix->cols = cols;
    matrix->nnz = nnz;
    if ( ld ) {
        matrix->ld = ld;
    } else {
        matrix->ld = rows;
    }
    matrix->store_type = storetype;
    matrix->data_type = datatype;
    matrix->symmetry = MESS_GENERAL;

    return 0;
}


/**
 * @brief Wrapper for @ref mess_matrix_realloc.
 * @param[in,out] matrix    input/output matrix to be resized
 * @param[in] rows          input new number of rows
 * @param[in] cols          input new number of cols
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_resize function is only a wrapper for \ref mess_matrix_realloc
 * to get a mathematical correct name for it. The name realloc tends to be too
 * computer science specific. \n
 * If the matrix is a sparse one, the number of non zeros entries is untouched. Otherwise
 * it behaves exactly like the @ref mess_matrix_realloc function and returns the typical error codes.
 *
 */
int mess_matrix_resize ( mess_matrix matrix, mess_int_t rows, mess_int_t cols)
{
    if ( MESS_IS_DENSE(matrix)) {
        return mess_matrix_realloc(matrix, rows, cols, rows*cols);
    } else {
        return mess_matrix_realloc(matrix, rows, cols, matrix->nnz);
    }
}    /* -----  end of function mess_matrix_resize  ----- */

/**
 * @brief Reallocate or resize a \f$ (m \times n) \f$ matrix.
 * @param[in,out] matrix    input/output matrix to be resized
 * @param[in] rows           input new number of rows
 * @param[in] cols          input new number of cols
 * @param[in] nnz          input new number of nonzero elements
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_realloc function reallocates data structures for
 * a matrix. \n
 * In case of a dense matrix the matrix is resized correctly  in sense of cutting of the left upper
 * \f$ (1:rows, 1:cols) \f$ block. \n
 * If the matrix is a compressed sparse matrix only the values array and the
 * corresponding index array is reallocated. Changing the dimension
 * is not possible yet. This might change in the future.\n
 * In case of coordinate matrices it acts like the dense one and ignores
 * the \f$ nnz \f$ argument if the new matrix is smaller than the original one.
 *
 */
int mess_matrix_realloc( mess_matrix matrix,
        mess_int_t rows,
        mess_int_t cols,
        mess_int_t nnz)
{
    MSG_FNAME ( __func__ ) ;
    mess_int_t ld = 0;
    mess_check_nullpointer(matrix);
    if ( rows < 0 || cols < 0 || nnz < 0 ){
        MSG_ERROR("wrong arguments: rows = " MESS_PRINTF_INT ", cols = " MESS_PRINTF_INT ", nnz = " MESS_PRINTF_INT "\n", rows, cols, nnz);
        return MESS_ERROR_ARGUMENTS;
    }
    ld = rows;
    switch (matrix->store_type) {
        case MESS_CSC:
            if ( matrix->rows != rows || matrix->cols != cols) {
                MSG_ERROR("can not resize CSC matrix\n");
                return MESS_ERROR_NOSUPPORT;
            }
            mess_try_realloc(matrix->rowptr, mess_int_t *, sizeof(mess_int_t) * nnz);
            if ( MESS_IS_REAL(matrix)) mess_try_realloc (matrix->values, double *, sizeof(double) * nnz);
            if ( MESS_IS_COMPLEX(matrix)) mess_try_realloc (matrix->values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * nnz);
            break;

        case MESS_CSR:
            if ( matrix->rows != rows || matrix->cols != cols) {
                MSG_ERROR("can not resize CSC matrix\n");
                return MESS_ERROR_NOSUPPORT;
            }
            mess_try_realloc(matrix->colptr, mess_int_t *, sizeof(mess_int_t) * nnz);
            if ( MESS_IS_REAL(matrix)) mess_try_realloc (matrix->values, double *, sizeof(double) * nnz);
            if ( MESS_IS_COMPLEX(matrix)) mess_try_realloc (matrix->values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * nnz);
            break;

        case MESS_DENSE:
            matrix->rowptr = NULL;
            matrix->colptr = NULL;
            ld = matrix->ld;
            nnz = rows*cols;
            if ( matrix->rows == rows) {
                mess_int_t i;
                /*-----------------------------------------------------------------------------
                 *  row size is the same
                 *-----------------------------------------------------------------------------*/
                if ( matrix->data_type == MESS_REAL ) {
                    mess_try_realloc(matrix->values, double*, sizeof(double) * cols*ld);
                    if ( cols > matrix->cols) {
                        for (i = matrix->cols*matrix->ld; i < cols*ld; i++) matrix->values[i] = 0;
                    }
                }
                else if ( matrix->data_type == MESS_COMPLEX ){
                    mess_try_realloc(matrix->values_cpx, mess_double_cpx_t*, sizeof (mess_double_cpx_t) * (cols*ld));
                    if ( cols > matrix->cols) {
                        for (i = matrix->cols*matrix->ld; i < cols*ld; i++) matrix->values_cpx[i] = 0;
                    }
                }
            } else if ( rows < matrix->rows) {
                if (matrix->data_type == MESS_REAL) {
                    ld = __select_ld(rows, sizeof(double));
                } else {
                    ld = __select_ld(rows, sizeof(mess_double_cpx_t));
                }

                /*-----------------------------------------------------------------------------
                 *  rows < matrix->rows
                 *-----------------------------------------------------------------------------*/
                mess_int_t i, j,mn;
                if (MESS_IS_REAL(matrix) ){
                    mn = MESS_MIN(cols, matrix->cols);
                    for (i=1;i<mn;i++){
                        for (j=0; j < rows; j++){
                            matrix->values[i*ld+j] = matrix->values[i*matrix->ld+j];
                        }
                    }
                    mess_try_realloc(matrix->values, double*, sizeof(double)*ld*cols);

                }
                else if (MESS_IS_COMPLEX(matrix) ){
                    mn = MESS_MIN(cols, matrix->cols);
                    for (i=1;i<mn;i++){
                        for (j=0; j < rows; j++){
                            matrix->values_cpx[i*ld+j] = matrix->values_cpx[i*matrix->ld+j];
                        }
                    }
                    mess_try_realloc(matrix->values_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*ld*cols);
                }
            } else {
                if (matrix->data_type == MESS_REAL) {
                    ld = __select_ld(rows, sizeof(double));
                } else {
                    ld = __select_ld(rows, sizeof(mess_double_cpx_t));
                }
                /*-----------------------------------------------------------------------------
                 *  rows > matrix->rows
                 *  we need an in place algorithm here
                 *-----------------------------------------------------------------------------*/
                mess_int_t i, j,mn;
                if (MESS_IS_REAL(matrix) ){
                    double *help = NULL;
                    mess_try_alloc ( help, double *, sizeof(double)*ld*cols);
                    for ( i = 0; i < ld*cols ; i++ ) { help[i] = 0.0; };
                    mn= MESS_MIN(cols, matrix->cols);
                    for (i=0;i<mn;i++){
                        for (j=0; j < matrix->rows; j++){
                            help[i*ld+j] = matrix->values[i*matrix->ld+j];
                        }
                    }
                    mess_free(matrix->values);
                    matrix->values = help;
                }
                else if (MESS_IS_COMPLEX(matrix) ){
                    mess_double_cpx_t *help = NULL;
                    mess_try_alloc ( help, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*ld*cols);
                    for ( i = 0; i < ld*cols ; i++ ) { help[i] = 0.0; };
                    mn = MESS_MIN(cols, matrix->cols);
                    for (i=0;i<mn;i++){
                        for (j=0; j < matrix->rows; j++){
                            help[i*ld+j] = matrix->values_cpx[i*matrix->ld+j];
                        }
                    }
                    mess_free(matrix->values_cpx);
                    matrix->values_cpx = help;
                }
            }

            break;
        case MESS_COORD:{
                            mess_int_t i,k,new_nnz;
                            mess_int_t *rowptr, *colptr;
                            new_nnz = 0;
                            for ( i = 0; i < matrix->nnz; i++ ) {
                                if ( matrix->rowptr[i] < rows && matrix->colptr[i] < cols ) {
                                    new_nnz ++;
                                }
                            }
                            if ( new_nnz ) {
                                nnz = new_nnz;
                            }
                            mess_try_alloc(rowptr, mess_int_t *, sizeof(mess_int_t) * nnz);
                            mess_try_alloc(colptr, mess_int_t *, sizeof(mess_int_t) * nnz);
                            if (MESS_IS_REAL(matrix)){
                                double *values;
                                k = 0;
                                mess_try_alloc(values, double *, sizeof(double)*nnz);
                                for ( i = 0; i < matrix->nnz; i++ ) {
                                    if ( matrix->rowptr[i] < rows && matrix->colptr[i] < cols ) {
                                        rowptr[k] = matrix->rowptr[i];
                                        colptr[k] = matrix->colptr[i];
                                        values[k] = matrix->values[i];
                                        k++;
                                    }
                                }
                                mess_free(matrix->rowptr);
                                mess_free(matrix->colptr);
                                mess_free(matrix->values);
                                matrix->rowptr = rowptr;
                                matrix->colptr = colptr;
                                matrix->values = values;
                            } else {
                                // complex
                                mess_double_cpx_t  *values;
                                k = 0;
                                mess_try_alloc(values, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*nnz);
                                for ( i = 0; i < matrix->nnz; i++ ) {
                                    if ( matrix->rowptr[i] < rows && matrix->colptr[i] < cols ) {
                                        rowptr[k] = matrix->rowptr[i];
                                        colptr[k] = matrix->colptr[i];
                                        values[k] = matrix->values_cpx[i];
                                        k++;
                                    }
                                }
                                mess_free(matrix->rowptr);
                                mess_free(matrix->colptr);
                                mess_free(matrix->values_cpx);
                                matrix->rowptr = rowptr;
                                matrix->colptr = colptr;
                                matrix->values_cpx = values;
                            }
                        }
                        break;
        default:
                        MSG_ERROR("storage type not supported\n");
                        return MESS_ERROR_STORAGETYPE;
    }

    matrix->rows = rows;
    matrix->cols = cols;
    matrix->ld = ld;/* necessary when the matrix row dimension is changed */
    matrix->nnz = nnz;

    return 0;
}

/**
 *
 * @brief Resets a @ref mess_matrix object.
 * @param[in,out] matrix  input/output matrix object to reset
 * @return zeros on success or a non zero error value.
 *
 * The @ref mess_matrix_reset function resets a \ref mess_matrix object that it behaves like a newly initialized one.
 * In contrast to @ref mess_matrix_clear it does only free the internal data and not the surrounding structure.
 */
int mess_matrix_reset ( mess_matrix matrix  )
{
    MSG_FNAME(__func__);

    mess_check_nullpointer(matrix);

    mess_free(matrix->values);
    mess_free(matrix->values_cpx);
    mess_free(matrix->colptr);
    mess_free(matrix->rowptr);
    matrix->data_type = MESS_REAL;
    matrix->nnz = 0;
    matrix->cols = 0;
    matrix->rows = 0;
    matrix->store_type =MESS_UNKNOWN;
    matrix->symmetry = MESS_GENERAL;
    return 0;
}    /* -----  end of function mess_matrix_reset  ----- */





