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
 * @file lib/matrix/cimport.c
 * @brief Import matrices from @c C and @c Fortran arrays.
 * @author @koehlerm
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>

/**
 * @brief Create a  @ref MESS_DENSE @ref mess_matrix from a @c Fortran array.
 * @param[out] mat      output matrix to copy the data to
 * @param[in] rows      input number of rows of the matrix
 * @param[in] cols      input number of columns of the matrix
 * @param[in] ld        input leading dimension of the matrix
 * @param[in] realv     input array with real values
 * @param[in] complexv  input array with complex values
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_dense_from_farray function copies a @c Fortran-like (column-major)
 * matrix to a @ref mess_matrix in @ref MESS_DENSE format. \n
 * If the provided leading dimension is smaller than the number of rows, the number of rows is used.
 * Depending whether @p realv  or @p complexv  is given a real matrix or a complex matrix is created.
 * If both are given  the complex data is used. \n
 * The leading dimension of the @ref mess_matrix is determined by \ref mess_matrix_alloc.
 *
 * @sa mess_matrix_csr
 * @sa mess_matrix_coord
 * @sa mess_matrix_dense_from_carray
 *
 */
int mess_matrix_dense_from_farray ( mess_matrix mat, mess_int_t rows, mess_int_t cols, mess_int_t ld,
        double *realv, mess_double_cpx_t *complexv )
{
    MSG_FNAME(__func__);
    int ret;
    mess_datatype_t dt;
    /*-----------------------------------------------------------------------------
     * Check Inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(mat);
    mess_check_positive(rows);
    mess_check_positive(cols);

    if ( !realv && !complexv) {
        MSG_ERROR("realv xor complexv must be given");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( realv ) dt = MESS_REAL;
    if ( complexv ) dt = MESS_COMPLEX;
    if ( ld < rows ) ld = rows;


    ret = mess_matrix_alloc(mat, rows, cols, rows*cols, MESS_DENSE, dt);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    /*-----------------------------------------------------------------------------
     *  real matrix
     *-----------------------------------------------------------------------------*/
    if ( dt == MESS_REAL) {
        F77_GLOBAL(dlacpy,DLACPY)("A",&rows,&cols,realv, &ld, mat->values, &mat->ld);
    }
    /*-----------------------------------------------------------------------------
     *  complex matrix
     *-----------------------------------------------------------------------------*/
    else {
        F77_GLOBAL(zlacpy,ZLACPY)("A",&rows,&cols,complexv, &ld, mat->values_cpx, &mat->ld);
    }
    return 0;
}       /* -----  end of function mess_matrix_dense_from_farray  ----- */

/**
 * @brief Create a @ref MESS_DENSE @ref mess_matrix from a @c C \f$ 2 \f$-dimensional array.
 * @param[out] mat          output matrix to copy the data to
 * @param[in] rows          input number of rows of the matrix
 * @param[in] cols          input number of columns of the matrix
 * @param[in] realv         input \f$ 2 \f$ D-array with real values
 * @param[in] complexv      input \f$ 2 \f$ D-array with complex values
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_dense_from_carray function copies a @c C \f$ 2 \f$-dimensional array(row-major,
 * double** pointer) matrix to a @ref mess_matrix in @ref MESS_DENSE format.\n
 * Depending whether @p realv  or @p complexv  is given a real matrix or a complex matrix is created.
 * If both are given the complex data is used. \n
 * The leading dimension of the @ref mess_matrix is determined by \ref mess_matrix_alloc.
 *
 * @sa mess_matrix_csc
 * @sa mess_matrix_csr
 * @sa mess_matrix_coord
 * @sa mess_matrix_dense_from_farray
 */
int mess_matrix_dense_from_carray ( mess_matrix mat, mess_int_t rows, mess_int_t cols,
        double **realv, mess_double_cpx_t **complexv )
{
    MSG_FNAME(__func__);
    int ret;
    mess_datatype_t dt;
    mess_int_t i,j;

    /*-----------------------------------------------------------------------------
     * Check Inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(mat);
    mess_check_positive(rows);
    mess_check_positive(cols);

    if ( !realv && !complexv) {
        MSG_ERROR("realv xor complexv must be given");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( realv ) dt = MESS_REAL;
    if ( complexv ) dt = MESS_COMPLEX;

    ret = mess_matrix_alloc(mat, rows, cols, rows*cols, MESS_DENSE, dt);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    /*-----------------------------------------------------------------------------
     *  real matrix
     *-----------------------------------------------------------------------------*/
    if ( dt == MESS_REAL) {
        // memcpy ( mat->values, realv, rows*cols*sizeof(double));
        for ( j = 0 ; j<cols; j++ ) {
            for ( i = 0; i<rows; i++ ) {
                mat->values[i+j*mat->ld] = realv[i][j];
            }
        }
    }
    /*-----------------------------------------------------------------------------
     *  complex matrix
     *-----------------------------------------------------------------------------*/
    else {
        // memcpy( mat->values_cpx, complexv, rows*cols*sizeof(mess_double_cpx_t));
        for ( j = 0 ; j<cols; j++ ) {
            for ( i = 0; i<rows; i++ ) {
                mat->values_cpx[i+j*mat->ld] = complexv[i][j];
            }
        }
    }
    return 0;
}       /* -----  end of function mess_matrix_dense_from_farray  ----- */



/**
 * @brief Create a @ref MESS_CSR @ref mess_matrix from given arrays.
 * @param[out]  matrix      output matrix which should be created from the given data
 * @param[in]   rows        input number of rows
 * @param[in]   cols        input number of columns
 * @param[in]   rowptr      input array containing row pointers
 * @param[in]   colptr      input array containing column pointers
 * @param[in]   values      input array containing real values
 * @param[in]   values_cpx  input array containing complex values
 * @return zero on success or a non zero error code
 *
 * The @ref mess_matrix_csr function creates a @ref MESS_CSR @ref mess_matrix out of the given arrays. \n
 * Depending on the first element in the row pointer array @p rowptr  it detects if zero or one based indexing
 * is used. \n
 * The last entry of the row pointer array @p rowptr[rows]  contains the number of non zero elements of the matrix.
 * If the @p values  array is @c NULL, the @p values_cpx  array is used and the generated matrix is complex.
 * If both value pointer are @c NULL an error is returned. If both are not @c NULL the complex data is used.
 *
 * @sa mess_matrix_csc
 * @sa mess_matrix_coord
 * @sa mess_matrix_dense_from_farray
 * @sa mess_matrix_dense_from_carray
 *
 */
int mess_matrix_csr ( mess_matrix matrix, mess_int_t rows, mess_int_t cols, mess_int_t * rowptr, mess_int_t * colptr, double * values, mess_double_cpx_t *values_cpx )
{
    MSG_FNAME(__func__);
    mess_int_t nnz;
    mess_int_t offset = 0;
    mess_int_t i;
    int ret = 0;
    int cpx = 0;


    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(rowptr);
    mess_check_nullpointer(colptr);
    mess_check_nonnegative(rows);
    mess_check_nonnegative(cols);
    if ( values == NULL && values_cpx == NULL ) {
        MSG_ERROR("One of values or values_cpx must be given.\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( values_cpx != NULL) {
        cpx = 1;
    }

    nnz = rowptr[rows];
    if (rowptr[0]>0) {
        offset = 1;
        MSG_INFO("Use one-based indexing.\n");
    }

    ret = mess_matrix_alloc(matrix, rows, cols, nnz, MESS_CSR, (cpx)?MESS_COMPLEX: MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    for (i = 0; i < rows; i++) {
        matrix->rowptr[i] = rowptr[i]-offset;
    }
    matrix->rowptr[rows] = nnz;

    if ( cpx ) {
        for (i = 0; i < nnz; i++) {
            matrix->colptr[i]=colptr[i]-offset;
            matrix->values_cpx[i]= values_cpx[i];
        }
    } else {
        for (i = 0; i < nnz; i++) {
            matrix->colptr[i]=colptr[i]-offset;
            matrix->values[i]= values[i];
        }
    }
    return 0;
}       /* -----  end of function mess_matrix_csr  ----- */

/**
 * @brief Create a @ref MESS_CSC @ref mess_matrix from given arrays.
 * @param[out]  matrix      output matrix which should be created from the given data
 * @param[in]   rows        input number of rows
 * @param[in]   cols        input number of columns
 * @param[in]   rowptr      input array containing row pointers
 * @param[in]   colptr      input array containing column pointers
 * @param[in]   values      input array containing real values
 * @param[in]   values_cpx  input array containing complex values
 * @return zero on success or a non zero error code
 *
 * The @ref mess_matrix_csc function creates a @ref MESS_CSC @ref matrix out of the given arrays.\n
 * Depending on the first element in the column pointer array @p colptr  it detects if zero or one based indexing
 * is used. \n
 * The last entry of the column pointer array @p colptr[cols]  contains the number of non zero elements of the
 * matrix. If the @p values array is @c NULL, the @p values_cpx array is used and the generated matrix is complex.
 * If both value pointer are @c NULL an error is returned. If both are not @c NULL the complex data is used.
 *
 * @sa mess_matrix_csr
 * @sa mess_matrix_coord
 * @sa mess_matrix_dense_from_farray
 * @sa mess_matrix_dense_from_carray
 *
 */
int mess_matrix_csc ( mess_matrix matrix, mess_int_t rows, mess_int_t cols, mess_int_t * rowptr, mess_int_t * colptr, double * values, mess_double_cpx_t *values_cpx )
{
    MSG_FNAME(__func__);
    mess_int_t nnz;
    mess_int_t offset = 0;
    mess_int_t i;
    int ret = 0;
    int cpx = 0;


    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(rowptr);
    mess_check_nullpointer(colptr);
    mess_check_nonnegative(rows);
    mess_check_nonnegative(cols);
    if ( values == NULL && values_cpx == NULL ) {
        MSG_ERROR("One of values or values_cpx must be given.\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( values_cpx != NULL) {
        cpx = 1;
    }

    nnz = colptr[cols];
    if (colptr[0]>0) {
        offset = 1;
        MSG_INFO("Use one-based indexing.\n");
    }

    ret = mess_matrix_alloc(matrix, rows, cols, nnz, MESS_CSC, (cpx)?MESS_COMPLEX: MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    for (i = 0; i < cols; i++) {
        matrix->colptr[i] = colptr[i]-offset;
    }
    matrix->colptr[cols] = nnz;

    if ( cpx ) {
        for (i = 0; i < nnz; i++) {
            matrix->rowptr[i]=rowptr[i]-offset;
            matrix->values_cpx[i]= values_cpx[i];
        }
    } else {
        for (i = 0; i < nnz; i++) {
            matrix->rowptr[i]=rowptr[i]-offset;
            matrix->values[i]= values[i];
        }
    }
    return 0;
}       /* -----  end of function mess_matrix_csc  ----- */


/**
 * @brief Create a @ref MESS_COORD @ref mess_matrix from given arrays.
 * @param[out]  matrix      output matrix which should be created from the given data
 * @param[in]   rows        input number of rows
 * @param[in]   cols        input number of columns
 * @param[in]   nnz         input number of non zero elements
 * @param[in]   onebased    input flag to identify one based indexing
 * @param[in]   rowptr      input array containing row pointers
 * @param[in]   colptr      input array containing column pointers
 * @param[in]   values      input array containing real values
 * @param[in]   values_cpx  input array containing complex values
 * @return zero on success or a non zero error code
 *
 * The @ref mess_matrix_coord function creates a @ref MESS_COORD @ref mess_matrix out of the given arrays. \n
 * If the @p values  array is @c NULL, the @p values_cpx  array is used and the generated matrix is complex.
 * If both value pointer are @c NULL an error is returned. If both are not @c NULL the complex data is used.
 *
 * @sa mess_matrix_csr
 * @sa mess_matrix_csc
 * @sa mess_matrix_dense_from_farray
 * @sa mess_matrix_dense_from_carray
 *
 */
int mess_matrix_coord ( mess_matrix matrix, mess_int_t rows, mess_int_t cols, mess_int_t nnz, mess_int_t onebased, mess_int_t * rowptr, mess_int_t * colptr, double * values, mess_double_cpx_t *values_cpx )
{
    MSG_FNAME(__func__);
    mess_int_t offset = 0;
    mess_int_t i;
    int ret = 0;
    int cpx = 0;


    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(rowptr);
    mess_check_nullpointer(colptr);
    mess_check_nonnegative(rows);
    mess_check_nonnegative(cols);
    if ( values == NULL && values_cpx == NULL ) {
        MSG_ERROR("One of values or values_cpx must be given.\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( values_cpx != NULL) {
        cpx = 1;
    }

    if (onebased){
        offset = 1;
        MSG_INFO("Use one-based indexing.\n");
    }

    ret = mess_matrix_alloc(matrix, rows, cols, nnz, MESS_COORD, (cpx)?MESS_COMPLEX: MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    for (i = 0; i < nnz; i++) {
        matrix->rowptr[i] = rowptr[i]-offset;
        matrix->colptr[i] = colptr[i]-offset;
    }
    if ( cpx ) {
        for (i = 0; i < nnz; i++) {
            matrix->values_cpx[i]= values_cpx[i];
        }
    } else {
        for (i = 0; i < nnz; i++) {
            matrix->values[i]= values[i];
        }
    }
    return 0;
}       /* -----  end of function mess_matrix_csr  ----- */

