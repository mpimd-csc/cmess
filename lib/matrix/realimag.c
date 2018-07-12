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
 * @file lib/matrix/realimag.c
 * @brief Handle real and imaginary parts of complex matrices.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Get the real part of a matrix.
 * @param[in] matrix input matrix
 * @param[out] realpart matrix with the real part
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_realpart function copies the real part of a given matrix to
 * the output matrix. \n
 * If the input is real the matrix is copied.
 *
 */
int mess_matrix_realpart ( mess_matrix matrix, mess_matrix realpart )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i,j;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(realpart);

    if ( MESS_IS_REAL(matrix)) {
        ret = mess_matrix_copy(matrix, realpart);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    } else {
        MESS_MATRIX_RESET(realpart);
        ret = mess_matrix_alloc (realpart, matrix->rows, matrix->cols, matrix->nnz, matrix->store_type, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        if ( MESS_IS_DENSE ( matrix )) {
            for ( j = 0; j < matrix->cols; j++ ) {
                for ( i = 0; i< matrix->rows; i++ ) {
                    realpart->values[i+j*realpart->ld] = creal(matrix->values_cpx[i+j*matrix->ld]);
                }
            }
        } else {
            if ( MESS_IS_CSR( matrix )) {
                memcpy(realpart->rowptr, matrix->rowptr, (matrix->rows+1)*sizeof(mess_int_t));
                memcpy(realpart->colptr, matrix->colptr, (matrix->nnz)*sizeof(mess_int_t));
            } else   if ( MESS_IS_CSC( matrix )) {
                memcpy(realpart->rowptr, matrix->rowptr, (matrix->nnz)*sizeof(mess_int_t));
                memcpy(realpart->colptr, matrix->colptr, (matrix->cols+1)*sizeof(mess_int_t));
            } else   if ( MESS_IS_COORD( matrix )) {
                memcpy(realpart->rowptr, matrix->rowptr, (matrix->nnz)*sizeof(mess_int_t));
                memcpy(realpart->colptr, matrix->colptr, (matrix->nnz)*sizeof(mess_int_t));
            } else {
                MSG_ERROR("unknown/unsupported storage type\n");
                return MESS_ERROR_STORAGETYPE;
            }
            for ( i =0 ; i < matrix->nnz; i++) {
                realpart->values[i]=creal(matrix->values_cpx[i]);
            }
            ret = mess_matrix_eliminate_zeros( realpart );        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eliminate_zeros);
        }
        //ret = mess_matrix_copy(matrix, realpart);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_copy);
        //ret = mess_matrix_toreal(realpart);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_toreal);;
    }
    return 0 ;
}       /* -----  end of function mess_matrix_realpart  ----- */

/**
 * @brief Get the imaginary part of a matrix.
 * @param[in] matrix input matrix
 * @param[out] imagpart  matrix with the imaginary part
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_imagpart function copies the imaginary part of a given matrix to
 * the output matrix. \n
 * If the input is real the output is a zero matrix.
 *
 */
int mess_matrix_imagpart ( mess_matrix matrix, mess_matrix imagpart )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i,j;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(imagpart);

    MESS_MATRIX_RESET(imagpart);
    ret = mess_matrix_alloc (imagpart, matrix->rows, matrix->cols, matrix->nnz, matrix->store_type, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    if ( MESS_IS_REAL(matrix)) {
        ret = mess_matrix_zeros(imagpart);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);
    } else {
        if ( MESS_IS_DENSE ( matrix )) {
            for ( j = 0; j < matrix->cols; j++ ) {
                for ( i = 0; i< matrix->rows; i++ ) {
                    imagpart->values[i+j*imagpart->ld] = cimag(matrix->values_cpx[i+j*matrix->ld]);
                }
            }
        } else {
            if ( MESS_IS_CSR( matrix )) {
                memcpy(imagpart->rowptr, matrix->rowptr, (matrix->rows+1)*sizeof(mess_int_t));
                memcpy(imagpart->colptr, matrix->colptr, (matrix->nnz)*sizeof(mess_int_t));
            } else   if ( MESS_IS_CSC( matrix )) {
                memcpy(imagpart->rowptr, matrix->rowptr, (matrix->nnz)*sizeof(mess_int_t));
                memcpy(imagpart->colptr, matrix->colptr, (matrix->cols+1)*sizeof(mess_int_t));
            } else   if ( MESS_IS_COORD( matrix )) {
                memcpy(imagpart->rowptr, matrix->rowptr, (matrix->nnz)*sizeof(mess_int_t));
                memcpy(imagpart->colptr, matrix->colptr, (matrix->nnz)*sizeof(mess_int_t));
            } else {
                MSG_ERROR("unknown/unsupported storage type\n");
                return MESS_ERROR_STORAGETYPE;
            }
            for ( i =0 ; i < matrix->nnz; i++) {
                imagpart->values[i]=cimag(matrix->values_cpx[i]);
            }
            mess_matrix_eliminate_zeros( imagpart );
        }
    }
    return 0 ;
}       /* -----  end of function mess_matrix_imagpart  ----- */


/**
 * @brief Conjugate a matrix.
 * @param[in,out] A input/output matrix \f$A\f$ to conjugate
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_conj function congujates a complex matrix, i.e. it computes
 * \f[ A\leftarrow\bar{A} . \f]
 *
 */
int mess_matrix_conj ( mess_matrix A )
{
    MSG_FNAME(__func__);
    mess_int_t i,j ;

    mess_check_nullpointer(A);
    if (!MESS_IS_COMPLEX(A)) return 0;
    if ( MESS_IS_DENSE(A)) {
        for ( j =0 ; j < A->cols; j++) {
            for ( i = 0;i < A->rows ; i++ ) {
                A->values_cpx[i+j*A->ld] = conj(A->values_cpx[i+j*A->ld]);
            }
        }
    } else {
        for (i = 0; i < A->nnz; i++){
            A->values_cpx[i] = conj(A->values_cpx[i]);
        }
    }
    return 0;
}       /* -----  end of function mess_matrix_conj  ----- */



/**
 * @brief Create a complex matrix from real and imaginary parts.
 * @param[in] xr    input real part \f$xr\f$
 * @param[in] xc    input complex part \f$xc\f$
 * @param[out] x    output complex matrix \f$x\f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_complex_from_parts function creates a complex matrix from
 * separate real and imaginary parts. \n
 * Mathematically this means
 * \f[ x = xr +xc \cdot i. \f]
 * @attention If the matrices \f$ xr \f$ and \f$ xc \f$ are sparse they need to have the same pattern structure otherwise
 * the result will be completely wrong.
 *
 * @see mess_matrix_addc
 * @see mess_matrix_scale
 */
int  mess_matrix_complex_from_parts ( mess_matrix xr, mess_matrix xc, mess_matrix x  )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i,j;

    /*-----------------------------------------------------------------------------
     *  check inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(xr);
    mess_check_nullpointer(xc);
    mess_check_real(xc);
    mess_check_real(xr);
    mess_check_same_size(xc,xr);

    if ( xc->store_type != xr->store_type){
        MSG_ERROR("Both matrices must have the same storage type\n");
        return MESS_ERROR_STORAGETYPE;
    }
    if (mess_matrix_need_alloc(x, xr->rows, xr->cols, xr->nnz, xr->store_type, MESS_COMPLEX)) {
        ret = mess_matrix_alloc(x, xr->rows, xr->cols, xr->nnz, xr->store_type, MESS_COMPLEX);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    }

    if ( MESS_IS_DENSE(xr)) {
        for ( j = 0 ;  j < xr->cols; j++) {
            for ( i = 0; i< xr->rows ; i++) {
                x->values_cpx[i+j*x->ld] = xr->values[i+j*xr->ld] + xc->values[i+j*xc->ld]*I;
            }
        }
    } else {
        ret = mess_matrix_copy(xr,x);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_addc(1.0*I,xc,1.0,x);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_addc);
    }
    return 0;
}       /* -----  end of function mess_matrix_complex_from_parts  ----- */


