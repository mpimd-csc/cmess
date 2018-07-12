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
 * @file lib/matrix/scal.c
 * @brief Scale a matrix.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/eigenvalue.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>


/**
 * @brief Scale a real matrix by a scalar.
 * @param[in] alpha     input scaling factor \f$\alpha\f$
 * @param[in,out] A     matrix to scale \f$A\f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_scale function scales a real matrix by \f$ alpha \f$:
 * \f[ A \leftarrow \alpha A.\f]
 *
 * @sa mess_matrix_scalec
 */
int mess_matrix_scale(double alpha, mess_matrix A){
    MSG_FNAME(__func__);
    mess_int_t one = 1;
    mess_int_t N;
    mess_check_nullpointer(A);
    mess_check_real_or_complex(A);

    if ( MESS_IS_COMPLEX(A)) return(mess_matrix_scalec(alpha, A));
    mess_check_real(A);

    if (MESS_IS_DENSE(A)){
        int i ;
        N = A->rows;
        for ( i = 0; i < A->cols; i++){
            F77_GLOBAL(dscal,DSCAL)(&N,&alpha, &A->values[i*A->ld],&one);
        }
    } else {
        N = A->nnz;
        F77_GLOBAL(dscal,DSCAL)(&N, &alpha, A->values, &one);
    }

    return (0);
}

/**
 * @brief Scale a complex matrix by a scalar.
 * @param[in] alpha input scaling factor
 * @param[in,out] A     matrix to scale
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_scalec function scales a complex matrix by \f$ alpha \f$:
 * \f[ A \leftarrow \alpha A. \f]
 *
 * @sa mess_matrix_scale
 */
int mess_matrix_scalec(mess_double_cpx_t alpha, mess_matrix A){
    MSG_FNAME(__func__);
    mess_int_t one = 1;
    mess_int_t N;
    mess_check_nullpointer(A);
    int ret = 0;


    ret = mess_matrix_tocomplex(A);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
    if (MESS_IS_DENSE(A)){
        int i ;
        N = A->rows;
        for ( i = 0; i < A->cols; i++){
            F77_GLOBAL(zscal,ZSCAL)(&N,&alpha, &A->values_cpx[i*A->ld],&one);
        }
    } else {
        N = A->nnz;
        F77_GLOBAL(zscal,ZSCAL)(&N, &alpha, A->values_cpx, &one);
    }

    return(0) ;
}


/**
 * @brief Scale rows of a matrix.
 * @param[in] r         input @ref mess_vector with scaling vectors
 * @param[in,out] A     matrix to scale
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_rowscalem function scales \f$ i-th\f$ row  of a matrix by  the \f$  i-th \f$
 * entry of @p r.
 * @ref mess_matrix.rows should be equal to @ref mess_vector.dim.
 *
 * @sa mess_matrix_scale
 * @sa mess_matrix_colscalem
 */
int mess_matrix_rowscalem(mess_vector r, mess_matrix A){
    MSG_FNAME(__func__);
    mess_int_t ret =0, i = 0, j =0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(r);
    mess_check_real_or_complex(A);
    mess_check_real_or_complex(r);
    mess_check_same_rowsdim(A,r);

    /*-----------------------------------------------------------------------------
     *  perform operation
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(A) && MESS_IS_COMPLEX(r)){
        ret = mess_matrix_tocomplex(A);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
    }

    if ( MESS_IS_DENSE (A) ){
        if(MESS_IS_REAL(A) && MESS_IS_REAL(r)){
            for ( i = 0; i < A->rows; i++){
                for( j = 0; j< A->cols; j++){
                    A->values[i + j*A->ld] *= r->values[i];
                }
            }
        }else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(r)){
            for ( i = 0; i < A->rows; i++){
                for( j = 0; j< A->cols; j++){
                    A->values_cpx[i + j*A->ld] *= r->values[i];
                }
            }
        }else{
            for ( i = 0; i < A->rows; i++){
                for( j = 0; j< A->cols; j++){
                    A->values_cpx[i + j*A->ld] *= r->values_cpx[i];
                }
            }
        }
    } else if ( MESS_IS_CSC(A) ) {
        if(MESS_IS_REAL(A) && MESS_IS_REAL(r)){
            for ( i = 0 ; i < A->cols; i++ ) {
                for ( j = A->colptr[i]; j < A->colptr[i+1]; j++ ) {
                    A->values[j] *= r->values[A->rowptr[j]];
                }
            }

        }else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(r)){
            for ( i = 0 ; i < A->cols; i++ ) {
                for ( j = A->colptr[i]; j < A->colptr[i+1]; j++ ) {
                    A->values_cpx[j] *= r->values[A->rowptr[j]];
                }
            }

        }else{
            for ( i = 0 ; i < A->cols; i++ ) {
                for ( j = A->colptr[i]; j < A->colptr[i+1]; j++ ) {
                    A->values_cpx[j] *= r->values_cpx[A->rowptr[j]];
                }
            }

        }
    }  else if ( MESS_IS_CSR(A)) {
        if(MESS_IS_REAL(A) && MESS_IS_REAL(r)){
            for ( i = 0; i < A->rows; ++i){
                for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                    A->values[j] *=  r->values[i];
                }
            }

        }else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(r)){
            for ( i = 0; i < A->rows; ++i){
                for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                    A->values_cpx[j] *=  r->values[i];
                }
            }

        }else{
            for ( i = 0; i < A->rows; ++i){
                for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                    A->values_cpx[j] *=  r->values_cpx[i];
                }
            }
        }
    } else if (MESS_IS_COORD(A)){
        if(MESS_IS_REAL(A) && MESS_IS_REAL(r)){
            for(i=0;i<A->nnz;++i){
                A->values[i] *= r->values[A->rowptr[i]];
            }
        }else if ( MESS_IS_COMPLEX(A) && MESS_IS_REAL(r)){
            for(i=0;i<A->nnz;++i){
                A->values_cpx[i] *= r->values[A->rowptr[i]];
            }
        }else{
            for(i=0;i<A->nnz;++i){
                A->values_cpx[i] *= r->values_cpx[A->rowptr[i]];
            }
        }
    }else {
        MSG_ERROR("Unsupported Storage type: %s \n", mess_storage_t_str(A->store_type));
        return(MESS_ERROR_STORAGETYPE);
    }

    return(0);
}


/**
 * @brief Scale cols of a matrix.
 * @param[in] c         input @ref mess_vector with scaling vectors
 * @param[in,out] A     matrix to scale
 * @return zero on success or a non zero error code otherwise
 * @author @mbehr
 *
 * The @ref mess_matrix_colscalem function scales \f$i-th\f$ col of a matrix by  the \f$  i-th \f$
 * entry of @p c.
 * @ref mess_matrix.cols should be equal to @ref mess_vector.dim.
 *
 * @sa mess_matrix_scale
 * @sa mess_matrix_rowscalem
 */
int mess_matrix_colscalem(mess_vector c, mess_matrix A){
    MSG_FNAME(__func__);
    mess_int_t ret =0, i = 0, j =0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(c);
    mess_check_real_or_complex(A);
    mess_check_real_or_complex(c);
    mess_check_same_colsdim(A,c);

    /*-----------------------------------------------------------------------------
     *  perform operation
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(A) && MESS_IS_COMPLEX(c)){
        ret = mess_matrix_tocomplex(A);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
    }

    if ( MESS_IS_DENSE (A) ){
        if(MESS_IS_REAL(A) && MESS_IS_REAL(c)){
            for ( i = 0; i < A->rows; i++){
                for( j = 0; j< A->cols; j++){
                    A->values[i + j*A->ld] *= c->values[j];
                }
            }
        }else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(c)){
            for ( i = 0; i < A->rows; i++){
                for( j = 0; j< A->cols; j++){
                    A->values_cpx[i + j*A->ld] *= c->values[j];
                }
            }
        }else{
            for ( i = 0; i < A->rows; i++){
                for( j = 0; j< A->cols; j++){
                    A->values_cpx[i + j*A->ld] *= c->values_cpx[j];
                }
            }
        }
    } else if ( MESS_IS_CSC(A) ) {
        if(MESS_IS_REAL(A) && MESS_IS_REAL(c)){
            for ( i = 0 ; i < A->cols; i++ ) {
                for ( j = A->colptr[i]; j < A->colptr[i+1]; j++ ) {
                    A->values[j] *= c->values[i];
                }
            }

        }else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(c)){
            for ( i = 0 ; i < A->cols; i++ ) {
                for ( j = A->colptr[i]; j < A->colptr[i+1]; j++ ) {
                    A->values_cpx[j] *= c->values[i];
                }
            }

        }else{
            for ( i = 0 ; i < A->cols; i++ ) {
                for ( j = A->colptr[i]; j < A->colptr[i+1]; j++ ) {
                    A->values_cpx[j] *= c->values_cpx[i];
                }
            }
        }
    }  else if ( MESS_IS_CSR(A)) {
        if(MESS_IS_REAL(A) && MESS_IS_REAL(c)){
            for ( i = 0; i < A->rows; ++i){
                for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                    A->values[j] *=  c->values[A->colptr[j]];
                }
            }
        }else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(c)){
            for ( i = 0; i < A->rows; ++i){
                for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                    A->values_cpx[j] *=  c->values[A->colptr[j]];
                }
            }
        }else{
            for ( i = 0; i < A->rows; ++i){
                for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                    A->values_cpx[j] *=  c->values_cpx[A->colptr[j]];
                }
            }
        }
    } else if (MESS_IS_COORD(A)){
        if(MESS_IS_REAL(A) && MESS_IS_REAL(c)){
            for(i=0;i<A->nnz;++i){
                A->values[i] *= c->values[A->colptr[i]];
            }
        }else if ( MESS_IS_COMPLEX(A) && MESS_IS_REAL(c)){
            for(i=0;i<A->nnz;++i){
                A->values_cpx[i] *= c->values[A->colptr[i]];
            }
        }else{
            for(i=0;i<A->nnz;++i){
                A->values_cpx[i] *= c->values_cpx[A->colptr[i]];
            }
        }

    } else {
        MSG_ERROR("Unsupported Storage type: %s \n", mess_storage_t_str(A->store_type));
        return(MESS_ERROR_STORAGETYPE);
    }

    return(0);
}


