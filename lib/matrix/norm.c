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
 * @file lib/matrix/norm.c
 * @brief Compute various norms of a matrix.
 * @author @koehlerm
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>



/**
 * @brief Compute the \f$ 1\f$-norm of a matrix.
 * @param[in] A input matrix \f$A\f$
 * @param[out] nrm output \f$1\f$-norm
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_norm1 function calculates the \f$ 1 \f$-norm of a matrix:
 * \f[nrm=\Vert A \Vert_1.\f]
 * It is also know as the column-sum-norm.
 *
 * @sa mess_matrix_norminf
 * @sa mess_matrix_norm2
 * @sa mess_matrix_normf
 *
 */
int mess_matrix_norm1(mess_matrix A , double *nrm) {
    MSG_FNAME(__func__);
    double max = 0;
    mess_int_t i=0, j=0;
    int ret = 0;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(nrm);
    mess_check_real_or_complex(A);


    switch (A->store_type) {
        /*-----------------------------------------------------------------------------
         *  CSR Storage
         *-----------------------------------------------------------------------------*/
        case MESS_CSR: {
                           double *tmp;
                           mess_try_alloc(tmp, double *, sizeof(double) * A->cols);
                           for ( i = 0; i < A->cols; i++) tmp[i] = 0.0;
                           if ( MESS_IS_REAL(A)) {
                               for (i = 0; i < A->rows; i++){
                                   for (j=A->rowptr[i] ; j < A->rowptr[i+1]; j++){
                                       tmp[A->colptr[j]] += fabs(A->values[j]);
                                       if ( (MESS_IS_SYMMETRIC (A)|| MESS_IS_SKEWSYMMETRIC(A))&& A->colptr[j] != i) {
                                           tmp[i] += fabs(A->values[j]);
                                       }
                                   }
                               }
                           }else if ( MESS_IS_COMPLEX(A)){
                               for (i = 0; i < A->rows; i++){
                                   for (j=A->rowptr[i] ; j < A->rowptr[i+1]; j++){
                                       tmp[A->colptr[j]] += cabs(A->values_cpx[j]);
                                       if ( (MESS_IS_HERMITIAN (A))&& A->colptr[j] != i) {
                                           tmp[i] += cabs(A->values_cpx[j]);
                                       }
                                   }
                               }

                           }
                           max=tmp[0];
                           for (i = 1; i < A->cols; i++){
                               if (tmp[i] > max ) max= tmp[i];
                           }
                           mess_free(tmp);
                           break;
                       }

                       /*-----------------------------------------------------------------------------
                        *  CSC Storage
                        *-----------------------------------------------------------------------------*/
        case MESS_CSC: {
                           double *tmp;
                           mess_try_alloc(tmp, double *, sizeof(double) * A->cols);
                           for ( i = 0; i < A->cols; i++) tmp[i] = 0.0;
                           if ( MESS_IS_REAL(A)){
                               for ( j = 0 ; j < A->cols; j++){
                                   for ( i = A->colptr[j]; i < A->colptr[j+1]; i++){
                                       tmp[j] += fabs(A->values[i]);
                                       if ( (MESS_IS_SYMMETRIC (A)|| MESS_IS_SKEWSYMMETRIC(A))&& A->rowptr[i] != j) {
                                           tmp[A->rowptr[i]] += fabs(A->values[i]);
                                       }
                                   }
                               }
                           } else if ( MESS_IS_COMPLEX(A)){
                               for ( j = 0 ; j < A->cols; j++){

                                   for ( i = A->colptr[j]; i < A->colptr[j+1]; i++){
                                       tmp[j] += cabs(A->values_cpx[i]);
                                       if ( (MESS_IS_HERMITIAN (A))&& A->rowptr[i] != j) {
                                           tmp[A->rowptr[i]] += cabs(A->values_cpx[i]);
                                       }
                                   }
                               }
                           }
                           max=tmp[0];
                           for (i = 1; i < A->cols; i++){
                               if (tmp[i] > max ) max= tmp[i];
                           }
                           mess_free(tmp);
                           break;
                       }

                       /*-----------------------------------------------------------------------------
                        *  coordinate, using converting to CSR
                        *-----------------------------------------------------------------------------*/
        case MESS_COORD:{
                            int r;
                            mess_matrix tmp;
                            MESS_MATRIX_CHECKFORMAT(A, tmp, r, MESS_CSR);
                            if ( r == 0 ){
                                ret =  mess_matrix_norm1(tmp, &max);
                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm1);
                                mess_matrix_clear(&tmp);
                            }else{
                                MSG_ERROR("error convert COORD to CSR\n");
                                ret=MESS_ERROR_CONVERT;
                                mess_matrix_clear(&tmp);
                            }

                            break;
                        }

                        /*-----------------------------------------------------------------------------
                         *  DENSE
                         *-----------------------------------------------------------------------------*/
        case MESS_DENSE: {
                             if(MESS_IS_REAL(A)){
                                 max = F77_GLOBAL(dlange,DLANGE)("1",&(A->rows),&(A->cols),A->values,&(A->ld),NULL);
                             }else{
                                 max = F77_GLOBAL(zlange,ZLANGE)("1",&(A->rows),&(A->cols),A->values_cpx,&(A->ld),NULL);
                             }
                             break;
                         }

        default:
                         MSG_ERROR("unkown/unsupported storage type: %s\n", mess_storage_t_str(A->store_type));
                         return(MESS_ERROR_STORAGETYPE);
    }
    *nrm = max;
    return(ret);
}

/**
 * @brief Compute the \f$ \infty \f$-norm of a matrix.
 * @param[in] A input matrix \f$A\f$
 * @param[out] nrm output \f$ \infty\f$-norm
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_norminf function calculates the \f$ \infty \f$-norm of a matrix:
 * \f[nrm=\Vert A \Vert_\infty . \f]
 * It is also known as row-sum-norm.
 *
 * @sa mess_matrix_norm1
 * @sa mess_matrix_norm2
 * @sa mess_matrix_normf
 */
int mess_matrix_norminf(mess_matrix A, double *nrm){
    MSG_FNAME(__func__);
    double max=0;
    mess_int_t i, j;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(nrm);
    mess_check_real_or_complex(A);

    switch (A->store_type) {
        /*-----------------------------------------------------------------------------
         * CSR Storage
         *-----------------------------------------------------------------------------*/
        case MESS_CSR: {
                           double *tmp;
                           mess_try_alloc(tmp,double *, sizeof(double) * A->rows);
                           for ( i = 0; i < A->rows; i++) tmp[i] = 0.0;
                           if (MESS_IS_REAL(A)){
                               for ( i = 0; i < A->rows; i++){
                                   for (j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                                       tmp[i] += fabs(A->values[j]);
                                       if ( (MESS_IS_SYMMETRIC (A)|| MESS_IS_SKEWSYMMETRIC(A))&& A->colptr[j] != i) {
                                           tmp[A->colptr[j]] += fabs(A->values[j]);
                                       }
                                   }
                               }
                           } else if (MESS_IS_COMPLEX(A)){
                               for ( i = 0; i < A->rows; i++){
                                   for (j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                                       tmp[i] += cabs(A->values_cpx[j]);
                                       if ( (MESS_IS_HERMITIAN(A))&& A->colptr[j] != i) {
                                           tmp[A->colptr[j]] += cabs(A->values_cpx[j]);
                                       }
                                   }
                               }
                           }
                           max=tmp[0];
                           for (i = 1; i < A->rows; i++){
                               if (tmp[i] > max ) max= tmp[i];
                           }
                           mess_free(tmp);
                           break;
                       }

                       /*-----------------------------------------------------------------------------
                        *  CSC Storage
                        *-----------------------------------------------------------------------------*/
        case MESS_CSC: {
                           double *tmp;
                           mess_try_alloc(tmp,double *, sizeof(double) * A->rows);
                           for ( i = 0; i < A->rows; i++) tmp[i] = 0.0;
                           if ( MESS_IS_REAL(A)) {
                               for (i = 0; i < A->cols; i++){
                                   for (j=A->colptr[i] ; j < A->colptr[i+1]; j++){
                                       tmp[A->rowptr[j]] += fabs(A->values[j]);
                                       if ( (MESS_IS_SYMMETRIC (A)|| MESS_IS_SKEWSYMMETRIC(A))&& A->rowptr[j] != i) {
                                           tmp[i] += fabs(A->values[j]);
                                       }
                                   }
                               }
                           } else if (MESS_IS_COMPLEX(A)){
                               for (i = 0; i < A->cols; i++){
                                   for (j=A->colptr[i] ; j < A->colptr[i+1]; j++){
                                       tmp[A->rowptr[j]] += cabs(A->values_cpx[j]);
                                       if ( (MESS_IS_HERMITIAN (A))&& A->rowptr[j] != i) {
                                           tmp[i] += cabs(A->values_cpx[j]);
                                       }
                                   }
                               }
                           }
                           max=tmp[0];
                           for (i = 1; i < A->rows; i++){
                               if (tmp[i] > max ) max= tmp[i];
                           }
                           mess_free(tmp);
                           break;
                       }
                       /*-----------------------------------------------------------------------------
                        *  COORD
                        *-----------------------------------------------------------------------------*/
        case MESS_COORD:{
                            int r;
                            mess_matrix tmp;
                            MESS_MATRIX_CHECKFORMAT(A, tmp, r, MESS_CSR);
                            if ( r == 0 ){
                                ret = mess_matrix_norminf(tmp, &max);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_norminf);
                                mess_matrix_clear(&tmp);
                            }else{
                                MSG_ERROR("error convert COORD to CSR\n");
                                ret = MESS_ERROR_CONVERT;
                                mess_matrix_clear(&tmp);
                            }
                            break;
                        }

                        /*-----------------------------------------------------------------------------
                         *  DENSE
                         *-----------------------------------------------------------------------------*/
        case MESS_DENSE: {
                             double *work;
                             mess_try_alloc(work, double *, sizeof(double) * A->rows);
                             if(MESS_IS_REAL(A)){
                                 max = F77_GLOBAL(dlange,DLANGE)("I",&(A->rows),&(A->cols),A->values,&(A->ld),work);
                             }else{
                                 max = F77_GLOBAL(zlange,ZLANGE)("I",&(A->rows),&(A->cols),A->values_cpx,&(A->ld),work);
                             }
                             mess_free(work);
                             break;
                         }

        default:
                         MSG_ERROR("unkown/unsupported storage type: %s\n", mess_storage_t_str(A->store_type));
                         return( MESS_ERROR_STORAGETYPE);
    }
    *nrm  = max;
    return(ret);
}

/**
 * @brief Compute the square of the Frobenius-norm of a matrix.
 * @param[in] A input matrix \f$A\f$
 * @param[out] nrm  output squared Frobenius-norm
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_normf2 function calculates the square of Frobenius-norm of a matrix:
 * \f[ nrm= \sum\limits_{i,j=1}^{m,n} |A_{i,j}|^2  = \Vert A \Vert_F^2 . \f]
 * It is used by \ref mess_matrix_normf to compute the Frobenius-norm.
 *
 * @sa mess_matrix_norm1
 * @sa mess_matrix_norm2
 * @sa mess_matrix_norminf
 */
int mess_matrix_normf2(mess_matrix A, double *nrm){
    MSG_FNAME(__func__);
    double norm = 0.0;
    mess_int_t i=0, j=0;
    int ret =0 ;
    mess_check_nullpointer(A);
    mess_check_nullpointer(nrm);
    mess_check_real_or_complex(A);


    switch ( A->store_type){
        case MESS_CSR: {
                           if (MESS_IS_REAL(A)) {
                               for ( i = 0 ; i < A->rows; i++){
                                   for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                                       norm += A->values[j]*A->values[j];
                                       if ( (MESS_IS_SYMMETRIC (A)|| MESS_IS_SKEWSYMMETRIC(A))&& A->colptr[j] != i) {
                                           norm += A->values[j]*A->values[j];
                                       }
                                   }
                               }
                           } else {
                               for ( i = 0 ; i < A->rows; i++){
                                   for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                                       norm += cabs(A->values_cpx[j])*cabs(A->values_cpx[j]);
                                       if ( (MESS_IS_SYMMETRIC (A)|| MESS_IS_SKEWSYMMETRIC(A))&& A->colptr[j] != i) {
                                           norm += cabs(A->values_cpx[j])*cabs(A->values_cpx[j]);
                                       }
                                   }
                               }
                           }
                           break;
                       }
        case MESS_CSC: {
                           if ( MESS_IS_REAL(A)){
                               for ( i = 0 ; i < A->cols; i++){
                                   for ( j = A->colptr[i]; j < A->colptr[i+1]; j++){
                                       norm += A->values[j]*A->values[j];
                                       if ( (MESS_IS_SYMMETRIC (A)|| MESS_IS_SKEWSYMMETRIC(A))&& A->rowptr[j] != i) {
                                           norm += A->values[j]*A->values[j];
                                       }
                                   }
                               }
                           }else {
                               for ( i = 0 ; i < A->cols; i++){
                                   for ( j = A->colptr[i]; j < A->colptr[i+1]; j++){
                                       norm += cabs(A->values_cpx[j])*cabs(A->values_cpx[j]);
                                       if ( (MESS_IS_SYMMETRIC (A)|| MESS_IS_SKEWSYMMETRIC(A))&& A->rowptr[j] != i) {
                                           norm += cabs(A->values_cpx[j])*cabs(A->values_cpx[j]);
                                       }
                                   }
                               }
                           }
                           break;
                       }
        case MESS_COORD:{
                            int r;
                            mess_matrix tmp;
                            MESS_MATRIX_CHECKFORMAT(A, tmp, r, MESS_CSR);
                            if ( r == 0 ){
                                ret = mess_matrix_normf2(tmp, &norm);   FUNCTION_FAILURE_HANDLE(ret, (ret !=0), mess_matrix_normf2);
                                mess_matrix_clear(&tmp);
                            }else{
                                MSG_ERROR("error convert COORD to CSR\n");
                                return(MESS_ERROR_CONVERT);
                            }
                            break;
                        }
        case MESS_DENSE:{
                            if(MESS_IS_REAL(A)){
                                norm = F77_GLOBAL(dlange,DLANGE)("F",&(A->rows),&(A->cols),A->values,&(A->ld),NULL);
                            }else{
                                norm = F77_GLOBAL(zlange,ZLANGE)("F",&(A->rows),&(A->cols),A->values_cpx,&(A->ld),NULL);
                            }
                            norm *=norm;
                            break;
                        }
        default:
                        MSG_ERROR("unkown/unsupported storage type: %s\n", mess_storage_t_str(A->store_type));
                        return(MESS_ERROR_STORAGETYPE);
    }
    *nrm = norm;
    return(0);
}


/**
 * @brief Compute the square of the Frobenius-norm of a \f$op1(A)*op2(B)\f$.
 * @param[in] op1 operation of \f$A\f$
 * @param[in] A input matrix \f$A\f$
 * @param[in] op2 operation of \f$B\f$
 * @param[in] B input matrix \f$B\f$
 * @param[out] nrm      output squared Frobenius-norm
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_mulnormf2 function calculates:
 * \f[ nrm= \Vert op(A)*op(B) \Vert_F^2 . \f]
 *
 * @sa mess_matrix_norm1
 * @sa mess_matrix_norm2
 * @sa mess_matrix_norminf
 */
int mess_matrix_mulnormf2(mess_operation_t op1, mess_matrix A, mess_operation_t op2, mess_matrix B, double *nrm){
    MSG_FNAME(__func__);
    mess_int_t ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_real_or_complex(A);
    mess_check_real_or_complex(B);
    mess_check_operation_type(op1);
    mess_check_operation_type(op2);
    mess_check_nullpointer(nrm);

    /*-----------------------------------------------------------------------------
     *  perform operations
     *-----------------------------------------------------------------------------*/
    mess_matrix tmp;
    ret = mess_matrix_init(&tmp);                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_multiply(op1,A,op2,B,tmp);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
    ret = mess_matrix_normf2(tmp,nrm);                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf2);
    ret = mess_matrix_clear(&tmp);                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);

    return 0;
}



/**
 * @brief Compute the Frobenius-norm of a matrix.
 * @param[in] A input matrix \f$A\f$
 * @param[out] nrm output Frobenius-norm
 *
 * The @ref mess_matrix_normf function calculates the Frobenius-norm of a matrix:
 * \f[ nrm= \Vert A \Vert_F . \f]
 * It uses the @ref mess_matrix_normf2 function internally.
 *
 * @sa mess_matrix_normf2
 * @sa mess_matrix_norm1
 * @sa mess_matrix_norm2
 * @sa mess_matrix_norminf
 */
int mess_matrix_normf(mess_matrix A, double *nrm ){
    MSG_FNAME(__func__);
    int ret =0;

    mess_check_nullpointer(A);
    mess_check_nullpointer(nrm);
    mess_check_real_or_complex(A);

    *nrm = 0.0;
    ret = mess_matrix_normf2(A, nrm); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_normf2);
    *nrm=sqrt(*nrm);
    return(0);
}

struct norm2data {
    mess_matrix C;
    mess_vector x1;
};

static int norm2mvp(void *data, mess_operation_t op,  mess_vector x, mess_vector y) {
    MSG_FNAME(__func__);

    struct norm2data * d = (struct norm2data*) data;
    int ret = 0;
    ret = mess_matrix_mvp(MESS_OP_NONE,d->C, x, d->x1);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->C,  d->x1,y);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(y);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    // mess_vector_printshort(y);
    return(0);
}



/**
 * @brief Compute the \f$ 2 \f$-norm of a matrix.
 * @param[in] A input matrix \f$ A\f$
 * @param[out] nrm   output \f$2\f$-norm
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_norm2 function calculates the \f$ 2 \f$-norm of a matrix, i.e. the largest eigenvalue of
 * \f$ A^HA \f$.
 * Depending on the storage format of the matrix the SVD or a iterative scheme will be used. \n
 * Even if the matrix is dense and larger than \f$ 100 \f$ rows or columns the iterative scheme is used, too.
 *
 * @sa mess_matrix_norm1
 * @sa mess_matrix_normf
 * @sa mess_matrix_norminf
 * @sa mess_eigen_svd
 * @sa mess_eigen_arnoldi_template_nrm
 */
int mess_matrix_norm2(mess_matrix A, double *nrm){
    MSG_FNAME(__func__);
    struct norm2data dat;
    mess_mvpcall mvpcall;
    mess_vector sv = NULL;
    int ret ;
    mess_int_t it;

    mess_check_nullpointer(A);
    mess_check_nullpointer(nrm);
    mess_check_real_or_complex(A);

    if ( MESS_MIN(A->rows, A->cols) < 500 || MESS_IS_DENSE(A)  ) {
        mess_vector svd;
        mess_vector_init(&svd);
        mess_vector_alloc(svd, MESS_MIN(A->rows, A->cols), MESS_REAL);
        ret = mess_eigen_svd (A, svd, NULL, NULL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_svd);
        *nrm = svd->values[0];
        mess_vector_clear(&svd);
    } else {
        ret = mess_vector_init(&(dat.x1));                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc((dat.x1), A->rows, A->data_type);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        dat.C=A;
        ret = mess_mvpcall_operator(&mvpcall, A->cols, A->data_type, norm2mvp, &dat);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_mvpcall_operator);

        ret = mess_vector_init(&sv);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(sv, A->cols, A->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_rand(sv);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_rand);
        it = 50 ;
        if ( it > A->cols) it = A->cols;
        ret = mess_eigen_arnoldi_template_nrm(mvpcall, it, sv, nrm); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);

        *nrm = sqrt (*nrm);
        mess_vector_clear(&(dat.x1));
        mess_vector_clear(&sv);
        mess_mvpcall_clear(&mvpcall);
    }

    return(0);
}

struct norminv2data{
    mess_vector x1;
    mess_direct sol;
};

static int norminv2mvp(void *data, mess_operation_t op, mess_vector x, mess_vector y) {
    MSG_FNAME(__func__);
    struct norminv2data * d = (struct norminv2data*) data;
    int ret = 0;
    ret = mess_direct_solve(MESS_OP_NONE,d->sol,x,d->x1);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solve);
    ret = mess_direct_solve(MESS_OP_HERMITIAN,d->sol,d->x1,y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solve);
    ret = mess_vector_toreal_nowarn(y);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    // mess_vector_printshort(y);
    return (0);
}



/**
 * @brief Compute the \f$ 2 \f$-norm of the inverse of a matrix.
 * @param[in] A input matrix \f$A\f$
 * @param[out] nrm output \f$ 2\f$-norm
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_norm2inv function computes the \f$ 2 \f$-norm of \f$ A^{-1}\f$. \n
 * If the matrix is small SVD is used to compute the norm. Otherwise a LU decomposition of the matrix is computed
 * and the Arnoldi process is used to estimate norm of the inverse.
 *
 * @sa mess_matrix_norm2
 * @sa mess_eigen_svd
 * @sa mess_eigen_arnoldi_template_nrm
 * @sa mess_direct_lu
 *
 */
int mess_matrix_norm2inv(mess_matrix A, double *nrm){
    MSG_FNAME(__func__);
    struct norminv2data dat;
    mess_mvpcall mvpcall;
    mess_vector sv = NULL;
    int ret ;
    mess_int_t it;

    mess_check_nullpointer(A);
    mess_check_nullpointer(nrm);
    mess_check_real_or_complex(A);

    if ( MESS_MIN(A->rows, A->cols) < 100 ) {
        mess_vector svd;
        ret = mess_vector_init(&svd);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(svd, MESS_MIN(A->rows, A->cols), MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_eigen_svd (A, svd, NULL, NULL);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_svd);
        *nrm = 1.0/svd->values[svd->dim-1];
        mess_vector_clear(&svd);
    } else {
        mess_direct sol;
        ret = mess_direct_init(&sol);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
        ret = mess_direct_lu(A, sol);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);
        ret = mess_vector_init(&(dat.x1));                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc((dat.x1), A->rows, A->data_type);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        dat.sol=sol;
        ret = mess_mvpcall_operator(&mvpcall, A->cols, A->data_type, norminv2mvp, &dat);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_mvpcall_operator);

        ret = mess_vector_init(&sv);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(sv, A->cols, A->data_type);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_rand(sv);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_rand);
        it = 50 ;
        if ( it > A->cols) it = A->cols;

        ret = mess_eigen_arnoldi_template_nrm(mvpcall, it, sv, nrm); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);

        *nrm = sqrt (*nrm);
        mess_vector_clear(&(dat.x1));
        mess_vector_clear(&sv);
        mess_direct_clear(&sol);
        mess_mvpcall_clear(&mvpcall);
    }

    return(0);
}


/**
 * @brief Compute the norm of a matrix.
 * @param[in] A         input matrix \f$A\f$
 * @param[in] nrm_t     input @ref mess_norm_t the desired type of norm
 * @param[out] nrm      output  \f$ \Vert A \Vert \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_norm function calculates the norm of a matrix.
 * Supported norms are:
 *
 * @ref MESS_2_NORM, @ref MESS_FROBENIUS_NORM, @ref MESS_1_NORM, @ref MESS_INF_NORM.
 *
 * @sa mess_matrix_norm1
 * @sa mess_matrix_norm2
 * @sa mess_matrix_normf
 * @sa mess_matrix_norminf
 */
int mess_matrix_norm(mess_matrix A , mess_norm_t nrm_t, double *nrm){
    MSG_FNAME(__func__);
    int ret;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(nrm);


    /*-----------------------------------------------------------------------------
     *  compute norm
     *-----------------------------------------------------------------------------*/
    switch (nrm_t) {
        case MESS_2_NORM:
            ret = mess_matrix_norm2(A, nrm);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm2);
            break;
        case MESS_FROBENIUS_NORM:
            ret = mess_matrix_normf(A, nrm);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_normf);
            break;
        case MESS_1_NORM:
            ret = mess_matrix_norm1(A, nrm);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm1);
            break;
        case MESS_INF_NORM:
            ret = mess_matrix_norminf(A, nrm);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norminf);
            break;
        default:
            MSG_ERROR("unkown/unsupported norm type: %s\n", mess_norm_t_str(nrm_t));
            return MESS_ERROR_NOSUPPORT;
    }
    return(ret);
}

