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
 * @file lib/eigen/svd.c
 * @brief Compute the singular value decompositon (SVD).
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/eigenvalue.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>


/**
 * @brief Compute the singular value decompositon (SVD) of a matrix.
 * @param[in] A    input matrix
 * @param[out] S vector containing singular values of \f$ A \f$
 * @param[out] U left singular vectors
 * @param[out] V right singular vectors
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_svd function computes the singular value decompositon (SVD) of \f$ A \f$: \n
 * \f[ A= U S V^T . \f]
 * If \f$ U \f$ or \f$ V \f$ is @c NULL the corresponding singular vector are not computed. \n
 * To compute the economy size SVD of a matrix please use \ref mess_eigen_svd_econ.
 *
 */
int mess_eigen_svd ( mess_matrix A, mess_vector S, mess_matrix U, mess_matrix V )
{
    MSG_FNAME(__func__);
    mess_matrix work, VT;
    int leftV =0;
    int rightV = 0;
    int ret = 0;
    mess_int_t M,N,LDA,LDU,LDVT,LWORK;
    mess_int_t INFO = 0;

    char JOBU[4];
    char JOBV[4];

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(S);

    if ( U != NULL){
        leftV=1;
        strncpy(JOBU,"A",3);
    } else {
        strncpy(JOBU,"N",3);
    }
    if ( V != NULL){
        rightV=1;
        strncpy(JOBV,"A",3);
    } else {
        strncpy(JOBV,"N",3);
    }

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&work);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&VT);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_convert(A,work, MESS_DENSE);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);


    M = A->rows;
    N = A->cols;
    LDA = MESS_MAX(1,work->ld);
    ret = mess_vector_resize(S,MESS_MIN(M,N));      FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_resize);
    ret = mess_vector_toreal(S);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_toreal);

    LDU = 1;
    LDVT = 1;
    if (MESS_IS_REAL(A)) {
        double *uvalues, *vtvalues;
        double *workfield;
        double ws;
        if ( leftV == 1){
            ret = mess_matrix_alloc(U,A->rows,A->rows,A->rows*A->rows,MESS_DENSE, MESS_REAL);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
            uvalues = U->values;
            LDU = U->ld;
        } else {
            uvalues = NULL;
        }
        if ( rightV ==1 ){
            ret =mess_matrix_alloc(VT,N,N,N*N,MESS_DENSE, MESS_REAL);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
            vtvalues=VT->values;
            LDVT = VT->ld;
        } else {
            vtvalues = NULL;
        }

        /*-----------------------------------------------------------------------------
         *  Workspace Query
         *-----------------------------------------------------------------------------*/
        LWORK = -1;
        F77_GLOBAL(dgesvd,DGESVD)(JOBU, JOBV,&M,&N,work->values, &LDA,S->values, uvalues, &LDU, vtvalues, &LDVT,&ws, &LWORK, &INFO);
        LWORK = nearbyint(ws+1);
        mess_try_alloc(workfield, double *, sizeof(double)*LWORK);

        /*-----------------------------------------------------------------------------
         *  SVD
         *-----------------------------------------------------------------------------*/
        F77_GLOBAL(dgesvd,DGESVD)(JOBU, JOBV,&M,&N,work->values, &LDA,S->values, uvalues, &LDU, vtvalues, &LDVT,workfield, &LWORK, &INFO);
        if ( INFO !=0){
            MSG_ERROR("DGESVD returned with error " MESS_PRINTF_INT "\n", INFO);
            return MESS_ERROR_LAPACK;
        }

        /*------------------------------------------------------------------------
         * finalize
         *-----------------------------------------------------------------------------*/
        if ( rightV == 1) {
            ret = mess_matrix_ctranspose(VT,V);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_ctranspose);
        }

        mess_matrix_clear(&work);
        mess_matrix_clear(&VT);
        mess_free(workfield);
    } else {
        double *rworkfield =  NULL;
        mess_double_cpx_t *uvalues, *vtvalues, *workfield;
        mess_double_cpx_t ws = 0;

        if ( leftV == 1){
            ret = mess_matrix_alloc(U,A->rows,A->rows,A->rows*A->rows,MESS_DENSE, MESS_COMPLEX);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
            uvalues = U->values_cpx;
            LDU=U->ld;
        } else {
            uvalues = NULL;
        }
        if ( rightV ==1 ){
            ret =mess_matrix_alloc(VT,N,N,N*N,MESS_DENSE, MESS_COMPLEX);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
            vtvalues=VT->values_cpx;
            LDVT=VT->ld;
        } else {
            vtvalues = NULL;
        }

        /*-----------------------------------------------------------------------------
         *  Workspace Query
         *-----------------------------------------------------------------------------*/
        LWORK = -1;
        F77_GLOBAL(zgesvd,ZGESVD)(JOBU, JOBV,&M,&N,work->values_cpx, &LDA,S->values, uvalues, &LDU, vtvalues, &LDVT,&ws, &LWORK, rworkfield, &INFO);
        LWORK =  nearbyint(creal(ws)+1);
        mess_try_alloc(workfield, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*LWORK);
        mess_try_alloc(rworkfield, double * , sizeof(double) * 5* MESS_MIN(M,N));
        /*-----------------------------------------------------------------------------
         *  SVD
         *-----------------------------------------------------------------------------*/
        F77_GLOBAL(zgesvd,ZGESVD)(JOBU, JOBV,&M,&N,work->values_cpx, &LDA,S->values, uvalues, &LDU, vtvalues, &LDVT,workfield, &LWORK, rworkfield, &INFO);
        if ( INFO !=0){
            MSG_ERROR("ZGESVD returned with error " MESS_PRINTF_INT "\n", INFO);
            return MESS_ERROR_LAPACK;
        }
        /*------------------------------------------------------------------------
         * finalize
         *-----------------------------------------------------------------------------*/
        if ( rightV == 1) {
            ret = mess_matrix_ctranspose(VT,V);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_ctranspose);
        }

        mess_matrix_clear(&work);
        mess_matrix_clear(&VT);
        mess_free(workfield);
        mess_free(rworkfield);
    }
    return 0;
}       /* -----  end of function mess_eigen_svd  ----- */


/**
 * @brief Compute the singular value decompositon (SVD) of a matrix.
 * @param[in] A   input matrix
 * @param[out] S vector containing singular values of \f$ A \f$
 * @param[out] U left singular vectors
 * @param[out] V right singular vectors
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_svd_econ function computes the economy size singular value decompositon (SVD) of \f$ A \f$:
 * \f[ A = U S V^T .\f]
 * If \f$ U \f$ or \f$ V \f$ is @c NULL the corresponding singular vector are not computed. \n
 * This function is similar to @matlab  \c svd(X,0) or \c svd(X,'econ'). \n
 * To compute
 * the normal SVD please use \ref mess_eigen_svd.
 *
 */
int mess_eigen_svd_econ( mess_matrix A, mess_vector S, mess_matrix U, mess_matrix V )
{
    MSG_FNAME(__func__);
    mess_matrix work, VT;
    int leftV =0;
    int rightV = 0;
    int ret = 0;
    mess_int_t M,N,LDA,LDU,LDVT,LWORK;
    mess_int_t INFO = 0;

    char JOBU[4];
    char JOBV[4];

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(S);

    if ( U != NULL){
        leftV=1;
        strncpy(JOBU,"S",3);
    } else {
        strncpy(JOBU,"N",3);
    }
    if ( V != NULL){
        rightV=1;
        strncpy(JOBV,"S",3);
    } else {
        strncpy(JOBV,"N",3);
    }

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&work);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&VT);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_convert(A,work, MESS_DENSE);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);


    M = A->rows;
    N = A->cols;
    LDA = MESS_MAX(1,work->ld);
    ret = mess_vector_resize(S,MESS_MIN(M,N));          FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_resize);
    ret = mess_vector_toreal(S);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_toreal);

    LDU = 1;
    LDVT=1;
    if (MESS_IS_REAL(work)) {
        double *workfield;
        double *uvalues, *vtvalues;
        double ws;

        if ( leftV == 1){
            ret = mess_matrix_alloc(U,A->rows,MESS_MIN(M,N),A->rows*MESS_MIN(M,N),MESS_DENSE, MESS_REAL);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
            uvalues = U->values;
            LDU = U->ld;
        } else {
            uvalues = NULL;
        }
        LDVT = MESS_MIN(M,N);
        if ( rightV ==1 ){
            ret =mess_matrix_alloc(VT,MESS_MIN(N,M),N,MESS_MIN(M,N)*N,MESS_DENSE, MESS_REAL);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
            vtvalues=VT->values;
            LDVT=VT->ld;
        } else {
            vtvalues = NULL;
        }

        /*-----------------------------------------------------------------------------
         *  Workspace Query
         *-----------------------------------------------------------------------------*/
        LWORK = -1;
        F77_GLOBAL(dgesvd,DGESVD)(JOBU, JOBV,&M,&N,work->values, &LDA,S->values, uvalues, &LDU, vtvalues, &LDVT,&ws, &LWORK, &INFO);
        LWORK = nearbyint(ws+1);
        mess_try_alloc(workfield, double *, sizeof(double)*LWORK);

        /*-----------------------------------------------------------------------------
         *  SVD
         *-----------------------------------------------------------------------------*/
        F77_GLOBAL(dgesvd,DGESVD)(JOBU, JOBV,&M,&N,work->values, &LDA,S->values, uvalues, &LDU, vtvalues, &LDVT,workfield, &LWORK, &INFO);
        if ( INFO !=0){
            MSG_ERROR("DGESVD returned with error " MESS_PRINTF_INT "\n", INFO);
            return MESS_ERROR_LAPACK;
        }

        /*------------------------------------------------------------------------
         * finalize
         *-----------------------------------------------------------------------------*/
        if ( rightV == 1) {
            ret = mess_matrix_ctranspose(VT,V);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_ctranspose);
        }

        mess_matrix_clear(&work);
        mess_matrix_clear(&VT);
        mess_free(workfield);
    } else {
        double *rworkfield = NULL;
        mess_double_cpx_t *uvalues, *vtvalues, *workfield;
        mess_double_cpx_t ws = 0;

        if ( leftV == 1){
            ret = mess_matrix_alloc(U,A->rows,MESS_MIN(M,N),MESS_MIN(M,N)*A->rows,MESS_DENSE, MESS_COMPLEX);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
            uvalues = U->values_cpx;
            LDU = U->ld;
        } else {
            uvalues = NULL;
        }
        LDVT = MESS_MIN(M,N);
        if ( rightV ==1 ){
            ret =mess_matrix_alloc(VT,MESS_MIN(N,M),N,MESS_MIN(N,M)*N,MESS_DENSE, MESS_COMPLEX);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
            vtvalues=VT->values_cpx;
            LDVT = VT->ld;
        } else {
            vtvalues = NULL;
        }

        /*-----------------------------------------------------------------------------
         *  Workspace Query
         *-----------------------------------------------------------------------------*/
        LWORK = -1;
        F77_GLOBAL(zgesvd,ZGESVD)(JOBU, JOBV,&M,&N,work->values_cpx, &LDA,S->values, uvalues, &LDU, vtvalues, &LDVT,&ws, &LWORK, rworkfield, &INFO);
        LWORK = nearbyint(creal(ws)+1);
        mess_try_alloc(workfield, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*LWORK);
        mess_try_alloc(rworkfield, double * , sizeof(double) * 5* MESS_MIN(M,N));

        /*-----------------------------------------------------------------------------
         *  SVD
         *-----------------------------------------------------------------------------*/
        F77_GLOBAL(zgesvd,ZGESVD)(JOBU, JOBV,&M,&N,work->values_cpx, &LDA,S->values, uvalues, &LDU, vtvalues, &LDVT,workfield, &LWORK, rworkfield, &INFO);
        if ( INFO !=0){
            MSG_ERROR("ZGESVD returned with error " MESS_PRINTF_INT "\n", INFO);
            return MESS_ERROR_LAPACK;
        }

        /*------------------------------------------------------------------------
         * finalize
         *-----------------------------------------------------------------------------*/
        if ( rightV == 1) {
            ret = mess_matrix_ctranspose(VT,V);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_ctranspose);
        }

        mess_matrix_clear(&work);
        mess_matrix_clear(&VT);
        mess_free(workfield);
        mess_free(rworkfield);
    }
    return 0;
}       /* -----  end of function mess_eigen_svd  ----- */


