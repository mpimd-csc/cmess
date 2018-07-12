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
 * @file lib/direct/dgeqp3.c
 * @brief Interface for @c DGEQP3.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/eigenvalue.h"
#include "mess/error_macro.h"
#include <complex.h>
#include <math.h>
#include "blas_defs.h"


/**
 * @brief Interface for @lapack @c DGEQP3.
 * @param[in] A         input matrix
 * @param[out] Q        output orthogonal matrix
 * @param[out] R        output  upper triangular matrix
 * @param[out] perm     output permutation
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_dgeqp3 function provides a frontend to the @c DGEQP3 routine from
 * @lapack.
 * If a parameter points to @c NULL, \f$ R \f$  is not computed.
 *
 */
int mess_direct_dgeqp3 ( mess_matrix A, mess_matrix Q, mess_matrix R, mess_int_t *perm)
{
    MSG_FNAME(__func__);
    mess_matrix wA;
    mess_int_t LDA;
    mess_int_t M,N;
    double  *work = NULL;
    double *TAU = NULL ;
    mess_int_t *JPVT = NULL ;
    mess_int_t LWORK;
    mess_int_t info = 0;
    int ret = 0;
    mess_int_t i,j;
    double ws;
    // MSG_ERROR("eps: %lg\n", eps);
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_real(A);

    ret = mess_matrix_init(&wA); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_convert(A,wA,MESS_DENSE); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);

    LDA=A->ld;
    M = A->rows;
    N = A->cols;

    /*-----------------------------------------------------------------------------
     *  WORK SPACE QUERY
     *-----------------------------------------------------------------------------*/
    LWORK = -1;
    F77_GLOBAL(dgeqp3,DGEQP3)(&M,&N,wA->values,&LDA,JPVT,TAU, &ws, &LWORK, &info);
    // MSG_INFO("Workspace query (DGEQP3): " MESS_PRINTF_INT "\n", (mess_int_t) ws);



    // LWORK = 2*N+(N+1)*32;
    LWORK = (mess_int_t) 2*ws +2;
    mess_try_alloc(TAU, double *, sizeof(double)* MESS_MIN(M, N));
    mess_try_alloc(JPVT, mess_int_t *, sizeof(mess_int_t)*N);
    mess_try_alloc(work, double *, sizeof(double)*LWORK);

    for ( i = 0 ; i < N; i++){
        JPVT[i]=0;
    }


    /*-----------------------------------------------------------------------------
     *  run DGEQP3
     *-----------------------------------------------------------------------------*/
    F77_GLOBAL(dgeqp3,DGEQP3)(&M,&N,wA->values,&LDA,JPVT,TAU, work, &LWORK, &info);
    if ( info !=0){
        MSG_ERROR("LAPACK DGEQP3 returned with error: " MESS_PRINTF_INT "\n", info);
        mess_free(work); mess_free(JPVT); mess_free(TAU);
        mess_matrix_clear(&wA);
        return MESS_ERROR_LAPACK;
    }
    mess_free(work);


    //  for ( i = 0 ; i < N; i++){
    //          MSG_PRINT("JPVT["MESS_PRINTF_INT"]="MESS_PRINTF_INT"\n",i,JPVT[i]);
    //  }


    /*-----------------------------------------------------------------------------
     *  extract R
     *-----------------------------------------------------------------------------*/
    if (R!=NULL) {
        ret= mess_matrix_alloc(R, MESS_MIN(M, N),N, MESS_MIN(M, N)*N, MESS_DENSE,MESS_REAL);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        for (i=0;i<R->rows;i++){
            for (j=i;j < R->cols; j++){
                R->values[i+j*R->ld]=wA->values[i+j*wA->ld];
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  fix the permutation
     *-----------------------------------------------------------------------------*/
    if ( perm != NULL){
        for ( i = 0 ; i < N; i++){
            perm[i] = (mess_int_t) JPVT[i]-1;
        }
    }

    /*-----------------------------------------------------------------------------
     *  extract Q
     *-----------------------------------------------------------------------------*/
    if ( Q!=NULL){
        ret = mess_matrix_alloc(Q,M,M,M*M, MESS_DENSE, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        for (i  =0 ; i < M; i++){ Q->values[i*Q->ld+i] = 1.0; }
        mess_int_t K = MESS_MIN(M, N);
        LWORK = -1;
        F77_GLOBAL(dormqr,DORMQR)("L","N", &M, &M, &K, wA->values, &LDA,TAU,Q->values,&M,&ws, &LWORK,&info);
        LWORK = (mess_int_t) ws +1;
        // MSG_INFO("DORMQR workspace = " MESS_PRINTF_INT "\n", LWORK );
        mess_try_alloc(work, double *, sizeof(double) *(LWORK));
        F77_GLOBAL(dormqr,DORMQR)("L","N", &M, &M, &K, wA->values, &LDA,TAU,Q->values,&M,work, &LWORK,&info);

        if ( info!=0){
            MSG_ERROR("LAPACK DORMQR returned with error: " MESS_PRINTF_INT "\n", info);
            mess_free(work); mess_free(JPVT); mess_free(TAU);
            mess_matrix_clear(&wA);
            return MESS_ERROR_LAPACK;
        }

    }
    MESS_CLEAR_POINTERS(work,JPVT,TAU);
    mess_matrix_clear(&wA);
    return 0;
}       /* -----  end of function mess_direct_dgeqp3  ----- */


