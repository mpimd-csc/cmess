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
 * @file lib/matrix/dynorm.c
 * @brief Norm computation of a dyadic matrix product.
 * @author @koehlerm
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/misc.h"
#include "mess/eigenvalue.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>


/**
 * @brief Compute the \f$ 2 \f$-norm of \f$ A^HA \f$ or \f$ AA^H \f$.
 * @param[in] A     input matrix \f$A\f$
 * @param[out] nrm  output \f$2\f$-norm
 * @return zero on success or a non zero error code
 *
 * The @ref mess_matrix_dynorm2 function computes the \f$ 2 \f$-norm of \f$ A^HA \f$ or \f$ AA^H \f$  by using the equality:
 * \f[ \Vert A^HA \Vert_2 = \Vert AA^H \Vert_2 \f]
 * and exploiting the smaller sized product.
 *
 * @see mess_matrix_dynormf
 */
int mess_matrix_dynorm2 ( mess_matrix A, double *nrm )
{
    MSG_FNAME(__func__);
    int ret = 0;
    double *w, work, *workspace;
    mess_int_t lwork=-1, info=0,conv;
    mess_matrix mwork, tmp;

    /*-----------------------------------------------------------------------------
     *  Check Input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(nrm);
    mess_check_real_or_complex(A);

    /*-----------------------------------------------------------------------------
     *  compute small matrix
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&tmp);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init );
    if ( A->rows < A->cols) {
        ret = mess_matrix_multiply(MESS_OP_NONE, A, MESS_OP_HERMITIAN, A, tmp); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    } else {
        ret = mess_matrix_multiply(MESS_OP_HERMITIAN, A, MESS_OP_NONE, A, tmp); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    }

    /*-----------------------------------------------------------------------------
     *  compute eigenvalues
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(w, double*, sizeof(double)*tmp->rows);
    mess_matrix_toreal(tmp);
    MESS_MATRIX_CHECKFORMAT(tmp, mwork, conv, MESS_DENSE);
    F77_GLOBAL(dsyev,DSYEV)("N", "U", &(mwork->rows), mwork->values, &(mwork->ld), w, &work, &lwork, &info);
    lwork = nearbyint(work+1);
    mess_try_alloc(workspace, double*, sizeof(double)*(lwork));
    F77_GLOBAL(dsyev,DSYEV)("N", "U", &(mwork->rows), mwork->values, &(mwork->ld), w, workspace, &lwork, &info);
    mess_free(workspace);
    if(conv==0){
        MESS_CLEAR_MATRICES(&mwork);
    }

    if ( info != 0 ) {
        MESS_CLEAR_MATRICES(&tmp);
        mess_free(w);
        MSG_ERROR("An error occured in DSYEV: " MESS_PRINTF_INT "\n", info);
        return MESS_ERROR_LAPACK;
    }
    //get largest one
    *nrm = fabs(w[tmp->rows-1]);
    MESS_CLEAR_MATRICES(&tmp);
    mess_free(w);
    MESS_CLEAR_MATRICES(&tmp);
    return 0;
}       /* -----  end of function mess_norm_dynorm2  ----- */

/**
 * @brief Compute the squared Frobenius-norm of \f$ A^HA \f$ or \f$ AA^H \f$.
 * @param[in] A     input matrix  \f$A\f$
 * @param[out] nrm  output  Frobenius-norm
 *
 * @return zero on success or a non zero error code
 *
 * The @ref mess_matrix_dynormf2 function computes the squared Frobenius-norm of \f$ A^HA \f$ or \f$ AA^H \f$  by using the equality:
 * \f[ \Vert A^HA \Vert_F^2 = \Vert AA^H \Vert_F^2 \f]
 * and exploiting the smaller sized product.
 *
 * @see mess_matrix_dynorm2
 * @see mess_matrix_dynormf
 */
int mess_matrix_dynormf2 ( mess_matrix A, double *nrm )
{
    MSG_FNAME(__func__);
    mess_matrix tmp;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  Check Input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(nrm);
    mess_check_real_or_complex(A);

    /* Form the small matrix  */
    ret = mess_matrix_init(&tmp);                                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    if ( A->rows < A->cols) {
        ret = mess_matrix_multiply(MESS_OP_NONE, A, MESS_OP_HERMITIAN, A, tmp);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    } else {
        ret = mess_matrix_multiply(MESS_OP_HERMITIAN, A, MESS_OP_NONE, A, tmp);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    }
    ret = mess_matrix_normf2(tmp, nrm);                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_normf2);
    mess_matrix_clear(&tmp);

    return 0;
}       /* -----  end of function mess_norm_dynormf2  ----- */


/**
 * @brief Compute the Frobenius-norm of \f$ A^HA \f$ or \f$ AA^H \f$.
 * @param[in] A     input matrix  \f$A\f$
 * @param[out] nrm  output  Frobenius-norm
 *
 * @return zero on success or a non zero error code
 *
 * The @ref mess_matrix_dynormf function computes the  Frobenius-norm of \f$ A^HA \f$ or \f$ AA^H \f$  by using the equality:
 * \f[ \Vert A^HA \Vert_F = \Vert AA^H \Vert_F \f]
 * and exploiting the smaller sized product.
 *
 * @see mess_matrix_dynorm2
 * @see mess_matrix_dynormf2
 */
int mess_matrix_dynormf ( mess_matrix A, double *nrm )
{
    MSG_FNAME(__func__);
    mess_int_t ret =0;
    ret = mess_matrix_dynormf2(A,nrm);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_dynormf2);
    *nrm = sqrt(*nrm);


    return 0;
}       /* -----  end of function mess_norm_dynormf  ----- */



/**
 * @brief Compute the \f$2\f$-norm of \f$ WW^T-KK^T \f$,\f$ WW^H-KK^T \f$,\f$ WW^T-KK^H \f$  or \f$ WW^H-KK^H \f$.
 * @param[in] W     input matrix  \f$W\f$
 * @param[in] K     input matrix  \f$K\f$
 * @param[out] nrm  output  Frobenius-norm
 *
 * @return zero on success or a non zero error code
 *
 * The @ref mess_matrix_indefinite_dynorm2 function computes the \f$2\f$-norm of \f$ WW^T-KK^T \f$,\f$ WW^H-KK^T \f$,\f$ WW^T-KK^H \f$  or \f$ WW^H-KK^H \f$
 * by using the following facts:
 *
 * \f[
 * \begin{aligned}
 *  WW^T - KK^T                 &\phantom{:}=   \begin{bmatrix} W & K \end{bmatrix}
 *                                              \begin{bmatrix} I_m & 0 \\ 0 &  -I_p \end{bmatrix}
 *                                              \begin{bmatrix} W^T \\ K^T \end{bmatrix} \\
 *                              &:=U D U^T \\
 *  \sigma(UDU^T)\setminus\{0\} &\phantom{:}=   \sigma(U^TUD)\setminus\{0\} \\
 *  \|WW^T-KK^T\|_2 &\phantom{:}= max\{|\lambda|\ |\ \lambda \in \sigma(WW^T-KK)\}=max\{|\lambda|\ |\ \lambda \in \sigma (U^TUD)\}
 *  \end{aligned}
 *  \f]
 *
 * Analogue for the case \f$(\cdot)^H\f$ and the mixed case \f$(\cdot)^T,(\cdot)^H\f$. We also assume that the number of rows is smaller
 * than the number of columns for \f$W\f$ and \f$ K \f$.
 * For a real/complex matrix \f$W\f$/\f$K\f$ the operation \f$(\cdot)^T\f$/\f$(\cdot)^H\f$ is performed.
 *
 *
 * @see mess_matrix_dynorm2
 * @see mess_matrix_indefinite_dynorm2
 * @see mess_matrix_indefinite_dynormf
 *
 */
int mess_matrix_indefinite_dynorm2 ( mess_matrix W, mess_matrix K, double *nrm){
    MSG_FNAME(__func__);
    mess_int_t ret = 0;
    mess_matrix Wt,Kt;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(W); mess_check_nullpointer(K);
    mess_check_real_or_complex(W); mess_check_real_or_complex(K);

    MESS_INIT_MATRICES(&Wt,&Kt);

    if(W->rows < W->cols){
        MSG_INFO("W->rows: " MESS_PRINTF_INT " >= W->cols: " MESS_PRINTF_INT".\n Take Hermitian of W.\n", W->rows, W->cols);
        ret = mess_matrix_ctranspose(W,Wt);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    }else{
        ret = mess_matrix_copy(W,Wt);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    }

    if(K->rows < K->cols){
        MSG_INFO("K->rows: " MESS_PRINTF_INT " >= K->cols: " MESS_PRINTF_INT" \n Take Hermitian of K.\n", K->rows, K->cols);
        ret = mess_matrix_ctranspose(K,Kt);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    }else{
        ret = mess_matrix_copy(K,Kt);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    }


    /*-----------------------------------------------------------------------------
     *  compte eigenvalues of U^TUD
     *-----------------------------------------------------------------------------*/
    mess_matrix U,D,UTU;
    MESS_INIT_MATRICES(&U,&D,&UTU);
    ret = mess_matrix_cat(Wt,Kt,NULL,NULL,MESS_DENSE,U);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);
    ret = mess_matrix_eye(D,U->cols,U->cols,MESS_DENSE);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_eye);
    mess_int_t i = 0;
    for(i=W->cols;i<D->rows;i++){D->values[i+D->ld*i]*=-1;}
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, U, MESS_OP_NONE, U, UTU);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE,UTU,MESS_OP_NONE,D,U);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);

    mess_vector ev;
    ret = mess_vector_init(&ev);                                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
    ret = mess_vector_alloc(ev,U->rows,MESS_COMPLEX);                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_alloc);
    ret = mess_eigen_eig(U,ev,NULL);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_eigen_eig);
    ret = mess_vector_norminf(ev,nrm);                                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norminf);

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&U,&D,&UTU,&Wt,&Kt);
    MESS_CLEAR_VECTORS(&ev);

    return ret;
}

/**
 * @brief Compute the square of the Frobenius-norm of \f$ WW^T-KK^T \f$,\f$ WW^H-KK^T \f$,\f$ WW^T-KK^H \f$  or \f$ WW^H-KK^H \f$.
 * @param[in] W     input matrix  \f$W\f$
 * @param[in] K     input matrix  \f$K\f$
 * @param[out] nrm  output  Frobenius-norm
 *
 * @return zero on success or a non zero error code
 *
 * The @ref mess_matrix_indefinite_dynormf2 function computes the square of the Frobenius-norm
 * of \f$ WW^T-KK^T \f$,\f$ WW^H-KK^T \f$,\f$ WW^T-KK^H \f$  or \f$ WW^H-KK^H \f$
 * by using the following facts:
 *
 * \f[
 * \begin{aligned}
 *  WW^T - KK^T                 &\phantom{:}=   \begin{bmatrix} W & K \end{bmatrix}
 *                                              \begin{bmatrix} I_m & 0 \\ 0 &  -I_p \end{bmatrix}
 *                                              \begin{bmatrix} W^T \\ K^T \end{bmatrix} \\
 *                              &:=U D U^T \\
 *  \sigma(UDU^T)\setminus\{0\} &\phantom{:}=   \sigma(U^TUD)\setminus\{0\} \\
 *  \|WW^T-KK^T\|_F^2 &\phantom{:}= \sum\limits_{\lambda_i\in\sigma(UDU^T) }\lambda_i^2 = \sum\limits_{\lambda_i\in\sigma(U^TUD)}\lambda_i^2
 *  \end{aligned}
 *  \f]
 *
 * Analogue for the case \f$(\cdot)^H\f$ and the mixed case \f$(\cdot)^T,(\cdot)^H\f$. We also assume that the number of rows is smaller
 * than the number of columns for \f$W\f$ and \f$ K \f$.
 * For a real/complex matrix \f$W\f$/\f$K\f$ the operation \f$(\cdot)^T\f$/\f$(\cdot)^H\f$ is performed.
 *
 *
 * @see mess_matrix_dynorm2
 * @see mess_matrix_indefinite_dynorm2
 * @see mess_matrix_indefinite_dynormf
 */
int mess_matrix_indefinite_dynormf2 ( mess_matrix W, mess_matrix K, double *nrm){
    MSG_FNAME(__func__);
    mess_int_t ret = 0;
    mess_matrix Wt,Kt;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(W); mess_check_nullpointer(K);
    mess_check_real_or_complex(W); mess_check_real_or_complex(K);

    MESS_INIT_MATRICES(&Wt,&Kt);

    if(W->rows < W->cols){
        MSG_INFO("W->rows: " MESS_PRINTF_INT " >= W->cols: " MESS_PRINTF_INT".\n Take Hermitian of W.\n", W->rows, W->cols);
        ret = mess_matrix_ctranspose(W,Wt);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    }else{
        ret = mess_matrix_copy(W,Wt);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    }

    if(K->rows < K->cols){
        MSG_INFO("K->rows: " MESS_PRINTF_INT " >= K->cols: " MESS_PRINTF_INT".\n Take Hermitian of K.\n", K->rows, K->cols);
        ret = mess_matrix_ctranspose(K,Kt);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    }else{
        ret = mess_matrix_copy(K,Kt);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    }


    /*-----------------------------------------------------------------------------
     *  compte eigenvalues of U^TUD
     *-----------------------------------------------------------------------------*/
    mess_matrix U,D,UTU;
    MESS_INIT_MATRICES(&U,&D,&UTU);
    ret = mess_matrix_cat(Wt,Kt,NULL,NULL,MESS_DENSE,U);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);
    ret = mess_matrix_eye(D,U->cols,U->cols,MESS_DENSE);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_eye);
    mess_int_t i = 0;
    for(i=W->cols;i<D->rows;i++){D->values[i+D->ld*i]*=-1;}
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, U, MESS_OP_NONE, U, UTU);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE,UTU,MESS_OP_NONE,D,U);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);

    mess_vector ev;
    ret = mess_vector_init(&ev);                                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
    ret = mess_vector_alloc(ev,U->rows,MESS_COMPLEX);                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
    ret = mess_eigen_eig(U,ev,NULL);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_eigen_eig);
    ret = mess_vector_norm2(ev,nrm);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm);
    (*nrm)*=(*nrm);

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&U,&D,&UTU,&Wt,&Kt);
    MESS_CLEAR_VECTORS(&ev);

    return ret;
}
/**
 * @brief Compute the Frobenius-norm of \f$ WW^T-KK^T \f$,\f$ WW^H-KK^T \f$,\f$ WW^T-KK^H \f$  or \f$ WW^H-KK^H \f$.
 * @param[in] W     input matrix  \f$W\f$
 * @param[in] K     input matrix  \f$K\f$
 * @param[out] nrm  output  Frobenius-norm
 *
 * @return zero on success or a non zero error code
 *
 * The @ref mess_matrix_indefinite_dynormf function computes the square of the Frobenius-norm
 * of \f$ WW^T-KK^T \f$,\f$ WW^H-KK^T \f$,\f$ WW^T-KK^H \f$  or \f$ WW^H-KK^H \f$.
 * It directly calls @ref mess_matrix_indefinite_dynormf2 and takes the root.
 *
 * @see mess_matrix_dynorm2
 * @see mess_matrix_indefinite_dynorm2
 * @see mess_matrix_indefinite_dynormf
 */
int mess_matrix_indefinite_dynormf ( mess_matrix W, mess_matrix K, double *nrm){
    MSG_FNAME(__func__);
    mess_int_t ret=0;
    ret = mess_matrix_indefinite_dynormf2 ( W, K, nrm); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_indefinite_dynormf);
    *nrm = sqrt(*nrm);
    return 0;
}

