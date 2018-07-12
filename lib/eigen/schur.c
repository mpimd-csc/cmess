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
 * @file lib/eigen/schur.c
 * @brief Compute the Schur decomposition of a matrix.
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

#define COMPLEX_ZGEES_IF_REAL 1


/**
 * @brief Wrapper around both complex Schur decomposition routines.
 * @param[in]  A   input matrix to decompose
 * @param[out] T  Schur form of \f$ A \f$
 * @param[out] U  orthogonal transformation (@c NULL if not needed)
 * @param[out] EV vector containing eigenvalues of \f$ A \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_schur_complex function is a wrapper around the
 * two complex Schur decomposition subroutines \ref mess_eigen_schur_complex_zgees
 * and \ref mess_eigen_schur_complex_post. \n
 *
 */
int mess_eigen_schur_complex ( mess_matrix A, mess_matrix T, mess_matrix U, mess_vector EV ){
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(A);
    mess_check_nullpointer(T);

    if ( MESS_IS_COMPLEX(A)) {
        ret = mess_eigen_schur_complex_zgees(A,T,U,EV);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_schur_complex_zgees);
    } else {
        if ( COMPLEX_ZGEES_IF_REAL ){
            ret = mess_eigen_schur_complex_zgees(A,T,U,EV);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_schur_complex_zgees);
        } else {
            ret = mess_eigen_schur_complex_post(A,T,U,EV);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_schur_complex_post);
        }
    }
    return 0;
}       /* -----  end of function mess_eigen_schur_complex  ----- */

/**
 * @brief Compute the complex Schur decomposition of a matrix.
 * @param[in]  A         input matrix to decompose
 * @param[out] T    Schur form of \f$ A \f$
 * @param[out] U    orthogonal transformation (NULL if not needed)
 * @param[out] EV        vector containing eigenvalues of \f$ A \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_schur_complex_zgees function computes the complex Schur decomposition
 * of a \f$ A \f$:
 * \f[
 * T=U^H A U
 * \f]
 * using @lapack's ZGEES subroutine.
 *
 */
int mess_eigen_schur_complex_zgees ( mess_matrix A, mess_matrix T, mess_matrix U, mess_vector EV )
{
    MSG_FNAME(__func__);
    int withU=0;
    int withEV=0;
    int ret = 0;
    mess_int_t lwork=0, info=0, lda=0, n=0, sdim=0, ldu;
    double *rwork=NULL;
    mess_double_cpx_t *work=NULL, *W=NULL;
    mess_double_cpx_t ws;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_square(A);
    if ( !MESS_IS_DENSE(A)){
        MSG_ERROR("only dense matrices allowed.\n");
        return MESS_ERROR_STORAGETYPE;
    }
    mess_check_nullpointer(T);
    if ( U != NULL) {
        MSG_INFO("computing orthogonal transformation\n");
        withU=1;
    }
    if ( EV!=NULL){
        MSG_INFO("with seperate eigenvalues.\n");
        withEV=1;
    }

    ret = mess_matrix_copy(A,T);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_tocomplex(T);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);


    ldu = 1;
    if ( withU==1){
        ret=mess_matrix_alloc(U, A->rows, A->rows, A->rows*A->rows, MESS_DENSE, MESS_COMPLEX);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ldu = U->ld;
    }

    n=A->rows;
    lda=T->ld;
    mess_try_alloc(W,mess_double_cpx_t *, sizeof(mess_double_cpx_t)*n);


    /*-----------------------------------------------------------------------------
     *   work space query
     *-----------------------------------------------------------------------------*/
    lwork = -1;
    ws = 0;
    if ( withU )
        zgees_("V","N",NULL,&n,T->values_cpx, &lda, &sdim,W,U->values_cpx,&ldu,&ws, &lwork,rwork,NULL, &info);
    else
        zgees_("N","N",NULL,&n,T->values_cpx, &lda, &sdim,W,NULL,&ldu,&ws, &lwork,rwork,NULL, &info);

    lwork=nearbyint(creal(ws)+1);
    mess_try_alloc(work, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*lwork);
    mess_try_alloc(rwork, double*, sizeof(double) *n);

    /*-----------------------------------------------------------------------------
     *  call LAPACK
     *-----------------------------------------------------------------------------*/
    if ( withU )
        zgees_("V","N",NULL,&n,T->values_cpx, &lda, &sdim,W,U->values_cpx,&ldu,work, &lwork,rwork,NULL, &info);
    else
        zgees_("N","N",NULL,&n,T->values_cpx, &lda, &sdim,W,NULL,&ldu,work, &lwork,rwork,NULL, &info);

    if ( info != 0){
        MSG_ERROR("ZGEES returned with error: " MESS_PRINTF_INT "\n", info);
        ret = MESS_ERROR_LAPACK;
    } else {
        if (withEV){
            mess_int_t i;
            ret = mess_vector_tocomplex(EV);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            ret = mess_vector_resize(EV, n);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_resize);
            for (i=0; i < n ; i++){
                EV->values_cpx[i]=W[i];
            }
        }
        ret =0;
    }
    /*-----------------------------------------------------------------------------
     * clean up
     *-----------------------------------------------------------------------------*/
    mess_free(W);
    mess_free(work);
    mess_free(rwork);

    return ret;
}       /* -----  end of function mess_eigen_schur_complex  ----- */


/*
 * computes U * Ue^T
 * where Ue is [ ue11 ue12;
 *               ue21 ue22]
 *
 *
 */
static int applyUe_right(mess_int_t n, mess_int_t ld, mess_double_cpx_t *U, mess_double_cpx_t ue11, mess_double_cpx_t ue21, mess_double_cpx_t ue12,
        mess_double_cpx_t ue22){
    mess_int_t i;
    mess_double_cpx_t t1,t2;

    for ( i = 0 ; i < n ; i++) {
        t1 = U[i];
        t2 = U[i+ld];
        U[i] = ue11 * t1 + ue12 *t2;
        U[i+ld] = ue21 * t1 + ue22 *t2;
    }
    return 0;
}

/*
 * computes Ue * U
 * where Ue is [ ue11 ue12;
 *       ue21 ue22 ]
 *
 */
static int applyUe_left(mess_int_t n,mess_int_t ld,  mess_double_cpx_t *U, mess_double_cpx_t ue11, mess_double_cpx_t ue12, mess_double_cpx_t ue21,
        mess_double_cpx_t ue22){
    mess_int_t i;
    mess_double_cpx_t t1, t2;
    for ( i = 0 ; i< n ; i++) {
        t1=U[i*ld];
        t2=U[i*ld+1];
        U[i*ld]= ue11*t1 + ue12*t2;
        U[i*ld+1]= ue21 *t1 + ue22*t2;
    }
    return 0;

}

/**
 * @brief Compute the complex Schur decomposition of a matrix using post processing.
 * @param[in] A           input matrix to decompose
 * @param[out] T    Schur form of \f$ A \f$
 * @param[out] U    orthogonal transformation (NULL if not needed)
 * @param[out] EV   vector containing eigenvalues of \f$ A \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_schur_complex_post function computes the complex Schur decomposition
 * of \f$ A \f$:
 * \f[ T= U^H A U \f]
 * using @lapack's DGEES subroutine and a postprocessing step.
 *
 */
int mess_eigen_schur_complex_post( mess_matrix A, mess_matrix T, mess_matrix U, mess_vector EV )
{
    MSG_FNAME(__func__);
    int withU=0;
    int withEV=0;
    int ret = 0;
    mess_int_t lwork=0, info=0, lda=0, n=0, sdim=0, ldu= 1;
    double *work=NULL, *WR=NULL, *WI=NULL;
    mess_int_t k;
    mess_double_cpx_t  b,c,x,y;
    double ws=0;
    // double Ue[4];


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_square(A);
    mess_check_dense(A);
    mess_check_real(A);
    mess_check_nullpointer(T);
    if ( U != NULL) {
        MSG_INFO("computing orthogonal transformation\n");
        withU=1;
    }
    if ( EV!=NULL){
        MSG_INFO("with seperate eigenvalues.\n");
        withEV=1;
    }

    ret = mess_matrix_copy(A,T);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);


    if ( withU==1){
        ret=mess_matrix_alloc(U, A->rows, A->rows, A->rows*A->rows, MESS_DENSE, MESS_REAL);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ldu = U->ld;
    }

    n=A->rows;
    lda=T->ld;
    mess_try_alloc(WR,double *, sizeof(double)*n);
    mess_try_alloc(WI,double *, sizeof(double)*n);

    /*-----------------------------------------------------------------------------
     *  workspace query
     *-----------------------------------------------------------------------------*/

    lwork=-1;
    if (withU)
        dgees_("V","N",NULL,&n,T->values, &lda, &sdim,WR,WI,U->values,&ldu,&ws, &lwork,NULL, &info);
    else
        dgees_("N","N",NULL,&n,T->values, &lda, &sdim,WR,WI,NULL,&ldu,&ws, &lwork,NULL, &info);

    lwork = nearbyint(ws+1);
    mess_try_alloc(work, double *, sizeof(double )*lwork);

    /*-----------------------------------------------------------------------------
     *  call LAPACK
     *-----------------------------------------------------------------------------*/
    if (withU)
        dgees_("V","N",NULL,&n,T->values, &lda, &sdim,WR,WI,U->values,&ldu,work, &lwork,NULL, &info);
    else
        dgees_("N","N",NULL,&n,T->values, &lda, &sdim,WR,WI,NULL,&ldu,work, &lwork,NULL, &info);

    if ( info != 0){
        mess_free(WR);
        mess_free(WI);
        mess_free(work);
        MSG_ERROR("DGEES returned with error: " MESS_PRINTF_INT "\n", info);
        ret = MESS_ERROR_LAPACK;
        return ret;
    }
    /*-----------------------------------------------------------------------------
     *  post processing
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_tocomplex(T);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
    if (withU){
        ret = mess_matrix_tocomplex(U);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
    }

    for ( k =0 ; k < n-1; k++){
        if ( WI[k] != 0.0 ) {
            // a = T->values_cpx[k*n+k];
            b = T->values_cpx[(k+1)*T->ld+k];
            c = T->values_cpx[k*T->ld+k+1];
            x = 1.0/ csqrt( 1 + cabs(b/c));
            y = (csqrt(b*c)/c)*x;
            applyUe_right(n,T->ld, &T->values_cpx[T->ld*k], y,x,conj(x), -conj(y));
            applyUe_left(n, T->ld, &T->values_cpx[k], conj(y),conj(x),x, -y);
            if ( withU ) {
                applyUe_right(n,U->ld, &U->values_cpx[U->ld*k], y,x,conj(x), -conj(y));
            }
            k++;
        }
    }


    if (withEV){
        mess_int_t i;
        ret = mess_vector_tocomplex(EV);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_resize(EV, n);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_resize);
        for (i=0; i < n ; i++){
            EV->values_cpx[i]=WR[i]+WI[i]*I;
        }
    }

    /*-----------------------------------------------------------------------------
     * clean up
     *-----------------------------------------------------------------------------*/
    mess_free(WR);
    mess_free(WI);
    mess_free(work);

    return ret;
}       /* -----  end of function mess_eigen_schur_complex  ----- */




/**
 * @brief Compute the real Schur decomposition of a matrix.
 * @param[in] A     input matrix to decompose
 * @param[out] T    Schur form of \f$ A \f$
 * @param[out] U    orthogonal transformation (@c NULL if not needed)
 * @param[out] EV   vector containing eigenvalues of \f$ A \f$ (@c NULL if not needed)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_schur function computes the real Schur decomposition
 * of \f$ A \f$
 *
 * \f[ T= U^T A U \f]
 * \f[ I_n= U^T U \f]
 *
 * using @lapack's DGEES subroutine.
 *
 */
int mess_eigen_schur( mess_matrix A, mess_matrix T, mess_matrix U, mess_vector EV )
{
    MSG_FNAME(__func__);
    int withU=0;
    int withEV=0;
    int ret = 0;
    mess_int_t lwork=0, info=0, lda=0, n=0, sdim=0,ldu =1;
    double *work=NULL, *WR=NULL, *WI=NULL;
    double ws=0;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_square(A);
    mess_check_dense(A);
    mess_check_real(A);
    mess_check_nullpointer(T);
    if ( U != NULL) {
        MSG_INFO("computing orthogonal transformation\n");
        withU=1;
    }
    if ( EV!=NULL){
        MSG_INFO("with seperate eigenvalues.\n");
        withEV=1;
    }

    /*-----------------------------------------------------------------------------
     *  prepare lapack call
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_copy(A,T); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);

    if ( withU==1){
        ret=mess_matrix_alloc(U, A->rows, A->rows, A->rows*A->rows, MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ldu = U->ld;
    }

    n=A->rows;
    lda=T->ld;
    mess_try_alloc(WR,double *, sizeof(double)*n);
    mess_try_alloc(WI,double *, sizeof(double)*n);

    /*-----------------------------------------------------------------------------
     *  workspace query
     *-----------------------------------------------------------------------------*/
    lwork=-1;
    if (withU)
        dgees_("V","N",NULL,&n,T->values, &lda, &sdim,WR,WI,U->values,&ldu,&ws, &lwork,NULL, &info);
    else
        dgees_("N","N",NULL,&n,T->values, &lda, &sdim,WR,WI,NULL,&ldu,&ws, &lwork,NULL, &info);

    lwork = nearbyint(ws+1);
    mess_try_alloc(work, double *, sizeof(double )*lwork);

    /*-----------------------------------------------------------------------------
     *  call LAPACK
     *-----------------------------------------------------------------------------*/
    if (withU)
        dgees_("V","N",NULL,&n,T->values, &lda, &sdim,WR,WI,U->values,&ldu,work, &lwork,NULL, &info);
    else
        dgees_("N","N",NULL,&n,T->values, &lda, &sdim,WR,WI,NULL,&ldu,work, &lwork,NULL, &info);

    if ( info != 0){
        MSG_ERROR("DGEES returned with error: " MESS_PRINTF_INT "\n", info);
        ret = MESS_ERROR_LAPACK;
    } else {
        if (withEV){
            mess_int_t i;
            ret = mess_vector_tocomplex(EV);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            ret = mess_vector_resize(EV, n);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_resize);
            for (i=0; i < n ; i++){
                EV->values_cpx[i]=WR[i]+WI[i]*I;
            }
        }
    }


    /*-----------------------------------------------------------------------------
     * clean up
     *-----------------------------------------------------------------------------*/
    mess_free(WR);
    mess_free(WI);
    mess_free(work);

    return ret;
}       /* -----  end of function mess_eigen_schur  ----- */



/**
 * @brief Compute the real generalized Schur decomposition of a matrix pair \f$ (A,B) \f$.
 * @param[in] A      input matrix to decompose
 * @param[in] B          input matrix to decompose
 * @param[out] S    Schur form of \f$ A \f$
 * @param[out] T    Schur form of \f$ B \f$
 * @param[in,out] U     left orthogonal transformation
 * @param[in,out] V     right orthogonal transformation
 * @param[in,out] EVA   generalized eigenvalues of \f$ A \f$
 * @param[in,out] EVB   generalized eigenvalues of \f$ B \f$
 * @param[in,out] EV      vector containing eigenvalues of the generalized eigenvalue problem
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_gschur function computes the real generalized Schur decomposition of \f$ A \f$ and \f$ B \f$:
 * \f[
 * \begin{array}{ccc}
 * S &=& U^T A V, \\
 * T &=& U^T B V, \\
 * \end{array}
 * \f]
 * where \f$ S \f$ is an upper quasi-triangular matrix and \f$ T \f$ is an upper triangular matrix
 * and \c EVA and \c EVB using @lapack's DGGES subroutine. \n
 * Additionally this function computes the eigenvalues \c EV of the generalized eigenvalue problem :
 * \f[ A x = \lambda B x. \f]
 *
 */
int mess_eigen_gschur( mess_matrix A, mess_matrix B, mess_matrix S , mess_matrix T, mess_matrix U, mess_matrix V, mess_vector EVA, mess_vector EVB, mess_vector EV )
{
    MSG_FNAME(__func__);
    int withU=0, withV=0;
    int withEVA=0, withEVB=0, withEV=0;
    int ret = 0;
    mess_int_t lwork=0, info=0, lda=0, n=0, sdim=0, ldb = 1 , ldu = 1, ldv=1;
    double *work=NULL;
    double *alphar, *alphai, *beta;
    char JOBVSL[10], JOBVSR[10];
    double *uval=NULL, *vval=NULL;
    double ws;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_square(A);
    mess_check_square(B);
    mess_check_same_size(A,B);
    mess_check_real(A);
    mess_check_real(B);
    mess_check_nullpointer(S);
    mess_check_nullpointer(T);

    if ( U != NULL) {
        withU=1;
    }
    if ( V != NULL) {
        withV=1;
    }

    if ( EVA!=NULL){
        withEVA=1;
    }
    if ( EVB!=NULL){
        withEVB=1;
    }
    if ( EV!=NULL){
        withEV=1;
    }

    ret = mess_matrix_convert(A,S,MESS_DENSE); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    ret = mess_matrix_convert(B,T,MESS_DENSE); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);

    if ( withU==1){
        ret=mess_matrix_alloc(U, A->rows, A->rows, A->rows*A->rows, MESS_DENSE, MESS_REAL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        snprintf(JOBVSL,9,"V");
        uval = U->values;
        ldu = U->ld;
    } else {
        snprintf(JOBVSL,9,"N");
    }
    if ( withV==1){
        ret=mess_matrix_alloc(V, A->rows, A->rows, A->rows*A->rows, MESS_DENSE, MESS_REAL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        snprintf(JOBVSR,9,"V");
        vval= V->values;
        ldv = V->ld;
    } else snprintf(JOBVSR, 9, "N");

    n=A->rows;
    lda=S->ld;
    ldb=T->ld;
    mess_try_alloc(alphar,double *, sizeof(double)*n);
    mess_try_alloc(alphai,double *, sizeof(double)*n);
    mess_try_alloc(beta,double *, sizeof(double)*n);


    /*-----------------------------------------------------------------------------
     *  workspace query
     *-----------------------------------------------------------------------------*/
    ws = 0;
    lwork = -1;
    dgges_(JOBVSL, JOBVSR,"N",NULL,&n, S->values, &lda, T->values, &ldb, &sdim, alphar, alphai, beta, uval, &ldu, vval, &ldv,&ws, &lwork, NULL, &info);

    lwork=nearbyint(ws+1);
    mess_try_alloc(work, double *, sizeof(double )*lwork);

    /*-----------------------------------------------------------------------------
     *  call LAPACK
     *-----------------------------------------------------------------------------*/
    dgges_(JOBVSL, JOBVSR,"N",NULL,&n, S->values, &lda, T->values, &ldb, &sdim, alphar, alphai, beta, uval, &ldu, vval, &ldv,work, &lwork, NULL, &info);
    if ( info != 0){
        MSG_ERROR("DGEES returned with error: " MESS_PRINTF_INT "\n", info);
        ret = MESS_ERROR_LAPACK;
    } else {
        if (withEVA){
            mess_int_t i;
            ret = mess_vector_tocomplex(EVA); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            ret = mess_vector_resize(EVA, n); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_resize);
            for (i=0; i < n ; i++){
                EVA->values_cpx[i]=alphar[i]+alphai[i]*I;
            }
        }
        if(withEVB) {
            mess_int_t i;
            ret= mess_vector_toreal(EVB); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
            ret = mess_vector_resize(EVB, n); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
            for (i=0; i<n; i++){
                EVB->values[i] = beta[i];
            }
        }
        if (withEV){
            mess_int_t i;
            ret = mess_vector_tocomplex(EV); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            ret = mess_vector_resize(EV, n); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_resize);
            for (i=0; i < n ; i++){
                if (alphai[i] != 0 ) {
                    EV->values_cpx[i]=(alphar[i]+alphai[i]*I)/beta[i];
                    EV->values_cpx[i+1]= conj ( EV->values_cpx[i] );
                    i = i +1;
                } else {
                    EV->values_cpx[i]=(alphar[i])/beta[i];
                }
            }
        }

        ret =0;
    }


    /*-----------------------------------------------------------------------------
     * clean up
     *-----------------------------------------------------------------------------*/
    mess_free(alphar);
    mess_free(alphai);
    mess_free(beta);
    mess_free(work);

    return ret;
}

/**
 * @brief Compute the generalized Schur decomposition of a matrix pair \f$ (A,B) \f$.
 * @param[in] A   input matrix to decompose
 * @param[in] B     input matrix to decompose
 * @param[out] S    Schur form of \f$ A \f$
 * @param[out] T    Schur form of \f$ B \f$
 * @param[out] U    left orthogonal transformation
 * @param[out] V    right orthogonal transformation
 * @param[out] EVA  generalized eigenvalues of \f$ A \f$
 * @param[out] EVB  generalized eigenvalues of \f$ B  \f$
 * @param[out] EV    vector containing eigenvalues of the generalized eigenvalue problem
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref mess_eigen_gschur_complex function computes the generalized Schur decomposition of \f$ A \f$ and \f$ B \f$:
 * \f[
 * \begin{array}{ccc}
 * S &=& U^T A V, \\
 * T &=& U^T B V, \\
 * \end{array}
 * \f]
 * where \f$ S \f$ and \f$ T \f$ are upper triangular matrices and \c EVA and \c EVB using
 * @lapack's ZGGES subroutine. \n
 * Additionally this function computes the eigenvalues \c EV of the generalized eigenvalue problem :
 * \f[ A x = \lambda B x. \f]
 *
 */
int mess_eigen_gschur_complex( mess_matrix A, mess_matrix B, mess_matrix S , mess_matrix T, mess_matrix U, mess_matrix V, mess_vector EVA, mess_vector EVB, mess_vector EV  )
{
    MSG_FNAME(__func__);
    int withU=0, withV=0;
    int withEVA=0, withEVB=0, withEV=0;
    int ret = 0;
    mess_int_t lwork=0, info=0, lda=0, n=0, sdim=0, ldb =1, ldu=1, ldv=1;
    mess_double_cpx_t *work=NULL;
    mess_double_cpx_t *alpha, *beta;
    double *rwork;
    char JOBVSL[10], JOBVSR[10];
    mess_double_cpx_t *uval=NULL, *vval=NULL;
    mess_double_cpx_t ws;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_square(A);
    mess_check_square(B);
    mess_check_same_size(A,B);
    mess_check_nullpointer(S);
    mess_check_nullpointer(T);

    if ( U != NULL) {
        MSG_INFO("computing orthogonal transformation\n"); withU=1;
    }
    if ( V != NULL) {
        MSG_INFO("computing orthogonal transformation\n"); withV=1;
    }

    if ( EVA!=NULL){
        MSG_INFO("with seperate eigenvalues.\n"); withEVA=1;
    }
    if ( EVB!=NULL){
        MSG_INFO("with seperate eigenvalues.\n"); withEVB=1;
    }
    if ( EV!=NULL){
        MSG_INFO("with combined eigenvalues.\n"); withEV=1;
    }

    ret = mess_matrix_convert(A,S,MESS_DENSE); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    ret = mess_matrix_convert(B,T,MESS_DENSE); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    ret = mess_matrix_tocomplex(S); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
    ret = mess_matrix_tocomplex(T); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);

    if ( withU==1){
        ret=mess_matrix_alloc(U, A->rows, A->rows, A->rows*A->rows, MESS_DENSE, MESS_COMPLEX);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        snprintf(JOBVSL,9,"V");
        uval = U->values_cpx;
        ldu = U->ld;
    } else {
        snprintf(JOBVSL,9,"N");
    }
    if ( withV==1){
        ret=mess_matrix_alloc(V, A->rows, A->rows, A->rows*A->rows, MESS_DENSE, MESS_COMPLEX);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        snprintf(JOBVSR,9,"V");
        vval= V->values_cpx;
        ldv=V->ld;
    } else snprintf(JOBVSR, 9, "N");

    n=A->rows;
    lda=S->ld;
    ldb=T->ld;
    mess_try_alloc(alpha,mess_double_cpx_t *, sizeof(mess_double_cpx_t )*n);
    mess_try_alloc(beta,mess_double_cpx_t *, sizeof(mess_double_cpx_t)*n);
    mess_try_alloc(rwork, double *, sizeof(double) * 8 * n);

    /*-----------------------------------------------------------------------------
     *  work space query
     *-----------------------------------------------------------------------------*/
    lwork = -1;
    ws = 0.0;
    zgges_(JOBVSL, JOBVSR,"N",NULL,&n, S->values_cpx, &lda, T->values_cpx, &ldb, &sdim, alpha, beta, uval, &ldu, vval, &ldv,&ws, &lwork,rwork, NULL, &info);

    lwork=nearbyint(creal(ws)+1);
    mess_try_alloc(work, mess_double_cpx_t *, sizeof(mess_double_cpx_t )*lwork);

    /*-----------------------------------------------------------------------------
     *  call LAPACK
     *-----------------------------------------------------------------------------*/
    zgges_(JOBVSL, JOBVSR,"N",NULL,&n, S->values_cpx, &lda, T->values_cpx, &ldb, &sdim, alpha, beta, uval, &ldu, vval, &ldv,work, &lwork,rwork, NULL, &info);
    if ( info != 0){
        MSG_ERROR("ZGEES returned with error: " MESS_PRINTF_INT "\n", info);
        ret = MESS_ERROR_LAPACK;
    } else {
        if (withEVA){
            mess_int_t i;
            ret = mess_vector_tocomplex(EVA); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            ret = mess_vector_resize(EVA, n); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_resize);
            for (i=0; i < n ; i++){
                EVA->values_cpx[i]=alpha[i];
            }
        }
        if(withEVB) {
            mess_int_t i;
            ret= mess_vector_tocomplex(EVB); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            ret = mess_vector_resize(EVB, n); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
            for (i=0; i<n; i++){
                EVB->values_cpx[i] = beta[i];
            }
        }
        if ( withEV) {
            mess_int_t i;
            ret = mess_vector_tocomplex(EV); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            ret = mess_vector_resize(EV, n); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_resize);
            for (i=0; i < n ; i++){
                EV->values_cpx[i]=alpha[i]/beta[i];
            }

        }
        ret =0;
    }

    /*-----------------------------------------------------------------------------
     * clean up
     *-----------------------------------------------------------------------------*/
    mess_free(alpha);
    mess_free(rwork);
    mess_free(beta);
    mess_free(work);

    return ret;
}       /* -----  end of function mess_eigen_schur_complex  ----- */

