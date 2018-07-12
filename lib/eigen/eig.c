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
 * @file lib/eigen/eig.c
 * @brief Wrapper around various eigenvalue computations from @lapack.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>


/**
 * @brief Compute all eigenvalues and eigenvectors.
 * @param[in] A            input matrix
 * @param[out] ew   vector containing eigenvalues
 * @param[out] ev   matrix containing eigenvectors (if not wanted set it to @c NULL)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_eig function solves the eigenvalue problem
 *\f[ Ax = \lambda x. \f]
 * Therefore it determines all eigenvalues \f$ \lambda \f$ stored in \c ew and if wanted all
 * eigenvectors \f$ x \f$  stored in \c ev with xGEEV from @lapack.
 *
 */
int mess_eigen_eig(mess_matrix A, mess_vector ew, mess_matrix ev){
    MSG_FNAME(__func__);
    int cev = 0;    // compute (right) eigenvectors
    mess_int_t one = 1;
    mess_int_t info =0;
    mess_int_t lwork;
    mess_matrix intA;
    int ret = 0;
    mess_int_t n, i, j;
    int cpx = 0 ;


    /*-----------------------------------------------------------------------------
     *  check input and prepare
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_square(A);
    mess_check_nullpointer(ew);

    if ( ev != NULL) {
        cev =1;
        MESS_MATRIX_RESET(ev);
    }
    n = A->rows;
    if ( ew->dim != n) {
        ret = mess_vector_resize(ew, n);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }

    ret = mess_matrix_init(&intA);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_convert(A, intA, MESS_DENSE);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);

    if (MESS_IS_REAL(intA)) {

        double *wr, *wi;
        double dummy = 0.0;
        double *work;
        double ws;
        /*-----------------------------------------------------------------------------
         *  Real Matrix
         *-----------------------------------------------------------------------------*/
        if ( cev == 1) {
            ret = mess_matrix_alloc(ev, A->rows, A->rows, A->rows*A->rows, MESS_DENSE, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        }
        mess_try_alloc(wr, double *, sizeof(double) * A->rows);
        mess_try_alloc(wi, double *, sizeof(double) * A->rows);

        /*-----------------------------------------------------------------------------
         *  Workspace Query
         *-----------------------------------------------------------------------------*/
        lwork = -1;
        if (cev == 1) {
            F77_GLOBAL(dgeev,DGEEV)("N","V", &(A->rows), intA->values, &(A->ld), wr, wi, &dummy,&one, ev->values, &(ev->ld), &ws, &lwork, &info);
        } else {
            F77_GLOBAL(dgeev,DGEEV)("N","N", &(A->rows), intA->values, &(A->ld), wr, wi, &dummy,&one, &dummy, &one, &ws, &lwork, &info);
        }
        lwork=nearbyint(ws+1);
        mess_try_alloc(work, double*, sizeof(double)*lwork);


        /*-----------------------------------------------------------------------------
         *  compute eigenvalues
         *-----------------------------------------------------------------------------*/
        if (cev == 1) {
            F77_GLOBAL(dgeev,DGEEV)("N","V", &(A->rows), intA->values, &(A->ld), wr, wi, &dummy,&one, ev->values, &(ev->ld), work, &lwork, &info);
        } else {
            F77_GLOBAL(dgeev,DGEEV)("N","N", &(A->rows), intA->values, &(A->ld), wr, wi, &dummy,&one, &dummy, &one, work, &lwork, &info);
        }

        if ( info != 0 ){
            MSG_ERROR("dgeev returned with error = " MESS_PRINTF_INT "\n", info);
            return(MESS_ERROR_LAPACK);
        }

        /*-----------------------------------------------------------------------------
         *  postprocessing
         *-----------------------------------------------------------------------------*/
        for ( i = 0; i < n ; i++){
            if ( wi[i] !=0) {
                cpx = 1;
                break;
            }
        }
        if (cpx) {
            if ( cev) {
                ret = mess_matrix_tocomplex(ev);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
            }
            ret = mess_vector_tocomplex(ew);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);

            for (i = 0; i < n; i++){
                if ( wi[i] == 0.0) {
                    ew->values_cpx[i] = wr[i];
                }else {
                    ew->values_cpx[i] = wr[i] + I * wi[i];
                    ew->values_cpx[i+1] = wr[i+1] + I * wi[i+1];
                    if ( cev) {
                        for ( j = 0; j < n; j++){
                            //ev->values_cpx[n*i+j] += ev->values_cpx[n*(i+1)+j] * I;
                            //ev->values_cpx[n*(i+1)+j] = conj(ev->values_cpx[n*i+j]);
                            ev->values_cpx[(ev->ld)*i+j] += ev->values_cpx[(ev->ld)*(i+1)+j] * I;
                            ev->values_cpx[(ev->ld)*(i+1)+j] = conj(ev->values_cpx[(ev->ld)*i+j]);
                        }
                    }
                    i++;
                }
            }
        } else {
            /*-----------------------------------------------------------------------------
             *  only real eigenvalues
             *-----------------------------------------------------------------------------*/
            if ( MESS_IS_REAL(ew) ){
                for (i = 0; i < n ; i++){
                    ew->values[i] = wr[i];
                }
            } else {
                for (i = 0; i < n ; i++){
                    ew->values_cpx[i] = wr[i];
                }
            }
        }

        mess_free( wr );
        mess_free( wi );
        mess_free( work);
    } else {
        double *rwork;
        mess_double_cpx_t * work;
        mess_double_cpx_t *wr;
        mess_double_cpx_t dummy = 1.0;
        mess_double_cpx_t ws = 0;
        /*-----------------------------------------------------------------------------
         *  Complex Matrix
         *-----------------------------------------------------------------------------*/
        if ( cev == 1) {
            ret =   mess_matrix_alloc(ev, A->rows, A->rows, A->rows*A->rows, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        }
        mess_try_alloc(wr, mess_double_cpx_t*, sizeof(mess_double_cpx_t) * A->rows);
        mess_try_alloc(rwork, double *, sizeof(double) * A->rows * 2);

        /*-----------------------------------------------------------------------------
         *  Workspace Query
         *-----------------------------------------------------------------------------*/
        lwork = -1;
        if (cev == 1) {
            F77_GLOBAL(zgeev,ZGEEV)("N","V", &(A->rows), intA->values_cpx, &(A->ld), wr, &dummy,&one, ev->values_cpx, &(ev->ld), &ws, &lwork, rwork,  &info);
        } else {
            F77_GLOBAL(zgeev,ZGEEV)("N","N", &(A->rows), intA->values_cpx, &(A->ld), wr, &dummy,&one, &dummy, &one, &ws, &lwork, rwork, &info);
        }
        lwork=nearbyint(creal(ws)+1);
        mess_try_alloc(work, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*lwork);

        /*-----------------------------------------------------------------------------
         *  compute eigenvalues
         *-----------------------------------------------------------------------------*/
        if (cev == 1) {
            F77_GLOBAL(zgeev,ZGEEV)("N","V", &(A->rows), intA->values_cpx, &(A->ld), wr, &dummy,&one, ev->values_cpx, &(ev->ld), work, &lwork, rwork,  &info);
        } else {
            F77_GLOBAL(zgeev,ZGEEV)("N","N", &(A->rows), intA->values_cpx, &(A->ld), wr, &dummy,&one, &dummy, &one, work, &lwork, rwork, &info);
        }
        if ( info != 0 ){
            MSG_ERROR("zgeev returned with error = " MESS_PRINTF_INT "\n", info);
            return(MESS_ERROR_LAPACK);
        }
        /*-----------------------------------------------------------------------------
         *  postprocessing
         *-----------------------------------------------------------------------------*/
        ret = mess_vector_from_farray(ew, n, NULL, wr);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_from_farray);
        //for (i = 0; i < n ; i++){
        //    ew->values_cpx[i] = wr[i];
        //}

        mess_free( wr );
        mess_free( rwork );
        mess_free( work);
    }

    /*-----------------------------------------------------------------------------
     *  clear if new memory was alloced for intA
     *-----------------------------------------------------------------------------*/
    mess_matrix_clear(&intA);

    return(info) ;
}


/**
 * @brief Compute all eigenvalues and right eigenvectors of the generalized eigenvalue problem.
 * @param[in] A   input matrix
 * @param[in] B   input matrix
 * @param[out] ew   vector containing eigenvalues of the generalized eigenvalue problem (NULL if not wanted)
 * @param[out] alpha       in combination with \c \f$\beta\f$  a rational pair for the eigenvalues (NULL if not wanted)
 * @param[out] beta      in combination with \c \f$\alpha\f$  a rational pair for the eigenvalues (NULL if not wanted)
 * @param[out] ev   matrix containing right eigenvectors ( if not wanted set it to @c NULL)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_eigg function computes all eigenvalues \f$ \lambda \f$ stored in \c ew and if wanted all
 * right eigenvectors \f$ x \f$ stored in \c ev of a given real matrix pair \f$(A,B)\f$, i.e. this function solves
 * \f[ Ax= \lambda B x \f]
 *  with xGGEV from @lapack.
 *
 */
int mess_eigen_eigg(mess_matrix A, mess_matrix B, mess_vector ew, mess_vector alpha, mess_vector beta, mess_matrix ev){
    MSG_FNAME(__func__);
    int cev = 0;    // compute (right) eigenvectors
    double *alphar, *alphai;
    double *betar;
    mess_int_t one = 1;
    mess_int_t info =0;
    mess_int_t lwork;
    mess_matrix intA = NULL, intB = NULL;
    int ret = 0;
    mess_int_t n, i, j;
    int cpx = 0 ;
    int returnalpha = 0;
    int returnbeta = 0 ;
    int returnew   = 0;


    /*-----------------------------------------------------------------------------
     *  check input and prepare
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_square(A);
    mess_check_square(B);
    mess_check_nullpointer(ew);
    mess_check_same_size(A,B);

    n= A->rows;
    if ( ev != NULL) {
        cev =1;
        MESS_MATRIX_RESET(ev);
    }
    if ( alpha != NULL) {
        ret = mess_vector_resize(alpha, A->rows); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        returnalpha =  1;
    }
    if ( beta != NULL){
        ret = mess_vector_resize(beta, A->rows); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        returnbeta = 1;
    }
    if ( ew  != NULL ) {
        ret = mess_vector_resize(ew, n); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        returnew = 1;
    }


    ret = mess_matrix_init(&intA);                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&intB);                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_convert(A, intA, MESS_DENSE);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    ret = mess_matrix_convert(B, intB, MESS_DENSE);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);

    if (MESS_IS_REAL(intA) && MESS_IS_REAL(intB)) {
        double ws = 0;
        double *work;
        double dummy = 0.0;
        /*-----------------------------------------------------------------------------
         *  REAL CASE
         *-----------------------------------------------------------------------------*/
        if ( cev == 1) {
            ret = mess_matrix_alloc(ev, A->rows, A->rows, A->rows*A->rows, MESS_DENSE, MESS_REAL);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        }

        mess_try_alloc(alphar, double *, sizeof(double) * A->rows);
        mess_try_alloc(alphai, double *, sizeof(double) * A->rows);
        mess_try_alloc(betar, double *, sizeof(double) * A->rows);

        /*-----------------------------------------------------------------------------
         *  Workspace Query
         *-----------------------------------------------------------------------------*/
        lwork = -1;
        if (cev == 1) {
            F77_GLOBAL(dggev,DGGEV)("N","V", &(A->rows), intA->values, &(A->ld), intB->values, &(intB->ld), alphar, alphai,betar, &dummy,&one, ev->values, &(ev->ld), &ws, &lwork, &info);
        } else {
            F77_GLOBAL(dggev,DGGEV)("N","N", &(A->rows), intA->values, &(A->ld), intB->values, &(intB->ld), alphar, alphai,betar,  &dummy,&one, &dummy, &one, &ws, &lwork, &info);
        }

        lwork=nearbyint(ws+1);
        mess_try_alloc(work, double*, sizeof(double)*lwork);
        /*-----------------------------------------------------------------------------
         *  compute eigenvalues
         *-----------------------------------------------------------------------------*/
        if (cev == 1) {
            F77_GLOBAL(dggev,DGGEV)("N","V", &(A->rows), intA->values, &(A->ld), intB->values, &(intB->ld), alphar, alphai,betar, &dummy,&one, ev->values, &(ev->ld), work, &lwork, &info);
        } else {
            F77_GLOBAL(dggev,DGGEV)("N","N", &(A->rows), intA->values, &(A->ld), intB->values, &(intB->ld), alphar, alphai,betar,  &dummy,&one, &dummy, &one, work, &lwork, &info);
        }

        if ( info != 0 ){
            mess_free(alphai);
            mess_free(alphar);
            mess_free(betar);
            MSG_ERROR("dggev returned with error = " MESS_PRINTF_INT "\n", info);
            return(MESS_ERROR_LAPACK);
        }

        /*-----------------------------------------------------------------------------
         *  postprocessing
         *-----------------------------------------------------------------------------*/
        for ( i = 0; i < n ; i++){
            if ( alphai[i] !=0) cpx = 1;
        }
        if (cpx) {
            if ( cev) {
                ret = mess_matrix_tocomplex(ev);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
            }
            if ( returnew ) {
                ret = mess_vector_tocomplex(ew);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            }
            if (returnalpha) {
                ret = mess_vector_tocomplex(alpha); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            }
            if ( returnbeta ) {
                ret = mess_vector_toreal_nowarn(beta); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
            }


            int skipnext =0;
            for (i = 0; i < n; i++){
                if ( returnew ){
                    if ( alphai[i] == 0.0) {
                        ew->values_cpx[i] = (alphar[i])/betar[i];
                    } else if ( skipnext == 0)  {
                        ew->values_cpx[i] = (alphar[i] + I * alphai[i])/betar[i];
                        ew->values_cpx[i+1] =  (alphar[i+1] + I * alphai[i+1])/betar[i+1];
                        if ( cev) {
                            for ( j = 0; j < n; j++){
                                ev->values_cpx[(ev->ld)*i+j] += ev->values_cpx[(ev->ld)*(i+1)+j] * I;
                                ev->values_cpx[(ev->ld)*(i+1)+j] = conj(ev->values_cpx[(ev->ld)*i+j]);
                            }
                        }
                        skipnext = 1;
                    } else {
                        skipnext = 0;
                    }
                }
                if (returnalpha) alpha->values_cpx[i]=alphar[i]+alphai[i]*I;
                if (returnbeta)  beta->values[i] = betar[i];
            }
        } else {
            /*-----------------------------------------------------------------------------
             *  only real eigenvalues
             *-----------------------------------------------------------------------------*/
            if ( MESS_IS_REAL(ew) ){
                for (i = 0; i < n ; i++){
                    ew->values[i] = alphar[i]/betar[i];
                }
            } else {
                for (i = 0; i < n ; i++){
                    ew->values_cpx[i] = alphar[i]/betar[i];
                }
            }
            if ( returnalpha ){
                ret = mess_vector_toreal_nowarn(alpha); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
                for (i = 0; i < n ; i++){
                    alpha->values[i] = alphar[i];
                }

            }
            if ( returnbeta ){
                ret = mess_vector_toreal_nowarn(beta); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
                for (i = 0; i < n ; i++){
                    beta->values[i] = betar[i];
                }

            }

        }

        mess_free( alphar );
        mess_free( alphai );
        mess_free( betar);
        mess_free( work);
    } else {
        /*-----------------------------------------------------------------------------
         *  Complex Case
         *-----------------------------------------------------------------------------*/
        mess_double_cpx_t *work, *alphac, *betac;
        mess_double_cpx_t dummy = 1.0;
        double *rwork;
        mess_double_cpx_t ws = 0;

        ret = mess_matrix_tocomplex(intA);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        ret = mess_matrix_tocomplex(intB);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        if ( cev == 1) {
            ret = mess_matrix_alloc(ev, A->rows, A->rows, A->rows*A->rows, MESS_DENSE, MESS_COMPLEX);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        }

        mess_try_alloc(alphac, mess_double_cpx_t *, sizeof(mess_double_cpx_t ) * A->rows);
        mess_try_alloc(betac, mess_double_cpx_t*, sizeof(mess_double_cpx_t) * A->rows);
        mess_try_alloc(rwork, double *, sizeof(double) * A->rows * 8);

        /*-----------------------------------------------------------------------------
         *  Workspace Query
         *-----------------------------------------------------------------------------*/
        lwork = -1;
        if (cev == 1) {
            F77_GLOBAL(zggev,ZGGEV)("N","V", &(A->rows), intA->values_cpx, &(A->ld), intB->values_cpx, &(intB->ld), alphac, betac, &dummy,&one, ev->values_cpx, &(ev->ld), &ws, &lwork, rwork,  &info);
        } else {
            F77_GLOBAL(zggev,ZGGEV)("N","N", &(A->rows), intA->values_cpx, &(A->ld), intB->values_cpx, &(intB->ld), alphac, betac, &dummy,&one, &dummy, &one, &ws, &lwork, rwork, &info);
        }
        lwork=nearbyint(creal(ws)+1);
        mess_try_alloc(work, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*lwork);

        /*-----------------------------------------------------------------------------
         *  compute eigenvalues
         *-----------------------------------------------------------------------------*/
        if (cev == 1) {
            F77_GLOBAL(zggev,ZGGEV)("N","V", &(A->rows), intA->values_cpx, &(A->ld), intB->values_cpx, &(intB->ld), alphac, betac, &dummy,&one, ev->values_cpx, &(ev->ld), work, &lwork, rwork,  &info);
        } else {
            F77_GLOBAL(zggev,ZGGEV)("N","N", &(A->rows), intA->values_cpx, &(A->ld), intB->values_cpx, &(intB->ld), alphac, betac, &dummy,&one, &dummy, &one, work, &lwork, rwork, &info);
        }
        if ( info != 0 ){
            MSG_ERROR("dggev returned with error = " MESS_PRINTF_INT "\n", info);
            return(MESS_ERROR_LAPACK);
        }

        /*-----------------------------------------------------------------------------
         *  postprocessing
         *-----------------------------------------------------------------------------*/
        if ( returnew ) {
            ret = mess_vector_tocomplex(ew);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        }
        if (returnalpha) {
            ret = mess_vector_tocomplex(alpha); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        }
        if ( returnbeta ) {
            ret = mess_vector_toreal_nowarn(beta); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        }

        for (i = 0; i < n ; i++){
            ew->values_cpx[i] = alphac[i]/betac[i];
        }
        if ( returnalpha ){
            for (i = 0; i < n ; i++){
                alpha->values_cpx[i] = alphac[i];
            }
        }
        if ( returnbeta ){
            for (i = 0; i < n ; i++){
                beta->values_cpx[i] = betac[i];
            }
        }
        mess_free( alphac );
        mess_free( betac);
        mess_free( work);
        mess_free( rwork);
    }
    mess_matrix_clear(&intA);
    mess_matrix_clear(&intB);

    return(info);
}
