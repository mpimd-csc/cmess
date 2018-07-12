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
 * @file lib/eigen/sign.c
 * @brief Implementation of a basic sign function iteration.
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

//double sign function routines
void F77_GLOBAL(dgesgnnone,DGESGNNONE)(mess_int_t *N,double * A,mess_int_t *LDA, mess_int_t *MAXIT,double * TOL,double *WORK,mess_int_t *IWORK,double *WORK2,mess_int_t *LWORK2,mess_int_t * VERBOSE,mess_int_t *INFO);
void F77_GLOBAL(dgesgnfro,DGESGNFRO)(mess_int_t *N,double * A,mess_int_t *LDA, mess_int_t *MAXIT,double * TOL,double *WORK,mess_int_t *IWORK,double *WORK2,mess_int_t *LWORK2,mess_int_t *VERBOSE,mess_int_t *INFO);
void F77_GLOBAL(dgesgndet,DGESGNDET)(mess_int_t *N,double * A,mess_int_t *LDA, mess_int_t *MAXIT,double * TOL,double *WORK,mess_int_t *IWORK,double *WORK2,mess_int_t *LWORK2,mess_int_t *VERBOSE,mess_int_t *INFO);

//complex sign function routines
void F77_GLOBAL(zgesgnnone,ZGESGNNONE)(mess_int_t *N,mess_double_cpx_t * A,mess_int_t *LDA, mess_int_t *MAXIT,double * TOL,mess_double_cpx_t *WORK,mess_int_t *IWORK,
        mess_double_cpx_t *WORK2, mess_int_t *LWORK2, mess_int_t *VERBOSE,mess_int_t *INFO);
void F77_GLOBAL(zgesgndet,ZGESGNDET)(mess_int_t *N,mess_double_cpx_t * A,mess_int_t *LDA, mess_int_t *MAXIT,double * TOL,mess_double_cpx_t *WORK,mess_int_t *IWORK,
        mess_double_cpx_t *WORK2, mess_int_t *LWORK2, mess_int_t *VERBOSE,mess_int_t *INFO);
void F77_GLOBAL(zgesgnfro,ZGESGNFRO)(mess_int_t *N,mess_double_cpx_t * A,mess_int_t *LDA, mess_int_t *MAXIT,double * TOL,mess_double_cpx_t *WORK,mess_int_t *IWORK,
        mess_double_cpx_t *WORK2, mess_int_t *LWORK2, mess_int_t *VERBOSE,mess_int_t *INFO);

//real gsign function routines
void F77_GLOBAL(dgegsgnnone,DGEGSGNNONE)(mess_int_t *N,double * A,mess_int_t *LDA, double *B,mess_int_t *LD, mess_int_t* MAXIT,double *TOL,double *WORK,mess_int_t *IWORK,mess_int_t * VERBOSE, mess_int_t *INFO);
void F77_GLOBAL(dgegsgndet,DGEGSGNDET)(mess_int_t *N,double * A,mess_int_t *LDA, double *B,mess_int_t *LD, mess_int_t* MAXIT,double *TOL,double *WORK,mess_int_t *IWORK,mess_int_t * VERBOSE, mess_int_t *INFO);

//complex gsign function routines
void F77_GLOBAL(zgegsgnnone,ZGEGSGNNONE)(mess_int_t *N,mess_double_cpx_t * A,mess_int_t *LDA, mess_double_cpx_t *B,mess_int_t *LD, mess_int_t* MAXIT,double *TOL,mess_double_cpx_t *WORK,mess_int_t *IWORK,mess_int_t * VERBOSE, mess_int_t *INFO);
void F77_GLOBAL(zgegsgndet,ZGEGSGNDET)(mess_int_t *N,mess_double_cpx_t * A,mess_int_t *LDA, mess_double_cpx_t *B,mess_int_t *LD, mess_int_t* MAXIT,double *TOL,mess_double_cpx_t *WORK,mess_int_t *IWORK, mess_int_t * VERBOSE, mess_int_t *INFO);


/**
 * @brief Compute the sign of a matrix.
 * @param[in] A  input matrix
 * @param[out] Z \f$ Z = sign(A) \f$
 * @param[in] scale scaling method
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_sign function compute the sign of a matrix using the
 * Newton sign function iteration.
 *
 *
 * @see mess_sign_scale_t
 *
 */
int  mess_eigen_sign(mess_matrix A, mess_matrix Z, mess_sign_scale_t scale){
    MSG_FNAME(__func__);
    mess_int_t lwork2, *iwork;
    double *work, *work2;
    mess_double_cpx_t *work_cpx, *work2_cpx;
    double tol;
    mess_int_t n;
    int ret = 0;
    mess_int_t info = 0;
    mess_int_t maxit;
    mess_int_t verbose = mess_error_level >2;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(Z);
    mess_check_square(A);
    MESS_MATRIX_RESET(Z);
    ret = mess_matrix_convert(A,Z,MESS_DENSE);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);

    /*-----------------------------------------------------------------------------
     *  real/complex sign function
     *-----------------------------------------------------------------------------*/
    n = A->rows;
    maxit = 40;
    tol = mess_eps()*n*n;
    lwork2 = A->rows * 64;

    if(MESS_IS_REAL(A)){
        mess_try_alloc(work, double *, sizeof(double) * 2 *n *n);
        mess_try_alloc(work2, double *, sizeof(double) * lwork2);
        mess_try_alloc(iwork, mess_int_t *, sizeof(mess_int_t) * n);

        if(scale==MESS_SIGN_SCALE_NONE){
            F77_GLOBAL(dgesgnnone,DGESGNNONE)(&n, Z->values, &(Z->ld), &maxit, &tol, work, iwork, work2, &lwork2,&verbose, &info);
            MSG_INFO("DGESGNNONE returned with info = " MESS_PRINTF_INT ", maxit = " MESS_PRINTF_INT ", tol = %lg\n", info, maxit, tol );
        }else if (scale==MESS_SIGN_SCALE_FRO){
            F77_GLOBAL(dgesgnfro,DGESGNFRO)(&n, Z->values, &(Z->ld), &maxit, &tol, work, iwork, work2, &lwork2,&verbose, &info);
            MSG_INFO("DGESGNFRO returned with info = " MESS_PRINTF_INT ", maxit = " MESS_PRINTF_INT ", tol = %lg\n", info, maxit, tol );
        }else if (scale==MESS_SIGN_SCALE_DET){
            F77_GLOBAL(dgesgndet,DGESGNDET)(&n, Z->values, &(Z->ld), &maxit, &tol, work, iwork, work2, &lwork2,&verbose, &info);
            MSG_INFO("DGESGNDET returned with info = " MESS_PRINTF_INT ", maxit = " MESS_PRINTF_INT ", tol = %lg\n", info, maxit, tol );
        }else{
            MSG_ERROR("Unknown scaling method!\n");
            return MESS_ERROR_ARGUMENTS;
        }
        mess_free(work2);
        mess_free(work);
        mess_free(iwork);

        if ( info != 0 ) {
            MSG_ERROR("dgesgn_ returned with error " MESS_PRINTF_INT "\n", info);
            return MESS_ERROR_LAPACK;
        }
    }else{
        mess_try_alloc(work_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * 2 *n *n);
        mess_try_alloc(work2_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t) * lwork2);
        mess_try_alloc(iwork, mess_int_t *, sizeof(mess_int_t) * n);

        if(scale==MESS_SIGN_SCALE_NONE){
            F77_GLOBAL(zgesgnnone,ZGESGNNONE)(&n, Z->values_cpx, &(Z->ld), &maxit, &tol, work_cpx, iwork, work2_cpx, &lwork2,&verbose, &info);
            MSG_INFO("ZGESGNNONE returned with info = " MESS_PRINTF_INT ", maxit = " MESS_PRINTF_INT ", tol = %lg\n", info, maxit, tol );
        }else if (scale==MESS_SIGN_SCALE_FRO){
            F77_GLOBAL(zgesgnfro,ZGESGNFRO)(&n, Z->values_cpx, &(Z->ld), &maxit, &tol, work_cpx, iwork, work2_cpx, &lwork2,&verbose, &info);
            MSG_INFO("ZGESGNFRO returned with info = " MESS_PRINTF_INT ", maxit = " MESS_PRINTF_INT ", tol = %lg\n", info, maxit, tol );
        }else if (scale==MESS_SIGN_SCALE_DET){
            F77_GLOBAL(zgesgndet,ZGESGNDET)(&n, Z->values_cpx, &(Z->ld), &maxit, &tol, work_cpx, iwork, work2_cpx, &lwork2,&verbose, &info);
            MSG_INFO("ZGESGNDET returned with info = " MESS_PRINTF_INT ", maxit = " MESS_PRINTF_INT ", tol = %lg\n", info, maxit, tol );
        }else{
            MSG_ERROR("Unknown scaling method!\n");
            return MESS_ERROR_ARGUMENTS;
        }
        mess_free(work2_cpx);
        mess_free(work_cpx);
        mess_free(iwork);

        if ( info != 0 ) {
            MSG_ERROR("zgesgn_ returned with error " MESS_PRINTF_INT "\n", info);
            return MESS_ERROR_LAPACK;
        }

    }
    return 0;
}       /* -----  end of function mess_eigen_sign  ----- */


/**
 * @brief Compute the sign of a matrix pencil.
 * @param[in] A  input matrix
 * @param[in] B      input matrix
 * @param[out] Z \f$ Z = sign(A,B)\f$
 * @param[in] scale scaling method
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_gsign function computes the sign of a matrix pencil \f$(A,B)\f$ using the
 * Newton sign function iteration.
 * @attention In the case of @ref MESS_SIGN_SCALE_FRO we switch automatically to @ref MESS_SIGN_SCALE_DET,
 * because the Frobenius Norm scaling \f$ c_k=\sqrt(\frac{||B^{-1}A||_F}{||A^{-1}B ||_F}) \f$ would
 * introduce an additional solve.
 *
 * @see mess_sign_scale_t
 *
 */
int  mess_eigen_gsign ( mess_matrix A, mess_matrix B, mess_matrix Z, mess_sign_scale_t scale ){
    MSG_FNAME(__func__);
    double *work, tol;
    mess_double_cpx_t *work_cpx;
    mess_int_t ret, maxit ,n, info, VERBOSE= (mess_error_level>2), *iwork;
    mess_matrix workB;
    //double ts =0 , te =0 ;
    //double flops = 0;

    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(Z);
    mess_check_square(A);
    mess_check_square(B);
    mess_check_same_size(A,B);
    mess_check_real_or_complex(A);
    mess_check_real_or_complex(B);
    mess_check_dense(A); mess_check_dense(B);

    MESS_MATRIX_RESET(Z);
    MESS_INIT_MATRICES(&workB);
    ret = mess_matrix_convert(A,Z,MESS_DENSE);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    ret = mess_matrix_convert(B,workB,MESS_DENSE);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert);

    if(MESS_IS_COMPLEX(Z) || MESS_IS_COMPLEX(workB)){
        ret = mess_matrix_tocomplex(Z);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
        ret = mess_matrix_tocomplex(workB);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
        //note the mixed datatype case MESS_REAL/MESS_COMPLEX is now impossible
    }

    /*-----------------------------------------------------------------------------
     *  real/complex gsign function
     *-----------------------------------------------------------------------------*/
    n = A->rows;
    maxit = 50;
    tol = mess_eps()*n*n;

    if(MESS_IS_REAL(Z)){
        mess_try_alloc(work, double *, sizeof(double) * 3 *n *n);
        mess_try_alloc(iwork, mess_int_t *, sizeof(mess_int_t) * n);
        //ts = mess_wtime();
        if(scale == MESS_SIGN_SCALE_NONE){
            F77_GLOBAL(dgegsgnnone,DGEGSGNNONE)(&n, Z->values, &(Z->ld), workB->values, &n, &maxit, &tol, work, iwork, &VERBOSE, &info);
            MSG_INFO("DGEGSGNNONE returned with info = " MESS_PRINTF_INT ", maxit = " MESS_PRINTF_INT ", tol = %lg\n", info, maxit, tol);
        }else if (scale == MESS_SIGN_SCALE_FRO){
            MSG_INFO("No support for frobenius scaling for generalized sign function. Automatically switch to determinant scaling. \n");
            scale = MESS_SIGN_SCALE_DET;
            F77_GLOBAL(dgegsgndet,DGEGSGNDET)(&n, Z->values, &(Z->ld), workB->values, &n, &maxit, &tol, work, iwork, &VERBOSE, &info);
            MSG_INFO("DGEGSGNDET returned with info = " MESS_PRINTF_INT ", maxit = " MESS_PRINTF_INT ", tol = %lg\n", info, maxit, tol);
        }else if ( scale == MESS_SIGN_SCALE_DET){
            F77_GLOBAL(dgegsgndet,DGEGSGNDET)(&n, Z->values, &(Z->ld), workB->values, &n, &maxit, &tol, work, iwork, &VERBOSE, &info);
            MSG_INFO("DGEGSGNDET returned with info = " MESS_PRINTF_INT ", maxit = " MESS_PRINTF_INT ", tol = %lg\n", info, maxit, tol);
        }else{
            MSG_ERROR("Unknown scaling method!\n");
            return MESS_ERROR_ARGUMENTS;
        }
        mess_free(work);
        mess_free(iwork);
    }else{
        mess_try_alloc(work_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * 3 *n *n);
        mess_try_alloc(iwork, mess_int_t *, sizeof(mess_int_t) * n);

        if(scale == MESS_SIGN_SCALE_NONE){
            F77_GLOBAL(zgegsgnnone,ZGEGSGNNONE)(&n, Z->values_cpx, &(Z->ld), workB->values_cpx, &n, &maxit, &tol, work_cpx, iwork, &VERBOSE, &info);
            MSG_INFO("ZGEGSGNNONE returned with info = " MESS_PRINTF_INT ", maxit = " MESS_PRINTF_INT ", tol = %lg\n", info, maxit, tol);
        }else if (scale == MESS_SIGN_SCALE_FRO){
            MSG_INFO("No support for frobenius scaling for generalized sign function. Automatically switch to determinant scaling. \n");
            scale = MESS_SIGN_SCALE_DET;
            F77_GLOBAL(zgegsgndet,ZGEGSGNDET)(&n, Z->values_cpx, &(Z->ld), workB->values_cpx, &n, &maxit, &tol, work_cpx, iwork, &VERBOSE, &info);
            MSG_INFO("ZGEGSGNDET returned with info = " MESS_PRINTF_INT ", maxit = " MESS_PRINTF_INT ", tol = %lg\n", info, maxit, tol);
        }else if ( scale == MESS_SIGN_SCALE_DET){
            F77_GLOBAL(zgegsgndet,ZGEGSGNDET)(&n, Z->values_cpx, &(Z->ld), workB->values_cpx, &n, &maxit, &tol, work_cpx, iwork, &VERBOSE, &info);
            MSG_INFO("ZGEGSGNDET returned with info = " MESS_PRINTF_INT ", maxit = " MESS_PRINTF_INT ", tol = %lg\n", info, maxit, tol);
        }else{
            MSG_ERROR("Unknown scaling method!\n");
            return MESS_ERROR_ARGUMENTS;
        }

        mess_free(work_cpx);
        mess_free(iwork);
    }

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&workB);

    /*-----------------------------------------------------------------------------
     *  return
     *-----------------------------------------------------------------------------*/
    if ( info != 0 ) {
        MSG_ERROR("dgesgn_ returned with error " MESS_PRINTF_INT "\n", info);
        return MESS_ERROR_LAPACK;
    }
    return 0;
}       /* -----  end of function mess_eigen_sign  ----- */

