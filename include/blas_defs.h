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
 * @file include/blas_defs.h
 * @brief Fortran name mangling for @lapack and @blas.
 * @author @koehlerm
 */

#ifndef BLAS_DEFS_H
#define BLAS_DEFS_H

#ifdef __cplusplus
extern "C" {
#endif
#include "fortran_translate.h"
#include <stdbool.h>




/**
 * @brief Datatype to represent Fortran Logical.
 *
 * The @ref MESS_LAPACK_LOGICAL represents the @c LOGICAL fortran type.
 *
 **/
//typedef mess_int_t MESS_LAPACK_LOGICAL;
typedef mess_int_t MESS_LAPACK_LOGICAL;

    /*-----------------------------------------------------------------------------
     *  BLAS
     *-----------------------------------------------------------------------------*/
#ifdef ZDOTC_MKL
    void F77_GLOBAL(zdotc,ZDOTC)(mess_double_cpx_t *ret, mess_int_t *N, mess_double_cpx_t* X, mess_int_t *incx, mess_double_cpx_t *y, mess_int_t *incy);
    void F77_GLOBAL(zdotu,ZDOTU)(mess_double_cpx_t *ret, mess_int_t *N, mess_double_cpx_t* X, mess_int_t *incx, mess_double_cpx_t *y, mess_int_t *incy);
#else
    mess_double_cpx_t F77_GLOBAL(zdotc,ZDOTC)(mess_int_t *N, mess_double_cpx_t* X, mess_int_t *incx, mess_double_cpx_t *y, mess_int_t *incy);
    mess_double_cpx_t F77_GLOBAL(zdotu,ZDOTU)(mess_int_t *N, mess_double_cpx_t* X, mess_int_t *incx, mess_double_cpx_t *y, mess_int_t *incy);
#endif


    double F77_GLOBAL(ddot,DDOT)(mess_int_t *N, double *X, mess_int_t *incx, double *y, mess_int_t * incy);
    double F77_GLOBAL(dnrm2,DNRM2)(mess_int_t *N, double * X, mess_int_t *incx);
    double F77_GLOBAL(dznrm2,DZNRM2)(mess_int_t *N, mess_double_cpx_t * X, mess_int_t *incx);
    void F77_GLOBAL(daxpy,DAXPY)(mess_int_t *N, double *alpha, double *x, mess_int_t *incx, double *y, mess_int_t * incy);
    void F77_GLOBAL(dcopy,DCOPY)(mess_int_t *N, double *x, mess_int_t *incx, double *y, mess_int_t *incy);
    void F77_GLOBAL(dgemm,DGEMM)(char *transA, char *transB, mess_int_t *m, mess_int_t *n,  mess_int_t *k, double *alpha, double *A, mess_int_t *lda, double *B, mess_int_t *ldb,   double *beta, double *C, mess_int_t *ldc);
    void F77_GLOBAL(dgemv,DGEMV)(char *trans, mess_int_t *m, mess_int_t *n, double *alpha, double *a, mess_int_t *lda, double *x, mess_int_t *incx,double *beta, double *y, mess_int_t *incy);
    void F77_GLOBAL(dger,DGER)(mess_int_t *M,mess_int_t *N,double *alpha, double* x, mess_int_t *incx, double *y,mess_int_t *incy,double *a, mess_int_t *lda  );
    void F77_GLOBAL(dscal,DSCAL)(mess_int_t *N, double * alpha, double *X, mess_int_t * incx);
    void F77_GLOBAL(dswap,DSWAP)(mess_int_t *N, double * DX, mess_int_t * incx, double *DY, mess_int_t * incy);
    void F77_GLOBAL(dtrmm,DTRMM)(char *side, char *uplo, char *transA, char * diag, mess_int_t* m, mess_int_t *n, double *alpha, double *A, mess_int_t* lda,double* B, mess_int_t* ldb);
    void F77_GLOBAL(dtrsm,DTRSM)(char *side, char *uplo, char *transA, char * diag, mess_int_t* m, mess_int_t *n, double *alpha, double *A, mess_int_t* lda,double* B, mess_int_t* ldb);
    void F77_GLOBAL(dtrevc,DTREVC)(char *SIDE, char *HOWMNY, MESS_LAPACK_LOGICAL *SELECT, mess_int_t *N, double *T, mess_int_t *LDT, double *VL, mess_int_t *LDVL, double *VR, mess_int_t *LDVR, mess_int_t *MM, mess_int_t *M, double *WORK, mess_int_t *INFO);
    void F77_GLOBAL(zaxpy,ZAXPY)(mess_int_t *N, mess_double_cpx_t *alpha, mess_double_cpx_t *X, mess_int_t *incx, mess_double_cpx_t *y, mess_int_t *incy);
    void F77_GLOBAL(zcopy,ZCOPY)(mess_int_t *N, mess_double_cpx_t *x, mess_int_t *incx, mess_double_cpx_t *y, mess_int_t *incy);
    void F77_GLOBAL(zgemm,ZGEMM)(char *transA, char *transB, mess_int_t *m, mess_int_t *n,  mess_int_t *k, mess_double_cpx_t *alpha, mess_double_cpx_t *A, mess_int_t *lda, mess_double_cpx_t *B, mess_int_t *ldb,   mess_double_cpx_t *beta, mess_double_cpx_t *C, mess_int_t *ldc);
    void F77_GLOBAL(zgemv,ZGEMV)(char *trans, mess_int_t *m, mess_int_t *n, mess_double_cpx_t *alpha, mess_double_cpx_t *A, mess_int_t *lda, mess_double_cpx_t *x, mess_int_t *incx,mess_double_cpx_t *beta, mess_double_cpx_t *y, mess_int_t *incy);
    void F77_GLOBAL(zgeru,ZGERU)(mess_int_t *M,mess_int_t *N,mess_double_cpx_t *alpha, mess_double_cpx_t* x, mess_int_t *incx, mess_double_cpx_t *y,mess_int_t *incy,mess_double_cpx_t *a, mess_int_t *lda  );
    void F77_GLOBAL(zscal,ZSCAL)(mess_int_t *N, mess_double_cpx_t * alpha, mess_double_cpx_t *X, mess_int_t * incx);
    void F77_GLOBAL(zswap,ZSWAP)(mess_int_t *N, mess_double_cpx_t * DX, mess_int_t * incx, mess_double_cpx_t *DY, mess_int_t * incy);
    void F77_GLOBAL(ztrevc,ZTREVC)(char *SIDE, char *HOWMNY, MESS_LAPACK_LOGICAL *SELECT, mess_int_t *N, mess_double_cpx_t *T, mess_int_t *LDT, mess_double_cpx_t *VL, mess_int_t *LDVL, mess_double_cpx_t *VR, mess_int_t *LDVR, mess_int_t *MM, mess_int_t* M,  mess_double_cpx_t* WORK, double *RWORK, mess_int_t *INFO);

    /*----------------------------------------------------------------------------
     *  Extend BLAS
     *  ----------------------------------------------------------------------------*/
    void F77_GLOBAL(dgeadd,DGEADD) (mess_int_t *M, mess_int_t *N, double *alpha, double *A, mess_int_t *LDA, double *beta, double *B, mess_int_t *LDB);
    void F77_GLOBAL(zgeadd,ZGEADD) (mess_int_t *M, mess_int_t *N, mess_double_cpx_t *alpha, mess_double_cpx_t  *A, mess_int_t *LDA, mess_double_cpx_t *beta, mess_double_cpx_t *B, mess_int_t *LDB);
    void F77_GLOBAL(dzgeadd,DZGEADD) (mess_int_t *M, mess_int_t *N, mess_double_cpx_t *alpha, double *A, mess_int_t *LDA, mess_double_cpx_t *beta, mess_double_cpx_t *B, mess_int_t *LDB);
    void F77_GLOBAL(dzgemmx,DZGEMMX)(char *transA, char *transB, mess_int_t *m, mess_int_t *n,  mess_int_t *k, mess_double_cpx_t *alpha, double *A, mess_int_t *lda, mess_double_cpx_t *B, mess_int_t *ldb,   mess_double_cpx_t *beta, mess_double_cpx_t *C, mess_int_t *ldc);
    void F77_GLOBAL(zdgemmx,ZDGEMMX)(char *transA, char *transB, mess_int_t *m, mess_int_t *n,  mess_int_t *k, mess_double_cpx_t *alpha, mess_double_cpx_t *A, mess_int_t *lda, double *B, mess_int_t *ldb,   mess_double_cpx_t *beta, mess_double_cpx_t *C, mess_int_t *ldc);
    mess_double_cpx_t F77_GLOBAL(dzdotu,DZDOTU)(mess_int_t *N, double *ZX, mess_int_t *INCX, mess_double_cpx_t *ZY, mess_int_t *INCY);
    mess_double_cpx_t F77_GLOBAL(zddotc,ZDDOTC)(mess_int_t *N, mess_double_cpx_t *ZX, mess_int_t *INCX, double *ZY, mess_int_t *INCY);
    mess_double_cpx_t F77_GLOBAL(zddotu,ZDDOTU)(mess_int_t *N, mess_double_cpx_t *ZX, mess_int_t *INCX, double *ZY, mess_int_t *INCY);

    /*-----------------------------------------------------------------------------
     *  LAPACK
     *-----------------------------------------------------------------------------*/
    double F77_GLOBAL(dlange,DLANGE)(char* NORM, mess_int_t* M, mess_int_t* N, double* A, mess_int_t* LDA, double* WORK);
    double F77_GLOBAL(zlange,ZLANGE)(char* NORM, mess_int_t* M, mess_int_t* N, mess_double_cpx_t* A, mess_int_t* LDA, double* WORK);
    void F77_GLOBAL(dgbtrf,DGBTRF)(mess_int_t * M, mess_int_t * N, mess_int_t * KL, mess_int_t * KU, double * AB, mess_int_t * LDAB, mess_int_t * IPIV, mess_int_t *INFO);
    void F77_GLOBAL(dgbtrs,DGBTRS)(char *TRANS, mess_int_t * N, mess_int_t * KL, mess_int_t * KU, mess_int_t *NRHS, double * AB, mess_int_t * LDAB, mess_int_t * IPIV, double * B, mess_int_t * LDB, mess_int_t *INFO);
    void F77_GLOBAL(dgees,DGEES)(char *JOBVS, char *SORT, MESS_LAPACK_LOGICAL *SELECT, mess_int_t *N, double *A, mess_int_t *SDA, mess_int_t *SDIM, double *WR, double *WI, double *VS, mess_int_t *LVDS, double *work, mess_int_t *WORK, mess_int_t * BWORK, mess_int_t *INFO);
    void F77_GLOBAL(dgeev,DGEEV)(const char *JOBVL, const char *JOBVR, mess_int_t *N, double *A, mess_int_t *LDA,double *WR, double *WI, double *VL, mess_int_t *LDVL, double *VR, mess_int_t *LDVR, double *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dgelqf,DGELQF)(mess_int_t *M, mess_int_t *N, double *A, mess_int_t* LDA, double *TAU, double *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dgels,DGELS)(char* TRANS, mess_int_t *M, mess_int_t *N, mess_int_t *NRHS, double *A, mess_int_t *LDA, double *B, mess_int_t *LDB, double *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dgeqp3,DGEQP3)(mess_int_t *M, mess_int_t* N, double *A, mess_int_t *LDA, mess_int_t *JPVT, double * TAU, double *WORK, mess_int_t *LWORK, mess_int_t *info);
    void F77_GLOBAL(dgeqrf,DGEQRF)(mess_int_t *M, mess_int_t *N, double *A, mess_int_t *LDA, double *TAU, double *WORK, mess_int_t * LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dgesvd,DGESVD)(char *JOBU, char * JOBVT, mess_int_t *M, mess_int_t *N, double *A, mess_int_t *LDA, double *S, double *U, mess_int_t *LDU,double *VT, mess_int_t* LDVT, double *WORK, mess_int_t*LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dgetrf,DGETRF)(mess_int_t *M, mess_int_t *N, double *A, mess_int_t *LDA,  mess_int_t *ipiv, mess_int_t *info);
    void F77_GLOBAL(dgetri,DGETRI)(mess_int_t *N, double *A, mess_int_t *LDA, mess_int_t *IPIV, double *WORK, mess_int_t *LWORK, mess_int_t*INFO);
    void F77_GLOBAL(dgetrs,DGETRS)(char * TRANS, mess_int_t *N, mess_int_t *NRHS, double *A, mess_int_t *LDA, mess_int_t *IPIV, double *B, mess_int_t *LDB, mess_int_t *INFO);
    void F77_GLOBAL(dgges,DGGES)(char *JOBVSL, char * JOBVSR, char * SORT, MESS_LAPACK_LOGICAL *SELECT, mess_int_t * N, double *A, mess_int_t *LDA, double *B, mess_int_t *LDB, mess_int_t *SDIM, double *ALPHAR, double* ALPHAI, double *BETA, double *VSL, mess_int_t *LDVSL, double * VSR, mess_int_t *LDVSR, double *work, mess_int_t *LWORK, mess_int_t *BWORK, mess_int_t *INFO);
    void F77_GLOBAL(dggev,DGGEV)(char *JOBVL, char * JOBVR, mess_int_t * N, double *A, mess_int_t * LDA, double *B, mess_int_t *LDB, double *ALPHAR, double *ALPHAI, double *BETA, double *VL,
            mess_int_t *LDVL, double *VR, mess_int_t *LDVR, double *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dlacpy,DLACPY)(char * uplo, mess_int_t *M, mess_int_t *N, double *A, mess_int_t *LDA, double *B, mess_int_t *LDB);
    void F77_GLOBAL(dlapmt,DLAPMT)(mess_int_t *FORWARD, mess_int_t *M, mess_int_t *N, double* x, mess_int_t *LDX,mess_int_t *K);
    void F77_GLOBAL(dlarnv,DLARNV)( mess_int_t *type, mess_int_t *seed, mess_int_t *n2, double *values );
    void F77_GLOBAL(dhseqr,DHSEQR)(char* JOB, char* COMPZ, mess_int_t *N, mess_int_t *ILO, mess_int_t *IHI, double *H, mess_int_t *LDH, double *WR, double *WI, double *Z, mess_int_t *LDZ, double *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dhsein,DHSEIN)(char* SIDE, char* EIGSRC, char* INITV, MESS_LAPACK_LOGICAL *SELECT, mess_int_t *N, double *H, mess_int_t *LDH, double* WR, double* WI, double* VL, mess_int_t *LDVL, double* VR, mess_int_t* LDVR, mess_int_t *MM, mess_int_t *M, double *WORK, mess_int_t *IFAILL, mess_int_t *IFAILR, mess_int_t *INFO);

    void F77_GLOBAL(dorglq,DORGLQ)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A, mess_int_t *LDA, double *TAU, double *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dorgqr,DORGQR)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A, mess_int_t *LDA, double *TAU, double *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dorgqr,DORGQR)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A, mess_int_t *LDA, double *TAU, double *work, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dormlq,DORMLQ)(char* SIDE, char* TRANS, mess_int_t* M, mess_int_t* N, mess_int_t* K, double* A, mess_int_t* LDA, double* TAU, double* C,  mess_int_t *LDC, double* WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dormqr,DORMQR)(char *SIDE, char* Trans, mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A, mess_int_t *LDA, double *TAU, double *C, mess_int_t *LDC, double *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dormqr,DORMQR)(char *SIDE, char* Trans, mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A, mess_int_t *LDA, double *TAU, double *C, mess_int_t *LDC, double *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(dpotrf,DPOTRF)(const char * UPLO, mess_int_t *N , double *A, mess_int_t *LDA, mess_int_t *INFO);
    void F77_GLOBAL(dpotri,DPOTRI)(const char * UPLO, mess_int_t *N, double *A, mess_int_t *LDA, mess_int_t *INFO);
    void F77_GLOBAL(dpotrs,DPOTRS)(const char * uplo, mess_int_t *n, mess_int_t * nrhs, double *a, mess_int_t *lda, double *b, mess_int_t *ldb, mess_int_t *info);
    void F77_GLOBAL(dstev,DSTEV)(const char *JOBZ, mess_int_t *N, double *D, double* E,double *Z, mess_int_t *LDZ, double *WORK, mess_int_t *INFO);
    void F77_GLOBAL(dsyev,DSYEV)(char* JOBZ, char* UPLO, mess_int_t* N, double *A, mess_int_t* LDA, double* W, double* WORK, mess_int_t* LWORK, mess_int_t* INFO);
    void F77_GLOBAL(dtgevc,DTGEVC)( char *side, char *how, MESS_LAPACK_LOGICAL* SELECT, mess_int_t *N , double *A, mess_int_t *LDA, double *B, mess_int_t *LDB, double *VL, mess_int_t *LDVL, double *VR, mess_int_t *LDVR, mess_int_t *MM, mess_int_t *M, double *WORK, mess_int_t *INFO);
    void F77_GLOBAL(dtgsyl,DTGSYL)(char *TRANS, mess_int_t *ijob, mess_int_t *M, mess_int_t *N, double *A, mess_int_t *LDA, double *B, mess_int_t *LDB, double *C, mess_int_t *LDC, double *D, mess_int_t *LDD, double *E, mess_int_t *LDE, double *F, mess_int_t *LDF, double *SCALE, double *DIFF, double *work, mess_int_t *LDWORK, mess_int_t *Iwork, mess_int_t * info);
    void F77_GLOBAL(dtrsm,DTRSM)(char* SIDE, char* UPLO, char* TRANSA, char * UNIT, mess_int_t *M, mess_int_t *N, double *alpha, double *A, mess_int_t* LDA, double *B, mess_int_t *ldb);
    void F77_GLOBAL(dtrsyl,DTRSYL)(char *transA, char *transB, mess_int_t *ISGN, mess_int_t * M, mess_int_t * N, double *A, mess_int_t *LDA, double *B, mess_int_t *LDB, double *C, mess_int_t *LDC, double *scale, mess_int_t *info);
    void F77_GLOBAL(dtrtrs,DTRTRS)(char* UPLO, char *trans, char* diag, mess_int_t *N, mess_int_t *NRHS, double *A, mess_int_t *LDA, double *B, mess_int_t *LDB, mess_int_t * INFO);
    void F77_GLOBAL(dzgbtrs,DZGBTRS)(char *TRANS, mess_int_t * N, mess_int_t * KL, mess_int_t * KU, mess_int_t *NRHS, double * AB, mess_int_t * LDAB, mess_int_t * IPIV, mess_double_cpx_t * B, mess_int_t * LDB, mess_int_t *INFO);
    void F77_GLOBAL(dzgetrs,DZGETRS)(char * TRANS, mess_int_t *N, mess_int_t *NRHS, double *A, mess_int_t *LDA, mess_int_t *IPIV, mess_double_cpx_t *B, mess_int_t *LDB, mess_int_t *INFO);
    void F77_GLOBAL(dzpotrs,DZPOTRS)(const char * uplo, mess_int_t *n, mess_int_t * nrhs, double *a, mess_int_t *lda, mess_double_cpx_t  *b, mess_int_t *ldb, mess_int_t *info);
    void F77_GLOBAL(zgbtrf,ZGBTRF)(mess_int_t * M, mess_int_t * N, mess_int_t * KL, mess_int_t * KU, mess_double_cpx_t * AB, mess_int_t * LDAB, mess_int_t * IPIV, mess_int_t *INFO);
    void F77_GLOBAL(zgbtrs,ZGBTRS)(char *TRANS, mess_int_t * N, mess_int_t * KL, mess_int_t * KU, mess_int_t *NRHS, mess_double_cpx_t * AB, mess_int_t * LDAB, mess_int_t * IPIV, mess_double_cpx_t * B, mess_int_t * LDB, mess_int_t *INFO);
    void F77_GLOBAL(zgees,ZGEES)(char * JOBVS, char *SORT, MESS_LAPACK_LOGICAL *SELECT, mess_int_t *N, mess_double_cpx_t *A, mess_int_t* LDA, mess_int_t *SDIM, mess_double_cpx_t* W, mess_double_cpx_t *VS, mess_int_t *LDVS, mess_double_cpx_t* WORK, mess_int_t *LWORK, double *RWORK, mess_int_t* BWORK, mess_int_t *INFO);
    void F77_GLOBAL(zgeev,ZGEEV)(const char *JOBVL, const char *JOBVR, mess_int_t *N, mess_double_cpx_t *A,mess_int_t *LDA, mess_double_cpx_t *W, mess_double_cpx_t *VL, mess_int_t *LDVL,mess_double_cpx_t *VR, mess_int_t *LDVR, mess_double_cpx_t *WORK, mess_int_t *LWORK,double* RWORK,mess_int_t *INFO );
    void F77_GLOBAL(zgelqf,ZGELQF)(mess_int_t *M, mess_int_t *N, mess_double_cpx_t *A, mess_int_t* LDA, mess_double_cpx_t *TAU, mess_double_cpx_t *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(zgels,ZGELS)(char* TRANS, mess_int_t *M, mess_int_t *N, mess_int_t *NRHS, mess_double_cpx_t *A, mess_int_t *LDA, mess_double_cpx_t* B, mess_int_t *LDB, mess_double_cpx_t *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(zgeqrf,ZGEQRF)(mess_int_t *M, mess_int_t *N, mess_double_cpx_t *A, mess_int_t *LDA, mess_double_cpx_t *TAU, mess_double_cpx_t *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(zgesvd,ZGESVD)(char* JOBU, char * JOBVT, mess_int_t *M, mess_int_t *N, mess_double_cpx_t *A, mess_int_t *LDA, double *S, mess_double_cpx_t *U, mess_int_t *LDU, mess_double_cpx_t *VT, mess_int_t *LDVT, mess_double_cpx_t *WORK, mess_int_t *LWORK, double *RWORK, mess_int_t *INFO);
    void F77_GLOBAL(zgetrf,ZGETRF)(mess_int_t *M, mess_int_t *N, mess_double_cpx_t *A, mess_int_t *LDA,  mess_int_t *ipiv, mess_int_t *info);
    void F77_GLOBAL(zgetri,ZGETRI)(mess_int_t *N, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *IPIV, mess_double_cpx_t *work, mess_int_t *lwork, mess_int_t *info);
    void F77_GLOBAL(zgetrs,ZGETRS)(char * TRANS, mess_int_t *N, mess_int_t *NRHS, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *IPIV, mess_double_cpx_t *B, mess_int_t *LDB, mess_int_t *INFO);
    void F77_GLOBAL(zgges,ZGGES)(char *JOBVSL, char *JOBVSR, char* SORT, MESS_LAPACK_LOGICAL *SELECT, mess_int_t *N, mess_double_cpx_t *A, mess_int_t *LDA, mess_double_cpx_t *B , mess_int_t *LDB, mess_int_t *SDIM, mess_double_cpx_t * alpha, mess_double_cpx_t * beta, mess_double_cpx_t *VSL, mess_int_t *LDVSL, mess_double_cpx_t *VSR, mess_int_t *LDVSR, mess_double_cpx_t *work, mess_int_t *lwork, double *rwork, mess_int_t *BWORK, mess_int_t *INFO);
    void F77_GLOBAL(zggev,ZGGEV)(const char *JOBL, const char *JOBR, mess_int_t *N, mess_double_cpx_t *A, mess_int_t *LDA, mess_double_cpx_t *B, mess_int_t *LDB,mess_double_cpx_t *alpha, mess_double_cpx_t *beta, mess_double_cpx_t *VL, mess_int_t *LDVL, mess_double_cpx_t *VR, mess_int_t *LDVR, mess_double_cpx_t *work, mess_int_t *lwork, double *rwork, mess_int_t *info);
    void F77_GLOBAL(zhseqr,ZHSEQR)(char* JOB, char* COMPZ, mess_int_t *N, mess_int_t *ILO, mess_int_t *IHI, mess_double_cpx_t *H, mess_int_t *LDH, mess_double_cpx_t *W, mess_double_cpx_t *Z, mess_int_t *LDZ, mess_double_cpx_t *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(zhsein,zHSEIN)(char* SIDE, char* EIGSRC, char* INITV, MESS_LAPACK_LOGICAL *SELECT, mess_int_t *N, mess_double_cpx_t *H, mess_int_t *LDH, mess_double_cpx_t* W, mess_double_cpx_t* VL, mess_int_t *LDVL, mess_double_cpx_t* VR, mess_int_t* LDVR, mess_int_t *MM, mess_int_t *M, mess_double_cpx_t *WORK, mess_int_t *IFAILL, mess_int_t *IFAILR, mess_int_t *INFO);

    void F77_GLOBAL(zlacpy,ZLACPY)(char * uplo, mess_int_t *M, mess_int_t *N, mess_double_cpx_t  *A, mess_int_t *LDA, mess_double_cpx_t *B, mess_int_t *LDB);
    void F77_GLOBAL(zlapmt,ZLAPMT)(mess_int_t *FORWARD, mess_int_t *M, mess_int_t *N, mess_double_cpx_t* x, mess_int_t *LDX,mess_int_t *K);
    void F77_GLOBAL(zpotrf,ZPOTRF)(const char * UPLO, mess_int_t *N , mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *INFO);
    void F77_GLOBAL(zpotri,ZPOTRI)(const char * UPLO, mess_int_t *N, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *INFO);
    void F77_GLOBAL(zpotrs,ZPOTRS)(const char * UPLO, mess_int_t *N, mess_int_t * NRHS, mess_double_cpx_t *A, mess_int_t *LDA, mess_double_cpx_t *B, mess_int_t *LDB, mess_int_t *INFO);
    void F77_GLOBAL(ztrsm,ZTRSM)(char* SIDE, char* UPLO, char* TRANSA, char * UNIT, mess_int_t *M, mess_int_t *N, mess_double_cpx_t *alpha, mess_double_cpx_t *A, mess_int_t* LDA, mess_double_cpx_t *B, mess_int_t *ldb);
    void F77_GLOBAL(ztrsyl,ZTRSYL)(char *transA, char *transB, mess_int_t *ISGN, mess_int_t * M, mess_int_t * N, mess_double_cpx_t *A, mess_int_t *LDA, mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC, double *scale, mess_int_t *info);
    void F77_GLOBAL(ztgsyl,ZTGSYL)(char *TRANS, mess_int_t *ijob, mess_int_t *M, mess_int_t *N, mess_double_cpx_t *A, mess_int_t *LDA, mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC, mess_double_cpx_t *D, mess_int_t *LDD, mess_double_cpx_t *E, mess_int_t *LDE, mess_double_cpx_t *F, mess_int_t *LDF, double *SCALE, double *DIFF, mess_double_cpx_t *work, mess_int_t *LDWORK, mess_int_t *Iwork, mess_int_t * info);
    void F77_GLOBAL(ztrtrs,ZTRTRS)(char* UPLO, char *TRANS, char* diag, mess_int_t *N, mess_int_t *NRHS, mess_double_cpx_t *A, mess_int_t *LDA, mess_double_cpx_t *B, mess_int_t *LDB, mess_int_t *INFO);
    void F77_GLOBAL(zunglq,ZUNGLQ)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_double_cpx_t *TAU, mess_double_cpx_t *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(zungqr,ZUNGQR)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_double_cpx_t *TAU, mess_double_cpx_t *WORK, mess_int_t *LWORK, mess_int_t *INFO);
    void F77_GLOBAL(zunmlq,ZUNMLQ)(char *SIDE, char *TRANS, mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA,  mess_double_cpx_t *TAU, mess_double_cpx_t *C, mess_int_t *LDC, mess_double_cpx_t *WORK, mess_int_t *LWORK, mess_int_t *INFO);
   void F77_GLOBAL(zunmqr,ZUNMQR)(char *SIDE, char* TRANS, mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_double_cpx_t *TAU,  mess_double_cpx_t *C, mess_int_t *LDC, mess_double_cpx_t *work, mess_int_t *lwork, mess_int_t *INFO);


   /*-----------------------------------------------------------------------------
     *  glyap3 stuff
     *-----------------------------------------------------------------------------*/
    void F77_GLOBAL(dgglyp,DGGLYP)(char *DICO, char *JOB, char* FACT, char *TRANS, char *UPLO, mess_int_t *ISOLVE, mess_int_t *NB,
            mess_int_t *N, double *A, mess_int_t *LDA,
            double *E, mess_int_t *LDE, double *Q, mess_int_t *LDQ, double * Z, mess_int_t *LDZ, double *X, mess_int_t *LDX,
            double *SCALE, double *sep, double *ferr, double *ALPHAR, double *ALPHAI, double *BETA,
            mess_int_t *IWORK, double *DWORK, mess_int_t *LDWORK, mess_int_t *INFO);

    void F77_GLOBAL(dgelyp,DGELYP)(char *DICO, char* JOB, char *FACT, char*TRANSA, mess_int_t *N, double *A, mess_int_t *LDA, double* U, mess_int_t *LDU, double *C, mess_int_t *LDC,
            double *SCALE, double *SEP, double *FERR, double *WR, double *WI, mess_int_t *IWORK, double *DWORK, mess_int_t *LDWORK, mess_int_t* INFO);

#ifdef __cplusplus
};
#endif

#endif /* end of include guard: BLAS_DEFS_H */
