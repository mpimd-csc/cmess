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
 * @file lib/direct/multisolver/multilu_umfpack.c
 * @brief @umfpack multisolver interface.
 * @author @koehlerm
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <umfpack.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/stat.h>

#define IS_REAL(X)  (fabs(cimag((X))) < eps)
#define IS_COMPLEX(X) (fabs(cimag((X))) >=  eps)

#ifdef MESS64
#define umfpack_dl_solve            umfpack_dl_solve
#define umfpack_zl_solve            umfpack_zl_solve
#define umfpack_zl_load_numeric     umfpack_zl_load_numeric
#define umfpack_dl_load_numeric     umfpack_dl_load_numeric
#define umfpack_zl_defaults         umfpack_zl_defaults
#define umfpack_dl_defaults         umfpack_dl_defaults
#define umfpack_dl_save_numeric     umfpack_dl_save_numeric
#define umfpack_zl_save_numeric     umfpack_zl_save_numeric
#define umfpack_zl_free_numeric     umfpack_zl_free_numeric
#define umfpack_dl_free_numeric     umfpack_dl_free_numeric
#define umfpack_dl_symbolic         umfpack_dl_symbolic
#define umfpack_dl_numeric          umfpack_dl_numeric
#define umfpack_zl_symbolic         umfpack_zl_symbolic
#define umfpack_zl_numeric          umfpack_zl_numeric
#define umfpack_dl_report_info      umfpack_dl_report_info
#define umfpack_dl_report_numeric   umfpack_dl_report_numeric
#define umfpack_dl_free_symbolic    umfpack_dl_free_symbolic
#define umfpack_zl_report_info      umfpack_zl_report_info
#define umfpack_zl_report_numeric   umfpack_zl_report_numeric
#define umfpack_zl_free_symbolic    umfpack_zl_free_symbolic
#else
#define umfpack_dl_solve            umfpack_di_solve
#define umfpack_zl_solve            umfpack_zi_solve
#define umfpack_zl_load_numeric     umfpack_zi_load_numeric
#define umfpack_dl_load_numeric     umfpack_di_load_numeric
#define umfpack_zl_defaults         umfpack_zi_defaults
#define umfpack_dl_defaults         umfpack_di_defaults
#define umfpack_dl_save_numeric     umfpack_di_save_numeric
#define umfpack_zl_save_numeric     umfpack_zi_save_numeric
#define umfpack_zl_free_numeric     umfpack_zi_free_numeric
#define umfpack_dl_free_numeric     umfpack_di_free_numeric
#define umfpack_dl_symbolic         umfpack_di_symbolic
#define umfpack_dl_numeric          umfpack_di_numeric
#define umfpack_zl_symbolic         umfpack_zi_symbolic
#define umfpack_zl_numeric          umfpack_zi_numeric
#define umfpack_dl_report_info      umfpack_di_report_info
#define umfpack_dl_report_numeric   umfpack_di_report_numeric
#define umfpack_dl_free_symbolic    umfpack_di_free_symbolic
#define umfpack_zl_report_info      umfpack_zi_report_info
#define umfpack_zl_report_numeric   umfpack_zi_report_numeric
#define umfpack_zl_free_symbolic    umfpack_zi_free_symbolic
#endif


struct multi_umfpack {
    void **numeric;
    unsigned short cpx;
    mess_int_t rows, cols;
    mess_int_t nlu;
    double **Control;
    double **Info;
    unsigned short *datatypes;
};


/**
 * @brief Provide the solve function for \f$ A(n)x=b \f$.
 * @param[in] data   input the solver data
 * @param[in] n  input index of the matrix
 * @param[in] b  input right hand side
 * @param[out] x solution
 * @return zero on success or a non-zero error value
 *
 * The @ref multiumf_solve function solve the system
 * \f[A(n)x=b\f]
 * for a single right-hand-side.
 *
 */
static int multiumf_solve(void *data, mess_int_t n, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct multi_umfpack *d = (struct multi_umfpack *) data;
    int ret = 0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( n < 0 || n>=d->nlu){
        MSG_ERROR("n is out of range ( = " MESS_PRINTF_INT " )\n", n);
        return MESS_ERROR_ARGUMENTS;
    }
    if ( x->dim != b->dim){
        MSG_WARN("resize x from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", x->dim, b->dim);
        mess_vector_resize(x, b->dim);
    }
    if (b->dim != d->rows ) {
        MSG_ERROR("The dimension of the right hand side don't fit.\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  real solver and real right hand side
     *-----------------------------------------------------------------------------*/
    if ( d->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)) {
        ret = mess_vector_toreal_nowarn(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        ret = umfpack_dl_solve (UMFPACK_A, NULL, NULL, NULL, x->values, b->values, d->numeric[n], d->Control[n], d->Info[n]) ;
        if ( ret != UMFPACK_OK) {
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_dl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }
    }

    /*-----------------------------------------------------------------------------
     *  real solver complex rhs
     *-----------------------------------------------------------------------------*/
    else if ( d->datatypes[n] == MESS_REAL && MESS_IS_COMPLEX(b) ){
        mess_vector tmp;
        mess_vector xr,xc;
        ret = mess_vector_tocomplex(x);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_resize(x, b->dim);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        MESS_INIT_VECTORS(&tmp,&xr,&xc);
        ret = mess_vector_alloc(tmp, b->dim, MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_alloc(xr, b->dim, MESS_REAL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_alloc(xc, b->dim, MESS_REAL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);

        // solve realpart
        ret = mess_vector_realpart(b, tmp);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_realpart);
        ret =  umfpack_dl_solve (UMFPACK_A, NULL, NULL, NULL, xr->values, tmp->values, d->numeric[n], d->Control[n], d->Info[n]) ;
        if ( ret != UMFPACK_OK) {
            mess_vector_clear(&xr);
            mess_vector_clear(&xc);
            mess_vector_clear(&tmp);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_dl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }
        // solve imag part
        ret = mess_vector_imagpart(b, tmp);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_imagpart);
        ret =  umfpack_dl_solve (UMFPACK_A, NULL, NULL, NULL, xc->values, tmp->values, d->numeric[n], d->Control[n], d->Info[n]) ;
        if ( ret != UMFPACK_OK) {
            mess_vector_clear(&xr);
            mess_vector_clear(&xc);
            mess_vector_clear(&tmp);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_dl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }
        // combine
        ret = mess_vector_complex_from_parts(xr, xc, x);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_complex_from_parts);
        mess_vector_clear(&xr);
        mess_vector_clear(&xc);
        mess_vector_clear(&tmp);
    }

    /*-----------------------------------------------------------------------------
     *  complex solver and real rhs
     *-----------------------------------------------------------------------------*/
    else if ( d->datatypes[n] == MESS_COMPLEX && MESS_IS_REAL(b)) {
        mess_vector bimag;
        ret = mess_vector_tocomplex(x);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_resize(x, b->dim);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        ret = mess_vector_init(&bimag);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(bimag, b->dim, MESS_COMPLEX);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_copy_tocomplex(b, bimag);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy_tocomplex);
        ret = umfpack_zl_solve (UMFPACK_A, NULL, NULL, NULL, NULL, (double *)x->values_cpx, NULL, (double*) bimag->values_cpx,NULL, d->numeric[n], d->Control[n], d->Info[n]) ;
        if ( ret != UMFPACK_OK) {
            mess_vector_clear(&bimag);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_zl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }

        mess_vector_clear(&bimag);
    }

    /*-----------------------------------------------------------------------------
     *  complex solver and complex rhs
     *-----------------------------------------------------------------------------*/
    else {
        ret = mess_vector_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
        ret = mess_vector_resize(x,b->dim); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        ret = umfpack_zl_solve (UMFPACK_A, NULL, NULL, NULL, NULL, (double *)x->values_cpx, NULL, (double*) b->values_cpx,NULL, d->numeric[n], d->Control[n], d->Info[n]) ;
        if ( ret != UMFPACK_OK) {
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_zl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }

    }
    return 0;
}

/**
 * @brief Provide the solve function for \f$ A(n)X=B \f$.
 * @param[in] data   input the solver data
 * @param[in] n  input index of the matrix
 * @param[in] b  input right hand sides
 * @param[out] x solution
 * @return zero on success or a non-zero error value
 *
 * The @ref multiumf_solvem function solve the system
 * \f[A(n)X=B \f]
 * for a multiple right-hand-sides.
 *
 */
static int multiumf_solvem(void *data, mess_int_t n, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    mess_matrix work;
    int conv = -1;
    mess_int_t i;
    int ret =0;
    struct multi_umfpack *d = (struct multi_umfpack *) data;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( n < 0 || n>=d->nlu){
        MSG_ERROR("n is out of range ( = " MESS_PRINTF_INT " )\n", n);
        return MESS_ERROR_ARGUMENTS;
    }
    if ( d->rows != b->rows) {
        MSG_ERROR("b don't have the right number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }

    MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
    MESS_MATRIX_RESET(x);

    /*-----------------------------------------------------------------------------
     *  real solver and real right hand side
     *-----------------------------------------------------------------------------*/
    if ( d->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)) {
        ret =  mess_matrix_alloc(x, d->rows, work->cols,d->rows*work->cols, MESS_DENSE, MESS_REAL ); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret =0;
#ifdef _OPENMP
#pragma omp parallel for private (i)  default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret += umfpack_dl_solve (UMFPACK_A, NULL, NULL, NULL, x->values+(i*x->ld), work->values+(i*work->ld), d->numeric[n], d->Control[n], d->Info[n]) ;
        }
        if ( ret != 0 ) {
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_dl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }
    }

    /*-----------------------------------------------------------------------------
     *  real solver complex rhs
     *-----------------------------------------------------------------------------*/
    else if ( d->datatypes[n] == MESS_REAL && MESS_IS_COMPLEX(b) ){
        mess_matrix tmp;
        mess_matrix xr,xc;
        ret = mess_matrix_init(&tmp);                                                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&xr);                                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&xc);                                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_alloc(xr, d->rows, work->cols,d->rows*work->cols, MESS_DENSE, MESS_REAL );    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(xc, d->rows, work->cols,d->rows*work->cols, MESS_DENSE, MESS_REAL );    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

        // solve realpart
        ret = mess_matrix_realpart(work, tmp);                                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realpart);

#ifdef _OPENMP
#pragma omp parallel for private (i)  default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret += umfpack_dl_solve (UMFPACK_A, NULL, NULL, NULL, xr->values+(i*xr->ld), tmp->values+(i*work->ld), d->numeric[n], d->Control[n], d->Info[n]) ;
        }
        if ( ret != 0 ) {
            mess_matrix_clear(&xr);
            mess_matrix_clear(&xc);
            mess_matrix_clear(&tmp);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_dl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }

        // solve imag part
        ret = mess_matrix_imagpart(work, tmp);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_imagpart);
#ifdef _OPENMP
#pragma omp parallel for private (i)  default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret += umfpack_dl_solve (UMFPACK_A, NULL, NULL, NULL, xc->values+(i*xc->ld), tmp->values+(i*work->ld), d->numeric[n], d->Control[n], d->Info[n]) ;
        }
        if ( ret != 0 ) {
            mess_matrix_clear(&xr);
            mess_matrix_clear(&xc);
            mess_matrix_clear(&tmp);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_dl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }
        ret = mess_matrix_complex_from_parts(xr, xc, x);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_complex_from_parts);
        mess_matrix_clear(&xr);
        mess_matrix_clear(&xc);
        mess_matrix_clear(&tmp);
    }

    /*-----------------------------------------------------------------------------
     *  complex solver and real rhs
     *-----------------------------------------------------------------------------*/
    else if ( d->datatypes[n] == MESS_COMPLEX && MESS_IS_REAL(b)) {

        mess_matrix bimag;
        ret =  mess_matrix_alloc(x, d->rows, work->cols,d->rows*work->cols, MESS_DENSE, MESS_COMPLEX );     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret =  mess_matrix_init(&bimag);                                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret =  mess_matrix_copy(work, bimag);                                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        ret =  mess_matrix_tocomplex(bimag);                                                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        ret = 0;
#ifdef _OPENMP
#pragma omp parallel for private (i)  default(shared)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret =  umfpack_zl_solve (UMFPACK_A, NULL, NULL, NULL, NULL, (double *) (x->values_cpx+(i*x->ld)),NULL,  (double*) (bimag->values_cpx+(i*bimag->ld)),NULL, d->numeric[n], d->Control[n], d->Info[n]) ;
        }
        if ( ret != UMFPACK_OK) {
            mess_matrix_clear(&bimag);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_zl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }
        mess_matrix_clear(&bimag);
    }

    /*-----------------------------------------------------------------------------
     *  complex solver and complex rhs
     *-----------------------------------------------------------------------------*/
    else {
        ret = mess_matrix_alloc(x, d->rows, work->cols,d->rows*work->cols, MESS_DENSE, MESS_COMPLEX );      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
        ret = 0 ;
#ifdef _OPENMP
#pragma omp parallel for private (i)  default(shared)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret = umfpack_zl_solve (UMFPACK_A, NULL, NULL, NULL, NULL, (double *) (x->values_cpx+(i*x->ld)),NULL,  (double*) (work->values_cpx+(i*work->ld)),NULL, d->numeric[n], d->Control[n], d->Info[n]) ;
        }
        if ( ret != UMFPACK_OK) {
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_zl_report_info(d->Control[n], d->Info[n]);
            abort ();
            return MESS_ERROR_GENERAL;
        }
    }

    if ( conv == 0) {
        mess_matrix_clear(&work);
    }
    return 0;
}

/**
 * @brief Provide the solve function for \f$ A(n)^Tx=b \f$.
 * @param[in] data   input the solver data
 * @param[in] n  input index of the matrix
 * @param[in] b  input right hand side
 * @param[out] x solution
 * @return zero on success or a non-zero error value
 *
 * The @ref multiumf_solvet function solve the system
 * \f[A(n)^Tx=b\f]
 * for a single right-hand-side.
 *
 */
static int multiumf_solvet(void *data, mess_int_t n, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct multi_umfpack *d = (struct multi_umfpack *) data;
    int ret = 0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( n < 0 || n>=d->nlu){
        MSG_ERROR("n is out of range ( = " MESS_PRINTF_INT " )\n", n);
        return MESS_ERROR_ARGUMENTS;
    }
    if ( x->dim != b->dim){
        MSG_WARN("resize x from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", x->dim, b->dim);
        mess_vector_resize(x, b->dim);
    }
    if (b->dim != d->rows ) {
        MSG_ERROR("The dimension of the right hand side don't fit.\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  real solver and real right hand side
     *-----------------------------------------------------------------------------*/
    if ( d->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)) {
        ret = mess_vector_toreal_nowarn(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        ret = mess_vector_resize (x, b->dim);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        ret = umfpack_dl_solve (UMFPACK_At, NULL, NULL, NULL, x->values, b->values, d->numeric[n], d->Control[n], d->Info[n]) ;
        if ( ret != UMFPACK_OK) {
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_dl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }

    }

    /*-----------------------------------------------------------------------------
     *  real solver complex rhs
     *-----------------------------------------------------------------------------*/
    else if ( d->datatypes[n] == MESS_REAL && MESS_IS_COMPLEX(b) ){
        mess_vector tmp;
        mess_vector xr,xc;
        ret = mess_vector_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_resize(x, b->dim); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);

        MESS_INIT_VECTORS(&tmp,&xr,&xc);
        ret = mess_vector_alloc(tmp, b->dim, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(xr, b->dim, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(xc, b->dim, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);

        // solve realpart
        ret = mess_vector_realpart(b, tmp); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_realpart);
        ret = umfpack_dl_solve (UMFPACK_At, NULL, NULL, NULL, xr->values, tmp->values, d->numeric[n], d->Control[n], d->Info[n]) ;
        if ( ret != UMFPACK_OK) {
            mess_vector_clear(&xr);
            mess_vector_clear(&xc);
            mess_vector_clear(&tmp);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_dl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }


        // solve imag part
        ret = mess_vector_imagpart(b, tmp);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_imagpart);
        ret =  umfpack_dl_solve (UMFPACK_At, NULL, NULL, NULL, xc->values, tmp->values, d->numeric[n], d->Control[n], d->Info[n]) ;
        if ( ret != UMFPACK_OK) {
            mess_vector_clear(&xr);
            mess_vector_clear(&xc);
            mess_vector_clear(&tmp);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_dl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }


        // combine
        ret = mess_vector_complex_from_parts(xr, xc, x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_complex_from_parts);
        mess_vector_clear(&xr);
        mess_vector_clear(&xc);
        mess_vector_clear(&tmp);
    }

    /*-----------------------------------------------------------------------------
     *  complex solver and real rhs
     *-----------------------------------------------------------------------------*/
    else if ( d->datatypes[n] == MESS_COMPLEX && MESS_IS_REAL(b)) {
        mess_vector bimag;
        ret = mess_vector_tocomplex(x);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_resize(x, b->dim);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        mess_vector_init(&bimag);
        ret = mess_vector_alloc(bimag, b->dim, MESS_COMPLEX);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_copy_tocomplex(b, bimag);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy_tocomplex);
        ret = umfpack_zl_solve (UMFPACK_Aat, NULL, NULL, NULL, NULL, (double *)x->values_cpx, NULL, (double*) bimag->values_cpx,NULL, d->numeric[n], d->Control[n], d->Info[n]) ;
        if ( ret != UMFPACK_OK) {
            mess_vector_clear(&bimag);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_zl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }

        mess_vector_clear(&bimag);
    }

    /*-----------------------------------------------------------------------------
     *  complex solver and complex rhs
     *-----------------------------------------------------------------------------*/
    else {
        ret = mess_vector_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
        ret = mess_vector_resize(x,b->dim); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        ret =  umfpack_zl_solve (UMFPACK_Aat, NULL, NULL, NULL, NULL, (double *)x->values_cpx, NULL, (double*) b->values_cpx,NULL, d->numeric[n], d->Control[n], d->Info[n]) ;
        if ( ret != UMFPACK_OK) {
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_zl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }

    }
    return 0;
}

/**
 * @brief Provide the solve function for \f$ A(n)^Hx=b \f$.
 * @param[in] data   input the solver data
 * @param[in] n  input index of the matrix
 * @param[in] b  input right hand side
 * @param[out] x solution
 * @return zero on success or a non-zero error value
 *
 * The @ref multiumf_solvet function solve the system
 * \f[ A(n)^Hx=b \f]
 * for a single right-hand-side.
 *
 */
static int multiumf_solveh(void *data, mess_int_t n, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct multi_umfpack *d = (struct multi_umfpack *) data;
    int ret = 0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( n < 0 || n>=d->nlu){
        MSG_ERROR("n is out of range ( = " MESS_PRINTF_INT " )\n", n);
        return MESS_ERROR_ARGUMENTS;
    }
    if ( x->dim != b->dim){
        MSG_WARN("resize x from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", x->dim, b->dim);
        mess_vector_resize(x, b->dim);
    }
    if (b->dim != d->rows ) {
        MSG_ERROR("The dimension of the right hand side don't fit.\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  real solver and real right hand side
     *-----------------------------------------------------------------------------*/
    if ( d->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)) {
        return multiumf_solvet(data, n, b, x);
    }

    /*-----------------------------------------------------------------------------
     *  real solver complex rhs
     *-----------------------------------------------------------------------------*/
    else if ( d->datatypes[n] == MESS_REAL && MESS_IS_COMPLEX(b) ){
        return multiumf_solvet(data, n, b, x);
    }

    /*-----------------------------------------------------------------------------
     *  complex solver and real rhs
     *-----------------------------------------------------------------------------*/
    else if ( d->datatypes[n] == MESS_COMPLEX && MESS_IS_REAL(b)) {
        mess_vector bimag;
        ret = mess_vector_tocomplex(x);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_resize(x, b->dim);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        ret = mess_vector_init(&bimag);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(bimag, b->dim, MESS_COMPLEX);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_copy_tocomplex(b, bimag);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy_tocomplex);
        ret = umfpack_zl_solve (UMFPACK_At, NULL, NULL, NULL, NULL, (double *)x->values_cpx, NULL, (double*) bimag->values_cpx,NULL, d->numeric[n], d->Control[n], d->Info[n]) ;
        if ( ret != UMFPACK_OK) {
            mess_vector_clear(&bimag);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_zl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }

        mess_vector_clear(&bimag);
    }

    /*-----------------------------------------------------------------------------
     *  complex solver and complex rhs
     *-----------------------------------------------------------------------------*/
    else {
        ret = mess_vector_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
        ret = mess_vector_resize(x,b->dim); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        ret =  umfpack_zl_solve (UMFPACK_At, NULL, NULL, NULL, NULL, (double *)x->values_cpx, NULL, (double*) b->values_cpx,NULL, d->numeric[n], d->Control[n], d->Info[n]) ;
        if ( ret != UMFPACK_OK) {
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_zl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }

    }
    return 0;
}

/**
 * @brief Provide the solve function for \f$ A(n)^TX=B \f$.
 * @param[in] data   input the solver data
 * @param[in] n  input index of the matrix
 * @param[in] b  input right hand sides
 * @param[out] x solution
 * @return zero on success or a non-zero error value
 *
 * The @ref multiumf_solvemt function solve the system
 * \f[ A(n)^T X=B \f]
 * for a multiple right-hand-sides.
 *
 */
static int multiumf_solvemt(void *data, mess_int_t n, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    mess_matrix work;
    int conv = -1;
    mess_int_t i;
    int ret =0;
    struct multi_umfpack *d = (struct multi_umfpack *) data;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( n < 0 || n>=d->nlu){
        MSG_ERROR("n is out of range ( = " MESS_PRINTF_INT " )\n", n);
        return MESS_ERROR_ARGUMENTS;
    }
    if ( d->rows != b->rows) {
        MSG_ERROR("b don't have the right number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }

    MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
    MESS_MATRIX_RESET(x);

    /*-----------------------------------------------------------------------------
     *  real solver and real right hand side
     *-----------------------------------------------------------------------------*/
    if ( d->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)) {
        ret =  mess_matrix_alloc(x, d->rows, work->cols,d->rows*work->cols, MESS_DENSE, MESS_REAL ); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret =0;
#ifdef _OPENMP
#pragma omp parallel for private (i)  default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret += umfpack_dl_solve (UMFPACK_Aat, NULL, NULL, NULL, x->values+(i*x->ld), work->values+(i*work->ld), d->numeric[n], d->Control[n], d->Info[n]) ;
        }
        if ( ret != 0 ) {
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_dl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }
    }

    /*-----------------------------------------------------------------------------
     *  real solver complex rhs
     *-----------------------------------------------------------------------------*/
    else if ( d->datatypes[n] == MESS_REAL && MESS_IS_COMPLEX(b) ){
        mess_matrix tmp;
        mess_matrix xr,xc;
        ret = mess_matrix_init(&tmp);                                                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&xr);                                                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&xc);                                                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret =  mess_matrix_alloc(xr, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, MESS_REAL ); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret =  mess_matrix_alloc(xc, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, MESS_REAL ); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

        // solve realpart
        ret = mess_matrix_realpart(work, tmp);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realpart);

#ifdef _OPENMP
#pragma omp parallel for private (i)  default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret += umfpack_dl_solve (UMFPACK_Aat, NULL, NULL, NULL, xr->values+(i*xr->ld), tmp->values+(i*tmp->ld), d->numeric[n], d->Control[n], d->Info[n]) ;
        }
        if ( ret != 0 ) {
            mess_matrix_clear(&xr);
            mess_matrix_clear(&xc);
            mess_matrix_clear(&tmp);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_dl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }

        // solve imag part
        ret = mess_matrix_imagpart(work, tmp);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_imagpart);
#ifdef _OPENMP
#pragma omp parallel for private (i)  default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret += umfpack_dl_solve (UMFPACK_A, NULL, NULL, NULL, xc->values+(i*xc->ld), tmp->values+(i*tmp->ld), d->numeric[n], d->Control[n], d->Info[n]) ;
        }
        if ( ret != 0 ) {
            mess_matrix_clear(&xr);
            mess_matrix_clear(&xc);
            mess_matrix_clear(&tmp);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_dl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }

        ret = mess_matrix_complex_from_parts(xr, xc, x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_complex_from_parts);

        mess_matrix_clear(&xr);
        mess_matrix_clear(&xc);
        mess_matrix_clear(&tmp);
    }

    /*-----------------------------------------------------------------------------
     *  complex solver and real rhs
     *-----------------------------------------------------------------------------*/
    else if ( d->datatypes[n] == MESS_COMPLEX && MESS_IS_REAL(b)) {
        mess_matrix bimag;
        ret =  mess_matrix_alloc(x, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, MESS_COMPLEX ); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret =  mess_matrix_init(&bimag); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret =  mess_matrix_copy(work, bimag); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        ret =  mess_matrix_tocomplex(bimag); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        ret = 0;
#ifdef _OPENMP
#pragma omp parallel for private (i)  default(shared)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret =  umfpack_zl_solve (UMFPACK_Aat, NULL, NULL, NULL, NULL, (double *) (x->values_cpx+(i*x->ld)),NULL,  (double*) (bimag->values_cpx+(i*bimag->ld)),NULL, d->numeric[n], d->Control[n], d->Info[n]) ;
        }
        if ( ret != UMFPACK_OK) {
            mess_matrix_clear(&bimag);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_zl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }
        mess_matrix_clear(&bimag);
        // mess_matrix_conj(x);
    }

    /*-----------------------------------------------------------------------------
     *  complex solver and complex rhs
     *-----------------------------------------------------------------------------*/
    else {
        ret = mess_matrix_alloc(x, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, MESS_COMPLEX );        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
        ret = 0 ;
#ifdef _OPENMP
#pragma omp parallel for private (i)  default(shared)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret = umfpack_zl_solve (UMFPACK_Aat, NULL, NULL, NULL, NULL, (double *) (x->values_cpx+(i*x->ld)),NULL,  (double*) (work->values_cpx+(i*work->ld)),NULL, d->numeric[n], d->Control[n], d->Info[n]) ;
        }
        if ( ret != UMFPACK_OK) {
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_zl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }
    }

    if ( conv == 0) {
        mess_matrix_clear(&work);
    }
    return 0;
}

/**
 * @brief Provide the solve function for \f$ A(n)^H X=B \f$.
 * @param[in] data   input the solver data
 * @param[in] n  input index of the matrix
 * @param[in] b  input right hand sides
 * @param[out] x solution
 * @return zero on success or a non-zero error value
 *
 * The @ref multiumf_solvemh function solve the system
 * \f[ A(n)^H X=B \f]
 * for a multiple right-hand-sides.
 *
 */
static int multiumf_solvemh(void *data, mess_int_t n, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    mess_matrix work;
    int conv = -1;
    mess_int_t i;
    int ret =0;
    struct multi_umfpack *d = (struct multi_umfpack *) data;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( n < 0 || n>=d->nlu){
        MSG_ERROR("n is out of range ( = " MESS_PRINTF_INT " )\n", n);
        return MESS_ERROR_ARGUMENTS;
    }
    if ( d->rows != b->rows) {
        MSG_ERROR("b don't have the right number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  real solver and real right hand side
     *-----------------------------------------------------------------------------*/
    if ( d->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)) {
        return multiumf_solvemt(data, n, b,x);
    }

    /*-----------------------------------------------------------------------------
     *  real solver complex rhs
     *-----------------------------------------------------------------------------*/
    else if ( d->datatypes[n] == MESS_REAL && MESS_IS_COMPLEX(b) ){
        return multiumf_solvemt(data, n, b,x);
    }

    /*-----------------------------------------------------------------------------
     *  complex solver and real rhs
     *-----------------------------------------------------------------------------*/
    else if ( d->datatypes[n] == MESS_COMPLEX && MESS_IS_REAL(b)) {
        mess_matrix bimag;
        MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
        MESS_MATRIX_RESET(x);

        ret =  mess_matrix_alloc(x, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, MESS_COMPLEX ); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret =  mess_matrix_init(&bimag); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret =  mess_matrix_copy(work, bimag); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        ret =  mess_matrix_tocomplex(bimag); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        ret = 0;
#ifdef _OPENMP
#pragma omp parallel for private (i)  default(shared)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret =  umfpack_zl_solve (UMFPACK_At, NULL, NULL, NULL, NULL, (double *) (x->values_cpx+(i*x->ld)),NULL,  (double*) (bimag->values_cpx+(i*bimag->ld)),NULL, d->numeric[n], d->Control[n], d->Info[n]) ;
        }
        if ( ret != UMFPACK_OK) {
            mess_matrix_clear(&bimag);
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_zl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }
        mess_matrix_clear(&bimag);
    }

    /*-----------------------------------------------------------------------------
     *  complex solver and complex rhs
     *-----------------------------------------------------------------------------*/
    else {
        MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
        MESS_MATRIX_RESET(x);
        ret = mess_matrix_alloc(x, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, MESS_COMPLEX );            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
        ret = 0 ;
#ifdef _OPENMP
#pragma omp parallel for private (i)  default(shared)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret = umfpack_zl_solve (UMFPACK_At, NULL, NULL, NULL, NULL, (double *) (x->values_cpx+(i*x->ld)),NULL,  (double*) (work->values_cpx+(i*work->ld)),NULL, d->numeric[n], d->Control[n], d->Info[n]) ;
        }
        if ( ret != UMFPACK_OK) {
            MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
            umfpack_zl_report_info(d->Control[n], d->Info[n]);
            return MESS_ERROR_GENERAL;
        }
    }

    if ( conv == 0) {
        mess_matrix_clear(&work);
    }
    return 0;
}


/**
 * @brief Clear an umfpack multisolver.
 * @param[in] data  input pointer to the data structure
 *@return always zero
 *
 * The @ref multiumf_clear function clears an umfpack multisolver object.
 *
 */
static int multiumf_clear(void *data){
    struct multi_umfpack *mlu = (struct multi_umfpack *) data;
    mess_int_t i;
    if ( mlu == NULL) return 0;
    for ( i = 0; i < mlu->nlu; i++){
        if ( mlu->numeric[i] != NULL ) {
            if ( mlu->datatypes[i] == MESS_COMPLEX ){
                umfpack_zl_free_numeric(&(mlu->numeric[i]));
            }else{
                umfpack_dl_free_numeric(&(mlu->numeric[i]));
            }
        }
        mess_free(mlu->Control[i]);
        mess_free(mlu->Info[i]);
    }
    mess_free(mlu->numeric);
    mess_free(mlu->Control);
    mess_free(mlu->Info);
    mess_free(mlu->datatypes);
    mess_free(mlu);
    return 0;
}


/**
 * @brief Create a multisolver on base of @umfpack.
 * @param[in] matrix  input system matrix
 * @param[in] shiftsl  input left shifts \f$ s \f$, in front of \f$ A \f$, @c NULL is not given
 * @param[in] shiftsr  input right shifts \f$ p \f$, in front of \f$ E \f$
 * @param[out] mlu  multidirect solver to create
 * @param[in] base   input existing solver to get the pattern if it is possible
 * @param[in] shiftmatrix_input matrix to shift if it points to @c NULL, the identity matrix is used
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_multidirect_create_umfpack function creates a multisolver from @umfpack.
 *
 */
int mess_multidirect_create_umfpack(mess_matrix matrix, mess_vector shiftsl, mess_vector shiftsr , mess_multidirect mlu, mess_direct base, mess_matrix shiftmatrix_input) {
    MSG_FNAME(__func__);
    mess_matrix work = NULL;
    mess_int_t *dpos;
    mess_int_t i, j;
    struct multi_umfpack *data;
    mess_int_t nshifts = 0;
    int use_gernalized=0;
    mess_matrix shiftmatrix;
    int ret = 0;
    double eps = mess_eps();
    mess_vector shiftsl_real = NULL, shiftsl_complex= NULL;
    mess_vector shiftsr_real = NULL, shiftsr_complex= NULL;

    int use_complex  = 0;
    int havelshift =0 ;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer (matrix);
    mess_check_nullpointer (shiftsr);
    mess_check_nullpointer (mlu);
    mess_check_square(matrix);


    nshifts = shiftsr->dim;
    if ( nshifts <=0 ) {
        MSG_ERROR("nshifts hasn't a useful value ( = " MESS_PRINTF_INT " )\n", nshifts);
        return MESS_ERROR_ARGUMENTS;
    }
    if ( shiftmatrix_input != NULL ) {
        MSG_INFO("using mass matrix to shift\n");
        mess_check_same_size(shiftmatrix_input, matrix);
        use_gernalized = 1;
    }
    if (shiftsl != NULL ) {
        if ( shiftsl->dim != shiftsr->dim) {
            MSG_ERROR("The left and right shift vector must have the same dimension. shiftsl->dim=" MESS_PRINTF_INT ", \t shiftsr->dim=" MESS_PRINTF_INT "\n", shiftsl->dim, shiftsr->dim);
            return MESS_ERROR_DIMENSION;
        }
        havelshift = 1;
    }


    mess_try_alloc(data, struct multi_umfpack*, sizeof (struct multi_umfpack));


    /*-----------------------------------------------------------------------------
     *  create input data
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&work);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    if ( MESS_IS_CSC(matrix)) {
        ret = mess_matrix_copy(matrix, work);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    } else {
        ret = mess_matrix_convert(matrix, work, MESS_CSC);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    }
    if ( use_gernalized ) {
        ret = mess_matrix_init(&shiftmatrix);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        if ( MESS_IS_CSC(shiftmatrix_input)) {
            ret = mess_matrix_copy(shiftmatrix_input, shiftmatrix);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        } else {
            ret = mess_matrix_convert(shiftmatrix_input, shiftmatrix, MESS_CSC);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
        }
    }
    // real and complex shifts
    MESS_INIT_VECTORS(&shiftsr_real,&shiftsr_complex);
    ret = mess_vector_alloc(shiftsr_real, nshifts, MESS_REAL);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_vector_alloc);
    ret = mess_vector_alloc(shiftsr_complex, nshifts, MESS_COMPLEX);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_vector_alloc);
    ret = mess_vector_copy(shiftsr, shiftsr_real);                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy);
    ret = mess_vector_copy(shiftsr, shiftsr_complex);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy);
    ret = mess_vector_toreal_nowarn(shiftsr_real);                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal_nowarn);
    ret = mess_vector_tocomplex(shiftsr_complex);                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
    if ( havelshift ) {
        MESS_INIT_VECTORS(&shiftsl_real,&shiftsl_complex);
        ret = mess_vector_alloc(shiftsl_real, nshifts, MESS_REAL);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_vector_alloc);
        ret = mess_vector_alloc(shiftsl_complex, nshifts, MESS_COMPLEX);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_vector_alloc);
        ret = mess_vector_copy(shiftsl, shiftsl_real);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy);
        ret = mess_vector_copy(shiftsl, shiftsl_complex);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy);
        ret = mess_vector_toreal_nowarn(shiftsl_real);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal_nowarn);
        ret = mess_vector_tocomplex(shiftsl_complex);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
    }

    /*-----------------------------------------------------------------------------
     *  select datatype
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data->datatypes, unsigned short*, sizeof(unsigned short)*nshifts);
    if (MESS_IS_REAL(work) /*  && MESS_IS_REAL(shiftmatrix) */ ){
        if ( havelshift ) {
            for (i = 0 ; i < nshifts; i++)  {
                if ( IS_REAL(shiftsr_complex->values_cpx[i]) && IS_REAL(shiftsl_complex->values_cpx[i])){
                    MSG_INFO("shiftr[" MESS_PRINTF_INT "] = %lg + %lgi, shiftl[" MESS_PRINTF_INT "] = %lg = %lgi are real\n",i, creal(shiftsr_complex->values_cpx[i]), cimag(shiftsr_complex->values_cpx[i]),i , creal(shiftsl_complex->values_cpx[i]), cimag(shiftsl_complex->values_cpx[i]));
                    data->datatypes[i]=MESS_REAL;
                } else {
                    MSG_INFO("shiftr[" MESS_PRINTF_INT "] = %lg + %lgi, shiftl[" MESS_PRINTF_INT "] = %lg = %lgi are complex\n",i, creal(shiftsr_complex->values_cpx[i]), cimag(shiftsr_complex->values_cpx[i]), i, creal(shiftsl_complex->values_cpx[i]), cimag(shiftsl_complex->values_cpx[i]));
                    data->datatypes[i]=MESS_COMPLEX;
                }
            }
        } else {
            for (i = 0 ; i < nshifts; i++)  {
                if ( IS_REAL(shiftsr_complex->values_cpx[i])){
                    MSG_INFO("shift[" MESS_PRINTF_INT "] = %lg + %lgi is real\n",i, creal(shiftsr_complex->values_cpx[i]), cimag(shiftsr_complex->values_cpx[i]));
                    data->datatypes[i]=MESS_REAL;
                } else {
                    MSG_INFO("shift[" MESS_PRINTF_INT "] = %lg + %lgi is complex\n",i, creal(shiftsr_complex->values_cpx[i]), cimag(shiftsr_complex->values_cpx[i]));
                    data->datatypes[i]=MESS_COMPLEX;
                }
            }
        }
        use_complex = 0;
    } else {
        MSG_INFO("all shifts are complex\n");
        for (i=0;i<nshifts;i++) data->datatypes[i]=MESS_COMPLEX;
        use_complex = 1;
    }

    data->cpx = use_complex;
    /*-----------------------------------------------------------------------------
     *  check if all diagonal elements exist
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_sort(work);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);

    mess_try_alloc(dpos, mess_int_t *, sizeof ( mess_int_t) * work->rows);
    for (i = 0; i < work->rows; i++) {
        dpos[i]=-1;
    }
    ret =   mess_matrix_diagpos(work , dpos);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_diagpos);

    /*-----------------------------------------------------------------------------
     *  blow up the matrices
     *-----------------------------------------------------------------------------*/
    if ( !use_gernalized ){
        int problems = 0;
        for ( i = 0; i < work->rows; i++) {
            if ( dpos[i] < 0) {
                MSG_WARN("diagonal entry in line " MESS_PRINTF_INT " doesn't exist.\n", i);
                problems ++;
            }
        }
        if ( problems > 0 ) {
            MSG_WARN("some diagonal entries don't exists. creating identity mass matrix.\n");
            ret = mess_matrix_init(&shiftmatrix);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
            if ( use_complex == 1) {
                ret = mess_matrix_eyec(shiftmatrix, work->rows, work->cols, MESS_CSC);
            } else {
                ret = mess_matrix_eye(shiftmatrix, work->rows, work->cols, MESS_CSC);
            }
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
            use_gernalized = 1;
        }
    }

    if (use_gernalized ) {
        ret = mess_matrix_add(0.0,shiftmatrix, 1, work);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_add(0,work, 1, shiftmatrix);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_sort(shiftmatrix);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
        ret = mess_matrix_sort(work);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
        ret = mess_matrix_diagpos(work , dpos);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_diagpos);
    }

    mess_try_alloc(data->numeric, void **,sizeof ( void *) * nshifts);


    /*-----------------------------------------------------------------------------
     *  configure umfpack
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data->Control, double**, sizeof ( double *) * nshifts);
    mess_try_alloc(data->Info, double ** , sizeof(double*)*nshifts);
    for ( i = 0; i < nshifts; i++) {
        mess_try_alloc(data->Control[i], double *, sizeof(double) * UMFPACK_CONTROL);
        mess_try_alloc(data->Info[i], double *, sizeof(double) * UMFPACK_INFO);
        if ( data->datatypes[i] == MESS_COMPLEX ){
            umfpack_zl_defaults(data->Control[i]);
        }else {
            umfpack_dl_defaults(data->Control[i]); // load defaults for UMFPACK
        }
        data->Control[i][UMFPACK_IRSTEP] = 0;   // disable iterative refinement
        //no pivoting
        //data->Control[i][UMFPACK_PIVOT_TOLERANCE] = 0.0;
    }




    /*-----------------------------------------------------------------------------
     *  Normal variant
     *-----------------------------------------------------------------------------*/

    // #pragma omp parallel for private(j,i) default(shared)
    for ( i = 0 ; i < nshifts; i++ ) {
        if (data->datatypes[i]==MESS_REAL) {
            void *SymbolicI;
            /*-----------------------------------------------------------------------------
             *  real system matrix and real shifts
             *-----------------------------------------------------------------------------*/
            double *vworkI;
            mess_try_alloc(vworkI, double *, work->nnz*sizeof(double));
            if (use_gernalized == 0) {
                memcpy(vworkI, work->values, work->nnz * sizeof(double));
                if (havelshift) {
                    for ( j =0; j < work->nnz; j++) vworkI[j] *= shiftsl_real->values[i];
                }
                for ( j = 0 ; j < work->rows; j++) {
                    if ( dpos[j] >= 0 ) {
                        vworkI[dpos[j]] += shiftsr_real->values[i];
                    }
                }

                // for ( j = 0 ; j < work->rows; j++)  vworkI[dpos[j]] += shiftsr_real->values[i];
            } else {
                if ( havelshift ) {
                    for ( j = 0; j < work->nnz; j++) vworkI[j] = shiftsl_real->values[i]*work->values [j] + shiftsr_real->values[i]*shiftmatrix->values[j];
                } else {
                    for ( j = 0; j < work->nnz; j++) vworkI[j] = work->values [j] + shiftsr_real->values[i]*shiftmatrix->values[j];
                }
            }

            ret = umfpack_dl_symbolic (work->rows, work->cols, work->colptr, work->rowptr, vworkI, &SymbolicI, data->Control[i], data->Info[i]) ;
            if ( ret != UMFPACK_OK) {
                MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
                umfpack_dl_report_info(data->Control[i], data->Info[i]);
                return MESS_ERROR_UMFPACK;
            }

            ret = umfpack_dl_numeric (work->colptr, work->rowptr, vworkI, SymbolicI, &(data->numeric[i]), data->Control[i], data->Info[i]) ;
            if ( ret != UMFPACK_OK) {
                MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
                umfpack_dl_report_info(data->Control[i], data->Info[i]);
                return MESS_ERROR_UMFPACK;
            }
            umfpack_dl_free_symbolic(&SymbolicI);
            mess_free(vworkI);
        }
        /*-----------------------------------------------------------------------------
         *  real system matrix and complex shifts
         *-----------------------------------------------------------------------------*/
        else if ( data->datatypes[i] == MESS_COMPLEX && use_complex == 0){
            mess_double_cpx_t *vwork;
            register mess_double_cpx_t s = shiftsr_complex->values_cpx[i];
            register mess_double_cpx_t t;
            void *Symbolic;
            mess_try_alloc(vwork, mess_double_cpx_t *, work->nnz*sizeof(mess_double_cpx_t));

            //printf("use_gernalized: %d\n", use_gernalized );
            if (use_gernalized == 0) {
                if (havelshift) {
                    t = shiftsl_complex->values_cpx[i];
                    for ( j =0; j < work->nnz; j++) vwork[j] = t*work->values[j];
                } else {
                    for ( j = 0; j < work->nnz; j++) vwork[j] = work->values[j];
                }
                for ( j = 0 ; j < work->rows; j++) {
                    if ( dpos[j] >= 0 ) {
                        vwork[dpos[j]] += s;
                    }
                    /* else {
                       MSG_PRINT("Zero entry on diagonal: %d\n",i);
                       } */
                }
            } else {
                if ( havelshift) {
                    t = shiftsl_complex->values_cpx[i];
                    for ( j = 0; j < work->nnz; j++) vwork[j] = t * work->values [j] + s *shiftmatrix->values[j];
                } else {
                    for ( j = 0; j < work->nnz; j++) vwork[j] = work->values [j] + s *shiftmatrix->values[j];
                }
            }


            ret = umfpack_zl_symbolic (work->rows, work->cols, work->colptr, work->rowptr, (double *)vwork,NULL,  &Symbolic, data->Control[i], data->Info[i]) ;
            if ( ret != UMFPACK_OK) {
                MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
                umfpack_zl_report_info(data->Control[i], data->Info[i]);
                return MESS_ERROR_UMFPACK;
            }
            ret = umfpack_zl_numeric (work->colptr, work->rowptr, (double *)vwork,NULL, Symbolic, &(data->numeric[i]), data->Control[i], data->Info[i]) ;
            if ( ret != UMFPACK_OK) {
                MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
                umfpack_zl_report_info(data->Control[i], data->Info[i]);
                return MESS_ERROR_UMFPACK;
            }
            umfpack_zl_free_symbolic(&Symbolic);
            mess_free(vwork);
        }
        /*-----------------------------------------------------------------------------
         *  complex system matrix
         *-----------------------------------------------------------------------------*/
        else {
            mess_double_cpx_t *vwork;
            register mess_double_cpx_t s = shiftsr_complex->values_cpx[i];
            register mess_double_cpx_t t;
            void *Symbolic;
            mess_try_alloc(vwork, mess_double_cpx_t *, work->nnz* sizeof(mess_double_cpx_t));

            if (use_gernalized == 0) {
                memcpy(vwork, work->values_cpx, work->nnz * sizeof(mess_double_cpx_t));
                if ( havelshift ) {
                    t =  shiftsl_complex->values_cpx[i];
                    for ( j = 0 ; j < work->nnz ; j++ ) vwork[j] *= t;
                }
                for ( j = 0 ; j < work->rows; j++)  vwork[dpos[j]] += s;
            } else {
                if ( havelshift ) {
                    t = shiftsl_complex->values_cpx[i];
                    for ( j = 0; j < work->nnz; j++) vwork[j] = t*work->values_cpx [j] + s *shiftmatrix->values_cpx[j];
                } else {
                    for ( j = 0; j < work->nnz; j++) vwork[j] = work->values_cpx [j] + s*shiftmatrix->values_cpx[j];
                }
            }


            ret = umfpack_zl_symbolic (work->rows, work->cols, work->colptr, work->rowptr, (double *)vwork,NULL,  &Symbolic, data->Control[i], data->Info[i]) ;
            if ( ret != UMFPACK_OK) {
                MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
                umfpack_zl_report_info(data->Control[i], data->Info[i]);
                return MESS_ERROR_UMFPACK;
            }
            ret = umfpack_zl_numeric (work->colptr, work->rowptr, (double *)vwork,NULL, Symbolic,  &(data->numeric[i]), data->Control[i], data->Info[i]) ;
            if ( ret != UMFPACK_OK) {
                MSG_ERROR("UMFPACK caused and error, err = %d\n", ret);
                umfpack_zl_report_info(data->Control[i], data->Info[i]);
                return MESS_ERROR_UMFPACK;
            }
            umfpack_zl_free_symbolic(&Symbolic);
            mess_free(vwork);
        }
    }
    mess_free(dpos);
    mlu->data_type = MESS_COMPLEX;
    mlu->indx = nshifts;
    data->nlu = nshifts;
    data->rows = mlu->rows = work->rows;
    data->cols = mlu->cols = work->cols;
    mlu->data = data;
    mlu->solve = multiumf_solve;
    mlu->solvet = multiumf_solvet;
    mlu->solveh = multiumf_solveh;
    mlu->solvem = multiumf_solvem;
    mlu->solvemt = multiumf_solvemt;
    mlu->solvemh = multiumf_solvemh;
    mlu->clear = multiumf_clear;

    mess_try_alloc(mlu->name, char*, sizeof(char)*20);
    strcpy(mlu->name, "UMFPACK MULTILU");

    mess_matrix_clear(&work);
    if (use_gernalized ) mess_matrix_clear(&shiftmatrix);
    mess_vector_clear(&shiftsr_real);
    mess_vector_clear(&shiftsr_complex);
    if ( havelshift ) {
        mess_vector_clear(&shiftsl_real);
        mess_vector_clear(&shiftsl_complex);
    }

    return 0;
}


