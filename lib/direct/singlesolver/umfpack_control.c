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
 * @file lib/direct/singlesolver/umfpack_control.c
 * @brief Interface to @umfpack.
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#ifdef MESS_HAVE_UMFPACK
#include <umfpack.h>

#ifdef MESS64
#define umfpack_dl_solve    umfpack_dl_solve
#define umfpack_zl_solve    umfpack_zl_solve
#define umfpack_zl_load_numeric umfpack_zl_load_numeric
#define umfpack_dl_load_numeric umfpack_dl_load_numeric
#define umfpack_zl_defaults     umfpack_zl_defaults
#define umfpack_dl_defaults     umfpack_dl_defaults
#define umfpack_dl_save_numeric umfpack_dl_save_numeric
#define umfpack_zl_save_numeric umfpack_zl_save_numeric
#define umfpack_zl_free_numeric umfpack_zl_free_numeric
#define umfpack_dl_free_numeric umfpack_dl_free_numeric
#define umfpack_dl_symbolic umfpack_dl_symbolic
#define umfpack_dl_numeric  umfpack_dl_numeric
#define umfpack_zl_symbolic     umfpack_zl_symbolic
#define umfpack_zl_numeric  umfpack_zl_numeric
#define umfpack_dl_report_info  umfpack_dl_report_info
#define umfpack_dl_free_symbolic umfpack_dl_free_symbolic
#define umfpack_zl_report_info   umfpack_zl_report_info
#define umfpack_zl_free_symbolic umfpack_zl_free_symbolic
#define umfpack_zl_get_lunz  umfpack_zl_get_lunz
#define umfpack_zl_get_numeric   umfpack_zl_get_numeric
#define umfpack_dl_get_lunz      umfpack_dl_get_lunz
#define umfpack_dl_get_numeric   umfpack_dl_get_numeric
#define umfpack_dl_get_determinant umfpack_dl_get_determinant
#define umfpack_zl_get_determinant umfpack_zl_get_determinant
#else
#define umfpack_dl_solve    umfpack_di_solve
#define umfpack_zl_solve    umfpack_zi_solve
#define umfpack_zl_load_numeric umfpack_zi_load_numeric
#define umfpack_dl_load_numeric umfpack_di_load_numeric
#define umfpack_zl_defaults     umfpack_zi_defaults
#define umfpack_dl_defaults     umfpack_di_defaults
#define umfpack_dl_save_numeric umfpack_di_save_numeric
#define umfpack_zl_save_numeric umfpack_zi_save_numeric
#define umfpack_zl_free_numeric umfpack_zi_free_numeric
#define umfpack_dl_free_numeric umfpack_di_free_numeric
#define umfpack_dl_symbolic umfpack_di_symbolic
#define umfpack_dl_numeric  umfpack_di_numeric
#define umfpack_zl_symbolic     umfpack_zi_symbolic
#define umfpack_zl_numeric  umfpack_zi_numeric
#define umfpack_dl_report_info  umfpack_di_report_info
#define umfpack_dl_free_symbolic umfpack_di_free_symbolic
#define umfpack_zl_report_info   umfpack_zi_report_info
#define umfpack_zl_free_symbolic umfpack_zi_free_symbolic
#define umfpack_zl_get_lunz  umfpack_zi_get_lunz
#define umfpack_zl_get_numeric   umfpack_zi_get_numeric
#define umfpack_dl_get_lunz      umfpack_di_get_lunz
#define umfpack_dl_get_numeric   umfpack_di_get_numeric
#define umfpack_dl_get_determinant umfpack_di_get_determinant
#define umfpack_zl_get_determinant umfpack_zi_get_determinant
#endif


#define UMFPACK_CONTROL_LOG2_10 3.321928094887362347870319429489390175864831393024

typedef struct umfpack_control_solver {
    void *Numeric;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    unsigned short cpx;
    mess_matrix work;
    mess_int_t dim;
}umfpack_control_solver;


/**
 * @internal
 * @brief Get the L factor of a matrix.
 * @param[in] data      input pointer to the data object
 * @param[in,out] L     input/output the L factor
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref umfpack_control_getL function gets the L factor of a matrix decomposed with @umfpack. L is in CSR format.
 *
 * @attention Internal use only.
 */
static int umfpack_control_getL(void *data, mess_matrix L ){
    MSG_FNAME(__func__);
    umfpack_control_solver *sol = (umfpack_control_solver *) data;
    mess_int_t lnz=0, unz = 0,rows =0 , cols = 0, nzdiag = 0;
    mess_int_t reci = 0 ;
    int ret = 0;
    mess_check_nullpointer(L);
    mess_check_nullpointer(data);

    if ( sol ->cpx) {
        // get L directly
        ret = umfpack_zl_get_lunz(&lnz, &unz, &rows, &cols, &nzdiag, sol->Numeric);             FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_zl_get_lunz);
        ret = mess_matrix_alloc(L, rows, cols, lnz, MESS_CSR, MESS_COMPLEX);            FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), mess_matrix_alloc);
        //          Lrow,   Lcol        LValues,   Ucol,  Urow, Uval, p,   q,    diag, recip, Scale, NUmeric)
        ret = umfpack_zl_get_numeric(L->rowptr, L->colptr,(double *) L->values_cpx, NULL,
                NULL, NULL, NULL, NULL,
                NULL,NULL, NULL, NULL,
                &reci, NULL, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_zl_get_numeric);
    } else {
        // get L directly
        ret = umfpack_dl_get_lunz(&lnz, &unz, &rows, &cols, &nzdiag, sol->Numeric);     FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_dl_get_lunz);
        ret = mess_matrix_alloc(L, rows, cols, lnz, MESS_CSR, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        //          Lrow,   Lcol        LValues,   Ucol,  Urow, Uval, p,   q,    diag, recip, Scale, NUmeric)
        ret = umfpack_dl_get_numeric(L->rowptr, L->colptr, L->values, NULL, NULL, NULL, NULL, NULL,NULL, &reci, NULL, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_dl_get_numeric);
    }
    return 0;
}



/**
 * @internal
 * @brief Get the U factor of a matrix.
 * @param[in] data  input pointer to the data object
 * @param[out] U output U matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref umfpack_control_getU function gets the U factor of a matrix decomposed with @umfpack. U is in CSR format.
 *
 * @attention Internal use only.
 */
static int umfpack_control_getU(void *data, mess_matrix U ){
    MSG_FNAME(__func__);
    umfpack_control_solver *sol = (umfpack_control_solver *) data;
    mess_matrix tmp;
    mess_int_t lnz=0, unz = 0,rows =0 , cols = 0, nzdiag = 0;
    mess_int_t reci = 0 ;
    int ret = 0;
    mess_check_nullpointer(U);
    mess_check_nullpointer(data);
    if ( sol-> cpx ) {
        ret = mess_matrix_init(&tmp);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = umfpack_zl_get_lunz(&lnz, &unz, &rows, &cols, &nzdiag, sol->Numeric); FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_zl_get_lunz);
        ret = mess_matrix_alloc(tmp, rows, cols, unz, MESS_CSC, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = umfpack_zl_get_numeric(NULL, NULL, NULL, NULL,
                tmp->colptr, tmp->rowptr, (double *) tmp->values_cpx, NULL,
                NULL, NULL, NULL,NULL,
                &reci, NULL, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_zl_get_numeric);
        ret = mess_matrix_convert(tmp, U, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
        mess_matrix_clear(&tmp);
    } else {
        ret = mess_matrix_init(&tmp);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = umfpack_dl_get_lunz(&lnz, &unz, &rows, &cols, &nzdiag, sol->Numeric); FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_dl_get_lunz);
        ret = mess_matrix_alloc(tmp, rows, cols, unz, MESS_CSC, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = umfpack_dl_get_numeric(NULL, NULL, NULL, tmp->colptr, tmp->rowptr, tmp->values, NULL, NULL,NULL, &reci, NULL, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_dl_get_numeric);
        ret = mess_matrix_convert(tmp, U, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
        mess_matrix_clear(&tmp);
    }
    return 0;

}


/**
 * @internal
 * @brief Get the row permutation of the matrix.
 * @param[in] data  input pointer to the data object
 * @param[in,out] p output permutation
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref umfpack_control_getpermp function gets the row perumation \f$ P \f$ for \f$ PRAQ=LU \f$.
 *
 * @attention Internal use only.
 */
static int umfpack_control_getpermp(void *data, mess_int_t *p){
    MSG_FNAME(__func__);
    umfpack_control_solver *sol = (umfpack_control_solver *) data;
    mess_int_t lnz=0, unz = 0,rows =0 , cols = 0, nzdiag = 0;
    mess_int_t reci = 0 ;
    int ret = 0;
    mess_check_nullpointer(p);
    mess_check_nullpointer(data);
    if (sol->cpx) {
        ret = umfpack_zl_get_lunz(&lnz, &unz, &rows, &cols, &nzdiag, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_zl_get_lunz);
        ret = umfpack_zl_get_numeric(NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL,
                p, NULL, NULL, NULL,
                &reci, NULL, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_zl_get_numeric);

    } else {
        ret = umfpack_dl_get_lunz(&lnz, &unz, &rows, &cols, &nzdiag, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_dl_get_lunz);
        ret = umfpack_dl_get_numeric(NULL, NULL, NULL, NULL, NULL, NULL,  p, NULL, NULL, &reci, NULL, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_dl_get_numeric);
    }
    return 0;
}

/**
 * @internal
 * @brief Get the column permutation of a matrix.
 * @param[in] data  input pointer to the data object
 * @param[in,out] q output permuataion
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref umfpack_control_getpermq function gets the colume permutation \f$ Q \f$ for \f$ PRAQ=LU \f$.
 *
 * @attention Internal use only.
 */
static int umfpack_control_getpermq(void *data, mess_int_t *q){
    MSG_FNAME(__func__);
    umfpack_control_solver *sol = (umfpack_control_solver *) data;
    mess_int_t lnz=0, unz = 0,rows =0 , cols = 0, nzdiag = 0;
    mess_int_t reci = 0 ;
    int ret =0;
    mess_check_nullpointer(q);
    mess_check_nullpointer(data);
    if (sol->cpx) {
        ret = umfpack_zl_get_lunz(&lnz, &unz, &rows, &cols, &nzdiag, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_zl_get_lunz);
        ret = umfpack_zl_get_numeric(NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL,
                NULL, q, NULL, NULL,
                &reci, NULL, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_zl_get_numeric);

    } else {
        ret = umfpack_dl_get_lunz(&lnz, &unz, &rows, &cols, &nzdiag, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_dl_get_lunz);
        ret = umfpack_dl_get_numeric(NULL, NULL, NULL, NULL, NULL, NULL,  NULL,q, NULL, &reci, NULL, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_dl_get_numeric);
    }
    return 0;
}

/**
 * @internal
 * @brief Get the row scaling of a matrix.
 * @param[in] data  input pointer to the data object
 * @param[in,out] r output permuataion
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref umfpack_control_getscalerow function gets the row scaling diagonal entries of \f$ R \f$ for \f$ PRAQ=LU \f$.
 *
 * @attention Internal use only.
 */
static int umfpack_control_getscalerow(void *data, mess_vector r){
    MSG_FNAME(__func__);
    umfpack_control_solver *sol = (umfpack_control_solver *) data;
    mess_int_t do_recip = 0, ret =0, i = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(r);


    /*-----------------------------------------------------------------------------
     *  get row scaling factors
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_toreal_nowarn(r);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal_nowarn);
    ret = mess_vector_resize(r,sol->dim);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);

    if (sol->cpx) {
        ret = umfpack_zl_get_numeric(NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL,
                &do_recip, r->values, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_zl_get_numeric);
    } else{
        ret = umfpack_dl_get_numeric(NULL, NULL, NULL,
                NULL, NULL, NULL,
                NULL, NULL, NULL,
                &do_recip, r->values, sol->Numeric);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=UMFPACK_OK), umfpack_dl_get_numeric);
    }

    /*-----------------------------------------------------------------------------
     *  check if reciprocal values are given
     *-----------------------------------------------------------------------------*/
    if(!do_recip){
        for(i=0;i<r->dim;++i){
            r->values[i]=1.0/r->values[i];
        }
    }


    return 0;
}


/**
 * @brief Solve \f$ Ax=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref umfpack_control_solve function solves \f$ Ax=b \f$.
 *
 */
static int umfpack_control_solve(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    umfpack_control_solver * sol = (umfpack_control_solver *)data;
    int ret = 0;
    mess_datatype_t dt;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->dim != sol->dim) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->dim, sol->dim);
        return MESS_ERROR_DIMENSION;
    }
    if ( x->dim != sol->dim){
        ret = mess_vector_resize(x, sol->dim); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }
    dt = b->data_type;

    if ( sol->cpx ) {
        ret = mess_vector_tocomplex(x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_tocomplex(b);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);

        ret= umfpack_zl_solve (UMFPACK_A, sol->work->colptr, sol->work->rowptr, (double *) sol->work->values_cpx, NULL,
                (double *) x->values_cpx, NULL, (double *)  b->values_cpx,NULL,  sol->Numeric, sol->Control, sol->Info);
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }
    } else if (!sol->cpx && MESS_IS_COMPLEX(b)){
        // MSG_WARN("Applying a real solver to a complex righthand side is inefficient. \n");
        mess_vector br,xr;
        mess_vector bi,xi;
        MESS_INIT_VECTORS(&br,&xr,&bi,&xi);
        ret = mess_vector_alloc(xr, sol->dim, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(br, sol->dim, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(xi, sol->dim, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(bi, sol->dim, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_realpart(b, br);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_realpart);
        ret = mess_vector_imagpart(b, bi);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_imagpart);
        ret= umfpack_dl_solve (UMFPACK_A, sol->work->colptr, sol->work->rowptr, sol->work->values, xr->values, br->values, sol->Numeric, sol->Control, sol->Info) ;
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }
        ret= umfpack_dl_solve (UMFPACK_A, sol->work->colptr, sol->work->rowptr, sol->work->values, xi->values, bi->values, sol->Numeric, sol->Control, sol->Info) ;
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }
        ret = mess_vector_complex_from_parts(xr,xi,x);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_complex_from_parts);
        mess_vector_clear(&xr);
        mess_vector_clear(&br);
        mess_vector_clear(&bi);
        mess_vector_clear(&xi);
    } else {
        ret = mess_vector_toreal_nowarn(x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        ret= umfpack_dl_solve (UMFPACK_A, sol->work->colptr, sol->work->rowptr, sol->work->values, x->values, b->values, sol->Numeric, sol->Control, sol->Info) ;
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }
    }

    if ( dt != b->data_type) {
        ret = mess_vector_totype(b, dt);    FUNCTION_FAILURE_HANDLE(ret, (ret !=0), mess_vector_totype);
    }

    return 0;
}

/**
 * @brief Solve \f$ A^Tx=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref umfpack_control_solvet function solves \f$ A^Tx=b \f$.
 *
 */
static int umfpack_control_solvet(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    umfpack_control_solver * sol = (umfpack_control_solver *)data;
    int ret = 0;
    mess_datatype_t dt;
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->dim != sol->dim) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->dim, sol->dim);
        return MESS_ERROR_DIMENSION;
    }
    if ( x->dim != sol->dim){
        ret = mess_vector_resize(x, sol->dim); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }
    dt = b->data_type;

    if ( sol->cpx ) {
        ret = mess_vector_tocomplex(x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_tocomplex(b);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret= umfpack_zl_solve (UMFPACK_Aat, sol->work->colptr, sol->work->rowptr, (double *) sol->work->values_cpx, NULL,
                (double *) x->values_cpx, NULL, (double *)  b->values_cpx,NULL,  sol->Numeric, sol->Control, sol->Info) ;
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_zl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

    } else if (!sol->cpx && MESS_IS_COMPLEX(b)){
        // MSG_WARN("Applying a real solver to a complex righthand side is inefficient. \n");
        mess_vector br,xr;
        mess_vector bi,xi;
        MESS_INIT_VECTORS(&br,&bi,&xi,&xr);
        ret = mess_vector_alloc(xr, sol->dim, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(br, sol->dim, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(xi, sol->dim, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(bi, sol->dim, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_realpart(b, br);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_realpart);
        ret = mess_vector_imagpart(b, bi);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_imagpart);
        ret= umfpack_dl_solve (UMFPACK_Aat, sol->work->colptr, sol->work->rowptr, sol->work->values, xr->values, br->values, sol->Numeric, sol->Control, sol->Info) ;
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

        ret= umfpack_dl_solve (UMFPACK_Aat, sol->work->colptr, sol->work->rowptr, sol->work->values, xi->values, bi->values, sol->Numeric, sol->Control, sol->Info) ;
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

        ret = mess_vector_complex_from_parts(xr,xi,x);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_complex_from_parts);
        mess_vector_clear(&xr);
        mess_vector_clear(&br);
        mess_vector_clear(&bi);
        mess_vector_clear(&xi);
    } else {
        ret = mess_vector_toreal_nowarn(x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        ret= umfpack_dl_solve (UMFPACK_Aat, sol->work->colptr, sol->work->rowptr, sol->work->values, x->values, b->values, sol->Numeric, sol->Control, sol->Info) ;
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

    }

    if ( dt != b->data_type) {
        ret = mess_vector_totype(b, dt);    FUNCTION_FAILURE_HANDLE(ret, (ret !=0), mess_vector_totype);
    }

    return 0;

}

/**
 * @brief Solve \f$ A^Hx=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref umfpack_control_solveh function solves \f$ A^Hx=b \f$.
 *
 */
static int umfpack_control_solveh(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    umfpack_control_solver * sol = (umfpack_control_solver *)data;
    int ret = 0;
    mess_datatype_t dt;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->dim != sol->dim) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->dim, sol->dim);
        return MESS_ERROR_DIMENSION;
    }
    if ( x->dim != sol->dim){
        ret = mess_vector_resize(x, sol->dim); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }
    dt = b->data_type;

    if ( sol->cpx) {
        ret = mess_vector_tocomplex(x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_tocomplex(b);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret= umfpack_zl_solve (UMFPACK_At, sol->work->colptr, sol->work->rowptr, (double *) sol->work->values_cpx, NULL,
                (double *) x->values_cpx, NULL, (double *)  b->values_cpx,NULL,  sol->Numeric, sol->Control, sol->Info) ;
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_zl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

    } else if ( ! sol->cpx && MESS_IS_COMPLEX(b)){
        // MSG_WARN("Applying a real solver to a complex righthand side is inefficient. \n");
        mess_vector br,xr;
        mess_vector bi,xi;
        MESS_INIT_VECTORS(&br,&bi,&xr,&xi);
        ret = mess_vector_alloc(xr, sol->dim, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(br, sol->dim, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(xi, sol->dim, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(bi, sol->dim, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_realpart(b, br);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_realpart);
        ret = mess_vector_imagpart(b, bi);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_imagpart);
        ret= umfpack_dl_solve (UMFPACK_At, sol->work->colptr, sol->work->rowptr, sol->work->values, xr->values, br->values, sol->Numeric, sol->Control, sol->Info) ;
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

        ret= umfpack_dl_solve (UMFPACK_At, sol->work->colptr, sol->work->rowptr, sol->work->values, xi->values, bi->values, sol->Numeric, sol->Control, sol->Info) ;
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

        ret = mess_vector_complex_from_parts(xr,xi,x);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_complex_from_parts);
        mess_vector_clear(&xr);
        mess_vector_clear(&br);
        mess_vector_clear(&bi);
        mess_vector_clear(&xi);

    } else {
        ret = mess_vector_toreal_nowarn(x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        ret= umfpack_dl_solve (UMFPACK_At, sol->work->colptr, sol->work->rowptr, sol->work->values, x->values, b->values, sol->Numeric, sol->Control, sol->Info) ;
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

    }
    if ( dt != b->data_type) {
        ret = mess_vector_totype(b, dt);    FUNCTION_FAILURE_HANDLE(ret, (ret !=0), mess_vector_totype);
    }

    return 0;
}

/**
 * @brief Solve\f$ AX=B \f$ (matrix version)
 * @param [in] data pointer to the data object
 * @param [in] b right hand side B
 * @param [in,out] x solution X
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref umfpack_control_solvem function solves \f$ AX=B \f$ with \f$ B \f$ and \f$ X \f$ matrices.
 *
 */
static int umfpack_control_solvem(void *data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    int conv;
    mess_matrix work;
    mess_int_t i;
    int ret = 0;
    umfpack_control_solver * sol = (umfpack_control_solver *)data;
    mess_datatype_t dt;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->rows != sol->dim) {
        MSG_ERROR("b has the wrong dimension (b->rows = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->rows, sol->dim);
        return MESS_ERROR_DIMENSION;
    }


    dt = b->data_type;
    MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
    MESS_MATRIX_RESET(x);
    ret = mess_matrix_alloc(x, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, (sol->cpx)?MESS_COMPLEX:MESS_REAL );
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    if ( sol->cpx) {
        ret = mess_matrix_tocomplex(work);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
#ifdef _OPENMP
#pragma omp parallel for private (i) default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret= umfpack_zl_solve (UMFPACK_A, sol->work->colptr, sol->work->rowptr, (double *)sol->work->values_cpx, NULL,
                    (double *) (x->values_cpx+i*x->ld), NULL,
                    (double *)  (work->values_cpx+i*work->ld),NULL,  sol->Numeric, sol->Control, sol->Info) ;

        }
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_zl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }
    } else if ( sol->cpx == 0 && MESS_IS_COMPLEX(b)){
        mess_matrix xr, xi, br, bi;
        ret = mess_matrix_init(&xr);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&xi);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&br);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&bi);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_realpart(work,br);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realpart);
        ret = mess_matrix_imagpart(work,bi);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_imagpart);

        ret = mess_matrix_alloc(xr, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(xi, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
#ifdef _OPENMP
#pragma omp parallel for private (i) default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret= umfpack_dl_solve (UMFPACK_A, sol->work->colptr, sol->work->rowptr, sol->work->values,
                    (xr->values+i*xr->ld), (br->values+i*br->ld), sol->Numeric, sol->Control, sol->Info) ;

            ret= umfpack_dl_solve (UMFPACK_A, sol->work->colptr, sol->work->rowptr, sol->work->values,
                    (xi->values+i*xi->ld), (bi->values+i*bi->ld), sol->Numeric, sol->Control, sol->Info) ;

        }
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }
        ret = mess_matrix_complex_from_parts(xr,xi,x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_complex_from_parts);
        mess_matrix_clear(&bi);
        mess_matrix_clear(&br);
        mess_matrix_clear(&xi);
        mess_matrix_clear(&xr);
    } else {
#ifdef _OPENMP
#pragma omp parallel for private (i) default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret= umfpack_dl_solve (UMFPACK_A, sol->work->colptr, sol->work->rowptr, sol->work->values,
                    (x->values+i*x->ld), (work->values+i*work->ld), sol->Numeric, sol->Control, sol->Info) ;

        }
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }
    }

    if ( conv == 0) {
        mess_matrix_clear(&work);
    }
    if ( dt != b->data_type) {
        ret = mess_matrix_totype(b, dt);    FUNCTION_FAILURE_HANDLE(ret, (ret !=0), mess_matrix_totype);
    }

    return ret;
}

/**
 * @brief Solve \f$ A^TX=B \f$ (matrix version)
 * @param [in] data pointer to the data object
 * @param [in] b right hand side B
 * @param [in,out] x solution X
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref umfpack_control_solvemt function solves \f$ A^T X=B \f$ with \f$ B \f$ and \f$ X \f$ matrices.
 *
 */
static int umfpack_control_solvemt(void *data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    int conv;
    mess_matrix work;
    mess_int_t i;
    int ret = 0;
    umfpack_control_solver * sol = (umfpack_control_solver *)data;
    mess_datatype_t dt;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->rows != sol->dim) {
        MSG_ERROR("b has the wrong dimension (b->rows = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->rows, sol->dim);
        return MESS_ERROR_DIMENSION;
    }

    dt = b->data_type;
    dt = b->data_type;
    MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
    MESS_MATRIX_RESET(x);
    ret = mess_matrix_alloc(x, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, (sol->cpx)?MESS_COMPLEX:MESS_REAL  );
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);


    if ( sol->cpx) {
        ret = mess_matrix_tocomplex(work);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
#ifdef _OPENMP
#pragma omp parallel for private (i) default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret= umfpack_zl_solve (UMFPACK_Aat, sol->work->colptr, sol->work->rowptr, (double*)sol->work->values_cpx, NULL, (double *) (x->values_cpx+i*x->ld), NULL,
                    (double *)  (work->values_cpx+i*work->ld),NULL,  sol->Numeric, sol->Control, sol->Info) ;

        }
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_zl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }
    } else if ( sol->cpx == 0 && MESS_IS_COMPLEX(b)){
        mess_matrix xr, xi, br, bi;
        ret = mess_matrix_init(&xr);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&xi);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&br);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&bi);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_realpart(work,br);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realpart);
        ret = mess_matrix_imagpart(work,bi);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_imagpart);
        ret = mess_matrix_alloc(xi, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, MESS_REAL );      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(xr, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, MESS_REAL );      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);


#ifdef _OPENMP
#pragma omp parallel for private (i) default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret= umfpack_dl_solve (UMFPACK_At, sol->work->colptr, sol->work->rowptr, sol->work->values,
                    (xr->values+i*xr->ld), (br->values+i*br->ld), sol->Numeric, sol->Control, sol->Info) ;

            ret= umfpack_dl_solve (UMFPACK_At, sol->work->colptr, sol->work->rowptr, sol->work->values,
                    (xi->values+i*xi->ld), (bi->values+i*bi->ld), sol->Numeric, sol->Control, sol->Info) ;

        }
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }
        ret = mess_matrix_complex_from_parts(xr,xi,x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_complex_from_parts);
        mess_matrix_clear(&bi);
        mess_matrix_clear(&br);
        mess_matrix_clear(&xi);
        mess_matrix_clear(&xr);

    } else {
#ifdef _OPENMP
#pragma omp parallel for private (i) default(shared) reduction(+:ret)
#endif
        for  ( i = 0; i < work->cols; i++) {
            ret= umfpack_dl_solve (UMFPACK_At, sol->work->colptr, sol->work->rowptr, sol->work->values,
                    (x->values+i*x->ld), (work->values+i*work->ld), sol->Numeric, sol->Control, sol->Info) ;

        }
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }
    }

    if ( conv == 0) {
        mess_matrix_clear(&work);
    }

    if ( dt != b->data_type) {
        ret = mess_matrix_totype(b, dt);    FUNCTION_FAILURE_HANDLE(ret, (ret !=0), mess_matrix_totype);
    }


    return ret;
}

/**
 * @brief Solve \f$ A^HX=B \f$ (matrix version)
 * @param [in] data pointer to the data object
 * @param [in] b right hand side B
 * @param [in,out] x solution X
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref umfpack_control_solvemh function solves \f$ A^H X=B \f$ with \f$ B \f$ and \f$ X \f$ matrices.
 *
 */
static int umfpack_control_solvemh(void *data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    int conv;
    mess_matrix work;
    mess_int_t i;
    int ret = 0;
    umfpack_control_solver * sol = (umfpack_control_solver *)data;
    mess_datatype_t dt;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->rows != sol->dim) {
        MSG_ERROR("b has the wrong dimension (b->rows = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->rows, sol->dim);
        return MESS_ERROR_DIMENSION;
    }
    if ( !sol->cpx){
        return umfpack_control_solvemt(data, b,x);
    }


    dt = b->data_type;
    MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
    MESS_MATRIX_RESET(x);
    ret = mess_matrix_alloc(x, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, (sol->cpx)?MESS_COMPLEX:MESS_REAL );   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_tocomplex(work);                                                                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);

#ifdef _OPENMP
#pragma omp parallel for private (i) default(shared) reduction(+:ret)
#endif
    for ( i = 0; i < work->cols; i++) {
        ret= umfpack_zl_solve (UMFPACK_At, sol->work->colptr, sol->work->rowptr, (double*)sol->work->values_cpx, NULL,
                (double *) (x->values_cpx+i*x->ld), NULL, (double *)  (work->values_cpx+i*work->ld),NULL,  sol->Numeric, sol->Control, sol->Info) ;

    }
    if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_zl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

    if ( conv == 0) {
        mess_matrix_clear(&work);
    }

    if ( dt != b->data_type) {
        ret = mess_matrix_totype(b, dt);    FUNCTION_FAILURE_HANDLE(ret, (ret !=0), mess_matrix_totype);
    }

    return ret;
}


/**
 * @brief Clear an umfpack_control_solver  object.
 * @param [in,out] solver pointer to the data object
 * @return always zero
 *
 * The @ref umfpack_control_clear function clears an @umfpack object.
 *
 */
static int umfpack_control_clear(void *solver){
    // MSG_FNAME(__func__);
    umfpack_control_solver * sol = (umfpack_control_solver*) solver;
    if ( sol != NULL) {
        if ( sol->Numeric != NULL)  {
            if ( sol->cpx )
                umfpack_zl_free_numeric(&(sol->Numeric));
            else
                umfpack_dl_free_numeric(&(sol->Numeric));
        }
        mess_matrix_clear(&(sol->work));
        mess_free( sol );
    }
    return 0;
}


/**
 * @brief Compute the inverse of a matrix
 * @param [in] data  pointer to the internal data.
 * @param [in,out] inv   output inverse
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref umfpack_control_inverse function computes \f$ A^{-1} \f$.
 *
 */
static int  umfpack_control_inverse ( void *data, mess_matrix inv ){
    MSG_FNAME(__func__);
    umfpack_control_solver * sol = (umfpack_control_solver*) data;
    mess_matrix eye;
    int ret = 0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(inv);
    MESS_MATRIX_RESET(inv);

    if ( sol-> cpx) {
        ret = mess_matrix_init(&eye);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_eyec(eye, sol->dim, sol->dim, MESS_DENSE);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
        ret = umfpack_control_solvem(data, eye, inv);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), umfpack_control_solvem);
        mess_matrix_clear(&eye);
    } else {
        ret = mess_matrix_init(&eye);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_eye(eye, sol->dim, sol->dim, MESS_DENSE);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
        ret = umfpack_control_solvem(data, eye, inv);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), umfpack_control_solvem);
        mess_matrix_clear(&eye);
    }

    return 0;
}       /* -----  end of function umfpack_inverse  ----- */

/**
 * @brief Compute the determinant of a matrix using @umfpack.
 * @param[in] data  input pointer to the internal data
 * @param[out] m output mantissa of determinant
 * @param[out] e output exponent of determinant
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref umfpack_control_det function computes \f$ \det( A ) \f$ using @umfpack.
 *
 */
static int umfpack_control_det( void *data, double * m, double *e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(data);
    mess_check_nullpointer(m);
    mess_check_nullpointer(e);
    int ret ;
    umfpack_control_solver * sol = (umfpack_control_solver*) data;
    ret = umfpack_dl_get_determinant (m,e,sol->Numeric,sol->Info);
    //*e *= log2(10);
    *e *= UMFPACK_CONTROL_LOG2_10;
    if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); }

    return 0;
}

/**
 * @brief Compute the determinant of a matrix using @umfpack.
 * @param[in] data  input pointer to the internal data
 * @param[out] mr output mantissa of real part of determinant
 * @param[out] mi output mantissa of imaginary part of determinant
 * @param[out] e output exponent of determinant
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref umfpack_control_detc function computes \f$ \det( A ) \f$ using @umfpack.
 *
 */
static int umfpack_control_detc( void *data, double * mr, double * mi, double *e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(data);
    mess_check_nullpointer(mr);
    mess_check_nullpointer(mi);
    mess_check_nullpointer(e);
    int ret;
    umfpack_control_solver * sol = (umfpack_control_solver*) data;
    ret = umfpack_zl_get_determinant (mr,mi,e,sol->Numeric,sol->Info);
    //*e *= log2(10);
    *e *= UMFPACK_CONTROL_LOG2_10;
    if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_zl_report_info(sol->Control, sol->Info); }
    return 0;
}






/**
 * @brief Generate a direct linear solver for standard linear systems \f$ Ax=b \f$ with @umfpack.
 * @param[in] matrix matrix to decompose
 * @param[in,out] solver direct solver
 * @param[in] Control umfpack control array
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref  mess_direct_create_umfpack_control function factorizes a matrix with @umfpack.
 * The Control argument provides the ability to control @umfpack solver options.
 *
 */
int mess_direct_create_umfpack_control(mess_matrix matrix, mess_direct solver, double* Control){
    MSG_FNAME(__func__);
    umfpack_control_solver * sol;
    int ret;
    void *Symbolic;
    mess_int_t i;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_square(matrix);
    mess_check_real_or_complex(matrix);

    /*-----------------------------------------------------------------------------
     *  init solver structure and umfpack controls
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(sol, umfpack_control_solver*, sizeof(umfpack_control_solver));
    MESS_INIT_MATRICES(&(sol->work));
    ret = mess_matrix_convert(matrix,(sol->work),MESS_CSC);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert);

    sol->cpx = MESS_IS_COMPLEX(sol->work);
    for(i=0;i<UMFPACK_CONTROL;++i){sol->Control[i]=Control[i];}

    // switch off stability scaling otherwise getLU is not working  sol->Control[UMFPACK_SCALE] = 0;

    //sol->Control[UMFPACK_SCALE] = 0;

    sol->dim = sol->work->rows;

    /*-----------------------------------------------------------------------------
     *  compute numeric and symbolic factorization
     *-----------------------------------------------------------------------------*/

    if ( sol->cpx ) {
        ret = umfpack_zl_symbolic (sol->work->rows, sol->work->cols, sol->work->colptr, sol->work->rowptr, (double *)  sol->work->values_cpx, NULL,
                &(Symbolic), sol->Control, sol->Info);
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_zl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

        ret = umfpack_zl_numeric (sol->work->colptr, sol->work->rowptr, (double *) sol->work->values_cpx, NULL,
                Symbolic, &(sol->Numeric), sol->Control, sol->Info);
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_zl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

        solver->data_type = MESS_COMPLEX;
    } else {
        ret = umfpack_dl_symbolic (sol->work->rows, sol->work->cols, sol->work->colptr,sol->work->rowptr, sol->work->values,
                &(Symbolic), sol->Control, sol->Info);
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

        ret = umfpack_dl_numeric (sol->work->colptr,sol->work->rowptr, sol->work->values, Symbolic, &(sol->Numeric), sol->Control, sol->Info) ;
        if (ret != UMFPACK_OK ) {MSG_ERROR("umfpack returned with %d\n",ret);umfpack_dl_report_info(sol->Control, sol->Info); return MESS_ERROR_UMFPACK; }

        solver->data_type = MESS_REAL;
    }

    if(sol->cpx){
        umfpack_zl_free_symbolic(&Symbolic);
    }else{
        umfpack_dl_free_symbolic(&Symbolic);
    }


    /*-----------------------------------------------------------------------------
     *  set solver functions
     *-----------------------------------------------------------------------------*/
    solver->rows =sol->work->rows;
    solver->cols =sol->work->cols;
    solver->data = (void *)sol;
    solver->solve = umfpack_control_solve;
    solver->solvet = umfpack_control_solvet;
    solver->solveh = umfpack_control_solveh;
    solver->solvem = umfpack_control_solvem;
    solver->solvemt = umfpack_control_solvemt;
    solver->solvemh = umfpack_control_solvemh;
    solver->clear = umfpack_control_clear;
    solver->getL = umfpack_control_getL;
    solver->getU = umfpack_control_getU;
    solver->getpermp = umfpack_control_getpermp;
    solver->getpermq = umfpack_control_getpermq;
    solver->getscalerow = umfpack_control_getscalerow;
    solver->inverse = umfpack_control_inverse;
    solver->det = umfpack_control_det;
    solver->detc = umfpack_control_detc;

    //SET_SOLVERNAME(solver->name, __func__);
    solver->name = NULL;
    return 0;
}


#endif


