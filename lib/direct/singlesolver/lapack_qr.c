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
 * @file lib/direct/singlesolver/lapack_qr.c
 * @brief   Generate a @lapack based solver based on QR decomposition.
 * @author @koehlerm
 * @author @mbehr
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/misc.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#ifdef _OPENMP_H
#include <omp.h>
#endif


/**
 * @internal
 * @brief Internal data structure.
 * Internal data structure for QR/LQ decompostitions.
 * @attention Internal use only.
 */
struct lapackqr_solver {
    double *val;                    /**< real value pointer of a matrix in fortran indexing style  */
    mess_double_cpx_t  *val_cpx;       /**< complex value pointer of a matrix in fortran indexing style. */
    double *tau;                    /**< real scalar factors of elementary reflectors of QR/LQ decomposition. */
    mess_double_cpx_t *tau_cpx;        /**< complex scalar factors of elementary reflectors of QR/LQ decomposition.*/
    mess_int_t rows;                /**< number of rows of the matrix decomposed to.*/
    mess_int_t cols;                /**< number of cols of the matrix decomposed to.*/
    unsigned short cpx;             /**< indicator for complex valued matrices*/
};


/**
 * @internal
 * @brief Clear a @lapack QR solver.
 * @param[in,out] data pointer to the data structure
 * @return always zero
 *
 * The @ref lapackqr_clear function clear a lapack-qr solver.
 *
 * @attention Internal use only.
 */
static int lapackqr_clear(void *data){
    struct lapackqr_solver *sol = (struct lapackqr_solver *) data;
    if ( sol->val != NULL)          mess_free(sol->val);
    if ( sol->val_cpx != NULL)      mess_free(sol->val_cpx);
    if ( sol->tau != NULL)          mess_free(sol->tau);
    if ( sol->tau_cpx != NULL)      mess_free(sol->tau_cpx);
    mess_free(sol);
    return 0;
}


/*-----------------------------------------------------------------------------
 *  OVERDETERMINED SYSTEM QR decomposition
 *-----------------------------------------------------------------------------*/

/**
 * @internal
 * @brief Solve \f$ Ax=(QR)x=b \f$.
 * @param[in] data   input data of the solver
 * @param[in] b      input right hand side
 * @param[in,out] x  input/output solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref lapackqr_solve_over function solves the overdetermined system
 * \f[ Ax=b \Leftrightarrow QRx=b.\f]
 *
 * @attention Internal use only.
 */
static int lapackqr_solve_over(void * data, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_int_t info = 0;
    mess_int_t one = 1;
    mess_int_t lwork = 0;
    double *work = NULL;
    mess_double_cpx_t *work_cpx = NULL;
    double workspace;
    mess_double_cpx_t workspace_cpx;


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    if(b->dim!=sol->rows){
        MSG_ERROR("b has wrong dimension "MESS_PRINTF_INT"!="MESS_PRINTF_INT"\n",b->dim,sol->rows);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/

    if (sol->cpx == 0) {
        if (  MESS_IS_REAL(b)) {
            mess_vector_copy(b,x);
            mess_vector_toreal(x);

            //workspace query and real work
            lwork = -1;
            F77_GLOBAL(dormqr,DORMQR)("L","T",&(sol->rows), &one, &(sol->cols), sol->val, &(sol->rows), sol->tau, x->values, &(x->dim), &workspace, &lwork, &info);
            lwork = nearbyint(workspace+1);

            mess_try_alloc(work, double *, sizeof(double)*lwork);

            // solve real system
            F77_GLOBAL(dormqr,DORMQR)("L","T",&(sol->rows), &one, &(sol->cols), sol->val, &(sol->rows), sol->tau, x->values, &(x->dim), work, &lwork, &info);
            F77_GLOBAL(dtrtrs,DTRTRS)("U","N","N", &(sol->cols), &one, sol->val, &(sol->rows), x->values, &(x->dim), &info);
            mess_free(work);
        } else {
            mess_vector xr,xi;
            MESS_INIT_VECTORS(&xr,&xi);
            mess_vector_alloc(xr, b->dim, MESS_REAL);
            mess_vector_alloc(xi, b->dim, MESS_REAL);
            mess_vector_realpart(b,xr);
            mess_vector_imagpart(b,xi);

            //workspace query and real work
            lwork = -1;
            F77_GLOBAL(dormqr,DORMQR)("L","T",&(sol->rows), &one, &(sol->cols), sol->val, &(sol->rows), sol->tau, xr->values, &(xr->dim), &workspace, &lwork, &info);
            lwork = nearbyint(workspace+1);

            mess_try_alloc(work, double *, sizeof(double)*lwork);

            // solve real system
            F77_GLOBAL(dormqr,DORMQR)("L","T",&(sol->rows), &one, &(sol->cols), sol->val, &(sol->rows), sol->tau, xr->values, &(xr->dim), work, &lwork, &info);
            F77_GLOBAL(dtrtrs,DTRTRS)("U","N","N", &(sol->cols), &one, sol->val, &(sol->rows), xr->values, &(xr->dim), &info);

            F77_GLOBAL(dormqr,DORMQR)("L","T",&(sol->rows), &one, &(sol->cols), sol->val, &(sol->rows), sol->tau, xi->values, &(xi->dim), work, &lwork, &info);
            F77_GLOBAL(dtrtrs,DTRTRS)("U","N","N", &(sol->cols), &one, sol->val, &(sol->rows), xi->values, &(xi->dim), &info);

            mess_free(work);

            mess_vector_complex_from_parts(xr,xi,x);
            mess_vector_clear(&xr);
            mess_vector_clear(&xi);
        }
    } else {
        mess_vector_copy(b,x);
        mess_vector_tocomplex(x);

        //workspace query and real work
        lwork = -1;
        F77_GLOBAL(zunmqr,ZUNMQR)("L","C",&(sol->rows), &one, &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, x->values_cpx, &(x->dim), &workspace_cpx, &lwork, &info);
        lwork = workspace_cpx;

        mess_try_alloc(work_cpx, mess_double_cpx_t  *, sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zunmqr,ZUNMQR)("L","C",&(sol->rows), &one, &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, x->values_cpx, &(x->dim), work_cpx, &lwork, &info);
        F77_GLOBAL(ztrtrs,ZTRTRS)("U","N","N", &(sol->cols), &one, sol->val_cpx, &(sol->rows), x->values_cpx, &(x->dim), &info);
        mess_free(work_cpx);
    }
    mess_vector_resize(x, sol->cols);
    if ( info < 0) {
        MSG_ERROR("Error calling DORMQR/DTRTRS ZUNMQR/ZTRTRS. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }
    return (int)(info);
}

/**
 * @internal
 * @brief Solve \f$ A^Tx=(QR)^Tx=b \f$.
 * @param[in] data   input data of the solver
 * @param[in] b      input right hand side
 * @param[in,out] x  input/output solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref lapackqr_solvet_over function solves the overdetermined system
 * \f[ A^Tx=b \Leftrightarrow (QR)^Tx=b.\f]
 *
 * @attention Internal use only.
 */
static int lapackqr_solvet_over(void * data, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_int_t info = 0;
    mess_int_t one = 1;
    mess_int_t lwork = 0;
    mess_int_t ret = 0;
    double *work = NULL;
    mess_double_cpx_t *work_cpx = NULL;
    double workspace;
    mess_double_cpx_t workspace_cpx;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    if(b->dim!=sol->cols){
        MSG_ERROR("b has wrong dimension "MESS_PRINTF_INT"!="MESS_PRINTF_INT"\n",b->dim,sol->cols);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/

    if (sol->cpx == 0) {
        if (MESS_IS_REAL(b)) {
            mess_vector_copy(b,x);
            mess_vector_toreal(x);

            // solve system R^Ty = b
            F77_GLOBAL(dtrtrs,DTRTRS)("U","T","N", &(sol->cols), &one, sol->val, &(sol->rows), x->values, &(x->dim), &info);

            // Number of rows > Number of cols for over case, Number of cols == b->dim, therefore we append some zeros at the end of x
            ret = mess_vector_resize(x, sol->rows);

            // workspace query and real work for operation x<-Qy
            lwork = -1;
            F77_GLOBAL(dormqr,DORMQR)("L","N",&(sol->rows), &one, &(sol->cols), sol->val, &(sol->rows), sol->tau, x->values, &(x->dim), &workspace, &lwork, &info);
            lwork = nearbyint(workspace+1);
            mess_try_alloc(work, double *, sizeof(double)*lwork);

            // perform multiplication
            F77_GLOBAL(dormqr,DORMQR)("L","N",&(sol->rows), &one, &(sol->cols), sol->val, &(sol->rows), sol->tau, x->values, &(x->dim), work, &lwork, &info);
            mess_free(work);

        }else{
            mess_vector xr,xi;
            MESS_INIT_VECTORS(&xr,&xi);
            mess_vector_alloc(xr,b->dim,MESS_REAL);
            mess_vector_alloc(xi,b->dim,MESS_REAL);
            mess_vector_realpart(b,xr);
            mess_vector_imagpart(b,xi);

            // solve R^T xr = Re(b) and R^T xc = Im(b)
            F77_GLOBAL(dtrtrs,DTRTRS)("U","T","N", &(sol->cols), &one, sol->val, &(sol->rows), xr->values, &(xr->dim), &info);
            F77_GLOBAL(dtrtrs,DTRTRS)("U","T","N", &(sol->cols), &one, sol->val, &(sol->rows), xi->values, &(xi->dim), &info);

            // Number of rows > Number of cols for over case, Number of cols == b->dim, therefore we append some zeros at the end of x
            ret = mess_vector_resize(xr, sol->rows);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);
            ret = mess_vector_resize(xi, sol->rows);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);

            // workspace query and real work for operation x<-Qy
            lwork = -1;
            F77_GLOBAL(dormqr,DORMQR)("L","N",&(sol->rows), &one, &(sol->cols), sol->val, &(sol->rows), sol->tau, xr->values, &(xr->dim), &workspace, &lwork, &info);
            lwork = nearbyint(workspace+1);
            mess_try_alloc(work, double *, sizeof(double)*lwork);

            // perform multiplication
            F77_GLOBAL(dormqr,DORMQR)("L","N",&(sol->rows), &one, &(sol->cols), sol->val, &(sol->rows), sol->tau, xr->values, &(xr->dim), work, &lwork, &info);
            F77_GLOBAL(dormqr,DORMQR)("L","N",&(sol->rows), &one, &(sol->cols), sol->val, &(sol->rows), sol->tau, xi->values, &(xi->dim), work, &lwork, &info);
            mess_free(work);

            mess_vector_complex_from_parts(xr,xi,x);
            mess_vector_clear(&xr);
            mess_vector_clear(&xi);

        }
    } else {
        mess_vector_copy(b,x);
        mess_vector_tocomplex(x);

        // solve R^T y = b
        F77_GLOBAL(ztrtrs,ZTRTRS)("U","T","N", &(sol->cols), &one, sol->val_cpx, &(sol->rows), x->values_cpx, &(x->dim), &info);

        // Number of rows > Number of cols for over case, Number of cols == b->dim, therefore we append some zeros at the end of x
        ret = mess_vector_resize(x, sol->rows);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);

        //workspace query and real work for operation x<-Qy
        ret = mess_vector_conj(x);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_conj);
        lwork = -1;
        F77_GLOBAL(zunmqr,ZUNMQR)("L", "N", &(sol->rows), &one, &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, x->values_cpx, &(x->dim), &workspace_cpx, &lwork, &info);
        lwork = workspace_cpx;

        //perform multiplication
        mess_try_alloc(work_cpx, mess_double_cpx_t  *, sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zunmqr,ZUNMQR)("L","N",&(sol->rows), &one, &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, x->values_cpx, &(x->dim), work_cpx, &lwork, &info);
        mess_free(work_cpx);
        ret = mess_vector_conj(x);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_conj);

    }

    if ( info < 0) {
        MSG_ERROR("Error calling DORMQR/DTRTRS ZUNMQR/ZTRTRS. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }
    return (int)(info);
}

/**
 * @internal
 * @brief Solve \f$ A^Hx=(QR)^Hx=b \f$.
 * @param[in] data   input data of the solver
 * @param[in] b      input right hand side
 * @param[in,out] x  input/output solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref lapackqr_solvet_over function solves the overdetermined system
 * \f[ A^Tx=b \Leftrightarrow (QR)^Hx=b.\f]
 *
 * @attention Internal use only.
 */
static int lapackqr_solveh_over(void * data, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_int_t info = 0;
    mess_int_t one = 1;
    mess_int_t lwork = 0;
    mess_int_t ret = 0;
    mess_double_cpx_t *work_cpx = NULL;
    mess_double_cpx_t workspace_cpx;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    if(b->dim!=sol->cols){
        MSG_ERROR("b has wrong dimension "MESS_PRINTF_INT"!="MESS_PRINTF_INT"\n",b->dim,sol->cols);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/

    if (sol->cpx == 0) {
            ret = lapackqr_solvet_over(data,b,x);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),lapackqr_solvet_over);
    } else {
        mess_vector_copy(b,x);
        mess_vector_tocomplex(x);

        // solve R^H y = b
        F77_GLOBAL(ztrtrs,ZTRTRS)("U","C","N", &(sol->cols), &one, sol->val_cpx, &(sol->rows), x->values_cpx, &(x->dim), &info);

        // Number of rows > Number of cols for over case, Number of cols == b->dim, therefore we append some zeros at the end of x
        ret = mess_vector_resize(x, sol->rows);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);

        //workspace query and real work for operation x<-Qy
        lwork = -1;
        F77_GLOBAL(zunmqr,ZUNMQR)("L", "N", &(sol->rows), &one, &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, x->values_cpx, &(x->dim), &workspace_cpx, &lwork, &info);
        lwork = workspace_cpx;

        //perform multiplication
        mess_try_alloc(work_cpx, mess_double_cpx_t  *, sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zunmqr,ZUNMQR)("L","N",&(sol->rows), &one, &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, x->values_cpx, &(x->dim), work_cpx, &lwork, &info);
        mess_free(work_cpx);
    }

    if ( info < 0) {
        MSG_ERROR("Error calling DORMQR/DTRTRS ZUNMQR/ZTRTRS. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }
    return (int)(info);
}

/**
 * @internal
 * @brief Solve \f$ AX=(QR)X=B \f$ (matrix version).
 * @param[in] data   input data of the solver
 * @param[in] b      input right hand side
 * @param[in,out] x     solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref lapackqr_solvem_over function solves the overdetermined system
 *\f[ Ax=b \Leftrightarrow QRx=b,\f]
 * where \f$ b \f$ and \f$ x \f$ are matrices.
 * @attention Internal use only.
 */
static int lapackqr_solvem_over(void * data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_int_t info = 0;
    mess_int_t lwork = 0;
    double *work = NULL;
    mess_double_cpx_t *work_cpx = NULL;
    double workspace;
    mess_double_cpx_t workspace_cpx;
    mess_int_t ret = 0 ;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    if(b->rows!=sol->rows){
        MSG_ERROR("b has wrong number of rows "MESS_PRINTF_INT"!="MESS_PRINTF_INT"\n",b->rows,sol->rows);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_convert(b,x,MESS_DENSE);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert);

    if (sol->cpx == 0) {
        if ( MESS_IS_REAL(b)) {

            //workspace query and real work
            lwork = -1;
            F77_GLOBAL(dormqr,DORMQR)("L","T",&(x->rows), &(x->cols), &(sol->cols), sol->val, &(sol->rows), sol->tau, x->values, &(x->ld), &workspace, &lwork, &info);
            lwork = nearbyint(workspace);

            mess_try_alloc(work, double *, sizeof(double)*lwork);
            F77_GLOBAL(dormqr,DORMQR)("L","T",&(x->rows), &(x->cols), &(sol->cols), sol->val, &(sol->rows), sol->tau, x->values, &(x->ld), work, &lwork, &info);
            F77_GLOBAL(dtrtrs,DTRTRS)("U","N","N", &(sol->cols), &(x->cols), sol->val, &(sol->rows), x->values, &(x->ld), &info);
            mess_free(work);
        } else {
            mess_matrix br, bi;
            ret = mess_matrix_init(&br);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
            ret = mess_matrix_init(&bi);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
            ret = mess_matrix_realpart(x,br);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realpart);
            ret = mess_matrix_imagpart(x,bi);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_imagpart);

            //workspace query and real work
            lwork = -1;
            F77_GLOBAL(dormqr,DORMQR)("L","T",&(x->rows), &(x->cols), &(sol->cols), sol->val, &(sol->rows), sol->tau, br->values, &(br->ld),&workspace, &lwork, &info);
            lwork = nearbyint(workspace+1);

            mess_try_alloc(work, double *, sizeof(double)*lwork);
            F77_GLOBAL(dormqr,DORMQR)("L","T",&(x->rows), &(x->cols), &(sol->cols), sol->val, &(sol->rows), sol->tau, br->values, &(br->ld), work, &lwork, &info);
            F77_GLOBAL(dtrtrs,DTRTRS)("U","N","N", &(sol->cols), &(x->cols), sol->val, &(sol->rows), br->values, &(br->ld), &info);
            F77_GLOBAL(dormqr,DORMQR)("L","T",&(x->rows), &(x->cols), &(sol->cols), sol->val, &(sol->rows), sol->tau, bi->values, &(bi->ld), work, &lwork, &info);
            F77_GLOBAL(dtrtrs,DTRTRS)("U","N","N", &(sol->cols), &(x->cols), sol->val, &(sol->rows), bi->values, &(br->ld), &info);
            ret = mess_matrix_complex_from_parts(br,bi, x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_complex_from_parts);
            mess_free(work);
            mess_matrix_clear(&br);
            mess_matrix_clear(&bi);
            ret = mess_matrix_resize(x, sol->cols, x->cols);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_resize);
        }
    } else {
        mess_matrix_tocomplex(x);

        //workspace query and real work
        lwork=-1;
        F77_GLOBAL(zunmqr,ZUNMQR)("L","C",&(x->rows), &(x->cols), &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, x->values_cpx, &(x->ld), &workspace_cpx, &lwork, &info);
        lwork = nearbyint(creal(workspace_cpx)+1);

        mess_try_alloc(work_cpx, mess_double_cpx_t  *, sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zunmqr,ZUNMQR)("L","C",&(x->rows), &(x->cols), &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, x->values_cpx, &(x->ld), work_cpx, &lwork, &info);
        F77_GLOBAL(ztrtrs,ZTRTRS)("U","N","N", &(sol->cols), &(x->cols), sol->val_cpx, &(sol->rows), x->values_cpx, &(x->ld), &info);
        mess_free(work_cpx);
    }

    mess_matrix_resize(x, sol->cols, x->cols);
    if ( info < 0) {
        MSG_ERROR("Error calling DORMQR/DTRTRS or ZUNMQR/ZTRTRS. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }
    return (int)(info);
}

/**
 * @internal
 * @brief Solve \f$ A^Tx=(QR)^Tx=b (matrix version)\f$.
 * @param[in] data   input data of the solver
 * @param[in] b      input right hand side
 * @param[in,out] x  input/output solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref lapackqr_solvemt_over function solves the overdetermined system
 * \f[ A^Tx=b \Leftrightarrow (QR)^Tx=b.\f]
 *
 * @attention Internal use only.
 */
static int lapackqr_solvemt_over(void * data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_int_t info = 0;
    mess_int_t lwork = 0;
    mess_int_t ret = 0;
    double *work = NULL;
    mess_double_cpx_t *work_cpx = NULL;
    double workspace;
    mess_double_cpx_t workspace_cpx;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    if(b->rows!=sol->cols){
        MSG_ERROR("b has wrong dimension "MESS_PRINTF_INT"!="MESS_PRINTF_INT"\n",b->rows,sol->cols);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    if (sol->cpx == 0) {
        if (MESS_IS_REAL(b)) {
            ret = mess_matrix_copy(b,x);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            ret = mess_matrix_toreal(x);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_toreal);

            // solve system R^Ty = b
            F77_GLOBAL(dtrtrs,DTRTRS)("U","T","N", &(sol->cols), &(x->cols), sol->val, &(sol->rows), x->values, &(x->ld), &info);

            // Number of rows > Number of cols for over case, Number of cols == b->dim, therefore we append some zeros at the end of x
            ret = mess_matrix_resize(x, sol->rows, x->cols);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_resize);

            // workspace query and real work for operation x<-Qy
            lwork = -1;
            F77_GLOBAL(dormqr,DORMQR)("L","N",&(sol->rows), &(x->cols), &(sol->cols), sol->val, &(sol->rows), sol->tau, x->values, &(x->ld), &workspace, &lwork, &info);
            lwork = nearbyint(workspace+1);
            mess_try_alloc(work, double *, sizeof(double)*lwork);

            // perform multiplication
            F77_GLOBAL(dormqr,DORMQR)("L","N",&(sol->rows), &(x->cols), &(sol->cols), sol->val, &(sol->rows), sol->tau, x->values, &(x->ld), work, &lwork, &info);
            mess_free(work);

        }else{
            mess_matrix xr,xi;
            MESS_INIT_MATRICES(&xr,&xi);
            ret = mess_matrix_realpart(b,xr);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_realpart);
            ret = mess_matrix_imagpart(b,xi);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_imagpart);


            // solve R^T xr = Re(b) and R^T xc = Im(b)
            F77_GLOBAL(dtrtrs,DTRTRS)("U","T","N", &(sol->cols), &(xr->cols), sol->val, &(sol->rows), xr->values, &(xr->ld), &info);
            F77_GLOBAL(dtrtrs,DTRTRS)("U","T","N", &(sol->cols), &(xr->cols), sol->val, &(sol->rows), xi->values, &(xi->ld), &info);

            // Number of rows > Number of cols for over case, Number of cols == b->dim, therefore we append some zeros at the end of x
            ret = mess_matrix_resize(xr, sol->rows, xr->cols);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_resize);
            ret = mess_matrix_resize(xi, sol->rows, xi->cols);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_resize);

            // workspace query and real work for operation x<-Qy
            lwork = -1;
            F77_GLOBAL(dormqr,DORMQR)("L","N",&(sol->rows), &(xr->cols), &(sol->cols), sol->val, &(sol->rows), sol->tau, xr->values, &(xr->ld), &workspace, &lwork, &info);
            lwork = nearbyint(workspace+1);
            mess_try_alloc(work, double *, sizeof(double)*lwork);

            // perform multiplication
            F77_GLOBAL(dormqr,DORMQR)("L","N",&(sol->rows), &(xr->cols), &(sol->cols), sol->val, &(sol->rows), sol->tau, xr->values, &(xr->ld), work, &lwork, &info);
            F77_GLOBAL(dormqr,DORMQR)("L","N",&(sol->rows), &(xi->cols), &(sol->cols), sol->val, &(sol->rows), sol->tau, xi->values, &(xi->ld), work, &lwork, &info);
            mess_free(work);

            ret = mess_matrix_complex_from_parts(xr,xi,x);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_complex_from_parts);
            MESS_CLEAR_MATRICES(&xr,&xi);

        }
    } else {
        ret = mess_matrix_copy(b,x);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_tocomplex(x);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);

        // solve R^T y = b
        F77_GLOBAL(ztrtrs,ZTRTRS)("U","T","N", &(sol->cols), &(x->cols), sol->val_cpx, &(sol->rows), x->values_cpx, &(x->ld), &info);

        // Number of rows > Number of cols for over case, Number of cols == b->dim, therefore we append some zeros at the end of x
        ret = mess_matrix_resize(x, sol->rows, x->cols);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_resize);

        //workspace query and real work for operation x<-Qy
        ret = mess_matrix_conj(x);                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_conj);
        lwork = -1;
        F77_GLOBAL(zunmqr,ZUNMQR)("L", "N", &(sol->rows), &(x->cols), &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, x->values_cpx, &(x->ld), &workspace_cpx, &lwork, &info);
        lwork = workspace_cpx;

        //perform multiplication
        mess_try_alloc(work_cpx, mess_double_cpx_t  *, sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zunmqr,ZUNMQR)("L","N",&(sol->rows), &(x->cols), &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, x->values_cpx, &(x->ld), work_cpx, &lwork, &info);
        mess_free(work_cpx);
        ret = mess_matrix_conj(x);                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_conj);

    }

    if ( info < 0) {
        MSG_ERROR("Error calling DORMQR/DTRTRS ZUNMQR/ZTRTRS. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }
    return (int)(info);
}

/**
 * @internal
 * @brief Solve \f$ A^Hx=(QR)^Hx=b (matrix version)\f$.
 * @param[in] data   input data of the solver
 * @param[in] b      input right hand side
 * @param[in,out] x  input/output solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref lapackqr_solvemh_over function solves the overdetermined system
 * \f[ A^Tx=b \Leftrightarrow (QR)^Hx=b.\f]
 *
 * @attention Internal use only.
 */
static int lapackqr_solvemh_over(void * data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_int_t info = 0;
    mess_int_t lwork = 0;
    mess_int_t ret = 0;
    mess_double_cpx_t *work_cpx = NULL;
    mess_double_cpx_t workspace_cpx;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    if(b->rows!=sol->cols){
        MSG_ERROR("b has wrong dimension "MESS_PRINTF_INT"!="MESS_PRINTF_INT"\n",b->rows,sol->cols);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/

    if (sol->cpx == 0) {
            ret = lapackqr_solvemt_over(data,b,x);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),lapackqr_solvemt_over);
    } else {
        ret = mess_matrix_copy(b,x);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_tocomplex(x);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);

        // solve R^H y = b
        F77_GLOBAL(ztrtrs,ZTRTRS)("U","C","N", &(sol->cols), &(x->cols), sol->val_cpx, &(sol->rows), x->values_cpx, &(x->ld), &info);

        // Number of rows > Number of cols for over case, Number of cols == b->dim, therefore we append some zeros at the end of x
        ret = mess_matrix_resize(x, sol->rows, x->cols);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_resize);

        //workspace query and real work for operation x<-Qy
        lwork = -1;
        F77_GLOBAL(zunmqr,ZUNMQR)("L", "N", &(sol->rows), &(x->cols), &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, x->values_cpx, &(x->ld), &workspace_cpx, &lwork, &info);
        lwork = workspace_cpx;

        //perform multiplication
        mess_try_alloc(work_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zunmqr,ZUNMQR)("L","N",&(sol->rows), &(x->cols), &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, x->values_cpx, &(x->ld), work_cpx, &lwork, &info);
        mess_free(work_cpx);
    }

    if ( info < 0) {
        MSG_ERROR("Error calling DORMQR/DTRTRS ZUNMQR/ZTRTRS. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }
    return (int)(info);
}

/**
 * @brief Get the L (the first factor) of a solver.
 * @param[in] data  input solver data
 * @param[in,out] L     output L factor
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapackqr_getL_over function gets the \f$ Q \f$ factor of \f$ A =QR \f$. Output is in @ref MESS_DENSE
 * format.
 *
 */
static int lapackqr_getL_over(void *data, mess_matrix L ){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_check_nullpointer(sol);
    mess_check_nullpointer(L);
    int ret=0;
    mess_int_t info=0;
    double * work=NULL;
    mess_double_cpx_t *work_cpx = NULL;
    mess_int_t lwork;
    double workspace;
    mess_double_cpx_t workspace_cpx;

    if(sol->rows<sol->cols){
        MSG_ERROR("For underdetermined systems this case is not implemented\n");
        return MESS_ERROR_MISSING;
    }

    ret = mess_matrix_reset(L);                                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_reset);
    ret = mess_matrix_alloc(L,sol->rows,sol->cols,sol->rows*sol->cols, MESS_DENSE,sol->cpx? MESS_COMPLEX:MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);

    if ( MESS_IS_REAL(L)) {
        F77_GLOBAL(dlacpy,DLACPY)("All", &(sol->rows), &(sol->cols), sol->val, &(sol->rows), L->values, &(L->ld));

        //work space query and real work
        lwork = -1;
        F77_GLOBAL(dorgqr,DORGQR)(&(sol->rows),&(sol->cols),&(sol->cols),L->values,&(L->ld),sol->tau, &workspace,&lwork,&info);
        lwork = nearbyint(workspace+1);

        mess_try_alloc(work, double*, sizeof(double)* lwork);
        F77_GLOBAL(dorgqr,DORGQR)(&(sol->rows),&(sol->cols),&(sol->cols),L->values,&(L->ld),sol->tau, work,&lwork,&info);
        mess_free(work);
    } else {
        F77_GLOBAL(zlacpy,ZLACPY)("All", &(sol->rows), &(sol->cols), sol->val_cpx, &(sol->rows), L->values_cpx, &(L->ld));

        //work space query and real work
        lwork = -1;
        F77_GLOBAL(zungqr,ZUNGQR)(&(sol->rows),&(sol->cols),&(sol->cols),L->values_cpx,&(L->ld),sol->tau_cpx,&workspace_cpx,&lwork,&info);
        lwork = nearbyint(creal(workspace_cpx)+1);

        mess_try_alloc(work_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)* lwork);
        F77_GLOBAL(zungqr,ZUNGQR)(&(sol->rows),&(sol->cols),&(sol->cols),L->values_cpx,&(L->ld),sol->tau_cpx,work_cpx,&lwork,&info);
        mess_free(work_cpx);
    }

    return ret;
}

/**
 * @brief Get the U (the second factor) of a solver.
 * @param[in] data  input solver data
 * @param[in,out] U     output U factor
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapackqr_getU_over function gets the \f$ R \f$ factor of \f$ A =QR \f$. Output is in @ref MESS_DENSE
 * format.
 *
 */
static int lapackqr_getU_over(void *data, mess_matrix U ){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_check_nullpointer(sol);
    mess_check_nullpointer(U);
    mess_int_t cols;
    cols = sol->cols;
    int ret=0;

    if(sol->rows<sol->cols){
        MSG_ERROR("For underdetermined systems this case is not implemented\n");
        return MESS_ERROR_MISSING;
    }

    ret = mess_matrix_reset(U);                                                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_reset);
    ret = mess_matrix_alloc(U,cols,cols,cols*cols,MESS_DENSE,sol->cpx? MESS_COMPLEX:MESS_REAL); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);

    if(MESS_IS_COMPLEX(U)){
        F77_GLOBAL(zlacpy,ZLACPY)("U",&(sol->cols),&(sol->cols),sol->val_cpx,&(sol->rows),U->values_cpx,&(U->ld));
    }else{
        F77_GLOBAL(dlacpy,DLACPY)("U",&(sol->cols),&(sol->cols),sol->val,&(sol->rows),U->values,&(U->ld));
    }
    return ret;
}

/**
 * @brief Construct a solver based on QR-decomposition from @lapack for overdetermined systems.
 * @param[in] matrix        input system matrix
 * @param[out] solver       output
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_create_lapack_qr_over function generates a solver based on the QR-decomposition
 * from @lapack for overdetermined systems, i.e. the number of rows is greater equals than the number of
 * cols.
 **/
static int mess_direct_create_lapack_qr_over(mess_matrix matrix, mess_direct solver){
    MSG_FNAME(__func__);
    mess_int_t info = 0;
    mess_int_t lwork = 0;
    mess_int_t ret = 0;
    double *work=NULL, workspace;
    mess_double_cpx_t *work_cpx = NULL, workspace_cpx;
    struct lapackqr_solver *sol;
    mess_int_t conv = 0;
    mess_matrix mwork;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_real_or_complex(matrix);
    if(matrix->rows<matrix->cols){
        MSG_ERROR("Input matrix has less rows than cols. Please provide an overdetermined system.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  create solver
     *-----------------------------------------------------------------------------*/
    MESS_MATRIX_CHECKFORMAT(matrix,mwork,conv,MESS_DENSE);
    mess_try_alloc(sol, struct lapackqr_solver*, sizeof(struct lapackqr_solver));
    sol->rows = matrix->rows;
    sol->cols = matrix->cols;

    if ( MESS_IS_REAL(matrix)) {
        MSG_INFO("real\n");
        sol->val_cpx = NULL;
        sol->tau_cpx = NULL;
        mess_try_alloc(sol->tau, double *, sizeof(double)*(MESS_MIN(sol->rows, sol->cols)));
        mess_try_alloc(sol->val, double *, sizeof(double)*sol->rows*sol->cols);
        F77_GLOBAL(dlacpy,DLACPY)("All", &matrix->rows, &matrix->cols, mwork->values, &mwork->ld, sol->val, &sol->rows);

        /*-----------------------------------------------------------------------------
         *   work space query and decomposition
         *-----------------------------------------------------------------------------*/
        lwork = -1;
        F77_GLOBAL(dgeqrf,DGEQRF)(&(sol->rows),&(sol->cols), sol->val, &(sol->rows), sol->tau, &workspace, &lwork, &info);
        lwork = nearbyint(workspace+1);
        MSG_INFO("Workspace:" MESS_PRINTF_INT"\n",lwork);

        mess_try_alloc(work, double*, sizeof(double)* lwork);
        F77_GLOBAL(dgeqrf,DGEQRF)(&(sol->rows),&(sol->cols), sol->val, &(sol->rows), sol->tau, work, &lwork, &info);
        mess_free(work);
        sol->cpx =0;
    } else {
        MSG_INFO("complex\n");
        sol->val = NULL;
        sol->tau = NULL;
        mess_try_alloc(sol->tau_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*(MESS_MIN(sol->rows, sol->cols)));
        mess_try_alloc(sol->val_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*sol->rows*sol->cols);
        F77_GLOBAL(zlacpy,ZLACPY)("All", &matrix->rows, &matrix->cols, mwork->values_cpx, &mwork->ld, sol->val_cpx, &sol->rows);

        /*-----------------------------------------------------------------------------
         *   work space query and decomposition
         *-----------------------------------------------------------------------------*/
        lwork = -1;
        F77_GLOBAL(zgeqrf,ZGEQRF)(&(sol->rows),&(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, &workspace_cpx, &lwork, &info);
        lwork = nearbyint(creal(workspace_cpx)+1);

        mess_try_alloc(work_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)* lwork);
        F77_GLOBAL(zgeqrf,ZGEQRF)(&(sol->rows),&(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, work_cpx, &lwork, &info);
        mess_free(work_cpx);
        sol->cpx =1;
    }

    solver->data_type = matrix->data_type;
    solver->data    = (void *) sol;
    solver->clear   = lapackqr_clear;
    solver->solve   = lapackqr_solve_over;
    solver->solvem  = lapackqr_solvem_over;
    solver->solvet  = lapackqr_solvet_over;
    solver->solvemt = lapackqr_solvemt_over;
    solver->solveh  = lapackqr_solveh_over;
    solver->solvemh = lapackqr_solvemh_over;

    solver->getL = lapackqr_getL_over;
    solver->getU = lapackqr_getU_over;
    solver->getpermp = NULL;
    solver->getpermq = NULL;
    SET_SOLVERNAME(solver->name, __func__);

    if ( conv == 0 ) {
        mess_matrix_clear(&mwork);
    }
    if ( info != 0 ) {
        MSG_ERROR("An error occured in DGEQRF/ZGEQRF: " MESS_PRINTF_INT "\n", info);
        ret = (int)(info);
    }

    return ret;
}


/*-----------------------------------------------------------------------------
 *  UNDERDETERMINED SYSTEMS LQ decomposition
 *-----------------------------------------------------------------------------*/

/**
 * @internal
 * @brief Solve \f$ Ax=(LQ)x=b \f$.
 * @param[in] data   input data of the solver
 * @param[in] b      input right hand side
 * @param[in,out] x  input/output solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref lapackqr_solve_under function solves the underdetermined system
 * \f[ Ax=b \Leftrightarrow LQx=b.\f]
 *
 * @attention Internal use only.
 */
static int lapackqr_solve_under(void * data, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_int_t info = 0;
    mess_int_t one = 1;
    mess_int_t lwork = 0;
    double *work = NULL;
    mess_double_cpx_t *work_cpx = NULL;
    double workspace;
    mess_double_cpx_t workspace_cpx;


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    if(b->dim!=sol->rows){
        MSG_ERROR("b has wrong dimension "MESS_PRINTF_INT"!="MESS_PRINTF_INT"\n",b->dim,sol->rows);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/

    if(sol->cpx){
        mess_vector_copy(b,x);
        mess_vector_tocomplex(x);
        mess_vector_resize(x,sol->cols);

        //solve with L
        F77_GLOBAL(ztrtrs,ZTRTRS)("L","N","N",&(sol->rows),&one,sol->val_cpx,&(sol->rows),x->values_cpx,&(x->dim),&info);

        //workspace query and real work
        lwork = -1;
        F77_GLOBAL(zunmlq,ZUNMLQ)("L","C",&(sol->cols),&one,&(sol->rows),sol->val_cpx,&(sol->rows),sol->tau_cpx,x->values_cpx,&(x->dim),&(workspace_cpx),&lwork,&(info));
        lwork =  nearbyint(creal(workspace_cpx)+1);
        mess_try_alloc(work_cpx,mess_double_cpx_t*, sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zunmlq,ZUNMLQ)("L","C",&(sol->cols),&one,&(sol->rows),sol->val_cpx,&(sol->rows),sol->tau_cpx,x->values_cpx,&(x->dim),work_cpx,&lwork,&(info));
        mess_free(work_cpx);
    }else{
        if(MESS_IS_REAL(b)){
            mess_vector_copy(b,x);
            mess_vector_resize(x,sol->cols);

            //solve with L
            F77_GLOBAL(dtrtrs,DTRTRS)("L","N","N",&(sol->rows), &(one), sol->val, &(sol->rows), x->values,&(x->dim),&info);

            //workspace query and real work
            lwork = -1;
            F77_GLOBAL(dormlq,DORMLQ)("L","T",&(sol->cols),&one,&(sol->rows),sol->val,&(sol->rows),sol->tau, x->values,&(x->dim),&(workspace),&lwork,&info);

            lwork = nearbyint(workspace+1);
            mess_try_alloc(work,double*,sizeof(double)*lwork);
            F77_GLOBAL(dormlq,DORMLQ)("L","T",&(sol->cols),&one,&(sol->rows),sol->val,&(sol->rows),sol->tau, x->values,&(x->dim),work,&lwork,&info);
            mess_free(work);
        }else{
            mess_vector xr,xc;
            MESS_INIT_VECTORS(&xr,&xc);
            mess_vector_alloc(xr,b->dim,MESS_REAL);
            mess_vector_alloc(xc,b->dim,MESS_REAL);
            mess_vector_realpart(b,xr); mess_vector_resize(xr,sol->cols);
            mess_vector_imagpart(b,xc); mess_vector_resize(xc,sol->cols);

            //solve with L
            F77_GLOBAL(dtrtrs,DTRTRS)("L","N","N",&sol->rows, &one, sol->val, &(sol->rows), xr->values,&(xr->dim),&info);
            F77_GLOBAL(dtrtrs,DTRTRS)("L","N","N",&sol->rows, &one, sol->val, &(sol->rows), xc->values,&(xc->dim),&info);

            //workspace query and real work
            lwork = -1;
            F77_GLOBAL(dormlq,DORMLQ)("L","T",&(sol->cols),&one,&(sol->rows),sol->val,&(sol->rows),sol->tau, xr->values,&(xr->dim),&(workspace),&lwork,&info);
            lwork = nearbyint(workspace+1);
            mess_try_alloc(work,double*,sizeof(double)*lwork);
            F77_GLOBAL(dormlq,DORMLQ)("L","T",&(sol->cols),&one,&(sol->rows),sol->val,&(sol->rows),sol->tau, xr->values,&(xr->dim),work,&lwork,&info);
            F77_GLOBAL(dormlq,DORMLQ)("L","T",&(sol->cols),&one,&(sol->rows),sol->val,&(sol->rows),sol->tau, xc->values,&(xc->dim),work,&lwork,&info);
            mess_free(work);
            mess_vector_complex_from_parts(xr,xc,x);
            mess_vector_clear(&xr);
            mess_vector_clear(&xc);
        }
    }

    if ( info < 0) {
        MSG_ERROR("Error calling DORMLQ/DTRTRS or ZUNMLQ/ZTRTRS. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }
    return (int)(info);
}


/**
 * @internal
 * @brief Solve \f$ Ax=(LQ)^Tx=b \f$.
 * @param[in] data   input data of the solver
 * @param[in] b      input right hand side
 * @param[in,out] x  input/output solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref lapackqr_solve_under function solves the underdetermined system
 * \f[ Ax=b \Leftrightarrow (LQ)^Tx=b.\f]
 *
 * @attention Internal use only.
 */
static int lapackqr_solvet_under(void * data, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_int_t info = 0;
    mess_int_t one = 1;
    mess_int_t lwork = 0;
    mess_int_t ret =0;
    double *work = NULL;
    mess_double_cpx_t *work_cpx = NULL;
    double workspace;
    mess_double_cpx_t workspace_cpx;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    if(b->dim!=sol->cols){
        MSG_ERROR("b has wrong dimension "MESS_PRINTF_INT"!="MESS_PRINTF_INT"\n",b->dim,sol->cols);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/

    if(sol->cpx){
        mess_vector_copy(b,x);
        mess_vector_tocomplex(x);

        //workspace query and real work
        ret = mess_vector_conj(x);                                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_conj);
        lwork = -1;
        F77_GLOBAL(zunmlq,ZUNMLQ)("L","N",&(sol->cols),&one,&(sol->rows),sol->val_cpx,&(sol->rows),sol->tau_cpx,x->values_cpx,&(x->dim),&(workspace_cpx),&lwork,&(info));
        lwork =  nearbyint(creal(workspace_cpx)+1);
        mess_try_alloc(work_cpx,mess_double_cpx_t*, sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zunmlq,ZUNMLQ)("L","N",&(sol->cols),&one,&(sol->rows),sol->val_cpx,&(sol->rows),sol->tau_cpx,x->values_cpx,&(x->dim),work_cpx,&lwork,&(info));
        mess_free(work_cpx);
        ret = mess_vector_conj(x);                                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_conj);

        //solve with L
        F77_GLOBAL(ztrtrs,ZTRTRS)("L","T","N",&(sol->rows),&one,sol->val_cpx,&(sol->rows),x->values_cpx,&(x->dim),&info);
        ret = mess_vector_resize(x,sol->rows);

   }else{
        if(MESS_IS_REAL(b)){
            mess_vector_copy(b,x);
            //workspace query and real work
            lwork = -1;
            F77_GLOBAL(dormlq,DORMLQ)("L","N",&(sol->cols),&one,&(sol->rows),sol->val,&(sol->rows),sol->tau, x->values,&(x->dim),&(workspace),&lwork,&info);
            lwork = nearbyint(workspace+1);
            mess_try_alloc(work,double*,sizeof(double)*lwork);
            F77_GLOBAL(dormlq,DORMLQ)("L","N",&(sol->cols),&one,&(sol->rows),sol->val,&(sol->rows),sol->tau, x->values,&(x->dim),work,&lwork,&info);
            mess_free(work);

            //solve with L
            ret = mess_vector_resize(x, sol->rows);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);
            F77_GLOBAL(dtrtrs,DTRTRS)("L","T","N",&(sol->rows),&(one),sol->val,&(sol->rows),x->values,&(x->dim),&info);

        }else{
            mess_vector xr,xc;
            MESS_INIT_VECTORS(&xr,&xc);
            mess_vector_alloc(xr,b->dim,MESS_REAL);
            mess_vector_alloc(xc,b->dim,MESS_REAL);
            mess_vector_realpart(b,xr);
            mess_vector_imagpart(b,xc);

            //workspace query and real work
            lwork = -1;
            F77_GLOBAL(dormlq,DORMLQ)("L","N",&(sol->cols),&one,&(sol->rows),sol->val,&(sol->rows),sol->tau, xr->values,&(xr->dim),&(workspace),&lwork,&info);
            lwork = nearbyint(workspace+1);
            mess_try_alloc(work,double*,sizeof(double)*lwork);
            F77_GLOBAL(dormlq,DORMLQ)("L","N",&(sol->cols),&one,&(sol->rows),sol->val,&(sol->rows),sol->tau, xr->values,&(xr->dim),work,&lwork,&info);
            F77_GLOBAL(dormlq,DORMLQ)("L","N",&(sol->cols),&one,&(sol->rows),sol->val,&(sol->rows),sol->tau, xc->values,&(xc->dim),work,&lwork,&info);
            mess_free(work);

            //solve with L
            F77_GLOBAL(dtrtrs,DTRTRS)("L","T","N",&(sol->rows),&(one),sol->val,&(sol->rows),xr->values,&(xr->dim),&info);
            F77_GLOBAL(dtrtrs,DTRTRS)("L","T","N",&(sol->rows),&(one),sol->val,&(sol->rows),xc->values,&(xc->dim),&info);

            mess_vector_complex_from_parts(xr,xc,x);
            mess_vector_clear(&xr);
            mess_vector_clear(&xc);
            mess_vector_resize(x,sol->rows);

        }
    }

    if ( info < 0) {
        MSG_ERROR("Error calling DORMLQ/DTRTRS or ZUNMLQ/ZTRTRS. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }
    return (int)(info);
}


/**
 * @internal
 * @brief Solve \f$ Ax=(LQ)^Hx=b \f$.
 * @param[in] data   input data of the solver
 * @param[in] b      input right hand side
 * @param[in,out] x  input/output solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref lapackqr_solve_under function solves the underdetermined system
 * \f[ Ax=b \Leftrightarrow (LQ)^Hx=b.\f]
 *
 * @attention Internal use only.
 */
static int lapackqr_solveh_under(void * data, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_int_t info = 0;
    mess_int_t one = 1;
    mess_int_t lwork = 0;
    mess_int_t ret =0;
    mess_double_cpx_t *work_cpx = NULL;
    mess_double_cpx_t workspace_cpx;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    if(b->dim!=sol->cols){
        MSG_ERROR("b has wrong dimension "MESS_PRINTF_INT"!="MESS_PRINTF_INT"\n",b->dim,sol->cols);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/

    if(sol->cpx){
        mess_vector_copy(b,x);
        mess_vector_tocomplex(x);
        printf("here\n") ;

        //workspace query and real work
        lwork = -1;
        F77_GLOBAL(zunmlq,ZUNMLQ)("L","N",&(sol->cols),&one,&(sol->rows),sol->val_cpx,&(sol->rows),sol->tau_cpx,x->values_cpx,&(x->dim),&(workspace_cpx),&lwork,&(info));
        lwork =  nearbyint(creal(workspace_cpx)+1);
        mess_try_alloc(work_cpx,mess_double_cpx_t*, sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zunmlq,ZUNMLQ)("L","N",&(sol->cols),&one,&(sol->rows),sol->val_cpx,&(sol->rows),sol->tau_cpx,x->values_cpx,&(x->dim),work_cpx,&lwork,&(info));
        mess_free(work_cpx);

        //solve with L
        F77_GLOBAL(ztrtrs,ZTRTRS)("L","C","N",&(sol->rows),&one,sol->val_cpx,&(sol->rows),x->values_cpx,&(x->dim),&info);
        ret = mess_vector_resize(x,sol->rows);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);

   }else{
       return lapackqr_solvet_under(data,b,x);
    }

    if ( info < 0) {
        MSG_ERROR("Error calling DORMLQ/DTRTRS or ZUNMLQ/ZTRTRS. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }
    return (int)(info);
}


/**
 * @internal
 * @brief Solve \f$ Ax=(LQ)x=b \f$ (matrix version).
 * @param[in] data   input data of the solver
 * @param[in] b      input right hand side
 * @param[in,out] x  input/output solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref lapackqr_solvem_under function solves the underdetermined system
 * \f[ Ax=b \Leftrightarrow LQx=b.\f]
 *
 * @attention Internal use only.
 */
static int lapackqr_solvem_under(void * data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_int_t info = 0;
    mess_int_t lwork = 0;
    double *work = NULL;
    mess_double_cpx_t *work_cpx = NULL;
    double workspace;
    mess_double_cpx_t workspace_cpx;


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    if(b->rows!=sol->rows){
        MSG_ERROR("b has wrong number of rows "MESS_PRINTF_INT"!="MESS_PRINTF_INT"\n",b->rows,sol->rows);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/

    if(sol->cpx){
        mess_matrix_copy(b,x);
        mess_matrix_tocomplex(x);
        mess_matrix_resize(x,sol->cols,x->cols);

        //solve with L
        F77_GLOBAL(ztrtrs,ZTRTRS)("L","N","N",&(sol->rows),&(x->cols),sol->val_cpx,&(sol->rows),x->values_cpx,&(x->ld),&info);

        //workspace query and real work
        lwork = -1;
        F77_GLOBAL(zunmlq,ZUNMLQ)("L","C",&(sol->cols),&(x->cols),&(sol->rows),sol->val_cpx,&(sol->rows),sol->tau_cpx,x->values_cpx,&(x->ld),&(workspace_cpx),&lwork,&(info));
        lwork =  nearbyint(creal(workspace_cpx)+1);
        mess_try_alloc(work_cpx,mess_double_cpx_t*, sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zunmlq,ZUNMLQ)("L","C",&(sol->cols),&(x->cols),&(sol->rows),sol->val_cpx,&(sol->rows),sol->tau_cpx,x->values_cpx,&(x->ld),work_cpx,&lwork,&(info));
        mess_free(work_cpx);
    }else{
        if(MESS_IS_REAL(b)){
            mess_matrix_copy(b,x);
            mess_matrix_resize(x,sol->cols,x->cols);

            //solve with L
            F77_GLOBAL(dtrtrs,DTRTRS)("L","N","N",&(sol->rows), &(x->cols), sol->val, &(sol->rows), x->values,&(x->ld),&info);

            //workspace query and real work
            lwork = -1;
            F77_GLOBAL(dormlq,DORMLQ)("L","T",&(sol->cols),&(x->cols),&(sol->rows),sol->val,&(sol->rows),sol->tau, x->values,&(x->ld),&(workspace),&lwork,&info);

            lwork = nearbyint(workspace+1);
            mess_try_alloc(work,double*,sizeof(double)*lwork);
            F77_GLOBAL(dormlq,DORMLQ)("L","T",&(sol->cols),&(x->cols),&(sol->rows),sol->val,&(sol->rows),sol->tau, x->values,&(x->ld),work,&lwork,&info);
            mess_free(work);
        }else{
            mess_matrix xr,xc;
            mess_matrix_init(&xr);
            mess_matrix_init(&xc);
            mess_matrix_realpart(b,xr); mess_matrix_resize(xr,sol->cols,xr->cols);
            mess_matrix_imagpart(b,xc); mess_matrix_resize(xc,sol->cols,xc->cols);

            //solve with L
            F77_GLOBAL(dtrtrs,DTRTRS)("L","N","N",&sol->rows, &(xr->cols), sol->val, &(sol->rows), xr->values,&(xr->ld),&info);
            F77_GLOBAL(dtrtrs,DTRTRS)("L","N","N",&sol->rows, &(xc->cols), sol->val, &(sol->rows), xc->values,&(xc->ld),&info);

            //workspace query and real work
            lwork = -1;
            F77_GLOBAL(dormlq,DORMLQ)("L","T",&(sol->cols),&(xr->cols),&(sol->rows),sol->val,&(sol->rows),sol->tau, xr->values,&(xr->ld),&(workspace),&lwork,&info);
            lwork = nearbyint(workspace+1);
            mess_try_alloc(work,double*,sizeof(double)*lwork);
            F77_GLOBAL(dormlq,DORMLQ)("L","T",&(sol->cols),&(xr->cols),&(sol->rows),sol->val,&(sol->rows),sol->tau, xr->values,&(xr->ld),work,&lwork,&info);
            F77_GLOBAL(dormlq,DORMLQ)("L","T",&(sol->cols),&(xc->cols),&(sol->rows),sol->val,&(sol->rows),sol->tau, xc->values,&(xc->ld),work,&lwork,&info);
            mess_free(work);
            mess_matrix_complex_from_parts(xr,xc,x);
            mess_matrix_clear(&xr);
            mess_matrix_clear(&xc);
        }
    }
    //mess_vector_resize(x, sol->cols);

    if ( info < 0) {
        MSG_ERROR("Error calling DORMLQ/DTRTRS or ZUNMLQ/ZTRTRS. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }
    return (int)(info);
}


/**
 * @internal
 * @brief Solve \f$ Ax=(LQ)^Tx=b (matrix version) \f$.
 * @param[in] data   input data of the solver
 * @param[in] b      input right hand side
 * @param[in,out] x  input/output solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref lapackqr_solve_under function solves the underdetermined system
 * \f[ Ax=b \Leftrightarrow (LQ)^Tx=b.\f]
 *
 * @attention Internal use only.
 */
static int lapackqr_solvemt_under(void * data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_int_t info = 0;
    mess_int_t lwork = 0;
    mess_int_t ret =0;
    double *work = NULL;
    mess_double_cpx_t *work_cpx = NULL;
    double workspace;
    mess_double_cpx_t workspace_cpx;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    if(b->rows!=sol->cols){
        MSG_ERROR("b has wrong dimension "MESS_PRINTF_INT"!="MESS_PRINTF_INT"\n",b->rows,sol->cols);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/

    if(sol->cpx){
        mess_matrix_copy(b,x);
        mess_matrix_tocomplex(x);

        //workspace query and real work
        ret = mess_matrix_conj(x);                                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_conj);
        lwork = -1;
        F77_GLOBAL(zunmlq,ZUNMLQ)("L","N",&(sol->cols),&(x->cols),&(sol->rows),sol->val_cpx,&(sol->rows),sol->tau_cpx,x->values_cpx,&(x->ld),&(workspace_cpx),&lwork,&(info));
        lwork =  nearbyint(creal(workspace_cpx)+1);
        mess_try_alloc(work_cpx,mess_double_cpx_t*, sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zunmlq,ZUNMLQ)("L","N",&(sol->cols),&(x->cols),&(sol->rows),sol->val_cpx,&(sol->rows),sol->tau_cpx,x->values_cpx,&(x->ld),work_cpx,&lwork,&(info));
        mess_free(work_cpx);
        ret = mess_matrix_conj(x);                                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_conj);

        //solve with L
        F77_GLOBAL(ztrtrs,ZTRTRS)("L","T","N",&(sol->rows),&(x->cols),sol->val_cpx,&(sol->rows),x->values_cpx,&(x->ld),&info);
        ret = mess_matrix_resize(x,sol->rows,x->cols);

   }else{
        if(MESS_IS_REAL(b)){
            mess_matrix_copy(b,x);

            //workspace query and real work
            lwork = -1;
            F77_GLOBAL(dormlq,DORMLQ)("L","N",&(sol->cols),&(x->cols),&(sol->rows),sol->val,&(sol->rows),sol->tau, x->values,&(x->ld),&(workspace),&lwork,&info);
            lwork = nearbyint(workspace+1);
            mess_try_alloc(work,double*,sizeof(double)*lwork);
            F77_GLOBAL(dormlq,DORMLQ)("L","N",&(sol->cols),&(x->cols),&(sol->rows),sol->val,&(sol->rows),sol->tau, x->values,&(x->ld),work,&lwork,&info);
            mess_free(work);

            //solve with L
            ret = mess_matrix_resize(x, sol->rows,x->cols);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_resize);
            F77_GLOBAL(dtrtrs,DTRTRS)("L","T","N",&(sol->rows),&(x->cols),sol->val,&(sol->rows),x->values,&(x->ld),&info);

        }else{
            mess_matrix xr,xc;
            MESS_INIT_MATRICES(&xr,&xc);
            mess_matrix_realpart(b,xr);
            mess_matrix_imagpart(b,xc);

            //workspace query and real work
            lwork = -1;
            F77_GLOBAL(dormlq,DORMLQ)("L","N",&(sol->cols),&(xr->cols),&(sol->rows),sol->val,&(sol->rows),sol->tau, xr->values,&(xr->ld),&(workspace),&lwork,&info);
            lwork = nearbyint(workspace+1);
            mess_try_alloc(work,double*,sizeof(double)*lwork);
            F77_GLOBAL(dormlq,DORMLQ)("L","N",&(sol->cols),&(xr->cols),&(sol->rows),sol->val,&(sol->rows),sol->tau, xr->values,&(xr->ld),work,&lwork,&info);
            F77_GLOBAL(dormlq,DORMLQ)("L","N",&(sol->cols),&(xr->cols),&(sol->rows),sol->val,&(sol->rows),sol->tau, xc->values,&(xc->ld),work,&lwork,&info);
            mess_free(work);

            //solve with L
            F77_GLOBAL(dtrtrs,DTRTRS)("L","T","N",&(sol->rows),&(xr->cols),sol->val,&(sol->rows),xr->values,&(xr->ld),&info);
            F77_GLOBAL(dtrtrs,DTRTRS)("L","T","N",&(sol->rows),&(xr->cols),sol->val,&(sol->rows),xc->values,&(xc->ld),&info);

            mess_matrix_complex_from_parts(xr,xc,x);
            MESS_CLEAR_MATRICES(&xr,&xc);
            mess_matrix_resize(x,sol->rows,x->cols);

        }
    }

    if ( info < 0) {
        MSG_ERROR("Error calling DORMLQ/DTRTRS or ZUNMLQ/ZTRTRS. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }
    return (int)(info);
}


/**
 * @internal
 * @brief Solve \f$ Ax=(LQ)^Hx=b \f$.
 * @param[in] data   input data of the solver
 * @param[in] b      input right hand side
 * @param[in,out] x  input/output solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref lapackqr_solve_under function solves the underdetermined system
 * \f[ Ax=b \Leftrightarrow (LQ)^Hx=b.\f]
 *
 * @attention Internal use only.
 */
static int lapackqr_solvemh_under(void * data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_int_t info = 0;
    mess_int_t lwork = 0;
    mess_int_t ret =0;
    mess_double_cpx_t *work_cpx = NULL;
    mess_double_cpx_t workspace_cpx;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    if(b->rows!=sol->cols){
        MSG_ERROR("b has wrong dimension "MESS_PRINTF_INT"!="MESS_PRINTF_INT"\n",b->rows,sol->cols);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/

    if(sol->cpx){
        mess_matrix_copy(b,x);
        mess_matrix_tocomplex(x);

        //workspace query and real work
        lwork = -1;
        F77_GLOBAL(zunmlq,ZUNMLQ)("L","N",&(sol->cols),&(x->cols),&(sol->rows),sol->val_cpx,&(sol->rows),sol->tau_cpx,x->values_cpx,&(x->ld),&(workspace_cpx),&lwork,&(info));
        lwork =  nearbyint(creal(workspace_cpx)+1);
        mess_try_alloc(work_cpx,mess_double_cpx_t*, sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zunmlq,ZUNMLQ)("L","N",&(sol->cols),&(x->cols),&(sol->rows),sol->val_cpx,&(sol->rows),sol->tau_cpx,x->values_cpx,&(x->ld),work_cpx,&lwork,&(info));
        mess_free(work_cpx);

        //solve with L
        F77_GLOBAL(ztrtrs,ZTRTRS)("L","C","N",&(sol->rows),&(x->cols),sol->val_cpx,&(sol->rows),x->values_cpx,&(x->ld),&info);
        ret = mess_matrix_resize(x,sol->rows,x->cols);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_resize);

   }else{
       return lapackqr_solvemt_under(data,b,x);
    }

    if ( info < 0) {
        MSG_ERROR("Error calling DORMLQ/DTRTRS or ZUNMLQ/ZTRTRS. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }
    return (int)(info);
}


/**
 * @internal
 * @brief Get the L (the first factor) of a solver.
 * @param[in] data  input solver data
 * @param[in,out] L     output L factor
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapackqr_getL_under function gets the \f$ L \f$ factor of \f$ A =LQ \f$. Output is in @ref MESS_DENSE
 * format.
 * @attention Internal use only.
 */
static int lapackqr_getL_under(void *data, mess_matrix L ){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_check_nullpointer(sol);
    mess_check_nullpointer(L);
    mess_int_t rows = sol->rows;
    int ret=0;


    ret = mess_matrix_reset(L);                                                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_reset);
    ret = mess_matrix_alloc(L,rows,rows,rows*rows,MESS_DENSE,sol->cpx? MESS_COMPLEX:MESS_REAL); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);

    if(MESS_IS_COMPLEX(L)){
        F77_GLOBAL(zlacpy,ZLACPY)("L",&(sol->rows),&(sol->rows),sol->val_cpx,&(sol->rows),L->values_cpx,&(L->ld));
    }else{
        F77_GLOBAL(dlacpy,DLACPY)("L",&(sol->rows),&(sol->rows),sol->val,&(sol->rows),L->values,&(L->ld));
    }
    return ret;
}


/**
 * @internal
 * @brief Get the U (the second factor) of a solver.
 * @param[in] data  input solver data
 * @param[in,out] U     output U factor
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapackqr_getU_under function gets the \f$ Q \f$ factor of \f$ A =LQ \f$. Output is in @ref MESS_DENSE
 * format.
 * @attention Internal use only.
 *
 */
static int lapackqr_getU_under(void *data, mess_matrix U ){
    MSG_FNAME(__func__);
    struct lapackqr_solver * sol = (struct lapackqr_solver*) data;
    mess_check_nullpointer(sol);
    mess_check_nullpointer(U);
    int ret=0;
    mess_int_t info=0;
    double * work=NULL;
    mess_double_cpx_t *work_cpx = NULL;
    mess_int_t lwork;
    double workspace;
    mess_double_cpx_t workspace_cpx;


    ret = mess_matrix_reset(U);                                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_reset);
    ret = mess_matrix_alloc(U,sol->rows,sol->cols,sol->rows*sol->cols, MESS_DENSE,sol->cpx? MESS_COMPLEX:MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);

    if ( MESS_IS_REAL(U)) {
        F77_GLOBAL(dlacpy,DLACPY)("All", &(sol->rows), &(sol->cols), sol->val, &(sol->rows), U->values, &(U->ld));

        //work space query and real work
        lwork = -1;
        F77_GLOBAL(dorglq,DORGLQ)(&(sol->rows),&(sol->cols),&(sol->rows),U->values,&(U->ld),sol->tau, &workspace,&lwork,&info);
        lwork = nearbyint(workspace+1);

        mess_try_alloc(work, double*, sizeof(double)* lwork);
        F77_GLOBAL(dorglq,DORGLQ)(&(sol->rows),&(sol->cols),&(sol->rows),U->values,&(U->ld),sol->tau, work,&lwork,&info);
        mess_free(work);
    } else {
        F77_GLOBAL(zlacpy,ZLACPY)("All", &(sol->rows), &(sol->cols), sol->val_cpx, &(sol->rows), U->values_cpx, &(U->ld));

        //work space query and real work
        lwork = -1;
        F77_GLOBAL(zunglq,ZUNGLQ)(&(sol->rows),&(sol->cols),&(sol->rows),U->values_cpx,&(U->ld),sol->tau_cpx,&workspace_cpx,&lwork,&info);
        lwork = nearbyint(creal(workspace_cpx)+1);

        mess_try_alloc(work_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)* lwork);
        F77_GLOBAL(zunglq,ZUNGLQ)(&(sol->rows),&(sol->cols),&(sol->rows),U->values_cpx,&(U->ld),sol->tau_cpx,work_cpx,&lwork,&info);
        mess_free(work_cpx);
    }

    return ret;
}


/**
 * @internal
 * @brief Construct a solver based on LQ-decomposition from @lapack for underdetermined Systems.
 * @param[in] matrix        input system matrix
 * @param[out] solver       output
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_create_lapack_qr_under function generates a solver based on the LQ-decomposition
 * from @lapack for underdetermined systems, i.e. the number of cols is greater than the number of
 * rows.
 * @brief Internal use only.
 *
 **/
static int mess_direct_create_lapack_qr_under(mess_matrix matrix, mess_direct solver){
    MSG_FNAME(__func__);
    mess_int_t info = 0;
    mess_int_t lwork = 0;
    mess_int_t ret = 0;
    double *work=NULL, workspace;
    mess_double_cpx_t *work_cpx = NULL, workspace_cpx;
    struct lapackqr_solver *sol;
    mess_int_t conv = 0;
    mess_matrix mwork;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_real_or_complex(matrix);
    if(matrix->rows>=matrix->cols){
        MSG_ERROR("Input matrix has not less rows than cols. Please provide an underdetermined system.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  create solver
     *-----------------------------------------------------------------------------*/
    MESS_MATRIX_CHECKFORMAT(matrix,mwork,conv,MESS_DENSE);
    mess_try_alloc(sol, struct lapackqr_solver*, sizeof(struct lapackqr_solver));
    sol->rows = matrix->rows;
    sol->cols = matrix->cols;

    if(MESS_IS_REAL(matrix)){
        MSG_INFO("real\n");
        sol->val_cpx = NULL;
        sol->tau_cpx = NULL;
        mess_try_alloc(sol->tau, double*, sizeof(double)*(MESS_MIN(sol->rows,sol->cols)));
        mess_try_alloc(sol->val, double*, sizeof(double)*sol->rows*sol->cols);
        F77_GLOBAL(dlacpy,DLACPY)("All",&(matrix->rows),&(matrix->cols),mwork->values,&(mwork->ld),sol->val,&(sol->rows));

        /*-----------------------------------------------------------------------------
         *  workspace query and decomposition
         *-----------------------------------------------------------------------------*/
        lwork = -1;
        F77_GLOBAL(dgelqf,DGELQF)(&(sol->rows), &(sol->cols), sol->val, &(sol->rows), sol->tau, &workspace, &lwork, &info);
        lwork = nearbyint(workspace+1);
        MSG_INFO("Workspace:" MESS_PRINTF_INT"\n",lwork);
        mess_try_alloc(work,double*,sizeof(double)*lwork);
        F77_GLOBAL(dgelqf,DGELQF)(&(sol->rows),&(sol->cols),sol->val,&(sol->rows),sol->tau,work,&lwork,&info);
        mess_free(work);
        sol->cpx=0;
    }else{
        MSG_INFO("complex\n");
        sol->val = NULL;
        sol->tau = NULL;
        mess_try_alloc(sol->tau_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*(MESS_MIN(sol->rows,sol->cols)));
        mess_try_alloc(sol->val_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*sol->rows*sol->cols);
        F77_GLOBAL(zlacpy,ZLACPY)("All",&(matrix->rows),&(matrix->cols),mwork->values_cpx,&(mwork->ld),sol->val_cpx,&(sol->rows));

        /*-----------------------------------------------------------------------------
         *  workspace query and decomposition
         *-----------------------------------------------------------------------------*/
        lwork = -1;
        F77_GLOBAL(zgelqf,ZGELQF)(&(sol->rows), &(sol->cols), sol->val_cpx, &(sol->rows), sol->tau_cpx, &workspace_cpx, &lwork, &info);
        lwork = nearbyint(creal(workspace_cpx)+1);
        MSG_INFO("Workspace:" MESS_PRINTF_INT"\n",lwork);
        mess_try_alloc(work_cpx,mess_double_cpx_t*,sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zgelqf,ZGELQF)(&(sol->rows),&(sol->cols),sol->val_cpx,&(sol->rows),sol->tau_cpx,work_cpx,&lwork,&info);
        mess_free(work_cpx);
        sol->cpx=1;
    }

    solver->data_type = matrix->data_type;
    solver->data = (void *) sol;
    solver->clear = lapackqr_clear;
    solver->solve = lapackqr_solve_under;
    solver->solvet = lapackqr_solvet_under;
    solver->solveh = lapackqr_solveh_under;
    solver->solvem = lapackqr_solvem_under;
    solver->solvemt = lapackqr_solvemt_under;
    solver->solvemh = lapackqr_solvemh_under;
    solver->getL = lapackqr_getL_under;
    solver->getU = lapackqr_getU_under;
    solver->getpermp = NULL;
    solver->getpermq = NULL;
    SET_SOLVERNAME(solver->name, __func__);

    if ( conv == 0 ) {
        mess_matrix_clear(&mwork);
    }
    if ( info != 0 ) {
        MSG_ERROR("An error occured in DGELQF/ZGELQF: " MESS_PRINTF_INT "\n", info);
        ret = (int)(info);
    }

    return ret;
}


/*-----------------------------------------------------------------------------
 *  Callable function for QR / LQ decomposition
 *-----------------------------------------------------------------------------*/
/**
 * @brief   Construct a solver based on QR or LQ decomposition from @lapack.
 * @param[in]   matrix   input system matrix
 * @param[in,out]   solver  generated solver
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_create_lapack_qr function generates a solver based on the QR or LQ decomposition from
 * @lapack. If the number of rows is greater equals the number of cols, we perform a QR decomposition,
 * otherwise we perform a LQ decomposition.
 * In both cases it is necessary, that @c matrix has full rank.
 *
 * @attention Please keep in mind, that a sparse matrix is converted to a dense one.
 *
 */
int mess_direct_create_lapack_qr (mess_matrix matrix, mess_direct solver){
    MSG_FNAME(__func__);

    /*-----------------------------------------------------------------------------
     *  check input arguements
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_real_or_complex(matrix);
    //mess_check_dense(matrix);

    /*-----------------------------------------------------------------------------
     *  create solver
     *-----------------------------------------------------------------------------*/
    if(matrix->rows<matrix->cols){
        return mess_direct_create_lapack_qr_under(matrix,solver);
    }else{
        return mess_direct_create_lapack_qr_over(matrix,solver);
    }
}


