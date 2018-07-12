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
 * @file lib/direct/singlesolver/cholesky.c
 * @brief Cholesky decomposition based on DPOTRF.
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
#include "mess/error_macro.h"
#include "blas_defs.h"
#ifdef _OPENMP_H
#include <omp.h>
#endif




/**
 * @brief Internal structure for the @lapack based solver.
 *
 */
struct cholesky_solver {
    double *val;                /**< the real lu decomposed matrix  */
    mess_double_cpx_t  *val_cpx;   /**< the complex lu decomposed matrix  */
    mess_int_t n;               /**< dimension of the system */
    unsigned short cpx;         /**< flag real/complex */
};



/**
 * @brief Clear the cholesky solver.
 * @param[in,out] data  pointer to the internal data
 * @return always zero
 *
 * The @ref chol_clear function clears a @lapack sovler.
 *
 */
static int chol_clear(void *data){
    struct cholesky_solver *sol = (struct cholesky_solver *) data;
    if ( sol->val != NULL)          mess_free(sol->val);
    if ( sol->val_cpx != NULL)      mess_free(sol->val_cpx);
    mess_free(sol);
    return 0;
}

/**
 * @brief Get the L factor of the cholesky solver (\f$ A=LL^T \f$)
 * @param[in] data   input pointer to the internal data
 * @param[in,out] L the matrix L
 *
 * The @ref chol_getL function returns the L factor such that \f$ A=LL^T \f$.
 *
 */
static int chol_getL(void *data,mess_matrix L){
    MSG_FNAME(__func__);
    int ret = 0;
    struct cholesky_solver *sol = (struct cholesky_solver *) data;

    mess_check_nullpointer(data);
    mess_check_nullpointer(L);
    MESS_MATRIX_RESET(L);
    if ( sol->cpx ) {
        ret = mess_matrix_alloc(L,sol->n,sol->n, sol->n*sol->n, MESS_DENSE, MESS_COMPLEX);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        F77_GLOBAL(zlacpy,ZLACPY)("L", &(sol->n), &(sol->n), sol->val_cpx, &(sol->n), L->values_cpx , &(L->ld));
    } else {
        ret = mess_matrix_alloc(L,sol->n,sol->n, sol->n*sol->n, MESS_DENSE, MESS_REAL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        F77_GLOBAL(dlacpy,DLACPY)("L", &(sol->n), &(sol->n), sol->val, &(sol->n), L->values, &(L->ld));
    }

    return 0;
}

/**
 * @brief Get the U factor (second factor) of the cholesky solver (\f$ A=LL^T \f$)
 * @param[in] data  input pointer to the internal data
 * @param[out] U    output matrix contains \f$L^T\f$
 *
 * The @ref chol_getU function sets @p U to \f$L^T\f$.
 *
 */
static int chol_getU(void *data,mess_matrix U){
    MSG_FNAME(__func__);
    mess_int_t ret;
    mess_matrix L;
    ret = mess_matrix_init(&L);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = chol_getL(data,L);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),chol_getL);
    ret = mess_matrix_ctranspose(L,U);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    ret = mess_matrix_clear(&L);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    return ret;
}


/**
 * @brief Solve \f$ Ax=b \f$ using @lapack.
 * @param[in] data  input pointer to the internal data structure
 * @param[in] b      input right hand side
 * @param[in,out] x     solution
 *
 * The @ref chol_solve function solves \f$ Ax=b \f$ using @lapack.
 *
 */
static int chol_solve(void * data, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct cholesky_solver * sol = (struct cholesky_solver*) data;
    mess_int_t info = 0;
    mess_int_t one = 1;
    int ret = 0;

    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    ret = mess_vector_copy(b,x);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);

    if (sol->cpx == 0 && MESS_IS_REAL(x)) {
        // solve real system
        F77_GLOBAL(dpotrs,DPOTRS)("L", &sol->n, &one, sol->val, &sol->n, x->values, &sol->n, &info);
    } else if ( sol->cpx == 0 && MESS_IS_COMPLEX(x)){
        F77_GLOBAL(dzpotrs,DZPOTRS)("L", &sol->n, &one, sol->val, &sol->n, x->values_cpx, &sol->n, &info);
    } else {
        // solve complex system
        ret = mess_vector_tocomplex(x);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(zpotrs,ZPOTRS)("L", &sol->n, &one, sol->val_cpx, &sol->n, x->values_cpx, &sol->n, &info);
    }
    if ( info < 0) {
        MSG_ERROR("error calling ZPOTRS/DPOTRS. Invalid argument: " MESS_PRINTF_INT "\n", -info);
    }

    return 0;
}

/**
 * @brief Solve \f$ AX=B \f$ using @lapack (matrix version)
 * @param[in] data  input pointer to the internal data structure
 * @param[in] b      input right hand side
 * @param[in,out] x     solution
 *
 * The @ref chol_solvem function solves \f$ AX=B \f$ using @lapack, where \f$ B \f$ and \f$ X \f$ are matrices.
 *
 */
static int chol_solvem(void * data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct cholesky_solver * sol = (struct cholesky_solver*) data;
    mess_int_t info = 0;
    int ret =0;
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);

    ret = mess_matrix_copy(b,x);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);

    if (sol->cpx == 0 && MESS_IS_REAL(x)) {
        // solve real system
        F77_GLOBAL(dpotrs,DPOTRS)("L", &sol->n, &(x->cols), sol->val, &sol->n, x->values, &x->ld, &info);
    } else if ( sol->cpx == 0 && MESS_IS_COMPLEX(x)){
        F77_GLOBAL(dzpotrs,DZPOTRS)("L", &sol->n, &(x->cols), sol->val, &sol->n, x->values_cpx, &x->ld, &info);

        // F77_GLOBAL(dzgetrs,DZGETRS)("N",&(sol->n), &(x->cols), sol->val, &(sol->n), sol->ipiv, x->values_cpx, &(sol->n),&info);
    } else {
        // solve complex system
        ret = mess_matrix_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        F77_GLOBAL(zpotrs,ZPOTRS)("L", &sol->n, &(x->cols), sol->val_cpx, &sol->n, x->values_cpx, &x->ld, &info);
    }
    if ( info < 0) {
        MSG_ERROR("error calling ZPOTRS/DPOTRS. Invalid argument: " MESS_PRINTF_INT "\n", -info);
    }

    return 0;

}

/**
 * @brief Compute the inverse of a matrix using @lapack.
 * @param[in] data  input pointer to the internal data
 * @param[in,out] inv output inverse
 *
 * The @ref lapack_inverse function computes \f$ A^{-1} \f$ using @lapack.
 *
 */
static int chol_inverse( void *data, mess_matrix inv) {
    MSG_FNAME(__func__);
    struct cholesky_solver * sol = (struct cholesky_solver*) data;
    mess_int_t info = 0;
    mess_check_nullpointer(sol);
    int ret =  0;
    mess_int_t n, j ,i ;
    n = sol->n ;

    if (sol->cpx == 0) {
        ret = mess_matrix_alloc(inv, n,n,n*n, MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
        F77_GLOBAL(dlacpy,DLACPY)("All", &n, &n, sol->val, &sol->n, inv->values, &inv->ld);
        F77_GLOBAL(dpotri,DPOTRI)("L", &(sol->n), inv->values, &(inv->ld), &info);
        for ( j = 0; j<n; j++) {
            for (i=0; i<j; i++) {
                inv->values[i+j*inv->ld] = inv->values[j+i*inv->ld];
            }
        }
    } else {
        ret = mess_matrix_alloc(inv, n,n,n*n, MESS_DENSE, MESS_COMPLEX);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
        F77_GLOBAL(zlacpy,ZLACPY)("All", &n, &n, sol->val_cpx, &sol->n, inv->values_cpx, &inv->ld);
        F77_GLOBAL(zpotri,ZPOTRI)("L", &(sol->n), inv->values_cpx, &(inv->ld), &info);
        for ( j = 0; j<n; j++) {
            for (i=0; i<j; i++) {
                inv->values_cpx[i+j*inv->ld] = conj(inv->values_cpx[j+i*inv->ld]);
            }
        }
    }
    if ( info < 0) {
        MSG_ERROR("error calling DPOTRI/ZPOTRI. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }

    return 0;
}

/**
 * @brief Generate a dense @lapack based Cholesky solver.
 * @param[in] matrix    symmetric positive definite input matrix
 * @param[out] solver   generated solver
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_direct_create_cholesky function generates a Cholesky solver for dense matrices based on @lapack.\n
 * It only works for symmetric positive definite matrices. If the matrix does not fulfill this property an error
 * is thrown.
 *
 *
 */
int mess_direct_create_cholesky (mess_matrix matrix, mess_direct solver){
    MSG_FNAME(__func__);
    mess_int_t info = 0;
    // int ret = 0;
    struct cholesky_solver * sol ;
    mess_matrix mwork;
    int conv = 0;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_square(matrix);
    mess_check_real_or_complex(matrix);

    MESS_MATRIX_CHECKFORMAT(matrix, mwork, conv, MESS_DENSE);

    mess_try_alloc(sol, struct cholesky_solver *, sizeof(struct cholesky_solver));
    sol->n = matrix->rows;

    if ( MESS_IS_COMPLEX ( matrix )) {
        sol->cpx = 1;
        mess_try_alloc(sol->val_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * matrix->rows * matrix->cols);
        sol->val=NULL;
        F77_GLOBAL(zlacpy,ZLACPY)("All", &(sol->n), &(sol->n), mwork->values_cpx, &mwork->ld, sol->val_cpx, &(sol->n));
        F77_GLOBAL(zpotrf,ZPOTRF)("L", &(sol->n), sol->val_cpx, &(sol->n), &info);

        if (info < 0 ) {
            MSG_ERROR("Argument Error in ZPOTRF\n");
            mess_free(sol->val);
            mess_free(sol);
            return MESS_ERROR_LAPACK;
        }
        if (info > 0 ) {
            MSG_ERROR("Input matrix is not SPD\n");
            mess_free(sol->val);
            mess_free(sol);
            return MESS_ERROR_SYMMETRIC;
        }

    } else {
        sol->cpx = 0;
        mess_try_alloc(sol->val, double *, sizeof(double) * matrix->rows * matrix->cols);
        sol->val_cpx=NULL;
        F77_GLOBAL(dlacpy,DLACPY)("All", &(sol->n), &(sol->n), mwork->values, &mwork->ld, sol->val, &(sol->n));
        F77_GLOBAL(dpotrf,DPOTRF)("L", &(sol->n), sol->val, &(sol->n), &info);

        if (info < 0 ) {
            MSG_ERROR("Argument Error in DPOTRF\n");
            mess_free(sol->val);
            mess_free(sol);
            return MESS_ERROR_LAPACK;
        }
        if (info > 0 ) {
            MSG_ERROR("Input matrix is not SPD info = %d\n", (int) info);
            mess_free(sol->val);
            mess_free(sol);
            return MESS_ERROR_SYMMETRIC;
        }
    }

    solver->data = (void *) sol;
    solver->clear = chol_clear;
    solver->solve = chol_solve;
    solver->solvet = chol_solve;
    solver->solveh = chol_solve;
    solver->solvem = chol_solvem;
    solver->solvemt = chol_solvem;
    solver->solvemh = chol_solvem;
    solver->inverse = chol_inverse;
    solver->getL = chol_getL;
    solver->getU = chol_getU;
    solver->getpermp = NULL;
    solver->getpermq = NULL;
    solver->rows = matrix->rows;
    solver->cols = matrix->cols;
    SET_SOLVERNAME(solver->name, __func__);

    if ( conv == 0 ) {
        mess_matrix_clear(&mwork);
    }

    return 0;
}

