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
 * @file lib/direct/singlesolver/banded.c
 * @brief Interface to @lapack or @umfpack based banded solver.
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"

#ifdef MESS_HAVE_UMFPACK
#include <umfpack.h>
#ifdef MESS64
#define umfpack_zl_defaults     umfpack_zl_defaults
#define umfpack_dl_defaults     umfpack_dl_defaults
#else
#define umfpack_zl_defaults     umfpack_zi_defaults
#define umfpack_dl_defaults     umfpack_di_defaults
#endif
#endif




#define BANDWIDHT_SWITCH_THRESHOLD 0.3 // if BANDWIDTH_SWITCH_THRESHOLD < (kl +ku+1)/n -> switch to umfpack solver

#define IDX2F(i,j,ld) (((((j)-1))*(ld))+ (((i)-1)))

#define LAPACK_ERROR(INFO)                                                                                      \
    if ( (INFO) < 0) {                                                                                          \
        MSG_ERROR("error calling DGBTRS/DZGBTRS/ZGBTRS. Invalid argument: " MESS_PRINTF_INT "\n", -(INFO));     \
        return (int) (INFO);                                                                                    \
    }


struct banded_lapack {
    double *AB_real;
    mess_double_cpx_t * AB_cpx;
    mess_int_t n;
    mess_int_t nband;
    mess_int_t kl;
    mess_int_t ku;
    mess_int_t info;
    mess_int_t * perm;
    mess_int_t * invperm;
    mess_int_t * ipiv;
};

/**
 * @brief Solve \f$ Ax=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref banded_lapack_solve function solves \f$ Ax=b \f$.
 *
 */
static int banded_lapack_solve(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    struct banded_lapack * sol = (struct banded_lapack *)data;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->dim != sol->n) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->dim, sol->n);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    mess_int_t nrhs=1;
    if(sol->AB_real && MESS_IS_REAL(b)){
        ret = mess_vector_toreal_nowarn(x);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal_nowarn);
        ret = mess_vector_perm(b,sol->perm,x);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm);
        F77_GLOBAL(dgbtrs,DGBTRS)("N", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_real, &(sol->nband), sol->ipiv, x->values, &(sol->n), &(sol->info));
    }else if(sol->AB_real){
        ret = mess_vector_tocomplex(x);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        ret = mess_vector_perm(b,sol->perm,x);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm);
        F77_GLOBAL(dzgbtrs,DZGBTRS)("N", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_real, &(sol->nband), sol->ipiv, x->values_cpx, &(sol->n), &(sol->info));
    }else if(sol->AB_cpx){
        ret = mess_vector_tocomplex(x);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        ret = mess_vector_perm(b,sol->perm,x);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm);
        F77_GLOBAL(zgbtrs,ZGBTRS)("N", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_cpx, &(sol->nband), sol->ipiv, x->values_cpx, &(sol->n), &(sol->info));
    }
    ret = mess_vector_perm_inplace(x,sol->invperm);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm_inplace);
    LAPACK_ERROR(sol->info);

    return ret;
}

/**
 * @brief Solve \f$ A^Tx=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref banded_lapack_solvet function solves \f$ A^Tx=b \f$.
 *
 */
static int banded_lapack_solvet(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    struct banded_lapack * sol = (struct banded_lapack *)data;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->dim != sol->n) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->dim, sol->n);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    mess_int_t nrhs=1;
    if(sol->AB_real && MESS_IS_REAL(b)){
        ret = mess_vector_toreal_nowarn(x);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal_nowarn);
        ret = mess_vector_perm(b,sol->perm,x);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm);
        F77_GLOBAL(dgbtrs,DGBTRS)("T", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_real, &(sol->nband), sol->ipiv, x->values, &(sol->n), &(sol->info));
    }else if(sol->AB_real){
        ret = mess_vector_tocomplex(x);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        ret = mess_vector_perm(b,sol->perm,x);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm);
        F77_GLOBAL(dzgbtrs,DZGBTRS)("T", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_real, &(sol->nband), sol->ipiv, x->values_cpx, &(sol->n), &(sol->info));
    }else if(sol->AB_cpx){
        ret = mess_vector_tocomplex(x);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        ret = mess_vector_perm(b,sol->perm,x);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm);
        F77_GLOBAL(zgbtrs,ZGBTRS)("T", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_cpx, &(sol->nband), sol->ipiv, x->values_cpx, &(sol->n), &(sol->info));
    }
    ret = mess_vector_perm_inplace(x,sol->invperm); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mss_vector_perm_inplace);

    LAPACK_ERROR(sol->info);

    return ret;
}

/**
 * @brief Solve \f$ A^Hx=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref banded_lapack_solveh function solves \f$ A^Hx=b \f$.
 *
 */
static int banded_lapack_solveh(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    struct banded_lapack * sol = (struct banded_lapack *)data;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->dim != sol->n) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->dim, sol->n);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    mess_int_t nrhs=1;
    if(sol->AB_real && MESS_IS_REAL(b)){
        ret = mess_vector_toreal_nowarn(x);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal_nowarn);
        ret = mess_vector_perm(b,sol->perm,x);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm);
        F77_GLOBAL(dgbtrs,DGBTRS)("T", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_real, &(sol->nband), sol->ipiv, x->values, &(sol->n), &(sol->info));
    }else if(sol->AB_real){
        ret = mess_vector_tocomplex(x);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        ret = mess_vector_perm(b,sol->perm,x);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm);
        F77_GLOBAL(dzgbtrs,DZGBTRS)("T", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_real, &(sol->nband), sol->ipiv, x->values_cpx, &(sol->n), &(sol->info));
    }else if(sol->AB_cpx){
        ret = mess_vector_tocomplex(x);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        ret = mess_vector_perm(b,sol->perm,x);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm);
        F77_GLOBAL(zgbtrs,ZGBTRS)("C", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_cpx, &(sol->nband), sol->ipiv, x->values_cpx, &(sol->n), &(sol->info));
    }
    ret = mess_vector_perm_inplace(x,sol->invperm);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm_inplace);

    LAPACK_ERROR(sol->info);

    return ret;
}

/**
 * @brief Solve \f$ AX=B \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref banded_lapack_solvem function solves \f$ AX=B \f$.
 *
 */
static int banded_lapack_solvem(void*data, mess_matrix b, mess_matrix x) {
    MSG_FNAME(__func__);
    struct banded_lapack * sol = (struct banded_lapack *)data;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->rows != sol->n) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->rows, sol->n);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    mess_int_t nrhs=b->cols;
    if(sol->AB_real && MESS_IS_REAL(b)){
        ret = mess_matrix_copy(b,x);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_perm(x,sol->perm,NULL);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
        F77_GLOBAL(dgbtrs,DGBTRS)("N", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_real, &(sol->nband), sol->ipiv, x->values, &(sol->n), &(sol->info));
    }else if(sol->AB_real){
        ret = mess_matrix_copy(b,x);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_perm(x,sol->perm,NULL);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
        F77_GLOBAL(dzgbtrs,DZGBTRS)("N", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_real, &(sol->nband), sol->ipiv, x->values_cpx, &(sol->n), &(sol->info));
    }else if(sol->AB_cpx){
        ret = mess_matrix_copy(b,x);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_perm(x,sol->perm,NULL);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
        F77_GLOBAL(zgbtrs,ZGBTRS)("N", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_cpx, &(sol->nband), sol->ipiv, x->values_cpx, &(sol->n), &(sol->info));
    }
    ret = mess_matrix_perm(x,sol->invperm,NULL);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
    LAPACK_ERROR(sol->info);

    return ret;
}


/**
 * @brief Solve \f$ A^TX=B \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref banded_lapack_solvemt function solves \f$ A^TX=B \f$.
 *
 */
static int banded_lapack_solvemt(void*data, mess_matrix b, mess_matrix x) {
    MSG_FNAME(__func__);
    struct banded_lapack * sol = (struct banded_lapack *)data;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->rows != sol->n) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->rows, sol->n);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    mess_int_t nrhs=b->cols;
    if(sol->AB_real && MESS_IS_REAL(b)){
        ret = mess_matrix_copy(b,x);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_perm(x,sol->perm,NULL);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
        F77_GLOBAL(dgbtrs,DGBTRS)("T", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_real, &(sol->nband), sol->ipiv, x->values, &(sol->n), &(sol->info));
    }else if(sol->AB_real){
        ret = mess_matrix_copy(b,x);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_perm(x,sol->perm,NULL);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
        F77_GLOBAL(dzgbtrs,DZGBTRS)("T", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_real, &(sol->nband), sol->ipiv, x->values_cpx, &(sol->n), &(sol->info));
    }else if(sol->AB_cpx){
        ret = mess_matrix_copy(b,x);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_perm(x,sol->perm,NULL);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
        F77_GLOBAL(zgbtrs,ZGBTRS)("T", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_cpx, &(sol->nband), sol->ipiv, x->values_cpx, &(sol->n), &(sol->info));
    }
    ret = mess_matrix_perm(x,sol->invperm,NULL);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);

    LAPACK_ERROR(sol->info);

    return ret;
}


/**
 * @brief Solve \f$ A^HX=B \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref banded_lapack_solvemh function solves \f$ A^HX=B \f$.
 *
 */
static int banded_lapack_solvemh(void*data, mess_matrix b, mess_matrix x) {
    MSG_FNAME(__func__);
    struct banded_lapack * sol = (struct banded_lapack *)data;
    int ret=0;
    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->rows != sol->n) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->rows, sol->n);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    mess_int_t nrhs=b->cols;
    if(sol->AB_real && MESS_IS_REAL(b)){
        ret = mess_matrix_copy(b,x);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_perm(x,sol->perm,NULL);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
        F77_GLOBAL(dgbtrs,DGBTRS)("T", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_real, &(sol->nband), sol->ipiv, x->values, &(sol->n), &(sol->info));
    }else if(sol->AB_real){
        ret = mess_matrix_copy(b,x);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_perm(x,sol->perm,NULL);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
        F77_GLOBAL(dzgbtrs,DZGBTRS)("T", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_real, &(sol->nband), sol->ipiv, x->values_cpx, &(sol->n), &(sol->info));
    }else if(sol->AB_cpx){
        ret = mess_matrix_copy(b,x);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_perm(x,sol->perm,NULL);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
        F77_GLOBAL(zgbtrs,ZGBTRS)("C", &(sol->n), &(sol->kl), &(sol->ku), &nrhs, sol->AB_cpx, &(sol->nband), sol->ipiv, x->values_cpx, &(sol->n), &(sol->info));
    }
    ret = mess_matrix_perm(x,sol->invperm,NULL);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);

    LAPACK_ERROR(sol->info);

    return ret;
}

/**
 * @brief Clear the banded solver
 * @param [in,out] data     pointer to the internal data
 * @return always zero
 *
 * The @ref banded_lapack_clear function clears a banded solver.
 *
 */
static int banded_lapack_clear(void *data){
    struct banded_lapack *sol = (struct banded_lapack *) data;
    if ( sol->AB_real)    mess_free(sol->AB_real);
    if ( sol->AB_cpx)    mess_free(sol->AB_cpx);
    mess_free(sol->ipiv);
    mess_free(sol->perm);
    mess_free(sol->invperm);
    mess_free(sol);
    return 0;
}


#if defined(MESS_HAVE_UMFPACK) && defined(UMFPACK_ORDERING_NONE)
struct banded_umfpack{
    mess_direct umfpack_solver;
    mess_int_t* perm;
    mess_int_t* invperm;
    mess_int_t n;
};

/**
 * @brief Clear the banded solver.
 * @param [in,out] data     pointer to the internal data
 * @return always zero
 *
 * The @ref banded_umfpack_clear function clears a banded solver.
 *
 */
static int banded_umfpack_clear(void* data){
    struct banded_umfpack *sol = (struct banded_umfpack *) data;
    mess_direct_clear(&(sol->umfpack_solver));
    mess_free(sol->perm);
    mess_free(sol->invperm);
    mess_free(sol);

    return 0;
}

/**
 * @brief Solve \f$ Ax=b \f$.
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref banded_umfpack_solve function solves \f$ Ax=b \f$.
 *
 */
static int banded_umfpack_solve(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    struct banded_umfpack * sol = (struct banded_umfpack *)data;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->dim != sol->n) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->dim, sol->n);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_perm_inplace(b,sol->perm);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm_inplace);
    ret = mess_direct_solve(MESS_OP_NONE,sol->umfpack_solver,b,x);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solve);
    ret = mess_vector_perm_inplace(b,sol->invperm);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm_inplace);
    ret = mess_vector_perm_inplace(x,sol->invperm);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm_inplace);
    return ret;
}

/**
 * @brief Solve \f$ A^Tx=b \f$.
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref banded_umfpack_solve function solves \f$ A^Tx=b \f$.
 *
 */
static int banded_umfpack_solvet(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    struct banded_umfpack * sol = (struct banded_umfpack *)data;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->dim != sol->n) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->dim, sol->n);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_perm_inplace(b,sol->perm);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm_inplace);
    ret = mess_direct_solve(MESS_OP_TRANSPOSE,sol->umfpack_solver,b,x);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solve);
    ret = mess_vector_perm_inplace(b,sol->invperm);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm_inplace);
    ret = mess_vector_perm_inplace(x,sol->invperm);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm_inplace);
    return ret;
}

/**
 * @brief Solve \f$ A^Hx=b \f$.
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref banded_umfpack_solve function solves \f$ A^Hx=b \f$.
 *
 */
static int banded_umfpack_solveh(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    struct banded_umfpack * sol = (struct banded_umfpack *)data;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->dim != sol->n) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->dim, sol->n);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
        MSG_INFO("Build banded solver based on LAPACK");
     *  solve system
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_perm_inplace(b,sol->perm);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm_inplace);
    ret = mess_direct_solve(MESS_OP_HERMITIAN,sol->umfpack_solver,b,x);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solve);
    ret = mess_vector_perm_inplace(b,sol->invperm);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm_inplace);
    ret = mess_vector_perm_inplace(x,sol->invperm);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_perm_inplace);
    return ret;
}

/**
 * @brief Solve \f$ AX=B \f$ (matrix version).
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref banded_umfpack_solve function solves \f$ AX=B \f$.
 *
 */
static int banded_umfpack_solvem(void*data, mess_matrix b, mess_matrix x) {
    MSG_FNAME(__func__);
    struct banded_umfpack * sol = (struct banded_umfpack *)data;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->rows!= sol->n) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->rows, sol->n);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_perm(b,sol->perm,NULL);                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
        MSG_INFO("Build banded solver based on LAPACK.\n");
    ret = mess_direct_solvem(MESS_OP_NONE,sol->umfpack_solver,b,x);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solve);
    ret = mess_matrix_perm(b,sol->invperm,NULL);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
    ret = mess_matrix_perm(x,sol->invperm,NULL);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
    return ret;
}

/**
 * @brief Solve \f$ A^TX=B \f$ (matrix version).
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref banded_umfpack_solve function solves \f$ A^TX=B \f$.
 *
 */
static int banded_umfpack_solvemt(void*data, mess_matrix b, mess_matrix x) {
    MSG_FNAME(__func__);
    struct banded_umfpack * sol = (struct banded_umfpack *)data;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->rows!= sol->n) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->rows, sol->n);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_perm(b,sol->perm,NULL);                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
    ret = mess_direct_solvem(MESS_OP_TRANSPOSE,sol->umfpack_solver,b,x);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solve);
    ret = mess_matrix_perm(b,sol->invperm,NULL);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
    ret = mess_matrix_perm(x,sol->invperm,NULL);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
    return ret;
}

/**
 * @brief Solve \f$ A^HX=B \f$ (matrix version).
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref banded_umfpack_solve function solves \f$ A^HX=B \f$.
 *
 */
static int banded_umfpack_solvemh(void*data, mess_matrix b, mess_matrix x) {
    MSG_FNAME(__func__);
    struct banded_umfpack * sol = (struct banded_umfpack *)data;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->rows!= sol->n) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->rows, sol->n);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_perm(b,sol->perm,NULL);                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
    ret = mess_direct_solvem(MESS_OP_HERMITIAN,sol->umfpack_solver,b,x);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solve);
    ret = mess_matrix_perm(b,sol->invperm,NULL);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
    ret = mess_matrix_perm(x,sol->invperm,NULL);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
    return ret;
}
#endif

/**
 * @brief Generate a direct linear banded solver for standard linear systems \f$ Ax=b \f$ with LAPACK
 * @param[in] matrix matrix to decompose
 * @param[out] solver output solver
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref mess_direct_create_banded function uses bandwidth reducing reordering
 * and @lapack banded solvers
 * or @umfpack solver for solving linear systems. If bandwidth after reordering to large
 * we switch to @umfpack based solver (if available).
 *
 */
int mess_direct_create_banded(mess_matrix matrix, mess_direct solver){
    MSG_FNAME(__func__);
    mess_int_t i, j, kl, ku, nband, *perm, *invperm;
    int ret;
    mess_matrix temp;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_square(matrix);
    mess_check_real_or_complex(matrix);

    /*-----------------------------------------------------------------------------
     * reorder matrix and compute inverse permutation
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&temp);                                                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_convert(matrix,temp,MESS_CSR);                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert);
    mess_try_alloc(perm, mess_int_t*, sizeof(mess_int_t)*(matrix->rows));
    mess_try_alloc(invperm, mess_int_t*, sizeof(mess_int_t)*(matrix->rows));
    ret = mess_matrix_reorder_rcm(temp,perm);                                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_reorder_rcm);
    ret = mess_matrix_perm(temp,perm,perm);                                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_perm);
    for(i=0;i<matrix->rows;++i){invperm[perm[i]]=i;}
    ret = mess_matrix_bandwidth(temp,&(kl),&(ku));                                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_bandwidth);
    nband = 1+2*kl+ku;
    MSG_INFO("dim=" MESS_PRINTF_INT" kl=" MESS_PRINTF_INT" ku=" MESS_PRINTF_INT"\n",matrix->rows,kl,ku);

#if defined(MESS_HAVE_UMFPACK) && defined(UMFPACK_ORDERING_NONE)
    //decide umpfack lapack solver
    if(BANDWIDHT_SWITCH_THRESHOLD < nband/(double)temp->rows){
#endif
        //lapack solver
        MSG_INFO("Build banded solver based on LAPACK.\n");
        struct banded_lapack *lapack_sol;
        mess_try_alloc(lapack_sol, struct banded_lapack*, sizeof(struct banded_lapack));
        lapack_sol->AB_real=NULL;
        lapack_sol->AB_cpx=NULL;
        lapack_sol->n = matrix->rows;
        lapack_sol->nband=nband;
        lapack_sol->kl=kl;
        lapack_sol->ku=ku;
        lapack_sol->perm=perm;
        lapack_sol->invperm=invperm;
        mess_try_alloc(lapack_sol->ipiv,mess_int_t*,sizeof(mess_int_t)*lapack_sol->n);

        /*-----------------------------------------------------------------------------
         *  Bandwidth convert and lu decomposition
         *-----------------------------------------------------------------------------*/
        if(MESS_IS_REAL(temp)){
            mess_try_alloc(lapack_sol->AB_real,double*,sizeof(double)*(lapack_sol->nband)*(lapack_sol->n));
            for(i=0;i<(lapack_sol->n)*(lapack_sol->nband);++i){lapack_sol->AB_real[i]=0;}
            for(i=0;i<temp->rows;++i){
                for(j=temp->rowptr[i];j<temp->rowptr[i+1];++j){
                    lapack_sol->AB_real[IDX2F(lapack_sol->kl + lapack_sol->ku + 1 + (i+1)-(temp->colptr[j]+1),temp->colptr[j]+1,lapack_sol->nband)]=temp->values[j];
                }
            }
            F77_GLOBAL(dgbtrf,DGBTRF)(&(lapack_sol->n),&(lapack_sol->n),&(lapack_sol->kl),&(lapack_sol->ku),lapack_sol->AB_real,&(lapack_sol->nband),lapack_sol->ipiv,&(lapack_sol->info));
        }else{
            mess_try_alloc(lapack_sol->AB_cpx,mess_double_cpx_t*,sizeof(mess_double_cpx_t)*(lapack_sol->nband)*(lapack_sol->n));
            for(i=0;i<(lapack_sol->n)*(lapack_sol->nband);++i){lapack_sol->AB_cpx[i]=0;}
            for(i=0;i<temp->rows;++i){
                for(j=temp->rowptr[i];j<temp->rowptr[i+1];++j){
                    lapack_sol->AB_cpx[IDX2F(lapack_sol->kl + lapack_sol->ku + 1 + (i+1)-(temp->colptr[j]+1),temp->colptr[j]+1,lapack_sol->nband)]=temp->values_cpx[j];
                }
            }
            F77_GLOBAL(zgbtrf,ZGBTRF)(&(lapack_sol->n),&(lapack_sol->n),&(lapack_sol->kl),&(lapack_sol->ku),lapack_sol->AB_cpx,&(lapack_sol->nband),lapack_sol->ipiv,&(lapack_sol->info));
        }
        LAPACK_ERROR(lapack_sol->info);
        /*-----------------------------------------------------------------------------
         *  set pointers
         *-----------------------------------------------------------------------------*/
        solver->data = (void *) lapack_sol;
        solver->clear = banded_lapack_clear;
        solver->solve = banded_lapack_solve;
        solver->solvet = banded_lapack_solvet;
        solver->solveh = banded_lapack_solveh;
        solver->solvem = banded_lapack_solvem;
        solver->solvemt = banded_lapack_solvemt;
        solver->solvemh = banded_lapack_solvemh;
        solver->inverse = NULL;
        solver->det  = NULL;
        solver->detc  = NULL;
        solver->getL = NULL;
        solver->getU = NULL;
        solver->getpermp = NULL;
        solver->getpermq = NULL;
        solver->rows = matrix->rows;
        solver->cols = matrix->cols;
        solver->data_type = temp->data_type;
        SET_SOLVERNAME(solver->name, __func__);

#if defined(MESS_HAVE_UMFPACK) && defined(UMFPACK_ORDERING_NONE)
    }else{
        MSG_INFO("Build banded solver based on UMFPACK.\n");
        struct banded_umfpack * sol;
        mess_try_alloc(sol,struct banded_umfpack*, sizeof(struct banded_umfpack));
        sol->n = matrix->rows;
        sol->perm = perm;
        sol->invperm = invperm;
        ret = mess_direct_init(&(sol->umfpack_solver));                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
        double Control[UMFPACK_CONTROL];
        if ( MESS_IS_COMPLEX(matrix)) {
            umfpack_zl_defaults(Control);
        } else {
            umfpack_dl_defaults(Control);
        }
        Control[UMFPACK_ORDERING]=UMFPACK_ORDERING_NONE;
        ret = mess_direct_create_umfpack_control(temp,sol->umfpack_solver,Control);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_create_umfpack_control);

        /*-----------------------------------------------------------------------------
         *  set pointers
         *-----------------------------------------------------------------------------*/
        solver->data = (void *) sol;
        solver->clear = banded_umfpack_clear;
        solver->solve = banded_umfpack_solve;
        solver->solvet = banded_umfpack_solvet;
        solver->solveh = banded_umfpack_solveh;
        solver->solvem = banded_umfpack_solvem;
        solver->solvemt = banded_umfpack_solvemt;
        solver->solvemh = banded_umfpack_solvemh;
        solver->inverse = NULL;
        solver->det  = NULL;
        solver->detc  = NULL;
        solver->getL = NULL;
        solver->getU = NULL;
        solver->getpermp = NULL;
        solver->getpermq = NULL;
        solver->rows = matrix->rows;
        solver->cols = matrix->cols;
        solver->data_type = temp->data_type;
        SET_SOLVERNAME(solver->name, __func__);
    }
#endif

    /*-----------------------------------------------------------------------------
     *  clear memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&temp);

    return ret;
}



