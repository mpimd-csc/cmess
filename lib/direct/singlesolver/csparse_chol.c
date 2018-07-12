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
 * @file lib/direct/singlesolver/csparse_chol.c
 * @brief Interface for Tim Davis @cxsparse  (Cholesky decomposition).
 * @author @koehlerm
 *
 *
 *
 * This file provide an interface for @cxsparse. You need a working
 * @suitesparse  to use it.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#ifdef MESS_HAVE_CSPARSE
#include <cs.h>
#include "mess/interface_csparse.h"

/** Internal structure to save the factorization. */
struct csparse_solver {
    cs_dln *LU;     /**< numerical part */
    cs_dls *SLU;    /**< symbolical part */
    mess_int_t dim;     /**< dimension  */
    mess_direct_levelset Llevels;
};

/** Internal structure to save the factorization. */
struct csparse_solver_complex {
    cs_cln *LU;     /**< numerical part */
    cs_cls *SLU;    /**< symbolical part */
    mess_int_t dim;     /**< dimension  */
};


/* |+* type of reordering  */
#define CSPARSE_ORDER 2

static mess_int_t *perm_inv (mess_int_t const *p,mess_int_t n)
{
    MSG_FNAME(__func__);
    mess_int_t k, *pinv ;
    if (!p) return (NULL) ;
    mess_try_alloc2(pinv, mess_int_t*, n*sizeof (mess_int_t)) ;
    if (!pinv) return (NULL) ;
    for (k = 0 ; k < n ; k++) pinv [p [k]] = k ;
    return (pinv) ;
}

static int cs_cl_lttsolve (const cs_cl  *L, mess_double_cpx_t *x)
{
    long p, j, n, *Lp, *Li ;
    mess_double_cpx_t *Lx ;
    if (!CS_CSC (L) || !x) return (0) ;                     /* check inputs */
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (j = n-1 ; j >= 0 ; j--)
    {
        for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
        {
            x [j] -= (Lx [p]) * x [Li [p]] ;
        }
        x [j] /= (Lx [Lp [j]]) ;
    }
    return (1) ;
}

static int cs_dlc_lsolve (const cs_dl *L, mess_double_cpx_t *x)
{
    long p, j, n, *Lp, *Li ;
    double *Lx ;
    if (!CS_CSC (L) || !x) return (0) ;                     /* check inputs */
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (j = 0 ; j < n ; j++)
    {
        x [j] /= Lx [Lp [j]] ;
        for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
        {
            x [Li [p]] -= Lx [p] * x [j] ;
        }
    }
    return (1) ;
}


static int cs_dlc_ltsolve (const cs_dl  *L, mess_double_cpx_t *x)
{
    long  p, j, n, *Lp, *Li ;
    double *Lx ;
    if (!CS_CSC (L) || !x) return (0) ;                     /* check inputs */
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (j = n-1 ; j >= 0 ; j--)
    {
        for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
        {
            x [j] -= (Lx [p]) * x [Li [p]] ;
        }
        x [j] /= (Lx [Lp [j]]) ;
    }
    return (1) ;
}

/**
 * @brief Get the L (the first factor) of a solver.
 * @param[in] data  input solver data
 * @param[in,out] L     output L factor
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_getL function gets the L factor from the cs_dln object. Output is in CSR
 * format.
 *
 */
static int csparse_getL(void *data, mess_matrix L ){
    MSG_FNAME(__func__);
    struct csparse_solver *sol = (struct csparse_solver *) data;
    mess_matrix tmp =NULL ;
    MSG_INFO("get L\n");
    int ret;
    mess_check_nullpointer(data);
    mess_check_nullpointer(L);
    MESS_MATRIX_RESET(L);

    ret = mess_matrix_init(&tmp);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_from_csparse_dl(sol->LU->L,  tmp);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_from_csparse);
    ret = mess_matrix_convert ( tmp, L, MESS_CSR);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    mess_matrix_clear(&tmp);

    return 0;
}

/**
 * @brief Get the L (the first factor) of a solver (complex).
 * @param[in] data  input solver data
 * @param[in,out] L     output L factor
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_getL_complex function gets the L factor from the cs_cln object. Output is in CSR
 * format.
 *
 */
static int csparse_getL_complex(void *data, mess_matrix L ){
    MSG_FNAME(__func__);
    struct csparse_solver_complex *sol = (struct csparse_solver_complex *) data;
    mess_matrix tmp =NULL ;
    MSG_INFO("get L\n");
    int ret;
    mess_check_nullpointer(data);
    mess_check_nullpointer(L);
    MESS_MATRIX_RESET(L);

    ret = mess_matrix_init(&tmp);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_from_csparse_cl(sol->LU->L,  tmp);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_from_csparse);
    ret = mess_matrix_convert ( tmp, L, MESS_CSR);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    mess_matrix_clear(&tmp);

    return 0;
}


/**
 * @brief Get the U (the second factor) of a solver.
 * @param[in] data  input solver data
 * @param[in,out] U     output U factor
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_getU function gets the U factor from the cs_dln object. Output is in CSR
 * format.
 *
 */
static int csparse_getU(void *data, mess_matrix U ){
    MSG_FNAME(__func__);
    struct csparse_solver *sol = (struct csparse_solver *) data;
    mess_matrix tmp = NULL, tmp2 = NULL;
    int ret;
    mess_check_nullpointer(data);
    mess_check_nullpointer(U);
    MESS_MATRIX_RESET(U);
    MSG_INFO("get U\n");

    ret = mess_matrix_init(&tmp);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&tmp2);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_from_csparse_dl(sol->LU->L,  tmp);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_from_csparse);
    ret = mess_matrix_convert ( tmp, tmp2, MESS_CSR);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    ret = mess_matrix_ctranspose(tmp2, U);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
    mess_matrix_clear(&tmp);
    mess_matrix_clear(&tmp2);
    return 0;
}

/**
 * @brief Get the U (the second factor) of a solver (complex).
 * @param[in] data  input solver data
 * @param[in,out] U     output U factor
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_getU_complex function gets the U factor from the cs_cln object. Output is in CSR
 * format.
 *
 */
static int csparse_getU_complex(void *data, mess_matrix U ){
    MSG_FNAME(__func__);
    struct csparse_solver_complex *sol = (struct csparse_solver_complex *) data;
    mess_matrix tmp = NULL, tmp2=NULL;
    int ret;
    mess_check_nullpointer(data);
    mess_check_nullpointer(U);
    MESS_MATRIX_RESET(U);
    MSG_INFO("get U\n");

    ret = mess_matrix_init(&tmp);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&tmp2);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_from_csparse_cl(sol->LU->U,  tmp);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_from_csparse);
    ret = mess_matrix_ctranspose(tmp, tmp2);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
    ret = mess_matrix_convert ( tmp2, U, MESS_CSR);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    mess_matrix_clear(&tmp);
    mess_matrix_clear(&tmp2);
    return 0;
}

/**
 * @brief Get the row permutation (real).
 * @param[in] data  input solver data
 * @param[in,out] p output permutation
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_getpermp function gets the row permutation of the decomposition.
 *
 */
static int csparse_getpermp(void *data, mess_int_t *p){
    MSG_FNAME(__func__);
    struct csparse_solver *sol = (struct csparse_solver *) data;
    mess_int_t i=0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(p);

    mess_int_t *tmp;
    mess_int_t *t2;

    mess_try_alloc(tmp, mess_int_t*, sizeof(mess_int_t) * sol->dim);

    for (i = 0 ; i < sol-> dim; i++) {
        tmp [i ] = (sol->LU->pinv == NULL)?i:(sol->LU->pinv[i]);
    }
    t2 = perm_inv(tmp, sol->dim);

    for ( i=0;i<sol->dim ; i++ ) {
        // MSG_PRINT(" [ " MESS_PRINTF_INT " ] = " MESS_PRINTF_INT "\n", i,tmp[i]);
        p[i] = t2[i];
    }
    mess_free(tmp);
    mess_free(t2);

    return 0;
}

/**
 * @brief Get the row permutation (complex).
 * @param[in] data  input solver data
 * @param[in,out] p output permutation
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_getpermp_complex function gets the row permutation of the decomposition.
 *
 */
static int csparse_getpermp_complex(void *data, mess_int_t *p){
    MSG_FNAME(__func__);
    struct csparse_solver_complex *sol = (struct csparse_solver_complex *) data;
    mess_int_t i=0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(p);

    mess_int_t *tmp;
    mess_int_t *t2;

    mess_try_alloc(tmp, mess_int_t*, sizeof(mess_int_t) * sol->dim);

    for (i = 0 ; i < sol-> dim; i++) {
        tmp [i ] = (sol->LU->pinv == NULL)?i:(sol->LU->pinv[i]);
    }
    t2 = perm_inv(tmp, sol->dim);

    for ( i=0;i<sol->dim ; i++ ) {
        // MSG_PRINT(" [ " MESS_PRINTF_INT " ] = " MESS_PRINTF_INT "\n", i,tmp[i]);
        p[i] = t2[i];
    }
    mess_free(tmp);
    mess_free(t2);

    return 0;
}


/* x(p) = b, for dense vectors x and b; p=null denotes identity */
static int cs_dl_ipvecsplit (const long *p, const mess_double_cpx_t *b, double *x1, double *x2, mess_int_t n)
{
    int k ;
    if (!x1 || !x2 || !b) return (0) ;                              /* check inputs */
    for (k = 0 ; k < n ; k++) {
        x1 [p ? p [k] : k] = creal(b [k]) ;
        x2 [p ? p [k] : k] = cimag(b [k]) ;
    }
    return (1) ;
}

/* x(p) = b, for dense vectors x and b; p=null denotes identity */
static int cs_dl_pveccombine (const long *p,double *x1, double *x2, mess_double_cpx_t *b, mess_int_t n)
{
    int k ;
    if (!x1 || !x2 || !b) return (0) ;                              /* check inputs */
    for (k = 0 ; k < n ; k++) {
        b[ k ] = x1[p ? p[k]:k]+x2[p ? p[k]:k]*I;
    }
    return (1) ;
}




/**
 * @brief Solve \f$ Ax=b \f$
 * @param[in] data  input solver data
 * @param[in] b  input right hand side
 * @param[in,out] x solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_solve function solves the system \f$ Ax=b \f$ with \f$ x \f$ and \f$ b \f$ vectors.
 *
 */
static int csparse_solve(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    struct csparse_solver * sol = (struct csparse_solver *)data;
    mess_int_t n = 0;
    double *tmp = NULL, *tmp2=NULL;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);

    n = sol->dim;
    mess_vector_resize(x,n);
    if ( MESS_IS_REAL(b)){
        mess_vector_toreal_nowarn(x);
        mess_try_alloc(tmp, double *, sizeof(double) * n);
        cs_dl_ipvec (sol->SLU->pinv, b->values, tmp, n) ;       /* x = b(p) */
        cs_dl_lsolve (sol->LU->L, tmp) ;               /* x = L\x */
        cs_dl_ltsolve (sol->LU->L, tmp) ;               /* x = U\x */
        cs_dl_pvec (sol->SLU->pinv, tmp, x->values, n) ;          /* b(q) = x */
        mess_free(tmp);
    } else {
        mess_vector_tocomplex(x);
        mess_try_alloc(tmp, double *, sizeof(double) * n);
        mess_try_alloc(tmp2, double *, sizeof(double) * n);
        cs_dl_ipvecsplit (sol->SLU->pinv, b->values_cpx, tmp,tmp2, n) ;       /* x = b(p) */
        cs_dl_lsolve (sol->LU->L, tmp) ;               /* x = L\x */
        cs_dl_ltsolve (sol->LU->L, tmp) ;               /* x = U\x */
        cs_dl_lsolve (sol->LU->L, tmp2) ;        /* x = L\x */
        cs_dl_ltsolve (sol->LU->L, tmp2) ;               /* x = U\x */

        cs_dl_pveccombine (sol->SLU->pinv, tmp,tmp2, x->values_cpx, n) ;          /* b(q) = x */
        mess_free(tmp);
        mess_free(tmp2);

    }
    return 0;
}



/**
 * @brief Solve \f$ Ax=b \f$ (complex)
 * @param[in] data  input solver data
 * @param[in] b  input right hand side
 * @param[in,out] x solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_solve_complex function solves the system \f$ Ax=b \f$ with \f$ x \f$ and \f$ b \f$ vectors.
 *
 */
static int csparse_solve_complex(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    struct csparse_solver_complex * sol = (struct csparse_solver_complex *)data;
    mess_int_t n = 0;
    mess_double_cpx_t *tmp = 0;
    int ret = 0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    ret = mess_vector_tocomplex(b);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
    ret = mess_vector_tocomplex(x);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);

    n = sol->dim;
    mess_try_alloc(tmp, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * n);

    cs_cl_ipvec (sol->SLU->pinv, b->values_cpx, tmp, n) ;       /* x = b(p) */
    cs_cl_lsolve (sol->LU->L, tmp) ;               /* x = L\x */
    cs_cl_usolve (sol->LU->L, tmp) ;               /* x = U\x */
    cs_cl_pvec (sol->SLU->q, tmp, x->values_cpx, n) ;          /* b(q) = x */
    mess_free(tmp);
    return 0;
}



/**
 * @brief Solve \f$ A^Tx=b \f$
 * @param[in] data  input solver data
 * @param[in] b  input right hand side
 * @param[in,out] x solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_solvet function solves the system \f$ A^Tx=b \f$ with \f$ x \f$ and \f$ b \f$ vectors.
 *
 */
static int csparse_solvet(void*data, mess_vector b, mess_vector x) {
    return csparse_solve(data,b,x);
}

/**
 * @brief Solve \f$ A^Tx=b \f$ (complex)
 * @param[in] data  input solver data
 * @param[in] b  input right hand side
 * @param[in,out] x solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_solvet_complex function solves the system \f$ A^Tx=b \f$ with \f$ x \f$ and \f$ b \f$ vectors.
 *
 */
static int csparse_solvet_complex(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    struct csparse_solver_complex * sol = (struct csparse_solver_complex *)data;
    mess_int_t  n = 0;
    int ret = 0;
    mess_double_cpx_t *tmp = 0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    ret = mess_vector_tocomplex(b);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
    ret = mess_vector_tocomplex(x);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);

    n = sol->dim;
    mess_try_alloc(tmp, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * n);

    cs_cl_ipvec (sol->SLU->pinv, b->values_cpx, tmp, n) ;       /* x = b(Q') */
    cs_cl_lttsolve (sol->LU->L, tmp) ;               /* x = L\x */
    cs_cl_lsolve (sol->LU->L, tmp) ;               /* x = U\x */
    cs_cl_pvec (sol->SLU->pinv, tmp, x->values_cpx, n) ;          /* b(p') = x */
    mess_free(tmp);
    return 0;
}

/**
 * @brief Solve \f$ A^H x=b \f$ (complex)
 * @param[in] data  input solver data
 * @param[in] b  input right hand side
 * @param[in,out] x solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_solve function solves the system \f$ A^H x=b \f$ with \f$ x \f$ and \f$ b \f$ vectors.
 *
 */
static int csparse_solveh_complex(void*data, mess_vector b, mess_vector x) {
    return csparse_solve_complex(data, b, x);
}



/**
 * @brief Solve \f$ AX=B \f$ (matrix version)
 * @param[in] data  input solver data
 * @param[in] b  input right hand sides
 * @param[in,out] x solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_solvem function solves the system \f$ AX=B \f$ with \f$ X \f$  and \f$ B \f$ matrices.
 *
 */
static int csparse_solvem(void*data, mess_matrix b, mess_matrix x) {
    MSG_FNAME(__func__);
    struct csparse_solver * sol = (struct csparse_solver *)data;
    mess_int_t n = 0;
    mess_int_t i=0;
    double *tmp = 0;
    int conv = -1;
    mess_matrix work;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);

    n = sol->dim;
    MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
    int ret = mess_matrix_alloc(x, b->rows, b->cols, b->rows * b->cols, MESS_DENSE, b->data_type);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    if ( MESS_IS_REAL(b)){
#ifdef _OPENMP
#pragma omp parallel for private (i, tmp) default(shared)
#endif
        for ( i = 0; i<b->cols ; i++ ) {
            mess_try_alloc2( tmp,  double *,sizeof(double) * n);
            cs_dl_ipvec (sol->SLU->pinv, work->values+(i*work->ld), tmp, n) ;       /* x = b(p) */
            cs_dl_lsolve (sol->LU->L, tmp) ;               /* x = L\x */
            cs_dl_ltsolve (sol->LU->L, tmp) ;               /* x = U\x */
            cs_dl_pvec (sol->SLU->pinv, tmp, x->values+(i*x->ld), n) ;          /* b(q) = x */
            mess_free(tmp);
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for private (i, tmp) default(shared)
#endif
        for ( i = 0; i<b->cols ; i++ ) {
            mess_double_cpx_t *tmpc ;
            mess_try_alloc2(tmpc, mess_double_cpx_t *,sizeof(mess_double_cpx_t) * n);
            cs_cl_ipvec (sol->SLU->pinv, work->values_cpx+(i*work->ld), tmpc, n) ;       /* x = b(p) */
            cs_dlc_lsolve (sol->LU->L, tmpc) ;               /* x = L\x */
            cs_dlc_ltsolve (sol->LU->L, tmpc) ;               /* x = U\x */
            cs_cl_pvec (sol->SLU->pinv, tmpc, x->values_cpx+(i*x->ld), n) ;          /* b(q) = x */
            mess_free(tmpc);
        }

    }


    if ( conv == 0) {
        mess_matrix_clear(&work);
    }
    return 0;
}

/**
 * @brief Solve AX=B (complex matrix version).
 * @param[in] data  input solver data
 * @param[in] b  input right hand sides
 * @param[in,out] x solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_solvem_complex function solves the system \f$ AX=B \f$ with \f$ X \f$ and \f$ B \f$ matrices.
 *
 */
static int csparse_solvem_complex(void*data, mess_matrix b, mess_matrix x) {
    MSG_FNAME(__func__);
    struct csparse_solver_complex * sol = (struct csparse_solver_complex *)data;
    mess_int_t n = 0;
    int ret = 0;
    mess_int_t i=0;
    mess_double_cpx_t *tmp = 0;
    int conv = -1;
    mess_matrix work;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    ret = mess_matrix_tocomplex(b);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_tocomplex);

    n = sol->dim;
    MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
    ret = mess_matrix_alloc(x, b->rows, b->cols, b->rows * b->cols, MESS_DENSE, MESS_COMPLEX);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
#ifdef _OPENMP
#pragma omp parallel for private (i, tmp) default(shared)
#endif
    for ( i = 0; i<b->cols ; i++ ) {
        mess_try_alloc2(tmp , mess_double_cpx_t *,sizeof(mess_double_cpx_t) * n);
        cs_cl_ipvec (sol->SLU->pinv, work->values_cpx+(i*work->ld), tmp, n) ;       /* x = b(p) */
        cs_cl_lsolve (sol->LU->L, tmp) ;               /* x = L\x */
        cs_cl_usolve (sol->LU->L, tmp) ;               /* x = U\x */
        cs_cl_pvec (sol->SLU->pinv, tmp, x->values_cpx+(i*x->ld), n) ;          /* b(q) = x */
        mess_free(tmp);
    }

    if ( conv == 0) {
        mess_matrix_clear(&work);
    }
    return 0;
}

/**
 * @brief Solve \f$ A^T X=B \f$ (matrix version)
 * @param[in] data  input solver data
 * @param[in] b  input right hand sides
 * @param[in,out] x solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_solvemt function solves the system \f$ A^TX=B \f$ with \f$ X \f$ and \f$ B \f$ matrices.
 *
 */
static int csparse_solvemt(void*data, mess_matrix b, mess_matrix x) {
    return csparse_solvem(data,b,x);
}

/**
 * @brief Solve \f$ A^TX=B \f$ (complex matrix version)
 * @param[in] data  input solver data
 * @param[in] b  input right hand sides
 * @param[in,out] x solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_solvemt_complex function solves the system \f$ A^T X=B \f$ with \f$ X \f$ and \f$ B \f$ matrices.
 *
 */
static int csparse_solvemt_complex(void*data, mess_matrix b, mess_matrix x) {
    MSG_FNAME(__func__);
    struct csparse_solver_complex * sol = (struct csparse_solver_complex *)data;
    mess_int_t n = 0;
    mess_int_t i=0;
    int ret = 0;
    mess_double_cpx_t *tmp = NULL;
    int conv = -1;
    mess_matrix work;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    ret  = mess_matrix_tocomplex(b);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);

    n = sol->dim;
    MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
    ret  = mess_matrix_alloc(x, b->rows, b->cols, b->rows * b->cols, MESS_DENSE, MESS_COMPLEX);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
#ifdef _OPENMP
#pragma omp parallel for private (i, tmp) default(shared)
#endif
    for ( i = 0; i<b->cols ; i++ ) {
        mess_try_alloc2(tmp, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * n);
        cs_cl_ipvec (sol->SLU->pinv, work->values_cpx+(i*work->ld), tmp, n) ;       /* x = b(p) */
        cs_cl_lttsolve (sol->LU->L, tmp) ;               /* x = L\x */
        cs_cl_lsolve (sol->LU->L, tmp) ;               /* x = U\x */
        cs_cl_pvec (sol->LU->pinv, tmp, x->values_cpx+(i*x->ld), n) ;          /* b(q) = x */
        mess_free(tmp);
    }

    if ( conv == 0) {
        mess_matrix_clear(&work);
    }
    return 0;

}


/**
 * @brief Solve \f$ A^HX=B \f$ (complex matrix version)
 * @param[in] data  input solver data
 * @param[in] b  input right hand sides
 * @param[in,out] x solution
 *
 * The @ref csparse_solvemt_complex function solves the system \f$ A^T X=B \f$ with \f$ X \f$ and \f$ B \f$ matrices.
 *
 */
static int csparse_solvemh_complex(void*data, mess_matrix b, mess_matrix x) {
    return csparse_solvem_complex(data, b, x ) ;
}

/**
 * @brief Clear a @csparse solver.
 * @param[in,out] solver solver data
 * @return always zero
 *
 * The @ref csparse_clear function is the clean up function.
 *
 */
static int csparse_clear(void *solver){
    struct csparse_solver * sol = (struct csparse_solver*) solver;
    if ( sol != NULL) {
        cs_dl_sfree(sol->SLU);
        cs_dl_nfree(sol->LU);
        mess_free( sol );
    }
    return 0;
}

/**
 * @brief Clear a @csparse solver (complex).
 * @param[in,out] solver solver data
 * @return always zero
 *
 * The @ref csparse_clear_complex function is the clean up function.
 *
 */
static int csparse_clear_complex(void *solver){
    struct csparse_solver_complex * sol = (struct csparse_solver_complex*) solver;
    if ( sol != NULL) {
        cs_cl_sfree(sol->SLU);
        cs_cl_nfree(sol->LU);
        mess_free( sol );
    }
    return 0;
}


/**
 * @brief Compute the inverse.
 * @param[in] data   input pointer to the internal data structure
 * @param[in,out] inv   output inverse
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_inverse function computes \f$ A^{-1} \f$.
 *
 */
static int csparse_inverse ( void *data, mess_matrix inv )
{
    MSG_FNAME(__func__);
    struct csparse_solver * sol = (struct csparse_solver*) data;
    mess_matrix eye;
    int ret = 0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(inv);
    MESS_MATRIX_RESET(inv);

    ret = mess_matrix_init(&eye);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_eye(eye, sol->dim, sol->dim, MESS_DENSE);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
    ret = csparse_solvem(data, eye, inv);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), csparse_solvem);
    mess_matrix_clear(&eye);
    return 0;
}       /* -----  end of function csparse_inverse  ----- */

/**
 * @brief Compute the inverse (complex).
 * @param[in] data   input pointer to the internal data structure
 * @param[in,out] inv   output inverse
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref csparse_inverse_complex function computes \f$ A^{-1} \f$
 *
 */
static int csparse_inverse_complex ( void *data, mess_matrix inv )
{
    MSG_FNAME(__func__);
    struct csparse_solver * sol = (struct csparse_solver*) data;
    mess_matrix eye;
    int ret = 0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(inv);
    MESS_MATRIX_RESET(inv);

    ret = mess_matrix_init(&eye);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_eyec(eye, sol->dim, sol->dim, MESS_DENSE);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
    ret = csparse_solvem_complex(data, eye, inv);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), csparse_solvem_complex);
    mess_matrix_clear(&eye);
    return 0;
}       /* -----  end of function csparse_inverse  ----- */



/**
 * @brief Generate a CXSparse Cholesky solver.
 * @param[in] matrix  input system matrix
 * @param[in,out] solver output solver
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_create_csparse_cholesky function generates a LU-factorization based solver. \n
 * It uses @csparse[1] to reorder and factorize the system matrix. \n
 * The reordering can be modified by using following runtime configurations:
 * <ul>
 * <li>  mess.csparse.order, \f$ 0 \f$=no reordering, \f$ 1,2 \f$= different AMD reordering, default \f$ 2 \f$
 * <li>  mess.csparse.tol, tolerance in pivoting, default \f$ 1.0 \f$
 * </ul>
 * @cite Dav06.
 *
 */
int mess_direct_create_csparse_cholesky(mess_matrix matrix, mess_direct solver)
{
    MSG_FNAME(__func__);
    int ret = 0;
    int order = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_square(matrix);
    mess_check_real_or_complex(matrix);
    if ( MESS_IS_COMPLEX(matrix)){
        struct csparse_solver_complex * data;
        cs_cl *A;

        mess_try_alloc(data, struct csparse_solver_complex*, sizeof(struct csparse_solver_complex));
        order = CSPARSE_ORDER;

        /*-----------------------------------------------------------------------------
         *  factorize matrix
         *-----------------------------------------------------------------------------*/
        ret  =mess_matrix_to_csparse_cl( matrix , &A);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_to_csparse_complex);
        data->SLU = cs_cl_schol(order, A);
        data->LU  = cs_cl_chol(A, data->SLU);
        if ( data->SLU==NULL || data->LU == NULL){
            MSG_ERROR("csparse returned with error\n");
            return MESS_ERROR_MISC;
        }

        /*-----------------------------------------------------------------------------
         *  construct the solver structure
         *-----------------------------------------------------------------------------*/

        data->dim = matrix->rows;
        solver->rows = matrix->rows;
        solver->cols = matrix->cols;
        solver->data = (void *)data;
        solver->solve = csparse_solve_complex;
        solver->solvet = csparse_solvet_complex;
        solver->solveh = csparse_solveh_complex;
        solver->solvem = csparse_solvem_complex;
        solver->solvemt = csparse_solvemt_complex;
        solver->solvemh = csparse_solvemh_complex;
        solver->clear = csparse_clear_complex;
        solver->getL = csparse_getL_complex;
        solver->getU = csparse_getU_complex;
        solver->getpermp = csparse_getpermp_complex;
        solver->inverse = csparse_inverse_complex;
        solver->data_type = MESS_COMPLEX;

        cs_cl_spfree(A);
        //cs_cl_free(A);

    } else {
        struct csparse_solver * data;
        cs_dl *A;

        mess_try_alloc(data, struct csparse_solver*, sizeof(struct csparse_solver));

        order = CSPARSE_ORDER;

        /*-----------------------------------------------------------------------------
         *  factorize matrix
         *-----------------------------------------------------------------------------*/
        ret  =mess_matrix_to_csparse_dl ( matrix , &A);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_to_csparse);

        data->SLU = cs_dl_schol(order, A);
        data->LU  = cs_dl_chol(A, data->SLU);
        if ( data->SLU==NULL || data->LU == NULL){
            MSG_ERROR("csparse returned with error\n");
            return MESS_ERROR_MISC;
        }
        /*-----------------------------------------------------------------------------
         *  construct the solver structure
         *-----------------------------------------------------------------------------*/

        data->dim = matrix->rows;
        solver->rows = matrix->rows;
        solver->cols = matrix->cols;
        solver->data = (void *)data;
        solver->solve = csparse_solve;
        solver->solvet = csparse_solvet;
        solver->solveh = csparse_solvet;
        solver->solvem = csparse_solvem;
        solver->solvemt = csparse_solvemt;
        solver->solvemh = csparse_solvemt;
        solver->clear = csparse_clear;
        solver->getL = csparse_getL;
        solver->getU = csparse_getU;
        solver->getpermp = csparse_getpermp;
        solver->inverse = csparse_inverse;
        solver->data_type = MESS_REAL;

        cs_dl_spfree(A);
        //  cs_dl_free(A);
    }
    SET_SOLVERNAME(solver->name, __func__);
    return 0;
}

#else

#endif


