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
 * @file lib/direct/multisolver/multilu_adavanced.c
 * @brief Multisolver implementation with mess lu.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>
#include <math.h>

#ifdef MESS_HAVE_UMFPACK
#include <umfpack.h>
#endif

#define WORK_BACKGROUND 0

#include <pthread.h>
#include <unistd.h>
#include <sched.h>

#if  _OPENMP
#include <omp.h>
#endif

#ifdef MESS_DEBUG
#define DTINIT(X) unsigned short X##_dt = (X)->data_type;
#define DTCHECK(X) assert(X##_dt == (X)->data_type);
#else
#define DTINIT(X)
#define DTCHECK(X)
#endif

#define IS_REAL(X)  (fabs(cimag((X))) < 1.1*eps)
#define IS_COMPLEX(X) (fabs(cimag((X))) >=  1.1*eps)


/**
 * @brief Internal structure to store a multilu.
 *
 * */
struct multilu {
        // unsigned short object;   /*< object id */
        double **lvalues;       /*< values for the L factors */
        double **uvalues;       /*< values for the U factors */
        mess_double_cpx_t **lvalues_cpx; /*< values for the L factors (complex) */
        mess_double_cpx_t **uvalues_cpx; /*< values for the L factors (complex) */
        unsigned short *datatypes; /*< data type for the varios solvers */
        mess_int_t nlu;     /*< number of LU decompositions */
        mess_int_t *lind;       /*< index pointer  for L */
        mess_int_t *lp;     /*< row/column start for L */
        mess_int_t *uind;       /*< index pointer  for U */
        mess_int_t *up;     /*< row/column start for U */
        mess_int_t lnnz;    /*< number of non zeros elements in L*/
        mess_int_t unnz;    /*< number of non zeros elements in U*/
        mess_int_t *p;      /*< row permutation */
        mess_int_t *pinv;       /*< inverse row permutation */
        mess_int_t *q;      /*< column permutation */
        mess_int_t *qinv;       /*< row permutation */
        mess_int_t rows;        /*< number of rows */
        mess_int_t cols;        /*< number of columns */
        unsigned short lstore_type;
        unsigned short ustore_type;
        unsigned short  store_type;  /*< in which way are the pattern / values stored */

        unsigned short background ;
        pthread_cond_t *sig_ready;
        pthread_mutex_t *mutex_ready;
        mess_threadpool solverpool;
        mess_int_t *ready;
        mess_int_t *workrowptr;
        mess_int_t *workcolptr;
};


/*-----------------------------------------------------------------------------
 *  include other functions.
 *-----------------------------------------------------------------------------*/


/**
 * @brief Get the size of the data element.
 * @param[in] data  input pointer to the data
 * @return size of data, otherwise zero
 *
 * The @ref multilu_memsize function gets the size of data.
 *
 */
static size_t multilu_memsize ( void * data )
{
    struct multilu * mlu = (struct multilu *) data;
    size_t s = 0;
    mess_int_t i;
    if ( mlu == NULL) return 0;

    s = (sizeof (mess_int_t)) * mlu->lnnz;
    s+= (sizeof (mess_int_t)) * mlu->unnz;
    s+= (sizeof(mess_int_t)) * 2 * (mlu->rows+1);
    for (i = 0 ; i < mlu->nlu; i++){
        if ( mlu->datatypes[i] == MESS_REAL) {
            s+= (sizeof(double) * mlu->lnnz);
            s+= (sizeof(double) * mlu->unnz);
        } else {
            s+= (sizeof(double) * mlu->lnnz);
            s+= (sizeof(double) * mlu->unnz);
        }
    }
    if ( mlu->p) s+=sizeof(mess_int_t)*mlu->rows;
    if ( mlu->pinv) s+=sizeof(mess_int_t)*mlu->rows;
    if ( mlu->q) s+=sizeof(mess_int_t)*mlu->rows;
    if ( mlu->qinv) s+=sizeof(mess_int_t)*mlu->rows;
    return s;
}

/**
 * @brief Invert a permutation.
 * @param[in] p  input permutation
 * @param[in] n  input length of the permutation
 * @return inverse permutation
 *
 * The @ref pinv function returns the inverted permution of p.
 *
 */
static mess_int_t *pinv (mess_int_t const *p,mess_int_t n)
{
    MSG_FNAME(__func__);
    mess_int_t k, *pinv ;
    if (!p) return (NULL) ;
    mess_try_alloc2(pinv, mess_int_t *, n*sizeof (mess_int_t));
    if (!pinv) return (NULL) ;
    for (k = 0 ; k < n ; k++) pinv [p [k]] = k ;
    return (pinv) ;
}



/**
 * @brief Initialize a "single-pattern-multi-value" LU decompositon.
 * @param[out] mlu  the spmv to initialize
 * @return always zero
 *
 *
 * The @ref multilu_init function initializes a single-patter-multi-value
 * LU decompositon.
 *
 */
static int multilu_init(struct multilu *mlu){
    MSG_FNAME(__func__);

    mess_check_nullpointer(mlu);
    mlu->lvalues = mlu->uvalues = NULL;
    mlu->lvalues_cpx = mlu->uvalues_cpx = NULL;
    mlu->datatypes = NULL;
    mlu->nlu = 0;
    mlu->rows = 0;
    mlu->cols = 0;
    mlu->lind = NULL;
    mlu->uind = NULL;
    mlu->lp = NULL;
    mlu->up = NULL;
    mlu->lnnz = 0;
    mlu->unnz = 0;
    mlu->store_type = 0;
    mlu->p  = NULL;
    mlu->pinv = NULL;
    mlu->q = NULL;
    mlu->qinv = NULL;
    mlu->ready = NULL;
    mlu->sig_ready = NULL;
    mlu->mutex_ready = NULL;
    mlu->background = 0;
    mlu->workcolptr = NULL;
    mlu->workrowptr = NULL;
    return 0;
}


/**
 * @brief Clean up a multilu.
 * @param[in] data  input pointer to the internal data structure
 * @return always zero
 *
 * The @ref multilu_clear function is a internal function of the multilu-solver to
 * clean up when @ref mess_multidirect_clear is called.
 *
 */
static int multilu_clear(void  *data){
    MSG_FNAME(__func__);
    struct multilu *mlu = (struct multilu*) data;
    mess_int_t i;
    mess_int_t r = 0;
    mess_check_nullpointer(mlu);
    if ( mlu->background ) {
        pthread_mutex_lock ( mlu->mutex_ready ) ;
        r = 0 ;
        for (i = 0; i < mlu->nlu; i++) {
            if ( mlu->ready[i] ) r++;
        }
        while ( r != mlu->nlu) {
            pthread_cond_wait(mlu->sig_ready, mlu->mutex_ready);
            r = 0 ;
            for (i = 0; i < mlu->nlu; i++) {
                if ( mlu->ready[i] ) r++;
            }
        }
        pthread_mutex_unlock(mlu->mutex_ready);

        MSG_INFO("clear background solver \n");
        mess_threadpool_clear(&(mlu->solverpool));
        MSG_INFO("threadpool cleared.\n");

        pthread_cond_destroy(mlu->sig_ready);
        pthread_mutex_destroy(mlu->mutex_ready);
        mess_free(mlu->sig_ready);
        mess_free(mlu->mutex_ready);
        mess_free(mlu->ready);
        mess_free(mlu->workrowptr);
        mess_free(mlu->workcolptr);
     }

    if (mlu->nlu > 0 ) {
        for ( i = 0; i < mlu->nlu; i++){
            if (mlu->datatypes[i] == MESS_REAL && mlu->lvalues !=NULL){
                mess_free(mlu->lvalues[i]);
                mess_free(mlu->uvalues[i]);
            }
            if (mlu->datatypes[i] == MESS_COMPLEX && mlu->lvalues_cpx !=NULL){
                mess_free(mlu->lvalues_cpx[i]);
                mess_free(mlu->uvalues_cpx[i]);
            }
        }
        // if ( mlu->data_type == MESS_REAL) {
        mess_free(mlu->lvalues);
        mess_free(mlu->uvalues);
        // }
        // if ( mlu->data_type == MESS_COMPLEX) {
        mess_free(mlu->lvalues_cpx);
        mess_free(mlu->uvalues_cpx);
        // }

    }
    mess_free(mlu->lind);
    mess_free(mlu->uind);
    mess_free(mlu->lp);
    mess_free(mlu->up);
    mess_free(mlu->p);
    mess_free(mlu->pinv);
    mess_free(mlu->q);
    mess_free(mlu->qinv);
    mess_free(mlu->datatypes);

    mess_free(mlu);
    mlu = NULL;
    return 0;
}


/**
 * @brief Get the L factor of a decomposition.
 * @param[in] data   input pointer to the internal data structure
 * @param[in] ind    input index to get
 * @param[out] L    output L factor
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref multilu_getL function gets the L factor from a decomposition.
 *
 */
static int  multilu_getL ( void *data, mess_int_t ind, mess_matrix L )
{
    MSG_FNAME(__func__);
    int ret = 0 ;
    struct multilu *mlu = ( struct multilu *) data;

    mess_check_nullpointer(data);
    mess_check_nullpointer(L);

    if ( ind < 0 || ind>=mlu->nlu){
        MSG_ERROR("ind is out of range ( = " MESS_PRINTF_INT " )\n", ind);
        return MESS_ERROR_ARGUMENTS;
    }


    ret = mess_matrix_alloc (L, mlu->rows, mlu->cols, mlu->lnnz, MESS_CSR, mlu->datatypes[ind]);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    memcpy(L->rowptr, mlu->lp, sizeof(mess_int_t) * (L->rows +1) ) ;
    memcpy(L->colptr, mlu->lind, sizeof(mess_int_t) * (L->nnz) ) ;
    if (mlu->datatypes[ind]== MESS_REAL){
        memcpy(L->values, mlu->lvalues[ind], sizeof(double) * (L->nnz) ) ;
    } else {
        memcpy(L->values_cpx, mlu->lvalues_cpx[ind], sizeof(double) * (L->nnz) ) ;
    }

    return 0;
}       /* -----  end of function multilu_getL  ----- */

/**
 * @brief Get the U factor of a decomposition.
 * @param[in] data   input pointer to the internal data structure
 * @param[in] ind    input index to get
 * @param[out] U    output U factor
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref multilu_getU function gets the U factor from a decomposition.
 *
 */
static int  multilu_getU ( void *data, mess_int_t ind, mess_matrix U  )
{
    MSG_FNAME(__func__);
    int ret = 0 ;
    struct multilu *mlu = ( struct multilu *) data;

    mess_check_nullpointer(data);
    mess_check_nullpointer(U);

    if ( ind < 0 || ind>=mlu->nlu){
        MSG_ERROR("ind is out of range ( = " MESS_PRINTF_INT " )\n", ind);
        return MESS_ERROR_ARGUMENTS;
    }

    ret = mess_matrix_alloc (U, mlu->rows, mlu->cols, mlu->lnnz, MESS_CSR, mlu->datatypes[ind]);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    memcpy(U->rowptr, mlu->up, sizeof(mess_int_t) * (U->rows +1) ) ;
    memcpy(U->colptr, mlu->uind, sizeof(mess_int_t) * (U->nnz) ) ;
    if (mlu->datatypes[ind]== MESS_REAL){
        memcpy(U->values, mlu->uvalues[ind], sizeof(double) * (U->nnz) ) ;
    } else {
        memcpy(U->values_cpx, mlu->uvalues_cpx[ind], sizeof(double) * (U->nnz) ) ;
    }

    return 0;
}       /* -----  end of function multilu_getU  ----- */


/**
 * @brief Get the row permutation from the solver.
 * @param[in] data  input pointer to internal data structure
 * @param[in] ind    input index of the solver
 * @param[out] p    output row permutation
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref multilu_getp function gets the row permutation of a solver.
 * p must be preallocated.
 *
 */
static int  multilu_getp ( void * data, mess_int_t ind, mess_int_t * p  )
{
    MSG_FNAME(__func__);
    struct multilu *mlu = ( struct multilu *) data;
    mess_check_nullpointer(data);
    mess_check_nullpointer(p);

    memcpy (p, mlu->p, sizeof(mess_int_t) * mlu->rows);
    return 0;
}       /* -----  end of function multilu_getp  ----- */

/**
 * @brief Get the column permutation from the solver.
 * @param[in] data  input pointer to internal data structure
 * @param[in] ind    input index of the solver
 * @param[out] q    output column permutation
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref multilu_getq function gets the row permutation of a solver.
 * q must be preallocated.
 *
 */
static int  multilu_getq ( void * data, mess_int_t ind, mess_int_t * q  )
{
    MSG_FNAME(__func__);
    struct multilu *mlu = ( struct multilu *) data;
    mess_check_nullpointer(data);
    mess_check_nullpointer(q);

    memcpy (q, mlu->q, sizeof(mess_int_t) * mlu->cols);
    return 0;
}       /* -----  end of function multilu_getp  ----- */


/**
 * @brief Solve the system \f$ A(n)x=b \f$.
 * @param[in] data  input pointer to a multilu
 * @param[in] n  input index of the decomposition
 * @param[in] b  input right hand side
 * @param[out] x solution vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref multilu_solve function solves the system
 * \f[ A(n)x=b, \f]
 * where \f$ b \f$ and \f$ x \f$ are vectors.
 *
 */
static int multilu_solve(void *data, mess_int_t n, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct multilu *mlu = ( struct multilu *) data;
    mess_vector y;
    int ret = 0 ;

    mess_check_nullpointer(mlu);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    // mess_check_real(x);
    // mess_check_real(b);
    // mess_check_real(mlu);

    if ( n < 0 || n>=mlu->nlu){
        MSG_ERROR("n is out of range ( = " MESS_PRINTF_INT " )\n", n);
        return MESS_ERROR_ARGUMENTS;
    }
    if ( x->dim != mlu->cols){
        MSG_WARN("resize x from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", x->dim, mlu->cols);
        ret = mess_vector_resize(x, mlu->cols);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }

#ifdef _OPENMP
    if ( mlu->background ) {
        // MSG_INFO("waiting to get the solver ready\n");
        pthread_mutex_lock ( mlu->mutex_ready ) ;
        while ( mlu->ready[n] == 0) {
            pthread_cond_wait(mlu->sig_ready, mlu->mutex_ready);
        }
        pthread_mutex_unlock(mlu->mutex_ready);
    }
#endif



    /*-----------------------------------------------------------------------------
     *  real solver and real right hand sides
     *-----------------------------------------------------------------------------*/
    if ( mlu->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)){
        ret = mess_vector_init(&y);                                                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(y, x->dim, MESS_REAL);                                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_perm(b, mlu->p, y);                                                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_perm);
        ret = mess_solver_lsolve_kernelcsr_real(mlu->rows, mlu->lvalues[n], mlu->lind, mlu->lp, y->values);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_lsolve_kernelcsr_real);
        ret = mess_solver_usolve_kernelcsr_real(mlu->rows, mlu->uvalues[n], mlu->uind, mlu->up, y->values);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_usolve_kernelcsr_real);
        ret = mess_vector_iperm(y, mlu->q, x);                                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_iperm);
        mess_vector_clear(&y);
    }

    /*-----------------------------------------------------------------------------
     *  real solver and complex rhs
     *-----------------------------------------------------------------------------*/
    else if (mlu->datatypes[n]==MESS_REAL && MESS_IS_COMPLEX(b)) {
        ret = mess_vector_init(&y);                                                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(y, x->dim, MESS_COMPLEX);                                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_perm(b, mlu->p, y);                                                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_perm);
        ret = mess_solver_lsolve_kernelcsr_real_complex(mlu->rows, mlu->lvalues[n], mlu->lind, mlu->lp, y->values_cpx); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_lsolve_kernelcsr_real_complex);
        ret = mess_solver_usolve_kernelcsr_real_complex(mlu->rows, mlu->uvalues[n], mlu->uind, mlu->up, y->values_cpx); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_usolve_kernelcsr_real_complex);
        ret = mess_vector_tocomplex(x);                                                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_iperm(y, mlu->q, x);                                                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_iperm);
        mess_vector_clear(&y);
    }
    /*-----------------------------------------------------------------------------
     * complex solver
     *-----------------------------------------------------------------------------*/
    else {
        ret = mess_vector_init(&y);                                                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(y, x->dim, MESS_COMPLEX);                                                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_perm(b, mlu->p, y);                                                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_perm);
        ret = mess_solver_lsolve_kernelcsr_complex(mlu->rows, mlu->lvalues_cpx[n], mlu->lind, mlu->lp, y->values_cpx);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_lsolve_kernelcsr_real_complex);
        ret = mess_solver_usolve_kernelcsr_complex(mlu->rows, mlu->uvalues_cpx[n], mlu->uind, mlu->up, y->values_cpx);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_usolve_kernelcsr_real_complex);
        ret = mess_vector_tocomplex(x);                                                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_iperm(y, mlu->q, x);                                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_iperm);
        mess_vector_clear(&y);
    }
    return 0;
}

/**
 * @brief Solve the system \f$ A(n)X=B \f$.
 * @param[in] data  input pointer to a multilu
 * @param[in] n  input index of the decomposition
 * @param[in] b  input right hand side
 * @param[out] x solution matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref multilu_solvem function solves the system
 * \f[ A(n)X=B, \f]
 * where \f$ B \f$ and \f$ X \f$ are matrices.
 *
 */
static int multilu_solvem(void *data, mess_int_t n, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct multilu *mlu = ( struct multilu *) data;
    mess_vector y;
    mess_matrix work;
    int conv = -1;
    mess_int_t i,j;
    int ret=0;
    unsigned short outdt = MESS_REAL;


    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    DTINIT(b);

    if ( n < 0 || n>=mlu->nlu){
        MSG_ERROR("n is out of range ( = " MESS_PRINTF_INT " )\n", n)
        return MESS_ERROR_ARGUMENTS;
    }
    if ( mlu->rows != b->rows) {
        MSG_ERROR("b don't have the right number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }
#ifdef _OPENMP
    if ( mlu->background ) {
        // MSG_INFO("waiting to get the solver ready\n");
        pthread_mutex_lock ( mlu->mutex_ready ) ;
        while ( mlu->ready[n] == 0) {
            pthread_cond_wait(mlu->sig_ready, mlu->mutex_ready);
        }
        pthread_mutex_unlock(mlu->mutex_ready);
    }
#endif



    MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
    MESS_MATRIX_RESET(x);
    DTCHECK(b);

    if ( mlu->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)){
        outdt  = MESS_REAL;
    } else  if ( mlu->datatypes[n] == MESS_COMPLEX || MESS_IS_COMPLEX(b)){
        outdt = MESS_COMPLEX;
    }
    ret = mess_matrix_alloc(x, mlu->rows, work->cols,mlu->rows*work->cols, MESS_DENSE, outdt );
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    /*-----------------------------------------------------------------------------
     *  real solver and real right hand sides
     *-----------------------------------------------------------------------------*/
    if ( mlu->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)){
        ret = 0 ;
        #ifdef _OPENMP
        #pragma omp parallel for private (y, j, i)  default(shared) reduction(+:ret)
        #endif
        for ( i = 0; i < work->cols; i++) {
            ret += mess_vector_init(&y);
            ret += mess_vector_alloc(y, work->rows, MESS_REAL);
            // perm LU->p
            for (j = 0 ; j < work->rows ; j++) { y->values[j] = work->values[(mlu->p ? mlu->p [j] : j)+i*work->ld] ;}
            ret += mess_solver_lsolve_kernelcsr_real(mlu->rows, mlu->lvalues[n], mlu->lind, mlu->lp,y->values);
            ret += mess_solver_usolve_kernelcsr_real(mlu->rows, mlu->uvalues[n], mlu->uind, mlu->up, y->values);
            // iperm LU->q
            for (j = 0 ; j < work->rows; j++){ x->values[(mlu->q ? mlu->q [j] : j)+i*work->ld] = y->values[j] ;}
            mess_vector_clear(&y);
        }
    }

    /*-----------------------------------------------------------------------------
     *  real solver and complex rhs
     *-----------------------------------------------------------------------------*/
    else if (mlu->datatypes[n]==MESS_REAL && MESS_IS_COMPLEX(b)) {
        ret = 0 ;
        #ifdef _OPENMP
        #pragma omp parallel for private (y, j, i)  default(shared) reduction(+:ret)
        #endif
        for ( i = 0; i < work->cols; i++) {
            ret += mess_vector_init(&y);
            ret += mess_vector_alloc(y, work->rows, MESS_COMPLEX);
            // perm LU->p
            for (j = 0 ; j < work->rows ; j++) { y->values_cpx[j] = work->values_cpx[(mlu->p ? mlu->p [j] : j)+i*work->ld] ;}
            ret += mess_solver_lsolve_kernelcsr_real_complex(mlu->rows, mlu->lvalues[n], mlu->lind, mlu->lp,y->values_cpx);
            ret += mess_solver_usolve_kernelcsr_real_complex(mlu->rows, mlu->uvalues[n], mlu->uind, mlu->up, y->values_cpx);
            // iperm LU->q
            for (j = 0 ; j < work->rows; j++){ x->values_cpx[(mlu->q ? mlu->q [j] : j)+i*work->ld] = y->values_cpx[j] ;}
            mess_vector_clear(&y);
        }
    }
    /*-----------------------------------------------------------------------------
     * complex solver
     *-----------------------------------------------------------------------------*/
    else {
        ret = 0 ;
        #ifdef _OPENMP
        #pragma omp parallel for private (y, j, i)  default(shared) reduction(+:ret)
        #endif
        for ( i = 0; i < work->cols; i++) {
            ret += mess_vector_init(&y);
            ret += mess_vector_alloc(y, work->rows, MESS_COMPLEX);
            // perm LU->p
            if ( MESS_IS_REAL(work)) {
                for (j = 0 ; j < work->rows ; j++) { y->values_cpx[j] = work->values[(mlu->p ? mlu->p [j] : j)+i*work->ld] ;}
            } else {
                for (j = 0 ; j < work->rows ; j++) { y->values_cpx[j] = work->values_cpx[(mlu->p ? mlu->p [j] : j)+i*work->ld] ;}
            }
            ret += mess_solver_lsolve_kernelcsr_complex(mlu->rows, mlu->lvalues_cpx[n], mlu->lind, mlu->lp,y->values_cpx);
            ret += mess_solver_usolve_kernelcsr_complex(mlu->rows, mlu->uvalues_cpx[n], mlu->uind, mlu->up, y->values_cpx);
            // iperm LU->q
            for (j = 0 ; j < work->rows; j++){ x->values_cpx[(mlu->q ? mlu->q [j] : j)+i*work->ld] = y->values_cpx[j] ;}
            mess_vector_clear(&y);
        }
    }
    if ( conv == 0) {
        mess_matrix_clear(&work);
    }

    if ( ret!=0) {
        MSG_ERROR("an error occured while solving\n");
        return MESS_ERROR_GENERAL;
    }

    return 0;
}


/**
 * @brief Solve the system \f$ A(n)^T x=b \f$.
 * @param[in] data  input pointer to a multilu
 * @param[in] n  input index of the decomposition
 * @param[in] b  input right hand side
 * @param[out] x solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref multilu_solvet function solves the system
 * \f[ A(n)^T x=b, \f]
 * where \f$ b \f$ and \f$ x \f$ are vectors.
 *
 */
static int multilu_solvet(void *data, mess_int_t n, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct multilu *mlu = ( struct multilu *) data;
    mess_vector y;
    int ret = 0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( n < 0 || n>=mlu->nlu){
        MSG_ERROR("n is out of range ( = " MESS_PRINTF_INT " )\n", n)
        return MESS_ERROR_ARGUMENTS;
    }
    if ( x->dim != mlu->cols){
        MSG_WARN("size x from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", x->dim, mlu->cols);
        mess_vector_resize(x, mlu->cols);
    }
#ifdef _OPENMP
    if ( mlu->background ) {
        // MSG_INFO("waiting to get the solver ready\n");
        pthread_mutex_lock ( mlu->mutex_ready ) ;
        while ( mlu->ready[n] == 0) {
            pthread_cond_wait(mlu->sig_ready, mlu->mutex_ready);
        }
        pthread_mutex_unlock(mlu->mutex_ready);
    }
#endif


    /*-----------------------------------------------------------------------------
     *  real solver and real right hand sides
     *-----------------------------------------------------------------------------*/
    if ( mlu->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)){
        ret = mess_vector_init(&y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(y, x->dim, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_perm(b, mlu->q, y);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_perm);
        ret = mess_solver_utsolve_kernelcsr_real(mlu->rows, mlu->uvalues[n], mlu->uind, mlu->up, y->values);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_lsolve_kernelcsr_real);
        ret = mess_solver_ltsolve_kernelcsr_real(mlu->rows, mlu->lvalues[n], mlu->lind, mlu->lp, y->values);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_usolve_kernelcsr_real);
        ret = mess_vector_iperm(y, mlu->p, x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_iperm);
        mess_vector_clear(&y);
    }

    /*-----------------------------------------------------------------------------
     *  real solver and complex rhs
     *-----------------------------------------------------------------------------*/
    else if (mlu->datatypes[n]==MESS_REAL && MESS_IS_COMPLEX(b)) {
        ret = mess_vector_init(&y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(y, x->dim, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_perm(b, mlu->q, y);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_perm);
        ret = mess_solver_utsolve_kernelcsr_real_complex(mlu->rows, mlu->uvalues[n], mlu->uind, mlu->up, y->values_cpx);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_lsolve_kernelcsr_real_complex);
        ret = mess_solver_ltsolve_kernelcsr_real_complex(mlu->rows, mlu->lvalues[n], mlu->lind, mlu->lp, y->values_cpx);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_usolve_kernelcsr_real_complex);
        ret = mess_vector_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        ret = mess_vector_iperm(y, mlu->p, x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_iperm);
        mess_vector_clear(&y);
    }
    /*-----------------------------------------------------------------------------
     * complex solver
     *-----------------------------------------------------------------------------*/
    else {
        // MSG_WARN("complex transpose solve\n");
        // ret = mess_vector_tocomplex(b); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        ret = mess_vector_init(&y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(y, b->dim, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_perm(b, mlu->q, y);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_perm);
        ret = mess_solver_utsolve_kernelcsr_complex(mlu->rows, mlu->uvalues_cpx[n], mlu->uind, mlu->up, y->values_cpx);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_lsolve_kernelcsr_real_complex);
        ret = mess_solver_ltsolve_kernelcsr_complex(mlu->rows, mlu->lvalues_cpx[n], mlu->lind, mlu->lp, y->values_cpx);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_usolve_kernelcsr_real_complex);
        ret = mess_vector_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        ret = mess_vector_iperm(y, mlu->p, x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_iperm);
        mess_vector_clear(&y);
    }

    return 0;
}

/**
 * @brief Solve the system \f$ A(n)^H x=b \f$.
 * @param[in] data  input pointer to a multilu
 * @param[in] n  input index of the decomposition
 * @param[in] b  input right hand side
 * @param[out] x solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref multilu_solveh function solves the system
 * \f[ A(n)^H x=b, \f]
 * where \f$ b \f$ and \f$ x \f$ are vectors.
 *
 */
static int multilu_solveh(void *data, mess_int_t n, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct multilu *mlu = ( struct multilu *) data;
    mess_vector y;
    int ret = 0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( n < 0 || n>=mlu->nlu){
        MSG_ERROR("n is out of range ( = " MESS_PRINTF_INT " )\n", n)
        return MESS_ERROR_ARGUMENTS;
    }
    if ( x->dim != mlu->cols){
        MSG_WARN("size x from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", x->dim, mlu->cols);
        mess_vector_resize(x, mlu->cols);
    }
#ifdef _OPENMP
    if ( mlu->background ) {
        // MSG_INFO("waiting to get the solver ready\n");
        pthread_mutex_lock ( mlu->mutex_ready ) ;
        while ( mlu->ready[n] == 0) {
            pthread_cond_wait(mlu->sig_ready, mlu->mutex_ready);
        }
        pthread_mutex_unlock(mlu->mutex_ready);
    }
#endif

    /*-----------------------------------------------------------------------------
     *  real solver and real right hand sides
     *-----------------------------------------------------------------------------*/
    if ( mlu->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)){
        return multilu_solvet(data,n,b,x);
    }

    /*-----------------------------------------------------------------------------
     *  real solver and complex rhs
     *-----------------------------------------------------------------------------*/
    else if (mlu->datatypes[n]==MESS_REAL && MESS_IS_COMPLEX(b)) {
        return multilu_solvet(data,n,b,x);
    }
    /*-----------------------------------------------------------------------------
     * complex solver
     *-----------------------------------------------------------------------------*/
    else {
        // MSG_WARN("complex hermitian solve\n");
        // ret = mess_vector_tocomplex(b); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        ret = mess_vector_init(&y);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(y, b->dim, MESS_COMPLEX);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_perm(b, mlu->q, y);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_perm);
        ret = mess_solver_uhsolve_kernelcsr_complex(mlu->rows, mlu->uvalues_cpx[n], mlu->uind, mlu->up, y->values_cpx);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_lsolve_kernelcsr_real_complex);
        ret = mess_solver_lhsolve_kernelcsr_complex(mlu->rows, mlu->lvalues_cpx[n], mlu->lind, mlu->lp, y->values_cpx);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_usolve_kernelcsr_real_complex);
        ret = mess_vector_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        ret = mess_vector_iperm(y, mlu->p, x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_iperm);
        mess_vector_clear(&y);
    }

    return 0;
}

/**
 * @brief Solve the system \f$ A(n)^T X=B \f$.
 * @param[in] data  input pointer to a multilu
 * @param[in] n  input index of the decomposition
 * @param[in] b  input right hand side
 * @param[out] x solution matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref multilu_solvemt function solves the system
 * \f[A(n)^T X=B,\f]
 * where \f$ B \f$ and \f$ X \f$ are matrices.
 *
 */
static int multilu_solvemt(void *data, mess_int_t n, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct multilu *mlu = ( struct multilu *) data;
    mess_vector y;
    mess_matrix work;
    int conv = -1;
    mess_int_t i,j;
    int ret = 0;
    unsigned short outdt = MESS_REAL;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( n < 0 || n>=mlu->nlu){
        MSG_ERROR("n is out of range ( = " MESS_PRINTF_INT " )\n", n)
        return MESS_ERROR_ARGUMENTS;
    }
    if ( mlu->rows != b->rows) {
        MSG_ERROR("b don't have the right number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }
#ifdef _OPENMP
    if ( mlu->background ) {
        // MSG_INFO("waiting to get the solver ready\n");
        pthread_mutex_lock ( mlu->mutex_ready ) ;
        while ( mlu->ready[n] == 0) {
            pthread_cond_wait(mlu->sig_ready, mlu->mutex_ready);
        }
        pthread_mutex_unlock(mlu->mutex_ready);
    }
#endif



    MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
    MESS_MATRIX_RESET(x);

    if ( mlu->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)){
        outdt  = MESS_REAL;
    } else  if ( mlu->datatypes[n] == MESS_COMPLEX || MESS_IS_COMPLEX(b)){
        outdt = MESS_COMPLEX;
    }
    ret = mess_matrix_alloc(x, mlu->rows, work->cols,mlu->rows*work->cols, MESS_DENSE, outdt );
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    /*-----------------------------------------------------------------------------
     *  real solver and real right hand sides
     *-----------------------------------------------------------------------------*/
    if ( mlu->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)){
        ret = 0 ;
        #ifdef _OPENMP
        #pragma omp parallel for private (y, j, i)  default(shared) reduction(+:ret)
        #endif
        for ( i = 0; i < work->cols; i++) {
            ret += mess_vector_init(&y);
            ret += mess_vector_alloc(y, work->rows, MESS_REAL);
            // perm LU->q
            for (j = 0 ; j < work->rows ; j++) { y->values[j] = work->values[(mlu->q ? mlu->q [j] : j)+i*work->ld] ;}
            ret += mess_solver_utsolve_kernelcsr_real(mlu->rows, mlu->uvalues[n], mlu->uind, mlu->up,y->values);
            ret += mess_solver_ltsolve_kernelcsr_real(mlu->rows, mlu->lvalues[n], mlu->lind, mlu->lp, y->values);
            // iperm LU->p
            for (j = 0 ; j < work->rows; j++){ x->values[(mlu->p ? mlu->p [j] : j)+i*work->ld] = y->values[j] ;}
            mess_vector_clear(&y);
        }
    }

    /*-----------------------------------------------------------------------------
     *  real solver and complex rhs
     *-----------------------------------------------------------------------------*/
    else if (mlu->datatypes[n]==MESS_REAL && MESS_IS_COMPLEX(b)) {
        ret = 0 ;
        #ifdef _OPENMP
        #pragma omp parallel for private (y, j, i)  default(shared) reduction(+:ret)
        #endif
        for ( i = 0; i < work->cols; i++) {
            ret += mess_vector_init(&y);
            ret += mess_vector_alloc(y, work->rows, MESS_COMPLEX);
            // perm LU->q
            for (j = 0 ; j < work->rows ; j++) { y->values_cpx[j] = work->values_cpx[(mlu->q ? mlu->q [j] : j)+i*work->ld] ;}
            ret += mess_solver_utsolve_kernelcsr_real_complex(mlu->rows, mlu->uvalues[n], mlu->uind, mlu->up, y->values_cpx);
            ret += mess_solver_ltsolve_kernelcsr_real_complex(mlu->rows, mlu->lvalues[n], mlu->lind, mlu->lp,y->values_cpx);
            // iperm LU->p
            for (j = 0 ; j < work->rows; j++){ x->values_cpx[(mlu->p ? mlu->p [j] : j)+i*work->ld] = y->values_cpx[j] ;}
            mess_vector_clear(&y);
        }
    }
    /*-----------------------------------------------------------------------------
     * complex solver
     *-----------------------------------------------------------------------------*/
    else {
        // ret = mess_matrix_tocomplex (work); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
        ret = 0 ;
        #ifdef _OPENMP
        #pragma omp parallel for private (y, j, i)  default(shared) reduction(+:ret)
        #endif
        for ( i = 0; i < work->cols; i++) {
            ret += mess_vector_init(&y);
            ret += mess_vector_alloc(y, work->rows, MESS_COMPLEX);
            // perm LU->q
            if (MESS_IS_COMPLEX(work)){
                for (j = 0 ; j < work->rows ; j++) { y->values_cpx[j] = work->values_cpx[(mlu->q ? mlu->q [j] : j)+i*work->ld] ;}
            } else {
                for (j = 0 ; j < work->rows ; j++) { y->values_cpx[j] = work->values[(mlu->q ? mlu->q [j] : j)+i*work->ld] ;}
            }
            ret += mess_solver_utsolve_kernelcsr_complex(mlu->rows, mlu->uvalues_cpx[n], mlu->uind, mlu->up, y->values_cpx);
            ret += mess_solver_ltsolve_kernelcsr_complex(mlu->rows, mlu->lvalues_cpx[n], mlu->lind, mlu->lp,y->values_cpx);
            // iperm LU->p
            for (j = 0 ; j < work->rows; j++){ x->values_cpx[(mlu->p ? mlu->p [j] : j)+i*work->ld] = y->values_cpx[j] ;}
            mess_vector_clear(&y);
        }
    }
    if ( conv == 0) {
        mess_matrix_clear(&work);
    }


    if ( ret!=0) {
        MSG_ERROR("an error occured while solving\n");
        return MESS_ERROR_GENERAL;
    }
    return 0;
}

/**
 * @brief Solve the system \f$ A(n)^H X=B \f$.
 * @param[in] data  input pointer to a multilu
 * @param[in] n  input index of the decomposition
 * @param[in] b  input right hand side
 * @param[out] x solution matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref multilu_solvemh function solves the system
 * \f[A(n)^H X=B,\f]
 * where \f$ B \f$ and \f$ X \f$ are matrices.
 *
 */
static int multilu_solvemh(void *data, mess_int_t n, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct multilu *mlu = ( struct multilu *) data;
    mess_vector y;
    mess_matrix work;
    int conv = -1;
    mess_int_t i,j;
    int ret = 0;
    unsigned short outdt = MESS_REAL;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( n < 0 || n>=mlu->nlu){
        MSG_ERROR("n is out of range ( = " MESS_PRINTF_INT " )\n", n)
        return MESS_ERROR_ARGUMENTS;
    }
    if ( mlu->rows != b->rows) {
        MSG_ERROR("b don't have the right number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }
#ifdef _OPENMP
    if ( mlu->background ) {
        // MSG_INFO("waiting to get the solver ready\n");
        pthread_mutex_lock ( mlu->mutex_ready ) ;
        while ( mlu->ready[n] == 0) {
            pthread_cond_wait(mlu->sig_ready, mlu->mutex_ready);
        }
        pthread_mutex_unlock(mlu->mutex_ready);
    }
#endif



    if ( mlu->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)){
        outdt  = MESS_REAL;
    } else  if ( mlu->datatypes[n] == MESS_COMPLEX || MESS_IS_COMPLEX(b)){
        outdt = MESS_COMPLEX;
    }
    /*-----------------------------------------------------------------------------
     *  real solver and real right hand sides
     *-----------------------------------------------------------------------------*/
    if ( mlu->datatypes[n] == MESS_REAL && MESS_IS_REAL(b)){
        return multilu_solvemt(data, n, b,x);
    }

    /*-----------------------------------------------------------------------------
     *  real solver and complex rhs
     *-----------------------------------------------------------------------------*/
    else if (mlu->datatypes[n]==MESS_REAL && MESS_IS_COMPLEX(b)) {
        return multilu_solvemt(data, n, b,x);
    }
    /*-----------------------------------------------------------------------------
     * complex solver
     *-----------------------------------------------------------------------------*/
    else {
        MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
        MESS_MATRIX_RESET(x);
        ret = mess_matrix_alloc(x, mlu->rows, work->cols,mlu->rows*work->cols, MESS_DENSE, outdt );
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

        // ret = mess_matrix_tocomplex (work); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
        ret = 0 ;
        #ifdef _OPENMP
        #pragma omp parallel for private (y, j, i)  default(shared) reduction(+:ret)
        #endif
        for ( i = 0; i < work->cols; i++) {
            ret += mess_vector_init(&y);
            ret += mess_vector_alloc(y, work->rows, MESS_COMPLEX);
            // perm LU->q
            if (MESS_IS_COMPLEX(work)){
                for (j = 0 ; j < work->rows ; j++) { y->values_cpx[j] = work->values_cpx[(mlu->q ? mlu->q [j] : j)+i*work->ld] ;}
            } else {
                for (j = 0 ; j < work->rows ; j++) { y->values_cpx[j] = work->values[(mlu->q ? mlu->q [j] : j)+i*work->ld] ;}
            }
            ret += mess_solver_uhsolve_kernelcsr_complex(mlu->rows, mlu->uvalues_cpx[n], mlu->uind, mlu->up, y->values_cpx);
            ret += mess_solver_lhsolve_kernelcsr_complex(mlu->rows, mlu->lvalues_cpx[n], mlu->lind, mlu->lp,y->values_cpx);
            // iperm LU->p
            for (j = 0 ; j < work->rows; j++){ x->values_cpx[(mlu->p ? mlu->p [j] : j)+i*work->ld] = y->values_cpx[j] ;}
            mess_vector_clear(&y);
        }
    }
    if ( conv == 0) {
        mess_matrix_clear(&work);
    }


    if ( ret!=0) {
        MSG_ERROR("an error occured while solving\n");
        return MESS_ERROR_GENERAL;
    }
    return 0;
}

typedef struct jobdata_st {
    mess_int_t rows;
    mess_int_t cols;
    void *workvalues;
    void *vwork;
    mess_int_t *workcolptr;
    mess_int_t *workrowptr;
    mess_int_t *lind;
    mess_int_t *lp;
    mess_int_t *uind;
    mess_int_t *up;
    mess_int_t *pperm;
    mess_int_t *qperm;
    void *lvalues;
    void *uvalues;
    void *lvalues_next;
    void *uvalues_next;
    mess_int_t index;
    int *ret;
    struct multilu *mlu;
    int pair;
} jobdata;

static void job_decompose_real(void *data){
    MSG_FNAME(__func__);
    jobdata * job = (jobdata *) data;
    MSG_INFO("Start shift " MESS_PRINTF_INT "\n", job->index);
    // showSchedParam(pthread_self());
    // setPriority(pthread_self(), 90);
    // showSchedParam(pthread_self());

    // usleep((job->index-1)*1000);
    mess_decomp_lureuse_kernelcsr(job->rows, job->cols,(double *) job->vwork,  job->workcolptr, job->workrowptr,
                        job->lind, job->lp,                         // L pattern
                        job->uind, job->up,                     // U pattern
                        job->pperm, job->qperm,                                 // permuatations
                        (double *)job->lvalues, (double *) job->uvalues);
    mess_free(job->vwork);

    if ( job->pair ) {
        mess_int_t lnz = job->lp[job->rows];
        mess_int_t unz = job->up[job->rows];
        double *lv, *lvn;
        double *uv, *uvn;
        lv = (double  *) job->lvalues;
        lvn = (double *) job->lvalues_next;
        uv = (double  *) job->uvalues;
        uvn = (double *) job->uvalues_next;

        memcpy(lvn , lv, sizeof(double) * lnz);
        memcpy(uvn, uv, sizeof(double)*unz);
    }


    if (job->mlu->background){
        pthread_mutex_lock(job->mlu->mutex_ready);
        job->mlu->ready[job->index] = 1;
        if (job->pair) job->mlu->ready[job->index+1] = 1;

        pthread_cond_broadcast(job->mlu->sig_ready);
        pthread_mutex_unlock(job->mlu->mutex_ready);


    }
    mess_free(job);
}

static void job_decompose_complex(void *data){
    MSG_FNAME(__func__);
    jobdata * job = (jobdata *) data;
    // usleep((job->index-1)*1000);
    MSG_INFO("Start shift " MESS_PRINTF_INT "\n", job->index);
    // showSchedParam(pthread_self());

    mess_decomp_lureuse_kernelcsr_complex(job->rows, job->cols,(mess_double_cpx_t*) job->vwork,  job->workcolptr, job->workrowptr,
                        job->lind, job->lp,                         // L pattern
                        job->uind, job->up,                     // U pattern
                        job->pperm, job->qperm,                                 // permuatations
                        (mess_double_cpx_t*)job->lvalues, (mess_double_cpx_t*) job->uvalues);
    mess_free(job->vwork);
    if ( job->pair ) {
        mess_int_t lnz = job->lp[job->rows];
        mess_int_t unz = job->up[job->rows];
        mess_int_t i;
        mess_double_cpx_t *lv, *lvn;
        mess_double_cpx_t *uv, *uvn;
        lv = (mess_double_cpx_t *) job->lvalues;
        lvn = (mess_double_cpx_t *) job->lvalues_next;
        uv = (mess_double_cpx_t *) job->uvalues;
        uvn = (mess_double_cpx_t *) job->uvalues_next;

        for (i=0; i<lnz; i++){
            lvn [i] = conj (lv[i]);
        }
        for (i=0; i<unz; i++){
            uvn [i] = conj (uv[i]);
        }
    }
    if (job->mlu->background){
        pthread_mutex_lock(job->mlu->mutex_ready);
        job->mlu->ready[job->index] = 1;
        if (job->pair) job->mlu->ready[job->index+1] = 1;
        pthread_cond_broadcast(job->mlu->sig_ready);
        pthread_mutex_unlock(job->mlu->mutex_ready);

    }
    mess_free(job);
}



/**
 * @brief Create a single-pattern-multi-value LU decomposition of \f$ (sA+pI) \f$ or \f$ (sA+pE) \f$.
 * @param[in] matrix  input system matrix
 * @param[in] shiftsl  input left shifts \f$ s \f$, in front of \f$ A \f$, @c NULL is not given
 * @param[in] shiftsr  input right shifts \f$ p \f$, in front of \f$ E \f$
 * @param[out] mlu  multidirect solver to create
 * @param[in] base_input    existing solver to get the pattern if it is possible
 * @param[in] shiftmatrix_input matrix to shift if it points to @c NULL, the identity matrix is used
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_multidirect_create_sparse_lu function create a single pattern multi values LU decomposition. \n
 * If you specify a base solver you have to ensure that the base solver is created from a matrix with the same pattern
 * as matrix. \n
 * If you do not want to use a base solver set it to @c NULL. \n
 * If you want to shift with the identity matrix set @p shiftmatrix_input to @c NULL.
 */
int mess_multidirect_create_sparse_lu(mess_matrix matrix, mess_vector shiftsl, mess_vector shiftsr , mess_multidirect mlu, mess_direct base_input, mess_matrix shiftmatrix_input) {
    MSG_FNAME(__func__);
    mess_matrix work = NULL;
    mess_matrix shiftmatrix = NULL;
    struct multilu * m;
    mess_int_t *dpos;
    mess_int_t nshifts;
    mess_int_t i, j;
    mess_int_t *lind, *lp, *uind, *up ;
    mess_int_t *pperm, *qperm;
    int use_base = 0;
    int use_gernalized = 0;
    int ret = 0;
    int clear_base = 0;
    int perm_before = 0;
    mess_direct base = NULL;
    mess_vector shiftsr_real, shiftsl_real;
    mess_vector shiftsr_complex, shiftsl_complex;
    int use_complex = 0;
    double eps = mess_eps();
    mess_int_t ii;
    int havelshift = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer (matrix);
    mess_check_nullpointer (shiftsr);
    mess_check_nullpointer (mlu);

    // DTINIT(matrix);
    nshifts = shiftsr->dim;
    if ( base_input != NULL && base_input->getL != NULL && base_input->getU != NULL ) {
        use_base  = 1;
    }
    if ( nshifts <=0 ) {
        MSG_ERROR("nshifts hasn't a useful value ( = " MESS_PRINTF_INT " )\n", nshifts);
        return MESS_ERROR_ARGUMENTS;
    }
    if ( shiftmatrix_input != NULL ) {
        MSG_INFO("using mass matrix to shift\n");
        if ( shiftmatrix_input->rows != matrix->rows || shiftmatrix_input->cols != matrix->cols) {
            MSG_ERROR("shiftmatrix and system matrix don't have the same dimension.\n");
            return MESS_ERROR_DIMENSION;
        }
        use_gernalized = 1;
    }
    if (shiftsl != NULL ) {
        if ( shiftsl->dim != shiftsr->dim) {
            MSG_ERROR("The left and right shift vector must have the same dimension. shiftsl->dim=" MESS_PRINTF_INT ", \t shiftsr->dim=" MESS_PRINTF_INT "\n", shiftsl->dim, shiftsr->dim);
            return MESS_ERROR_DIMENSION;
        }
        havelshift = 1;
    }

    mess_try_alloc(m, struct multilu*, sizeof(struct multilu));

    /*-----------------------------------------------------------------------------
     *  create input data
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&work);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    if ( MESS_IS_CSR(matrix)) {
        ret = mess_matrix_copy(matrix, work);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    } else {
        ret = mess_matrix_convert(matrix, work, MESS_CSR);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    }
    if ( use_gernalized ) {
        ret = mess_matrix_init(&shiftmatrix);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        if ( MESS_IS_CSR(shiftmatrix_input)) {
            ret = mess_matrix_copy(shiftmatrix_input, shiftmatrix);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        } else {
            ret = mess_matrix_convert(shiftmatrix_input, shiftmatrix, MESS_CSR);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
        }
    }
    // real and complex shifts
    MESS_INIT_VECTORS(&shiftsr_real, &shiftsr_complex);
    ret = mess_vector_alloc(shiftsr_real, nshifts, MESS_REAL);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_vector_alloc);
    ret = mess_vector_alloc(shiftsr_complex, nshifts, MESS_COMPLEX);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_vector_alloc);
    ret = mess_vector_copy(shiftsr, shiftsr_real);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy);
    ret = mess_vector_copy(shiftsr, shiftsr_complex);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy);
    ret = mess_vector_toreal_nowarn(shiftsr_real);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal_nowarn);
    ret = mess_vector_tocomplex(shiftsr_complex);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
    if (havelshift){
        MESS_INIT_VECTORS(&shiftsl_real,&shiftsl_complex);
        ret = mess_vector_alloc(shiftsl_real, nshifts, MESS_REAL);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_vector_alloc);
        ret = mess_vector_alloc(shiftsl_complex, nshifts, MESS_COMPLEX);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_vector_alloc);
        ret = mess_vector_copy(shiftsl, shiftsl_real);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy);
        ret = mess_vector_copy(shiftsl, shiftsl_complex);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy);
        ret = mess_vector_toreal_nowarn(shiftsl_real);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal_nowarn);
        ret = mess_vector_tocomplex(shiftsl_complex);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
    }

    ret  =  multilu_init(m);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), multilu_init);

    // mess_matrix_toreal(work);
    // DTINIT(work);
    // mess_matrix_printinfo(shiftmatrix_input);
    // mess_matrix_printinfo(shiftmatrix);

    /*-----------------------------------------------------------------------------
     *  select datatype
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(m->datatypes, unsigned short*, sizeof(unsigned short)*nshifts);
    // mess_matrix_printinfo(work);
    // mess_matrix_printinfo(shiftmatrix);
    if (MESS_IS_REAL(work) /* && MESS_IS_REAL(shiftmatrix)*/ ){
        if ( havelshift ) {
            for (i = 0 ; i < nshifts; i++)  {
                if ( IS_REAL(shiftsr_complex->values_cpx[i]) && IS_REAL(shiftsl_complex->values_cpx[i])){
                    MSG_INFO("shiftr[" MESS_PRINTF_INT "] = %lg + %lgi, shiftl[" MESS_PRINTF_INT "] = %lg = %lgi are real\n",i, creal(shiftsr_complex->values_cpx[i]), cimag(shiftsr_complex->values_cpx[i]),i , creal(shiftsl_complex->values_cpx[i]), cimag(shiftsl_complex->values_cpx[i]));
                    m->datatypes[i]=MESS_REAL;
                } else {
                    MSG_INFO("shiftr[" MESS_PRINTF_INT "] = %lg + %lgi, shiftl[" MESS_PRINTF_INT "] = %lg = %lgi are complex\n",i, creal(shiftsr_complex->values_cpx[i]), cimag(shiftsr_complex->values_cpx[i]), i, creal(shiftsl_complex->values_cpx[i]), cimag(shiftsl_complex->values_cpx[i]));
                    m->datatypes[i]=MESS_COMPLEX;
                }
            }

        } else {
            for (i = 0 ; i < nshifts; i++)  {
                if ( IS_REAL(shiftsr_complex->values_cpx[i])){
                    MSG_INFO("shift[" MESS_PRINTF_INT "] = %lg + %lgi is real\n",i, creal(shiftsr_complex->values_cpx[i]), cimag(shiftsr_complex->values_cpx[i]));
                    m->datatypes[i]=MESS_REAL;
                } else {
                    MSG_INFO("shift[" MESS_PRINTF_INT "] = %lg + %lgi is complex\n",i, creal(shiftsr_complex->values_cpx[i]), cimag(shiftsr_complex->values_cpx[i]));
                    m->datatypes[i]=MESS_COMPLEX;
                }
            }
        }
        use_complex = 0;
    } else {
        MSG_INFO("all shifts are complex");
        for (i=0;i<nshifts;i++) m->datatypes[i]=MESS_COMPLEX;
        use_complex = 1;
    }



    /*-----------------------------------------------------------------------------
     *  check if all diagonal elements exist
     *-----------------------------------------------------------------------------*/
    if ( !use_gernalized ){
        int problems = 0;
        mess_try_alloc(dpos, mess_int_t *, sizeof ( mess_int_t) * work->rows);
        ret =   mess_matrix_diagpos(work , dpos);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_diagpos);
        for ( i = 0; i < work->rows; i++) {
            if ( dpos[i] < 0) {
                // MSG_WARN("diagonal entry in line " MESS_PRINTF_INT " doesn't exist.\n", i);
                problems ++;
            }
        }
        if ( problems > 0 ) {
            MSG_WARN("some diagonal entries don't exists. creating identity mass matrix.\n");
            ret = mess_matrix_init(&shiftmatrix);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
            if ( use_complex == 1) {
                ret = mess_matrix_eyec(shiftmatrix, work->rows, work->cols, MESS_CSR);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_eyec);
            } else {
                ret = mess_matrix_eye(shiftmatrix, work->rows, work->cols, MESS_CSR);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_eye);
            }
            perm_before = 1;
            use_base = 0;
            use_gernalized = 1;
        }
        mess_free(dpos);
    }

    /*-----------------------------------------------------------------------------
     *  generate the pattern
     *-----------------------------------------------------------------------------*/
    if ( use_gernalized && mess_matrix_equalpattern(work,shiftmatrix) != 1) {
        use_base = 0;
        MSG_WARN("The pattern of System and Shiftmatrix aren't the same. Can not use base solver.\n");
    }
    /*-----------------------------------------------------------------------------
     *  use a given base solver, only in case of eye shifted
     *-----------------------------------------------------------------------------*/
    if ( use_base  ==  1 ){
        MSG_INFO("use the given base solver, name = %s \n", base_input->name);
        // printf("use base name = %s \n", base_input->name);
        // printf("base = %lx \n", base_input );
        base = base_input;
        clear_base = 0;
    }
    /*-----------------------------------------------------------------------------
     *  gernerate a new solver
     *-----------------------------------------------------------------------------*/
    else {
        MSG_INFO("create solver...\n");
        clear_base = 1;
        ret = mess_direct_init(&base);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init) ;

        if ( use_gernalized ){
            MSG_INFO("create generalized solver...\n");
            if ( use_complex == 0) {
                mess_matrix work2;
                // mess_matrix_toreal(work);

                ret = mess_matrix_init(&work2);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_copy(work, work2);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
                ret = mess_matrix_add(0.0,shiftmatrix, 1, work);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                if (havelshift) {
                    ret = mess_matrix_add(shiftsr_real->values[0],shiftmatrix, shiftsl_real->values[0], work2);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                } else {
                    ret = mess_matrix_add(shiftsr_real->values[0],shiftmatrix, 1 , work2);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                }
                ret = mess_matrix_add(0,work, 1, shiftmatrix);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                ret = mess_matrix_sort(work2);                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
                ret = mess_matrix_sort(work);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
                ret = mess_matrix_sort(shiftmatrix);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
                //#ifdef MESS_HAVE_CSPARSE
                // Don't use umfpack because the generated solver isn't useable.
                //ret = mess_direct_create_csparse_lu(work2, base); FUNCTION_FAILURE_HANDLE (ret, (ret!=0), mess_direct_create_csparse_lu);
                // ret = mess_direct_create_umfpack(work2, base); FUNCTION_FAILURE_HANDLE (ret, (ret!=0), mess_direct_create_umfpack);
                //#else
                ret =  mess_direct_create_sparse_lu(work2,base);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sparse_lu);
                //#endif
                mess_matrix_clear(&work2);

            } else {
                // ret = mess_matrix_addc(0.0,shiftmatrix, 1, work); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
                // ret = mess_matrix_addc(0,work, 1, shiftmatrix); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
                // ret = mess_matrix_sort(work);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
                // ret = mess_matrix_sort(shiftmatrix); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
                mess_matrix work2;
                // mess_matrix_toreal(work);
                // mess_matrix_toreal(shiftmatrix);

                ret = mess_matrix_init(&work2);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_copy(work, work2);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
                ret = mess_matrix_add(0.0,shiftmatrix, 1, work);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                if ( havelshift) {
                    ret = mess_matrix_addc(shiftsr_complex->values_cpx[0],shiftmatrix, shiftsl_complex->values_cpx[0], work2);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                } else {
                    ret = mess_matrix_addc(shiftsr_complex->values_cpx[0],shiftmatrix, 1 , work2);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                }
                ret = mess_matrix_add(0,work, 1, shiftmatrix);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                ret = mess_matrix_sort(work2);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
                ret = mess_matrix_sort(work);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
                ret = mess_matrix_sort(shiftmatrix);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);

                /* #ifdef MESS_HAVE_CSPARSE
                // Don't use umfpack because the generated solver isn't useable.
                ret = mess_direct_create_csparse_lu(work2, base); FUNCTION_FAILURE_HANDLE (ret, (ret!=0), mess_direct_create_csparse_lu);
                // ret = mess_direct_create_umfpack(work2, base); FUNCTION_FAILURE_HANDLE (ret, (ret!=0), mess_direct_create_umfpack);
                #else */
                ret = mess_direct_create_sparse_lu(work2, base);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sparse_lu);
                // #endif
                mess_matrix_clear(&work2);

            }
        } else {
            if ( use_complex == 1) {
                #ifdef MESS_HAVE_CSPARSE
                // ret = mess_direct_create_csparse_lu_complex(work, base);
                ret = mess_direct_create_umfpack(work, base);
                FUNCTION_FAILURE_HANDLE (ret, (ret!=0), mess_direct_umpack);
                #else
                ret = mess_direct_create_sparse_lu(work, base);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sparse_lu);
                #endif
            } else {
                #ifdef MESS_HAVE_UMFPACK
                // ret = mess_direct_create_csparse_lu(work, base);
                ret = mess_direct_create_umfpack(work, base);
                FUNCTION_FAILURE_HANDLE (ret, (ret!=0), mess_direct_umpack);
                #else
                ret = mess_direct_create_sparse_lu(work, base);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sparse_lu);
                #endif
            }
            perm_before = 0 ;
        }
    }

    /*-----------------------------------------------------------------------------
     *  Use a exsisting solver to get the pattern
     *-----------------------------------------------------------------------------*/
    mess_matrix L=NULL, U=NULL;
    ret = mess_matrix_init(&L); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&U); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

    // get L
    ret  = mess_direct_getL(base, L); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_getL);
    if ( !MESS_IS_CSR(L) ) {
        MSG_ERROR ("at the moment only CSR in base\n");
        return MESS_ERROR_STORAGETYPE;
    }
    m->lnnz = L->nnz; m->lind= L->colptr ; L->colptr = NULL;
    m->lp= L->rowptr ; L->rowptr = NULL;
    mess_matrix_clear(&L);

    // get U
    ret = mess_direct_getU(base,U); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_getU);
    if ( !MESS_IS_CSR(U) ) {
        MSG_ERROR ("at the moment only CSR in base\n");
        return MESS_ERROR_STORAGETYPE;
    }
    m->unnz = U->nnz; m->uind= U->colptr ; U->colptr = NULL;
    m->up= U->rowptr ; U->rowptr = NULL;
    mess_matrix_clear(&U);

    // get permutations
    mess_try_alloc ( m->p, mess_int_t *, sizeof(mess_int_t) * work->rows);
    mess_try_alloc ( m->q, mess_int_t *, sizeof(mess_int_t) * work->cols);
    ret = mess_direct_getpermp(base, m->p); FUNCTION_FAILURE_HANDLE(ret, (ret!=0) , mess_direct_getpermp);
    ret = mess_direct_getpermq(base, m->q); FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_direct_getpermq);
    m->pinv = pinv(m->p, work->rows);
    m->qinv = pinv(m->q, work->cols);

    m->store_type = MESS_CSR;
    m->nlu = 0;
    m->rows= work->rows;
    m->cols= work->cols;


    /*-----------------------------------------------------------------------------
     *  generate the solvers.
     *-----------------------------------------------------------------------------*/
    if ( !use_gernalized ){
        mess_try_alloc(dpos, mess_int_t *, sizeof ( mess_int_t) * work->rows);
        ret =   mess_matrix_diagpos(work , dpos); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_diagpos);
    }

    mess_try_alloc(m->lvalues,  double **, sizeof (double *) * nshifts);
    mess_try_alloc(m->uvalues,  double **, sizeof (double *) * nshifts);
    mess_try_alloc(m->lvalues_cpx,  mess_double_cpx_t **, sizeof (mess_double_cpx_t *) * nshifts);
    mess_try_alloc(m->uvalues_cpx,  mess_double_cpx_t **, sizeof (mess_double_cpx_t *) * nshifts);
    for ( i  = 0; i < nshifts ; i++){
        m->lvalues[i]=NULL;
        m->lvalues_cpx[i]=NULL;
        m->uvalues[i]=NULL;
        m->uvalues_cpx[i]=NULL;
    }

    m->nlu = nshifts;

    lind = m->lind;
    uind = m->uind;
    lp = m-> lp;
    up = m-> up;


    /*-----------------------------------------------------------------------------
     *  permute the matrix and the shift if we need it
     *-----------------------------------------------------------------------------*/
    /* for ( i = 0; i < 20; i++) {
        printf(" p [ %6ld ] = " MESS_PRINTF_INT " \t q[ %6ld ] = " MESS_PRINTF_INT "\n", i, m->p[i] , i,  m->q[i]);
    } */
    // perm_before = 1;
    if ( perm_before == 1){
        MSG_INFO("permute system\n");
        ret = mess_matrix_perm(work, m->p, m->q);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_perm);
        if (use_gernalized) {
            ret = mess_matrix_perm(shiftmatrix, m->p,m->q);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_perm);
        }
        pperm = NULL;
        qperm = NULL;
    } else {
        pperm = m->p;
        qperm = m->qinv;
    }

    m->background =  WORK_BACKGROUND;
#ifndef _OPENMP
    m->background = 0;
#else
    if ( m->background ) {
        int thr;
        MSG_INFO("preparing background mode\n");
        thr = omp_get_num_procs()-1;
        if (thr ==0 ) thr =1 ;
        // if (thr > 4 ) thr = 3 + thr/2-1;
        MSG_INFO("create a threadpool with %d jobs\n",thr);
        ret = mess_threadpool_init(&(m->solverpool),thr,2*nshifts); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_threadpool_init);
        mess_try_alloc(m->sig_ready, pthread_cond_t *, sizeof(pthread_cond_t));
        mess_try_alloc(m->mutex_ready, pthread_mutex_t *, sizeof(pthread_mutex_t));
        pthread_cond_init(m->sig_ready, NULL);
        pthread_mutex_init(m->mutex_ready, NULL);
        mess_try_alloc(m->ready, mess_int_t *, sizeof(mess_int_t)*nshifts);
        for ( i = 0 ; i < nshifts; i++) m->ready[i] = 0;

    }
#endif

    /*-----------------------------------------------------------------------------
     *  hold all in memory
     *-----------------------------------------------------------------------------*/
    if ( m->background == 0 ) {
        mess_int_t *mapping;
        mess_int_t mapcount = 0;
        mess_try_alloc(mapping, mess_int_t *, nshifts*sizeof(mess_int_t));
        for (i = 0; i < nshifts; i++){
            if ( m->datatypes[i] == MESS_REAL ){
                mess_try_calloc(m->lvalues [i] , double *, m->lnnz, sizeof (double));
                mess_try_calloc(m->uvalues [i] , double *, m->unnz, sizeof (double));
                m->lvalues_cpx[i] = NULL;
                m->uvalues_cpx[i] = NULL;
            } else {
                mess_try_calloc(m->lvalues_cpx [i] , mess_double_cpx_t *, m->lnnz, sizeof (mess_double_cpx_t));
                mess_try_calloc(m->uvalues_cpx [i] , mess_double_cpx_t *, m->unnz, sizeof (mess_double_cpx_t));
                m->lvalues[i] = NULL;
                m->uvalues[i] = NULL;
            }
        }
        if ( havelshift ) {
            for ( i =0; i < nshifts; i++) {
                mapping[mapcount++]=i;
            }
        } else {
            for (i=0; i <nshifts; i++){
                if ((i <nshifts-1) && m->datatypes[i] == MESS_COMPLEX && cabs(shiftsr_complex->values_cpx[i]-conj(shiftsr_complex->values_cpx[i+1])) < eps ) {
                    MSG_INFO(" " MESS_PRINTF_INT " complex pair. (diff = %lg)\n", i, cabs(shiftsr_complex->values_cpx[i]-conj(shiftsr_complex->values_cpx[i+1]))) ;
                    mapping[mapcount++]=i;
                    i++;
                } else if ((i < nshifts-1) && m->datatypes[i]== MESS_REAL && (fabs(shiftsr_real->values[i]-shiftsr_real->values[i+1]) < 2*eps )){
                    MSG_INFO(" " MESS_PRINTF_INT " real pair. \n",i) ;
                    mapping[mapcount++]=i;
                    i++;
                } else {
                    mapping[mapcount++]=i;
                }
            }
        }

#ifdef _OPENMP
        #pragma omp parallel for private(ii,i,j) default(shared)
#endif
        for ( ii = 0; ii < mapcount; ii++) {
            i = mapping[ii];
            jobdata *job;
            mess_try_alloc2(job, jobdata *, sizeof(jobdata));
            job->rows = work->rows;
            job->cols = work->cols;
            job->workrowptr = work->rowptr;
            job->workcolptr = work->colptr;
            job->lind = lind;
            job->lp = lp;
            job->uind = uind;
            job->up = up;
            job->pperm = pperm ;
            job->qperm = qperm ;
            job->ret = &ret;
            job->mlu = m ;
            job->index = i ;
            job->pair = 0;

            /*-----------------------------------------------------------------------------
             *  real system matrix and real shifts
             *-----------------------------------------------------------------------------*/
            if ( m->datatypes[i] == MESS_REAL) {
                double *vwork;
                mess_try_alloc2 ( vwork ,  double *, work->nnz*sizeof(double));

                if (use_gernalized == 0) {
                    memcpy(vwork, work->values, work->nnz * sizeof(double));
                    if (havelshift) {
                        for ( j =0; j < work->nnz; j++) vwork[j] *= shiftsl_real->values[i];
                    }
                    for ( j = 0 ; j < work->rows; j++)  vwork[dpos[j]] += shiftsr_real->values[i];
                } else {
                    if ( havelshift ) {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = shiftsl_real->values[i]*work->values [j] + shiftsr_real->values[i]*shiftmatrix->values[j];
                    } else {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = work->values [j] + shiftsr_real->values[i]*shiftmatrix->values[j];
                    }
                }
                if ( (i < nshifts-1) && (fabs(shiftsr_real->values[i]-shiftsr_real->values[i+1]) < 2*eps ) && havelshift == 0 ) {
                    job->pair = 1;
                } else {
                    job->pair = 0;
                }
                // mess_try_alloc(job, jobdata *, sizeof(jobdata));
                job->vwork = (void *) vwork;
                job->lvalues = (void *) m->lvalues[i];
                job->uvalues = (void *) m->uvalues[i];
                if ( i < nshifts -1 ){
                    job->lvalues_next = (void *) m->lvalues[i+1];
                    job->uvalues_next = (void *) m->uvalues[i+1];
                }
                job_decompose_real(job);
                /* @ref mess_decomp_lureuse_kernelcsr(work->rows, work->cols,   vwork,  work->colptr, work->rowptr,
                        lind, lp,                       // L pattern
                        uind, up,                       // U pattern
                        pperm, qperm,                                   // permuatations
                        m->lvalues[i], m->uvalues[i]);

                mess_free(vwork); */
            }
            /*-----------------------------------------------------------------------------
             *  real system matrix and complex shifts
             *-----------------------------------------------------------------------------*/
            else if ( m->datatypes[i] == MESS_COMPLEX && use_complex == 0){
                mess_double_cpx_t *vwork;
                mess_try_alloc2(vwork, mess_double_cpx_t *, work->nnz *sizeof(mess_double_cpx_t));
                if (use_gernalized == 0) {
                    if (havelshift) {
                        for ( j =0; j < work->nnz; j++) vwork[j] = shiftsl_complex->values_cpx[i]*work->values[j];
                    } else {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = work->values[j];
                    }
                    for ( j = 0 ; j < work->rows; j++)  vwork[dpos[j]] += shiftsr_complex->values_cpx[i];
                } else {
                    if ( havelshift) {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = shiftsl_complex->values_cpx[i]* work->values [j] + shiftsr_complex->values_cpx[i]*shiftmatrix->values[j];
                    } else {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = work->values [j] + shiftsr_complex->values_cpx[i]*shiftmatrix->values[j];
                    }
                }

                job->vwork = (void *) vwork;
                job->lvalues = (void *) m->lvalues_cpx[i];
                job->uvalues = (void *) m->uvalues_cpx[i];
                if ( i < nshifts - 1 ) {
                    job->lvalues_next = (void *) m->lvalues_cpx[i+1];
                    job->uvalues_next = (void *) m->uvalues_cpx[i+1];
                }
                if ((i < nshifts-1) && (cabs(shiftsr_complex->values_cpx[i]-conj(shiftsr_complex->values_cpx[i+1])) < eps )&&havelshift == 0 ) {
                    job->pair = 1;
                } else {
                    job->pair = 0;
                }
                job_decompose_complex(job);

                /* @ref mess_decomp_lureuse_kernelcsr_complex(work->rows, work->cols,   vwork,  work->colptr, work->rowptr,
                            lind, lp,                       // L pattern
                            uind, up,                       // U pattern
                            pperm, qperm,                                   // permuatations
                            m->lvalues_cpx[i], m->uvalues_cpx[i]);
                mess_free(vwork); */
            }

            /*-----------------------------------------------------------------------------
             *  complex system matrix
             *-----------------------------------------------------------------------------*/
            else {
                mess_double_cpx_t *vwork;
                mess_try_alloc2(vwork, mess_double_cpx_t *, work->nnz* sizeof(mess_double_cpx_t));
                if (use_gernalized == 0) {

                    memcpy(vwork, work->values_cpx, work->nnz * sizeof(mess_double_cpx_t));
                    if ( havelshift ) {
                        for ( j = 0 ; j < work->nnz ; j++ ) vwork[j] *= shiftsl_complex->values_cpx[i];
                    }
                    for ( j = 0 ; j < work->rows; j++)  vwork[dpos[j]] += shiftsr_complex->values_cpx[i];
                } else {
                    if ( havelshift ) {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = shiftsl_complex->values_cpx[j]*work->values_cpx [j] + shiftsr_complex->values_cpx[i]*shiftmatrix->values_cpx[j];
                    } else {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = work->values_cpx [j] + shiftsr_complex->values_cpx[i]*shiftmatrix->values_cpx[j];
                    }
                }
                job->vwork = (void *) vwork;
                job->lvalues = (void *) m->lvalues_cpx[i];
                job->uvalues = (void *) m->uvalues_cpx[i];
                if ( i < nshifts - 1) {
                    job->lvalues_next = (void *) m->lvalues_cpx[i+1];
                    job->uvalues_next = (void *) m->uvalues_cpx[i+1];
                }
                if ((i < nshifts -1)&& (cabs(shiftsr_complex->values_cpx[i]-conj(shiftsr_complex->values_cpx[i+1])) < eps) && havelshift == 0 ) {
                    job->pair = 1;
                } else {
                    job->pair = 0;
                }

                job_decompose_complex(job);

                /* @ref mess_decomp_lureuse_kernelcsr_complex(work->rows, work->cols,   vwork,  work->colptr, work->rowptr,
                            lind, lp,                       // L pattern
                            uind, up,                       // U pattern
                            pperm, qperm,                                   // permuatations
                            m->lvalues_cpx[i], m->uvalues_cpx[i]);
                mess_free(vwork); */

            }
            // mess_free(job);
        }
        mess_free(mapping);
    }
#ifdef _OPENMP
    /*-----------------------------------------------------------------------------
     *  back ground mode
     *-----------------------------------------------------------------------------*/
    else if (m->background) {
        MSG_INFO("decompose in background mode\n");

        for (i = 0; i < nshifts; i++){
            if ( m->datatypes[i] == MESS_REAL ){
                mess_try_calloc(m->lvalues [i] , double *, m->lnnz, sizeof (double));
                mess_try_calloc(m->uvalues [i] , double *, m->unnz, sizeof (double));
                m->lvalues_cpx[i] = NULL;
                m->uvalues_cpx[i] = NULL;
            } else {
                mess_try_calloc(m->lvalues_cpx [i] , mess_double_cpx_t *, m->lnnz, sizeof (mess_double_cpx_t));
                mess_try_calloc(m->uvalues_cpx [i] , mess_double_cpx_t *, m->unnz, sizeof (mess_double_cpx_t));
                m->lvalues[i] = NULL;
                m->uvalues[i] = NULL;
            }
        }


        for ( i = 0; i < nshifts; i++) {
            jobdata * job ;
            mess_try_alloc2( job, jobdata *, sizeof(jobdata));
            job->rows = work->rows;
            job->cols = work->cols;
            job->workrowptr = work->rowptr;
            job->workcolptr = work->colptr;
            job->lind = lind;
            job->lp = lp;
            job->uind = uind;
            job->up = up;
            job->pperm = pperm ;
            job->qperm = qperm ;
            job->ret = &ret;
            job->mlu = m ;
            job->index = i ;
            job->pair = 0 ;

            /*-----------------------------------------------------------------------------
             *  real system matrix and real shifts
             *-----------------------------------------------------------------------------*/
            if ( m->datatypes[i] == MESS_REAL) {
                double *vwork;
                mess_try_alloc2(vwork,  double *, work->nnz*sizeof(double));
                mess_int_t id;

                if (use_gernalized == 0) {
                    if ( havelshift ) {
                        for ( j = 0 ; j < work->nnz; j++) vwork[j] = shiftsl_real->values[i] * work->values[j];
                    } else {
                        memcpy(vwork, work->values, work->nnz * sizeof(double));
                    }
                    for ( j = 0 ; j < work->rows; j++)  vwork[dpos[j]] += shiftsr_real->values[i];
                } else {
                    if ( havelshift ) {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = shiftsl_real->values[i]*work->values [j] + shiftsr_real->values[i]*shiftmatrix->values[j];
                    } else {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = work->values [j] + shiftsr_real->values[i]*shiftmatrix->values[j];
                    }
                }
                job->vwork = (void *) vwork;
                job->lvalues = (void *) m->lvalues[i];
                job->uvalues = (void *) m->uvalues[i];
                job->lvalues_next = (void *) m->lvalues[i+1];
                job->uvalues_next = (void *) m->uvalues[i+1];

                // job_decompose_real(job);
                if ( (fabs(shiftsr_real->values[i]-shiftsr_real->values[i+1]) < 2*eps )&& havelshift == 0 ) {
                    job->pair = 1;
                    i++;
                } else {
                    job->pair = 0;
                }

                mess_threadpooljob tpjob;
                mess_threadpooljob_init(&tpjob);
                tpjob->worker = job_decompose_real;
                tpjob->data = job;
                mess_threadpool_insert(m->solverpool, tpjob,&id);


                /* @ref mess_decomp_lureuse_kernelcsr(work->rows, work->cols,   vwork,  work->colptr, work->rowptr,
                        lind, lp,                       // L pattern
                        uind, up,                       // U pattern
                        pperm, qperm,                                   // permuatations
                        m->lvalues[i], m->uvalues[i]);

                mess_free(vwork); */
            }
            /*-----------------------------------------------------------------------------
             *  real system matrix and complex shifts
             *-----------------------------------------------------------------------------*/
            else if ( m->datatypes[i] == MESS_COMPLEX && use_complex == 0){
                mess_int_t id;
                mess_double_cpx_t *vwork;
                mess_try_alloc2(vwork, mess_double_cpx_t *, work->nnz *sizeof(mess_double_cpx_t));
                if (use_gernalized == 0) {
                    if ( havelshift ) {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = shiftsl_complex->values_cpx[i] * work->values[j];
                    } else {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = work->values[j];
                    }
                    for ( j = 0 ; j < work->rows; j++)  vwork[dpos[j]] += shiftsr_complex->values_cpx[i];
                } else {
                    if ( havelshift) {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = shiftsl_complex->values_cpx[i]* work->values [j] + shiftsr_complex->values_cpx[i]*shiftmatrix->values[j];
                    } else {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = work->values [j] + shiftsr_complex->values_cpx[i]*shiftmatrix->values[j];
                    }
                }

                job->vwork = (void *) vwork;
                job->lvalues = (void *) m->lvalues_cpx[i];
                job->uvalues = (void *) m->uvalues_cpx[i];
                job->lvalues_next = (void *) m->lvalues_cpx[i+1];
                job->uvalues_next = (void *) m->uvalues_cpx[i+1];

                if ((cabs(shiftsr_complex->values_cpx[i]-conj(shiftsr_complex->values_cpx[i+1])) < eps ) && havelshift== 0 )  {
                    printf(" " MESS_PRINTF_INT " complex pair. (diff = %lg)\n", i, cabs(shiftsr_complex->values_cpx[i]-conj(shiftsr_complex->values_cpx[i+1]))) ;
                    job->pair = 1;
                    // skip the next parameter
                    i++;
                } else {
                    job->pair = 0;
                }

                // job_decompose_complex(job);
                mess_threadpooljob tpjob;
                mess_threadpooljob_init(&tpjob);
                tpjob->worker = job_decompose_complex;
                tpjob->data = job;
                mess_threadpool_insert(m->solverpool, tpjob,&id);

                /* @ref mess_decomp_lureuse_kernelcsr_complex(work->rows, work->cols,   vwork,  work->colptr, work->rowptr,
                            lind, lp,                       // L pattern
                            uind, up,                       // U pattern
                            pperm, qperm,                                   // permuatations
                            m->lvalues_cpx[i], m->uvalues_cpx[i]);
                mess_free(vwork); */
            }

            /*-----------------------------------------------------------------------------
             *  complex system matrix
             *-----------------------------------------------------------------------------*/
            else {
                mess_int_t id;
                mess_double_cpx_t *vwork;
                mess_try_alloc2(vwork, mess_double_cpx_t *, work->nnz* sizeof(mess_double_cpx_t));
                if (use_gernalized == 0) {
                    if ( havelshift ) {
                        for ( j = 0 ; j < work->nnz  ; j++ ) vwork[j] = shiftsl_complex->values_cpx[i]*work->values_cpx[j];
                    } else {
                        memcpy(vwork, work->values_cpx, work->nnz * sizeof(mess_double_cpx_t));
                    }
                    for ( j = 0 ; j < work->rows; j++)  vwork[dpos[j]] += shiftsr_complex->values_cpx[i];
                } else {
                    if (havelshift) {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = shiftsl_complex->values_cpx[i]*work->values_cpx [j] + shiftsr_complex->values_cpx[i]*shiftmatrix->values_cpx[j];
                    } else {
                        for ( j = 0; j < work->nnz; j++) vwork[j] = work->values_cpx [j] + shiftsr_complex->values_cpx[i]*shiftmatrix->values_cpx[j];
                    }
                }

                job->vwork = (void *) vwork;
                job->lvalues = (void *) m->lvalues_cpx[i];
                job->uvalues = (void *) m->uvalues_cpx[i];
                job->lvalues_next = (void *) m->lvalues_cpx[i+1];
                job->uvalues_next = (void *) m->uvalues_cpx[i+1];

                if ((cabs(shiftsr_complex->values_cpx[i]-conj(shiftsr_complex->values_cpx[i+1])) < eps ) && havelshift == 0 )  {
                    printf(" " MESS_PRINTF_INT " complex pair. (diff = %lg)\n", i, cabs(shiftsr_complex->values_cpx[i]-conj(shiftsr_complex->values_cpx[i+1]))) ;
                    job->pair = 1;
                    // skip the next parameter
                    i++;
                } else {
                    job->pair = 0;
                }

                // job_decompose_complex(job);
                mess_threadpooljob tpjob;
                mess_threadpooljob_init(&tpjob);
                tpjob->worker = job_decompose_complex;
                tpjob->data = job;
                mess_threadpool_insert(m->solverpool, tpjob,&id);


                /* @ref mess_decomp_lureuse_kernelcsr_complex(work->rows, work->cols,   vwork,  work->colptr, work->rowptr,
                            lind, lp,                       // L pattern
                            uind, up,                       // U pattern
                            pperm, qperm,                                   // permuatations
                            m->lvalues_cpx[i], m->uvalues_cpx[i]);
                mess_free(vwork); */

            }
        }

        // avoid freeing the memory for work
        m->workrowptr = work->rowptr;
        m->workcolptr = work->colptr;
        work->rowptr = NULL;
        work->colptr = NULL;


    }
#endif

    mlu->data = (void *) m;
    mlu->cols = matrix->cols;
    mlu->rows = matrix->rows;
    mlu->clear = multilu_clear;
    mlu->solve = multilu_solve;
    mlu->solvem = multilu_solvem;
    mlu->solvet = multilu_solvet;
    mlu->solveh = multilu_solveh;
    mlu->solvemt = multilu_solvemt;
    mlu->solvemh = multilu_solvemh;
    mlu->indx = nshifts;
    mlu->memsize = multilu_memsize;
    mlu->getL = multilu_getL;
    mlu->getU = multilu_getU;
    mlu->getp = multilu_getp;
    mlu->getq = multilu_getq;
    mess_try_alloc(mlu->name, char *, sizeof(char) * 20);
    strcpy(mlu->name, "MLU-mess_lu");


    mess_matrix_clear(&work);

    if (use_gernalized) {
        mess_matrix_clear(&shiftmatrix);
    }
    if ( clear_base == 1){
        mess_direct_clear(&base);
    }
    mess_vector_clear(&shiftsr_real);
    mess_vector_clear(&shiftsr_complex);
    if ( havelshift ) {
        mess_vector_clear(&shiftsl_real);
        mess_vector_clear(&shiftsl_complex);
    }

    if (!use_gernalized) mess_free(dpos);
    return 0;
}

