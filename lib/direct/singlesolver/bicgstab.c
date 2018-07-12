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
 * @file lib/direct/singlesolver/bicgstab.c
 * @brief Wrapper to use bicgstab as single solver.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"

#define ILU_LEVEL 3
#ifdef _OPENMP_H
#include <omp.h>
#endif

struct bicgstab_solver {
    mess_matrix A;
    mess_precond pre;
    mess_solver_options opt;
    mess_solver_status stat;
    mess_int_t *p;
    mess_int_t *q;
};




static int bicgstab_clear(void *data){
    struct bicgstab_solver *sol = (struct bicgstab_solver *) data;
    mess_matrix_clear(&(sol->A));
    mess_precond_clear(&(sol->pre));
    mess_solver_options_clear(&(sol->opt));
    mess_solver_status_clear(&(sol->stat));
    mess_free(sol->p);
    mess_free(sol->q);
    mess_free(sol);
    return 0;
}

static int bicgstab_solve(void * data, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct bicgstab_solver * sol = (struct bicgstab_solver*) data;
    mess_mvpcall mvpcall;
    int ret = 0;
    mess_vector bint;
    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);


    /*-----------------------------------------------------------------------------
     *  bicg stab can handle complex right hand sides.
     *-----------------------------------------------------------------------------*/
    // ret = mess_vector_toreal_nowarn(x);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_vector_init(&bint);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(bint, b->dim, b->data_type);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_perm(b, sol->p,bint);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_perm_inplace);
    ret = mess_vector_zeros(x);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_zeros);
    ret = mess_mvpcall_matrix(&mvpcall, MESS_OP_NONE, sol->A);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_matrix );
    ret = mess_solver_bicgstab(mvpcall, sol->pre, bint, x, sol->opt, sol->stat);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_bicgstab);
    ret = mess_vector_iperm_inplace(x, sol->q);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_iperm_inplace);
    mess_vector_clear(&bint);
    mess_mvpcall_clear(&mvpcall);
    return 0 ;
}

static int bicgstab_solvem (void * data, mess_matrix B, mess_matrix X){
    MSG_FNAME(__func__);
    struct bicgstab_solver * sol = (struct bicgstab_solver*) data;
    mess_mvpcall mvpcall;
    int ret = 0;
    mess_int_t i;
    mess_vector bint,b,x;
    mess_check_nullpointer(sol);
    mess_check_nullpointer(B);
    mess_check_nullpointer(X);
    mess_check_real_or_complex(B);
    ret = mess_matrix_alloc(X, B->rows, B->cols, 0, MESS_DENSE, B->data_type);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    MESS_INIT_VECTORS(&b,&x);
    ret = mess_vector_alloc(b, B->rows, B->data_type);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(x, B->rows, B->data_type);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_mvpcall_matrix(&mvpcall, MESS_OP_NONE, sol->A);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_matrix );


    /*-----------------------------------------------------------------------------
     *  bicg stab can handle complex right hand sides.
     *-----------------------------------------------------------------------------*/
    for ( i = 0 ; i < B->cols; i++) {
        ret = mess_matrix_getcol(B,i,b);                                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
        ret = mess_vector_init(&bint);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(bint, b->dim, b->data_type);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_perm(b, sol->p,bint);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_perm_inplace);
        ret = mess_vector_zeros(x);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_zeros);
        ret = mess_solver_bicgstab(mvpcall, sol->pre, bint, x, sol->opt, sol->stat);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_bicgstab);
        ret = mess_vector_iperm_inplace(x, sol->q);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_iperm_inplace);
        mess_vector_clear(&bint);
        ret = mess_matrix_setcol(X, i, x);                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);

    }
    mess_mvpcall_clear(&mvpcall);
    mess_vector_clear(&b);
    mess_vector_clear(&x);
    return 0 ;



}


/**
 * @brief Generate a bicgstab solver in order to use it as direct solver.
 * @param[in] matrix     input System matrix
 * @param[out] solver   generated solver based on BiCGSTAB
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_direct_create_bicgstab function generates a pseudo direct solver. \n
 * It relies internally on BiCGSTAB and wraps it in a way that it can be used like a direct solver.
 * It reorders the matrix using the AMD reordering and computes a \f$ ILU(3) \f$ preconditioner with AMD reordering as default.
 * The maximum number of iterations is at least \f$ 100 \f$ in real case or
 * \f$ 300 \f$ in complex case. \n
 *
 */
int mess_direct_create_bicgstab(mess_matrix matrix, mess_direct solver){
    MSG_FNAME(__func__);
    int ret = 0;
    struct bicgstab_solver * sol ;
    mess_int_t maxit;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_square(matrix);

    mess_try_alloc(sol, struct bicgstab_solver *, sizeof(struct bicgstab_solver));
    ret = mess_matrix_init(&(sol->A)); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_solver_options_init(&(sol->opt)); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_options_init);
    ret = mess_solver_status_init(&(sol->stat)); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_status_init);
    ret = mess_matrix_convert(matrix, sol->A, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);

    /*-----------------------------------------------------------------------------
     *  compute the reordering
     *-----------------------------------------------------------------------------*/
    if (1) {
        sol->p = NULL;
        sol->q =NULL;
        mess_try_alloc(sol->p, mess_int_t *, sizeof(mess_int_t)*(sol->A->rows));
        mess_try_alloc(sol->q, mess_int_t *, sizeof(mess_int_t)*(sol->A->rows));
        ret = mess_matrix_sort(sol->A); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);

        // select reordering
#ifdef MESS_HAVE_AMD
        ret = mess_matrix_reorder( MESS_REORDER_AMD, sol->A, sol->p, sol->q);
#else
        ret = mess_matrix_reorder( MESS_REORDER_RCM, sol->A, sol->p, sol->q);
#endif

        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_reorder);
    } else {
        sol->p = NULL;
        sol->q = NULL;
    }

    ret = mess_matrix_perm(sol->A, sol->p, sol->q);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_perm);
    ret = mess_matrix_sort(sol->A);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);

    ret = mess_precond_init(&(sol->pre));                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_precond_init);
    ret = mess_precond_iluk(sol->pre, sol->A, ILU_LEVEL);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_precond_iluk);

    maxit = MESS_MAX(5000,log(log(matrix->rows)) * sqrt(matrix->rows)) ;
    if ( maxit < 300 ) maxit = 300;

    solver->data_type = sol->A->data_type;
    sol->opt->maxit =  maxit;
    sol->opt->tol   = MESS_MAX(1e-13, mess_eps()*matrix->rows);

    if ( mess_error_level > 2) {
        printf("bicgstab solver options:\n");
        mess_solver_options_print(sol->opt);
    }
    // sol->pre = NULL;

    solver->data = (void *) sol;
    solver->clear = bicgstab_clear;
    solver->solve = bicgstab_solve;
    solver->solvet = NULL;
    solver->solvem = bicgstab_solvem;
    solver->solvemt = NULL;
    solver->getL = NULL;
    solver->getU = NULL;
    solver->getpermp = NULL;
    solver->getpermq = NULL;
    solver->rows = matrix->rows;
    solver->cols = matrix->cols;
    SET_SOLVERNAME(solver->name, __func__);

    return ret;
}

