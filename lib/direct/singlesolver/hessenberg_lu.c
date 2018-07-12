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
 * @file lib/direct/singlesolver/hessenberg_lu.c
 * @brief Generate a dense @lapack based solver for upper Hessenberg Matrices.
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



/**
 * @internal
 * @brief Internal structure for the @lapack based solver.
 *
 * @attention Internal use only.
 */
typedef struct hessenberg_solver {
    double *val;                    /**< the real lu decomposed matrix */
    mess_double_cpx_t *val_cpx;     /**< the complex lu decomposed matrix */
    mess_int_t *ipiv;               /**< row permutation */
    mess_int_t n;                   /**< dimension of the system */
    mess_datatype_t data_type;      /**< flag real/complex */
} hessenberg_solver;



/**
 * @internal
 * @brief Clear the internal structure @ref hessenberg_solver.
 * @param[in,out] data  pointer to the internal data
 * @return always zero
 *
 * The @ref hessenberg_clear function clears a @lapack solver.
 *
 * @attention Internal use only.
 */
static int hessenberg_clear(void *data){
    hessenberg_solver *sol = (hessenberg_solver *) data;
    if (sol->val)          mess_free(sol->val);
    if (sol->val_cpx)      mess_free(sol->val_cpx);
    mess_free(sol->ipiv);
    mess_free(sol);
    return 0;
}


/**
 * @internal
 * @brief Solve a linear system using the @ref hessenberg_solver
 * @param[in] data          input pointer to the internal data
 * @param[in] op            input determine which equation has to be solved
 * @param[in] nrhs          input number of right hand sides
 * @param[in,out] rhs       input/output points to the real right hand side
 * @param[in,out] rhs_cpx   input/output points to the complex right hand side
 * @param[in] ld_rhs        input leading dimension of the right hand side array
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref hessenberg_solvex function solves a linear system \f$ Ax =b\f$ using the
 * @ref hessenberg_solver. The right hand side data array is @p rhs or @p rhs_cpx.
 * If the @ref hessenberg_solver is @ref MESS_REAL @p rhs must be given and @p rhs_cpx must points to @c NULL or vice versa.
 * The solution is at exit written to @p rhs / @p rhs_cpx respectively.
 *
 * @attention Internal use only.
 *
 */
static int hessenberg_solvex(void *data, mess_operation_t op, mess_int_t nrhs, double *rhs, mess_double_cpx_t *rhs_cpx, mess_int_t ld_rhs){
    MSG_FNAME(__func__);
    hessenberg_solver *sol = (hessenberg_solver*) data;
    mess_int_t  i = 0, ip = 0;
    double one_d = 1.0, temp;
    mess_double_cpx_t  one_c = 1.0, temp_c;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_real_or_complex(sol);
    mess_check_operation_type(op);
    mess_check_positive(nrhs);

    if(!rhs && !rhs_cpx){
        MSG_ERROR("At least rhs or rhs_cpx must be given.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if(rhs && rhs_cpx){
        MSG_ERROR("Either rhs or rhs_cpx must be given.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if(rhs && MESS_IS_COMPLEX(sol)){
        MSG_ERROR("Complex solver and real data are given.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if(rhs_cpx && MESS_IS_REAL(sol)){
        MSG_ERROR("Real solver and complex data are given.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    mess_check_positive(ld_rhs);

    /*-----------------------------------------------------------------------------
     *  solve equation
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(sol)){

        /*-----------------------------------------------------------------------------
         *  real case
         *-----------------------------------------------------------------------------*/
        if(op==MESS_OP_NONE){
            for(i=0;i<(sol->n-1);++i){
                ip = sol->ipiv[i];
                if(ip != i){
                    F77_GLOBAL(dswap,DSWAP)(&nrhs, rhs+ip, &ld_rhs, rhs+i, &ld_rhs);
                }
                temp = -sol->val[i+1+i*sol->n];
                F77_GLOBAL(daxpy,DAXPY)(&nrhs, &temp, rhs+i, &ld_rhs, rhs+i+1, &ld_rhs);
            }
            F77_GLOBAL(dtrsm,DTRSM)("L", "U", "N", "N", &(sol->n), &nrhs, &one_d, sol->val, &(sol->n), rhs, &ld_rhs);
        }else{
            // MESS_OP_TRANSPOSE, MESS_OP_HERMITIAN is equal for real solver
            F77_GLOBAL(dtrsm,DTRSM)("L","U","T","N", &(sol->n), &nrhs, &one_d, sol->val, &(sol->n), rhs, &ld_rhs);

            for(i=(sol->n-2); i>=0; --i){
                temp = -sol->val[i+1+i*sol->n];
                F77_GLOBAL(daxpy,DAXPY)(&nrhs, &temp, rhs+i+1, &ld_rhs, rhs+i, &ld_rhs);
                ip = sol->ipiv[i];
                if(ip != i){
                    F77_GLOBAL(dswap,DSWAP)(&nrhs, rhs+ip, &ld_rhs, rhs+i, &ld_rhs);
                }
            }
        }
    }else{

        /*-----------------------------------------------------------------------------
         *  complex case
         *-----------------------------------------------------------------------------*/
        if(op==MESS_OP_NONE){
            for(i=0;i<(sol->n-1);++i){
                ip = sol->ipiv[i];
                if(ip != i){
                    F77_GLOBAL(zswap,ZSWAP)(&nrhs, rhs_cpx+ip, &ld_rhs, rhs_cpx+i, &ld_rhs);
                }
                temp_c = -sol->val_cpx[i+1+i*sol->n];
                F77_GLOBAL(zaxpy,ZAXPY)(&nrhs, &temp_c, rhs_cpx+i, &ld_rhs, rhs_cpx+i+1, &ld_rhs);
            }
            F77_GLOBAL(ztrsm,ZTRSM)("L","U","N","N", &(sol->n), &nrhs, &one_c, sol->val_cpx, &(sol->n), rhs_cpx, &ld_rhs);
        }else if(op == MESS_OP_TRANSPOSE){
            F77_GLOBAL(ztrsm,ZTRSM)("L","U","T","N", &(sol->n), &nrhs, &one_c, sol->val_cpx, &(sol->n), rhs_cpx, &ld_rhs);

            for(i=(sol->n-2); i>=0; --i){
                temp_c = -sol->val_cpx[i+1+i*sol->n];
                F77_GLOBAL(zaxpy,ZAXPY)(&nrhs, &temp_c, rhs_cpx+i+1, &ld_rhs, rhs_cpx+i, &ld_rhs);
                ip = sol->ipiv[i];
                if(ip != i){
                    F77_GLOBAL(zswap,ZSWAP)(&nrhs, rhs_cpx+ip, &ld_rhs, rhs_cpx+i, &ld_rhs);
                }
            }
        }else{
            F77_GLOBAL(ztrsm,ZTRSM)("L","U","C","N", &(sol->n), &nrhs, &one_c, sol->val_cpx, &(sol->n), rhs_cpx, &ld_rhs);
            for(i=(sol->n-2); i>=0; --i){
                temp_c = -conj(sol->val_cpx[i+1+i*sol->n]);
                F77_GLOBAL(zaxpy,ZAXPY)(&nrhs, &temp_c, rhs_cpx+i+1, &ld_rhs, rhs_cpx+i, &ld_rhs);
                ip = sol->ipiv[i];
                if(ip != i){
                    F77_GLOBAL(zswap,ZSWAP)(&nrhs, rhs_cpx+ip, &ld_rhs, rhs_cpx+i, &ld_rhs);
                }
            }
        }
    }

    return 0;
}


/**
 * @internal
 * @brief Solve a linear system using the @ref hessenberg_solver
 * @param[in] data          input pointer to the internal data
 * @param[in] op            input determine which equation has to be solved
 * @param[in] b             input right hand side
 * @param[out] x            output the solution
 * @return zero on success or a non-zero error value otherwise
 *
 * @attention Internal use only.
 *
 */
static int hessenberg_solve_all(void * data, mess_operation_t op, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t nrhs=1;
    hessenberg_solver *sol = (hessenberg_solver*) data;

    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(x);
    mess_check_nullpointer(b);
    mess_check_operation_type(op);
    if(b->dim != sol->n){
        MSG_ERROR("Number of rows of right hand side does not fit to solver.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  hessenberg_solvex call
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_copy(b,x);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);

    if(MESS_IS_COMPLEX(sol) && MESS_IS_REAL(x)){
        ret = mess_vector_tocomplex(x);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tcocomplex);
    }

    if(MESS_IS_COMPLEX(sol)){
        ret = hessenberg_solvex(data, op, nrhs, NULL, x->values_cpx, x->dim);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), hessenberg_solvex);
    }else{
        if(MESS_IS_REAL(x)){
            ret = hessenberg_solvex(data, op, nrhs, x->values, NULL, x->dim);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), hessenberg_solvex);
        }else{
            mess_vector xr, xi;
            MESS_INIT_VECTORS(&xr, &xi);
            ret = mess_vector_realpart(x, xr);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_realpart);
            ret = mess_vector_imagpart(x, xi);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_imagpart);

            //solve for real and imaginary part
            ret = hessenberg_solvex(data, op, nrhs, xr->values, NULL, xr->dim);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), hessenberg_solvex);
            ret = hessenberg_solvex(data, op, nrhs, xi->values, NULL, xi->dim);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), hessenberg_solvex);

            ret = mess_vector_complex_from_parts(xr, xi, x);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_complex_from_parts);

            MESS_CLEAR_VECTORS(&xr, &xi);
        }
    }
    return ret;
}

/**
 * @internal
 * @brief Solve a linear system using the @ref hessenberg_solver
 * @param[in] data          input pointer to the internal data
 * @param[in] op            input determine which equation has to be solved
 * @param[in] b             input right hand side
 * @param[out] x            output the solution
 * @return zero on success or a non-zero error value otherwise
 *
 * @attention Internal use only.
 *
 */
static int hessenberg_solvem_all(void * data, mess_operation_t op, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    int ret = 0;
    hessenberg_solver *sol = (hessenberg_solver*) data;

    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(x);
    mess_check_nullpointer(b);
    mess_check_operation_type(op);
    if(b->rows != sol->n){
        MSG_ERROR("Number of rows of right hand side does not fit to solver\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  hessenberg_solvex call
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_copy(b,x);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);

    if(MESS_IS_COMPLEX(sol) && MESS_IS_REAL(x)){
        ret = mess_matrix_tocomplex(x);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
    }

    if(MESS_IS_COMPLEX(sol)){
        ret = hessenberg_solvex(data, op, x->cols, NULL, x->values_cpx, x->ld);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), hessenberg_solvex);
    }else{
        if(MESS_IS_REAL(x)){
            ret = hessenberg_solvex(data, op, x->cols, x->values, NULL, x->ld);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), hessenberg_solvex);
        }else{
            mess_matrix xr, xi;
            MESS_INIT_MATRICES(&xr, &xi);
            ret = mess_matrix_realpart(x, xr);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realpart);
            ret = mess_matrix_imagpart(x, xi);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_imagpart);

            //solve for real and imaginary part
            ret = hessenberg_solvex(data, op, xr->cols, xr->values, NULL, xr->ld);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), hessenberg_solvex);
            ret = hessenberg_solvex(data, op, xi->cols, xi->values, NULL, xi->ld);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), hessenberg_solvex);

            ret = mess_matrix_complex_from_parts(xr, xi, x);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_complex_from_parts);

            MESS_CLEAR_MATRICES(&xr, &xi);
        }
    }
    return ret;
}



static int hessenberg_solve(void * data, mess_vector b, mess_vector x){
    return hessenberg_solve_all(data, MESS_OP_NONE, b, x);
}

static int hessenberg_solvet(void * data, mess_vector b, mess_vector x){
    return hessenberg_solve_all(data, MESS_OP_TRANSPOSE, b, x);
}

static int hessenberg_solveh(void * data, mess_vector b, mess_vector x){
    return hessenberg_solve_all(data, MESS_OP_HERMITIAN, b, x);
}


static int hessenberg_solvem(void * data, mess_matrix b, mess_matrix x){
    return hessenberg_solvem_all(data, MESS_OP_NONE, b, x);
}

static int hessenberg_solvemt(void * data, mess_matrix b, mess_matrix x){
    return hessenberg_solvem_all(data, MESS_OP_TRANSPOSE, b, x);
}

static int hessenberg_solvemh(void * data, mess_matrix b, mess_matrix x){
    return hessenberg_solvem_all(data, MESS_OP_HERMITIAN, b, x);
}




/**
 * @brief Generate a dense @lapack based LU-solver for upper hessenberg matrices.
 * @param[in] matrix    input upper hessenberg matrix to decompose
 * @param[out] solver   output generated solver
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_direct_create_hessenberg_lu function generates a LU solver for dense upper hessenberg matrices based on @lapack.
 * It is not checked that @p matrix is actually upper Hessenberg,
 *
 * See @cite GolV13 Chapter 4.3.4
 *
 */
int mess_direct_create_hessenberg_lu(mess_matrix matrix, mess_direct solver){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i=0, ip=0, temp;
    double temp_d;
    mess_double_cpx_t temp_c;
    hessenberg_solver *sol;

    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(solver);
    mess_check_nullpointer(matrix);
    mess_check_square(matrix);
    mess_check_real_or_complex(matrix);
    mess_check_dense(matrix);

    /*-----------------------------------------------------------------------------
     *  create lu factorization
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(sol, hessenberg_solver *, sizeof(hessenberg_solver));
    memset(sol, 0, sizeof(hessenberg_solver));
    sol->n = matrix->rows;
    mess_try_alloc(sol->ipiv, mess_int_t *, sizeof(mess_int_t)*(sol->n));
    sol->data_type = matrix->data_type;

    if(MESS_IS_REAL(matrix)){
        mess_try_alloc(sol->val, double *, sizeof(double)*sol->n*sol->n);
        F77_GLOBAL(dlacpy,DLACPY)("A", &(matrix->rows), &(matrix->cols), matrix->values, &(matrix->ld), sol->val, &(sol->n));
    }else{
        mess_try_alloc(sol->val_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*sol->n*sol->n);
        F77_GLOBAL(zlacpy,ZLACPY)("A", &(matrix->rows), &(matrix->cols), matrix->values_cpx, &(matrix->ld), sol->val_cpx, &(sol->n));
    }

    /*-----------------------------------------------------------------------------
     *  compute lu decomposition using algorithm see reference
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(sol)){
        //find largest element in absolute value for pivoting
        for(i=0;i<sol->n;++i){
            ip = i;
            if(i<(sol->n-1)){
                if(fabs(sol->val[i+1+i*(sol->n)]) >= fabs(sol->val[i+i*(sol->n)])){
                    ip = i+1;
                }
            }
            sol->ipiv[i] = ip;

            if(fabs(sol->val[ip+i*sol->n]) != 0){
                //swap rows
                if(ip != i){
                    temp = sol->n-i;
                    F77_GLOBAL(dswap,DSWAP)(&temp, sol->val+(i+i*sol->n), &(sol->n), sol->val+(ip+i*sol->n), &(sol->n));
                }

                if(i<(sol->n-1)){
                    sol->val[i+1+i*sol->n] /= sol->val[i+i*sol->n];
                }

            }else{
                MSG_ERROR("Zero pivot element during upper Hessenberg LU-factorization.\n");
                goto error;
                ret = MESS_ERROR_ARGUMENTS;
            }

            // update trailing matrix
            if(i<(sol->n-1)){
                temp = sol->n-i-1;
                temp_d = -sol->val[i+1+i*sol->n];
                F77_GLOBAL(daxpy,DAXPY)(&temp, &temp_d, sol->val+(i+(i+1)*sol->n), &(sol->n), sol->val+(i+1+(i+1)*sol->n), &(sol->n));
            }
        }
    }else{
        //find largest element in absolute value for pivoting
        for(i=0;i<sol->n;++i){
            ip = i;
            if(i<(sol->n-1)){
                if(cabs(sol->val_cpx[i+1+i*(sol->n)]) >= cabs(sol->val_cpx[i+i*(sol->n)])){
                    ip = i+1;
                }
            }
            sol->ipiv[i] = ip;

            if(cabs(sol->val_cpx[ip+i*sol->n]) != 0){
                //swap rows
                if(ip != i){
                    temp = sol->n-i;
                    F77_GLOBAL(zswap,ZSWAP)(&temp, sol->val_cpx+(i+i*sol->n), &(sol->n), sol->val_cpx+(ip+i*sol->n), &(sol->n));
                }

                if(i<(sol->n-1)){
                    sol->val_cpx[i+1+i*sol->n] /= sol->val_cpx[i+i*sol->n];
                }

            }else{
                MSG_ERROR("Zero pivot element during upper Hessenberg LU-factorization.\n");
                goto error;
                ret = MESS_ERROR_ARGUMENTS;
            }

            // update trailing matrix
            if(i<(sol->n-1)){
                temp = sol->n-i-1;
                temp_c = -sol->val_cpx[i+1+i*sol->n];
                F77_GLOBAL(zaxpy,ZAXPY)(&temp, &temp_c, sol->val_cpx+(i+(i+1)*sol->n), &(sol->n), sol->val_cpx+(i+1+(i+1)*sol->n), &(sol->n));
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  setup solver
     *-----------------------------------------------------------------------------*/
    solver->data        = (void *) sol;
    solver->clear       = hessenberg_clear;
    solver->solve       = hessenberg_solve;
    solver->solvet      = hessenberg_solvet;
    solver->solveh      = hessenberg_solveh;
    solver->solvem      = hessenberg_solvem;
    solver->solvemt     = hessenberg_solvemt;
    solver->solvemh     = hessenberg_solvemh;
    solver->inverse     = NULL;
    solver->det         = NULL;
    solver->detc        = NULL;
    solver->getL        = NULL;
    solver->getU        = NULL;
    solver->getpermp    = NULL;
    solver->getpermq    = NULL;
    solver->rows        = matrix->rows;
    solver->cols        = matrix->cols;
    solver->data_type   = matrix->data_type;
    SET_SOLVERNAME(solver->name, __func__);

    return 0;

    /*-----------------------------------------------------------------------------
     *  error handling
     *-----------------------------------------------------------------------------*/
error:
    hessenberg_clear(sol);
    return ret;
}

