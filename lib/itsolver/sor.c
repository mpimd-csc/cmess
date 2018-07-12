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
 * @file lib/itsolver/sor.c
 * @brief Solve \f$ Ax=b \f$ with SOR and similar algorithms.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>

#ifdef _OPENMP_H
#include <omp.h>
#endif

/**
 * @brief Solve a system of linear equations with successive overrelaxation.
 * @param[in] matrix   input matrix \f$ A \f$ of the equation
 * @param[in] pre   input preconditioning placeholder to fit in the API of the other iterative solvers (ignored for this method)
 * @param[in] b     input right hand side
 * @param[in,out] x     on input: intial solution, \n on output: computed solution
 * @param[in] opt    input options for the solver \n (details see \ref mess_solver_options)
 * @param[out] stat statistics about the solution process \n (details see \ref mess_solver_status)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_solver_sor function solves a system of linear equations
 * \f[ Ax=b \f]
 * with successive overrelaxation (SOR). \n
 * The overrelaxation can be set in the omega value in \b opt. If you use the stepdebug feature, the stepdebug
 * function is only evaluated every \f$ 5 \f$th step. The preconditioner is not used.
 *
 */
int mess_solver_sor(mess_matrix matrix, mess_precond pre, mess_vector b, mess_vector x, mess_solver_options opt,
        mess_solver_status stat)
{
    MSG_FNAME(__func__);
    mess_int_t dim;
    int stepdebug = 0;
    double tol = 0;
    double omega = 1.0;
    mess_int_t maxit = 0;
    int conv = -1;
    double s = 0;
    mess_matrix work;
    mess_int_t i = 0;
    mess_int_t j, k;
    double d = 1;
    double resid;
    double normb;
    double t;
    int ret = 0 ;
    double relresid;

    /*-----------------------------------------------------------------------------
     *  check inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_real(matrix);
    mess_check_square(matrix);
    mess_check_nullpointer(x);
    mess_check_nullpointer(b);

    if ( opt == NULL ) {
        MSG_WARN("no options specified, using defaults\n");
        tol = 1e-5;
        maxit = matrix->rows;
        stepdebug = 0;
        omega = 1.0;
    }else{
        tol = opt->tol;
        maxit = opt->maxit;
        if ( opt->stepdebug != NULL ) stepdebug = 1;
        omega =  opt->omega;
    }
    if ( stat == NULL ) {
        MSG_WARN("no statistics avaialable\n");
    }
    if ( x->dim != matrix->rows) {
        MSG_WARN("resize x from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", x->dim, matrix->rows);
        ret = mess_vector_resize(x, matrix->rows);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }
    if ( b->dim != matrix->cols) {
        MSG_ERROR("dimension of b mismatch. b->dim = " MESS_PRINTF_INT " \t matrix->cols = " MESS_PRINTF_INT "\n", b->dim, matrix->cols);
        return MESS_ERROR_DIMENSION;
    }
    ret = mess_vector_toreal(x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
    ret = mess_vector_toreal(b);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);


    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/

    dim = matrix->rows;

    MESS_MATRIX_CHECKFORMAT(matrix, work, conv, MESS_CSR);
    ret = mess_direct_res2(MESS_OP_NONE, work, x, b, &resid, &relresid);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_res2);
    ret = mess_vector_norm2 (b,&normb);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);

    for ( i = 1; i <= maxit; i++ ){
        for ( j = 0; j < dim; j++){
            s = 0;
            for ( k = work->rowptr[j] ; k < work->rowptr[j+1];k++){
                if (work->colptr[k]!= j){
                    s += work->values[k] * x->values [work->colptr[k]];
                } else {
                    d = work->values[k];
                }
            }

            s = (b->values[j]-s)/d;
            t = x->values[j];
            x->values[j] = t + omega*(s-t);
        }
        if ( i % 5 == 0){
            ret = mess_direct_res2(MESS_OP_NONE, work, x, b, &resid, &relresid);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_res2);
            MSG_INFO("i = " MESS_PRINTF_INT "  resid = %lg                           \r", i, resid/normb);

            if (stepdebug ){
                opt->stepdebug(resid, resid/normb, i, opt->aux_stepdebug);
            }
        }
        if ( resid/normb < tol){
            MSG_INFO("converged i = " MESS_PRINTF_INT "\n", i);
            break;
        }
    }
    ret = mess_direct_res2(MESS_OP_NONE, work, x, b, &resid, &relresid);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_res2);
    if (stat) {
        stat->it = i;
        stat->relres =resid/normb;
        stat->res = resid;
        stat->converged = (resid/normb < tol) ? 1: 0;
    }
    if ( conv == 0)
        mess_matrix_clear(&work);

    return 0;
}

/**
 * @brief Solve a system of linear equations with symmetric successive overrelaxation ).
 * @param[in] matrix  input matrix \f$ A \f$ of the equation
 * @param[in] pre   input preconditioning placeholder to fit in the API of the other iterative solvers (ignored for this method)
 * @param[in] b     input right hand side of the equation
 * @param[in,out] x  on input: intial solution,\n on output: computed solution
 * @param[in] opt    input options for the solver \n (details see \ref mess_solver_options)
 * @param[out] stat statistics about the solution process \n (details see \ref mess_solver_status)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_solver_ssor function solves a system of linear equations
 * \f[ Ax=b \f]
 * with symmetric successive overrelaxation (SSOR). \n
 * The overrelaxation can be set in the omega value in \b opt. If you use the stepdebug feature, the stepdebug
 * function is only evaluated every 5th step. The preconditioner is not used.
 *
 */
int mess_solver_ssor(mess_matrix matrix, mess_precond pre, mess_vector b, mess_vector x, mess_solver_options opt,
        mess_solver_status stat)
{
    MSG_FNAME(__func__);
    mess_int_t dim;
    int stepdebug = 0;
    double tol = 0;
    double omega = 1.0;
    mess_int_t maxit = 0;
    int conv = -1;
    double s = 0;
    mess_matrix work;
    mess_int_t i = 0;
    mess_int_t j, k;
    double d = 1;
    double resid;
    double normb;
    double t;
    int ret = 0 ;
    double relresid;

    /*-----------------------------------------------------------------------------
     *  check inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_real(matrix);
    mess_check_square(matrix);
    mess_check_nullpointer(x);
    mess_check_nullpointer(b);

    if ( opt == NULL ) {
        MSG_WARN("no options specified, using defaults\n");
        tol = 1e-5;
        maxit = matrix->rows;
        stepdebug = 0;
        omega = 1.0;
    }else{
        tol = opt->tol;
        maxit = opt->maxit;
        if ( opt->stepdebug != NULL ) stepdebug = 1;
        omega =  opt->omega;
    }
    if ( stat == NULL ) {
        MSG_WARN("no statistics avaialable\n");
    }
    if ( x->dim != matrix->rows) {
        MSG_WARN("resize x from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", x->dim, matrix->rows);
        ret = mess_vector_resize(x, matrix->rows);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }
    if ( b->dim != matrix->cols) {
        MSG_ERROR("dimension of b mismatch. b->dim = " MESS_PRINTF_INT " \t matrix->cols = " MESS_PRINTF_INT "\n", b->dim, matrix->cols);
        return MESS_ERROR_DIMENSION;
    }
    ret = mess_vector_toreal(x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
    ret = mess_vector_toreal(b);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);


    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/

    dim = matrix->rows;

    MESS_MATRIX_CHECKFORMAT(matrix, work, conv, MESS_CSR);
    ret = mess_direct_res2(MESS_OP_NONE, work, x, b, &resid, &relresid);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_res2);
    ret = mess_vector_norm2 (b,&normb);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);

    for ( i = 1; i <= maxit; i++ ){

        /*-----------------------------------------------------------------------------
         *  First Sweep
         *-----------------------------------------------------------------------------*/
        for ( j = 0; j < dim; j++){
            s = 0;
            for ( k = work->rowptr[j] ; k < work->rowptr[j+1];k++){
                if (work->colptr[k]< j){
                    s += work->values[k] * x->values[work->colptr[k]];
                } else if ( work->colptr[k] > j) {
                    s += work->values[k] * x->values[work->colptr[k]];
                } else {
                    d = work->values[k];
                }
            }

            s = (b->values[j]-s)/d;
            t = x->values[j];
            x->values[j] = t + omega*(s-t);
        }

        /*-----------------------------------------------------------------------------
         *  Second Sweep
         *-----------------------------------------------------------------------------*/
        for (j=dim-1; j >=0; j--) {
            s = 0;
            for ( k = work->rowptr[j] ; k < work->rowptr[j+1];k++){
                if (work->colptr[k] < j){
                    s += work->values[k] * x->values[work->colptr[k]];
                } else if (work->colptr[k] > j) {
                    s += work->values[k] * x->values[work->colptr[k]];
                } else {
                    d = work->values[k];
                }
            }

            s = (b->values[j]-s)/d;
            t = x->values[j];
            x->values[j] = t + omega*(s-t);
        }
        /*-----------------------------------------------------------------------------
         *  Finalize the step
         *-----------------------------------------------------------------------------*/
        if ( i % 5 == 0){
            ret = mess_direct_res2(MESS_OP_NONE, work, x, b, &resid, &relresid);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_res2);
            MSG_INFO("i = " MESS_PRINTF_INT "  resid = %lg                           \r", i, resid/normb);

            if (stepdebug ){
                opt->stepdebug(resid, resid/normb, i, opt->aux_stepdebug);
            }
        }
        if ( resid/normb < tol){
            MSG_INFO("converged i = " MESS_PRINTF_INT "\n", i);
            break;
        }
    }
    ret = mess_direct_res2(MESS_OP_NONE, work, x, b, &resid, &relresid);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_res2);
    if (stat) {
        stat->it = i;
        stat->relres =resid/normb;
        stat->res = resid;
        stat->converged = (resid/normb < tol) ? 1: 0;
    }
    if ( conv == 0)
        mess_matrix_clear(&work);

    return 0;
}

/**
 * @brief Solve a system of linear equations with the Gauss-Seidel method.
 * @param[in] matrix   input matrix \f$ A \f$ of the equation
 * @param[in] pre   input preconditioning placeholder to fit in the API of the other iterative solvers (ignored for this method)
 * @param[in] b     input right hand side
 * @param[in,out] x     on input: intial solution, \n
 *                       on output: computed solution
 * @param[in] opt    input options for the solver \n (details see \ref mess_solver_options)
 * @param[out] stat statistics about the solution process \n (details see \ref mess_solver_status)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_solver_gaussseidel function solves
 * \f[ Ax=b\f]
 * with Gauss-Seidel method. \n
 * It is only a wrapper around the \ref mess_solver_sor function setting the overrelaxation parameter \b omega to 1.
 *
 * \sa mess_solver_sor
 */
int mess_solver_gaussseidel(mess_matrix matrix, mess_precond pre, mess_vector b, mess_vector x, mess_solver_options opt,
        mess_solver_status stat)
{
    MSG_FNAME(__func__);
    double omega_old = 0 ;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_real(matrix);
    mess_check_square(matrix);
    mess_check_nullpointer(x);
    mess_check_nullpointer(b);
    mess_check_nullpointer(opt);

    omega_old = opt->omega;
    ret = mess_solver_sor(matrix, pre, b, x, opt, stat);
    opt->omega = omega_old;
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_sor);
    return 0;
}

