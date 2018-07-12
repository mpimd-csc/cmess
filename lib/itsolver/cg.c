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
 * @file lib/itsolver/cg.c
 * @brief Conjugate gradients itsolver.
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


/**
 * @brief Solve a system of linear equations using conjugate gradients.
 * @param[in] matrix input      matrix \f$ A \f$ of the equation
 * @param[in] pre input     preconditioner
 * @param[in] b input            right hand side
 * @param[in,out] x     solution
 * @param[in] opt input options
 * @param[out] stat output status
 * @return zero on success or a non-zero error value otherwise

 * The @ref mess_solver_cg function solves a system of linear equations
 * \f[ Ax=b \f]
 * using the conjugate gradients (CG) method. Therefore the matrix \f$ A \f$ need to be symmetric positive definite
 * otherwise the method will fail. \n
 * Additionally a preconditioner can be used to accelerate the convergence. The preconditioner
 * must not destroy symmetry and positive definiteness. A list of avialable preconditioners is listed in
 * \ref itsolver_pre.
 *
 */
int mess_solver_cg  (   mess_mvpcall matrix,
        mess_precond pre,
        mess_vector b,
        mess_vector x,
        mess_solver_options opt,
        mess_solver_status stat
        )
{
    MSG_FNAME(__func__);
    int precond = 0;
    int stepdebug = 0;
    double tol = 0;
    mess_int_t maxit = 0;
    int stats = 1;
    mess_int_t i;
    mess_int_t dim;
    double normb;
    double normr;
    double alpha, alpha2;
    double resid = 0;
    double vr;
    double lambda;
    mess_vector r, tmp, p, v, z;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(x);
    mess_check_nullpointer(b);
    mess_check_real(matrix);
    mess_check_positive(matrix->dim);

    if ( pre != NULL ){
        MSG_INFO("preconditioner used\n");
        precond = 1;
    }
    if ( opt == NULL ) {
        MSG_WARN("no options specified, using defaults\n");
        tol = 1e-5;
        maxit = matrix->dim;
        stepdebug = 0;
    }else{
        tol = opt->tol;
        maxit = opt->maxit;
        if ( opt->stepdebug != NULL ) stepdebug = 1;
    }
    if ( stat == NULL ) {
        stats = 0;
    }
    if ( x->dim != matrix->dim) {
        MSG_WARN("resize x from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", x->dim, matrix->dim);
        ret = mess_vector_resize(x, matrix->dim);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }
    if ( b->dim != matrix->dim) {
        MSG_ERROR("dimension of b mismatch. b->dim = " MESS_PRINTF_INT " \t matrix->cols = " MESS_PRINTF_INT "\n", b->dim, matrix->dim);
        return MESS_ERROR_DIMENSION;
    }
    ret = mess_vector_toreal(x);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
    ret = mess_vector_toreal(b);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
    if (stats) stat->num_mvp = 0;

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    dim = matrix->dim;
    MESS_INIT_VECTORS(&r,&tmp,&p,&v,&z);
    ret = mess_vector_alloc(r, dim, MESS_REAL);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(tmp, dim, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(p, dim, MESS_REAL);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(v, dim, MESS_REAL);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(z, dim, MESS_REAL);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_norm2(b,&normb);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_norm2);
    if ( normb == 0) normb = 1;

    // r = b-Ax
    mess_vector_copy(b, r);
    // mess_matrix_mvp(MESS_OP_NONE,matrix, x, tmp);
    ret =  mess_mvpcall_apply(matrix, MESS_OP_NONE,  x, tmp);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_apply);
    if (stats) stat->num_mvp++;


    mess_vector_axpy(-1.0, tmp, r);
    ret = mess_vector_norm2(r,&normr);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);

    if ( precond ) {
        pre->solve(pre,  opt, r, p);
        mess_vector_dot(r, p,&alpha);
    }else {
        mess_vector_copy(r, p);
        alpha = normr * normr;
    }

    for ( i = 1; i <= maxit; i++){
        if ( (resid = normr/normb) < tol){
            MSG_INFO("converged after " MESS_PRINTF_INT " iterations\n", i);
            break;
        }
        // v = A*p;
        // mess_matrix_mvp(MESS_OP_NONE, matrix, p, v);
        ret = mess_mvpcall_apply(matrix, MESS_OP_NONE, p,v); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_apply);
        if (stats) stat->num_mvp++;

        //vr = mess_vector_dot(v, p);
        mess_vector_dot(v,p,&vr);
        if ( vr <= 0.0){
            MSG_ERROR("matrix isn't positive definite\n");
            break;
        }
        lambda = alpha  / vr;
        mess_vector_axpy(lambda, p, x);
        mess_vector_axpy(-lambda, v, r);
        alpha2 = alpha;
        ret = mess_vector_norm2(r,&normr);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
        if (precond){
            pre->solve ( pre, opt, r, z);
            //alpha = mess_vector_dot(r, z);
            mess_vector_dot(r,z,&alpha);
            mess_vector_scale(alpha/alpha2, p);
            mess_vector_axpy(1.0, z, p);
        }else {
            alpha = normr*normr;
            mess_vector_scale(alpha/alpha2, p);
            mess_vector_axpy(1.0, r, p);
        }
        MSG_INFO("i = " MESS_PRINTF_INT "  resid = %lg                           \r", i, resid);
        if (stepdebug ){
            opt->stepdebug(resid * normb, resid, i, opt->aux_stepdebug);
        }
    }
    MSG_INFO("\n");

    if (stats){
        stat->converged = (resid  < tol) ? 1 : 0;
        stat->it = i;
        stat->relres = resid;
        ret = mess_vector_norm2(r,&(stat->res));       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    }

    mess_vector_clear(&r);
    mess_vector_clear(&v);
    mess_vector_clear(&tmp);
    mess_vector_clear(&p);
    mess_vector_clear(&z);
    return 0;
}

