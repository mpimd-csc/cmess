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
 * @file lib/itsolver/bicgstab.c
 * @brief Solve \f$ Ax=b \f$ with BiCGStab.
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
 * @brief Solve a system of linear equations with BICGSTAB algorithm.
 * @param[in] matrix input   maxtrix \f$ A \f$ of the equation
 * @param[in] pre input     left preconditioner for the equation \n (if not used, set to \c NULL)
 * @param[in] b input       right hand side
 * @param[in,out] x   on input: initial guess, <br >on output: approximated solution
 * @param[in] opt input object containing the parameters for the iteration, \n if it is set to \c NULL defaults are used
 * @param[out] stat output statitics about the convergence \n (if not used, set to \c NULL)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_solver_bicgstab function implements the BICGSTAB algorithm to solve
 * \f[ Ax = b . \f]
 * It can be accelerated using a left preconditioner defined by \c pre. It works for real and complex data. \n
 * Details about it are available in \cite Mei11.
 *
 */
int mess_solver_bicgstab(   mess_mvpcall matrix,
        mess_precond pre,
        mess_vector b,
        mess_vector x,
        mess_solver_options opt,
        mess_solver_status  stat)
{
    MSG_FNAME(__func__);
    double resid, normb;
    double rho_1 = 0, rho_2 = 0, alpha = 0, beta = 0, omega = 0, nrm=.0 ;
    mess_double_cpx_t crho_1=0, crho_2=0,calpha=0, cbeta=0,comega=0;
    mess_vector  p, phat, s, shat, t, v, r, rtilde, tmp;
    int precond = 0;
    int stepdebug = 0;
    double tol;
    mess_int_t maxit;
    mess_int_t dim;
    mess_int_t i;
    int cancel = 0;
    int ret = 0;
    mess_datatype_t data_type;


    /*-----------------------------------------------------------------------------
     *  check inputs
     *-----------------------------------------------------------------------------*/
    // mess_check_nullpointer(matrix);
    // mess_check_square(matrix);
    // mess_check_real(matrix);
    mess_check_nullpointer(x);
    mess_check_nullpointer(b);

    if ( pre != NULL ){
        // MSG_INFO("preconditioner used\n");
        precond = 1;
    }
    if ( opt == NULL ) {
        MSG_WARN("no options specified, using defaults\n");
        tol = 1e-7;
        maxit = sqrt(matrix->dim) * 5 ;
        stepdebug = 0;
    }else{
        tol = opt->tol;
        maxit = opt->maxit;
        if ( opt->stepdebug != NULL ) stepdebug = 1;
    }
    if ( stat == NULL ) {
        MSG_WARN("no statistics available\n");
    }

    if ( x->dim != matrix->dim) {
        MSG_WARN("resize x from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", x->dim, matrix->dim);
        ret = mess_vector_resize(x, matrix->dim); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }
    if ( b->dim != matrix->dim) {
        MSG_ERROR("dimension of b mismatch. b->dim = " MESS_PRINTF_INT " \t matrix->cols = " MESS_PRINTF_INT "\n", b->dim, matrix->dim);
        return MESS_ERROR_DIMENSION;
    }


    /*-----------------------------------------------------------------------------
     *  Datatype selection
     *-----------------------------------------------------------------------------*/
    if (matrix->data_type == MESS_COMPLEX || MESS_IS_COMPLEX(b) ){
        ret = mess_vector_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_tocomplex(b); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        // MSG_INFO("Using complex BICGSTAB\n");
        data_type = MESS_COMPLEX;
    } else {
        ret = mess_vector_toreal(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
        ret = mess_vector_toreal(b); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
        data_type = MESS_REAL;
    }


    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/

    dim = matrix->dim;
    MESS_INIT_VECTORS(&p,&phat,&s,&shat,&t,&v,&r,&rtilde,&tmp);
    ret = mess_vector_alloc(p, dim, data_type);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(phat, dim, data_type);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(s, dim, data_type);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(shat, dim, data_type);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(t, dim, data_type);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(v, dim, data_type);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(r, dim, data_type);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(rtilde, dim, data_type);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(tmp, dim, data_type);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);


    ret = mess_vector_norm2(b,&normb);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    // MSG_INFO("Normb: %lg\n", normb );
    // Vector r = b - A * x;
    ret = mess_vector_copy(b, r);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
    ret = mess_mvpcall_apply(matrix, MESS_OP_NONE, x, tmp); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), matrix.mvp);
    if ( stat) stat->num_mvp = 1;
    ret = mess_vector_axpy(-1.0, tmp, r);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // Vector rtilde = r;
    ret = mess_vector_copy ( r, rtilde);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);

    if (normb == 0.0) normb = 1;
    ret =mess_vector_norm2(r,&nrm); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    if (fabs(resid =( nrm / normb)) <= tol) {
        if ( stat ) {
            stat->converged = 1;
            stat->it = 0;
            stat->relres = resid;
            ret = mess_vector_norm2(r,&(stat->res));       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
        }

        // clean up
        mess_vector_clear(&p);
        mess_vector_clear(&phat);
        mess_vector_clear(&s);
        mess_vector_clear(&shat);
        mess_vector_clear(&t);
        mess_vector_clear(&v);
        mess_vector_clear(&r);
        mess_vector_clear(&rtilde);
        mess_vector_clear(&tmp);
        return 0;
    }

    for ( i = 1; i <= maxit; i++) {
        if ( data_type == MESS_COMPLEX ) {mess_vector_dotc(rtilde,r,&crho_1);}
        else { mess_vector_dot(rtilde,r,&(rho_1)); crho_1=rho_1; }

        // MSG_INFO(" i = %4ld \t rho1  = %15lg + %15lgi\n", i, creal(crho_1), cimag(crho_1));
        if (crho_1 == 0) {
            cancel = 1;
            break;
        }
        if (i == 1){
            ret = mess_vector_copy(r, p); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
        } else {
            if ( data_type == MESS_COMPLEX) {
                cbeta = (crho_1/crho_2) * (calpha/comega);
                // p = r + beta(0) * (p - omega(0) * v);
                ret  = mess_vector_copy(p, tmp); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
                mess_vector_axpyc(-comega, v, tmp); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
                mess_vector_copy(r, p); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
                mess_vector_axpyc(cbeta, tmp, p); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
            } else {
                beta = (rho_1/rho_2) * (alpha/omega); cbeta = beta;
                // p = r + beta(0) * (p - omega(0) * v);
                ret  = mess_vector_copy(p, tmp); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
                mess_vector_axpy(-omega, v, tmp); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
                mess_vector_copy(r, p); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
                mess_vector_axpy(beta, tmp, p); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
            }
            // MSG_INFO(" i = %4ld \t beta  = %15lg + %15lgi\n", i, creal(cbeta), cimag(cbeta));

        }
        if ( precond ){
            // apply preconditioner
            ret =  pre->solve(pre, opt, p, phat); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), pre->solve);
            ret = mess_mvpcall_apply(matrix, MESS_OP_NONE, phat, v); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_apply );
            if ( stat ) stat->num_mvp++;
        }else {
            // without preconditioner
            ret = mess_mvpcall_apply(matrix,MESS_OP_NONE, p, v); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_apply );
            if ( stat ) stat->num_mvp++;
        }
        if ( data_type == MESS_COMPLEX) {
            mess_double_cpx_t rtv;
            ret = mess_vector_dotc(rtilde, v,&rtv); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_dotc);
            // MSG_INFO(" i = %4ld \t rtv  = %15lg + %15lgi\n", i, creal(rtv), cimag(rtv));

            calpha = crho_1 / rtv;
            //  s = r - alpha(0) * v;
            ret = mess_vector_copy(r, s);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy);
            ret = mess_vector_axpyc(-calpha, v, s); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_axpyc);
        } else {
            //alpha = rho_1 / mess_vector_dot(rtilde, v);
            ret = mess_vector_dot(rtilde,v,&alpha);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_dot);
            alpha = rho_1 / alpha;
            //  s = r - alpha(0) * v;
            ret = mess_vector_copy(r, s);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy);
            ret = mess_vector_axpy(-alpha, v, s);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_axpy);
            calpha= alpha;
        }
        // MSG_INFO(" i = %4ld \t alpha = %15lg + %15lgi\n", i, creal(calpha), cimag(calpha));
        ret =mess_vector_norm2(s,&nrm); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
        if ((resid = (nrm/normb)) < tol) {
            if ( precond ) {
                if ( data_type == MESS_COMPLEX) mess_vector_axpyc(calpha, phat, x);
                else mess_vector_axpy(alpha, phat, x);
            } else {
                if ( data_type == MESS_COMPLEX) mess_vector_axpyc(calpha, p, x);
                else  mess_vector_axpy(alpha, p, x);
            }
            cancel = 2;
            break;
        }
        if ( precond ) {
            pre->solve(pre,opt, s, shat);
            ret = mess_mvpcall_apply(matrix, MESS_OP_NONE, shat, t); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_apply );
            if ( stat ) stat->num_mvp++;
            // shat = M.solve(s);
            // t = A * shat;
        } else {
            ret = mess_mvpcall_apply(matrix, MESS_OP_NONE , s, t); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_apply );
            if ( stat ) stat->num_mvp++;
        }


        if ( data_type==MESS_COMPLEX) {
            mess_double_cpx_t temp;
            mess_vector_dotc(t,s,&temp);
            mess_vector_dotc(t,t,&comega);
            comega = temp/comega;
            //comega = mess_vector_dotc(t,s) / mess_vector_dotc(t,t);
        }
        else {
            double temp1,temp2;
            ret = mess_vector_dot(t,s,&temp1);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_dot);
            ret = mess_vector_dot(t,t,&temp2);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_dot);
            omega = temp1/temp2;
            comega = omega;
            //omega = mess_vector_dot(t,s) / mess_vector_dot(t,t); comega = omega;
        }
        // MSG_INFO(" i = %4ld \t omega = %15lg + %15lgi\n", i, creal(comega), cimag(comega));

        if (precond) {
            //x += alpha(0) * phat + omega(0) * shat;
            if (data_type == MESS_COMPLEX) {
                mess_vector_axpyc(calpha, phat, x);
                mess_vector_axpyc(comega, shat, x);
            } else {
                mess_vector_axpy(alpha, phat, x);
                mess_vector_axpy(omega, shat, x);
            }
        } else {
            if ( data_type == MESS_COMPLEX) {
                mess_vector_axpyc(calpha, p, x);
                mess_vector_axpyc(comega, s, x);
            } else {
                mess_vector_axpy(alpha, p, x);
                mess_vector_axpy(omega, s, x);
            }
        }
        mess_vector_copy(s, r);
        if ( data_type == MESS_COMPLEX) { mess_vector_axpyc(-comega, t, r); }
        else { mess_vector_axpy(-omega, t, r); }

        rho_2 = rho_1;
        crho_2 = crho_1;
        ret =mess_vector_norm2(r,&nrm); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
        if ((resid = (nrm / normb)) < tol) {
            // MSG_WARN("resid = %lg\n", resid);
            cancel = 3;
            break;
        }
        if (comega == 0) {
            cancel = 4;
            break;
        }
        MSG_INFO("i = " MESS_PRINTF_INT "  resid = %lg                           \r", i, resid);
        if (stepdebug ){
            opt->stepdebug(resid * normb, resid, i, opt->aux_stepdebug);
        }
        // MSG_INFO("\n");
    }
    // MSG_INFO("\n");

    switch (cancel){
        case 1:     // rho_1 == 0
            ret = mess_vector_norm2(r,&nrm); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
            tol = nrm / normb;
            if (stat) {
                stat->converged = 0;
                stat->relres = tol;
                stat->it = i;
                ret = mess_vector_norm2(r,&(stat->res));       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
            }
            // MSG_WARN("cancel because rho_1 == 0\n");
            break;
        case 2:     // norm(s) < tolb
            if ( stat ) {
                stat->converged = 1;
                stat->relres = resid;
                stat->it = i;
                ret = mess_vector_norm2(r,&(stat->res));       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
            }
            // MSG_INFO("converged i = " MESS_PRINTF_INT " res=%lg, case 2 \n", i, resid);
            break;
        case 3:     // norm(r) < tolb
            if (stat) {
                stat->converged = 1;
                stat->relres = resid;
                stat->it = i;
                ret = mess_vector_norm2(r,&(stat->res));       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
            }
            // MSG_INFO("\nconverged i = " MESS_PRINTF_INT " res=%lg, case 3 \n", i, resid);            break;
        case 4:     // omega == 0
            ret = mess_vector_norm2(r,&nrm); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
            tol = nrm / normb;
            if ( stat ) {
                stat->converged = 0;
                stat->relres = tol;
                stat->it = i;
                ret = mess_vector_norm2(r,&(stat->res));       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
            }
            // MSG_INFO("converged i = " MESS_PRINTF_INT " res=%lg, case 4 \n", i, resid);

            break;
        default:
            if ( stat ) {
                stat->converged = (resid < tol) ? 1 : 0;
                stat->relres = resid;
                ret = mess_vector_norm2(r,&(stat->res));       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
                stat->it = maxit;
            }
            break;
    }


    mess_vector_clear(&p);
    mess_vector_clear(&phat);
    mess_vector_clear(&s);
    mess_vector_clear(&shat);
    mess_vector_clear(&t);
    mess_vector_clear(&v);
    mess_vector_clear(&r);
    mess_vector_clear(&rtilde);
    mess_vector_clear(&tmp);
    return 0;
}



