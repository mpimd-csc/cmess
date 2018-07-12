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
 * @file lib/dynsys/h2/h2_norm.c
 * @brief Compute the \f$ \mathcal{H}_2 \f$-norm of a LTI system.
 * @author @koehlerm
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"

#ifndef INFINITY
 #define INFINITY 1.0/0.0
#endif

/**
 * @brief Compute the \f$ \mathcal{H}_2 \f$-norm of a LTI System given as matrices.
 * @param[in] A      input system matrix
 * @param[in] B     input matrix
 * @param[in] C      input output matrix
 * @param[in] E  input mass matrix ( @c NULL if not need )
 * @param[out] norm pointer to the double value where the norm is put in
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_h2_norm_internal function computes the \f$ \mathcal{H}_2 \f$-norm of a LTI system using Lyapunov Equations,
 * i.e. the \f$ \mathcal{H}_2 \f$ norm is computed via
 * \f[ \Vert \Sigma \Vert_{ \mathcal{H}_2 } = \sqrt{ trace\left( C P C^T \right)}, \f]
 * where \f$ \Sigma \f$ denotes the given LTI system and \f$ P \f$ is the solution of the Lyapunov Equation
 * \f[ APE^T + EPA^T + BB^T = 0.  \f]
 * This function is internally used by \ref mess_h2_norm. \n
 * More details about the norm computation are available in \cite morAnt05 .
 *
 */
int mess_h2_norm_internal ( mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix E,  double *norm)
{
    MSG_FNAME(__func__);
    mess_int_t n, p,mc;
    double tr = 0;
    mess_options opt = NULL;
    mess_status stat = NULL;
    mess_matrix P = NULL, X=NULL, X2=NULL;
    double h2n = 0;
    int ret =0;
    unsigned short skiprest = 0;



    /*-----------------------------------------------------------------------------
     *  check inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);
    mess_check_nullpointer(norm);

    if ( E != NULL) {
        mess_check_square(E);
        if ( E->rows!=A->rows){
            MSG_ERROR("E must have the same dimension as A\n");
            return MESS_ERROR_DIMENSION;
        }
    }
    n=A->rows;
    p=C->rows;
    mc=C->cols;

    if ( B->rows != n) {
        MSG_ERROR("Input matrix B has the wrong number of rows\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( p > mc || C->cols !=n) {
        MSG_ERROR("Output matrix C has the wrong dimension.\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_options_init(&opt);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_options_init);
    opt->adi_res2_tol = mess_eps() * n * 100;
    opt->adi_res2c_tol = mess_eps() *sqrt(n);
    opt->adi_rel_change_tol = mess_eps()*sqrt(n);
    opt->adi_shifts_arp_p = 50;
    opt->adi_shifts_arp_m = 25;
    opt->adi_shifts_l0 = 16;
    opt->adi_maxit = 250;

    if (mess_error_level > 2) {
        opt->adi_output = 1;
        mess_options_print(opt);
    } else {
        opt->adi_output = 0;
    }
    ret = mess_status_init(&stat);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_status_init);

    ret = mess_matrix_init(&P);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

    /*-----------------------------------------------------------------------------
     *  solve the Lyapunov Equation AP+PA^T+BB^T = 0
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_DENSE(A)){
        MSG_WARN("using dense Lyapunov solver\n");
        mess_direct lyap;
        ret  = mess_direct_init(&lyap);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
        ret = mess_direct_create_generalized_lyapunovchol(A,E, lyap);   FUNCTION_FAILURE_HANDLE(ret,!( ret==0 || ret == MESS_ERROR_EIGENVALUES), mess_direct_create_generalized_lyapunovchol);
        if ( ret == MESS_ERROR_EIGENVALUES) {
            h2n=INFINITY;
            skiprest = 1;
        } else {
            ret  = mess_direct_solvem(MESS_OP_NONE,lyap, B, P);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
        }
        mess_direct_clear(&lyap);
    } else {
        MSG_INFO("using ADI-Lyapunov solver\n");
        mess_equation     eqn;
        int ok = 0;

        ret = mess_equation_init(&eqn);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_init);
        ret = mess_equation_lyap(eqn, opt, A, E, B);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_lyap);
        ret = mess_parameter(eqn, opt, stat);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_parameter);
        ret = mess_lrcfadi_adi(eqn, opt, stat, P);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_lrcfadi_adi);
        mess_equation_clear(&eqn);

        ok = stat->stop_res2 || stat->stop_res2c || stat->stop_rel || stat->stop_user;

        if ( !ok  ){
            h2n=INFINITY;
            skiprest = 1;
        }
        /*      else if ( stat -> unstable && !stat->stop_res2 ) {
            MSG_WARN("It seems that the system is unstable. Set H2-Norm to infinity.\n");
            h2n = INFINITY;
            skiprest = 1;
        }
          else if ( ! (stat->stop_rel || stat->stop_res2 || stat->stop_res2c)) {
            h2n=INFINITY;
            skiprest = 1;
            // mess_status_print(stat);
            // MSG_ERROR("The LRCF-ADI do not converge. restart with complex parameters.\n");
            |+ mess_matrix_clear(&P);
            mess_status_clear(&stat);
            mess_options_clear(&opt);
            *norm=1/0.0;
            return MESS_ERROR_CONVERGE; +|
        } */
    }
    if (skiprest == 0) {
        ret = mess_matrix_init(&X);                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&X2);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_multiply(MESS_OP_NONE, C, MESS_OP_NONE, P, X);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, X, MESS_OP_HERMITIAN, X, X2);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

        if ( MESS_IS_REAL(X2)) {
            ret = mess_matrix_trace(X2,&tr);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_trace);
            h2n =sqrt(tr);
        } else {
            mess_double_cpx_t tr2 = 0;
            ret = mess_matrix_tracec(X2,&tr2);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_tracec);
            h2n =sqrt(tr2);
        }
        mess_matrix_clear(&X);
        mess_matrix_clear(&X2);

    }
    *norm = h2n;

    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
    mess_matrix_clear(&P);
    mess_status_clear(&stat);

    mess_options_clear(&opt);

    return 0;
}


/**
 * @brief Compute the \f$ \mathcal{H}_2 \f$-norm of a LTI system.
 * @param[in] lti  input given LTI system
 * @param[out] norm pointer to the norm
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_h2_norm function computes the \f$ \mathcal{H}_2 \f$-norm of a LTI system
 * \f[\begin{array}{ccccc}
 *  E & \dot{x} &=& A x& + Bu \\
 *    &   y     &=& C x&
 *    \end{array}
 *  \f]
 * via
 * \f[ \Vert \Sigma \Vert_{ \mathcal{H}_2 } = \sqrt{ trace\left( C P C^T \right)}, \f]
 * where \f$ \Sigma \f$ denotes the given LTI system and \f$ P \f$ is the solution of the Lyapunov Equation
 * \f[ APE^T + EPA^T + BB^T = 0.  \f]
 * This function is only a wrapper around \ref mess_h2_norm_internal to work with the \ref mess_dynsys
 * object.
 *
 */
int  mess_h2_norm ( mess_dynsys lti, double *norm )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(lti);
    mess_check_nullpointer(norm);
    if (!(MESS_IS_DYNSYS_LTI(lti) || MESS_IS_DYNSYS_GLTI(lti))){
        MSG_ERROR("The type of the dynamic system is wrong.");
        return MESS_ERROR_DYNSYS;
    }
    ret = mess_h2_norm_internal(lti->A, lti->B, lti->C, lti->E, norm);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_h2_norm_internal);
    return 0;
}       /* -----  end of function mess_h2_norm  ----- */

