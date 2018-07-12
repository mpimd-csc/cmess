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
 * @file lib/easyfrontend/lyap.c
 * @brief Easy to use Lyapunov solver interface.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>

/**
 * @brief Frontend to compute a factorized solution of a (generalized) Lyapunov Equation.
 * @param[in] A   input \f$ A \f$ matrix of the Lyapunov Equation
 * @param[in] E   input \f$ E \f$ matrix of the Lyapunov Equation \n
 *                (NULL if \f$ E \f$ does not exist)
 * @param[in] B    input right hand side of the Lyapunov Equation
 * @param[out] Z Cholesky factor \f$ Z \f$ of the solution \f$ X \approx ZZ^T \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lyap function solves either
 * \f[
 *  AX +XA^T +BB^T =0
 * \f]
 * or
 * \f[
 *  AXE^T + EXA^T + BB^T = 0
 * \f]
 * depending on the existence of the matrix \f$ E \f$. \n
 * If the matrix object for \f$ E  \f$ is set to @c NULL the standard Lyapunov Equation is solved, otherwise the
 * generalized one is solved.\n
 * Depending on the storage type of the matrix \f$ A \f$ the solver decides whether a dense solver
 * \ref mess_direct_create_generalized_lyapunovchol or a sparse solver \ref mess_lrcfadi_adi is used. In both cases the solution is returned
 * as a Cholesky factor \f$ X \approx ZZ^T \f$.
 *
 * \attention
 * The \ref mess_lrcfadi_adi uses its default parameters set by \ref mess_options_init. \n
 * If these settings are not applicable to the given problem, the equation have to be solved by configuring the
 * ADI properly and call \ref mess_lrcfadi_adi separately.
 *
 */
int mess_lyap ( mess_matrix A, mess_matrix E, mess_matrix B, mess_matrix Z )
{
    int haveE = 0 ;
    MSG_FNAME(__func__);
    mess_matrix  Etmp= NULL ;
    int ret = 0 ;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_real(A);
    mess_check_real(B);
    mess_check_nullpointer(Z);

    if ( E != NULL ) {
        haveE=1;
        mess_check_same_size(A,E);
    }

    if ( MESS_IS_DENSE(A)) {
        /*-----------------------------------------------------------------------------
         *  solve a dense lyapunov eqn
         *-----------------------------------------------------------------------------*/
        if (haveE) {
            /*-----------------------------------------------------------------------------
             *  with mass matrix
             *-----------------------------------------------------------------------------*/
            if ( MESS_IS_DENSE(E)) {
                Etmp = E;
            } else {
                ret = mess_matrix_init(&Etmp);                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(E,Etmp,MESS_DENSE);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
            }

            mess_direct lyapsol;
            ret = mess_direct_init(&lyapsol);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
            ret = mess_direct_create_generalized_lyapunovchol(A, Etmp, lyapsol);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_generalize_lyapunovchol);
            ret = mess_direct_solvem(MESS_OP_NONE, lyapsol, B, Z);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
            mess_direct_clear(&lyapsol);

            if ( Etmp != E ) {
                mess_matrix_clear(&Etmp);
            }
        } else {
            /*-----------------------------------------------------------------------------
             *  standard case
             *-----------------------------------------------------------------------------*/
            mess_direct lyapsol;
            ret = mess_direct_init(&lyapsol);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
            ret = mess_direct_create_generalized_lyapunovchol(A,NULL,lyapsol);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_lyapunovchol);
            ret = mess_direct_solvem(MESS_OP_NONE, lyapsol, B, Z);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
            mess_direct_clear(&lyapsol);
        }

    } else {

        /*-----------------------------------------------------------------------------
         *  Solve with LRCFADI
         *-----------------------------------------------------------------------------*/
            int ok =0 ;
            mess_options opt;
            mess_status  stat;
            mess_equation     eqn;

            ret = mess_options_init(&opt);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_options_init);
            ret = mess_status_init(&stat);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_status_init);
            ret = mess_equation_init(&eqn);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_init);
            opt->adi_res2_tol = MESS_MIN(1e-13 * sqrt(A->rows), sqrt(mess_eps()));
            opt->adi_res2c_tol =MESS_MIN(1e-14 * sqrt(A->rows), sqrt(mess_eps()));
            opt->adi_rel_change_tol = MESS_MIN(mess_eps() * A->rows, sqrt(mess_eps()));

            ret = mess_equation_lyap(eqn, opt, A, E, B);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_lyap);

            ret = mess_parameter(eqn, opt, stat);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_parameter);
            ret = mess_lrcfadi_adi(eqn, opt, stat, Z);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_lrcfadi_adi);

            ok = stat->stop_res2 || stat->stop_res2c || stat->stop_rel || stat->stop_user;
            if ( !ok ) {
                MSG_ERROR("The ADI did not converge. Please check your matrix or change the solver parameters.\n");
                mess_status_print(stat);
                mess_equation_clear(&eqn);
                mess_status_clear(&stat);
                mess_options_clear(&opt);
                return MESS_ERROR_CONVERGE;
            }
            mess_equation_clear(&eqn);
            mess_status_clear(&stat);
            mess_options_clear(&opt);
    }

    return 0;
}       /* -----  end of function mess_easy_lyap  ----- */




