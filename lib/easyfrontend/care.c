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
 * @file lib/easyfrontend/care.c
 * @brief Easy to use Riccati solver interface.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


/**
 * @brief Frontend to compute a factorized solution of a (generalized) Riccati Equation.
 * @param[in] A  input matrix of the Riccati Equation \f$ A \f$
 * @param[in] E  input matrix of the Riccati Equation \f$ E \f$ \n
 *                (@c NULL if \f$ E \f$ does not exist)
 * @param[in] B   input matrix of the Riccati Equation  \f$ B \f$
 * @param[in] C   input matrix of the Riccati Equation  \f$ C \f$
 * @param[out] Z Cholesky factor \f$ Z \f$ of the solution \f$ X \approx ZZ^T \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_care function solves either
 * \f[
 *  A^T X + XA - XBB^TX + C^TC = 0
 * \f]
 * or
 * \f[
 *  A^T X E+E^T X A - E^T XBB^TX E + C^TC = 0
 * \f]
 * depending on the existence of the matrix \f$ E \f$. \n
 * If the matrix object for \f$ E  \f$ is set to @c NULL the standard Riccati Equation is solved, otherwise the
 * generalized one is solved.\n
 * Depending on the storage type of the matrix \f$ A \f$ the solver decides whether a dense solver
 * \ref mess_direct_care or a sparse solver \ref mess_lrcfadi_nm is used. In both cases the solution is returned
 * as a Cholesky factor \f$ X \approx ZZ^T \f$.
 *
 * \attention
 *
 * The \ref mess_lrcfadi_nm uses its default parameters set by \ref mess_options_init. \n
 * If these settings are not applicable to the given problem, the equation have to be solved by configuring the
 * ADI properly and call \ref mess_lrcfadi_adi separately.
 *
 */
int mess_care ( mess_matrix A, mess_matrix E, mess_matrix B, mess_matrix C, mess_matrix Z )
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
    mess_check_nullpointer(C);
    mess_check_real(A);
    mess_check_real(B);
    mess_check_real(C);
    mess_check_nullpointer(Z);

    mess_check_same_rows(A,B);
    mess_check_same_cols(A,C);

    if ( E != NULL ) {
        haveE=1;
        mess_check_same_size(A,E);
    }


    if ( MESS_IS_DENSE(A)) {
        mess_matrix ZZ;
        mess_matrix G;
        mess_matrix Q;

        ret = mess_matrix_init(&ZZ); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init );

        /*-----------------------------------------------------------------------------
         *  Solve a dense Riccati Equation
         *-----------------------------------------------------------------------------*/
        if (haveE) {
            /*-----------------------------------------------------------------------------
             *  with mass matrix
             *-----------------------------------------------------------------------------*/
            if ( MESS_IS_DENSE(E)) {
                Etmp = E;
            } else {
                ret = mess_matrix_init(&Etmp);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(E,Etmp,MESS_DENSE);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
            }
            ret = mess_matrix_init(&G);                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
            ret = mess_matrix_init(&Q);                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
            ret = mess_matrix_multiply(MESS_OP_TRANSPOSE, C, MESS_OP_NONE, C, Q);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_NONE, B, MESS_OP_TRANSPOSE, B, G);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

            ret = mess_direct_care(A, Etmp, Q, G, ZZ);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_care);
            ret = mess_direct_cholfactor(ZZ, Z);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_cholfactor);

            mess_matrix_clear(&G);
            mess_matrix_clear(&Q);

            if ( Etmp != E ) {
                mess_matrix_clear(&Etmp);
            }
        } else {
            /*-----------------------------------------------------------------------------
             *  standard case
             *-----------------------------------------------------------------------------*/
            ret = mess_matrix_init(&G);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
            ret = mess_matrix_init(&Q);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
            ret = mess_matrix_multiply(MESS_OP_TRANSPOSE, C, MESS_OP_NONE, C, Q);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_NONE, B, MESS_OP_TRANSPOSE, B, G);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

            ret = mess_direct_care(A, NULL, Q, G, ZZ);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_care);
            ret = mess_direct_cholfactor(ZZ, Z);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_cholfactor);
            mess_matrix_clear(&G);
            mess_matrix_clear(&Q);


        }
        mess_matrix_clear(&ZZ);

    } else {

        /*-----------------------------------------------------------------------------
         *  Solve with LRCFNM
         *-----------------------------------------------------------------------------*/

        int ok =0 ;
        mess_options opt;
        mess_status  stat;
        mess_equation     eqn;

        ret = mess_options_init(&opt);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_options_init);
        ret = mess_status_init(&stat);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_status_init);
        ret = mess_equation_init(&eqn);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_init);

        opt->adi_res2_tol =  MESS_MIN(1e-13 * sqrt(A->rows), sqrt(mess_eps()));
        opt->adi_res2c_tol = MESS_MIN(1e-14 * sqrt(A->rows), sqrt(mess_eps()));
        opt->adi_rel_change_tol = MESS_MIN(mess_eps() * A->rows, sqrt(mess_eps()));


        opt->nm_res2_tol = MESS_MIN(1e-12 * sqrt(A->rows), sqrt(mess_eps())/10.0);
        opt->type = MESS_OP_TRANSPOSE;
        opt->nm_gpStep = 0;
        opt->nm_maxit  = 30;

        ret = mess_equation_riccati(eqn, opt, A, E, B, C);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_riccati);

        ret = mess_lrcfadi_nm(eqn, opt, stat, Z);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_lrcfadi_nm);

        ok = stat->stop_res2 || stat->stop_res2c || stat->stop_rel || stat->stop_user;
        if ( !ok ) {
            MSG_ERROR("The ADI-NM did not converge. Please check your matrix or change the solver parameters.\n");
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
}




