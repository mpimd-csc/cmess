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
 * @addtogroup test_lrcfadi
 * @{
 *
 * @file tests/lrcf_adi/check_lradi_unsym.c
 * @brief Check the solution of standard Lyapunov Equations using the LRCF-ADI method.
 * @test
 * This function checks the @ref mess_lrcfadi_adi function defined in lradi.c that means it checks if a standard Lyapunov Equation
 * \f[ op(F) X +  X op(F)^T + op(G) op(G)^T = 0 \f]
 * is solved correctly using the low-rank Cholesky factor alternate direction implicit (LRCF-ADI) method.\n
 * Operation applied to matrices \f$ F \f$ and \f$ G \f$ can be
 * * \ref MESS_OP_NONE - \f$ op(C) = C \f$
 * * \ref MESS_OP_TRANSPOSE - \f$ op(C)= C^T \f$
 * This Lyapunov Equation is generated using the @ref mess_equation_lyap function defined in equation_lyap.c.
 * @}
 */



#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "../call_macro.h"

/**
 * @brief Determine values used in a FDM matrix.
 * @param[in] x      input space variable \f$ x \f$
 * @param[in] y          input space variable \f$ y \f$
 * @return 1 if \f$ x \f$ is between \f$ 0.1 \f$ and \f$ 0.4 \f$, otherwise 0
 *
 * The @ref fdmcol function returns \f$ 1 \f$ if the space variable \f$ x \f$ fullfills
 * \f[ 0.1 < x < 0.4, \f]
 * otherwise \f$ 0 \f$.
 */
double fdmcol(double x, double y) {
    if ( x > 0.1 && x < 0.4 ) return 1;
    return 0;
}

/**
 * @brief Determine values used in a FDM matrix.
 * @param[in] x      input space variable \f$ x \f$
 * @param[in] y          input space variable \f$ y \f$
 * @return \f$ 400x \f$
 *
 * The @ref fx function always returns \f$ 400x \f$.
 */
double fx(double x, double y ) {
    return 400.0*x;
}

/**
 * @brief Determine values used in a FDM matrix.
 * @param[in] x      input space variable \f$ x \f$
 * @param[in] y          input space variable \f$ y \f$
 * @return \f$ 10 \sin(y) \f$
 *
 * The @ref fy function always returns \f$ 10 \sin(y)\f$.
 */
double fy(double x, double y){
    return 10.0*sin(y);
}

int main ( int argc, char **argv){
    mess_version();
    mess_init();
    mess_error_level = 3;

    mess_matrix F, G, Z, G2;
    mess_options opt;
    mess_status stat;
    mess_int_t dim, paratype;
    mess_equation eqn;
    int ret, ok1 = 0, ok2 = 0;
    double res2, rel, adi_res2, res2_diffexp;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 6 ) {
        fprintf(stderr, "usage: %s dim paratype direct mdirect resmeth\n", argv[0]);
        return 1;
    }

    dim         = atoi(argv[1]);
    paratype    = atoi(argv[2]);

    /*-----------------------------------------------------------------------------
     *  generate matrix
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&F,&G,&Z,&G2);
    CALL(mess_matgen_fdmmatrix(F,dim, fx, fy, NULL));
    CALL(mess_matgen_fdmcolumn(G,dim, fdmcol));

    /*-----------------------------------------------------------------------------
     *  Setup options
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    opt->adi_shifts_paratype    = paratype;
    opt->adi_output             = 1;
    opt->type                   = MESS_OP_NONE;
    opt->residual_method        = atoi(argv[5]);
    CALL(mess_options_print(opt));

    //select direct solver
    mess_direct_lupackage_t direct = atoi(argv[3]);
    mess_direct_lu_select(direct);

    mess_multidirect_t mdirect = atoi(argv[4]);
    mess_multidirect_select(mdirect);


    /*-----------------------------------------------------------------------------
     *  Setup and Solve the Equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_init(&stat));
    CALL(mess_equation_init(&eqn));
    CALL(mess_equation_std_lyap(eqn, opt, F, G));
    CALL(mess_parameter(eqn,opt,stat));
    CALL(mess_lrcfadi_adi(eqn,opt,stat,Z));
    CALL(mess_matrix_dynorm2(G,&rel));
    CALL(mess_lrcfadi_residual(eqn, opt, Z, &res2));
    printf("res2 = %lg\n",res2/rel );
    CALL(mess_status_print(stat));
    ok1 = stat->stop_res2 || stat->stop_res2c || stat->stop_rel || stat->stop_user;


    /*-----------------------------------------------------------------------------
     * compare orders of residuals
     *-----------------------------------------------------------------------------*/
    adi_res2        = stat->res2_norm;
    res2_diffexp    = fabs(log10(adi_res2)-log10(res2));
    printf("2-Norm absolute residual ADI                    = %e\n", adi_res2);
    printf("2-Norm absolute residual mess_lrcfadi_residual  = %e\n", res2);

    if(res2_diffexp >= 1){
        printf("Difference in exponents                         = %e\n", res2_diffexp);
        return 1;
    }

    mess_matrix_clear(&Z);
    mess_status_clear(&stat);
    mess_options_clear(&opt);
    mess_equation_clear(&eqn);

    /*-----------------------------------------------------------------------------
     *  Setup options
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    opt->adi_shifts_paratype    = paratype;
    opt->adi_output             = 1;
    opt->type                   = MESS_OP_TRANSPOSE;
    opt->residual_method        = atoi(argv[5]);
    CALL(mess_options_print(opt));

    /*-----------------------------------------------------------------------------
     *  check if error occurs by wrong data
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_init(&stat));
    CALL(mess_equation_init(&eqn));

    if ( mess_equation_std_lyap(eqn, opt, F, G) == 0 ) {
        fprintf(stderr, "mess_equation_lyap does not detect wrong sized RHS.\n");
        return 1;
    }
    printf("Wrong size of right hand side detected correctly.\n");
    CALL(mess_matrix_ctranspose(G,G2));

    /*-----------------------------------------------------------------------------
     *  setup, solve and clear equation
     *-----------------------------------------------------------------------------*/
    CALL( mess_matrix_init(&Z));
    CALL(mess_equation_std_lyap(eqn, opt, F, G2));
    CALL(mess_parameter(eqn,opt,stat));
    CALL(mess_lrcfadi_adi(eqn,opt,stat,Z));
    CALL(mess_matrix_dynorm2(G2,&rel));
    CALL(mess_lrcfadi_residual(eqn, opt, Z, &res2));
    printf("res2: %lg\n", res2/rel);
    CALL(mess_status_print(stat));
    ok2 = stat->stop_res2 || stat->stop_res2c || stat->stop_rel || stat->stop_user;

    /*-----------------------------------------------------------------------------
     * compare orders of residuals
     *-----------------------------------------------------------------------------*/
    adi_res2        = stat->res2_norm;
    res2_diffexp    = fabs(log10(adi_res2)-log10(res2));
    printf("2-Norm absolute residual ADI                    = %e\n", adi_res2);
    printf("2-Norm absolute residual mess_lrcfadi_residual  = %e\n", res2);

    if(res2_diffexp >= 1){
        printf("Difference in exponents                         = %e\n", res2_diffexp);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  clear data and finisch
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&F,&G,&G2,&Z);
    mess_status_clear(&stat);
    mess_options_clear(&opt);
    mess_equation_clear(&eqn);

    if ( ok1 && ok2  ) return 0;
    return 1;
}

