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
 * @file tests/lrcf_adi/check_lrnm_unsym.c
 * @brief Check the solution of Riccati Equations using LRCF-NM.
 * @test
 * This function checks the @ref mess_lrcfadi_nm function defined in lrnm.c that means it checks if the Riccati Equations
 * \f[ A X + X A^T - X C^T C X + B B^T = 0 \f]
 * and
 * \f[ A^T X + X A - X B B^T X + C^T C = 0 \f]
 * is solved correctly using the low-rank-Cholesky factor Newton method (LRCF-NM).\n
 * These Riccati Equations are generated using the @ref mess_equation_riccati function defined in equation_lyap.c.
 *
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
 * @return \f$ 1 \f$ if \f$ x \f$ is between \f$ 0.4 \f$ and \f$ 0.7 \f$, otherwise 0
 *
 * The @ref fdmcol2 function returns \f$ 1 \f$ if the space variable \f$ x \f$ fullfills
 * \f[ 0.4 < x < 0.7, \f]
 * otherwise \f$ 0 \f$.
 */
double fdmcol2(double x, double y) {
    if ( x > 0.4 && x < 0.7 ) return 1;
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
    // return 400.0*x;
}

int main ( int argc, char **argv){
    mess_init();
    mess_error_level=3;
    mess_version();

    mess_matrix  A,B,C,Z, G;
    mess_options opt;
    mess_status stat1,stat2;
    mess_int_t dim ,paratype, linesearch;
    mess_equation eqn;
    int ret, ok1 = 0,  ok2 = 0;
    double res2, rel, lrnm_res2, res2_diffexp;
    unsigned int nm_gpstep, singleshift;


    /*-----------------------------------------------------------------------------
     * check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 9 ) {
        fprintf(stderr, "usage: %s dim paratype nm_gpstep singleshift direct mdirect linesearch resmeth \n", argv[0]);
        return 1;
    }

    dim         = atoi(argv[1]);
    paratype    = atoi(argv[2]);
    nm_gpstep   = atoi(argv[3]);
    singleshift = atoi(argv[4]);
    linesearch  = atoi(argv[7]);

    /*-----------------------------------------------------------------------------
     *  generate matrix
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&B,&C,&G,&Z);
    CALL(mess_matgen_fdmmatrix(A,dim, fx, fy, NULL));
    CALL(mess_matgen_fdmcolumn(B,dim, fdmcol));
    CALL(mess_matgen_fdmcolumn(G,dim, fdmcol2));
    CALL(mess_matrix_ctranspose(G,C));

    /*-----------------------------------------------------------------------------
     *  setup options
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    opt->type                   = MESS_OP_NONE;
    opt->adi_shifts_paratype    = paratype;
    opt->nm_singleshifts        = singleshift;
    opt->nm_gpStep              = nm_gpstep;
    opt->adi_output             = 1;
    opt->nm_output              = 1;
    opt->nm_linesearch          = linesearch;
    opt->residual_method        = atoi(argv[8]);
    CALL(mess_options_print(opt));

    //select direct solver
    mess_direct_lupackage_t direct = atoi(argv[5]);
    mess_direct_lu_select(direct);

    mess_multidirect_t mdirect = atoi(argv[6]);
    mess_multidirect_select(mdirect);


    /*-----------------------------------------------------------------------------
     *  solve riccati equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_init(&stat1));
    CALL(mess_equation_init(&eqn));
    CALL(mess_equation_std_riccati(eqn, opt, A,B, C));
    CALL(mess_parameter(eqn, opt, stat1));
    CALL(mess_lrcfadi_nm(eqn,opt,stat1,Z));
    CALL(mess_matrix_dynorm2(G,&rel));
    CALL(mess_lrcfadi_residual(eqn, opt, Z, &res2));
    printf("res2 = %lg res2/rel = %lg\n",res2, res2/rel);
    ok1 = stat1->stop_res2 || stat1->stop_rel || stat1->stop_user;
    CALL(mess_options_clear(&opt));
    CALL(mess_equation_clear(&eqn));

    /*-----------------------------------------------------------------------------
     * compare orders of residuals
     *-----------------------------------------------------------------------------*/
    lrnm_res2       = stat1->res2_norm;
    res2_diffexp    = fabs(log10(lrnm_res2)-log10(res2));
    printf("2-Norm absolute residual LRNM                   = %e\n", lrnm_res2);
    printf("2-Norm absolute residual mess_lrcfadi_residual  = %e\n", res2);

    if(res2_diffexp >= 1){
        printf("Difference in exponents                         = %e\n", res2_diffexp);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  setup options
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    opt->type                   = MESS_OP_TRANSPOSE;
    opt->adi_shifts_paratype    = paratype;
    opt->nm_gpStep              = nm_gpstep;
    opt->nm_singleshifts        = singleshift;
    opt->adi_output             = 1;
    opt->nm_output              = 1;
    opt->nm_linesearch          = linesearch;
    opt->residual_method        = atoi(argv[8]);
    CALL(mess_options_print(opt));

    /*-----------------------------------------------------------------------------
     *  solve riccati equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_init(&stat2));
    CALL(mess_equation_init(&eqn));
    CALL(mess_equation_std_riccati(eqn, opt, A,B, C));
    CALL(mess_parameter(eqn, opt, stat2));
    CALL(mess_lrcfadi_nm(eqn,opt,stat2,Z));
    CALL(mess_matrix_dynorm2(G,&rel));
    CALL(mess_lrcfadi_residual(eqn, opt, Z, &res2));
    printf("res2 = %lg res2/rel = %lg\n",res2, res2/rel);
    ok2 = stat2->stop_res2 || stat2->stop_rel || stat2->stop_user;

    /*-----------------------------------------------------------------------------
     * compare orders of residuals
     *-----------------------------------------------------------------------------*/
    lrnm_res2       = stat2->res2_norm;
    res2_diffexp    = fabs(log10(lrnm_res2)-log10(res2));
    printf("2-Norm absolute residual LRNM                   = %e\n", lrnm_res2);
    printf("2-Norm absolute residual mess_lrcfadi_residual  = %e\n", res2);

    if(res2_diffexp >= 1){
        printf("Difference in exponents                         = %e\n", res2_diffexp);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  print status
     *-----------------------------------------------------------------------------*/
    printf("Status - op = MESS_OP_NONE \n");
    CALL(mess_status_print(stat1));
    printf("\nStatus - op = MESS_OP_TRANSPOSE \n");
    CALL(mess_status_print(stat2));

    /*-----------------------------------------------------------------------------
     *  clear memory
     *-----------------------------------------------------------------------------*/
    mess_status_clear(&stat1);
    mess_status_clear(&stat2);
    mess_options_clear(&opt);
    mess_equation_clear(&eqn);
    MESS_CLEAR_MATRICES(&A,&B,&C,&G,&Z);
    mess_exit();
    if ( ok1 && ok2  ) return 0;
    return 1;
}

