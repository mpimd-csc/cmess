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
 * @file tests/lrcf_adi/check_lrnm_g.c
 * @brief Check the solution of a generalized Riccati Equation in transposed and nontransposed version using LRCF-NM.
 * @test
 * This function checks the @ref mess_lrcfadi_nm function defined in lrnm.c that means it checks if the Riccati Equations
 * \f[ AXE^T+EXA^T-EXC^TCXE^T+BB^T = 0 \f]
 * and
 * \f[ A^TXE+E^TXA-E^TXBB^TXE+C^TC = 0 \f]
 * are solved correctly using the low-rank-Cholesky factor Newton method (LRCF-NM).\n
 * Riccati Equations are constructed using the @ref mess_equation_griccati function defined in equation_glyap.c.
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


int main ( int argc, char **argv){

    mess_version();
    mess_init();
    mess_error_level = 3;
    mess_matrix A, E, B, C, Z;
    mess_options opt;
    mess_status stat1, stat2;
    mess_int_t paratype;
    mess_equation eqn;
    int ret, ok1 = 0, ok2 = 0, nm_gpstep, singleshifts, linesearch;
    double res2, rel, lrnm_res2, res2_diffexp;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 12 ) {
        fprintf(stderr, "usage: %s A.mtx E.mtx B.mtx C.mtx  paratype  nm_gpstep singleshift direct mdirect linesearch resmeth\n", argv[0]);
        return 1;
    }

    MESS_INIT_MATRICES(&A, &E, &B, &C, &Z);
    CALL(mess_matrix_read_formated(argv[1], A, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[2], E, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[3], B, MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[4], C, MESS_DENSE));

    paratype     = atoi(argv[5]);
    nm_gpstep    = atoi(argv[6]);
    singleshifts = atoi(argv[7]);

    mess_direct_lupackage_t direct = atoi(argv[8]);
    mess_direct_lu_select(direct);

    mess_multidirect_t mdirect = atoi(argv[9]);
    mess_multidirect_select(mdirect);

    linesearch  = atoi(argv[10]);

    /*-----------------------------------------------------------------------------
     *  Setup Options for MESS_OP_NONE
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    opt->adi_shifts_paratype    = paratype;
    opt->nm_gpStep              = nm_gpstep;
    opt->nm_singleshifts        = singleshifts;
    opt->nm_maxit               = 25;
    opt->adi_output             = 0;
    opt->nm_output              = 1;
    opt->nm_res2_tol            = 1e-7;
    opt->type                   = MESS_OP_NONE;
    opt->nm_linesearch          = linesearch;
    opt->residual_method        = atoi(argv[11]);
    CALL(mess_options_print(opt));

    /*-----------------------------------------------------------------------------
     *  Setup, Solve and Clear Riccati Equation for MESS_OP_NONE
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_init(&stat1));
    CALL(mess_equation_init(&eqn));
    CALL(mess_equation_griccati(eqn, opt, A, E, B, C));
    CALL(mess_parameter(eqn, opt, stat1));
    CALL(mess_lrcfadi_nm(eqn,opt,stat1,Z));
    CALL(mess_matrix_dynorm2(B,&rel));
    CALL(mess_lrcfadi_residual(eqn, opt, Z, &res2));
    printf("rel = %e \t res2 = %e \t %e\n",rel, res2, res2/rel );
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
     *  Solve the equation  op = TRANSPOSE
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    opt->adi_shifts_paratype    = paratype;
    opt->nm_gpStep              = nm_gpstep;
    opt->nm_singleshifts        = singleshifts;
    opt->nm_maxit               = 25;
    opt->adi_output             = 1;
    opt->nm_output              = 1;
    opt->nm_res2_tol            = 1e-7;
    opt->type                   = MESS_OP_TRANSPOSE;
    opt->nm_linesearch          = linesearch;
    opt->residual_method        = atoi(argv[11]);
    CALL(mess_options_print(opt));

    /*-----------------------------------------------------------------------------
     * Setup, Solve and Clear Riccati Equation for MESS_OP_TRANSPOSE
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_init(&stat2));
    CALL(mess_equation_init(&eqn));
    CALL(mess_equation_griccati(eqn, opt, A, E, B, C));
    CALL(mess_parameter(eqn, opt, stat2));
    CALL(mess_lrcfadi_nm(eqn,opt,stat2,Z));
    CALL(mess_matrix_dynorm2(C,&rel));
    CALL(mess_lrcfadi_residual(eqn, opt, Z, &res2));
    printf("rel = %lg \t res2 = %lg \t %lg\n",rel, res2, res2/rel );
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
     *  STATUS PRINT for MESS_OP_NONE and MESS_OP_TRANSPOSE
     *-----------------------------------------------------------------------------*/
    printf("Status - op = MESS_OP_NONE \n");
    CALL(mess_status_print(stat1));
    CALL(mess_status_printfull(stat1));

    printf("\nStatus - op = MESS_OP_TRANSPOSE \n");
    CALL(mess_status_print(stat2));
    CALL(mess_status_printfull(stat2));

    mess_status_clear(&stat1);
    mess_status_clear(&stat2);
    mess_options_clear(&opt);
    MESS_CLEAR_MATRICES(&A,&B,&C,&Z,&E);
    mess_equation_clear(&eqn);

    printf("ok1: %d ok2: %d\n", ok1, ok2);
    if ( ok1 && ok2  ) return 0;
    return 1;

}

