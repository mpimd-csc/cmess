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
 * @file tests/lrcf_adi/check_lradi_glyap.c
 * @brief Check the solution of a Lyapunov Equation using the LRCF-ADI method.
 * @test
 * This function checks the @ref mess_lrcfadi_adi function defined in lrcf_adi/lradi.c that means it checks if a Lyapunov Equation
 * \f[ A X E^T + E X A^T + B B^T = 0 \f]
 * is solved correctly using the low-rank Cholesky factor alternate direction implicit (LRCF-ADI) method.\n
 * This equation is generated using the @ref mess_equation_glyap function defined in equation_glyap.c.
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
    mess_init();
    mess_error_level=3;
    mess_version();

    mess_matrix  A, B, C, E, Z;
    mess_options opt;
    mess_status stat;
    mess_int_t paratype;
    mess_equation eqn;
    int ret, ok1 = 0, ok2 = 0;
    double res2, rel,adi_res2,res2_diffexp;

    /*-----------------------------------------------------------------------------
     *  check input args
     *-----------------------------------------------------------------------------*/
    if ( argc != 9) {
        fprintf(stderr, "usage: %s A.mtx E.mtx B.mtx C.mtx paratype direct mdirect resmethod \n", argv[0]);
        return 1;
    }

    paratype = atoi(argv[5]);

    /*-----------------------------------------------------------------------------
     *  Read matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&B,&C,&E,&Z);
    CALL(mess_matrix_read_formated(argv[1], A, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[2], E, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[3], B, MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[4], C, MESS_DENSE));

    /*-----------------------------------------------------------------------------
     *  setup options
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    opt->adi_shifts_paratype    = paratype;
    opt->type                   = MESS_OP_NONE;
    opt->residual_method        = atoi(argv[8]);
    CALL(mess_options_print(opt));

    mess_direct_lupackage_t direct = atoi(argv[6]);
    mess_direct_lu_select(direct);

    mess_multidirect_t mdirect = atoi(argv[7]);
    mess_multidirect_select(mdirect);

    /*-----------------------------------------------------------------------------
     *  setup equation and solve
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_init(&stat));
    CALL(mess_equation_init(&eqn));
    CALL(mess_equation_glyap(eqn, opt, A, E, B));
    CALL(mess_parameter(eqn, opt, stat));
    CALL(mess_lrcfadi_adi(eqn,opt,stat,Z));
    CALL(mess_matrix_dynorm2(B,&rel));
    CALL(mess_lrcfadi_residual(eqn, opt, Z, &res2));
    printf("res2 = %.9e \t res2/rel = %.9e\n",res2,res2/rel);
    CALL(mess_status_print(stat));
    CALL(mess_status_printfull(stat));
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


    /*-----------------------------------------------------------------------------
     *  reset memory
     *-----------------------------------------------------------------------------*/
    mess_matrix_clear(&Z);
    mess_status_clear(&stat);
    mess_options_clear(&opt);
    mess_equation_clear(&eqn);

    /*-----------------------------------------------------------------------------
     *  setup options
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    opt->adi_shifts_paratype    = paratype;
    opt->type                   = MESS_OP_TRANSPOSE;
    opt->residual_method        = atoi(argv[8]);
    CALL(mess_options_print(opt));

    /*-----------------------------------------------------------------------------
     *  setup equation and solve
     *-----------------------------------------------------------------------------*/
    CALL( mess_matrix_init(&Z));
    CALL(mess_status_init(&stat));
    CALL(mess_equation_init(&eqn));

    //intermediate test
    if ( mess_equation_glyap(eqn, opt, A, E, B) == 0 ) {
        fprintf(stderr, "mess_equation_glyap does not detect wrong sized RHS.\n");
        return 1;
    }
    CALL(mess_equation_glyap(eqn, opt, A, E, C));
    CALL(mess_parameter(eqn, opt, stat));
    CALL(mess_lrcfadi_adi(eqn,opt,stat,Z));
    CALL(mess_status_print(stat));
    CALL(mess_matrix_dynorm2(C,&rel));
    CALL(mess_lrcfadi_residual(eqn, opt, Z, &res2));
    printf("res2 = %.9e \t res2/rel = %.9e\n",res2,res2/rel );
    ok2 = stat->stop_res2 || stat->stop_res2c || stat->stop_rel || stat->stop_user;
    CALL(mess_status_print(stat));
    CALL(mess_status_printfull(stat));

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
     *  clear memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&B,&C,&E,&Z);
    mess_status_clear(&stat);
    mess_options_clear(&opt);
    mess_equation_clear(&eqn);

    printf("ok1: %d ok2: %d\n", ok1, ok2);
    if ( ok1 && ok2  ) return 0;
    return 1;
}

