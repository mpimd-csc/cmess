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
 * @file tests/lrcf_adi/check_lradi_dae1.c
 * @brief Check the solution of a Lyapunov Equation using the LRCF-ADI method and DAE 1 function handles.
 * @author @mbehr
 *
 * This function checks the @ref mess_lrcfadi_adi function for the DAE 1 function handles defined in lrcf_adi/equation_glyap_dae1.c.
 *
 * @test
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
    int ret;
    mess_int_t ok = 0;
    double res2,rel;
    mess_init();
    mess_error_level = 3;
    mess_version();

    mess_matrix E11,A11,A12,A21,A22,RHS,Z;
    mess_options opt;
    mess_status stat;
    mess_equation eqn;


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/

    if ( argc != 10) {
        fprintf(stderr, "usage: %s E11.mtx A11.mtx A12.mtx A21.mtx A22.mtx RHS.mtx paratype direct mdirect\n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  Init and read data
     *-----------------------------------------------------------------------------*/
    //init matrices
    MESS_INIT_MATRICES(&E11,&A11,&A12,&A21,&A22,&RHS,&Z);

    //read matrices
    CALL(mess_matrix_read_formated(argv[1], E11, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[2], A11, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[3], A12, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[4], A21, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[5], A22, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[6], RHS, MESS_DENSE));


    /*-----------------------------------------------------------------------------
     *  Setup options structure
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    //choose automatically the correct lyapunov equation, with help of shape of B
    opt->type                       = (RHS->cols<RHS->rows) ? MESS_OP_NONE: MESS_OP_TRANSPOSE;
    opt->adi_res2_tol               = 1e-6;
    opt->adi_maxit                  = 10000;
    opt->adi_shifts_paratype        = atoi(argv[7]);
    opt->adi_output                 = 1;
    CALL(mess_options_print(opt));

    //select direct solver
    mess_direct_lupackage_t direct = atoi(argv[8]);
    mess_direct_lu_select(direct);

    mess_multidirect_t mdirect = atoi(argv[9]);
    mess_multidirect_select(mdirect);


    /*-----------------------------------------------------------------------------
     *  Solve Lyapunov Equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_init(&stat));
    CALL(mess_equation_init(&eqn));
    CALL(mess_equation_glyap_dae1(eqn,opt,E11,A11,A12,A21,A22,RHS));

    /*-----------------------------------------------------------------------------
     *  compute shift parameters
     *-----------------------------------------------------------------------------*/
    CALL(mess_parameter(eqn, opt, stat));

    /*-----------------------------------------------------------------------------
     *  solve lyapunov equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_lrcfadi_adi(eqn,opt,stat,Z));

    /*-----------------------------------------------------------------------------
     *  calculate residual
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dynorm2(eqn->RHS,&rel));
    CALL(mess_lrcfadi_residual(eqn, opt, Z, &res2));
    printf("res2 = %e \t rel= %e \t res2/rel = %e\n",res2,rel,res2/rel );

    /*-----------------------------------------------------------------------------
     *  print status
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_print(stat));

    ok = stat->stop_res2 || stat->stop_res2c || stat->stop_rel || stat->stop_user;

    /*-----------------------------------------------------------------------------
     *  Clean Data
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&E11,&A11,&A12,&A21,&A22,&RHS,&Z);
    CALL(mess_status_clear(&stat));
    CALL(mess_options_clear(&opt));
    CALL(mess_equation_clear(&eqn));

    mess_exit();
    if ( ok ) return 0;
    return 1;
}

