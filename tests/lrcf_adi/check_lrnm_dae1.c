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
 * @file tests/lrcf_adi/check_lrnm_dae1.c
 * @brief Check the solution of a Lyapunov Equation using the LRCF-ADI method and DAE 1 function handles.
 * @author @mbehr
 *
 *
 * This function checks the @ref mess_lrcfadi_nm function for the DAE 1 function handles defined in lrcf_adi/equation_glyap_dae1.c.
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
    mess_error_level = 4;
    mess_version();

    mess_matrix E11,A11,A12,A21,A22,B,C,Z;
    mess_options opt;
    mess_status stat;
    mess_equation eqn;

    mess_version();

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/

    if ( argc != 15) {
        fprintf(stderr, "usage: %s E11.mtx A11.mtx A12.mtx A21.mtx A22.mtx B.mtx C.mtx optype paratype singleshifts gpstep direct mdirect linesearch\n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  Init and read data
     *-----------------------------------------------------------------------------*/
    //init matrices
    MESS_INIT_MATRICES(&E11,&A11,&A12,&A21,&A22,&B,&C,&Z);

    //read matrices
    CALL(mess_matrix_read_formated(argv[1], E11, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[2], A11, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[3], A12, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[4], A21, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[5], A22, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[6], B, MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[7], C, MESS_DENSE));


    /*-----------------------------------------------------------------------------
     *  Setup options structure
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    opt->type                       = (atoi(argv[8])==0) ? MESS_OP_NONE:MESS_OP_TRANSPOSE;
    opt->adi_res2_tol               = 1e-10;
    opt->adi_maxit                  = 1500;
    opt->adi_shifts_paratype        = atoi(argv[9]);
    opt->nm_singleshifts            = atoi(argv[10]);
    opt->nm_gpStep                  = atoi(argv[11]);
    opt->nm_linesearch              = atoi(argv[14]);
    opt->adi_output                 = 0;
    opt->nm_output                  = 1;
    opt->nm_gpStep                  = 0;
    opt->nm_maxit                   = 50;

    CALL(mess_options_print(opt));

    //select direct solver
    mess_direct_lupackage_t direct = atoi(argv[12]);
    mess_direct_lu_select(direct);

    mess_multidirect_t mdirect = atoi(argv[13]);
    mess_multidirect_select(mdirect);

    /*-----------------------------------------------------------------------------
     *  Solve Riccati Equation
     *-----------------------------------------------------------------------------*/
    //set controls in algebraic constraints to zero

    if(opt->type==MESS_OP_TRANSPOSE){
        CALL(mess_matrix_resize(B,E11->rows,B->cols));
        CALL(mess_matrix_resize(C,C->rows,E11->cols+A22->cols));
    }else{
        CALL(mess_matrix_resize(B,E11->rows + A22->rows,B->cols));
        CALL(mess_matrix_resize(C,C->rows,E11->cols));
    }

    CALL(mess_status_init(&stat));
    CALL(mess_equation_init(&eqn));
    CALL(mess_equation_griccati_dae1(eqn,opt,E11,A11,A12,A21,A22,B,C));

    /*-----------------------------------------------------------------------------
     *  compute shift parameters
     *-----------------------------------------------------------------------------*/
    CALL(mess_parameter(eqn, opt, stat));

    /*-----------------------------------------------------------------------------
     *  solve lyapunov equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_lrcfadi_nm(eqn,opt,stat,Z));

    /*-----------------------------------------------------------------------------
     *  calculate residual
     *-----------------------------------------------------------------------------*/
    if(opt->type==MESS_OP_NONE){
        CALL(mess_matrix_dynorm2(eqn->B,&rel));
    }else{
        CALL(mess_matrix_dynorm2(eqn->C,&rel));
    }
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
    MESS_CLEAR_MATRICES(&E11,&A11,&A12,&A21,&A22,&B,&C,&Z);
    CALL(mess_status_clear(&stat));
    CALL(mess_options_clear(&opt));
    CALL(mess_equation_clear(&eqn));

    mess_exit();
    if ( ok ) return 0;
    return 1;
}

