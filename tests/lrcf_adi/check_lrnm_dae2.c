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
 * @file tests/lrcf_adi/check_lrnm_dae2.c
 * @brief Check the solution of generalized Riccati Equation arising from a Hessenberg index 2 DAE system using LRCF-NM.
 * @test
 *
 * This function uses the @ref mess_equation_griccati_dae2 function defined in equation_glyap_dae2.c.
 *
 * Model description available under @cite BaeBSetal14.
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

    int ret, ok = 0;
    mess_int_t optype=1;
    double res2, rel, delta = -0.02;
    mess_matrix  M, A, G, B, C, Z, K0;
    mess_equation eqn;
    mess_options opt;
    mess_status stat;
    mess_int_t have_initialfeed=0;


    /*-----------------------------------------------------------------------------
     *  check number of input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 13 && argc != 12 ) {
        fprintf(stderr, "usage: %s M.mtx A.mtx G.mtx B.mtx C.mtx op memory_usage para nm_gp singleshifts linesearch [K0.mtx] \n", argv[0]);
        return 1;
    }


    /*-----------------------------------------------------------------------------
     *  init and read matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&M,&A,&G,&B,&C,&Z);
    CALL(mess_matrix_read_formated(argv[1], M, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[2], A, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[3], G, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[4], B, MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[5], C, MESS_DENSE));
    optype = (atoi(argv[6])==0) ? MESS_OP_NONE:MESS_OP_TRANSPOSE;

    if(argc==13){
        have_initialfeed=1;
        mess_matrix_init(&K0);
        CALL(mess_matrix_read_formated(argv[12], K0, MESS_DENSE));
    }

    /*-----------------------------------------------------------------------------
     *  build opt structure
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    opt->type                       = optype;
    opt->adi_shifts_arp_m           = 40;
    opt->adi_shifts_arp_p           = 40;
    opt->adi_shifts_l0              = 25;
    opt->nm_output                  = 1;
    opt->nm_maxit                   = 25;
    opt->adi_output                 = 1;
    opt->adi_res2_tol               = 1e-12;
    opt->nm_res2_tol                = 1e-7;
    opt->adi_maxit                  = 500;
//    opt->adi_maxit                  = 5;
//    opt->nm_maxit                   = 2;
    opt->memory_usage               = atoi(argv[7]);
    opt->adi_shifts_paratype        = atoi(argv[8]);
    opt->nm_gpStep                  = atoi(argv[9]);
    opt->nm_singleshifts            = atoi(argv[10]);
    opt->nm_linesearch              = atoi(argv[11]);

    CALL(mess_options_print(opt));

    /*-----------------------------------------------------------------------------
     *  build stat structure
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_init(&stat));

    /*-----------------------------------------------------------------------------
     *  build eqn structure
     *-----------------------------------------------------------------------------*/
    CALL(mess_equation_init(&eqn));
    CALL(mess_equation_griccati_dae2(eqn,opt,M,A,G,B,C,delta));

    //set if initial feedback is given
    if(have_initialfeed){
        opt->K0 = K0;
    }else{
        opt->K0 = NULL;
    }

    /*-----------------------------------------------------------------------------
     *  solve generalized algebraic riccatti eqn
     *-----------------------------------------------------------------------------*/
    CALL(mess_lrcfadi_nm(eqn,opt,stat,Z));
    if(optype==MESS_OP_NONE){
        CALL(mess_matrix_dynorm2(eqn->B,&rel));
    }else{
        CALL(mess_matrix_dynorm2(eqn->C,&rel));
    }

    CALL(mess_lrcfadi_residual(eqn, opt, Z, &res2));
    printf("rel = %lg \t res2 = %lg res2/rel = %lg\n",rel, res2, res2/rel );
    ok = stat->stop_res2 || stat->stop_rel || stat->stop_user;

    /*-----------------------------------------------------------------------------
     *  Clear Data
     *-----------------------------------------------------------------------------*/
    //MESS_CLEAR_MATRICES(&M,&A,&G,&B,&C,&Z,&K0);
    //K0 is cleared by opt
    MESS_CLEAR_MATRICES(&M,&A,&G,&B,&C,&Z);


    CALL(mess_equation_clear(&eqn));
    CALL(mess_options_clear(&opt));
    CALL(mess_status_clear(&stat));


    if ( ok ) return 0;
    return 1;

}


