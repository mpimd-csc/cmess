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
 * @file tests/lrcf_adi/check_lradi_dae2.c
 * @brief Check the solution of a Lyapunov Equation using the LRCF-ADI method and DAE 2 function handles.
 *
 * This function checks the @ref mess_lrcfadi_adi function defined in lrcf_adi/equation_glyap_dae2.c.
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
    mess_init();
    mess_error_level = 1;
    mess_version();

    int ret, ok1 = 0 ;
    mess_int_t  have_feed=0;
    double rel,res2, delta = -0.02;
    mess_matrix M, A, G, B, Z, K0=NULL;
    mess_options opt;
    mess_status stat;
    mess_equation eqn, eqn_stable;

    mess_version();

    if ( argc != 10 && argc != 9) {
        fprintf(stderr, "usage: %s M.mtx A.mtx G.mtx B/C.mtx memory_usage paratype direct mdirect [K0.mtx]\n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  Init and read data
     *-----------------------------------------------------------------------------*/
    //init matrices
    MESS_INIT_MATRICES(&M,&A,&G,&B,&Z);

    //read matrices
    CALL(mess_matrix_read_formated(argv[1], M, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[2], A, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[3], G, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[4], B, MESS_DENSE));
    if(argc==10){
        have_feed=1;
        CALL(mess_matrix_init(&K0));
        CALL(mess_matrix_read_formated(argv[9], K0, MESS_DENSE));
    }


    /*-----------------------------------------------------------------------------
     *  Setup options structure
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    //choose automatically the correct lyapunov equation, with help of shape of B
    opt->type                       = (B->cols<B->rows) ? MESS_OP_NONE: MESS_OP_TRANSPOSE;
    opt->adi_shifts_arp_m           = 40;
    opt->adi_shifts_arp_p           = 40;
    opt->adi_shifts_l0              = 25;
    opt->memory_usage               = atoi(argv[5]);
    opt->adi_shifts_paratype        = atoi(argv[6]);
    opt->adi_maxit                  = 1500;
    opt->adi_output                 = 1;

    CALL(mess_options_print(opt));

    //select direct solver
    mess_direct_lupackage_t direct = atoi(argv[7]);
    mess_direct_lu_select(direct);

    mess_multidirect_t mdirect = atoi(argv[8]);
    mess_multidirect_select(mdirect);

    /*-----------------------------------------------------------------------------
     *  Solve Lyapunov Equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_init(&stat));
    CALL(mess_equation_init(&eqn));
    CALL(mess_equation_init(&eqn_stable));
    CALL(mess_equation_glyap_dae2(eqn,opt,M,A,G,B,delta));

    if(have_feed){
        /*--------stabilize equation with given feedback----*/
        if(opt->type == MESS_OP_NONE){
            //build (A-BK) Operator for stabilization
            CALL(mess_equation_stable(eqn_stable,opt,eqn,eqn->B,K0));
        }else{
            // build (A-KC) Operator for stabilization
            CALL(mess_equation_stable(eqn_stable,opt,eqn,K0,eqn->B))
        }

        /*-------------compute shift parameters-------------*/
        CALL(mess_parameter(eqn_stable, opt, stat));

        /*-------------solve lyapunov equation--------------*/
        CALL(mess_lrcfadi_adi(eqn_stable,opt,stat,Z));

        /*-------------compute residual---------------------*/
        //the projection of the RHS is done by mess_lrcfadi_adi int init_rhs
        CALL(mess_matrix_dynorm2(eqn_stable->RHS,&rel));
        CALL(mess_lrcfadi_residual(eqn_stable, opt, Z, &res2));
        printf("res2 = %lg \t rel= %lg \t res2/rel = %lg\n",res2,rel,res2/rel );

    }else{
        //no initial feedback
        /*-------------compute shift parameters-------------*/
        CALL(mess_parameter(eqn, opt, stat));

        /*-------------solve lyapunov equation--------------*/
        CALL(mess_lrcfadi_adi(eqn,opt,stat,Z));

        /*-------------compute residual---------------------*/
        CALL(mess_matrix_dynorm2(eqn->RHS,&rel));
        CALL(mess_lrcfadi_residual(eqn, opt, Z, &res2));
        printf("res2 = %lg \t rel= %lg \t res2/rel = %lg\n",res2,rel,res2/rel );
    }

    /*-----------------------------------------------------------------------------
     *  print status
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_print(stat));

    ok1 = stat->stop_res2 || stat->stop_res2c || stat->stop_rel || stat->stop_user;

    /*-----------------------------------------------------------------------------
     *  Clean Data
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&M,&A,&G,&B,&Z);
    if(K0){
        CALL(mess_matrix_clear(&K0));
    }

    CALL(mess_status_clear(&stat));
    CALL(mess_options_clear(&opt));
    CALL(mess_equation_clear(&eqn_stable));
    CALL(mess_equation_clear(&eqn));

    mess_exit();
    if ( ok1 ) return 0;
    return 1;
}

