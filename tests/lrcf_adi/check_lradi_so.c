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
 * @file tests/lrcf_adi/check_lradi_so.c
 * @brief Check the solution of a Lyapunov Equation arising from a second order system.
 * @test
 * This function checks the @ref mess_lrcfadi_adi function defined in lradi.c that means it checks if a Lyapunov Equation
 * \f[ A X E^T + E X A^T + B B^T = 0 \f]
 * is solved correctly using the low-rank Cholesky factor alternate direction implicit (LRCF-ADI) method.\n
 * This Lyapunov Equation is generated from a second order system
 * \f[ M \ddot{x} + D \dot{x} + Kx = Bu \f]
 * either using the @ref mess_equation_glyap_so1 function defined in equation_glyap_so1.c, where
 * \f[
 *  \begin{array}{cc}
 *   E & = \left[
 *          \begin{array}{cc}
 *           -K &  0 \\  0 &  M
 *          \end{array}
 *         \right] \\
 *   A & = \left[
 *          \begin{array}{cc}
 *           0 & -K \\ -K & -D
 *          \end{array}
 *         \right]\\
 *   B & = \left[
 *          \begin{array}{c}
 *              0 \\ B
 *          \end{array}
 *        \right]
 *  \end{array}
 * \f]
 * or using the @ref mess_equation_glyap_so2 function defined in equation_glyap_so2.c, where
 * \f[
 * \begin{array}{cc}
 * E & = \left[
 *  \begin{array}{cc}
 *      -I &  0 \\  0 &  M
 *  \end{array}
 *   \right] \\
 * A &= \left[
 *  \begin{array}{cc}
 *      0 & I \\ -K & -D
 *  \end{array}
 *    \right]\\
 *
 * B & = \left[
 *  \begin{array}{c}
 *      0 \\ B
 *  \end{array}
 *  \right]
 * \end{array}.
 *  \f]
 *
 * @}
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "../call_macro.h"



int main ( int argc, char **argv){
    mess_init();
    mess_error_level = 3;
    mess_version();

    mess_matrix M, D, K, B, Z;
    mess_options opt;
    mess_equation eqn;
    mess_status stat;
    mess_int_t paratype, variant;
    mess_operation_t optype;
    int ret=0, ok=0;
    double rel=.0, res2=.0;

    /*-----------------------------------------------------------------------------
     *  check input args
     *-----------------------------------------------------------------------------*/
    if ( argc != 9) {
        printf("usage: %s <M.mtx> <D.mtx> <K.mtx> optype variant paratype direct mdirect\n", argv[0]);
        return 1;

    }
    optype      = (atoi(argv[4])==0)? MESS_OP_NONE:MESS_OP_TRANSPOSE;
    variant     = atoi(argv[5]);
    paratype    = atoi(argv[6]);

    /*-----------------------------------------------------------------------------
     *  Read matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&M,&D,&K,&B,&Z);
    CALL(mess_matrix_read_formated(argv[1], M, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[2], D, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[3], K, MESS_CSR));

    //set random entries to RHS input matrix B
    if(optype == MESS_OP_NONE){
        //CALL(mess_matrix_rand(B,2*M->rows,1,MESS_DENSE,1));
        mess_matrix_alloc(B,2*M->rows,1,2*M->rows,MESS_DENSE,MESS_REAL);
        CALL(mess_matrix_one_value(B,1));
    }else{
        CALL(mess_matrix_alloc(B,1,2*M->rows,2*M->rows,MESS_DENSE,MESS_REAL));
        CALL(mess_matrix_one_value(B,1));
    }

    /*-----------------------------------------------------------------------------
     *  Solve the Lyapunov Equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_status_init(&stat));
    CALL(mess_options_init(&opt));
    opt->type                   = optype;
    opt->adi_shifts_paratype    = paratype;
    opt->adi_output             = 1;

    CALL(mess_options_print(opt));
    CALL(mess_equation_init(&eqn));

    //select direct solver
    mess_direct_lupackage_t direct = atoi(argv[7]);
    mess_direct_lu_select(direct);

    mess_multidirect_t mdirect = atoi(argv[8]);
    mess_multidirect_select(mdirect);


    //depending on input use so1 or so2
    if(variant==1){
        printf("mess_equation_glyap_so1\n");
        CALL(mess_equation_glyap_so1(eqn, opt, M, D, K, B,1e-8,1e+8));
    }else{
        printf("mess_equation_glyap_so2\n");
        CALL(mess_equation_glyap_so2(eqn, opt, M, D, K, B,1e-8,1e+8));
    }

    CALL(mess_parameter(eqn, opt, stat));
    CALL(mess_lrcfadi_adi(eqn,opt,stat,Z));
    CALL(mess_matrix_dynorm2(B,&rel));
    CALL(mess_lrcfadi_residual(eqn,opt,Z, &res2));
    printf("res2 = %lg \t rel= %lg \t res2/rel = %lg\n",res2,rel,res2/rel );
    ok= stat->stop_res2 || stat->stop_res2c || stat->stop_rel || stat->stop_user;

    printf("Status - op = %s \n",mess_operation_t_str(optype));
    CALL(mess_status_print(stat));
    CALL(mess_options_print(opt));

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&M,&D,&K,&B,&Z);
    CALL(mess_equation_clear(&eqn));
    CALL(mess_options_clear(&opt));
    CALL(mess_status_clear(&stat));

    mess_exit();
    printf("ok: %d \n", ok);
    if ( ok ) return 0;
    return 1;
}


