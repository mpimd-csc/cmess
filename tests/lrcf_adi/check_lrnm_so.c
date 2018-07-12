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
 * @file tests/lrcf_adi/check_lrnm_so.c
 * @brief Check the solution of a Riccati Equation arising from a second order system.
 * @test
 * This function checks the @ref mess_lrcfadi_nm function defined in lrnm.c that means it checks if the Riccati Equation
 * \f[ AXE^T+EXA^T-EXC^TCXE^T+BB^T = 0 \f]
 * or
 * \f[ A^TXE+E^TXA-EXBB^TXE^T+C^TC = 0 \f]
 * is solved correctly using the low-rank-Cholesky factor Newton method (LRCF-NM).\n
 * This Riccati Equation is generated from a second order system
 * \f[\begin{array}{cccc}
 *    M \ddot{x} + D  \dot{x} + K & x &=& B u \\
 *    & y & = & C x
 *   \end{array} \f]
 * either using the @ref mess_equation_griccati_so1 function defined in equation_glyap_so1.c, where
 * \f[
 * \begin{array}{cc}
 *  E & = \left[
 *          \begin{array}{cc}
 *              -K &  0 \\  0 &  M
 *          \end{array}
 *        \right] \\
 *  A &=  \left[
 *          \begin{array}{cc}
 *              0 & -K \\ -K & -D
 *          \end{array}
 *       \right]\\
 *  B & = \left[
 *          \begin{array}{c}
 *              0 \\ B
 *          \end{array}
 *       \right]\\
 *  C & = \left[
 *          \begin{array}{cc}
 *              C & 0
 *          \end{array}
 *        \right]
 * \end{array}
 * \f]
 * or using the @ref mess_equation_griccati_so2 function defined in equation_glyap_so2.c, where
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
 *  \right]\\
 * C & = \left[
 *  \begin{array}{cc}
 *      C & 0
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
    mess_error_level = 1;
    mess_version();

    mess_matrix M, D, K, C, B, Z;
    mess_options opt;
    mess_equation eqn;
    mess_status stat;
    mess_int_t  optype, paratype, nm_gp, singleshifts, variant;
    mess_int_t colsB=1, rowsC=1, linesearch;
    int ret=0;
    double rel=.0, res2=.0;
    int ok=0;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/

    if ( argc != 11) {
        printf("usage: %s <M.mtx> <D.mtx> <K.mtx>  optype variant paratype nm_gp singleshifts direct linesearch\n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  Read matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&M,&D,&K,&C,&B,&Z);
    CALL(mess_matrix_read_formated(argv[1], M, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[2], D, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[3], K, MESS_CSR));
    CALL(mess_matrix_sort(M));
    CALL(mess_matrix_sort(D));
    CALL(mess_matrix_sort(K));

    /*-----------------------------------------------------------------------------
     *  setup options
     *-----------------------------------------------------------------------------*/
    optype          = atoi(argv[4])==0 ? MESS_OP_NONE: MESS_OP_TRANSPOSE;
    variant         = atoi(argv[5]);
    paratype        = atoi(argv[6]);
    nm_gp           = atoi(argv[7]);
    singleshifts    = atoi(argv[8]);
    linesearch      = atoi(argv[10]);


    //set all entries to input matrix B and output Matrix C to 1
    if(optype==MESS_OP_NONE){
        CALL(mess_matrix_alloc(B, 2*M->rows, colsB, colsB*2*M->rows, MESS_DENSE, MESS_REAL));
        CALL(mess_matrix_alloc(C, rowsC, M->rows, rowsC*M->rows, MESS_DENSE, MESS_REAL));
        CALL(mess_matrix_ones(B));
        CALL(mess_matrix_ones(C));
        mess_int_t i,j=0;
        for(i=0;i<B->rows;++i){
            for(j=0;j<B->cols;++j){
                if(i<(B->rows+1)/2){ B->values[i+j*(B->ld)]=0; }
            }
        }
        //CALL(mess_matrix_rand(B,2*M->rows,colsB,MESS_DENSE,0.5));
        //CALL(mess_matrix_rand(C,rowsC,M->rows,MESS_DENSE,0.5));
    }else{
        CALL(mess_matrix_alloc(B, M->rows, colsB,   colsB*M->rows, MESS_DENSE, MESS_REAL));
        CALL(mess_matrix_alloc(C, rowsC, 2*M->rows, rowsC*2*M->rows, MESS_DENSE, MESS_REAL));
        CALL(mess_matrix_ones(B));
        CALL(mess_matrix_ones(C));
        mess_int_t i,j=0;
        for(i=0;i<C->rows;++i){
            for(j=0;j<C->cols;++j){
                if(j<(C->cols+1)/2){ C->values[i+j*(C->ld)]=0; }
            }
        }
        //CALL(mess_matrix_rand(B,M->rows,colsB,MESS_DENSE,0.5));
        //CALL(mess_matrix_rand(C,rowsC,2*M->rows,MESS_DENSE,0.5));
    }

    /*-----------------------------------------------------------------------------
     *  Solve the Riccati Equation
     *-----------------------------------------------------------------------------*/
    CALL(mess_options_init(&opt));
    opt->adi_shifts_paratype    = paratype;
    opt->nm_gpStep              = nm_gp;
    opt->nm_singleshifts        = singleshifts;
    opt->adi_output             = 1;
    opt->nm_output              = 1;
    opt->type                   = optype;
    opt->nm_linesearch          = linesearch;
    CALL(mess_options_print(opt));
    CALL(mess_status_init(&stat));
    CALL(mess_equation_init(&eqn));

    //select direct solver
    mess_direct_lupackage_t direct = atoi(argv[9]);
    mess_direct_lu_select(direct);

    if(variant==1){
        CALL(mess_equation_griccati_so1(eqn, opt, M, D, K, B, C,1e-8,1e+8));
    }else{
        CALL(mess_equation_griccati_so2(eqn, opt, M, D, K, B, C,1e-8,1e+8));
    }

    CALL(mess_parameter(eqn, opt, stat));
    CALL(mess_lrcfadi_nm(eqn, opt, stat, Z ));
    CALL(mess_matrix_dynorm2((opt->type==MESS_OP_NONE)?B:C,&rel));
    CALL(mess_lrcfadi_residual(eqn,opt,Z, &res2));
    printf("rel = %lg \t res2 = %lg res2/rel = %lg\n",rel, res2, res2/rel );
    ok= stat->stop_res2 || stat->stop_rel || stat->stop_user;
    printf("Status - op = %s \n",mess_operation_t_str(optype));
    CALL(mess_status_print(stat));

    /*-----------------------------------------------------------------------------
     *  clear additional memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&Z,&M,&D,&K,&B,&C);
    CALL(mess_equation_clear(&eqn));
    CALL(mess_status_clear(&stat));
    CALL(mess_options_clear(&opt));

    mess_exit();
    printf("ok: %d\n", ok);
    if ( ok ) return 0;
    return 1;
}


