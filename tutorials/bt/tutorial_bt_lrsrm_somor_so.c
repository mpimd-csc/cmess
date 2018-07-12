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
 * @file tutorials/bt/tutorial_bt_lrsrm_somor_so.c
 * @brief  Demonstrate how to reduce a second order LTI system to another second order LTI system using a low-rank square root balanced truncation (BT-LRSRM-SOMOR-SO).
 * @author @koehlerm
 *
 * @sa mess_bt_lrsrm_somor_so_gram
 *
 * # Tutorial: Reduce a Second Order LTI to a Second Order LTI System Using a Low-rank Square Root Balanced Trunction Method
 *
 * This function demonstrates how to convert a second order linear time invariant (LTI) system
 * \f[
 * \begin{array}{ccccc}
 * M \ddot{x}(t) + G \dot{x}(t) + & K x(t) &=& B u(t) &  \\
 *  & y(t) &=& C_p x(t) & + C_v \dot{x}(t)
 * \end{array}
 * \f]
 * to a another reduced second order LTI system
 * using a low-rank square root balanced truncation.
 *
 *
 *
 * @snippet "tutorials/bt/tutorial_bt_lrsrm_somor_so.c" CODE
 *
 */

///@cond
///[CODE]
#include <stdlib.h>
#include <stdarg.h>
#include "mess/mess.h"

#define NSAMPLE 200

int main ( int argc, char **argv){
    mess_matrix M,K,G,B,Cv,Cp;
    mess_matrix *Go, *Gr;
    mess_matrix V,W;
    mess_options adiopt;
    mess_bt_options btopt;
    mess_bt_status  btstat;
    mess_dynsys sys, red;
    mess_int_t maxr = 0 ;
    double tol = 0;
    mess_vector omega,G1, G2;
    mess_int_t i;
    double err2, rel2, errinf, relinf;
    mess_vector diff;

    printf("mess - BT LRSRM\n");
    printf("================\n");

    if ( argc != 9 ) {
        printf("usage: %s M.mtx G.mtx K.mtx B.mtx Cp.mtx Cv.mtx tol maxr\n", argv[0]);
        return 0;
    }

    /*-----------------------------------------------------------------------------
     *  Init Data
     *-----------------------------------------------------------------------------*/

    mess_matrix_init(&M);
    mess_matrix_init(&G);
    mess_matrix_init(&K);
    mess_matrix_init(&B);
    mess_matrix_init(&Cv);
    mess_matrix_init(&Cp);
    mess_dynsys_init(&sys);
    mess_matrix_init(&V);
    mess_matrix_init(&W);
    mess_dynsys_init(&red);

    mess_set_errorlevel(1);


    /*-----------------------------------------------------------------------------
     *  Read Data
     *-----------------------------------------------------------------------------*/

    mess_matrix_read_formated(argv[1], M, MESS_CSR);
    mess_matrix_read_formated(argv[2], G, MESS_CSR);
    mess_matrix_read_formated(argv[3], K, MESS_CSR);
    mess_matrix_read_formated(argv[4], B, MESS_DENSE);
    mess_matrix_read_formated(argv[5], Cp, MESS_DENSE);
    mess_matrix_read_formated(argv[6], Cv, MESS_DENSE);
    tol = atof(argv[7]);
    maxr = atoi(argv[8]);

    printf("Dimension: " MESS_PRINTF_INT "\n", M->rows);
    printf("Inputs:    " MESS_PRINTF_INT "\n", B->cols);
    printf("Outputs:   " MESS_PRINTF_INT "\n\n", Cp->rows);

    mess_dynsys_2nd(sys,M,G,K,B,Cp,Cv);


    /*-----------------------------------------------------------------------------
     *  Evaluate original transferfunction
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&omega,&G1);
    mess_vector_alloc(omega,NSAMPLE,MESS_REAL);
    mess_vector_alloc(G1, NSAMPLE, MESS_REAL);
    Go = malloc (sizeof(mess_matrix)*NSAMPLE);
    for ( i = 0 ; i < NSAMPLE; i++){
        mess_matrix_init(&Go[i]);
    }
    mess_dynsys_evaltransfer(sys,-2,4,NSAMPLE,omega,G1,Go);

    /*-----------------------------------------------------------------------------
     *  setup BT
     *-----------------------------------------------------------------------------*/
    mess_bt_options_init(&btopt);
    mess_bt_status_init(&btstat);
    mess_options_init(&adiopt);
    adiopt->adi_output = 0 ;

    btopt->rdim =  maxr;
    btopt->tol = tol;
    btopt->so_type = MESS_BT_POSITION_VELOCITY;

    btopt->chooseorder = mess_bt_chooseorder_default;

    /*-----------------------------------------------------------------------------
     *  compute V and W
     *-----------------------------------------------------------------------------*/
    printf("\nBT options:\n");
    mess_bt_options_print(btopt);
    mess_bt_lrsrm_somor_so(sys,btopt,adiopt,  V,W, btstat);
    printf("\nStatus of mess_bt_lrsrm_somor:\n");
    mess_bt_status_print(btstat);
    printf("\n");

    /*-----------------------------------------------------------------------------
     *  project the system
     *-----------------------------------------------------------------------------*/
    mess_dynsys_project_2nd_to_2nd(sys, V,W,red);

    /*-----------------------------------------------------------------------------
     *  eval red. TFM
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&G2);
    mess_vector_alloc(G2, NSAMPLE, MESS_REAL);
    Gr = malloc (sizeof(mess_matrix)*NSAMPLE);
    for ( i = 0 ; i < NSAMPLE; i++){
        mess_matrix_init(&Gr[i]);
    }
    mess_dynsys_evaltransfer(red,-2,4,NSAMPLE,omega,G2,Gr);


    /*-----------------------------------------------------------------------------
     *  Compute the error
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&diff);
    mess_vector_alloc(diff, NSAMPLE, MESS_REAL);
    mess_dynsys_transfer_diff(NSAMPLE, Go, Gr,&err2, &rel2, &errinf, &relinf, diff);
    mess_vector_clear(&diff);

    printf("error discrete ||G1-G2||_2           = %lg\n", err2);
    printf("error discrete ||G1-G2||_2 rel.      = %lg\n", rel2);
    printf("error discrete ||G1-G2||_infty       = %lg\n", errinf);
    printf("error discrete ||G1-G2||_inftyi rel. = %lg\n", relinf);


    mess_vector_clear(&G2);
    mess_matrix_clear(&V);
    mess_matrix_clear(&W);
    mess_dynsys_clear(&red);

    for ( i = 0; i < NSAMPLE; i++){
        mess_matrix_clear(&Go[i]);
        mess_matrix_clear(&Gr[i]);
    }
    free(Go); Go=NULL;
    free(Gr); Gr=NULL;


    mess_vector_clear(&G1);
    mess_vector_clear(&omega);
    mess_dynsys_clear(&sys);
    mess_options_clear(&adiopt);
    mess_bt_options_clear(&btopt);
    mess_bt_status_clear(&btstat);
    return 0;
}
///[CODE]
///@endcond
