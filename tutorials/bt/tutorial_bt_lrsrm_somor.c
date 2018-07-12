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
 * @file tutorials/bt/tutorial_bt_lrsrm_somor.c
 * @brief Demonstrate how to reduce a second order LTI system to a first order LTI system using a low-rank square root balanced truncation (BT-LRSRM-SOMOR).
 * @author @koehlerm
 *
 * @sa mess_bt_lrsrm_somor
 *
 * # Tutorial: Low-rank Square Root Balanced Trunctation Second Order LTI to First Order LTI System
 *
 * This function demonstrates how to convert a second order linear time invariant (LTI) system
 * \f[
 * \begin{array}{ccccc}
 * M \ddot{x}(t) + G \dot{x}(t) + & K x(t) &=& B u(t) &  \\
 *  & y(t) &=& C_p x(t) & + C_v \dot{x}(t)
 * \end{array}
 * \f]
 * to a reduced first order LTI system
 * \f[
 * \begin{array}{ccc}
 * \dot{x} &= A x & + B u \\
 * y &= C x&
 *  \end{array}
 * \f]
 * using a low-rank square root balanced truncation.
 *
 *
 * @snippet "tutorials/bt/tutorial_bt_lrsrm_somor.c" CODE
 *
 */

///@cond
///[CODE]
#include <stdlib.h>
#include <stdarg.h>
#include "mess/mess.h"

int main ( int argc, char **argv){
    mess_matrix M,K,G,B,Cv,Cp;
    mess_matrix V,W;
    mess_options adiopt;
    mess_bt_options btopt;
    mess_bt_status  btstat;
    mess_dynsys sys, red;
    mess_int_t maxr = 0 ;
    double tol = 0, err1, err2;
    mess_vector omega,G1, G2;

    printf("Balanced Truncation Example for second order LTI Systems\n");
    printf("========================================================\n");

    mess_set_errorlevel(1);

    /*-----------------------------------------------------------------------------
     *  Read matrices
     *-----------------------------------------------------------------------------*/
    if ( argc != 9 ) {
        printf("usage: %s M.mtx G.mtx K.mtx B.mtx Cp.mtx Cv.mtx tol maxr\n", argv[0]);
        return 0;
    }
    mess_init();


    mess_matrix_init(&M);
    mess_matrix_init(&G);
    mess_matrix_init(&K);
    mess_matrix_init(&B);
    mess_matrix_init(&Cv);
    mess_matrix_init(&Cp);
    mess_dynsys_init(&sys);
    mess_dynsys_init(&red);

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


    /*-----------------------------------------------------------------------------
     *  Setup BT
     *-----------------------------------------------------------------------------*/
    mess_dynsys_2nd(sys,M,G,K,B,Cp,Cv);

    mess_options_init(&adiopt);
    mess_bt_options_init(&btopt);
    mess_bt_status_init(&btstat);
    btopt->rdim =  maxr;
    btopt->tol = tol;
    adiopt->adi_output = 0;
    adiopt->adi_maxit = 100;

    mess_matrix_init(&V);
    mess_matrix_init(&W);
    printf("BT options\n");
    mess_bt_options_print(btopt);


    /*-----------------------------------------------------------------------------
     *  Compute the Projection matrices
     *-----------------------------------------------------------------------------*/
    mess_bt_lrsrm_somor(sys,btopt, adiopt, V,W, btstat);

    printf("\nStatus of mess_bt_lrsrm_somor:\n");
    mess_bt_status_print(btstat);
    printf("\n");


    /*-----------------------------------------------------------------------------
     *  Compute the Reduced Order Model
     *-----------------------------------------------------------------------------*/
    mess_dynsys_project_2nd_to_1st(sys, V,W,red);


    /*-----------------------------------------------------------------------------
     *  Evaluate the Transferfunction
     *-----------------------------------------------------------------------------*/
#define NSAMPLE 200

    MESS_INIT_VECTORS(&omega,&G1,&G2);
    mess_vector_alloc(omega,NSAMPLE,MESS_REAL);
    mess_vector_alloc(G1, NSAMPLE, MESS_REAL);
    mess_vector_alloc(G2, NSAMPLE, MESS_REAL);

    mess_dynsys_evaltransfer(sys,0,5,NSAMPLE,omega,G1,NULL);
    mess_dynsys_evaltransfer(red,0,5,NSAMPLE,omega,G2,NULL);

    mess_vector_diffnorm(G1,G2,&err1);
    mess_vector_diffnorminf(G1,G2,&err2);

    printf("error discrete ||G1-G2||_2 = %lg\n", err1);
    printf("error discrete ||G1-G2||_infty = %lg\n", err2);


    /*-----------------------------------------------------------------------------
     *  Cleanup
     *-----------------------------------------------------------------------------*/
    mess_vector_clear(&G1);
    mess_vector_clear(&G2);
    mess_vector_clear(&omega);

    mess_matrix_clear(&V);
    mess_matrix_clear(&W);
    mess_dynsys_clear(&sys);
    mess_dynsys_clear(&red);
    mess_options_clear(&adiopt);
    mess_bt_options_clear(&btopt);
    mess_bt_status_clear(&btstat);
    mess_exit();
    return 0;
}
///[CODE]
///@endcond
