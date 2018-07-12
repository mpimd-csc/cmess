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
 *
 * @file tutorials/bt/tutorial_bt_lrsrm.c
 * @brief  Demonstrate how to reduce a first order LTI system using a low-rank square root balanced truncation (BT-LRSRM).
 * @author @koehlerm
 * @sa mess_bt_lrsrm
 *
 * # Tutorial: Low-rank Square Root Balanced Trunctation First Order LTI to First Order LTI System
 *
 * This function demonstrates how to convert a first order linear time invariant (LTI) system
 * \f[
 * \begin{array}{ccc}
 * \dot{x} &= A x & + B u \\
 * y &= C x&
 *  \end{array}
 * \f]
 * to another reduced first order LTI system
 * using a low-rank square root balanced truncation.
 *
 * @snippet "tutorials/bt/tutorial_bt_lrsrm.c" CODE
 *
 */

///@cond
///[CODE]
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/mess.h"


int main ( int argc, char **argv){
    mess_matrix A,B,C;
    mess_matrix V,W;
    mess_options adiopt;
    mess_bt_options btopt;
    mess_bt_status  btstat;
    mess_dynsys lti, red;
    mess_int_t maxr = 0 ;
    double tol = 0, err1, err2;
    mess_vector omega,G1, G2;


    printf("Balanced Truncation Example for LTI Systems\n");
    printf("===========================================\n");


    /*-----------------------------------------------------------------------------
     *  Read Matrices
     *-----------------------------------------------------------------------------*/
    if ( argc != 6 ) {
        printf("usage: %s A.mtx B.mtx C.mtx tol maxr\n", argv[0]);
        return 0;
    }
    mess_init();

    mess_matrix_init(&A);
    mess_matrix_init(&B);
    mess_matrix_init(&C);

    mess_set_errorlevel(1);

    mess_matrix_read_formated(argv[1], A, MESS_CSR);
    mess_matrix_read_formated(argv[2], B, MESS_DENSE);
    mess_matrix_read_formated(argv[3], C, MESS_DENSE);

    tol  = atof(argv[4]);
    maxr = atoi(argv[5]);

    printf("Dimension: " MESS_PRINTF_INT "\n", A->rows);
    printf("Inputs:    " MESS_PRINTF_INT "\n", B->cols);
    printf("Outputs:   " MESS_PRINTF_INT "\n\n", C->rows);

    /*-----------------------------------------------------------------------------
     *  Setup BT
     *-----------------------------------------------------------------------------*/

    mess_dynsys_init(&lti);
    mess_dynsys_lti_copy(lti, A, B, C);

    mess_options_init(&adiopt);
    mess_bt_options_init(&btopt);
    mess_bt_status_init(&btstat);
    btopt->rdim =  maxr;
    btopt->tol = tol;

    adiopt->adi_output = 0;
    mess_matrix_init(&V);
    mess_matrix_init(&W);
    printf("BT options\n");
    mess_bt_options_print(btopt);

    /*-----------------------------------------------------------------------------
     *  Compute the Projection Matrices
     *-----------------------------------------------------------------------------*/
    mess_bt_lrsrm(lti,btopt, adiopt, V,W, btstat);

    printf("\nStatus of mess_bt_lrsrm:\n");
    mess_bt_status_print(btstat);
    printf("\n");


    /*-----------------------------------------------------------------------------
     *  Compute the Reduced Order Model
     *-----------------------------------------------------------------------------*/
    mess_dynsys_init(&red);
    mess_dynsys_project_to_lti(lti, V,W,red);

    /*-----------------------------------------------------------------------------
     *  Evaluate the Transferfunctions
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&omega,&G1,&G2);
    mess_vector_alloc(omega,100,MESS_REAL);
    mess_vector_alloc(G1, 100, MESS_REAL);
    mess_vector_alloc(G2, 100, MESS_REAL);

    mess_dynsys_evaltransfer(lti,-1,6,100,omega,G1,NULL);
    mess_dynsys_evaltransfer(red,-1,6,100,omega,G2,NULL);

    mess_vector_diffnorm(G1,G2,&err1);
    mess_vector_diffnorminf(G1,G2,&err2);
    printf("error discrete ||G1-G2||_2 = %lg\n", err1);
    printf("error discrete ||G1-G2||_infty = %lg\n", err2);


    /*-----------------------------------------------------------------------------
     *  Clean up
     *-----------------------------------------------------------------------------*/

    mess_vector_clear(&G1);
    mess_vector_clear(&G2);
    mess_vector_clear(&omega);

    mess_matrix_clear(&A);
    mess_matrix_clear(&B);
    mess_matrix_clear(&C);
    mess_matrix_clear(&V);
    mess_matrix_clear(&W);
    mess_dynsys_clear(&lti);
    mess_dynsys_clear(&red);
    mess_options_clear(&adiopt);
    mess_bt_options_clear(&btopt);
    mess_bt_status_clear(&btstat);
    return 0;
}
///[CODE]
///@endcond



