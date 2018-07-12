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
 * @file tutorials/bt/tutorial_bt_gglrsrm.c
 * @brief Demonstrate how to reduce a first order GLTI system using a generalized low-rank square root balanced truncation (BT-GGLRSRM).
 * @author @koehlerm
 *
 * @sa mess_bt_gglrsrm
 *
 * # Tutorial: Reduce a First Order Generalized LTI System Using a Low-rank Square Root Balanced Trunction Method.
 *
 * This function demonstrates how to convert a first order generalized linear time invariant (GLTI) system
 * \f[
 * \begin{array}{ccc}
 * E \dot{x} &= A x & + B u \\
 * y &= C x&
 *  \end{array}
 * \f]
 * to another reduced first order GLTI system using a generalized low-rank square root balanced truncation.
 *
 * @snippet "tutorials/bt/tutorial_bt_gglrsrm.c" CODE
 *
 */

///@cond
///[CODE]
#include <stdlib.h>
#include "mess/mess.h"

#define NSAMPLE 25
#define S_LOW -2
#define S_HIGH 6


int main ( int argc, char **argv){
    mess_matrix A,B,C,E;
    mess_matrix V,W;
    mess_options adiopt;
    mess_bt_options btopt;
    mess_bt_status  btstat;
    mess_dynsys lti, red;
    mess_int_t maxr = 0 ;
    double tol = 0;
    mess_vector omega,G1, G2;

    printf("Balanced Truncation Example for generalized LTI Systems\n");
    printf("=======================================================\n");


    /*-----------------------------------------------------------------------------
     *  Read and init matrices
     *-----------------------------------------------------------------------------*/
    if ( argc != 7 ) {
        printf("usage: %s A.mtx B.mtx C.mtx E.mtx tol maxr\n", argv[0]);
        return 0;
    }

    MESS_INIT_MATRICES(&A,&B,&C,&E);

    mess_matrix_read_formated(argv[1], A, MESS_CSR);
    mess_matrix_read_formated(argv[2], B, MESS_DENSE);
    mess_matrix_read_formated(argv[3], C, MESS_DENSE);
    mess_matrix_read_formated(argv[4], E, MESS_CSR);
    tol  = atof(argv[5]);
    maxr = atoi(argv[6]);


    mess_dynsys_init(&lti);
    mess_set_errorlevel(1);


    printf("Dimension: " MESS_PRINTF_INT "\n", A->rows);
    printf("Inputs:    " MESS_PRINTF_INT "\n", B->cols);
    printf("Outputs:   " MESS_PRINTF_INT "\n", C->rows);

    /*-----------------------------------------------------------------------------
     *  Setup BT
     *-----------------------------------------------------------------------------*/
    mess_dynsys_glti_copy(lti, A, B, C, E);

    mess_options_init(&adiopt);
    mess_bt_options_init(&btopt);
    mess_bt_status_init(&btstat);
    printf("maxr_: " MESS_PRINTF_INT "\n", maxr);
    btopt->rdim =  maxr;
    btopt->tol = tol;

    adiopt->adi_output = 0;

    MESS_INIT_MATRICES(&V,&W);
    printf("ADI options\n");
    mess_options_print(adiopt);
    printf("\nBT options\n");
    mess_bt_options_print(btopt);


    /*-----------------------------------------------------------------------------
     *  Compute the Projection Matrices
     *-----------------------------------------------------------------------------*/
    mess_bt_gglrsrm(lti,btopt, adiopt, V,W, btstat);

    printf("\nStatus of mess_bt_lrsrm:\n");
    mess_bt_status_print(btstat);
    printf("\n");

    /*-----------------------------------------------------------------------------
     *  Compute the reduced order model
     *-----------------------------------------------------------------------------*/
    mess_dynsys_init(&red);
    mess_dynsys_project_to_glti(lti, V,W,red);

    /*-----------------------------------------------------------------------------
     *  Evaluate the Transferfunction
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&omega,&G1,&G2);
    mess_vector_alloc(omega,100,MESS_REAL);
    mess_vector_alloc(G1, 100, MESS_REAL);
    mess_vector_alloc(G2, 100, MESS_REAL);

    int i;
    mess_matrix *Go,*Gr;
    Gr = malloc (sizeof(mess_matrix)*NSAMPLE);
    Go = malloc (sizeof(mess_matrix)*NSAMPLE);
    for ( i = 0 ; i < NSAMPLE; i++){
        mess_matrix_init(&Gr[i]);
        mess_matrix_init(&Go[i]);
    }
    mess_dynsys_evaltransfer(lti,S_LOW,S_HIGH,NSAMPLE,omega,G1,Go);
    mess_dynsys_evaltransfer(red,S_LOW,S_HIGH,NSAMPLE,omega,G2,Gr);


    double err2, rel2, errinf, relinf;
    mess_vector diff;
    MESS_INIT_VECTORS(&diff);
    mess_vector_alloc(diff, NSAMPLE, MESS_REAL);
    mess_dynsys_transfer_diff(NSAMPLE, Go, Gr,&err2, &rel2, &errinf, &relinf, diff);
    for ( i = 0; i < NSAMPLE; i++){
        mess_matrix_clear(&Gr[i]);
        mess_matrix_clear(&(Go[i]));
    }
    free(Go); Go=NULL;
    free(Gr); Gr=NULL;

    printf("error discrete ||G1-G2||_2           = %lg\n", err2);
    printf("error discrete ||G1-G2||_2 rel.      = %lg\n", rel2);
    printf("error discrete ||G1-G2||_infty       = %lg\n", errinf);
    printf("error discrete ||G1-G2||_infty rel.  = %lg\n", relinf);

    if(err2 > 1e-4)  { printf("FAILED:error discrete ||G1-G2||_2 to large.\n"); return 1; }
    if(rel2 > 1e-8)  { printf("FAILED:error discrete ||G1-G2||_2_rel to large.\n"); return 1; }
    if(errinf > 1e-5){ printf("FAILED:error discrete ||G1-G2||_2_infty to large.\n"); return 1; }
    if(relinf > 1e-8){ printf("FAILED:error discrete ||G1-G2||_2_infty rel. to large.\n"); return 1; }

    /*-----------------------------------------------------------------------------
     *  Clean up
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_VECTORS(&G1,&G2,&omega,&diff);
    MESS_CLEAR_MATRICES(&A,&B,&C,&E,&V,&W);

    mess_dynsys_clear(&lti);
    mess_dynsys_clear(&red);
    mess_options_clear(&adiopt);
    mess_bt_options_clear(&btopt);
    mess_bt_status_clear(&btstat);
    return 0;
}
///[CODE]
///@endcond

