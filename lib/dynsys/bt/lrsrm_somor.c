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
 * @file lib/dynsys/bt/lrsrm_somor.c
 * @brief LR-SRM balanced truncation.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"

#ifdef _OPENMP
#include <omp.h>
#endif


/**
 * @brief Low-rank square root balanced truncation for second order LTI systems to first order LTI systems.
 * @param[in] sys   generalized input system
 * @param[in] btopt      input balanced truncation options
 * @param[in] adiopt     input options for the underlying ADI Solver
 * @param[out] V    right projection matrix
 * @param[out] W    left projection matrix
 * @param[out] status   status information
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_bt_lrsrm_somor function computes the left and the right projection matrix for
 * a balanced truncation reduced order model. \n
 * In order to get the reduced order model from the output, you have to call the
 * \ref mess_dynsys_project_2nd_to_1st function with the computed \f$ V \f$ and \f$ W \f$. \n
 * In the options structure the \b tol value is given to set the cut off
 * tolerance for the Hankel singular values. Also the \b rdim parameter is given to set the maximum
 * wanted dimension of the reduced oder model.\n
 * If you do not want to choose the reduced order in a adaptive way so can specify a \b chooseorder function
 * to decide it in your own way. The default choose order function is \ref mess_bt_chooseorder_default.
 *
 */
int mess_bt_lrsrm_somor(mess_dynsys sys, mess_bt_options btopt, mess_options adiopt, mess_matrix V, mess_matrix W, mess_bt_status status) {
    MSG_FNAME(__func__);
    mess_matrix M,G,K,B,Cp,Cv,ZB,ZC, ZCZB, C, B2;
    mess_matrix UC, UB,S, T;
    mess_vector SIGMA;
    double cuterror = 0;
    int ret = 0;
    mess_int_t k0 = 0;
    mess_bt_chooseorder_func chooseorder = mess_bt_chooseorder_default;
    mess_int_t maxr = 0;
    double tol = 0;
    double ts,te;
    double talls;
    mess_equation eqn;
    mess_operation_t old_op;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sys);
    mess_check_nullpointer(btopt);
    mess_check_nullpointer(adiopt);
    mess_check_nullpointer(W);
    mess_check_nullpointer(V);
    mess_check_nullpointer(status);

    if ( !MESS_IS_DYNSYS_2ND(sys) ){
        MSG_ERROR("Only for 2nd Order Systems\n");
        return MESS_ERROR_DYNSYS;
    }

    if ( btopt->chooseorder != NULL ) {
        MSG_INFO("using user supplied reduced order chooser.\n");
        chooseorder = btopt->chooseorder;
    }
    maxr = btopt->rdim;
    tol = btopt->tol;
    mess_check_positive(maxr);
    if ( tol < 0.0) {
        MSG_ERROR("The given tolerance must be at least zero.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    talls = mess_wtime();
    ret = mess_matrix_init(&ZB); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&ZC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&ZCZB); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);

    M=sys->M;
    G=sys->G;
    K=sys->K;
    B=sys->B;
    Cp=sys->Cp;
    Cv=sys->Cv;

    /*-----------------------------------------------------------------------------
     *  compute the solutions of
     *    AX+XA^T=-BB^T
     *  and
     *    A^TY+YA=-C^TC
     *-----------------------------------------------------------------------------*/
    if ( mess_error_level > 2 ){
        MSG_PRINT("%s - ADI Options\n", __func__);
        mess_options_print(adiopt);
    }
    ts = mess_wtime();

    old_op = adiopt->type;
    adiopt->type = MESS_OP_NONE;
    ret = mess_equation_init(&eqn);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_init);
    ret = mess_matrix_init(&B2);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_init );
    ret = mess_matrix_init(&C);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_init );
    ret = mess_matrix_alloc(C, B->rows, B->cols, B->nnz, MESS_DENSE,MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_zeros(C);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);
    ret = mess_matrix_cat(C, NULL, B, NULL, MESS_DENSE, B2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);

    ret = mess_equation_glyap_so2(eqn, adiopt, M,G,K, B2, 1e-8,1e8); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_glyap);
    ret = mess_lrcfadi_adi(eqn, adiopt, status->statB, ZB);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_lrcfadi_adi);
    mess_matrix_clear(&B2);

    adiopt->type = MESS_OP_TRANSPOSE;
    ret = mess_matrix_cat(Cp, Cv, NULL, NULL, MESS_DENSE, C);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);
    ret = mess_equation_set_rhs(eqn,adiopt, C);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_set_rhs);
    ret = mess_lrcfadi_adi(eqn, adiopt, status->statC, ZC);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_lrcfadi_adi);
    mess_matrix_clear(&C);
    mess_equation_clear(&eqn);
    adiopt->type = old_op;


    /* ret =mess_lrcfadi_lrcfadidual_2nd(M,G,K,B,Cp,Cv,ZB,ZC,adiopt, status->statB,status->statC);
       FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_lrcfadi_lrcfadidual_2nd); */
    te = mess_wtime();
    status->time_lyap=te-ts;

    /* if ( mess_error_level > 2) {
       MSG_PRINT("%s - ADI status for ZB:\n",__func__);
       mess_status_print(status->statB);
       MSG_PRINT("%s - ADI status for ZC:\n",__func__);
       mess_status_print(status->statC);
       } */

    // mess_matrix_write("ZB.mtx", ZB, MESS_FILE_TEXT);
    // mess_matrix_write("ZC.mtx", ZC, MESS_FILE_TEXT);



    /*-----------------------------------------------------------------------------
     *  compute the SVD
     *-----------------------------------------------------------------------------*/
    ts = mess_wtime();
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, ZC, MESS_OP_NONE, ZB, ZCZB);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

    ret = mess_matrix_init(&UC);                                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&UB);                                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&S);                                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);

    ret = mess_vector_init(&SIGMA);                                                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
    ret = mess_vector_alloc(SIGMA, (ZCZB->rows < ZCZB->cols)? ZCZB->rows:ZCZB->cols, MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_alloc);

    ret = mess_eigen_svd_econ(ZCZB,SIGMA, UC, UB);                                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_eigen_svd_econ);


    /*-----------------------------------------------------------------------------
     *  compute the reduced order
     *-----------------------------------------------------------------------------*/
    // mess_vector_printshort(SIGMA);
    if ( maxr  > SIGMA->dim) {  MSG_INFO("setting maximum reduced order to " MESS_PRINTF_INT "\n", SIGMA->dim); maxr = SIGMA->dim; }
    ret  = chooseorder(SIGMA, &tol, &maxr,2*M->rows,btopt->chooseorder_aux);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), chooseorder);
    if ( maxr  > SIGMA->dim) {  MSG_INFO("reduced order is too big, setting reduced order to " MESS_PRINTF_INT "\n", SIGMA->dim); maxr = SIGMA->dim; }
    MSG_INFO("maxr = " MESS_PRINTF_INT "\n", maxr);
    status->rdim=maxr;

    /*-----------------------------------------------------------------------------
     *  estimate the error
     *-----------------------------------------------------------------------------*/
    cuterror = 0;
    mess_int_t i = 0 ;
    for ( i = maxr; i < SIGMA->dim;i++ ) cuterror+=SIGMA->values[i];
    MSG_INFO("2*cuterror = %lg\n", 2*cuterror);
    status->esterror = 2 * cuterror;
    k0=maxr;


    /*-----------------------------------------------------------------------------
     *  calc the new sigma
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_alloc(S,k0,k0,k0*k0,MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    for (i = 0 ; i < k0; i++){
        S->values[i+i*k0]=1.0/sqrt(SIGMA->values[i]);
    }


    /*-----------------------------------------------------------------------------
     *  cut of UC and UB
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_resize(UC, UC->rows, k0);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_resize);
    ret = mess_matrix_resize(UB, UB->rows, k0);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_resize);


    /*-----------------------------------------------------------------------------
     *  compute V and W
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&T); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_multiply(MESS_OP_NONE, UB, MESS_OP_NONE, S, T);FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, ZB, MESS_OP_NONE, T, V);FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, UC, MESS_OP_NONE, S, T);FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, ZC, MESS_OP_NONE, T, W);FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    te = mess_wtime();
    status->time =  te -talls;
    status->time_VW = te - ts;


    mess_matrix_clear(&S);
    mess_matrix_clear(&ZB);
    mess_matrix_clear(&ZC);
    mess_matrix_clear(&ZCZB);
    mess_matrix_clear(&UC);
    mess_matrix_clear(&UB);
    mess_vector_clear(&SIGMA);
    mess_matrix_clear(&T);

    return 0;
}

