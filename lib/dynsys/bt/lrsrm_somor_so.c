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
 * @file lib/dynsys/bt/lrsrm_somor_so.c
 * @brief LR-SRM balanced truncation for second order LTI systems to second order LTI systems.
 * @author @koehlerm
 *
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
 * @brief Low-rank square root balanced truncation for second order LTI systems to second order LTI systems.
 * @param[in] sys   generalized input system
 * @param[in] btopt      input balanced truncation options
 * @param[in] adiopt     input options for the underlying ADI Solver
 * @param[out] V    right projection matrix
 * @param[out] W    left projection matrix
 * @param[out] status   status information
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_bt_lrsrm_somor_so function computes the left and the right projection matrix for
 * a balanced truncation reduced order model.\n
 * In order to get the reduced order model from the output, you have to call the
 * \ref mess_dynsys_project_2nd_to_2nd function with the computed \f$ V \f$ and \f$ W \f$. \n
 * In the options structure the \b tol value is given to set the cut off
 * tolerance for the Hankel singular values. Also the \b rdim parameter is given to set the maximum
 * wanted dimension of the reduced oder model. \n
 * If you do not want to choose the reduced order in a adaptive way so can specify a \b chooseorder function
 * to decide it in your own way. The default choose order function is \ref mess_bt_chooseorder_minreal.\n
 *
 */
int mess_bt_lrsrm_somor_so(mess_dynsys sys, mess_bt_options btopt, mess_options adiopt, mess_matrix V, mess_matrix W, mess_bt_status status) {
    MSG_FNAME(__func__);
    mess_matrix M,G,K,B,Cp,Cv,ZB,ZC, ZCZB, ZBr, ZCr, ZBM, C, B2;
    mess_matrix UC, UB,S, T;
    mess_vector SIGMA;
    // double cuterror = 0;
    int ret = 0;
    mess_int_t k0 = 0;
    mess_bt_chooseorder_func chooseorder = mess_bt_chooseorder_minreal;
    mess_int_t maxr = 0;
    double tol = 0;
    double ts,te;
    double talls;
    unsigned short meth;
    mess_int_t dim;
    mess_int_t i;
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
    meth = btopt->so_type;
    if ( meth < 1 || meth > 4 ) {
        MSG_ERROR("The type (%s) of the wanted second order system is unknown.\n", mess_bt_getsotypestr(meth));
        return MESS_ERROR_ARGUMENTS;
    }

    talls = mess_wtime();
    ret = mess_matrix_init(&ZB); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&ZC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&ZCZB); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&ZBr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&ZBM); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&ZCr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);

    M=sys->M;
    G=sys->G;
    K=sys->K;
    B=sys->B;
    Cp=sys->Cp;
    Cv=sys->Cv;
    dim  = M->rows;

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

    /*-----------------------------------------------------------------------------
     *  get the ZBr and ZCr from the wanted second order system.
     *-----------------------------------------------------------------------------*/
    switch(meth) {
        case MESS_BT_POSITION:
            MSG_INFO("get postion gramians\n");
            ret = mess_matrix_rowsub(ZB,0,dim-1,ZBr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            ret = mess_matrix_rowsub(ZC,0,dim-1,ZCr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            break;
        case MESS_BT_VELOCITY:
            MSG_INFO("get velocity gramians\n");
            ret = mess_matrix_rowsub(ZB,dim,2*dim-1,ZBr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            ret = mess_matrix_rowsub(ZC,dim,2*dim-1,ZCr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            break;
        case MESS_BT_POSITION_VELOCITY:
            MSG_INFO("get postion velocity gramians\n");
            ret = mess_matrix_rowsub(ZB,0,dim-1,ZBr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            ret = mess_matrix_rowsub(ZC,dim,2*dim-1,ZCr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            break;
        case MESS_BT_VELOCITY_POSITION:
            MSG_INFO("get velocity position gramians\n");
            ret = mess_matrix_rowsub(ZB,dim,2*dim-1,ZBr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            ret = mess_matrix_rowsub(ZC,0,dim-1,ZCr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            break;
        default:
            break;
            // normaly this should be catched before.
    }
    mess_matrix_clear(&ZB);
    mess_matrix_clear(&ZC);

    /*-----------------------------------------------------------------------------
     *  compute the SVD
     *-----------------------------------------------------------------------------*/
    ts = mess_wtime();
    ret = mess_matrix_multiply(MESS_OP_NONE, M, MESS_OP_NONE, ZBr, ZBM); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, ZCr, MESS_OP_NONE, ZBM, ZCZB); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    mess_matrix_clear(&ZBM);

    MESS_INIT_MATRICES(&UC,&UB,&S);
    ret = mess_vector_init(&SIGMA);                                                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
    ret = mess_vector_alloc(SIGMA, (ZCZB->rows < ZCZB->cols)? ZCZB->rows:ZCZB->cols, MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_init);

    ret = mess_eigen_svd_econ(ZCZB,SIGMA, UC, UB);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_eigen_svd_econ);


    /*-----------------------------------------------------------------------------
     *  compute the reduced order
     *-----------------------------------------------------------------------------*/
    // mess_vector_printshort(SIGMA);
    if ( maxr  > SIGMA->dim) {  MSG_INFO("setting maximum reduced order to " MESS_PRINTF_INT "\n", SIGMA->dim); maxr = SIGMA->dim; }
    ret  = chooseorder(SIGMA, &tol, &maxr,2*M->rows,btopt->chooseorder_aux);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), chooseorder);
    if ( maxr  > SIGMA->dim) {  MSG_INFO("reduced order is too big, setting reduced order to " MESS_PRINTF_INT "\n", SIGMA->dim); maxr = SIGMA->dim; }
    MSG_INFO("maxr = " MESS_PRINTF_INT "\n", maxr);
    status->rdim=maxr;

    /*-----------------------------------------------------------------------------
     *  estimate the error, not for SO MOR
     *-----------------------------------------------------------------------------*/
    /* cuterror = 0;
       mess_int_t i = 0 ;
       for ( i = maxr; i < SIGMA->dim;i++ ) cuterror+=SIGMA->values[i];
       MSG_INFO("2*cuterror = %lg\n", 2*cuterror);
       status->esterror = 2 * cuterror; */
    k0=maxr;


    /*-----------------------------------------------------------------------------
     *  calc the new sigma
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_alloc(S,k0,k0,k0*k0,MESS_DENSE, MESS_REAL);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    for (i = 0 ; i < k0; i++){
        S->values[i+i*k0]=1.0/sqrt(SIGMA->values[i]);
    }


    /*-----------------------------------------------------------------------------
     *  cut of UC and UB
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_resize(UC, UC->rows, k0);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_resize);
    ret = mess_matrix_resize(UB, UB->rows, k0);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_resize);


    /*-----------------------------------------------------------------------------
     *  compute V and W
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&T); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_multiply(MESS_OP_NONE, UB, MESS_OP_NONE, S, T);FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, ZBr, MESS_OP_NONE, T, V);FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, UC, MESS_OP_NONE, S, T);FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, ZCr, MESS_OP_NONE, T, W);FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    te = mess_wtime();
    status->time =  te -talls;
    status->time_VW = te - ts;


    mess_matrix_clear(&S);
    mess_matrix_clear(&ZBr);
    mess_matrix_clear(&ZCr);
    mess_matrix_clear(&ZCZB);
    mess_matrix_clear(&UC);
    mess_matrix_clear(&UB);
    mess_vector_clear(&SIGMA);
    mess_matrix_clear(&T);

    return 0;
}

/**
 * @brief Low-rank square root balanced truncation for second order LTI systems to second order LTI systems (Gramian version).
 * @param[in] sys   second order input system
 * @param[in] ZB     input controllability gramian
 * @param[in] ZC     input observability gramian
 * @param[in] btopt      input balanced truncation options
 * @param[out] V    right projection matrix
 * @param[out] W    left projection matrix
 * @param[out] status   status information
 * @return zero on succes or a non-zero error value
 *
 * The @ref mess_bt_lrsrm_somor_so function computes the left and the right projection matrix for
 * a balanced truncation reduced order model. In order to get the reduced order model
 * from the output, you have to call the \ref mess_dynsys_project_2nd_to_2nd function with the
 * computed \f$ V \f$ and \f$ W \f$. \n
 * In the options structure the \b tol value is given to set the cut off
 * tolerance for the Hankel singular values. Also the \b rdim parameter is given to set the maximum
 * wanted dimension of the reduced oder model. If you do not want to choose the reduced
 * order in a adaptive way so can specify a \b chooseorder function to decide it in your
 * own way. The default choose order function is \ref mess_bt_chooseorder_minreal.
 *
 */
int mess_bt_lrsrm_somor_so_gram(mess_dynsys sys, mess_matrix ZB, mess_matrix ZC, mess_bt_options btopt, mess_matrix V, mess_matrix W, mess_bt_status status) {
    MSG_FNAME(__func__);
    mess_matrix M,ZCZB, ZBr, ZCr, ZBM;
    mess_matrix UC, UB,S, T;
    mess_vector SIGMA;
    // double cuterror = 0;
    int ret = 0;
    mess_int_t k0 = 0;
    mess_bt_chooseorder_func chooseorder = mess_bt_chooseorder_minreal;
    mess_int_t maxr = 0;
    double tol = 0;
    double ts,te;
    double talls;
    unsigned short meth;
    mess_int_t dim;
    mess_int_t i;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sys);
    mess_check_nullpointer(ZB);
    mess_check_nullpointer(ZC);
    mess_check_nullpointer(btopt);
    mess_check_nullpointer(W);
    mess_check_nullpointer(V);
    mess_check_nullpointer(status);

    if ( !MESS_IS_DYNSYS_2ND(sys) ){
        MSG_ERROR("Only for 2nd Order Systems\n");
        return MESS_ERROR_DYNSYS;
    }
    if ( ZB->rows != sys->M->rows*2) {
        MSG_ERROR("ZB has the wrong number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( ZC->rows != sys->M->rows*2) {
        MSG_ERROR("ZC has the wrong number of rows.\n");
        return MESS_ERROR_DIMENSION;
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
    meth = btopt->so_type;
    if ( meth < 1 || meth > 4 ) {
        MSG_ERROR("The type (%s) of the wanted second order system is unknown.\n", mess_bt_getsotypestr(meth));
        return MESS_ERROR_ARGUMENTS;
    }

    talls = mess_wtime();
    ret = mess_matrix_init(&ZCZB); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&ZBr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&ZBM); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&ZCr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);

    M=sys->M;
    dim  = M->rows;

    // mess_matrix_write("ZB.mtx", ZB, MESS_FILE_TEXT);
    // mess_matrix_write("ZC.mtx", ZC, MESS_FILE_TEXT);



    /*-----------------------------------------------------------------------------
     *  get the ZBr and ZCr from the wanted second order system.
     *-----------------------------------------------------------------------------*/
    switch(meth) {
        case MESS_BT_POSITION:
            MSG_INFO("get postion gramians\n");
            ret = mess_matrix_rowsub(ZB,0,dim-1,ZBr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            ret = mess_matrix_rowsub(ZC,0,dim-1,ZCr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            break;
        case MESS_BT_VELOCITY:
            MSG_INFO("get velocity gramians\n");
            ret = mess_matrix_rowsub(ZB,dim,2*dim-1,ZBr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            ret = mess_matrix_rowsub(ZC,dim,2*dim-1,ZCr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            break;
        case MESS_BT_POSITION_VELOCITY:
            MSG_INFO("get postion velocity gramians\n");
            ret = mess_matrix_rowsub(ZB,0,dim-1,ZBr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            ret = mess_matrix_rowsub(ZC,dim,2*dim-1,ZCr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            break;
        case MESS_BT_VELOCITY_POSITION:
            MSG_INFO("get velocity position gramians\n");
            ret = mess_matrix_rowsub(ZB,dim,2*dim-1,ZBr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            ret = mess_matrix_rowsub(ZC,0,dim-1,ZCr); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            break;
        default:
            break;
            // normaly this should be catched before.
    }
    // mess_matrix_clear(&ZB);
    // mess_matrix_clear(&ZC);

    /*-----------------------------------------------------------------------------
     *  compute the SVD
     *-----------------------------------------------------------------------------*/
    ts = mess_wtime();
    ret = mess_matrix_multiply(MESS_OP_NONE, M, MESS_OP_NONE, ZBr, ZBM); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, ZCr, MESS_OP_NONE, ZBM, ZCZB); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    mess_matrix_clear(&ZBM);

    ret = mess_matrix_init(&UC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&UB); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&S); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);

    ret = mess_vector_init(&SIGMA);                                                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
    ret = mess_vector_alloc(SIGMA, (ZCZB->rows < ZCZB->cols)? ZCZB->rows:ZCZB->cols, MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_init);

    ret = mess_eigen_svd_econ(ZCZB,SIGMA, UC, UB);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_eigen_svd_econ);


    /*-----------------------------------------------------------------------------
     *  compute the reduced order
     *-----------------------------------------------------------------------------*/
    // mess_vector_printshort(SIGMA);
    if ( maxr  > SIGMA->dim) {  MSG_INFO("setting maximum reduced order to " MESS_PRINTF_INT "\n", SIGMA->dim); maxr = SIGMA->dim; }
    ret  = chooseorder(SIGMA, &tol, &maxr,2*M->rows,btopt->chooseorder_aux);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), "choose order failed.");
    if ( maxr  > SIGMA->dim) {  MSG_INFO("reduced order is too big, setting reduced order to " MESS_PRINTF_INT "\n", SIGMA->dim); maxr = SIGMA->dim; }
    MSG_INFO("maxr = " MESS_PRINTF_INT "\n", maxr);
    status->rdim=maxr;

    /*-----------------------------------------------------------------------------
     *  estimate the error, not for SO MOR
     *-----------------------------------------------------------------------------*/
    /* cuterror = 0;
       mess_int_t i = 0 ;
       for ( i = maxr; i < SIGMA->dim;i++ ) cuterror+=SIGMA->values[i];
       MSG_INFO("2*cuterror = %lg\n", 2*cuterror);
       status->esterror = 2 * cuterror; */
    k0=maxr;


    /*-----------------------------------------------------------------------------
     *  calc the new sigma
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_alloc(S,k0,k0,k0*k0,MESS_DENSE, MESS_REAL);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    for (i = 0 ; i < k0; i++){
        S->values[i+i*k0]=1.0/sqrt(SIGMA->values[i]);
    }


    /*-----------------------------------------------------------------------------
     *  cut of UC and UB
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_resize(UC, UC->rows, k0);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_resize);
    ret = mess_matrix_resize(UB, UB->rows, k0);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_resize);


    /*-----------------------------------------------------------------------------
     *  compute V and W
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&T); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_multiply(MESS_OP_NONE, UB, MESS_OP_NONE, S, T);FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, ZBr, MESS_OP_NONE, T, V);FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, UC, MESS_OP_NONE, S, T);FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, ZCr, MESS_OP_NONE, T, W);FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    te = mess_wtime();
    status->time =  te -talls;
    status->time_VW = te - ts;


    mess_matrix_clear(&S);
    mess_matrix_clear(&ZBr);
    mess_matrix_clear(&ZCr);
    mess_matrix_clear(&ZCZB);
    mess_matrix_clear(&UC);
    mess_matrix_clear(&UB);
    mess_vector_clear(&SIGMA);
    mess_matrix_clear(&T);

    return 0;
}

