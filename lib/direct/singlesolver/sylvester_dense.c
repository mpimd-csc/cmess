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
 * @file lib/direct/singlesolver/sylvester_dense.c
 * @brief Generate a solver for dense Sylvester eqautions based on  @lapack
 * @c  dtrsyl / @c ztrsyl / @c dtgsyl / @c ztgsyl.
 * @author @koehlerm
 * @author @dykstra
 * @author @mbehr
 *
 * This code implements a sparse-dense Sylvester equation solver for the following Sylvester equations:
 * <ul>
 * <li> Standard:  \f$AX+XH+M=0\f$, \f$A^TX+XH^T+M=0\f$, \f$A^HX+XH^H+M=0\f$
 * <li> Generalized Sylvester: \f$AXF+EXH+M=0\f$, \f$A^TXF^T+E^TXH^T+M=0\f$, \f$A^HXF^H+E^HXH^H+M=0\f$
 * </ul>
 * with \f$ A,F,E,H \f$ small and dense matrices.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"



/**
 * @internal
 *
 * @attention Internal use only.
 */
typedef struct _sylv_solver_d_st {
    mess_matrix Ahat;                   /** Schur Form of Matrix \f$A\f$ of Sylvester Equation.*/
    mess_matrix QA;                     /** Schur Transformation Matrix of Matrix \f$A\f$.*/
    mess_matrix Ehat;                   /** Schur Form of Matrix \f$E\f$ of Sylvester Equation.*/
    mess_matrix QE;                     /** Schur Transformation Matrix of Matrix \f$E\f$.*/
    mess_matrix Fhat;                   /** Schur Form of Matrix \f$F\f$ of Sylvester Equation.*/
    mess_matrix QF;                     /** Schur Transformation Matrix of Matrix \f$F\f$.*/
    mess_matrix Hhat;                   /** Schur Form of Matrix \f$H\f$ of Sylvester Equation.*/
    mess_matrix QH;                     /** Schur Transformation Matrix of Matrix \f$H\f$.*/
    mess_int_t rowsX;                   /** Rows of Solution \f$ X\f$ of Sylvester Equation.*/
    mess_int_t colsX;                   /** Columns of Solution \f$ X\f$ of Sylvester Equation.*/
    int isreal;                         /** @c 1 if all Matrices of Sylvester Equation are real, @c 0 otherwise.*/
} _sylv_solver_d;




/**
 * @internal
 * @brief Clean up handler for @ref mess_direct_create_sylvester_dense.
 * @param[in] solver  input pointer to the solver data
 * @return always zero
 *
 * The @ref sylvester_d_clear function is the clean up function.
 *
 * @attention Internal use only.
 */
static int sylvester_d_clear(void *solver){
    _sylv_solver_d * sol = (_sylv_solver_d*) solver;
    if ( sol != NULL) {
        mess_matrix_clear(&(sol->Ahat));
        mess_matrix_clear(&(sol->QA));
        mess_matrix_clear(&(sol->Ehat));
        mess_matrix_clear(&(sol->QE));
        mess_matrix_clear(&(sol->Fhat));
        mess_matrix_clear(&(sol->QF));
        mess_matrix_clear(&(sol->Hhat));
        mess_matrix_clear(&(sol->QH));
        mess_free( sol );
    }
    return 0;
}


/**
 * @internal
 * @brief Solve \f$ AX+XH+M=0 \f$.
 * @param[in] datain  input solver data
 * @param[in] M  input right hand sides
 * @param[out] X  output solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref sylvester_d_solvem_standard function solves dense Sylvester equations:
 *  <ul>
 *  <li> \f$AX+XH+M=0\f$
 *  </ul>
 *
 * @attention Internal use only.
 */
static int sylvester_d_solvem_standard(void *datain, mess_matrix M, mess_matrix X){
    MSG_FNAME(__func__);
    _sylv_solver_d * data = (_sylv_solver_d *)datain;
    int ret = 0;
    mess_matrix Mtilde, Xtilde, Tmp1, Tmp2;
    double scale1 = 1, scale2 = 1;
    mess_int_t info = 0, sgn=1;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(X);
    mess_check_nullpointer(M);
    mess_check_dense(M);
    mess_check_real_or_complex(M);

    if ( M->rows != data->rowsX) {
        MSG_ERROR("number of rows doesn't match.\n");
        return MESS_ERROR_DIMENSION;
    }

    if ( M->cols != data->colsX){
        MSG_ERROR("number of columns doesn't match.\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  init and alloc some additional matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&Mtilde,&Xtilde,&Tmp1,&Tmp2);

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, data->QA, MESS_OP_NONE, M, X);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, X, MESS_OP_NONE, data->QH, Mtilde);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);


    /*-----------------------------------------------------------------------------
     *  standard dense Sylvester equation A X + X H = -M
     *-----------------------------------------------------------------------------*/
    if(data->isreal){
        if(MESS_IS_REAL(Mtilde)){
            F77_GLOBAL(dtrsyl,DTRSYL)("N","N",&sgn, &(data->rowsX), &(data->colsX), data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld),
                                                                                    Mtilde->values, &(Mtilde->ld), &scale1, &info);
            ret = mess_matrix_scale(-1.0/scale1, Mtilde);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
        }else{
            ret = mess_matrix_realpart(Mtilde,Tmp1);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realpart);
            ret = mess_matrix_imagpart(Mtilde,Tmp2);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_imagpart);
            F77_GLOBAL(dtrsyl,DTRSYL)("N","N",&sgn, &(data->rowsX), &(data->colsX), data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld),
                                                                                    Tmp1->values, &(Tmp1->ld), &scale1, &info);
            F77_GLOBAL(dtrsyl,DTRSYL)("N","N",&sgn, &(data->rowsX), &(data->colsX), data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld),
                                                                                    Tmp2->values, &(Tmp2->ld), &scale2, &info);
            ret = mess_matrix_copy(Tmp2, Mtilde);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            ret = mess_matrix_addc(-1.0/scale1, Tmp1, -1.0*I / scale2, Mtilde);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        }
    }else{
        ret = mess_matrix_tocomplex(Mtilde);                                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        F77_GLOBAL(ztrsyl,ZTRSYL)("N","N",&sgn, &(data->rowsX), &(data->colsX), data->Ahat->values_cpx, &(data->Ahat->ld), data->Hhat->values_cpx, &(data->Hhat->ld),
                                                                                Mtilde->values_cpx, &(Mtilde->ld), &scale1, &info);
        ret = mess_matrix_scale(-1.0/scale1, Mtilde);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
    }

    /*-----------------------------------------------------------------------------
     *  check info value
     *-----------------------------------------------------------------------------*/
    if ( info < 0 ) {
        MSG_ERROR("DTRSYL / ZTRSYL / DTGSYL / ZTGSYL returned with error " MESS_PRINTF_INT " \n", info);
        ret = MESS_ERROR_LAPACK;
        goto clean;
    }

    if ( info == 1 ) {
        MSG_WARN("Pertubed values were used to solve the Sylvester equation");
    }

    if ( (scale1 != 1) || (scale2 != 1)){
        MSG_WARN("DTRSYL / ZTRSYL / DTGSYL / ZTGSYL returned scale1 = %e and scale2 = %e to avoid overflow.\n", scale1, scale2);
    }

    /*-----------------------------------------------------------------------------
     *  post processing
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_multiply(MESS_OP_NONE, data->QA, MESS_OP_NONE, Mtilde, Xtilde);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, Xtilde, MESS_OP_HERMITIAN, data->QH , X);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
clean:
    MESS_CLEAR_MATRICES(&Mtilde,&Xtilde,&Tmp1, &Tmp2);
    return ret;

}


/**
 * @internal
 * @brief Solve \f$ A^HX+XH^H+M=0 \f$.
 * @param[in] datain  input solver data
 * @param[in] M  input right hand sides
 * @param[out] X  output solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref sylvester_d_solvemh_standard function solves dense Sylvester equations:
 *  <ul>
 *  <li> \f$A^HX+XH^H+M=0\f$
 *  </ul>
 *
 * @attention Internal use only.
 */
static int sylvester_d_solvemh_standard(void *datain, mess_matrix M, mess_matrix X){

    MSG_FNAME(__func__);
    _sylv_solver_d * data = (_sylv_solver_d *)datain;
    int ret = 0;
    mess_matrix Mtilde, Xtilde, Tmp1, Tmp2;
    double scale1 = 1, scale2 = 1;
    mess_int_t info = 0, sgn=1;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(X);
    mess_check_nullpointer(M);
    mess_check_dense(M);
    mess_check_real_or_complex(M);

    if ( M->rows != data->rowsX) {
        MSG_ERROR("number of rows doesn't match.\n");
        return MESS_ERROR_DIMENSION;
    }

    if ( M->cols != data->colsX){
        MSG_ERROR("number of columns doesn't match.\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  init and alloc some additional matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&Mtilde,&Xtilde,&Tmp1,&Tmp2);

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, data->QA, MESS_OP_NONE, M, X);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, X, MESS_OP_NONE, data->QH, Mtilde);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);


    /*-----------------------------------------------------------------------------
     *  standard dense Sylvester equation A^H X + X H^H = -M
     *-----------------------------------------------------------------------------*/
    if(data->isreal){
        if(MESS_IS_REAL(Mtilde)){
            F77_GLOBAL(dtrsyl,DTRSYL)("C","C",&sgn, &(data->rowsX), &(data->colsX), data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Mtilde->values, &(Mtilde->ld), &scale1, &info);
            ret = mess_matrix_scale(-1.0/scale1, Mtilde);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
        }else{
            ret = mess_matrix_realpart(Mtilde,Tmp1);                                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realpart);
            ret = mess_matrix_imagpart(Mtilde,Tmp2);                                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_imagpart);
            F77_GLOBAL(dtrsyl,DTRSYL)("C","C",&sgn, &(data->rowsX), &(data->colsX), data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Tmp1->values, &(Tmp1->ld), &scale1, &info);
            F77_GLOBAL(dtrsyl,DTRSYL)("C","C",&sgn, &(data->rowsX), &(data->colsX), data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Tmp2->values, &(Tmp2->ld), &scale2, &info);
            ret = mess_matrix_copy(Tmp2, Mtilde);                                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            ret = mess_matrix_addc(-1.0/scale1, Tmp1, -1.0*I / scale2, Mtilde);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
        }
    }else{
        ret = mess_matrix_tocomplex(Mtilde);                                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        F77_GLOBAL(ztrsyl,ZTRSYL)("C","C",&sgn, &(data->rowsX), &(data->colsX), data->Ahat->values_cpx, &(data->Ahat->ld), data->Hhat->values_cpx, &(data->Hhat->ld), Mtilde->values_cpx, &(Mtilde->ld), &scale1, &info);
        ret = mess_matrix_scale(-1.0/scale1, Mtilde);                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
    }

    /*-----------------------------------------------------------------------------
     *  check info value
     *-----------------------------------------------------------------------------*/
    if ( info < 0 ) {
        MSG_ERROR("DTRSYL / ZTRSYL / DTGSYL / ZTGSYL returned with error " MESS_PRINTF_INT " \n", info);
        ret = MESS_ERROR_LAPACK;
        goto clean;
    }

    if ( info == 1 ) {
        MSG_WARN("Pertubed values were used to solve the Sylvester equation");
    }

    if ( (scale1 != 1) || (scale2 != 1)){
        MSG_WARN("DTRSYL / ZTRSYL / DTGSYL / ZTGSYL returned scale1 = %e and scale2 = %e to avoid overflow.\n", scale1, scale2);
    }

    /*-----------------------------------------------------------------------------
     *  post processing
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_multiply(MESS_OP_NONE, data->QA, MESS_OP_NONE, Mtilde, Xtilde);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, Xtilde, MESS_OP_HERMITIAN, data->QH , X);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
clean:
    MESS_CLEAR_MATRICES(&Mtilde, &Xtilde, &Tmp1, &Tmp2);

    return ret;

}


/**
 * @internal
 * @brief Solve \f$ A^TX+XH^T+M=0 \f$.
 * @param[in] datain  input solver data
 * @param[in] M  input right hand sides
 * @param[out] X  output solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref sylvester_d_solvemt_standard function solves dense Sylvester equations:
 *  <ul>
 *  <li> \f$A^TX+XH^T+M=0\f$
 *  </ul>
 *
 * @attention Internal use only.
 */
static int sylvester_d_solvemt_standard(void *datain, mess_matrix M, mess_matrix X){
    MSG_FNAME(__func__);
    int ret = 0;
    ret = mess_matrix_conj(M);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_conj);
    ret = sylvester_d_solvemh_standard(datain, M, X);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), sylvester_d_solvemh_standard);
    ret = mess_matrix_conj(X);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_conj);
    ret = mess_matrix_conj(M);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_conj);
    return ret;
}




/**
 * @internal
 * @brief Solve \f$ AXF+EXH+M=0 \f$.
 * @param[in] datain  input solver data
 * @param[in] M  input right hand sides
 * @param[out] X  output solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref sylvester_d_solvem_generalized function solves dense Sylvester equations:
 *  <ul>
 *  <li> \f$AXF+EXH+M=0\f$
 *  </ul>
 *
 * @attention Internal use only.
 */
static int sylvester_d_solvem_generalized(void * datain, mess_matrix M, mess_matrix X){
    MSG_FNAME(__func__);
    _sylv_solver_d * data = (_sylv_solver_d *)datain;
    int ret = 0;
    mess_matrix Mtilde, Xtilde, Tmp1, Tmp2, Tmp3;
    double scale1 = 1, scale2 = 1, diff;
    mess_int_t info = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(X);
    mess_check_nullpointer(M);
    mess_check_dense(M);
    mess_check_real_or_complex(M);

    if ( M->rows != data->rowsX) {
        MSG_ERROR("number of rows doesn't match.\n");
        return MESS_ERROR_DIMENSION;
    }

    if ( M->cols != data->colsX){
        MSG_ERROR("number of columns doesn't match.\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  init and alloc some additional matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&Mtilde,&Xtilde,&Tmp1,&Tmp2,&Tmp3);

    /*-----------------------------------------------------------------------------
     *  Idea:
     *  Solve A X F + E X H + M = 0.
     *  (Schur decomposition) <==>
     *  Ahat * (QE^H * X * QH) * Fhat + Ehat * (QE^H * X * QH) * Hhat + QA^H M QF = 0
     *  <==>
     *  Ahat * Y * Fhat + Ehat* Y * Hhat + Mtilde = 0
     *  <==>
     *  Ahat L + R Hhat = Mtilde
     *  Ehat L - R Fhat = 0.
     *  Y = - Ehat^(-1) R.
     *-----------------------------------------------------------------------------*/

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, data->QA, MESS_OP_NONE, M, X);                                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, X, MESS_OP_NONE, data->QF, Mtilde);                                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_alloc(Tmp3, M->rows, M->cols, M->rows*M->cols, MESS_DENSE, data->isreal?MESS_REAL:MESS_COMPLEX);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    /*-----------------------------------------------------------------------------
     *  generalized dense Sylvester equation A X F + E X H = -M
     *  we first do a workspace query and then the call
     *-----------------------------------------------------------------------------*/
    double workspace, *work=NULL;
    mess_double_cpx_t workspace_cpx, *work_cpx=NULL;
    mess_int_t lwork, *iwork=NULL, ijob=0;
    lwork = -1;
    mess_try_alloc(iwork, mess_int_t *, sizeof(mess_int_t)*(data->rowsX + data->colsX + 6));

    mess_matrix_scale(-1.0, data->Hhat);
    if(data->isreal){
        if(MESS_IS_REAL(Mtilde)){
            F77_GLOBAL(dtgsyl,DTGSYL)("N", &ijob, &(data->rowsX), &(data->colsX),   data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Mtilde->values, &(Mtilde->ld),
                    data->Ehat->values, &(data->Ehat->ld), data->Fhat->values, &(data->Fhat->ld), Tmp3->values, &(Tmp3->ld), &diff, &scale1, &workspace, &lwork, iwork, &info);
            lwork = nearbyint(workspace);
            mess_try_alloc(work, double*, sizeof(double)*lwork);
            F77_GLOBAL(dtgsyl,DTGSYL)("N", &ijob, &(data->rowsX), &(data->colsX),   data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Mtilde->values, &(Mtilde->ld),
                    data->Ehat->values, &(data->Ehat->ld), data->Fhat->values, &(data->Fhat->ld), Tmp3->values, &(Tmp3->ld), &diff, &scale1, work, &lwork, iwork, &info);

            ret = mess_matrix_scale(-1.0/scale1, Tmp3);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
        }else{
            ret = mess_matrix_realpart(Mtilde,Tmp1);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realpart);
            F77_GLOBAL(dtgsyl,DTGSYL)("N", &ijob,   &(data->rowsX), &(data->colsX), data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Tmp1->values, &(Tmp1->ld),
                    data->Ehat->values, &(data->Ehat->ld), data->Fhat->values, &(data->Fhat->ld), Tmp3->values, &(Tmp3->ld), &diff, &scale1, &workspace, &lwork, iwork, &info);
            lwork = nearbyint(workspace);
            mess_try_alloc(work, double*, sizeof(double)*lwork);

            F77_GLOBAL(dtgsyl,DTGSYL)("N", &ijob,   &(data->rowsX), &(data->colsX), data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Tmp1->values, &(Tmp1->ld),
                    data->Ehat->values, &(data->Ehat->ld), data->Fhat->values, &(data->Fhat->ld), Tmp3->values, &(Tmp3->ld), &diff, &scale1, work, &lwork, iwork, &info);

            ret = mess_matrix_copy(Tmp3, Tmp1);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            ret = mess_matrix_imagpart(Mtilde,Tmp2);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_imagpart);
            ret = mess_matrix_zeros(Tmp3);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);

            F77_GLOBAL(dtgsyl,DTGSYL)("N", &ijob,   &(data->rowsX), &(data->colsX), data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Tmp2->values, &(Tmp2->ld),
                    data->Ehat->values, &(data->Ehat->ld), data->Fhat->values, &(data->Fhat->ld), Tmp3->values, &(Tmp3->ld), &diff, &scale2, work, &lwork, iwork, &info);

            ret = mess_matrix_addc(-1.0/scale1, Tmp1, (-1.0*I)/scale2, Tmp3);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        }
    }else{
        ret = mess_matrix_tocomplex(Mtilde);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        F77_GLOBAL(ztgsyl,ZTGSYL)("N", &ijob,   &(data->rowsX), &(data->colsX), data->Ahat->values_cpx, &(data->Ahat->ld), data->Hhat->values_cpx, &(data->Hhat->ld), Mtilde->values_cpx, &(Mtilde->ld),
                data->Ehat->values_cpx, &(data->Ehat->ld), data->Fhat->values_cpx, &(data->Fhat->ld), Tmp3->values_cpx, &(Tmp3->ld), &diff, &scale1, &workspace_cpx, &lwork, iwork, &info);
        lwork = nearbyint(creal(workspace_cpx));
        mess_try_alloc(work_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*lwork);

        F77_GLOBAL(ztgsyl,ZTGSYL)("N", &ijob,   &(data->rowsX), &(data->colsX), data->Ahat->values_cpx, &(data->Ahat->ld), data->Hhat->values_cpx, &(data->Hhat->ld), Mtilde->values_cpx, &(Mtilde->ld),
                data->Ehat->values_cpx, &(data->Ehat->ld), data->Fhat->values_cpx, &(data->Fhat->ld), Tmp3->values_cpx, &(Tmp3->ld), &diff, &scale1, work_cpx, &lwork, iwork, &info);

        ret = mess_matrix_scale(-1.0/scale1, Tmp3);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
    }

    mess_matrix_scale(-1.0, data->Hhat);

    /*-----------------------------------------------------------------------------
     *  check info value
     *-----------------------------------------------------------------------------*/
    if ( info < 0 ) {
        MSG_ERROR("DTRSYL / ZTRSYL / DTGSYL / ZTGSYL returned with error " MESS_PRINTF_INT " \n", info);
        ret = MESS_ERROR_LAPACK;
        goto clean;
    }

    if ( info == 1 ) {
        MSG_WARN("Pertubed values were used to solve the Sylvester equation");
    }

    if ( (scale1 != 1) || (scale2 != 1)){
        MSG_WARN("DTRSYL / ZTRSYL / DTGSYL / ZTGSYL returned scale1 = %e and scale2 = %e to avoid overflow.\n", scale1, scale2);
    }

    /*-----------------------------------------------------------------------------
     *  post processing
     *-----------------------------------------------------------------------------*/
    // Backslash variant
    //ret = mess_matrix_backslashm(MESS_OP_NONE, data->Ehat, Tmp3, Mtilde);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_backslashm);

    // Hessenberg Solver variant
    mess_direct hessenberg;
    mess_direct_init(&hessenberg);
    mess_direct_create_hessenberg_lu(data->Ehat, hessenberg);
    mess_direct_solvem(MESS_OP_NONE, hessenberg, Tmp3, Mtilde);
    mess_direct_clear(&hessenberg);

    ret = mess_matrix_multiply(MESS_OP_NONE, data->QE, MESS_OP_NONE, Mtilde, Tmp3);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, Tmp3, MESS_OP_HERMITIAN, data->QH, X);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
clean:
    MESS_CLEAR_MATRICES(&Mtilde,&Xtilde,&Tmp1,&Tmp2,&Tmp3);
    if(work) mess_free(work);
    if(work_cpx) mess_free(work_cpx);
    if(iwork) mess_free(iwork);
    return ret;


}

/**
 * @internal
 * @brief Solve \f$ A^HXF^H+E^HXH^H+M=0 \f$.
 * @param[in] datain  input solver data
 * @param[in] M  input right hand sides
 * @param[out] X  output solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref sylvester_d_solvemh_generalized function solves dense Sylvester equations:
 *  <ul>
 *  <li> \f$A^HXF^H+E^HXH^H+M=0\f$
 *  </ul>
 *
 * @attention Internal use only.
 */
static int sylvester_d_solvemh_generalized(void * datain, mess_matrix M, mess_matrix X){
    MSG_FNAME(__func__);
    _sylv_solver_d * data = (_sylv_solver_d *)datain;
    int ret = 0;
    mess_matrix Mtilde, Xtilde, Tmp1, Tmp2, Tmp3;
    double scale1 = 1, scale2 = 1, diff;
    mess_int_t info = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(X);
    mess_check_nullpointer(M);
    mess_check_dense(M);
    mess_check_real_or_complex(M);

    if ( M->rows != data->rowsX) {
        MSG_ERROR("number of rows doesn't match.\n");
        return MESS_ERROR_DIMENSION;
    }

    if ( M->cols != data->colsX){
        MSG_ERROR("number of columns doesn't match.\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  init and alloc some additional matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&Mtilde,&Xtilde,&Tmp1,&Tmp2,&Tmp3);

    /*-----------------------------------------------------------------------------
     *  Idea:
     *  Solve A^H X F^H + E^H X H^H + M = 0.
     *  (Schur decomposition) <==>
     *  Ahat^H * (QA^H * X * QF) * Fhat^H + Ehat^H * (QA^H * X * QF) * Hhat^H + QE^H M QH = 0
     *  <==>
     *  Ahat^H * Y * Fhat^H + Ehat^H* Y * Hhat^H + Mtilde = 0
     *  <==>
     *  Ahat^H R - Hhat^H L = Mtilde * Hhat^{-1}
     *  Ehat^H R + Fhat^H L = 0.
     *  Y = - L.
     *-----------------------------------------------------------------------------*/


    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, data->QE, MESS_OP_NONE, M, X);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, X, MESS_OP_NONE, data->QH, Xtilde);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

    // Backslash Variant
    //ret = mess_matrix_ctranspose(Xtilde, Mtilde);                                                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
    //ret = mess_matrix_backslashm(MESS_OP_NONE, data->Hhat, Mtilde, Xtilde);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_backslashm);
    //ret = mess_matrix_ctranspose(Xtilde, Mtilde);                                                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);

    // Hessenberg Solver variant
    mess_direct hessenberg;
    ret = mess_direct_init(&hessenberg);                                                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
    ret = mess_direct_create_hessenberg_lu(data->Hhat, hessenberg);                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_hessenberg_lu);
    ret = mess_matrix_ctranspose(Xtilde, Mtilde);                                                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
    ret = mess_direct_solvem(MESS_OP_NONE, hessenberg, Mtilde, Xtilde);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
    ret = mess_matrix_ctranspose(Xtilde, Mtilde);                                                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
    mess_direct_clear(&hessenberg);



    ret = mess_matrix_alloc(Tmp3, M->rows, M->cols, M->rows*M->cols, MESS_DENSE, data->isreal?MESS_REAL:MESS_COMPLEX);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);


    /*-----------------------------------------------------------------------------
     *  generalized dense Sylvester equation A X F + E X H = -M
     *  we first do a workspace query and then the call
     *-----------------------------------------------------------------------------*/
    double workspace, *work=NULL;
    mess_double_cpx_t workspace_cpx, *work_cpx=NULL;
    mess_int_t lwork, *iwork=NULL, ijob=0;
    lwork = -1;
    mess_try_alloc(iwork, mess_int_t *, sizeof(mess_int_t)*(data->rowsX + data->colsX + 6));

    mess_matrix_scale(-1.0, data->Hhat);
    if(data->isreal){
        if(MESS_IS_REAL(Mtilde)){
            F77_GLOBAL(dtgsyl,DTGSYL)("T", &ijob, &(data->rowsX), &(data->colsX),   data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Mtilde->values, &(Mtilde->ld),
                    data->Ehat->values, &(data->Ehat->ld), data->Fhat->values, &(data->Fhat->ld), Tmp3->values, &(Tmp3->ld), &diff, &scale1, &workspace, &lwork, iwork, &info);
            lwork = nearbyint(workspace);
            mess_try_alloc(work, double*, sizeof(double)*lwork);
            F77_GLOBAL(dtgsyl,DTGSYL)("T", &ijob, &(data->rowsX), &(data->colsX),   data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Mtilde->values, &(Mtilde->ld),
                    data->Ehat->values, &(data->Ehat->ld), data->Fhat->values, &(data->Fhat->ld), Tmp3->values, &(Tmp3->ld), &diff, &scale1, work, &lwork, iwork, &info);

            ret = mess_matrix_scale(-1.0/scale1, Tmp3);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
        }else{
            ret = mess_matrix_realpart(Mtilde,Tmp1);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realpart);
            F77_GLOBAL(dtgsyl,DTGSYL)("T", &ijob,   &(data->rowsX), &(data->colsX), data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Tmp1->values, &(Tmp1->ld),
                    data->Ehat->values, &(data->Ehat->ld), data->Fhat->values, &(data->Fhat->ld), Tmp3->values, &(Tmp3->ld), &diff, &scale1, &workspace, &lwork, iwork, &info);
            lwork = nearbyint(workspace);
            mess_try_alloc(work, double*, sizeof(double)*lwork);

            F77_GLOBAL(dtgsyl,DTGSYL)("T", &ijob,   &(data->rowsX), &(data->colsX), data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Tmp1->values, &(Tmp1->ld),
                    data->Ehat->values, &(data->Ehat->ld), data->Fhat->values, &(data->Fhat->ld), Tmp3->values, &(Tmp3->ld), &diff, &scale1, work, &lwork, iwork, &info);

            ret = mess_matrix_copy(Tmp3, Tmp1);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            ret = mess_matrix_imagpart(Mtilde,Tmp2);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_imagpart);
            ret = mess_matrix_zeros(Tmp3);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);

            F77_GLOBAL(dtgsyl,DTGSYL)("T", &ijob,   &(data->rowsX), &(data->colsX), data->Ahat->values, &(data->Ahat->ld), data->Hhat->values, &(data->Hhat->ld), Tmp2->values, &(Tmp2->ld),
                    data->Ehat->values, &(data->Ehat->ld), data->Fhat->values, &(data->Fhat->ld), Tmp3->values, &(Tmp3->ld), &diff, &scale2, work, &lwork, iwork, &info);


            ret = mess_matrix_addc(-1.0/scale1, Tmp1, -1.0*I/scale2, Tmp3);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        }
    }else{
        ret = mess_matrix_tocomplex(Mtilde);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        F77_GLOBAL(ztgsyl,ZTGSYL)("C", &ijob,   &(data->rowsX), &(data->colsX), data->Ahat->values_cpx, &(data->Ahat->ld), data->Hhat->values_cpx, &(data->Hhat->ld), Mtilde->values_cpx, &(Mtilde->ld),
                data->Ehat->values_cpx, &(data->Ehat->ld), data->Fhat->values_cpx, &(data->Fhat->ld), Tmp3->values_cpx, &(Tmp3->ld), &diff, &scale1, &workspace_cpx, &lwork, iwork, &info);
        lwork = nearbyint(creal(workspace_cpx));
        mess_try_alloc(work_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*lwork);

        F77_GLOBAL(ztgsyl,ZTGSYL)("C", &ijob,   &(data->rowsX), &(data->colsX), data->Ahat->values_cpx, &(data->Ahat->ld), data->Hhat->values_cpx, &(data->Hhat->ld), Mtilde->values_cpx, &(Mtilde->ld),
                data->Ehat->values_cpx, &(data->Ehat->ld), data->Fhat->values_cpx, &(data->Fhat->ld), Tmp3->values_cpx, &(Tmp3->ld), &diff, &scale1, work_cpx, &lwork, iwork, &info);

        ret = mess_matrix_scale(-1.0/scale1, Tmp3);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
    }

    mess_matrix_scale(-1.0, data->Hhat);

    /*-----------------------------------------------------------------------------
     *  check info value
     *-----------------------------------------------------------------------------*/
    if ( info < 0 ) {
        MSG_ERROR("DTRSYL / ZTRSYL / DTGSYL / ZTGSYL returned with error " MESS_PRINTF_INT " \n", info);
        ret = MESS_ERROR_LAPACK;
        goto clean;
    }

    if ( info == 1 ) {
        MSG_WARN("Pertubed values were used to solve the Sylvester equation");
    }

    if ( (scale1 != 1) || (scale2 != 1)){
        MSG_WARN("DTRSYL / ZTRSYL / DTGSYL / ZTGSYL returned scale1 = %e and scale2 = %e to avoid overflow.\n", scale1, scale2);
    }

    /*-----------------------------------------------------------------------------
     *  post processing
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_multiply(MESS_OP_NONE, data->QA, MESS_OP_NONE, Tmp3, Mtilde);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, Mtilde, MESS_OP_HERMITIAN, data->QF, X);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);


    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
clean:
    MESS_CLEAR_MATRICES(&Mtilde,&Xtilde,&Tmp1,&Tmp2,&Tmp3);
    mess_free(iwork);
    if(work)mess_free(work);
    if(work_cpx)mess_free(work_cpx);
    return ret;


}


/**
 * @internal
 * @brief Solve \f$ A^TXF^T+E^TXH^T+M=0 \f$.
 * @param[in] datain  input solver data
 * @param[in] M  input right hand sides
 * @param[out] X  output solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref sylvester_d_solvemt_generalized function solves dense Sylvester equations:
 *  <ul>
 *  <li> \f$A^TXF^T+E^TXH^T+M=0\f$
 *  </ul>
 *
 * @attention Internal use only.
 */
static int sylvester_d_solvemt_generalized(void *datain, mess_matrix M, mess_matrix X){
    MSG_FNAME(__func__);
    int ret = 0;
    ret = mess_matrix_conj(M);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_conj);
    ret = sylvester_d_solvemh_generalized(datain, M, X);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), sylvester_d_solvemh_generalized);
    ret = mess_matrix_conj(X);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_conj);
    ret = mess_matrix_conj(M);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_conj);
    return ret;
}




/**
 * @brief Generate a solver for a dense Sylvester equation using @lapack.
 * @param[in] A  input matrix
 * @param[in] F  input matrix
 * @param[in] E  input matrix
 * @param[in] H  input matrix
 * @param[out] solver solver to create
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_direct_create_sylvester_dense function generates a Sylvester
 * equation solver using @lapack.
 * It solves
 * <ul>
 * <li> Standard:  \f$AX+XH+M=0\f$, \f$A^TX+XH^T+M=0\f$, \f$A^HX+XH^H+M=0\f$
 * <li> Generalized Sylvester: \f$AXF+EXH+M=0\f$, \f$A^TXF^T+E^TXH^T+M=0\f$, \f$A^HXF^H+E^HXH^H+M=0\f$
 * </ul>
 * with \f$ A,F,E,H \f$ small and dense matrices.
 * If a standard Sylvester equation solver is wanted use @p NULL and @p NULL for the matrices @p E and @p F.
 * @p F must be a regular matrix.
 *
 */
int mess_direct_create_sylvester_dense(mess_matrix A, mess_matrix F, mess_matrix E, mess_matrix H, mess_direct solver){
    MSG_FNAME(__func__);
    int ret =0;
    _sylv_solver_d* data;
    int generalized = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(solver);
    mess_check_nullpointer(A);
    mess_check_dense(A);
    mess_check_square(A);
    mess_check_real_or_complex(A);
    mess_check_nullpointer(H);
    mess_check_dense(H);
    mess_check_square(H);
    mess_check_real_or_complex(H);

    if((F && !E) || (!E && F)){
        MSG_ERROR("Both matrices E and F must be given.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if(E){
        mess_check_nullpointer(E);
        mess_check_dense(E);
        mess_check_square(E);
        mess_check_real_or_complex(E);
        mess_check_same_size(A,E);
        generalized = 1;
    }

    if(F){
        mess_check_nullpointer(F);
        mess_check_dense(F);
        mess_check_square(F);
        mess_check_real_or_complex(F);
        mess_check_same_size(F,H);
    }

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data, _sylv_solver_d*, sizeof(_sylv_solver_d));
    data->isreal = (MESS_IS_REAL(A) && MESS_IS_REAL(H) && (!(E&&F) || (MESS_IS_REAL(E) && MESS_IS_REAL(F)))) ? 1:0;
    MESS_INIT_MATRICES(&(data->Ahat), &(data->QA), &(data->Ehat), &(data->QE), &(data->Fhat), &(data->QF), &(data->Hhat), &(data->QH));
    data->rowsX=A->cols;
    data->colsX=H->rows;

    /*-----------------------------------------------------------------------------
     *  Standard Schur Decompositions
     *  A = QA * Ahat * QA^H
     *  H = QH * Hhat * QH^H
     *
     *  Generalized Schur Decompositions
     *  A = QA * Ahat * QE^H and E = QA * Ehat * QE^H
     *  H = QH * Hhat * QF^H and F = QH * Fhat * QF^H
     *-----------------------------------------------------------------------------*/
    if(!generalized){
        if(data->isreal){
            ret = mess_eigen_schur(A, data->Ahat, data->QA, NULL);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_schur);
            ret = mess_eigen_schur(H, data->Hhat, data->QH, NULL);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_schur);
        }else{
            ret = mess_eigen_schur_complex(A, data->Ahat, data->QA, NULL);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_schur_complex);
            ret = mess_eigen_schur_complex(H, data->Hhat, data->QH, NULL);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_schur_complex);
        }
    }else{
        if(data->isreal){
            ret = mess_eigen_gschur(A, E, data->Ahat, data->Ehat, data->QA, data->QE, NULL, NULL, NULL);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_gschur);
            ret = mess_eigen_gschur(H, F, data->Hhat, data->Fhat, data->QH, data->QF, NULL, NULL, NULL);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_gschur);
        }else{
            ret = mess_eigen_gschur_complex(A, E, data->Ahat, data->Ehat, data->QA, data->QE, NULL, NULL, NULL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_gschur_complex);
            ret = mess_eigen_gschur_complex(H, F, data->Hhat, data->Fhat, data->QH, data->QF, NULL, NULL, NULL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_gschur_complex);
        }
    }

    /*-----------------------------------------------------------------------------
     *  finalize solver
     *-----------------------------------------------------------------------------*/
    solver->data_type = data->isreal?MESS_REAL:MESS_COMPLEX;
    solver->rows = data->rowsX;
    solver->cols = data->colsX;
    if(!generalized){
        SET_SOLVERNAME(solver->name, "SYLV_DENSE_STANDARD");
    }else{
        SET_SOLVERNAME(solver->name, "SYLV_DENSE_GENERALIZED");
    }

    solver->data = (void *)data;
    solver->solve = NULL;
    solver->solvet = NULL;
    solver->solveh = NULL;

    /*-----------------------------------------------------------------------------
     *  set function pointers for standard and generalized solver
     *-----------------------------------------------------------------------------*/
    if(!generalized){
        solver->solvem = sylvester_d_solvem_standard;
        solver->solvemt = sylvester_d_solvemt_standard;
        solver->solvemh = sylvester_d_solvemh_standard;

    }else{
        solver->solvem = sylvester_d_solvem_generalized;
        solver->solvemt = sylvester_d_solvemt_generalized;
        solver->solvemh = sylvester_d_solvemh_generalized;

    }
    solver->clear = sylvester_d_clear;
    solver->getL = NULL;
    solver->getU = NULL;
    solver->getpermp = NULL;
    solver->getpermq = NULL;

    return 0;
}

