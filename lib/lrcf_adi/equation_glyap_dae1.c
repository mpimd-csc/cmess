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
 * @file lib/lrcf_adi/equation_glyap_dae1.c
 * @brief Function handles for Index 1 DAE System.
 * @author @mbehr
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @internal
 * @brief Contains matrices, solvers and additional data for a DAE 1 function handles.
 *
 * The @ref __glyap_dae1 structure contains all matrices, solvers and additional data for DAE 1 function handles.
 * @attention Internal use only.
 *
 */
typedef struct {
    //system matrices
    mess_matrix E11;                /**< System matrix \f$E_{11}\f$.*/
    mess_direct E11sol;         /**< Direct solver for \f$E_{11}\f$.*/
    mess_matrix fullE;          /**< System matrix \f$\begin{bmatrix}E_{11} & 0 \\ 0 & 0 \end{bmatrix}\f$.*/

    mess_matrix A11;            /**< System matrix \f$A_{11}\f$.*/
    mess_matrix A12;            /**< System matrix \f$A_{12}\f$.*/
    mess_matrix A21;            /**< System matrix \f$A_{21}\f$.*/
    mess_matrix A22;            /**< System matrix \f$A_{22}\f$.*/

    mess_direct A22sol;         /**< Direct solver for \f$A_{22}\f$.*/

    mess_matrix fullA;          /**< System block matrix \f$\begin{bmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix}\f$.*/
    mess_direct fullAsol;       /**< Direct solver fo system block matrix \f$\begin{bmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix}\f$.*/

    mess_int_t n1;          /**< Number of rows of \f$A_{11}\f$.*/
    mess_int_t n2;          /**< Number of rows of \f$A_{21}\f$.*/

    mess_multidirect Amsolver;      /**< Multidirect solver for the shifted system \f$\begin{bmatrix} A_{11}+pE_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix} \f$.*/

} __glyap_dae1;


/**
 * @internal
 * @brief Compute matrix-matrix multiplication with \f$ A_{11} - A_{12}A_{22}^{-1}A_{22}\f$.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref AX_apply function computes the matrix-matrix product, i.e. \n
 * \f[ in = op(A_{11} - A_{12}A_{22}^{-1}A_{22}) out. \f]
 *
 * @see mess_matrix_multiply
 * @see mess_matrix_add
 * @see mess_direct_solvem
 *
 * @attention Internal use only.
 */
static int AX_apply(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae1 *eqn = (__glyap_dae1 *) e->aux;
    mess_int_t ret=0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);
    mess_check_operation_type(op);

    if(in->rows == eqn->n1){
        /*-----------------------------------------------------------------------------
         *  this part is only coded for residual computation or shift parameter computation
         *
         *-----------------------------------------------------------------------------*/
        mess_matrix tmp1,tmp2; MESS_INIT_MATRICES(&tmp1,&tmp2);
        //op(A*PI^T)*x->out
        if(op==MESS_OP_NONE){
            ret = mess_matrix_multiply(MESS_OP_NONE,eqn->A21,MESS_OP_NONE,in,tmp1);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
            ret = mess_direct_solvem(MESS_OP_NONE,eqn->A22sol,tmp1,tmp2);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
            ret = mess_matrix_multiply(MESS_OP_NONE,eqn->A12,MESS_OP_NONE,tmp2,tmp1);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_NONE,eqn->A11,MESS_OP_NONE,in,out);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
            ret = mess_matrix_add(-1.0,tmp1,1.0,out);                                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
        }else{
            ret = mess_matrix_multiply(MESS_OP_TRANSPOSE,eqn->A12,MESS_OP_NONE,in,tmp1);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
            ret = mess_direct_solvem(MESS_OP_TRANSPOSE,eqn->A22sol,tmp1,tmp2);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
            ret = mess_matrix_multiply(MESS_OP_TRANSPOSE,eqn->A21,MESS_OP_NONE,tmp2,tmp1);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_TRANSPOSE,eqn->A11,MESS_OP_NONE,in,out);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
            ret = mess_matrix_add(-1.0,tmp1,1.0,out);                                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
        }
        MESS_CLEAR_MATRICES(&tmp1,&tmp2);
    }else{
        return MESS_ERROR_NOSUPPORT;
    }

    return 0;
}

/**
 * @internal
 * @brief Solve linear system with block matrix \f$ op(A_{11} - A_{12}A_{22}^{-1}A_{22})\f$.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref AINV_apply function solve a linear system \n
 * \f[
 *  op(A_{11} - A_{12}A_{22}^{-1}A_{22}) out = in
 * \f]
 * by using the equivalent block formulation
 * \f[
 *  op\left(\begin{bmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix}\right)
 *  \begin{bmatrix} out \\ \ast \end{bmatrix}
 *  =
 *  \begin{bmatrix} in \\ 0 \end{bmatrix}.
 *  \f]
 *
 *
 * @see mess_matrix_rowsub
 * @see mess_direct_solvem
 * @see mess_matrix_lift
 *
 * @attention Internal use only.
 */
static int AINV_apply(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae1 *eqn = (__glyap_dae1 *) e->aux;
    mess_int_t ret=0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);
    mess_check_operation_type(op);

    /*-----------------------------------------------------------------------------
     *  perfom operation
     *-----------------------------------------------------------------------------*/
    if(in->rows == eqn->n1){
        /*-----------------------------------------------------------------------------
         *  this part is only coded for shift parameter computation
         *-----------------------------------------------------------------------------*/
        mess_matrix tmp1,tmp2;
        MESS_INIT_MATRICES(&tmp1,&tmp2);
        ret = mess_matrix_lift(in,eqn->n2,tmp1);                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_lift);
        ret = mess_direct_solvem(op,eqn->fullAsol,tmp1,tmp2);                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
        ret = mess_matrix_rowsub(tmp2,0,eqn->n1-1,out);                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
        MESS_CLEAR_MATRICES(&tmp1,&tmp2);
    }else{
        return MESS_ERROR_NOSUPPORT;
    }

    return 0;
}

/**
 * @internal
 * @brief Compute matrix-matrix multiplication with \f$ E_{11}\f$.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref EX_apply function computes the matrix-matrix product, i.e. \n
 * \f[ out = op(E_{11})in. \f]
 *
 * @see mess_matrix_multiply
 *
 * @attention Internal use only.
 */
static int EX_apply(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae1 *eqn = (__glyap_dae1 *) e->aux;
    int ret=0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);
    mess_check_operation_type(op);
    mess_check_same_colsrows(eqn->E11,in);

    /*-----------------------------------------------------------------------------
     *  perform E*in -> out
     *-----------------------------------------------------------------------------*/
    if(in->rows == eqn->n1){
        ret = mess_matrix_multiply(op,eqn->E11,MESS_OP_NONE,in,out);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
    }else{
        return MESS_ERROR_NOSUPPORT;
    }

    return 0;
}

/**
 * @internal
 * @brief Solve linear system with block matrix \f$ op(E_{11})\f$.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref EINV_apply function solve a linear system \n
 * \f[
 *  op(E_{11} - A_{12}A_{22}^{-1}A_{22}) out = in.
 * \f]
 *
 * @see mess_direct_solvem
 *
 * @attention Internal use only.
 */
static int EINV_apply(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae1 *eqn = (__glyap_dae1 *) e->aux;
    mess_int_t ret=0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);
    mess_check_operation_type(op);

    /*-----------------------------------------------------------------------------
     *  perfom operation
     *-----------------------------------------------------------------------------*/
    if(in->rows == eqn->n1){
        ret = mess_direct_solvem(op,eqn->E11sol,in,out);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
    }else{
        return MESS_ERROR_NOSUPPORT;
    }

    return 0;
}

/**
 * @internal
 * @brief Generate a @ref mess_multidirect solve for the system matrix \f$ \begin{bmatrix} A_{11} + pE_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix} \f$.
 * @param[in] e             input @ref mess_equation structure
 * @param[in] parameters    input @ref mess_vector with shift parameters
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref ApEINV2_generate builds a @ref mess_multidirect solver for the matrices \n
 * \f[
 *  \begin{bmatrix} A_{11} + pE_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix}.
 * \f]
 *
 * @attention Internal use only.
 *
 */
static int ApEINV2_generate(mess_equation e, mess_vector parameters) {
    MSG_FNAME(__func__);
    mess_int_t ret = 0;
    mess_check_nullpointer(e);
    __glyap_dae1 *eqn = (__glyap_dae1 *) e->aux;

    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(parameters);

    /*-----------------------------------------------------------------------------
     *  perform decomposition for every shift
     *-----------------------------------------------------------------------------*/
    if ( e->ApEINV.to_clear) return 0;

    ret = mess_multidirect_init(&(eqn->Amsolver));                                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_multidirect_init);
    ret = mess_multidirect_create(eqn->fullA, NULL, parameters , eqn->Amsolver, NULL, eqn->fullE);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_multidirect_create);
    e->ApEINV.to_clear = 1;
    return 0;
}

/**
 * @internal
 * @brief Solve shifted linear block system with matrix \f$ \begin{bmatrix} A_{11} + pE_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix} \f$.
 * @param[in, out] e        input @ref mess_equation structure
 * @param[in]   op          input operation type
 * @param[in]   p           input shift parameter
 * @param[in]   idx_p       input index of parameter for shift
 * @param[in]   in          input matrix
 * @param[out]  out         output matrix
 * @return zero on succes or a non-zero error value otherwise
 *
 *  The @ref ApEINV2_apply function solves the system
 * \f[
 *  op((A_{11} - A_{12}A_{22}^{-1}A_{21}) + p E_{11}) out = in
 * \f]
 * via the equivalent block formulation
 *
 *
 * \f[
 *  op\left(\begin{bmatrix} A_{11} + pE_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix}\right)
 *  \begin{bmatrix} out \\ \ast \end{bmatrix}
 *  =
 *  \begin{bmatrix} in \\ 0 \end{bmatrix}.
 * \f]
 *
 *
 * @attention Internal use only.
 */
static int ApEINV2_apply(mess_equation e, mess_operation_t op, mess_double_cpx_t p, mess_int_t idx_p, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    __glyap_dae1 *eqn = (__glyap_dae1 *) e->aux;
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);

    mess_matrix temp1, temp2;
    MESS_INIT_MATRICES(&temp1,&temp2);
    ret = mess_matrix_lift(in,eqn->n2,temp1);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_lift);
    ret = mess_multidirect_solvem(op,eqn->Amsolver, idx_p,temp1,temp2);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_multidirect_solvem);
    ret = mess_matrix_rowsub(temp2,0,eqn->n1-1,out);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_rowsub);
    MESS_CLEAR_MATRICES(&temp1,&temp2);



    return 0;
}

/**
 * @internal
 * @brief Clear a @ref mess_multidirect solve for the system matrix \f$ \begin{bmatrix} A_{11} + pE_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix} \f$.
 * @param[in] e         input @ref mess_equation structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref ApEINV2_clear builds a @ref mess_multidirect solver for the matrices \n
 * \f[
 *  \begin{bmatrix} A_{11} + pE_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix}.
 * \f]
 * by using @ref mess_multidirect_clear the field @ref __glyap_dae1.Amsolver is cleared.
 *
 * @attention Internal use only.
 *
 */
static int ApEINV2_clear(mess_equation e){
    __glyap_dae1 *eqn = (__glyap_dae1 *) e->aux;
    MSG_FNAME(__func__);
    mess_check_nullpointer(eqn);
    if ( e->ApEINV.to_clear == 0) return 0;
    if ( eqn -> Amsolver != NULL ) {
        mess_multidirect_clear(&(eqn->Amsolver));
    }
    e->ApEINV.to_clear = 0;
    return 0;
}

/**
 * @internal
 * @brief Call the @ref mess_lrcfadi_parameter function.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] opt   input options
 * @param[in] stat  input status
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref parameter function calls directly @ref mess_lrcfadi_parameter.
 * If @ref mess_equation.parent of @p in is not @c NULL the equation was stabilized
 * by a feedback and @ref mess_lrcfadi_parameter is called
 * with the stabilized equation object @ref mess_equation.parent otherwise
 * @p e is simply used.
 *
 * @see mess_lrcfadi_parameter
 *
 * @attention Internal use only.
 * */
static int parameter(mess_equation e, mess_options opt, mess_status stat){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae1 *eqn = (__glyap_dae1 *) e->aux;
    mess_check_nullpointer(eqn);

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(opt);
    mess_check_nullpointer(e);
    mess_check_nullpointer(stat);

    /*-----------------------------------------------------------------------------
     *  adaptive shift parameters need no lift for matrix B and K
     *-----------------------------------------------------------------------------*/

    if(e->parent){
        //equation was stabilised by equation_stable
        //in this case we use the stabilized function handles to compute shift parameters
        return mess_lrcfadi_parameter(e->parent,opt,stat);
    }
    //equation was not stabilised by equation_stable
    return mess_lrcfadi_parameter(e,opt,stat);
}

/**
 * @internal
 * @brief Clear additional allocated field of @ref __glyap_dae1 structure.
 * @param[in, out] e         input @ref mess_equation structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref clear function, clears the following additional allocated data
 *
 * <center>
 * |          field                           |  mathematical description                                                                       |
 * |:----------------------------------------:|:-----------------------------------------------------------------------------------------------:|
 * | @ref __glyap_dae1.fullA                  | \f$ \begin{bmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix} \f$                        |
 * | @ref __glyap_dae1.fullE                  | \f$ \begin{bmatrix} E_{11} & 0 \\ 0 & 0 \end{bmatrix} \f$                       |
 * | @ref __glyap_dae1.A22sol                 | \f$ A_{22}^{-1} \f$                                             |
 * | @ref __glyap_dae1.fullAsol               | \f$ \left( \begin{bmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix}\right)^{-1} \f$     |
 * | @ref __glyap_dae1.E11sol                 | \f$ E_{11}^{-1} \f$                                         |
 * </center>.
 *
 * @attention Internal use only.
 * */
static int clear(mess_equation e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae1 *eqn = (__glyap_dae1 *) e->aux;
    mess_check_nullpointer(eqn);

    //clear additional matrices
    if(eqn->fullA) mess_matrix_clear(&(eqn->fullA));
    if(eqn->fullE) mess_matrix_clear(&(eqn->fullE));


    //clear additional direct solvers
    if(eqn->A22sol)     mess_direct_clear(&(eqn->A22sol));
    if(eqn->fullAsol)   mess_direct_clear(&(eqn->fullAsol));
    if(eqn->E11sol)     mess_direct_clear(&(eqn->E11sol));

    if(eqn)     mess_free(eqn);
    eqn = NULL;
    e->aux = NULL;
    return 0;
}


/**
 * @brief Generate an Equation object to solve a generalized Lyapunov Equation from Index 1 DAE.
 * @param[in,out]   e           equation object
 * @param[in]       opt         input options for ADI iteration
 * @param[in]       E11         input \f$ E_{11} \in \mathbb{R}^{n_1 \times n_1} \f$ matrix of the DAE system
 * @param[in]       A11         input \f$ A_{11} \in \mathbb{R}^{n_1 \times n_1} \f$ matrix of the DAE system
 * @param[in]       A12         input \f$ A_{12} \in \mathbb{R}^{n_1 \times n_2} \f$ matrix of the DAE system
 * @param[in]       A21         input \f$ A_{21} \in \mathbb{R}^{n_2 \times n_1} \f$ matrix of the DAE system
 * @param[in]       A22         input \f$ A_{22} \in \mathbb{R}^{n_2 \times n_2} \f$ matrix of the DAE system
 * @param[in]       rhs             input \f$ rhs \in \mathbb{R}^{n_1+n_2 \times p }\f$ or \f$ rhs \in \mathbb{R}^{p \times n_1 + n_2} \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_equation_glyap_dae1 function creates a mess_equation object from Index 1 DAE System
 *
 *  \f[
 *  \begin{aligned}
 *  \begin{bmatrix} E_{11} & 0 \\ 0 & 0  \end{bmatrix}
 *  \begin{bmatrix} \dot{x_1} \\ \dot{x_2}   \end{bmatrix}
 *  =
 *  \begin{bmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix}
 *  \begin{bmatrix} x_1 \\ x_2 \end{bmatrix}
 *  +
 *  \begin{bmatrix} B \\ 0 \end{bmatrix}
 *  u  \\
 *  \end{aligned}
 *  \f]
 * Preliminaries:
 *  @li \f$ E_{11}, A_{11} \in \mathbb{R}^{n_1 \times n_1}, A_{12} \in \mathbb{R}^{n_1 \times n_2}, A_{21} \in \mathbb{R}^{n_2 \times n_1}, A_{22} \in \mathbb{R}^{n_2 \times n_2} \f$.
 *  @li \f$ A_{22} \f$ is invertible.
 *
 *  Using that \f$ A_{22} \f$ is invertible leads to the ODE system
 *
 *  \f[
 *  \begin{aligned}
 *  E_{11} \frac{dx_1}{dt} =
 *  (A_{11} - A_{12}A_{22}^{-1}A_{21})x_1 + (B(1:n_1,:)- A_{12}A_{22}^{-1}B(n_1+1:end,:)) u
 *  \end{aligned}
 *  \f]
 *
 *
 * Depending on the @ref mess_operation_t defined in @ref mess_options.type we build an equation
 * structure to solve the following Lyapunov Equations:
 * If the operation type is @ref MESS_OP_NONE, we solve:
 * \f[
 * \begin{aligned}
 *  E_{11} X (A_{11}-A_{12}A_{22}^{-1}A_{21})^T + (A_{11}-A_{12}A_{22}^{-1}A_{21}) X E_{11}^T = -BB^T \\
 *  B := (rhs(1:n_1,:) - A_{12}A_{22}^{-1}rhs(n_1+1:end,:))
 * \end{aligned}
 * \f]
 *
 * otherwise
 * \f[
 * \begin{aligned}
 *  E_{11}^T X (A_{11}-A_{12}A_{22}^{-1}A_{21}) + (A_{11}-A_{12}A_{22}^{-1}A_{21})^T X E_{11} = -B^TB \\
 *  B := (rhs(:,1:n_1) - rhs(:,n_1+1:end) A_{22}^{-1}A_{21}).
 * \end{aligned}
 * \f]
 *
 *
 * \sa mess_lrcfadi_adi
 * \sa mess_equation_lyap
 * \sa mess_equation_st
 * \sa mess_equation_griccati_dae1
 *
 * References: \cite morGugSW13
 *
 */
int mess_equation_glyap_dae1(mess_equation e, mess_options opt , mess_matrix E11, mess_matrix A11, mess_matrix A12, mess_matrix A21, mess_matrix A22,  mess_matrix rhs) {
    MSG_FNAME(__func__);
    __glyap_dae1 *data ;
    int ret = 0;
    mess_operation_t optype;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(e);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(E11);    mess_check_real(E11); mess_check_square(E11);
    mess_check_nullpointer(A11);    mess_check_real(A11); mess_check_square(A11); mess_check_same_size(A11,E11);
    mess_check_nullpointer(A12);    mess_check_real(A12);
    mess_check_nullpointer(A21);    mess_check_real(A21); mess_check_same_colsrows(A12,A21); mess_check_same_colsrows(A21,A12);
    mess_check_nullpointer(A22);    mess_check_real(A22); mess_check_square(A22); mess_check_same_cols(A22,A12); mess_check_same_rows(A22,A21);
    mess_check_nullpointer(rhs);    mess_check_dense(rhs);

    optype = (opt->type==MESS_OP_NONE)?MESS_OP_NONE:MESS_OP_TRANSPOSE;

    if(optype==MESS_OP_NONE){
        if(rhs->rows != (E11->rows+A22->rows)){
            MSG_ERROR("rhs has wrong number of rows.\n");
            return MESS_ERROR_ARGUMENTS;
        }
        if ( rhs->cols >= rhs->rows) {
            MSG_WARN("Oversized right hand side factor for LRCF-ADI (op = %s, rows=%d, cols=%d)\n ", mess_operation_t_str(optype), (int) rhs->rows, (int) rhs->cols);
        }
    }else{
        if(rhs->cols != (E11->rows+A22->rows)){
            MSG_ERROR("rhs has wrong number of rows.\n");
            return MESS_ERROR_ARGUMENTS;
        }
        if ( rhs->rows >= rhs->cols) {
            MSG_WARN("Oversized right hand side factor for LRCF-ADI (op = %s, rows=%d, cols=%d)\n ", mess_operation_t_str(optype), (int) rhs->rows, (int) rhs->cols);
        }
    }

    /*-----------------------------------------------------------------------------
     * Set DEFAULT
     *-----------------------------------------------------------------------------*/
    if ( e->clear != NULL ) {
        e->clear(e->aux);
    }

    if ( e->clearRHS) {
        ret = mess_matrix_clear(&e->RHS);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
        e->RHS=NULL;
        e->clearRHS = 0 ;
    }

    if ( e->clearB) {
        ret = mess_matrix_clear(&e->B);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
        e->B=NULL;
        e->clearB = 0;
    }

    if ( e->clearC) {
        ret = mess_matrix_clear(&e->C);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
        e->C=NULL;
        e->clearC = 0;
    }


    /*-----------------------------------------------------------------------------
     *  build the __glyap_dae1 structure
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data,__glyap_dae1 *, sizeof(__glyap_dae1));
    data->E11           = E11;
    data->A11           = A11;
    data->A12           = A12;
    data->A21           = A21;
    data->A22           = A22;
    data->n1            = E11->rows;
    data->n2            = A22->rows;


    /*-----------------------------------------------------------------------------
     *  add direct solvers to data and matrices to data
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_init(&data->A22sol);                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
    ret = mess_direct_lu(data->A22,data->A22sol);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);

    ret = mess_direct_init(&data->E11sol);                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
    ret = mess_direct_lu(data->E11,data->E11sol);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);

    MESS_INIT_MATRICES(&(data->fullA));
    ret = mess_matrix_cat(A11,A12,A21,A22,MESS_CSR,data->fullA);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);
    ret = mess_direct_init(&(data->fullAsol));                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
    ret = mess_direct_lu(data->fullA,data->fullAsol);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);

    mess_matrix zeros;
    MESS_INIT_MATRICES(&zeros,&data->fullE);
    ret = mess_matrix_alloc(zeros,data->n2,data->n2,0,MESS_CSR,MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
    ret = mess_matrix_cat(E11,NULL,NULL,zeros,MESS_CSR,data->fullE);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);
    MESS_CLEAR_MATRICES(&zeros);

    /*-----------------------------------------------------------------------------
     *  build equation structure
     *-----------------------------------------------------------------------------*/
    //set dimension
    e->dim                  = data->A11->rows;

    //set DAE 1 structure
    e->aux                  = data;

    //add function pointers
    e->eqn_type             = MESS_EQN_GLYAP;
    e->clear                = clear;
    e->AX.apply             = AX_apply;
    e->EX.apply             = EX_apply;
    e->AINV.generate        = NULL;
    e->AINV.apply           = AINV_apply;
    e->AINV.clear           = NULL;
    e->EINV.generate        = NULL;
    e->EINV.apply           = EINV_apply;
    e->EINV.clear           = NULL;
    e->ApEINV.generate      = ApEINV2_generate;
    e->ApEINV.apply         = ApEINV2_apply;
    e->ApEINV.clear         = ApEINV2_clear;
    e->parameter            = parameter;
    e->init_rhs             = NULL;

    /*-----------------------------------------------------------------------------
     * compute reduced right hand side
     *-----------------------------------------------------------------------------*/
    mess_matrix myrhs,temp1,temp2;
    MESS_INIT_MATRICES(&myrhs,&temp1,&temp2);
    if(optype == MESS_OP_NONE){
        ret = mess_matrix_rowsub(rhs,0,data->n1-1,myrhs);                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
        ret = mess_matrix_rowsub(rhs,data->n1,data->n1+data->n2-1,temp1);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
        ret = mess_direct_solvem(MESS_OP_NONE,data->A22sol,temp1,temp2);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
        ret = mess_matrix_multiply(MESS_OP_NONE,data->A12,MESS_OP_NONE,temp2,temp1);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        ret = mess_matrix_add(-1.0,temp1,1.0,myrhs);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);

        //set the reduced right hand side
        e->B            = myrhs;
        e->clearB       = 1;

        e->RHS          = myrhs;
        e->clearRHS     = 0;        //B is cleared by user


    }else{
        mess_matrix temp;
        MESS_INIT_MATRICES(&temp);
        mess_matrix_ctranspose(rhs,temp);                                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
        ret = mess_matrix_rowsub(temp,0,data->n1-1,myrhs);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
        ret = mess_matrix_rowsub(temp,data->n1,data->n1+data->n2-1,temp1);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
        ret = mess_direct_solvem(MESS_OP_TRANSPOSE,data->A22sol,temp1,temp2);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
        ret = mess_matrix_multiply(MESS_OP_TRANSPOSE,data->A21,MESS_OP_NONE,temp2,temp1);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        ret = mess_matrix_add(-1.0,temp1,1.0,myrhs);                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
        //set the reduced right hand side
        e->B            = myrhs;
        e->clearB       = 1;        //PIB is cleared by mess_equation_clear

        e->RHS          = myrhs;    //projection is done in  mess_lradi by init_rhs
        e->clearRHS     = 0;        //B is cleared by user

        MESS_CLEAR_MATRICES(&temp);

    }

    /*-----------------------------------------------------------------------------
     *  clear additional memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&temp1,&temp2);
    return 0;
}

/**
 * @brief Generate an Equation object to solve a generalized Lyapunov Equation from Index 1 DAE.
 * @param[in,out]   e           equation object
 * @param[in]       opt         input options for ADI iteration
 * @param[in]       E11         input \f$ E_{11} \in \mathbb{R}^{n_1 \times n_1} \f$ matrix of the DAE system
 * @param[in]       A11         input \f$ A_{11} \in \mathbb{R}^{n_1 \times n_1} \f$ matrix of the DAE system
 * @param[in]       A12         input \f$ A_{12} \in \mathbb{R}^{n_1 \times n_2} \f$ matrix of the DAE system
 * @param[in]       A21         input \f$ A_{21} \in \mathbb{R}^{n_2 \times n_1} \f$ matrix of the DAE system
 * @param[in]       A22         input \f$ A_{22} \in \mathbb{R}^{n_2 \times n_2} \f$ matrix of the DAE system
 * @param[in]       B           input \f$ B \in \mathbb{R}^{n_1+n_2 \times p }\f$ or \f$ B \in \mathbb{R}^{ n_1 \times p} \f$
 * @param[in]       C           input \f$ C \in \mathbb{R}^{q \times n_1 }\f$ or \f$ C \in \mathbb{R}^{ q \times n_1 + n_2} \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_equation_glyap_dae1 function creates a mess_equation object from Index 1 DAE System.
 *
 *  \f[
 *  \begin{aligned}
 *  \begin{bmatrix} E_{11} & 0 \\ 0 & 0  \end{bmatrix}
 *  \begin{bmatrix} \frac{dx_1}{dt} \\ \frac{dx_2}{dt}   \end{bmatrix}
 *  &=
 *  \begin{bmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{bmatrix}
 *  \begin{bmatrix} x_1 \\ x_2 \end{bmatrix}
 *  +
 *  B u  \\
 *  y &= C \begin{bmatrix}  x_1 & x_2 \end{bmatrix}
 *  \end{aligned}
 *  \f]
 * Preliminaries:
 *  @li \f$ E_{11}, A_{11} \in \mathbb{R}^{n_1 \times n_1}, A_{12} \in \mathbb{R}^{n_1 \times n_2}, A_{21} \in \mathbb{R}^{n_2 \times n_1}, A_{22} \in \mathbb{R}^{n_2 \times n_2} \f$.
 *  @li \f$ A_{22} \f$ is invertible.
 *
 *  Using that \f$ A_{22} \f$ is invertible leads to the ODE system
 *
 *  \f[
 *  \begin{aligned}
 *  E_{11} \frac{dx_1}{dt}  &= (A_{11} - A_{12}A_{22}^{-1}A_{21})x_1 + (B(1:n_1,:)- A_{12}A_{22}^{-1}B(n_1+1:end,:)) u   \\
 *  y           &= C \begin{bmatrix} x_1 \\ -A_{22}^{-1}A_{21} x_1 \end{bmatrix} -  C \begin{bmatrix} 0 \\ A_{22}^{-1}B(n_1+1:end,:) u \end{bmatrix}
 *  \end{aligned}
 *  \f]
 *
 *
 *
 * Depending on the @ref mess_operation_t defined in @ref mess_options.type we build an equation
 * structure to solve the following Riccati Equations:
 * If the operation type is @ref MESS_OP_NONE, we solve:
 * \f[
 * \begin{aligned}
 *  E_{11} X (A_{11}-A_{12}A_{22}^{-1}A_{21})^T + (A_{11}-A_{12}A_{22}^{-1}A_{21}) X E_{11}^T - E_{11}X C^T C X E_{11}^T + \tilde{B}\tilde{B}^T = 0 \\
 *  \tilde{B} := (B(1:n_1,:) - A_{12}A_{22}^{-1}B(n_1+1:end,:))
 * \end{aligned}
 * \f]
 * otherwise
 * \f[
 * \begin{aligned}
 *  E_{11}^T X (A_{11}-A_{12}A_{22}^{-1}A_{21}) + (A_{11}-A_{12}A_{22}^{-1}A_{21})^T X E_{11} - E_{11}^TX B B^T X E_{11} + \tilde{C}^T\tilde{C} = 0 \\
 *  \tilde{C} := (C(:,1:n_1) - C(:,n_1+1:end) A_{22}^{-1}A_{21}).
 * \end{aligned}
 * \f]
 *
 * \sa mess_lrcfadi_nm
 * \sa mess_equation_lyap
 * \sa mess_equation_st
 * \sa mess_equation_glyap_dae1
 *
 * References: \cite morGugSW13
 *
 */
int mess_equation_griccati_dae1(mess_equation e, mess_options opt , mess_matrix E11, mess_matrix A11, mess_matrix A12, mess_matrix A21, mess_matrix A22,  mess_matrix B, mess_matrix C) {
    MSG_FNAME(__func__);
    __glyap_dae1 *data ;
    int ret = 0;
    mess_operation_t optype;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(e);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(E11);    mess_check_real(E11); mess_check_square(E11);
    mess_check_nullpointer(A11);    mess_check_real(A11); mess_check_square(A11); mess_check_same_size(A11,E11);
    mess_check_nullpointer(A12);    mess_check_real(A12);
    mess_check_nullpointer(A21);    mess_check_real(A21); mess_check_same_colsrows(A12,A21); mess_check_same_colsrows(A21,A12);
    mess_check_nullpointer(A22);    mess_check_real(A22); mess_check_square(A22); mess_check_same_cols(A22,A12); mess_check_same_rows(A22,A21);
    mess_check_nullpointer(B);      mess_check_dense(B);
    mess_check_nullpointer(C);      mess_check_dense(C);

    optype = (opt->type==MESS_OP_NONE)?MESS_OP_NONE:MESS_OP_TRANSPOSE;

    if(optype==MESS_OP_NONE){
        mess_check_same_cols(E11,C);
        if(B->rows != (E11->rows+A22->rows)){
            MSG_ERROR("B has wrong number of rows="MESS_PRINTF_INT", expected="MESS_PRINTF_INT".\n", B->rows,E11->rows+A22->rows);
            return MESS_ERROR_ARGUMENTS;
        }
    }else{
        mess_check_same_rows(E11,B);
        if(C->cols != (E11->rows+A22->rows)){
            MSG_ERROR("C has wrong number of cols="MESS_PRINTF_INT", expected="MESS_PRINTF_INT".\n",C->cols,E11->rows+A22->rows);
            return MESS_ERROR_ARGUMENTS;
        }
    }

    if(B->cols >= B->rows){
        MSG_WARN("Oversized factor (op = %s, rows="MESS_PRINTF_INT", cols="MESS_PRINTF_INT")\n ", mess_operation_t_str(optype), B->rows, B->cols);
    }

    if (C->rows >= C->cols){
        MSG_WARN("Oversized factor (op = %s, rows="MESS_PRINTF_INT", cols="MESS_PRINTF_INT")\n ", mess_operation_t_str(optype), C->rows, C->cols);
    }

    /*-----------------------------------------------------------------------------
     * Set DEFAULT
     *-----------------------------------------------------------------------------*/
    if ( e->clear != NULL ) {
        e->clear(e->aux);
    }

    if ( e->clearRHS) {
        ret = mess_matrix_clear(&e->RHS);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
        e->RHS=NULL;
        e->clearRHS = 0 ;
    }

    if ( e->clearB) {
        ret = mess_matrix_clear(&e->B);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
        e->B=NULL;
        e->clearB = 0;
    }

    if ( e->clearC) {
        ret = mess_matrix_clear(&e->C);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
        e->C=NULL;
        e->clearC = 0;
    }


    /*-----------------------------------------------------------------------------
     *  build the __glyap_dae1 structure
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data,__glyap_dae1 *, sizeof(__glyap_dae1));
    data->E11           = E11;
    data->A11           = A11;
    data->A12           = A12;
    data->A21           = A21;
    data->A22           = A22;
    data->n1            = E11->rows;
    data->n2            = A22->rows;


    /*-----------------------------------------------------------------------------
     *  add direct solvers to data and matrices to data
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_init(&data->A22sol);                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
    ret = mess_direct_lu(data->A22,data->A22sol);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);

    ret = mess_direct_init(&data->E11sol);                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
    ret = mess_direct_lu(data->E11,data->E11sol);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);

    MESS_INIT_MATRICES(&(data->fullA));
    ret = mess_matrix_cat(A11,A12,A21,A22,MESS_CSR,data->fullA);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);
    ret = mess_direct_init(&(data->fullAsol));                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
    ret = mess_direct_lu(data->fullA,data->fullAsol);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);

    mess_matrix zeros;
    mess_matrix_init(&zeros);
    mess_matrix_init(&data->fullE);
    ret = mess_matrix_alloc(zeros,data->n2,data->n2,0,MESS_CSR,MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
    ret = mess_matrix_cat(E11,NULL,NULL,zeros,MESS_CSR,data->fullE);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);

    /*-----------------------------------------------------------------------------
     *  build equation structure
     *-----------------------------------------------------------------------------*/
    //set dimension
    e->dim                  = data->A11->rows;

    //set DAE 1 structure
    e->aux                  = data;

    //add function pointers
    e->eqn_type             = MESS_EQN_GRICCATI;
    e->clear                = clear;
    e->AX.apply             = AX_apply;
    e->EX.apply             = EX_apply;
    e->AINV.generate        = NULL;
    e->AINV.apply           = AINV_apply;
    e->AINV.clear           = NULL;
    e->EINV.generate        = NULL;
    e->EINV.apply           = EINV_apply;
    e->EINV.clear           = NULL;
    e->ApEINV.generate      = ApEINV2_generate;
    e->ApEINV.apply         = ApEINV2_apply;
    e->ApEINV.clear         = ApEINV2_clear;
    e->parameter            = parameter;
    e->init_rhs             = NULL;

    /*-----------------------------------------------------------------------------
     * compute reduced right hand side
     *-----------------------------------------------------------------------------*/
    mess_matrix myrhs,temp1,temp2;
    MESS_INIT_MATRICES(&myrhs,&temp1,&temp2);
    if(optype == MESS_OP_NONE){
        ret = mess_matrix_rowsub(B,0,data->n1-1,myrhs);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
        ret = mess_matrix_rowsub(B,data->n1,data->n1+data->n2-1,temp1);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
        ret = mess_direct_solvem(MESS_OP_NONE,data->A22sol,temp1,temp2);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
        ret = mess_matrix_multiply(MESS_OP_NONE,data->A12,MESS_OP_NONE,temp2,temp1);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        ret = mess_matrix_add(-1.0,temp1,1.0,myrhs);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);

        //set the reduced right hand side
        e->B            = myrhs;
        e->clearB       = 1;
        e->C            = C;
        e->clearC       = 0;


    }else{
        mess_matrix temp;
        MESS_INIT_MATRICES(&temp);
        mess_matrix_ctranspose(C,temp);                                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
        ret = mess_matrix_rowsub(temp,0,data->n1-1,myrhs);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
        ret = mess_matrix_rowsub(temp,data->n1,data->n1+data->n2-1,temp1);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
        ret = mess_direct_solvem(MESS_OP_TRANSPOSE,data->A22sol,temp1,temp2);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
        ret = mess_matrix_multiply(MESS_OP_TRANSPOSE,data->A21,MESS_OP_NONE,temp2,temp1);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        ret = mess_matrix_add(-1.0,temp1,1.0,myrhs);                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
        ret = mess_matrix_ctranspose(myrhs,temp1);                                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
        ret = mess_matrix_copy(temp1,myrhs);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        MESS_CLEAR_MATRICES(&temp);

        //set the reduced right hand side
        e->C            = myrhs;
        e->clearC       = 1;
        e->B            = B;
        e->clearB       = 0;

    }

    /*-----------------------------------------------------------------------------
     *  clear additional memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&temp1,&temp2,&zeros);
    return 0;
}

