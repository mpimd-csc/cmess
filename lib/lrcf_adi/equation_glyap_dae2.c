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
 * @file lib/lrcf_adi/equation_glyap_dae2.c
 * @brief Function handles for Hessenberg Index 2 DAE System.
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
 * @brief Contains matrices, solvers and additional data for a DAE 2 function handles.
 *
 * The @ref __glyap_dae2 structure contains all matrices, solvers and additional data for DAE 2 function handles.
 * @attention Internal use only.
 *
 */
typedef struct {
    mess_matrix M;              /**< System matrix \f$ M \f$.*/
    mess_matrix A;              /**< System matrix \f$ A \f$.*/
    mess_matrix G;              /**< System matrix \f$ G \f$.*/
    double delta;               /**< Shift parameter to avoid infiite eigenvalues of the matrix pencil.*/
    mess_matrix fullA;          /**< System block matrix \f$ \begin{bmatrix} A & G \\  G^T & 0 \end{bmatrix} \f$.*/
    mess_direct fullAsolver;    /**< Direct solver for \f$ \begin{bmatrix} A & G \\  G^T & 0 \end{bmatrix} \f$.*/
    mess_matrix fullM;          /**< System block matrix \f$ \begin{bmatrix} M & \delta G \\ \delta G^T & 0 \end{bmatrix} \f$ for shift paramter computations.*/
    mess_direct fullMsolver;    /**< Direct solver for \f$ \begin{bmatrix} M & \delta G \\ \delta G^T & 0 \end{bmatrix} \f$ for shift paramter computations.*/
    mess_direct applyPIsolver;  /**< Direct solver for \f$ \begin{bmatrix} M & \delta G \\ \delta G^T & 0 \end{bmatrix} \f$. */
    mess_matrix MM;             /**< System matrix \f$ \begin{bmatrix} M & 0 \\ 0 & 0\end{bmatrix} \f$. */
    mess_direct * shiftsolvers; /**< Array of direct solvers for shifted system \f$ \begin{bmatrix} A + p M & G \\ G^T & 0 \end{bmatrix} \f$ .*/
    mess_int_t num_solvers;     /**< Number of direct solvers in @ref __glyap_dae2.shiftsolvers. */

} __glyap_dae2;



/**
 * @internal
 * @brief Apply the projector \f$ \Pi:=I - G(G^TM^{-1}G)^{-1}G^TM^{-1}\f$ to a @ref mess_matrix.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref mess_matrix_applyPI_dae2 function applies the projector \f$ \Pi \f$ to a matrix, i.e. \n
 * it computes
 * \f[ out = op(\Pi) in = op(I - G(G^TM^{-1}G)^{-1}G^TM^{-1})in \f]
 * by solving
 * \f[
 *  \begin{aligned}
 *      \begin{bmatrix} M & G \\ G^T & 0  \end{bmatrix}
 *      \begin{bmatrix} M^{-1}out \\ \ast \end{bmatrix}
 *      &=
 *      \begin{bmatrix} in \\ 0 \end{bmatrix} \\
 *      \text{for}  \\
 *      \Pi in &= out \\
 *      \text{or}
 *      \begin{bmatrix} M & G \\ G^T & 0 \end{bmatrix}
 *      \begin{bmatrix} out \\ \ast \end{bmatrix}
 *      =
 *      \begin{bmatrix} M in \\ 0 \end{bmatrix}
 *      \text{for}
 *      \Pi^T in &= out.
 *  \end{aligned}
 * \f]
 *
 * @see mess_matrix_multiply
 * @see mess_matrix_add
 * @see mess_direct_solvem
 *
 * @attention Internal use only.
 */
static int mess_matrix_applyPI_dae2(mess_equation e, mess_operation_t op,  mess_matrix in, mess_matrix out){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    int ret;
    mess_int_t nv,np;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_operation_type(op);
    mess_check_nullpointer(in);
    mess_check_real_or_complex(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(out);
    mess_check_nullpointer(eqn);
    mess_check_same_colsrows(eqn->M,in);

    /*-----------------------------------------------------------------------------
     *  perform operation
     *-----------------------------------------------------------------------------*/

    nv = eqn->M->cols;
    np = eqn->G->cols;

    mess_matrix Tmp1;
    ret = mess_matrix_init(&Tmp1);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);

    if(op==MESS_OP_NONE){
        ret = mess_matrix_lift(in,np,Tmp1);                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_lift);
        ret = mess_direct_solvem(MESS_OP_NONE,eqn->applyPIsolver,Tmp1, out);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
        ret = mess_matrix_rowsub(out,0,nv-1,Tmp1);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
        ret = mess_matrix_multiply(MESS_OP_NONE, eqn->M, MESS_OP_NONE, Tmp1, out);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
    }else{
        ret = mess_matrix_multiply(MESS_OP_NONE, eqn->M, MESS_OP_NONE, in, Tmp1);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        ret = mess_matrix_lift(Tmp1,np,out);                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_lift);
        ret = mess_direct_solvem(MESS_OP_NONE,eqn->applyPIsolver,out,Tmp1);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
        ret = mess_matrix_rowsub(Tmp1,0,nv-1,out);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
    }

    ret = mess_matrix_clear(&Tmp1);                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    return ret;
}

/**
 * @internal
 * @brief Apply the projector \f$ op(A\Pi^T) \f$ to a @ref mess_matrix.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref AX_apply_T function computes \f$ out = op(A\Pi^T) in\f$.
 * This function is only used in residual computation, therefore it is assumed that \f$ in \f$  \f$ n_v \f$ rows.
 *
 * @see mess_matrix_applyPI_dae2
 * @see AX_apply_N
 *
 * @attention Internal use only.
 */
static int AX_apply_T(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    int ret=0;
    mess_int_t nv = eqn->A->cols;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);
    mess_check_operation_type(op);

    if(in->rows == nv){
        /*-----------------------------------------------------------------------------
         *  this part is only coded for residual computation
         *-----------------------------------------------------------------------------*/
        mess_matrix tmp1; MESS_INIT_MATRICES(&tmp1);

        //op(A*PI^T)*x->out
        if(op==MESS_OP_NONE){
            ret = mess_matrix_applyPI_dae2(e,MESS_OP_TRANSPOSE,in,tmp1);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_applyPI_dae2);
            ret = mess_matrix_multiply(MESS_OP_NONE, eqn->A, MESS_OP_NONE,tmp1, out);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        }else{
            ret = mess_matrix_multiply(MESS_OP_TRANSPOSE, eqn->A,MESS_OP_NONE,in,tmp1);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
            ret = mess_matrix_applyPI_dae2(e,MESS_OP_NONE,tmp1, out);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_applyPI_dae2);
        }
        MESS_CLEAR_MATRICES(&tmp1);

    }else{
        return MESS_ERROR_NOSUPPORT;
    }

    return 0;
}

/**
 * @internal
 * @brief Apply the projector \f$ op(\Pi A) \f$ to a @ref mess_matrix.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref AX_apply_N function computes \f$ out = op(\Pi A) in \f$.
 * This function is only used in residual computation, therefore it is assumed that \f$ in \f$ has \f$ n_v \f$ rows.
 *
 * @see mess_matrix_applyPI_dae2
 * @see AX_apply_N
 *
 * @attention Internal use only.
 */
static int AX_apply_N(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    int ret=0;
    mess_int_t nv = eqn->A->cols;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);
    mess_check_operation_type(op);

    if(in->rows == nv){
        /*-----------------------------------------------------------------------------
         *  this part is only coded for residual computation
         *-----------------------------------------------------------------------------*/
        mess_matrix tmp1; MESS_INIT_MATRICES(&tmp1);

        //op(PI*A)*x->out
        if(op==MESS_OP_NONE){
            ret = mess_matrix_multiply(MESS_OP_NONE, eqn->A, MESS_OP_NONE, in, tmp1);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
            ret = mess_matrix_applyPI_dae2(e,MESS_OP_NONE,tmp1,out);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_applyPI_dae2);
        }else{
            ret = mess_matrix_applyPI_dae2(e,MESS_OP_TRANSPOSE,in, tmp1);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_applyPI_dae2);
            ret = mess_matrix_multiply(MESS_OP_TRANSPOSE, eqn->A,MESS_OP_NONE,tmp1,out);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        }
        MESS_CLEAR_MATRICES(&tmp1);

    }else{
        return MESS_ERROR_NOSUPPORT;
    }

    return 0;
}

/**
 * @internal
 * @brief Apply the projector \f$ op(\Pi^T A \Pi) \f$ or \f$ \begin{bmatrix} A & G \\  G^T & 0 \end{bmatrix} \f$ to a @ref mess_matrix.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref AX_apply_shifts function computes \f$ op(\Pi^T A \Pi) \f$ or \f$ \begin{bmatrix} A & G \\  G^T & 0 \end{bmatrix} \f$.
 * If the number of rows of \f$ in \f$ is \f$ n_v\f$ then we perform
 * \f[ out = op(\Pi^T A \Pi)in \f],
 * because it is assumed that the context of usage is in a shift parameter computation in the case of
 * @ref MESS_LRCFADI_PARA_ADAPTIVE_V or @ref MESS_LRCFADI_PARA_ADAPTIVE_Z.
 * If the number of rows of \f$ in \f$ is \f$ n_v + n_p \f$ we perform
 * \f$  out = \begin{bmatrix} A & G \\  G^T & 0 \end{bmatrix} in  \f$.
 * In that case it is assumed that the function is used for @ref MESS_LRCFADI_PARA_MINMAX or @ref MESS_LRCFADI_PARA_MINMAX_REAL.
 *
 * @see mess_matrix_applyPI_dae2
 *
 * @attention Internal use only.
 */
static int AX_apply_shifts(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    int ret=0;
    mess_int_t nv = eqn->A->cols, np = eqn->G->cols;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);
    mess_check_operation_type(op);

    /*-----------------------------------------------------------------------------
     *  perform shift parameter computation
     *-----------------------------------------------------------------------------*/

    if(in->rows == nv){
        /*-----------------------------------------------------------------------------
         *  This part is only coded for adaptive shift parameter strategy
         *-----------------------------------------------------------------------------*/
        mess_matrix temp, temp2;
        MESS_INIT_MATRICES(&temp,&temp2);
        ret = mess_matrix_applyPI_dae2(e,MESS_OP_TRANSPOSE,in,temp);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_applyPI_dae2);
        ret = mess_matrix_multiply(op,eqn->A,MESS_OP_NONE,temp,temp2);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        ret = mess_matrix_applyPI_dae2(e,MESS_OP_NONE,temp2,out);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_applyPI_dae2);
        MESS_CLEAR_MATRICES(&temp,&temp2);
    }else if(in->rows == nv+np){
        /*-----------------------------------------------------------------------------
         *  This part is only coded for shift parameter computation
         *  (not adaptive shift parameter strategy)
         *-----------------------------------------------------------------------------*/
        ret = mess_matrix_multiply(op, eqn->fullA, MESS_OP_NONE, in, out);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
    }else{
        MSG_ERROR("Unknown usage for this function handle.\n");
        return MESS_ERROR_NOSUPPORT;
    }

    return 0;
}

/**
 * @internal
 * @brief Multiplication with \f$ M \f$ and a @ref mess_matrix.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref EX_apply function computes \f$ out=op(M)in \f$.
 * It is assumed that \f$ in \f$ has \f$ n_v \f$ rows.
 *
 * @attention Internal use only.
 */
static int EX_apply(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    int ret=0;
    mess_int_t nv;
    nv = eqn->M->cols;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);
    mess_check_operation_type(op);

    /*-----------------------------------------------------------------------------
     *  perform E*in -> out
     *-----------------------------------------------------------------------------*/
    if(in->rows ==nv){
        ret = mess_matrix_multiply(op,eqn->M,MESS_OP_NONE,in,out);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
    }else{
        return MESS_ERROR_NOSUPPORT;
    }

    return 0;
}

/**
 * @internal
 * @brief Multiplication with \f$ M \f$  or \f$ \begin{bmatrix} M & \delta G \\ \delta G^T & 0 \end{bmatrix} \f$ and a @ref mess_matrix.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref EX_apply_shifts function computes \f$ out=op(M)in \f$
 * or \f$ out = op\left(\begin{bmatrix} M & \delta G \\ \delta G^T & 0 \end{bmatrix}\right) in\f$.
 * If the number of rows of \f$ in \f$ is \f$ n_v \f$ we perform the first operation.
 * If the number of rows of \f$ in \f$ is \f$ n_v + n_p \f$ we perform the second operation.
 * Therfore it is assumed that \f$ in \f$ has \f$ n_v \f$ or \f$ n_v + n_p \f$ rows.
 * This function is only used for shift parameter computations.
 *
 * @attention Internal use only.
 */
static int EX_apply_shifts(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    int ret=0;
    mess_int_t nv, np;
    nv = eqn->M->cols;
    np = eqn->G->cols;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);
    mess_check_operation_type(op);

    /*-----------------------------------------------------------------------------
     *  perform E*in -> out
     *-----------------------------------------------------------------------------*/
    if(in->rows ==nv){
        ret = mess_matrix_multiply(op,eqn->M,MESS_OP_NONE,in,out);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
    }else if (in->rows == nv +np){
        ret = mess_matrix_multiply(op,eqn->fullM,MESS_OP_NONE,in,out);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
    }else{
        return MESS_ERROR_NOSUPPORT;
    }

    return 0;
}


/**
 * @internal
 * @brief Apply \f$ \left(\begin{bmatrix} A & G \\  G^T & 0 \end{bmatrix}\right)^{-1} \f$ to a @ref mess_matrix.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref AINV_apply_shifts function computes \f$ out = op \left(\begin{bmatrix} A & G \\  G^T & 0 \end{bmatrix}\right)^{-1} in \f$.
 * The number of rows of \f$ in \f$ must be \f$ n_v +n_p\f$.
 * This function should only be used for @ref MESS_LRCFADI_PARA_MINMAX or @ref MESS_LRCFADI_PARA_MINMAX_REAL shift parameter computations.
 *
 * @attention Internal use only.
 */
static int AINV_apply_shifts(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    int ret=0;
    mess_int_t nv, np;
    nv = eqn->A->cols;
    np = eqn->G->cols;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);
    mess_check_operation_type(op);

    if(in->rows != nv+np ){
        return MESS_ERROR_NOSUPPORT;
    }

    /*-----------------------------------------------------------------------------
     *  perform A\in = fullA\in -> out
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_solvem(op,eqn->fullAsolver,in,out);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);

    return 0;
}

/**
 * @internal
 * @brief Apply \f$ \left(\begin{bmatrix} M & \delta G \\  \delta G^T & 0 \end{bmatrix}\right)^{-1} \f$ to a @ref mess_matrix.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref EINV_apply_shifts function computes \f$ out = op \left(\begin{bmatrix} M & \delta G \\  \delta G^T & 0 \end{bmatrix}\right)^{-1} in \f$.
 * The number of rows of \f$ in \f$ must be \f$ n_v +n_p\f$.
 * This function should only be used for @ref MESS_LRCFADI_PARA_MINMAX or @ref MESS_LRCFADI_PARA_MINMAX_REAL shift parameter computations.
 *
 * @attention Internal use only.
 */
static int EINV_apply_shifts(mess_equation e, mess_operation_t op, mess_matrix in , mess_matrix out){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);

    int ret=0;

    mess_int_t nv, np;
    nv = eqn->A->cols;
    np = eqn->G->cols;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);

    if(in->rows != nv+np ){
        return MESS_ERROR_NOSUPPORT;
    }

    /*-----------------------------------------------------------------------------
     *  perform fullE\in -> out
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_solvem(op,eqn->fullMsolver,in,out);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);

    return 0;
}


/**
 * @internal
 * @brief Generate a shifted solver for the system matrix \f$ \begin{bmatrix} A + p M & G \\ G^T & 0 \end{bmatrix} \f$ (high memory variant).
 * @param[in] e             input @ref mess_equation structure
 * @param[in] parameters    input @ref mess_vector of shift parameters
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref ApEINV_generate_memory1 builds a @ref mess_multidirect solver for the matrices \n
 * \f[
 *      \begin{bmatrix} A + p M & G \\ G^T & 0 \end{bmatrix}.
 * \f]
 * It uses for each shift parameter a @ref mess_direct instance.
 *
 * @attention Internal use only.
 */
static int ApEINV_generate_memory1(mess_equation e, mess_vector parameters) {
    MSG_FNAME(__func__);

    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(parameters);
    if ( e->ApEINV.to_clear) return 0;

    if(!eqn->shiftsolvers){
        int ret = 0;
        mess_matrix temp;
        mess_try_alloc(eqn->shiftsolvers, mess_direct*,sizeof(mess_direct)*(parameters->dim));
        mess_int_t i =0;
        ret = mess_matrix_init(&temp);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
        if(MESS_IS_REAL(parameters)){
            //#warning parallelize code
            //#ifdef _OPENMP
            //#pragma omp parallel for private (i)
            //#endif
            for(i=0;i<parameters->dim;++i){
                ret = mess_matrix_copy(eqn->fullM, temp);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
                MSG_INFO("DAE2 Multishiftsolver: %d\t %e\n",(int) i,parameters->values[i]);
                ret = mess_matrix_add(1.0,eqn->fullA,parameters->values[i],temp);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
                ret = mess_direct_init(&(eqn->shiftsolvers[i]));            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
                ret = mess_direct_lu(temp,eqn->shiftsolvers[i]);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);
            }
            eqn->num_solvers=parameters->dim;
        }else{
            //#ifdef _OPENMP
            //#pragma omp parallel for private (i)
            //#endif
            for(i=0;i<parameters->dim;++i){
                ret = mess_matrix_copy(eqn->fullM, temp);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
                MSG_INFO("DAE2 Multishiftsolver: %d\t %e + %ei\n",(int)i,creal(parameters->values_cpx[i]),cimag(parameters->values_cpx[i]));
                ret = mess_matrix_addc(1.0,eqn->fullA,parameters->values_cpx[i],temp);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mes_matrix_addc);
                ret = mess_direct_init(&(eqn->shiftsolvers[i]));            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
                ret = mess_direct_lu(temp,eqn->shiftsolvers[i]);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);
            }
            eqn->num_solvers=parameters->dim;
        }
        ret = mess_matrix_clear(&temp);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    }
    e->ApEINV.to_clear = 1;
    return 0;
}

/**
 * @internal
 * @brief Solve shifted linear block system with matrix \f$ \begin{bmatrix} A + pM & G \\ G^T & 0 \end{bmatrix} \f$ (high memory variant).
 * @param[in, out] e        input @ref mess_equation structure
 * @param[in]   op          input operation type
 * @param[in]   p           input shift parameter
 * @param[in]   idx_p       input index of parameter for shift
 * @param[in]   in          input matrix
 * @param[out]  out         output matrix
 * @return zero on succes or a non-zero error value otherwise
 *
 *  The @ref ApEINV_apply_memory1 function solves the system
 * \f[
 *      op(\Pi A + p M ) out = in
 * \f]
 * via the equivalent block reformulation
 *
 *
 * \f[
 *      op\left(\begin{bmatrix} A + pM & G \\ G^T & 0 \end{bmatrix} \right)
 *      \begin{bmatrix} out \\ \ast \end{bmatrix}
 *      =
 *      \begin{bmatrix} in \\ 0 \end{bmatrix}.
 * \f]
 * The formulation is mathematical equivalent if \f$ in \in image(\Pi)\f$.
 *
 * @attention Internal use only.
 */
static int ApEINV_apply_memory1(mess_equation e, mess_operation_t op, mess_double_cpx_t p, mess_int_t idx_p, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    int ret=0;

    mess_int_t nv, np;
    nv = eqn->A->cols;
    np = eqn->G->cols;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);

    if(in->rows != nv ){
        //System size for lrfcadi/lrnm computation is nv
        MSG_ERROR("mess_matrix in has wrong size. ApEINV_apply in DAE2 is only coded for lrcfadi/lrnm computation.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  (fullA+p*MM)\[in;zeros(np,:)]->[out,*], with * is a block with np rows and not needed
     *-----------------------------------------------------------------------------*/
    mess_matrix temp1, temp2;
    ret = mess_matrix_init(&temp1);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&temp2);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_lift(in,np,temp1);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_lift);
    ret = mess_direct_solvem(op, eqn->shiftsolvers[idx_p],temp1,temp2);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
    ret = mess_matrix_rowsub(temp2,0,nv-1,out);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);

    ret = mess_matrix_clear(&temp1);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret = mess_matrix_clear(&temp2);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);

    return 0;
}

/**
 * @internal
 * @brief Clear a array of direct solvers for the system matrix \f$ \begin{bmatrix} A + pM & G \\ G^T & 0 \end{bmatrix} \f$ (high memory variant).
 * @param[in] e         input @ref mess_equation structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref ApEINV_clear_memory1 clear generated solvers for the matrices \n
 * \f[
 *      \begin{bmatrix} A + pM & G \\ G^T & 0 \end{bmatrix}.
 * \f]
 * by using @ref mess_direct_clear the fields of @ref __glyap_dae2.shiftsolvers are cleared.
 *
 * @attention Internal use only.
 *
 */
static int ApEINV_clear_memory1(mess_equation e){
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    MSG_FNAME(__func__);
    int ret;
    mess_check_nullpointer(eqn);
    if ( e->ApEINV.to_clear == 0) return 0;
    if(eqn->shiftsolvers){
        mess_int_t i=0;
        for(i=0;i<eqn->num_solvers;++i){
            ret = mess_direct_clear(&(eqn->shiftsolvers[i]));       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_clear);
        }
    }
    mess_free(eqn->shiftsolvers);
    eqn->shiftsolvers=NULL;

    e->ApEINV.to_clear = 0;
    return 0;
}

/**
 * @internal
 * @brief Dummy for generating direct solvers (low memory variant).
 * @param[in] e             input @ref mess_equation structure
 * @param[in] parameters        input shift parameters
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref ApEINV_generate_memory0 is only a dummy function which sets the
 * field. The decomposition for solving the shifted system is directly computed in
 * @ref ApEINV_apply_memory0.
 *
 * @attention Internal use only.
 */
static int ApEINV_generate_memory0(mess_equation e, mess_vector parameters) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    if ( e->ApEINV.to_clear) return 0;

    e->ApEINV.to_clear = 1;
    return 0;
}

/**
 * @internal
 * @brief Solve shifted linear block system with matrix \f$ \begin{bmatrix} A + pM & G \\ G^T & 0 \end{bmatrix} \f$ (low memory variant).
 * @param[in, out] e        input @ref mess_equation structure
 * @param[in]   op          input operation type
 * @param[in]   p           input shift parameter
 * @param[in]   idx_p       input index of parameter for shift
 * @param[in]   in          input matrix
 * @param[out]  out         output matrix
 * @return zero on succes or a non-zero error value otherwise
 *
 *  The @ref ApEINV_apply_memory0 function solves the system
 * \f[
 *      op(\Pi A + p M ) out = in
 * \f]
 * via the equivalent block reformulation
 *
 *
 * \f[
 *      op\left(\begin{bmatrix} A + pM & G \\ G^T & 0 \end{bmatrix} \right)
 *      \begin{bmatrix} out \\ \ast \end{bmatrix}
 *      =
 *      \begin{bmatrix} in \\ 0 \end{bmatrix}.
 * \f]
 * The formulation is mathematical equivalent if \f$ in \in image(\Pi)\f$.
 * Needed @ref mess_direct instances to solve the linear system are generated and cleared in that function.
 *
 * @attention Internal use only.
 */
static int ApEINV_apply_memory0(mess_equation e, mess_operation_t op, mess_double_cpx_t p, mess_int_t idx_p, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    int ret=0;

    mess_int_t nv, np;
    nv = eqn->A->cols;
    np = eqn->G->cols;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);
    mess_check_operation_type(op);

    if(in->rows != nv ){
        //System size for lrfcadi/lrnm computation is nv
        MSG_ERROR("mess_matrix in has wrong size. ApEINV_apply in DAE2 is only coded for lrcfadi/lrnm computation.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  (fullA+p*MM)\[in;zeros(np,:)]->[out,*], with * is a block with np rows and not needed
     *-----------------------------------------------------------------------------*/

    mess_direct solver;
    mess_matrix temp1, temp2, full;
    ret = mess_matrix_init(&temp1);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&temp2);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&full);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_lift(in,np,temp1);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_lift);
    ret = mess_matrix_copy(eqn->MM,full);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    ret = mess_matrix_addc(1.0,eqn->fullA,p,full);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mes_matrix_addc);
    ret = mess_direct_init(&solver);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
    ret = mess_direct_lu(full,solver);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);
    ret = mess_direct_solvem(op,solver,temp1,temp2);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
    ret = mess_matrix_sub(temp2,0,nv-1,0,in->cols-1,out);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_sub);
    ret = mess_direct_clear(&solver);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_clear);
    ret = mess_matrix_clear(&full);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret = mess_matrix_clear(&temp1);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret = mess_matrix_clear(&temp2);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);

    return 0;
}

/**
 * @internal
 * @brief Dummy function (low memory variant).
 * @param[in] e         input @ref mess_equation structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref ApEINV_clear_memory0 is only a dummy function because all data are cleared directly in
 * @ref ApEINV_apply_memory0.
 *
 * @attention Internal use only.
 */
static int ApEINV_clear_memory0(mess_equation e){
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    MSG_FNAME(__func__);
    mess_check_nullpointer(eqn);
    if ( e->ApEINV.to_clear == 0) return 0;
    e->ApEINV.to_clear = 0;
    return 0;
}

/**
 * @internal
 * @brief Call the @ref mess_lrcfadi_parameter and exchange function pointers.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] opt       input options
 * @param[in] stat      input status
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref parameter function
 * performs the shift parameter computation for the dae2 function handles.
 * In the case of @ref MESS_LRCFADI_PARA_ADAPTIVE_V or @ref MESS_LRCFADI_PARA_ADAPTIVE_Z
 * function pointers
 * are exchanged and the dimension of the equation is changed. After that
 * the @ref mess_lrcfadi_parameter function is called with the stabilized
 * or original equation.
 *
 * In the case of @ref MESS_LRCFADI_PARA_MINMAX or @ref MESS_LRCFADI_PARA_MINMAX_REAL
 * function pointers are exchanged the matrix \f$ B \f$ and \f$ K \f$ are lifted,
 * that the \f$(1,1)\f$ block of the matrix if stabilized.
 *
 * @see AX_apply_shifts
 * @see EX_apply_shifts
 * @see mess_lrcfadi_parameter
 *
 * @attention Internal use only.
 */
static int parameter(mess_equation e, mess_options opt, mess_status stat){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    int ret;
    mess_int_t nv = eqn->M->cols, np = eqn->G->cols;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(opt);
    mess_check_nullpointer(e);
    mess_check_nullpointer(stat);

    /*-----------------------------------------------------------------------------
     *  adaptive shift parameters need no lift for matrix B and K
     *-----------------------------------------------------------------------------*/
    if(opt->adi_shifts_paratype == MESS_LRCFADI_PARA_ADAPTIVE_V || opt->adi_shifts_paratype ==MESS_LRCFADI_PARA_ADAPTIVE_Z){
        //exchange function pointers
        int(*AX)(mess_equation,mess_operation_t,mess_matrix,mess_matrix) = e->AX.apply;
        e->AX.apply         = AX_apply_shifts;
        int(*EX)(mess_equation,mess_operation_t,mess_matrix,mess_matrix) = e->EX.apply;
        e->EX.apply         = EX_apply_shifts;

        //make sure system size is nv
        if(e->parent) e->parent->dim = nv;
        e->dim = nv;

        //perform shift parameter computation
        if(e->parent){
            ret = mess_lrcfadi_parameter(e->parent, opt, stat);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_lrcfadi_parameter);
        }else{
            ret = mess_lrcfadi_parameter(e, opt,stat);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_lrcfadi_parameter);
        }

        //reset function handles
        e->AX.apply         = AX;
        e->EX.apply         = EX;
    }else if(e->parent && e->parent->B && e->parent->K){
        /*-----------------------------------------------------------------------------
         *  lift up equation structure to space np+nv for shift parameter computation
         *  - check if operator was stabilized
         *  - store memory adress of original matrices
         *  - copy and lift original matrices
         *  - perform shift parameter computation
         *  - undo changes
         *-----------------------------------------------------------------------------*/
        //operator was stabilized
        //store memory address of original matrices
        mess_matrix B_old = e->parent->B, K_old = e->parent->K;

        // copy matrices
        mess_matrix B, K, Blift, Klift;
        MESS_INIT_MATRICES(&B,&K,&Blift,&Klift);
        ret = mess_matrix_copy(B_old,B);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_copy(K_old,K);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);

        //lift matrices
        ret = mess_matrix_lift(B,np,Blift);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_lift);
        ret = mess_matrix_ctranspose(K,Klift);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
        ret = mess_matrix_lift(Klift,np,K);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_lift);
        ret = mess_matrix_ctranspose(K,Klift);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);

        //set lifted matrices
        e->parent->B = Blift;
        e->parent->K = Klift;

        //set function handles related to shift parameter computation
        int(*AX)(mess_equation,mess_operation_t,mess_matrix,mess_matrix) = e->AX.apply;
        e->AX.apply         = AX_apply_shifts;
        int(*EX)(mess_equation,mess_operation_t,mess_matrix,mess_matrix) = e->EX.apply;
        e->EX.apply         = EX_apply_shifts;

        //set new system size
        e->parent->dim = nv+np;

        //perform shift parameter computation
        ret = mess_lrcfadi_parameter(e->parent, opt, stat);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_lrcfadi_parameter);

        //reset system size
        e->parent->dim = nv;

        //reset function handles
        e->AX.apply         = AX;
        e->EX.apply         = EX;

        //set matrices back to original matrices
        e->parent->B = B_old;
        e->parent->K = K_old;

        //clear copies
        MESS_CLEAR_MATRICES(&Blift,&Klift,&B,&K);

    }else{
        //operator was not stabilized
        //no matrices for lift

        //exchange function pointers
        int(*AX)(mess_equation,mess_operation_t,mess_matrix,mess_matrix) = e->AX.apply;
        e->AX.apply         = AX_apply_shifts;
        int(*EX)(mess_equation,mess_operation_t,mess_matrix,mess_matrix) = e->EX.apply;
        e->EX.apply         = EX_apply_shifts;

        //set new system size
        e->dim              = nv + np;

        //compute shifts
        ret = mess_lrcfadi_parameter(e,opt,stat);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_lrcfadi_parameter);

        //reset system size
        e->dim              = nv;

        //set function pointers back
        e->AX.apply         = AX;
        e->EX.apply         = EX;
    }

    return 0;
}

/**
 * @internal
 * @brief Apply \f$ \Pi \f$ to the right hand side of an equation.
 * @param[in, out] e        input @ref mess_equation structure
 * @param[in] opt           input @ref mess_options structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref init_rhs function applies \f$ \Pi \f$ to the
 * @ref mess_matrix @ref mess_equation.RHS.
 * If the equation is stabilized, i.e. @ref mess_equation.parent
 * is not @c NULL then \f$ \Pi \f$ is applied to the parent equation,
 * otherwise to the original one.
 *
 * @see mess_matrix_applyPI_dae2
 *
 * @attention Internal use only.
 */
static int init_rhs(mess_equation e, mess_options opt){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);
    int ret;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(e);
    mess_check_nullpointer(opt);

    /*-----------------------------------------------------------------------------
     *  perform RHS = PI*RHS (projection) in parent equation if given
     *-----------------------------------------------------------------------------*/
    if(e->parent){
        //project RHS in parent equation structure
        mess_check_nullpointer(e->parent->RHS);
        ret = mess_matrix_applyPI_dae2(e,MESS_OP_NONE,e->parent->RHS,e->parent->RHS);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_applyPI);
    }else{
        //project RHS in original equation structure
        mess_check_nullpointer(e->RHS);
        ret = mess_matrix_applyPI_dae2(e,MESS_OP_NONE,e->RHS,e->RHS);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_applyPI);
    }
    return 0;
}

/**
 * @internal
 * @brief Clear additional allocated field of @ref __glyap_dae2 structure.
 * @param[in, out] e         input @ref mess_equation structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref clear function, clears the following additional allocated data
 *
 * <center>
 * |          field                           |  mathematical description                                                                       |
 * |:----------------------------------------:|:-----------------------------------------------------------------------------------------------:|
 * | @ref __glyap_dae2.fullA                  | \f$ \begin{bmatrix} A & G \\  G^{T} & 0 \end{bmatrix} \f$                                       |
 * | @ref __glyap_dae2.fullM                  | \f$ \begin{bmatrix} M & \delta G \\ \delta G^T & 0 \end{bmatrix} \f$                            |
 * | @ref __glyap_dae2.MM                     | \f$ \begin{bmatrix} M & 0 \\ 0 & 0\end{bmatrix}  \f$                                            |
 * | @ref __glyap_dae2.fullAsolver            | \f$ \left( \begin{bmatrix} A & G \\  G^T & 0 \end{bmatrix}\right)^{-1} \f$                      |
 * | @ref __glyap_dae2.fullMsolver            | \f$ \begin{bmatrix} M & \delta G \\ \delta G^T & 0 \end{bmatrix} \f$                            |
 * | @ref __glyap_dae2.applyPIsolver          | \f$ \left( \begin{bmatrix} M & \delta G \\ \delta G^T & 0 \end{bmatrix}\right)^{-1} \f$         |
 * </center>.
 *
 * @attention Internal use only.
 */
static int clear(mess_equation e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_dae2 *eqn = (__glyap_dae2 *) e->aux;
    mess_check_nullpointer(eqn);

    //clear additional matrices
    MESS_CLEAR_MATRICES(&(eqn->fullA),&(eqn->fullM),&(eqn->MM));

    //clear additional solvers
    MESS_CLEAR_DIRECTS(&(eqn->fullAsolver),&(eqn->fullMsolver),&(eqn->applyPIsolver));

    mess_free(eqn);
    eqn = NULL;
    e->aux = NULL;
    //e->getritzvalues = NULL;
    e->parameter = NULL;
    return 0;
}

/**
 * @brief Generate an Equation object to solve a generalized Lyapunov Equation from Hessenberg Index 2 DAE.
 * @param[in,out]   e        equation object
 * @param[in]       opt           input options for ADI iteration
 * @param[in]       M         input \f$ M \f$ matrix of the DAE system
 * @param[in]       A         input \f$ A \f$ matrix of the DAE system
 * @param[in]       G         input \f$ G \f$ matrix for the linear algebraic constraint
 * @param[in]       rhs       input \f$ B \f$ or \f$ C \f$ (unprojected) right hand side of the equation
 * @param[in]       delta     input \f$ \delta \f$ for the infinite eigenvalues of the matrix pencil
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_equation_glyap_dae2 function creates a mess_equation object from a Hessenberg Index 2 DAE System
 *
 *  \f[
 *  \begin{aligned}
 *  \begin{bmatrix} M & 0 \\ 0 & 0  \end{bmatrix}
 *  \begin{bmatrix} \frac{dz}{dt} \\ 0  \end{bmatrix}
 *  =
 *  \begin{bmatrix} A & G \\ G^T & 0 \end{bmatrix}
 *  \begin{bmatrix} z \\ p \end{bmatrix}
 *  +
 *  \begin{bmatrix} B \\ 0 \end{bmatrix}
 *  u
 *  \end{aligned}
 *  \f]
 *
 * Preliminaries:
 *  \f$ A \in \mathbb{R}^{n_v \times n_v}, B \in \mathbb{R}^{n_v \times p}, C \in \mathbb{R}^{n_a \times n_v} \f$. \n
 *  \f$ n_p < n_v \f$. \n
 *  \f$ \delta \in \mathbb{R}\f$ and \f$ \delta <0 \f$. \n
 *  \f$ M \in \mathbb{R}^{n_v \times n_v} \f$ must be symmetric and positive definite. \n
 *  \f$ G \in \mathbb{R}^{n_v \times n_p} \f$ and \f$ G \f$ has full rank \f$ n_p \f$. \n
 *
 *  The matrix pencil
 *  \f[
 *  \left(
 *  \begin{bmatrix} A & G \\ G^T & 0 \end{bmatrix}
 *  ,
 *  \begin{bmatrix} M & 0 \\ 0 & 0 \end{bmatrix}
 *  \right)
 *  \f]
 *  belongs to that system.
 *  To avoid infinite eigenvalues, we use
 *  \f[
 *  \left(
 *  \begin{bmatrix}A & G \\ G^{T} & 0 \end{bmatrix},
 *  \begin{bmatrix}M & \delta G \\ \delta G^{T} & 0 \end{bmatrix}
 *  \right)
 *  \f]
 *  as matrix pencil for shift parameter computation in \ref mess_lrcfadi_parameter.
 *  This transforms all infinite eigenvalues to \f$ \displaystyle \frac{1}{\delta} \f$. \n
 *  Depending on the @ref mess_options.type  component  a normal or transposed projected Lyapunov Equation is created. \n
 *  The discrete Leray Projector is applied implicitly. \n
 *  For more details look at \cite BaeBSetal14 \cite BaeBSetal15.
 *
 *  Depending on the @ref mess_options.memory_usage component of the \ref mess_options structure the memory usage is controlled. \n
 *   *  @ref mess_options.memory_usage  = 0:
 *  In every iteration a new solver is computed in ApEINV_apply.
 *   *  @ref mess_options.memory_usage  = 1:
 *  For every shift a solver is precomputed and reused in ApEINV_apply.
 *
 *   * @ref mess_options.type = MESS_OP_NONE:
 *      \f$ (\Pi A) X  M^T  +  M  X  (\Pi A)^T = -(\Pi B)(\Pi B)^T \f$
 *   * @ref mess_options.type = MESS_OP_TRANSPOSE:
 *      \f$ (\Pi A)^T X  M  +  M^T  X  (\Pi A) = -(C\Pi^T )^T(C\Pi^T ) \f$
 *
 * If your system is unstable use \ref mess_equation_stable to build the stabilized operator
 * \f$ (A-BK^T) \f$ or \f$ (A-KB^T) \f$, to ensure the proper application of \f$ \Pi \f$ also after stabilization
 * This function sets the \f$ B \f$ component of the given \ref mess_equation_st structure to
 *
 * @ref mess_options.type = MESS_OP_NONE:
* \f[ \Pi B \f]
*
* @ref mess_options.type = MESS_OP_TRANSPOSE:
* \f[ B \Pi^T. \f]
*
*
* \sa mess_lrcfadi_adi
* \sa mess_equation_lyap
* \sa mess_equation_st
* \sa mess_equation_stable
*
*/
int mess_equation_glyap_dae2(mess_equation e, mess_options opt , mess_matrix M, mess_matrix A, mess_matrix G, mess_matrix rhs, double delta) {
    MSG_FNAME(__func__);
    __glyap_dae2 *data ;
    int ret = 0;
    mess_operation_t optype;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(e);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(M);
    mess_check_nullpointer(A);
    mess_check_nullpointer(G);
    mess_check_nullpointer(rhs);

    mess_check_real(M); mess_check_square(M); mess_check_same_size(M,A);
    mess_check_real(A); mess_check_square(A);
    mess_check_real(G); mess_check_same_rows(A,G);
    mess_check_real(rhs); mess_check_dense(rhs);

    optype = opt->type;
    if (optype == MESS_OP_HERMITIAN) optype=MESS_OP_TRANSPOSE;

    if(optype==MESS_OP_NONE){
        mess_check_same_rows(G,rhs);
        if ( rhs->cols >= rhs->rows) {
            MSG_WARN("Oversized right hand side factor for LRCF-ADI (op = %s, rows=%d, cols=%d)\n ", mess_operation_t_str(optype), (int) rhs->rows, (int) rhs->cols);
        }
    }else{
        mess_check_same_colsrows(rhs,G);
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
     *  build the __glyap_dae2 structure
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data,__glyap_dae2 *, sizeof(__glyap_dae2));
    data->M             = M;
    data->A             = A;
    data->G             = G;
    data->delta         = delta;

    //add additional matrices for solving the saddle point problem
    mess_matrix Gt,Gcopy, applyPI;
    MESS_INIT_MATRICES(&Gt,&Gcopy,&data->fullA,&data->fullM,&data->MM,&applyPI);
    ret = mess_matrix_copy(G,Gcopy);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    ret = mess_matrix_ctranspose(G,Gt);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);

    //build fullA=[A,G;G^T,0]
    ret = mess_matrix_cat(A,G,Gt,NULL,MESS_CSR,data->fullA);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);

    //build applyPI=[M,G;G^T,0]
    ret = mess_matrix_cat(M,G,Gt,NULL,MESS_CSR,applyPI);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);

    //build fullE=[M,delta*G;delta*G^T,0]
    ret = mess_matrix_scale(data->delta,Gt);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
    ret = mess_matrix_scale(data->delta,Gcopy);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
    ret = mess_matrix_cat(data->M,Gcopy,Gt,NULL,MESS_CSR,data->fullM);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);

    //build MM=[M,0;0,0]
    ret = mess_matrix_zeros(Gcopy);                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_zeros);
    ret = mess_matrix_zeros(Gt);                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_zeros);
    ret = mess_matrix_cat(M,Gcopy,Gt,NULL,MESS_CSR,data->MM);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_zeros);


    //init solvers
    data->shiftsolvers  = NULL;
    MESS_INIT_DIRECTS(&data->fullAsolver,&data->applyPIsolver,&data->fullMsolver);
    ret = mess_direct_lu(data->fullA,data->fullAsolver);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);
    ret = mess_direct_lu(applyPI,data->applyPIsolver);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);
    ret = mess_direct_lu(data->fullM,data->fullMsolver);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);

    /*-----------------------------------------------------------------------------
     *  build equation structure
     *-----------------------------------------------------------------------------*/
    //set dimension
    e->dim                  = M->cols;

    //set DAE 2 structure
    e->aux                  = data;

    //add function pointers
    e->eqn_type             = MESS_EQN_GLYAP;
    e->clear                = clear;
    e->AX.apply             = (optype==MESS_OP_NONE) ? AX_apply_N : AX_apply_T;
    e->EX.apply             = EX_apply;
    e->AINV.generate        = NULL;
    e->AINV.apply           = AINV_apply_shifts;
    e->AINV.clear           = NULL;
    e->EINV.generate        = NULL;
    e->EINV.apply           = EINV_apply_shifts;
    e->EINV.clear           = NULL;
    if(opt->memory_usage==0){
        e->ApEINV.generate  = ApEINV_generate_memory0;
        e->ApEINV.apply     = ApEINV_apply_memory0;
        e->ApEINV.clear     = ApEINV_clear_memory0;
    }else{
        e->ApEINV.generate  = ApEINV_generate_memory1;
        e->ApEINV.apply     = ApEINV_apply_memory1;
        e->ApEINV.clear     = ApEINV_clear_memory1;
    }
    e->parameter            = parameter;
    e->init_rhs             = init_rhs;

    /*-----------------------------------------------------------------------------
     *  Set RHS Project the RHS
     * Warning
     * Use PIB instead of B for e->B, because in mess_lrcfadi_residual
     * operation PI(A-B*K)*x = PI*Ax-PI*B*K*x = Ax.apply(x)-PIB*(K*x) -> y
     * with the stabilized System has to be performed.
     * using e->B = B, would have the following effect
     * PI(A-B*K)*x != PI*Ax-B*K*x = Ax.apply(x)-B*(K*x) -> y
     * For the Transpose Lyapunov Equation the argumentation is similiar and we use
     * B*PI^T.
     * \sa equation_bk.c Ax_apply
     *-----------------------------------------------------------------------------*/
    mess_matrix PIB,Bt;
    if(optype == MESS_OP_NONE){
        // PI*B -> e->B
        ret = mess_matrix_init(&PIB);                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
        ret = mess_matrix_applyPI_dae2(e,MESS_OP_NONE,rhs,PIB);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_applyPI);
        //set the Projected B to e->B
        e->B            = PIB;
        e->clearB       = 1; //PIB is cleared by mess_equation_clear

        e->RHS          = rhs; //projection is done in  mess_lradi by init_rhs
        e->clearRHS     = 0; //B is cleared by user
    }else{
        //  C*PI^T -> e->B
        MESS_INIT_MATRICES(&PIB,&Bt,&(e->B));
        ret = mess_matrix_ctranspose(rhs,Bt);                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
        ret = mess_matrix_applyPI_dae2(e,MESS_OP_NONE,Bt,PIB);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_applyPI);
        ret = mess_matrix_ctranspose(PIB,e->B);                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
        e->clearB       = 1; //e->B is cleared by mess_equation_clear

        e->RHS          = Bt; //projection is done in mess_lradi by init_rhs
        e->clearRHS     = 1;  //cleared by mess_equation_clear
        ret = mess_matrix_clear(&PIB);                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    }

    /*-----------------------------------------------------------------------------
     *  clear additional memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&Gt,&Gcopy,&applyPI);
    return 0;
}

/**
 *
 * @brief Generate an equation object to a Riccati Equation arising from a Hessenberg Index 2 DAE.
 * @param[out] e         equation object to setup
 * @param[in]  opt        input options for the ADI iteration
 * @param[in]  M          input \f$ M \f$ matrix of the DAE system
 * @param[in]  A          input \f$ A \f$ matrix of the DAE system
 * @param[in]  G          input \f$ G \f$ matrix for the linear algebraic constraint
 * @param[in]  B          input right hand side \f$ B \f$ for the control
 * @param[in]  C          input observer matrix for the output y
 * @param[in] delta       input shifting parameter for the infinite eigenvalues of the matrix pencil
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_equation_griccati_dae2 function creates a mess_equation object from a Hessenberg Index 2 DAE System
 *
 *  \f[
 *  \begin{aligned}
 *  \begin{bmatrix} M & 0 \\ 0 & 0 \end{bmatrix}
 *  \begin{bmatrix} \frac{dz}{dt} \\ 0 \end{bmatrix}
 *  &=
 *  \begin{bmatrix} A & G \\ G^T & 0 \end{bmatrix}
 *  \begin{bmatrix} z \\ p \end{bmatrix} \\
 *  y &= C z
 *  \end{aligned}
 *  \f]
 *
 * Preliminaries:
 *  @li \f$ A \in \mathbb{R}^{n_v \times n_v}, B \in \mathbb{R}^{n_v \times n_r}, C \in \mathbb{R}^{n_a \times n_v} \f$
 *  @li \f$ n_v > n_p \f$
 *  @li matrix \f$ M \in \mathbb{R}^{n_v \times n_v} \f$ need to be positive definite.
 *  @li \f$ G \in \mathbb{R}^{n_v \times n_p} \f$ and \f$ G\f$ has column rank \f$ n_p \f$
 *
 *  The following matrix pencil belongs to that system
 *  \f[
 *  \left(
 *  \begin{bmatrix}A & G \\ G^{T} & 0 \end{bmatrix},
 *  \begin{bmatrix}M & 0 \\ 0 & 0 \end{bmatrix}
 *  \right).
 *  \f]
 *  To avoid infinite eigenvalues, we use
 *  \f[
 *  \left(
 *  \begin{bmatrix}A & G \\ G^{T} & 0 \end{bmatrix},
 *  \begin{bmatrix}M & \delta G \\ \delta G^{T} & 0 \end{bmatrix}
 * \right)
 *  \f]
 *  as matrix pencil for shift parameter computation in \ref mess_lrcfadi_parameter.
 *  This transforms all infinite eigenvalues to \f$ \displaystyle \frac{1}{\delta} \f$.
 *
 *  Depending on the component \ref mess_options.type  a normal or transposed \f$ \Pi^{T} \f$
 *  -projected Riccati Equation is created. \n
 *  The discrete Leray projector is applied implicitly.
 *  For more details look at \cite BaeBSetal14 \cite BaeBSetal15 .
 *
 *  Depending on memory_usage component of the \ref mess_options structure the memory usage is controlled.
 *  @ref mess_options.memory_usage = 0:
 *  In every iteration a new solver is computed in ApEINV_apply.
 *  @ref mess_options.memory_usage = 1:
 *  For every shift a solver is precomputed and reused in ApEINV_apply.
 *
 *  @ref mess_options.type = MESS_OP_NONE:
 *  \f[
 *  \Pi A \Pi^T X \Pi M^T \Pi^T + \Pi M \Pi^T X \Pi A^T \Pi^T - \Pi M \Pi^T X \Pi C^T C \Pi^T X \Pi M^T \Pi^T + \Pi B B^T \Pi^T = 0
 *  \f]
 *  @ref mess_options.type = MESS_OP_TRANSPOSE:
 *  \f[
 *  \Pi A^T \Pi^T X \Pi M \Pi^T + \Pi M^T \Pi^T X \Pi A \Pi^T - \Pi M^T \Pi^T X \Pi B B^T \Pi^T X \Pi M \Pi^T + \Pi C^T C \Pi^T = 0
 *  \f]
 *
 * If your system is unstable use the @p K0 component of \ref mess_options structure to provide an initial
 * stabilizing feedback for \ref mess_lrcfadi_nm. \n
* This function sets the \f$ B \f$ and \f$ C \f$ component of the given \ref mess_equation_st
* structure.
*
* @ref mess_options.type = MESS_OP_NONE:
*
* \f[
    * \begin{array}{ccc}
    *  B &=& \Pi B \\
           *  C &=& C
           *  \end{array}
           * \f]
           * @ref mess_options.type = MESS_OP_TRANSPOSE:
           *
           * \f[
               * \begin{array}{ccc}
               *  B &=& B \\
                      *  C &=& C \Pi^T.
                      *  \end{array}
                      * \f]
                      *
                      *
                      * \sa mess_equation_lyap
                      * \sa mess_equation_st
                      * \sa mess_lrcfadi_nm
                      *
    */
int mess_equation_griccati_dae2(mess_equation e, mess_options opt , mess_matrix M, mess_matrix A, mess_matrix G, mess_matrix B, mess_matrix C, double delta )
{
    MSG_FNAME(__func__);
    __glyap_dae2 *data ;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(e);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(M);
    mess_check_nullpointer(A);
    mess_check_nullpointer(G);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);

    mess_check_real(M); mess_check_square(M); mess_check_same_size(M,A);
    mess_check_real(A); mess_check_square(A);
    mess_check_real(G); mess_check_same_rows(A,G);
    mess_check_real(B); mess_check_same_rows(G,B);
    mess_check_real(C); mess_check_same_cols(A,C);
    mess_check_dense(C);
    mess_check_dense(B);
    mess_check_same_colsrows(C,B);

    /*-----------------------------------------------------------------------------
     * Set DEFAULT
     *-----------------------------------------------------------------------------*/
    if ( e->clear != NULL ) {
        e->clear(e->aux);
    }

    if ( e->clearRHS) {
        mess_matrix_clear(&e->RHS);
        e->RHS=NULL;
        e->clearRHS = 0 ;
    }

    if ( e->clearB) {
        mess_matrix_clear(&e->B);
        e->B=NULL;
        e->clearB = 0;
    }

    if ( e->clearC) {
        mess_matrix_clear(&e->C);
        e->C=NULL;
        e->clearC = 0;
    }

    /*-----------------------------------------------------------------------------
     *  Setup the equation
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data,__glyap_dae2 *, sizeof(__glyap_dae2));

    //add system matrices
    data->M             = M;
    data->A             = A;
    data->G             = G;
    data->delta         = delta;        //parameter for infinite shifts

    //add additional matrices for solving the saddle point problem

    //build fullA=[A,G;G^T,0]
    mess_matrix Gt,Gcopy, applyPI;
    MESS_INIT_MATRICES(&Gt,&Gcopy,&(data->fullA),&(applyPI),&(data->fullM),&(data->MM));
    ret = mess_matrix_copy(G,Gcopy);                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    ret = mess_matrix_ctranspose(G,Gt);                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    ret = mess_matrix_cat(A,G,Gt,NULL,MESS_CSR,data->fullA);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);

    //build applyPI=[M,G;G^T,0]
    ret = mess_matrix_cat(M,G,Gt,NULL,MESS_CSR,applyPI);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);

    //build fullM=[M,delta*G;delta*G^T,0]
    ret = mess_matrix_scale(data->delta,Gt);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
    ret = mess_matrix_scale(data->delta,Gcopy);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
    ret = mess_matrix_cat(data->M,Gcopy,Gt,NULL,MESS_CSR,data->fullM);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);

    //build MM=[M,0;0,0]
    ret = mess_matrix_zeros(Gcopy);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_zeros);
    ret = mess_matrix_zeros(Gt);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_zeros);
    ret = mess_matrix_cat(M,Gcopy,Gt,NULL,MESS_CSR,data->MM);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_zeros);

    //add solvers
    data->shiftsolvers  = NULL;
    MESS_INIT_DIRECTS(&(data->fullAsolver),&(data->applyPIsolver),&(data->fullMsolver));

    ret = mess_direct_lu(data->fullA,data->fullAsolver);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);
    ret = mess_direct_lu(applyPI,data->applyPIsolver);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);
    ret = mess_direct_lu(data->fullM,data->fullMsolver);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);

    e->dim                  = M->rows;
    e->aux                  = data;
    e->RHS                  = NULL;

    //add function pointers
    e->eqn_type             = MESS_EQN_GRICCATI;
    e->clear                = clear;
    e->AX.apply             = (opt->type==MESS_OP_NONE) ? AX_apply_N : AX_apply_T;
    e->EX.apply             = EX_apply;
    e->AINV.generate        = NULL;
    e->AINV.apply           = AINV_apply_shifts;
    e->AINV.clear           = NULL;
    e->EINV.generate        = NULL;
    e->EINV.apply           = EINV_apply_shifts;
    e->EINV.clear           = NULL;
    if(opt->memory_usage==0){
        e->ApEINV.generate      = ApEINV_generate_memory0;
        e->ApEINV.apply         = ApEINV_apply_memory0;
        e->ApEINV.clear         = ApEINV_clear_memory0;
    }else{
        e->ApEINV.generate      = ApEINV_generate_memory1;
        e->ApEINV.apply         = ApEINV_apply_memory1;
        e->ApEINV.clear         = ApEINV_clear_memory1;
    }
    e->ApEX.generate    = NULL;
    e->ApEX.apply       = NULL;
    e->ApEX.clear       = NULL;
    e->parameter        = parameter;
    e->init_rhs         = init_rhs;

    if(opt->type == MESS_OP_NONE){
        //Project B RHS
        //PI*B -> eqn->B
        mess_matrix PIB;
        ret = mess_matrix_init(&PIB);                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
        ret = mess_matrix_applyPI_dae2(e,MESS_OP_NONE,B,PIB);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_applyPI);
        e->B            = PIB;
        e->clearB       = 1;

        e->C            = C;
        e->clearC       = 0;
    }else{
        //Project C RHS
        //C*PI^T -> eqn->C
        mess_matrix Tmp1,Tmp2;
        MESS_INIT_MATRICES(&(Tmp1),&(Tmp2));
        ret = mess_matrix_ctranspose(C,Tmp1);                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
        ret = mess_matrix_applyPI_dae2(e,MESS_OP_NONE,Tmp1,Tmp2);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_applyPI);
        ret = mess_matrix_ctranspose(Tmp2,Tmp1);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
        e->C            = Tmp1;
        e->clearC       = 1;
        ret = mess_matrix_clear(&Tmp2);                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);

        e->B            = B;
        e->clearB       = 0;
    }

    /*-----------------------------------------------------------------------------
     *  clear memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&Gt,&Gcopy,&applyPI);

    return 0;
}





