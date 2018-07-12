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
 * @file lib/lrcf_adi/equation_glyap.c
 * @brief  Build an equation object (generalized Lyapunov or Riccati Equation).
 * @author @koehlerm
 *
 * Generate the @ref mess_equation object for the generalized Lyapunov and Riccati equatuion.
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
 * @brief Contains matrices, solvers and additional data for a generalized Lyapunov Equation.
 *
 * The @ref __lyap structure contains all matrices, solvers and additional data
 * for a generalized Lyapunov Equation
 * \f[AXE^T +EAF^T = -BB^T\f]
 * or
 * \f[A^TXE +E^TXA = -B^TB.\f]
 * @attention Internal use only.
 *
 */
typedef struct {
    mess_matrix A;                  /**< System matrix \f$A\f$ of Lyapunov Equation. */
    mess_matrix E;                  /**< System matrix \f$E\f$ of Lyapunov Equation. */
    mess_direct Asolver;            /**< Direct solver for system matrix \f$A\f$. */
    mess_direct Esolver;            /**< Direct solver for system matrix \f$E\f$. */
    mess_multidirect Amsolver;      /**< Multisolver for shifted matrix \f$ A+pE\f$. */
    mess_matrix ApEX_tmp;           /**< Temporary shifted matrix \f$A+pI\f$. */
} __lyap;




/**
 * @internal
 * @brief Compute matrix-matrix product with system matrix \f$A\f$.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to \f$A\f$
 * @param[in] x         input  matrix \f$x\f$
 * @param[out] y        output matrix \f$y\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref AX_apply function computes the matrix-matrix product, i.e.
 * \f[ y = op(A) x. \f]
 *
 * @see mess_matrix_multiply
 * @attention Internal use only.
 **/
static int AX_apply(mess_equation e, mess_operation_t op, mess_matrix x, mess_matrix y) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    switch(op){
        case MESS_OP_NONE:
            return mess_matrix_multiply(MESS_OP_NONE, eqn->A, MESS_OP_NONE, x, y);
        case MESS_OP_TRANSPOSE:
            return mess_matrix_multiply(MESS_OP_HERMITIAN, eqn->A, MESS_OP_NONE, x, y);
        case MESS_OP_HERMITIAN:
            return mess_matrix_multiply(MESS_OP_HERMITIAN, eqn->A, MESS_OP_NONE, x, y);
        default:
            return MESS_ERROR_ARGUMENTS;
    }
    return 0;
}

/**
 * @internal
 * @brief Compute matrix-matrix product with system matrix \f$E\f$.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to \f$E\f$
 * @param[in] x         input matrix \f$x\f$
 * @param[out] y        output matrix \f$y\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref EX_apply function computes the matrix-matrix product, i.e.
 * \f[ y = op(E) x. \f]
 *
 * @see mess_matrix_multiply
 * @attention Internal use only.
 **/
static int EX_apply(mess_equation e, mess_operation_t op, mess_matrix x, mess_matrix y) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    switch(op){
        case MESS_OP_NONE:
            return mess_matrix_multiply(MESS_OP_NONE, eqn->E, MESS_OP_NONE, x, y);
        case MESS_OP_TRANSPOSE:
            return mess_matrix_multiply(MESS_OP_HERMITIAN, eqn->E, MESS_OP_NONE, x, y);
        case MESS_OP_HERMITIAN:
            return mess_matrix_multiply(MESS_OP_HERMITIAN, eqn->E, MESS_OP_NONE, x, y);
        default:
            return MESS_ERROR_ARGUMENTS;
    }
    return 0;
}

static int AINV_generate(mess_equation e){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    int ret = 0;
    if ( e->AINV.to_clear) return 0;
    if ( !eqn->Asolver ) {
        ret = mess_direct_init(&(eqn->Asolver));    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
        ret = mess_direct_lu(eqn->A, eqn->Asolver); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);
    }

    e->AINV.to_clear = 1;
    return 0;
}

static int AINV_apply(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    return mess_direct_solvem(op, eqn->Asolver, in, out);
}

static int AINV_clear(mess_equation e){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    if (eqn->Asolver) mess_direct_clear(&(eqn->Asolver));
    eqn->Asolver = NULL;
    //e->Ainv.to_clear = 0;
    e->AINV.to_clear = 0;
    return 0;

}

static int EINV_generate(mess_equation e){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    int ret = 0;
    if ( e->EINV.to_clear) return 0;
    ret = mess_direct_init(&(eqn->Esolver));    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
    ret = mess_direct_lu(eqn->E, eqn->Esolver); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);
    e->EINV.to_clear = 1;
    return 0;
}

static int EINV_apply(mess_equation e, mess_operation_t op, mess_matrix in , mess_matrix out){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    return mess_direct_solvem(op, eqn->Esolver, in, out);
}

static int EINV_clear(mess_equation e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    if ( e->EINV.to_clear == 0 ) return 0;
    mess_direct_clear(&(eqn->Esolver));
    eqn->Esolver = NULL;
    e->EINV.to_clear = 0;
    return 0;
}

static int ApEINV_generate(mess_equation e, mess_vector parameters) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(parameters);
    if ( e->ApEINV.to_clear) return 0;

    ret = mess_multidirect_init(&(eqn->Amsolver));                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_multidirect_init);
    ret = mess_multidirect_create(eqn->A, NULL, parameters , eqn->Amsolver, NULL, eqn->E);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_multidirect_create);

    e->ApEINV.to_clear = 1;
    return 0;

}

static int ApEINV_apply(mess_equation e, mess_operation_t op, mess_double_cpx_t p, mess_int_t idx_p, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);

    ret = mess_multidirect_solvem(op,eqn->Amsolver, idx_p, in,out);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_multidirect_solvem);

    return 0;
}

static int ApEINV_clear(mess_equation e){
    __lyap *eqn = (__lyap *) e->aux;
    MSG_FNAME(__func__);
    mess_check_nullpointer(eqn);
    if ( e->ApEINV.to_clear == 0) return 0;
    if ( eqn -> Amsolver != NULL ) {
        mess_multidirect_clear(&(eqn->Amsolver));
    }
    e->ApEINV.to_clear = 0;
    return 0;
}

static int ApEX_generate(mess_equation e, mess_vector parameters){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    mess_check_nullpointer(parameters);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    if ( e->ApEX.to_clear) return 0;
    ret = mess_matrix_init(&(eqn->ApEX_tmp));   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    e->ApEX.to_clear = 1;
    return 0;
}

static int ApEX_apply(mess_equation e, mess_operation_t op, mess_double_cpx_t p, mess_int_t idx_p, mess_matrix in, mess_matrix out){
    MSG_FNAME(__func__);
    double alpha = 1;
    int ret = 0;
    mess_check_nullpointer(e);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);

    switch(op){
        case MESS_OP_NONE:
            ret = mess_matrix_multiply(MESS_OP_NONE, eqn->A, MESS_OP_NONE, in, eqn->ApEX_tmp);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mul);
            ret = mess_matrix_multiply(MESS_OP_NONE, eqn->E, MESS_OP_NONE, in, out);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mul);
            ret = mess_matrix_addc(1,eqn->ApEX_tmp, alpha*p, out);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
            break;
        case MESS_OP_TRANSPOSE:
        case MESS_OP_HERMITIAN:
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, eqn->A, MESS_OP_NONE, in, eqn->ApEX_tmp);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mul);
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, eqn->E, MESS_OP_NONE, in, out);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mul);
            ret = mess_matrix_addc(1,eqn->ApEX_tmp, alpha*conj(p), out);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);

            break;
        default:
            MSG_ERROR("Unknown Operation\n");
            return MESS_ERROR_ARGUMENTS;
    }
    return 0;
}

static int ApEX_clear(mess_equation e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    if ( e->ApEX.to_clear == 0) return 0;
    mess_matrix_clear(&(eqn->ApEX_tmp));
    e->ApEX.to_clear = 0;
    return 0;
}

static int clear(mess_equation e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    mess_free(eqn);
    eqn = NULL;
    e->aux = NULL;
    return 0;
}

static int parameter(mess_equation e, mess_options opt, mess_status stat){
    if(e->parent){
        //equation was stabilised by equation_stable
        //in this case we use the stabilized function handles to compute shift parameters
        return mess_lrcfadi_parameter(e->parent,opt,stat);
    }
    //equation was not stabilised by equation_stable
    return mess_lrcfadi_parameter(e,opt,stat);
}

/**
 * @brief Setup a mess_Equation object to solve a generalized Lyapunov Equation.
 * @param[in,out]  e    mess_equation object
 * @param[in]      opt  input @ref mess_options object
 * @param[in]      A    input system matrix of the equation
 * @param[in]      E    input mass matrix of the equation
 * @param[in]      B    input right hand side of the equation
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_equation_glyap function sets up a mess_equation object such that the Lyapunov
 * equation
 * \f[ op(A)Xop(E)^T+op(E) X op(A)^T +op(B)op(B)^T = 0 \f]
 * is represented. \n
 * Matrices need to be real. The operation is selected by the @ref mess_options.type component. \n
 *
 * If @ref mess_options.type attribute of @p opt is set to \ref MESS_OP_NONE the ADI process will solve the standard Lyapunov Equation,
 * i.e. \f$ op(A)=A \f$. \n
 * Otherwise if @ref mess_options.type is set to \ref MESS_OP_TRANSPOSE the transposed Lyapunov Equation with
 * \f$ op(A)=A^T \f$ will be solved. The resulting object can be used with the \ref mess_lrcfadi_adi solver.
 *
 * \sa mess_lrcfadi_adi
 * \sa mess_equation_lyap
 * \sa mess_equation_st
 *
 */
int mess_equation_glyap(mess_equation e, mess_options opt , mess_matrix A, mess_matrix E, mess_matrix B) {
    MSG_FNAME(__func__);
    __lyap *data ;
    int ret = 0;
    mess_operation_t op;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(E);
    mess_check_nullpointer(B);
    mess_check_nullpointer(e);
    mess_check_nullpointer(opt);
    mess_check_real(A);
    mess_check_real(B);
    mess_check_real(E);
    mess_check_square(A);
    mess_check_square(E);
    mess_check_same_size(A,E);


    op = opt->type;
    if (op == MESS_OP_HERMITIAN) op= MESS_OP_TRANSPOSE;

    if ( op == MESS_OP_NONE &&  A->rows != B->rows) {
        MSG_ERROR("A and B have the wrong number of rows. A->rows = " MESS_PRINTF_INT ", B->rows = " MESS_PRINTF_INT"\n", A->rows, B->rows);
        return MESS_ERROR_DIMENSION;
    }
    if ( op == MESS_OP_TRANSPOSE && A->rows != B->cols) {
        MSG_ERROR("A and B do not have fitting dimensions.  A->rows = " MESS_PRINTF_INT ", B->cols = " MESS_PRINTF_INT"\n", A->rows, B->cols);
        return MESS_ERROR_DIMENSION;
    }

    if ( op == MESS_OP_NONE ) {
        if ( B->cols >= B->rows) {
            MSG_WARN("Oversized right hand side factor for LRCF-ADI (op = %s, rows=" MESS_PRINTF_INT ", cols=" MESS_PRINTF_INT")\n ", mess_operation_t_str(op), B->rows, B->cols);
        }
    }
    if ( op == MESS_OP_TRANSPOSE ) {
        if ( B->rows >= B->cols) {
            MSG_WARN("Oversized right hand side factor for LRCF-ADI (op = %s, rows=" MESS_PRINTF_INT ", cols=" MESS_PRINTF_INT")\n ", mess_operation_t_str(op), B->rows, B->cols);
        }
    }

    /*-----------------------------------------------------------------------------
     * Set DEFAULT
     *-----------------------------------------------------------------------------*/
    if ( e->clear != NULL ) {
        e->clear(e);
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
    mess_try_alloc(data,__lyap *, sizeof(__lyap));
    data->A = A;
    data->E = E;
    data->Asolver = NULL;
    data->Amsolver = NULL;
    data->ApEX_tmp = NULL;


    e->dim = A->rows;
    e->aux = data;
    if ( op == MESS_OP_NONE) {
        e->clearRHS = 0;
        e->RHS = B;
    } else {
        ret = mess_matrix_init(&e->RHS);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        e->clearRHS = 1;
        ret = mess_matrix_ctranspose(B, e->RHS); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
    }
    e->eqn_type         = MESS_EQN_GLYAP;
    e->clear            = clear;
    e->AX.apply         = AX_apply;
    e->EX.apply         = EX_apply;
    e->AINV.generate    = AINV_generate;
    e->AINV.apply       = AINV_apply;
    e->AINV.clear       = AINV_clear;
    e->EINV.generate    = EINV_generate;
    e->EINV.apply       = EINV_apply;
    e->EINV.clear       = EINV_clear;
    e->ApEINV.generate  = ApEINV_generate;
    e->ApEINV.apply     = ApEINV_apply;
    e->ApEINV.clear     = ApEINV_clear;
    e->ApEX.generate    = ApEX_generate;
    e->ApEX.apply       = ApEX_apply;
    e->ApEX.clear       = ApEX_clear;
    e->parameter        = parameter;

    return 0;
}

/**
 * @brief Setup a mess_equation object to solve a generalized Riccati Equation.
 * @param[in,out]  e    mess_equation object
 * @param[in]      opt input    mess_options object
 * @param[in]      A    input   system matrix of the equation
 * @param[in]      E    input mass matrix of the equation
 * @param[in]      B    input factor of the matrix in the quadratic term
 * @param[in]      C    input right hand side
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_equation_griccati function sets up a mess_equation object such that the Riccati
 * equation
 * \f[ AXE^T+EXA^T-EXC^TCXE^T+BB^T = 0 \f]
 * or the transposed one
 * \f[ A^TXE+E^TXA-E^TXBB^TXE+C^TC = 0 \f]
 * is represented.\n
 * Depending on the @ref mess_options.type component in the \ref mess_options object the first or the second one is solved.\n
 * The first one corresponds to th that is equal to \ref MESS_OP_NONE and the second one
 * corresponds to  the @ref mess_options.type  that is equal to \ref MESS_OP_TRANSPOSE.
 *
 * \sa mess_lrcfadi_adi
 * \sa mess_equation_riccati
 * \sa mess_equation_st
 * \sa mess_lrcfadi_nm
 *
 */
int mess_equation_griccati(mess_equation e, mess_options opt , mess_matrix A, mess_matrix E, mess_matrix B, mess_matrix C) {
    MSG_FNAME(__func__);
    __lyap *data ;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(E);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);
    mess_check_nullpointer(e);
    mess_check_nullpointer(opt);
    mess_check_real(A);
    mess_check_real(B);
    mess_check_real(E);
    mess_check_square(A);
    mess_check_square(E);
    mess_check_same_size(A,E);
    if ( A->rows != B->rows) {
        MSG_ERROR("A and B have the wrong number of rows. A->rows = " MESS_PRINTF_INT ", B->rows = " MESS_PRINTF_INT"\n", A->rows, B->rows);
        return MESS_ERROR_DIMENSION;
    }
    if ( C->cols != A->rows) {
        MSG_ERROR("A and C have the wrong number of rows/columns. A->rows = " MESS_PRINTF_INT ", C->cols = " MESS_PRINTF_INT"\n", A->rows, C->cols);
        return MESS_ERROR_DIMENSION;
    }

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
    mess_try_alloc(data,__lyap *, sizeof(__lyap));
    data->A = A;
    data->E = E;
    data->Asolver = NULL;
    data->Amsolver = NULL;
    data->ApEX_tmp = NULL;


    e->dim = A->rows;
    e->aux = data;
    e->B                = B;
    e->C                = C;
    e->RHS              = NULL;
    e->eqn_type         = MESS_EQN_GRICCATI;
    e->clear            = clear;
    e->AX.apply         = AX_apply;
    e->EX.apply         = EX_apply;
    e->AINV.generate    = AINV_generate;
    e->AINV.apply       = AINV_apply;
    e->AINV.clear       = AINV_clear;
    e->EINV.generate    = EINV_generate;
    e->EINV.apply       = EINV_apply;
    e->EINV.clear       = EINV_clear;
    e->ApEINV.generate  = ApEINV_generate;
    e->ApEINV.apply     = ApEINV_apply;
    e->ApEINV.clear     = ApEINV_clear;
    e->ApEX.generate    = ApEX_generate ;
    e->ApEX.apply       = ApEX_apply;
    e->ApEX.clear       = ApEX_clear;
    e->parameter        = parameter;

    return 0;
}


