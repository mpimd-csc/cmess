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
 * @file lib/lrcf_adi/equation_lyap.c
 * @brief Build an equation object (Lyapunov or Riccati Equation).
 * @author @koehlerm
 *
 * Generate the @ref mess_equation object for the standard Lyapunov and Riccati Equations.
 *
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
 * @brief Contains matrices, solvers and additional data for a Lyapunov Equation.
 *
 * The @ref __lyap structure contains all matrices, solvers and additional data for a Lyapunov Equation
 * \f[AX +AF^T = -BB^T\f]
 * or
 * \f[A^TX +XA = -B^TB.\f]
 * @attention Internal use only.
 *
 */
typedef struct {
    mess_matrix A;                  /**< System matrix \f$A\f$ of Lyapunov Equation.*/
    mess_matrix B;                  /**< System matrix \f$B\f$ of Lyapunov Equation.*/
    mess_direct Asolver;            /**< Direct solver for system matrix \f$A.\f$ */
    mess_multidirect Amsolver;      /**< Multisolver for shifted  matrix \f$ A +pI.\f$ */
    mess_matrix ApEX_tmp;           /**< Temporary shifted \f$ A+pI.\f$ */
} __lyap;

/**
 * @internal
 * @brief Compute matrix-matrix multitplication with system matrix \f$A\f$.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to \f$A\f$
 * @param[in] X         input matrix \f$X\f$
 * @param[out] Y        output matrix \f$Y\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref __lyap_A_multiplym function computes the matrix-matrix product, i.e. \n
 * \f[ Y = op(A) X. \f]
 *
 * @see mess_matrix_multiply
 * @attention Internal use only.
 */
static int __lyap_A_multiplym(mess_equation e, mess_operation_t op, mess_matrix X, mess_matrix Y) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    switch(op){
        case MESS_OP_NONE:
            return mess_matrix_multiply(MESS_OP_NONE, eqn->A, MESS_OP_NONE, X, Y);
        case MESS_OP_TRANSPOSE:
            return mess_matrix_multiply(MESS_OP_HERMITIAN, eqn->A, MESS_OP_NONE, X, Y);
        case MESS_OP_HERMITIAN:
            return mess_matrix_multiply(MESS_OP_HERMITIAN, eqn->A, MESS_OP_NONE, X, Y);
        default:
            return MESS_ERROR_ARGUMENTS;
    }
    return 0;
}

/**
 * @internal
 * @brief Generate a direct solver for system matrix \f$A\f$.
 * @param[in, out] e         input @ref mess_equation structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref __lyap_AINV_generate function generates a @ref mess_direct_lu solver for system matrix \f$A\f$. \n
 * After function call the @ref __lyap.Asolver field holds a generated @ref mess_direct_lu solver.
 * If @ref __lyap.Asolver is already generated \f$0\f$ is return immediately.
 *
 * @see mess_direct_init
 * @see mess_direct_lu
 * @attention Internal use only.
 */
static int __lyap_AINV_generate(mess_equation e){
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

/**
 * @internal
 * @brief Solve system \f$op(A)X=Y\f$.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to \f$A\f$
 * @param[in] Y         input matrix right-hand side \f$Y\f$
 * @param[out] X        output matrix solution \f$X\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref __lyap_AINV_apply function solves \f$op(A)X=Y.\f$
 *
 * @see mess_direct_solvem
 * @attention Internal use only.
 */
static int __lyap_AINV_apply(mess_equation e, mess_operation_t op, mess_matrix Y, mess_matrix X){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    return mess_direct_solvem(op, eqn->Asolver, Y, X);

}

/**
 * @internal
 * @brief Clear a direct solver for system matrix \f$A\f$.
 * @param[in, out] e         input @ref mess_equation structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref __lyap_AINV_clear function clears a generated solver in  @ref __lyap.Asolver.
 *
 * @see mess_direct_clear
 * @attention Internal use only.
 */
static int __lyap_AINV_clear(mess_equation e){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    if (eqn->Asolver) mess_direct_clear(&(eqn->Asolver));
    eqn->Asolver = NULL;
    e->AINV.to_clear = 0;
    return 0;

}

/**
 * @internal
 * @brief
 * @param[in, out] e        input @ref mess_equation structure
 * @param[in] parameters    input shift parameter vector
 * @return zero on succes or a non-zero error value otherwise
 *
 * @attention Internal use only.
 * */
static int __lyap_ApE_generate(mess_equation e, mess_vector parameters) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(parameters);
    if ( e->ApEINV.to_clear) return 0;
    ret = mess_multidirect_init(&(eqn->Amsolver));                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_multidirect_init);
    ret = mess_multidirect_create(eqn->A, NULL, parameters , eqn->Amsolver, NULL, NULL);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_multidirect_create);

    e->ApEINV.to_clear = 1;
    return 0;

}

/**
 * @internal
 * @brief
 * @param[in, out] e        input @ref mess_equation structure
 * @param[in]   op          input operation type
 * @param[in]   p           input shift parameter
 * @param[in]   idx_p       input index of parameter for shift
 * @param[in]   in          input matrix
 * @param[out]  out         output matrix
 * @return zero on succes or a non-zero error value otherwise
 *
 * @attention Internal use only.
 */
static int __lyap_ApE_solve(mess_equation e, mess_operation_t op, mess_double_cpx_t p, mess_int_t idx_p, mess_matrix in, mess_matrix out) {
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

/**
 * @internal
 * @brief
 * @param[in, out] e         input @ref mess_equation structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * @attention Internal use only.
 * */
static int __lyap_ApE_clear(mess_equation e){
    __lyap *eqn = (__lyap *) e->aux;
    MSG_FNAME(__func__);
    mess_check_nullpointer(eqn);
    if ( eqn -> Amsolver != NULL ) {
        mess_multidirect_clear(&(eqn->Amsolver));
    }
    e->ApEINV.to_clear = 0;
    return 0;
}

/**
 * @internal
 * @brief
 * @param[in,out] e         input @ref mess_equation structure
 * @param[in] parameters    input vector of shift parameters
 * @return zero on succes or a non-zero error value otherwise
 *
 * @attention Internal use only.
 */
static int __lyap_ApEX_generate(mess_equation e, mess_vector parameters){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    mess_check_nullpointer(parameters);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    ret = mess_matrix_init(&(eqn->ApEX_tmp));   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    e->ApEX.to_clear = 1;
    return 0;
}

/**
 * @internal
 * @brief
 * @param[in, out] e        input @ref mess_equation structure
 * @param[in]   op          input operation type
 * @param[in]   p           input shift parameter
 * @param[in]   idx_p       input index of parameter for shift
 * @param[in]   in          input matrix
 * @param[out]  out         output matrix
 * @return zero on succes or a non-zero error value otherwise
 *
 * @attention Internal use only.
 */
static int __lyap_ApEX_apply(mess_equation e, mess_operation_t op,mess_double_cpx_t p, mess_int_t idx_p, mess_matrix in, mess_matrix out){
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
            ret = mess_matrix_copy(in, out);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            ret = mess_matrix_addc(1,eqn->ApEX_tmp, alpha*p, out);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
            break;
        case MESS_OP_TRANSPOSE:
        case MESS_OP_HERMITIAN:
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, eqn->A, MESS_OP_NONE, in, eqn->ApEX_tmp);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mul);
            ret = mess_matrix_copy(in, out);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            ret = mess_matrix_addc(1,eqn->ApEX_tmp, alpha*conj(p), out);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);

            break;
        default:
            MSG_ERROR("Unknown Operation\n");
            return MESS_ERROR_ARGUMENTS;
    }
    return 0;

}

/**
 * @internal
 * @brief
 * @param[in, out] e         input @ref mess_equation structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * @attention Internal use only.
 * */
static int __lyap_ApEX_clear(mess_equation e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    mess_matrix_clear(&(eqn->ApEX_tmp));
    e->ApEX.to_clear = 0;
    return 0;
}

/**
 * @internal
 * @brief
 * @param[in, out] e         input @ref mess_equation structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * @attention Internal use only.
 * */
static int __lyap_clear(mess_equation e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __lyap *eqn = (__lyap *) e->aux;
    mess_check_nullpointer(eqn);
    mess_free(eqn);
    eqn = NULL;
    e->aux = NULL;
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
    if(e->parent){
        //equation was stabilised by equation_stable
        //in this case we use the stabilized function handles to compute shift parameters
        return mess_lrcfadi_parameter(e->parent,opt,stat);
    }
    //equation was not stabilised by equation_stable
    return mess_lrcfadi_parameter(e,opt,stat);
}

/**
 * @brief Setup a mess_Equation object to solve a standard Lyapunov Equation.
 * @param[in,out]  e     @ref mess_equation object
 * @param[in]      opt   input @ref mess_options object
 * @param[in]      A     input system matrix of equation
 * @param[in]      B     input right hand side ofequation
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_equation_std_lyap function sets up a mess_equation object such that the Lyapunov
 * equation
 * \f[ op(A) X + X op(A)^T + op(B) op(B)^T = 0  \f]
 * is represented. Matrices need to be real. \n
 * The operation is selected by the @ref mess_options.type component of the  structure. \n
 * If the @ref mess_options.type is set to \ref MESS_OP_NONE  the ADI process will solve the standard Lyapunov Equation,
 * i.e. \f$ op(A)=A \f$. \n
 * Otherwise if @ref mess_options.type is set to \ref MESS_OP_TRANSPOSE the transposed Lyapunov Equation will be
 * solved, i.e.  \f$ op(A)=A^T \f$. \n
 * The resulting object can be used with the  \ref mess_lrcfadi_adi solver.
 *
 * @see mess_lrcfadi_adi
 * @see mess_equation_glyap
 * @see mess_equation_st
 * @see mess_lrcfadi_nm
 *
 */
int mess_equation_std_lyap(mess_equation e, mess_options opt, mess_matrix A, mess_matrix B) {
    MSG_FNAME(__func__);
    __lyap *data ;
    mess_operation_t op;
    int ret;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_real(A);
    mess_check_real(B);
    mess_check_square(A);
    mess_check_nullpointer(e);

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
        if ( B->cols > B->rows) {
            MSG_WARN("Oversized right hand side factor for LRCF-ADI (op = %s, rows=%d, cols=%d)\n ", mess_operation_t_str(op), (int) B->rows, (int) B->cols);
        }
    }
    if ( op == MESS_OP_TRANSPOSE ) {
        if ( B->rows > B->cols) {
            MSG_WARN("Oversized right hand side factor for LRCF-ADI (op = %s, rows=%d, cols=%d)\n ", mess_operation_t_str(op), (int) B->rows, (int) B->cols);
        }
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
    data->Asolver = NULL;
    data->Amsolver = NULL;
    data->ApEX_tmp = NULL;

    e->eqn_type = MESS_EQN_LYAP;
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

    e->clear            = __lyap_clear;
    e->AX.apply         = __lyap_A_multiplym;
    e->AINV.generate    = __lyap_AINV_generate;
    e->AINV.apply       = __lyap_AINV_apply;
    e->AINV.clear       = __lyap_AINV_clear;
    e->ApEINV.generate  = __lyap_ApE_generate;
    e->ApEINV.clear     = __lyap_ApE_clear;
    e->ApEINV.apply     = __lyap_ApE_solve;
    e->ApEX.generate    = __lyap_ApEX_generate;
    e->ApEX.apply       = __lyap_ApEX_apply;
    e->ApEX.clear       = __lyap_ApEX_clear;
    e->parameter        = parameter;

    return 0;
}


/**
 * @brief Setup a mess_equation object to solve a standard Riccati Equation.
 * @param[in,out]  e    mess_equation object
 * @param[in]      opt   input mess_options object
 * @param[in]      A     input system matrix of equation
 * @param[in]      B     input factor of the matrix in quadratic term
 * @param[in]      C     input right hand side
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_equation_std_riccati function sets up a mess_equation object such that the Riccati
 * equation
 * \f[ A X + X A^T - X C^T C X + B B^T = 0 \f]
 * or the transposed one
 * \f[ A^T X + X A - X B B^T X + C^T C = 0 \f]
 * is represented. \n
 * Depending on the @ref mess_options.type component the first or the second equation is
 * solved. \n
 * The first one corresponds to a @ref mess_options.type that is equal to \ref MESS_OP_NONE and the second one
 * corresponds to a @ref mess_options.type that is equal to \ref MESS_OP_TRANSPOSE.
 *
 * @see mess_lrcfadi_adi
 * @see mess_equation_griccati
 * @see mess_equation_st
 * @see mess_lrcfadi_nm
 *
 */
int mess_equation_std_riccati(mess_equation e, mess_options opt, mess_matrix A, mess_matrix B, mess_matrix C) {
    MSG_FNAME(__func__);
    __lyap *data ;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_real(A);
    mess_check_real(B);
    mess_check_square(A);
    if ( A->rows != B->rows) {
        MSG_ERROR("A and B have the wrong number of rows. A->rows = " MESS_PRINTF_INT ", B->rows = " MESS_PRINTF_INT"\n", A->rows, B->rows);
        return MESS_ERROR_DIMENSION;
    }
    mess_check_nullpointer(e);
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
    data->Asolver = NULL;
    data->Amsolver = NULL;
    data->ApEX_tmp = NULL;

    e->B = B;
    e->C = C;

    e->eqn_type = MESS_EQN_RICCATI;
    e->dim = A->rows;
    e->aux = data;
    e->RHS = NULL;
    e->clear            = __lyap_clear;
    e->AX.apply         = __lyap_A_multiplym;
    e->AINV.generate    = __lyap_AINV_generate;
    e->AINV.apply       = __lyap_AINV_apply;
    e->AINV.clear       = __lyap_AINV_clear;
    e->ApEINV.generate  = __lyap_ApE_generate;
    e->ApEINV.clear     = __lyap_ApE_clear;
    e->ApEINV.apply     = __lyap_ApE_solve;
    e->ApEX.generate    = __lyap_ApEX_generate;
    e->ApEX.apply       = __lyap_ApEX_apply;
    e->ApEX.clear       = __lyap_ApEX_clear;
    e->parameter        = parameter;

    return 0;
}




/**
 * @brief Setup a mess_Equation object to solve a matrix-valued Lyapunov Equation.
 * @param[in,out]  e    mess_equation object
 * @param[in]      opt   input mess_options object
 * @param[in]      A     input system matrix of equation
 * @param[in]      E     input mass matrix of the equation, @c NULL if it is the identity
 * @param[in]      B     input right hand side ofequation
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_equation_lyap function sets up a mess_equation object such that the Lyapunov
 * equation
 * \f[ op(A) X op(E)^T + op(E) X op(A)^T + op(B) op(B)^T = 0  \f]
 * is represented. Matrices need to be real. \n
 * The operation is selected by the @ref mess_options.type component. \n
 * i.e. \f$ op(A)=A \f$. \n
 * Otherwise if the @ref mess_options.type is set to \ref MESS_OP_TRANSPOSE the transposed Lyapunov Equation will be
 * solved, i.e.  \f$ op(A)=A^T \f$. \n
 * The resulting object can be used with the  \ref mess_lrcfadi_adi solver.
 *
 * Depending on if the E parameter is @c NULL or not the \ref mess_equation_lyap or \ref mess_equation_glyap function
 * is called to setup the equation object.
 *
 * @see mess_lrcfadi_adi
 * @see mess_equation_lyap
 * @see mess_equation_glyap
 * @see mess_equation_st
 * @see mess_lrcfadi_nm
 *
 */
int mess_equation_lyap(mess_equation e, mess_options opt, mess_matrix A, mess_matrix E, mess_matrix B)
{
    if ( E == NULL ) {
        return mess_equation_std_lyap(e, opt, A, B);
    } else {
        return mess_equation_glyap(e, opt, A, E, B);
    }
}

/**
 * @brief Setup a mess_equation object to solve a matrix-valued Riccati Equation.
 * @param[in,out]  e    mess_equation object
 * @param[in]      opt   input mess_options object
 * @param[in]      A     input system matrix of equation
 * @param[in]      E     input mass matrix of the equation, @c NULL if it is the identity
 * @param[in]      B     input factor of the matrix in quadratic term
 * @param[in]      C     input right hand side
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_equation_riccati function sets up a mess_equation object such that the Riccati
 * equation
 * \f[ A X E^T + E X A^T - E X C^T C X E^T + B B^T = 0 \f]
 * or the transposed one
 * \f[ A^T X E+ E^T X A - E^T X B B^T X E + C^T C = 0 \f]
 * is represented. \n
 * Depending on the @ref mess_options.type component the first or the second equation is
 * solved. \n
 * The first one corresponds to the @ref mess_options.type that is equal to \ref MESS_OP_NONE and the second one
 * corresponds to the @ref mess_options.type that is equal to \ref MESS_OP_TRANSPOSE.
 *
 * Depending on if the E parameter is @c NULL or not the \ref mess_equation_riccati  or \ref mess_equation_griccati function
 * is called to setup the equation object.
 *
 * @see mess_lrcfadi_adi
 * @see mess_equation_riccati
 * @see mess_equation_griccati
 * @see mess_equation_st
 * @see mess_lrcfadi_nm
 *
 */
int mess_equation_riccati(mess_equation e, mess_options opt, mess_matrix A, mess_matrix E, mess_matrix B, mess_matrix C)
{
    if ( E == NULL ) {
        return mess_equation_std_riccati(e, opt, A, B, C);
    } else {
        return mess_equation_griccati(e, opt, A, E, B, C);
    }
}

