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
 * @file lib/lrcf_adi/equation_glyap_so2.c
 * @brief Build an equation object (generalized Lyapunov or Riccati Equation) out of a second order system.
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
 * @brief Contains matrices, solvers and additional data for Second Order System function handles.
 *
 * The @ref __glyap_so2 structure contains all matrices, solvers and additional data for Second Order System function handles.
 * @attention Internal use only.
 *
 */
typedef struct {
    mess_matrix M;              /**< System matrix \f$ M \f$.*/
    mess_matrix D;              /**< System matrix \f$ D \f$.*/
    mess_matrix K;              /**< System matrix \f$ K \f$.*/

    mess_direct Ksolver;            /**< Direct solver fo \f$ K \f$.*/
    mess_direct Msolver;            /**< Direct solver fo \f$ M \f$.*/


    mess_direct* shiftsolver;       /**< Array of direct solvers for shifted system \f$ K-p*D+p^2*M \f$ .*/
    mess_int_t shifts;          /**< Number of direct solvers in @ref __glyap_so2.shiftsolver. */

    double upperbound;          /**< upperbound for absolute value of shift parameter.*/
    double lowerbound;          /**< lowerbound for absolute value of shift parameter.*/

} __glyap_so2;


/**
 * @internal
 * @brief Multiply with operator \f$ \begin{bmatrix} 0 & I \\ -K & -D \end{bmatrix} \f$ to a @ref mess_matrix.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref AX_apply function computes \f$ out = \begin{bmatrix} 0 & I \\ -K & -D \end{bmatrix}  in\f$.
 *
 * @attention Internal use only.
 */
static int AX_apply(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_so2 *eqn = (__glyap_so2 *) e->aux;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_operation_type(op);
    int ret =0;
    mess_int_t n =eqn->M->rows;

    //MSG_PRINT("AX_apply in->cols=%d \t in->rows=%d \t in->data_type=%d\n",in->cols,in->rows, in->data_type);
    /*-----------------------------------------------------------------------------
     *  perform [x2;-K*x1-D*x2]->out
     *-----------------------------------------------------------------------------*/
    mess_matrix f1, f2, x1, x2;
    MESS_INIT_MATRICES(&x1,&x2,&f1,&f2);

    //take upper part of in
    ret = mess_matrix_rowsub(in, 0,n-1, x1);                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);

    //take lower part of in
    ret = mess_matrix_rowsub(in, n, 2*n-1, x2);                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);

    if(op==MESS_OP_NONE){
        ret = mess_matrix_multiply(op, eqn->D, MESS_OP_NONE, x2, f1);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mul);
        ret = mess_matrix_multiply(op, eqn->K, MESS_OP_NONE, x1, f2);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mul);

        ret = mess_matrix_add(-1.0,f1,-1.0,f2);                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mul);
        ret = mess_matrix_cat(x2,NULL,f2,NULL,MESS_DENSE,out);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);

    }else{
        ret = mess_matrix_multiply(op, eqn->K, MESS_OP_NONE, x2, f2);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        ret = mess_matrix_scale(-1,f2);                                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
        ret = mess_matrix_multiply(op, eqn->D, MESS_OP_NONE, x2, f1);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        ret = mess_matrix_add(-1,f1,1,x1);                                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
        ret = mess_matrix_cat(f2,NULL,x1,NULL,MESS_DENSE,out);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);

    }

    /*-----------------------------------------------------------------------------
     *  clear additional memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&x1,&x2,&f1,&f2);

    return 0;
}

/**
 * @internal
 * @brief Multiply with operator \f$ \begin{bmatrix} I & 0 \\ 0 & M \end{bmatrix} \f$ to a @ref mess_matrix.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref EX_apply function computes \f$ out = \begin{bmatrix} I & 0 \\ 0 & M\end{bmatrix}  in\f$.
 *
 * @attention Internal use only.
 */
static int EX_apply(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_so2 *eqn = (__glyap_so2 *) e->aux;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_operation_type(op);
    int ret =0;
    mess_int_t n = eqn->M->rows;

    //MSG_PRINT("EX_apply in->cols=%d \t in->rows=%d \t in->data_type=%d\n",in->cols,in->rows, in->data_type);
    /*-----------------------------------------------------------------------------
     *  perform [x1;M*x2]->out
     *-----------------------------------------------------------------------------*/
    mess_matrix x1, x2, f2;

    ret = mess_matrix_init(&x1);                                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&x2);                                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&f2);                                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);

    //take upper part of in
    ret = mess_matrix_rowsub(in, 0,n-1, x1);                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);

    //take lower part of in
    ret = mess_matrix_rowsub(in, n, 2*n-1, x2);                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);

    ret = mess_matrix_multiply(op, eqn->M, MESS_OP_NONE, x2, f2);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mul);
    ret = mess_matrix_cat(x1,NULL,f2,NULL,MESS_DENSE,out);                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);


    /*-----------------------------------------------------------------------------
     *  clear additional memory
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_clear(&x1);                                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret = mess_matrix_clear(&x2);                                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret = mess_matrix_clear(&f2);                                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);

    return 0;
}

/**
 * @internal
 * @brief Solve with operator \f$ \begin{bmatrix} 0 & I \\ -K & -D \end{bmatrix} \f$ to a @ref mess_matrix.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref AINV_apply function solves the system \f$ \begin{bmatrix} 0 & I \\ -K & -D \end{bmatrix} out=in \f$
 *  via the equivalent formulation
 * \f[
 *  \begin{aligned}
 *      out(n+1:end,:)  &=  in(1:n,:) \\
 *      out(1:n,:)  &= -K^{-1}(in(n+1:end,:) +D out(n+1:end,:) )
 *  \end{aligned}.
 * \f]
 *
 * @attention Internal use only.
 */
static int AINV_apply(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out){

    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_so2 *eqn = (__glyap_so2 *) e->aux;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_operation_type(op);

    int ret =0;
    mess_int_t n = eqn->M->rows;

    //MSG_PRINT("AINV_apply in->cols=%d \t in->rows=%d \t in->data_type=%d\n",in->cols,in->rows, in->data_type);
    /*-----------------------------------------------------------------------------
     *  perform [ -K\(x2+ D*x1) ; x1 ]->out
     *          [ x2 +D^T(-K^T\x1) ;-K^T\x1]->out (MESS_OP_TRANSPOSE)
     *-----------------------------------------------------------------------------*/
    mess_matrix x1, x2, f1,f2;

    ret = mess_matrix_init(&x1);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&x2);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&f1);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&f2);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);

    ret = mess_matrix_rowsub(in,0,n-1,x1);                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
    ret = mess_matrix_rowsub(in,n,2*n-1,x2);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);

    if(op==MESS_OP_NONE){
        ret = mess_matrix_multiply(op, eqn->D, MESS_OP_NONE, x1, f2);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mul);
        ret = mess_matrix_add(1.0, x2, 1.0, f2);                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
        ret = mess_direct_solvem(op, eqn->Ksolver,f2, f1 );                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solve);
        ret = mess_matrix_scale(-1.0,f1);                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
        ret = mess_matrix_cat(f1, NULL, x1, NULL, MESS_DENSE, out);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);
    }else{
        ret = mess_direct_solvem(op,eqn->Ksolver,x1,f2);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
        ret = mess_matrix_scale(-1,f2);                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
        ret = mess_matrix_multiply(op, eqn->D, MESS_OP_NONE, f2, f1);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        ret = mess_matrix_add(1,f1,1,x2);                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
        ret = mess_matrix_cat(x2,NULL,f2,NULL,MESS_DENSE,out);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);
    }

    /*-----------------------------------------------------------------------------
     *  clear additional memory
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_clear(&x1);                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_clear);
    ret = mess_matrix_clear(&x2);                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_clear);
    ret = mess_matrix_clear(&f1);                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_clear);
    ret = mess_matrix_clear(&f2);                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_clear);

    return 0;


}

/**
 * @internal
 * @brief Solve with \f$ \begin{bmatrix} I & 0 \\ 0 & M \end{bmatrix} \f$.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] op        input operation applied to block matrix
 * @param[in] in        input matrix \f$in\f$
 * @param[out] out      output matrix \f$out\f$
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref EINV_apply function solves the system \f$ \begin{bmatrix} I & 0 \\ 0 & M\end{bmatrix}out=  in\f$ via
 *
 * \f[
 *  \begin{aligned}
 *      out(1:n,:) &= in(1:n,:)     \\
 *      out(n+1:end,:) &= M^{-1}(in(n+1:end,:))
 *  \end{aligned}
 *\f]
 *
 * @attention Internal use only.
 */
static int EINV_apply(mess_equation e, mess_operation_t op, mess_matrix in , mess_matrix out){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_so2 *eqn = (__glyap_so2 *) e->aux;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_operation_type(op);
    int ret = 0;
    mess_int_t n = eqn->M->rows, n2 = 2*n;

    //MSG_PRINT("EINV_apply in->cols=%d \t in->rows=%d \t in->data_type=%d\n",in->cols,in->rows, in->data_type);
    /*-----------------------------------------------------------------------------
     *  perform [x1;M\x2]->out
     *-----------------------------------------------------------------------------*/
    mess_matrix temp1, temp2, temp3;

    ret = mess_matrix_init(&temp1);                                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&temp2);                                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&temp3);                                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);

    ret = mess_matrix_rowsub(in, 0, n-1, temp1);                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
    ret = mess_matrix_rowsub(in, n, n2-1, temp2);                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
    ret = mess_direct_solvem(op,eqn->Msolver, temp2, temp3);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
    ret = mess_matrix_cat(temp1, NULL, temp3, NULL, MESS_DENSE, out);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);

    /*-----------------------------------------------------------------------------
     *  clear additional memory
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_clear(&temp1);                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret = mess_matrix_clear(&temp2);                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret = mess_matrix_clear(&temp3);                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);

    return 0;
}

/**
 * @internal
 * @brief Generate a shifted solver for the system matrix \f$ \begin{bmatrix} -pI & -I \\ -K & -D + pM \end{bmatrix}\f$.
 * @param[in] e             input @ref mess_equation structure
 * @param[in] parameters    input @ref mess_vector of shift parameters
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref ApEINV_generate builds a several direct solvers for the matrices \n
 * \f[
 *      \begin{bmatrix} -pI & -I \\ -K & -D + pM \end{bmatrix}.
 * \f]
 * It uses for each shift parameter a @ref mess_direct instance.
 *
 * @attention Internal use only.
 */
static int ApEINV_generate(mess_equation e, mess_vector parameters) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_so2 *eqn = (__glyap_so2 *) e->aux;
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(parameters);
    e->ApEINV.to_clear = 1;
    //return 0;
    int ret=0;

    if(!eqn->Ksolver){
        ret = mess_direct_init(&(eqn->Ksolver));                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
        ret = mess_direct_lu(eqn->K, eqn->Ksolver);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);
    }
    if(!eqn->shiftsolver){
        eqn->shifts = parameters->dim;
        mess_int_t i = 0;
        mess_matrix *temps;
        mess_try_alloc(eqn->shiftsolver,mess_direct*,sizeof(mess_direct)*parameters->dim);
        mess_try_alloc(temps,mess_matrix*,sizeof(mess_matrix)*parameters->dim);
        if(MESS_IS_REAL(parameters)){
            //#warning parallelize code
            for(i=0;i<parameters->dim;++i){
                ret = mess_matrix_init(&temps[i]);                                                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
                ret = mess_matrix_copy(eqn->M,temps[i]);                                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
                ret = mess_matrix_add(-1*parameters->values[i],eqn->D,pow(parameters->values[i],2),temps[i]);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
                ret = mess_matrix_add(1.0,eqn->K,1.0,temps[i]);                                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
                ret = mess_direct_init(&(eqn->shiftsolver[i]));                                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
                ret = mess_direct_lu(temps[i],eqn->shiftsolver[i]);                                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);
                ret = mess_matrix_clear(&temps[i]);                                                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
            }
        }else{
            //#warning parallelize code
            for(i=0;i<parameters->dim;++i){
                ret = mess_matrix_init(&temps[i]);                                                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
                ret = mess_matrix_copy(eqn->M,temps[i]);                                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
                ret = mess_matrix_addc(-1*parameters->values_cpx[i],eqn->D,cpow(parameters->values_cpx[i],2),temps[i]); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
                ret = mess_matrix_addc(1.0,eqn->K,1.0,temps[i]);                                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
                ret = mess_direct_init(&(eqn->shiftsolver[i]));                                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
                ret = mess_direct_lu(temps[i],eqn->shiftsolver[i]);                                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);
                ret = mess_matrix_clear(&temps[i]);                                                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
            }
        }
        mess_free(temps);
    }
    return ret;

    return 0;
}

/**
 * @internal
 * @brief Solves shifted linear block system matrix \f$ \begin{bmatrix} -pI &  I \\ -K & -D + pM \end{bmatrix}\f$.
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
 *      op\left(\begin{bmatrix} -pI & I \\ -K & -D + pM \end{bmatrix}\right) out = in
 * \f]
 * via
 *
 * \f[
 *  \begin{aligned}
 *       out(n+1:end,:) &= in(1:n,:) -p out(1:n,:)  \\
 *      (-K+pD-p^2M)out(1:n,:) &=  in(n+1:end,:) +D in(1:n,:)- p M in(1:n,:)
 *  \end{aligned}
 * \f]
 *
 * @attention Internal use only.
 */
static int ApEINV_apply (mess_equation e, mess_operation_t op, mess_double_cpx_t p, mess_int_t idx_p, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    __glyap_so2 *eqn = (__glyap_so2 *) e->aux;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_operation_type(op);
    mess_int_t n = eqn->M->cols;

    //MSG_PRINT("ApEINV_apply in->cols=%d \t in->rows=%d \t in->data_type=%d\n",in->cols,in->rows, in->data_type);
    /*-----------------------------------------------------------------------------
     *  perform [x1-p*( (K-pD+p^2M)\(p*x2+Kx1) );(K-pD+p^2M)\(p*x2+Kx1)]->out (MESS_OP_NONE)
     *          [1/p*x1+1/p*K^T( (K^T-pD^T+p^2M^T)\(p*x2-x1) );(K^T-pD^T+p^M^T)\(p*x2-x1)]->out (MESS_OP_TRANSPOSE)
     *-----------------------------------------------------------------------------*/
    mess_matrix x1,x2,f1,f2,rhs;
    ret = mess_matrix_init(&x1);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&x2);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&f1);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&f2);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&rhs);                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);

    //get submatrices of rhs
    ret = mess_matrix_rowsub(in,0,n-1,x1);                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);
    ret = mess_matrix_rowsub(in,n,2*n-1,x2);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rowsub);

    //compute rhs
    if(op==MESS_OP_NONE){
        ret = mess_matrix_multiply(op, eqn->K, MESS_OP_NONE, x1, rhs);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        ret = mess_matrix_addc(p,x2,1,rhs);                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_addc);
    }else{
        mess_matrix_copy(x1,rhs);                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        mess_matrix_addc(op==MESS_OP_HERMITIAN?conj(p):p,x2,-1,rhs);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_addc);
    }

    //build shiftsolver and solve
    ret = mess_direct_solvem(op,eqn->shiftsolver[idx_p],rhs,f2);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_solvem);

    //compute f1
    if(op==MESS_OP_NONE){
        ret = mess_matrix_addc(-1/p,f2,1/p,x1);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_addc);
    }else{
        ret = mess_matrix_multiply(op, eqn->K, MESS_OP_NONE, f2, rhs);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        ret = mess_matrix_addc(op==MESS_OP_HERMITIAN?1/conj(p):1/p,rhs,op==MESS_OP_HERMITIAN?1/conj(p):1/p,x1); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_addc);
    }
    //build result
    ret = mess_matrix_cat(x1,NULL,f2,NULL,MESS_DENSE,out);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);


    /*-----------------------------------------------------------------------------
     *  clear additional memory
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_clear(&f1);                                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret = mess_matrix_clear(&f2);                                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret = mess_matrix_clear(&x1);                                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret = mess_matrix_clear(&x2);                                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret = mess_matrix_clear(&rhs);                                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);

    return 0;
}

/**
 * @internal
 * @brief Clear a array of direct solvers.
 * @param[in] e         input @ref mess_equation structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref ApEINV_clear clear generated solvers by using @ref mess_direct_clear the fields of @ref __glyap_so2.shiftsolver are cleared.
 *
 * @attention Internal use only.
 *
 */
static int ApEINV_clear(mess_equation e){
    __glyap_so2 *eqn = (__glyap_so2 *) e->aux;
    MSG_FNAME(__func__);
    mess_check_nullpointer(eqn);
    int ret =0;
    if ( e->ApEINV.to_clear == 0) return 0;

    if(eqn->shiftsolver){
        mess_int_t i;
        for(i=0;i<eqn->shifts;++i){
            ret = mess_direct_clear(&(eqn->shiftsolver[i]));                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_clear);
        }
        mess_free(eqn->shiftsolver);
        eqn->shiftsolver=NULL;
        eqn->shifts = 0;
    }

    e->ApEINV.to_clear = 0;
    return 0;
}

/**
 * @internal
 * @brief Call the @ref mess_lrcfadi_parameter and filter shift parameters.
 * @param[in] e         input @ref mess_equation structure
 * @param[in] opt       input options
 * @param[in] stat      input status
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref parameter function
 * performs the shift parameter computation for the so1 function handles.
 * After calling @ref mess_lrcfadi_parameter shift parameters are removed which are
 * not fullfill the boundness criteria \f$ lowerbound < |p| < upperbound \f$.
 *
 * @see mess_lrcfadi_parameter
 *
 * @attention Internal use only.
 */
static int parameter(mess_equation e, mess_options opt, mess_status stat){
    MSG_FNAME(__func__);
    mess_int_t ret=0;
    mess_check_nullpointer(e);
    __glyap_so2 *eqn = (__glyap_so2 *) e->aux;
    mess_check_nullpointer(eqn);

    /*-----------------------------------------------------------------------------
     *  compute shift parameter, with stabilized or non stabilized version
     *-----------------------------------------------------------------------------*/

    if(e->parent){
        //equation was stabilised by equation_stable
        //in this case we use the stabilized function handles to compute shift parameters
        ret = mess_lrcfadi_parameter(e->parent,opt,stat);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_lrcfadi_parameter);
    }
    //equation was not stabilised by equation_stable
    ret = mess_lrcfadi_parameter(e,opt,stat);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_lrcfadi_parameter);

    if(opt->adi_shifts_p==NULL){
        MSG_INFO("No shifts found by mess_parameter.\n");
        return 0;
    }

    if(opt->adi_shifts_p==NULL){
        MSG_INFO("No shifts found by mess_parameter.\n");
        return 0;
    }

    /*-----------------------------------------------------------------------------
     *  Filter Ritz Values to small and to large ritz values
     *  instable ritz values are filtered by mess_lrcfadi_adi or mess_lrcfadi_nm
     *-----------------------------------------------------------------------------*/
    mess_int_t i,j=0;

    if(MESS_IS_REAL(opt->adi_shifts_p)){
        for(i=0;i<opt->adi_shifts_p->dim;++i){
            if( eqn->lowerbound<fabs(opt->adi_shifts_p->values[i]) && fabs(opt->adi_shifts_p->values[i])<eqn->upperbound){
                opt->adi_shifts_p->values[j]=opt->adi_shifts_p->values[i];
                ++j;
            }
        }
    }else{
        for(i=0;i<opt->adi_shifts_p->dim;++i){
            if(eqn->lowerbound<cabs(opt->adi_shifts_p->values_cpx[i]) && cabs(opt->adi_shifts_p->values_cpx[i])<eqn->upperbound) {
                opt->adi_shifts_p->values_cpx[j]=opt->adi_shifts_p->values_cpx[i];
                ++j;
            }
        }
    }

    MSG_INFO("Filtered out " MESS_PRINTF_INT " shift parameters \n",opt->adi_shifts_p->dim-j);
    ret = mess_vector_resize(opt->adi_shifts_p,j);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);
    if(!opt->adi_shifts_p->dim){
        MSG_ERROR("Filtered all shifts, no shift parameter available.\n");
        return MESS_ERROR_ARGUMENTS;
    }
    return 0;
}

/**
 * @internal
 * @brief Clear additional allocated field of @ref __glyap_so2 structure.
 * @param[in, out] e         input @ref mess_equation structure
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref clear function, clears the following additional allocated data
 *
 * <center>
 * |          field                           |  mathematical description                                                                       |
 * |:----------------------------------------:|:-----------------------------------------------------------------------------------------------:|
 * | @ref __glyap_so2.Ksolver                 | \f$ K^{-1} \f$                                                          |
 * | @ref __glyap_so2.Msolver                 | \f$ M^{-1} \f$                                                      |
 * </center>.
 *
 * @attention Internal use only.
 */
static int clear(mess_equation e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __glyap_so2 *eqn = (__glyap_so2 *) e->aux;
    mess_check_nullpointer(eqn);

    //clear solvers
    MESS_CLEAR_DIRECTS(&(eqn->Ksolver),&(eqn->Msolver));

    mess_free(eqn);
    eqn = NULL;
    e->aux = NULL;
    e->parameter = NULL;
    return 0;
}

/**
 * @brief Generate an equation object to a generalized Riccati Equation arising from a second order system.
 * @param[out] e equation object to setup
 * @param[in]  opt input options for the ADI iteration
 * @param[in]  M input   \f$ M \f$ matrix of the LTI system
 * @param[in]  D input   \f$ D \f$ matrix of the LTI system
 * @param[in]  K input   \f$ K \f$ matrix of the LTI system
 * @param[in]  B input   \f$ B \f$ matrix of the LTI system
 * @param[in]  C input   observer matrix of the LTI system
 * @param[in] lowerbound input lower bound for absolute value of shift parameter
 * @param[in] upperbound input upper bound for absolute value of shift parameter
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_equation_griccati_so2 function creates an equation object from a second order LTI system
 * \f[\begin{array}{cccc}
 *    M \ddot{x} + D  \dot{x} + K & x &=& B u \\
 *    & y & = & C x
 *   \end{array} \f]
 * where \f$ M,D,K \in \mathbb{R}^{n\times n} \f$. \n
 * Furthermore  \f$ K \f$ need to be invertible and \f$ M \f$ positive definite. \n
 * The system is implicitly transformed to its first order representation:
 *  \f[
 *   \underbrace{\left[
 *   \begin{array}{cc}
 *       I &  0 \\  0 &  M
 *   \end{array}
 *   \right]}_{\mathcal{E}}
 *   \left[
 *   \begin{array}{c}
 *      \dot{x} \\ \ddot{x}
 *   \end{array}
 *   \right] =
 *   \underbrace{ \left[
 *   \begin{array}{cc}
 *       0 &  I \\ -K & -D
 *   \end{array}
 *   \right]}_{\mathcal{A}}
 *   \left[
 *   \begin{array}{c}
 *       x  \\ \dot{x}
 *   \end{array}
 *   \right] +
 *   \underbrace{\left[
 *   \begin{array}{c}
 *       0 \\ B
 *   \end{array}
 *   \right]}_{\mathcal{B}}u
 *   \f]
 * and the Riccati solver \ref mess_lrcfadi_nm will then solve on
 * \f[ \mathcal{A}^T X \mathcal{E} + \mathcal{E}^T X \mathcal{A} - \mathcal{E}^T X \mathcal{B} \mathcal{B}^T X \mathcal{E} + \mathcal{C}^T \mathcal{C} = 0 \f]
 * or
 * \f[ \mathcal{A} X \mathcal{E}^T + \mathcal{E} X \mathcal{A}^T - \mathcal{E} X \mathcal{C}^T \mathcal{C} X \mathcal{E}^T + \mathcal{B} \mathcal{B}^T = 0 \f]
 *
 * If @ref mess_options.type  equals \ref MESS_OP_NONE
 * \f[
 * \begin{array}{cccc}
 *  \mathcal{B} & = & B & \in \mathbb{R}^{n\times t}, \\
 *  \mathcal{C} & = &
 *  \left[
 *  \begin{array}{cc}
 *      0 & C
 *  \end{array}
 *  \right]
 *  & \in \mathbb{R}^{s \times 2n},
 * \end{array}
 * \f]
 * otherwise
 * \f[
 * \begin{array}{cccc}
 * \mathcal{B} &=&
 * \left[
* \begin{array}{c}
* 0 \\ B
* \end{array}
* \right]
* & \in \mathbb{R}^{2n\times t}, \\
        * \mathcal{C} & = & C
        * & \in \mathbb{R}^{s \times n}
        * \end{array}
        * \f]
        *
        * The linearization of the second order system leads to quadratic shifts.
        * To avoid numerical problems one can defines lower and upper bounds for the absolute value of the shift.
        * Shifts that not fullfill the inequalitiy \f$ lowerbound < |p| < upperbound \f$ are automatically sorted out.
        * Common choices are \f$ lowerbound=10^{-5} \f$ and \f$ upperbound=10^{5} \f$.
        *
        *
        * \sa mess_lrcfadi_nm
        * \sa mess_equation_st
        *
        */
        int mess_equation_griccati_so2(mess_equation e, mess_options opt , mess_matrix M, mess_matrix D, mess_matrix K, mess_matrix B, mess_matrix C, double lowerbound, double upperbound){
            MSG_FNAME(__func__);
            __glyap_so2 *data;
            int ret = 0;

            /*-----------------------------------------------------------------------------
             *  check input
             *-----------------------------------------------------------------------------*/
            mess_check_nullpointer(M);
            mess_check_nullpointer(D);
            mess_check_nullpointer(K);
            mess_check_nullpointer(C);
            mess_check_nullpointer(B);
            mess_check_nullpointer(e);

            mess_check_real(M);
            mess_check_real(D);
            mess_check_real(K);
            mess_check_real(B);
            mess_check_real(C);
            mess_check_square(M);
            mess_check_square(D);
            mess_check_square(K);

            mess_check_same_size(M,D);
            mess_check_same_size(D,K);

            if(opt->type==MESS_OP_NONE){
                mess_check_same_cols(M,C);
                if (2* M->rows != B->rows ) {
                    MSG_ERROR("B has the wrong number of rows. Expected number of rows = " MESS_PRINTF_INT ",  but B->rows = " MESS_PRINTF_INT"\n", 2*M->cols, B->rows);
                    return MESS_ERROR_DIMENSION;
                }

            }else{
                mess_check_same_rows(M,B);
                if (2* M->rows != C->cols ) {
                    MSG_ERROR("C has the wrong number of cols. Expected number of cols = " MESS_PRINTF_INT ",  but C->cols = " MESS_PRINTF_INT"\n", 2*M->cols, C->cols);
                    return MESS_ERROR_DIMENSION;
                }
            }

            if(lowerbound>upperbound){
                MSG_ERROR("lowerbound=%e is higher than upperbound=%e\n",lowerbound,upperbound);
                return  MESS_ERROR_ARGUMENTS;
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
             *  Setup the equation
             *-----------------------------------------------------------------------------*/
            mess_try_alloc(data,__glyap_so2 *, sizeof(__glyap_so2));
            data->M             = M;
            data->D             = D,
                data->K             = K;
            data->lowerbound    = lowerbound;
            data->upperbound    = upperbound;
            data->shiftsolver   = NULL;
            MESS_INIT_DIRECTS(&(data->Ksolver),&(data->Msolver));
            ret = mess_direct_lu(data->K, data->Ksolver);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);
            ret = mess_direct_lu(data->M, data->Msolver);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);
            e->dim = 2*M->rows;
            e->aux = data;

            /*-----------------------------------------------------------------------------
             *  perform e->B = B        and e->C = [0,C]    (MESS_OP_NONE) or
             *          e->B = [0;B]    and e->C = C        (MESS_OP_TRANSPOSE)
             *-----------------------------------------------------------------------------*/
            mess_matrix zeros;
            ret = mess_matrix_init(&zeros);                                                                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
            if(opt->type==MESS_OP_NONE){
                ret = mess_matrix_alloc(zeros,C->rows,M->cols,C->rows*M->cols,MESS_DENSE,MESS_REAL);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
                ret = mess_matrix_init(&(e->C));                                                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
                ret = mess_matrix_cat(zeros,C,NULL,NULL,MESS_DENSE,e->C);                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);
                e->clearC = 1;
                e->B      = B;
                e->clearB = 0;
                ret = mess_matrix_clear(&zeros);                                                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear)
            }else{
                ret = mess_matrix_alloc(zeros,M->rows,B->cols,M->rows*B->cols,MESS_DENSE,MESS_REAL);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
                ret = mess_matrix_init(&(e->B));                                                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
                ret = mess_matrix_cat(zeros,NULL,B,NULL,MESS_DENSE,e->B);                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);
                e->clearB = 1;
                e->C      = C;
                e->clearC = 0;
                ret = mess_matrix_clear(&zeros);                                                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear)
            }

            e->RHS              = NULL;
            e->eqn_type         = MESS_EQN_GRICCATI;
            e->clear            = clear;
            e->EX.apply         = EX_apply;
            e->AX.generate      = NULL;
            e->AX.apply         = AX_apply;
            e->AX.clear         = NULL;
            e->AINV.generate    = NULL;//AINV_generate;
            e->AINV.apply       = AINV_apply;
            e->AINV.clear       = NULL;//AINV_clear;
            e->EINV.generate    = NULL;//EINV_generate;
            e->EINV.apply       = EINV_apply;
            e->EINV.clear       = NULL;//EINV_clear;
            e->ApEINV.generate  = ApEINV_generate;
            e->ApEINV.apply     = ApEINV_apply;
            e->ApEINV.clear     = ApEINV_clear;
            e->parameter        = parameter;

            return 0;
        }


/**
 * @brief Generate an Equation object to a generalized Lyapunov Equation arising from a second order system.
 * @param[out] e  equation object to setup
 * @param[in]  opt input  options for the ADI iteration
 * @param[in]  M input  \f$ M \f$ matrix of the LTI system
 * @param[in]  D input  \f$ D \f$ matrix of the LTI system
 * @param[in]  K input  \f$ K \f$ matrix if the LTI system
 * @param[in]  B input   right hand side \f$ B \f$ of the LTI system
 * @param[in] lowerbound input lower bound for absolute value of shift parameter
 * @param[in] upperbound input upper bound for absolute value of shift parameter
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_equation_glyap_so2 function creates an equation object from a second order LTI system
 * \f[
 *  \begin{array}{ccc}
 *    M\ddot{x} + D\dot{x} + K x &=& B u
 *  \end{array}
 * \f]
 * where \f$ M,D,K \in \mathbb{R}^{n\times n} \f$. \n
 * Furthermore \f$ K \f$ need to be invertible and \f$ M \f$ positive definite. \n
 * The system is implicitly transformed to its first order representation:
 *  \f[
 *    \underbrace{\left[
 *    \begin{array}{cc}
 *        I &  0 \\  0 &  M
 *    \end{array}
 *    \right]}_{\mathcal{E}}
 *    \left[
 *    \begin{array}{c}
 *       \dot{x} \\ \ddot{x}
 *    \end{array}\right] =
 *    \underbrace{\left[
 *    \begin{array}{cc}
 *       0 &  I \\ -K & -D
 *    \end{array}
 *    \right]}_{\mathcal{A}}
 *    \left[
 *    \begin{array}{c}
 *       x  \\ \dot{x}
 *    \end{array}
 *    \right] +
 *    \underbrace{ \left[
 *    \begin{array}{c}
 *       0 \\ B
 *    \end{array}
 *    \right]}_{\mathcal{B}}u
 *   \f]
 * and the \ref mess_lrcfadi_adi function will solve
 *
 * \f[ \mathcal{A} X \mathcal{E}^T + \mathcal{E} X \mathcal{A}^T + \mathcal{B} \mathcal{B}^T = 0\f]
 * or
 *\f[ \mathcal{A}^T X \mathcal{E} + \mathcal{E}^T X \mathcal{A} + \mathcal{B}^T \mathcal{B} = 0\f]
 *
 *
 * If @ref mess_options.type equals \ref MESS_OP_NONE the right hand side B must be passed in as
 * \f[ \mathcal{B} =
 *  \left[
 *  \begin{array}{c}
 *       0 \\ B
 *  \end{array}
 *  \right]
 *  \in \mathbb{R}^{2n\times t}, \f]
 * otherwise
 * \f[ \mathcal{B} =
 *  \left[
 *  \begin{array}{cc}
 *      0 & B^T
 *  \end{array}
 *   \right]\in \mathbb{R}^{t\times 2n} \f]
 *
 * The linearization of the second order system leads to quadratic shifts.
* To avoid numerical problems one can defines lower and upper bounds for the absolute value of the shift.
* Shifts that not fullfill the inequalitiy \f$ lowerbound < |p| < upperbound \f$ are automatically sorted out.
* Common choices are \f$ lowerbound=10^{-5} \f$ and \f$ upperbound=10^{5} \f$.
*
*
* \sa mess_lrcfadi_adi
* \sa mess_equation_lyap
*
*/
int mess_equation_glyap_so2(mess_equation e, mess_options opt , mess_matrix M, mess_matrix D, mess_matrix K, mess_matrix B, double lowerbound, double upperbound) {
    MSG_FNAME(__func__);
    __glyap_so2 *data ;
    int ret = 0;
    mess_int_t optype;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(e);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(M);
    mess_check_nullpointer(D);
    mess_check_nullpointer(K);
    mess_check_nullpointer(B);

    mess_check_real(M); mess_check_square(M); mess_check_same_size(M,D);
    mess_check_real(D); mess_check_square(D); mess_check_same_size(M,K);
    mess_check_real(K);
    mess_check_real(B); mess_check_dense(B);

    optype = opt->type;
    if (optype == MESS_OP_HERMITIAN) optype=MESS_OP_TRANSPOSE;

    if(optype==MESS_OP_NONE){
        if ( B->cols >= B->rows) {
            MSG_WARN("Oversized right hand side factor for LRCF-ADI (op = %s, rows=%d, cols=%d)\n ", mess_operation_t_str(optype), (int) B->rows, (int) B->cols);
        }
    }else{
        if ( B->rows >= B->cols) {
            MSG_WARN("Oversized right hand side factor for LRCF-ADI (op = %s, rows=%d, cols=%d)\n ", mess_operation_t_str(optype), (int) B->rows, (int) B->cols);
        }
    }

    if(lowerbound>upperbound){
        MSG_ERROR("lowerbound=%e is larger than upperbound=%e\n",lowerbound,upperbound);
        return  MESS_ERROR_ARGUMENTS;
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
     *  Setup the equation
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data,__glyap_so2 *, sizeof(__glyap_so2));
    data->M             = M;
    data->D             = D;
    data->K             = K;
    data->lowerbound    = lowerbound;
    data->upperbound    = upperbound;
    data->shiftsolver   = NULL;
    ret = mess_direct_init(&(data->Ksolver));               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
    ret = mess_direct_lu(data->K, data->Ksolver);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);
    ret = mess_direct_init(&(data->Msolver));               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
    ret = mess_direct_lu(data->M, data->Msolver);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);

    e->dim              = 2*M->rows;
    e->aux              = data;

    if(optype==MESS_OP_NONE){
        e->B        = B;
        e->clearB   = 0;
        e->RHS      = B;
        e->clearRHS = 0;
    }else{
        mess_matrix Bt;
        ret = mess_matrix_init(&Bt);                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_ctranspose(B,Bt);                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
        e->B        = Bt;
        e->clearB   = 1;
        e->RHS      = Bt;
        e->clearRHS = 0;
    }

    e->eqn_type         = MESS_EQN_GLYAP;
    e->clear            = clear;
    e->AX.apply         = AX_apply;
    e->EX.apply         = EX_apply;
    e->AINV.generate    = NULL;
    e->AINV.apply       = AINV_apply;
    e->AINV.clear       = NULL;
    e->EINV.generate    = NULL;
    e->EINV.apply       = EINV_apply;
    e->EINV.clear       = NULL;
    e->ApEINV.generate  = ApEINV_generate;
    e->ApEINV.apply     = ApEINV_apply;
    e->ApEINV.clear     = ApEINV_clear;
    e->parameter        = parameter;

    return 0;
}


