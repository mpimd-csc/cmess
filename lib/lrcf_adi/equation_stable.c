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
 * @file lib/lrcf_adi/equation_stable.c
 * @brief Stabilize a given Lyapunov Equation with matrices \f$ U \f$ and \f$ V\f$, i.e. build an operator \f$ (A-UV^T) \f$.
 * @author @mbehr
 */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "mess/mess.h"
#include "mess/error_macro.h"

#define IS_REAL(X) (cimag((X))==0.0)

typedef struct {
    mess_matrix ApEX_tmp, AX, KX,x2,z2,Khat,Bhat;
} __stable_lyap;


static int AX_generate(mess_equation e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    int ret = 0;
    if ( e ->AX.to_clear) return 0;
    MESS_CALL_IF_EXIST(e->child->AX.generate, e->child);
    e->AX.to_clear = 1;
    return 0;

}


static int AX_apply(mess_equation e, mess_operation_t op, mess_matrix in , mess_matrix out){
    MSG_FNAME(__func__);
    int ret = 0 ;

    mess_check_nullpointer(e);
    __stable_lyap *eqn = (__stable_lyap *) e->aux;
    mess_check_nullpointer(eqn);

    if (op == MESS_OP_NONE){
        if ( e-> B && e->K) {
            ret = mess_matrix_multiply(MESS_OP_NONE, e->K, MESS_OP_NONE , in, eqn->KX);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_NONE, e->B, MESS_OP_NONE, eqn->KX, out);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = e->child->AX.apply(e->child, MESS_OP_NONE, in, eqn->AX);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), e->child->AX.apply);
            ret = mess_matrix_add(1, eqn->AX, -1, out);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
        } else {
            ret = e->child->AX.apply(e->child, MESS_OP_NONE, in, out);                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), e->child->AX.apply);
        }
    } else {
        if ( e->B && e->K ) {
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, e->B, MESS_OP_NONE , in, eqn->KX);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, e->K, MESS_OP_NONE, eqn->KX, out);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = e->child->AX.apply(e->child, MESS_OP_HERMITIAN, in, eqn->AX);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), e->child->AX.apply);
            ret = mess_matrix_add(1, eqn->AX, -1, out);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
        } else {
            ret = e->child->AX.apply(e->child, MESS_OP_HERMITIAN, in, out);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), e->child->Ax.apply);
        }
    }
    return 0;
}


static int AX_clear(mess_equation e){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    if ( e ->AX.to_clear == 0 ) return 0;
    e->AX.to_clear = 0;
    return 0;
}


static int EX_generate(mess_equation e){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    int ret = 0;
    if ( e->EX.to_clear) return 0;
    MESS_CALL_IF_EXIST(e->child->EX.generate,e->child);
    e->EX.to_clear = 1;
    return 0;
}


static int EX_apply(mess_equation e, mess_operation_t op, mess_matrix  in, mess_matrix out){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    return e->child->EX.apply(e->child, op, in, out);
}


static int EX_clear(mess_equation e){
    MSG_FNAME(__func__);
    int ret =0;
    mess_check_nullpointer(e);
    if(e->EX.to_clear==0) return 0;
    MESS_CALL_IF_EXIST(e->child->EX.clear,e->child);
    e->EX.to_clear = 0;
    return 0;
}

static int EINV_generate(mess_equation e){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    if ( e->EINV.to_clear) return 0;
    MESS_CALL_IF_EXIST(e->child->EINV.generate,e->child);
    e->EINV.to_clear = 1;
    return 0;
}


static int EINV_apply(mess_equation e, mess_operation_t op, mess_matrix  in, mess_matrix out){
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    return e->child->EINV.apply(e->child, op, in, out);
}


static int EINV_clear(mess_equation e){
    MSG_FNAME(__func__);
    int ret =0;
    mess_check_nullpointer(e);
    if(e->EINV.to_clear==0) return 0;
    MESS_CALL_IF_EXIST(e->child->EINV.clear,e->child);
    e->EINV.to_clear = 0;
    return 0;
}

static int ApEINV_generate(mess_equation e, mess_vector parameters){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    mess_check_nullpointer(parameters);
    __stable_lyap *eqn = (__stable_lyap *) e->aux;
    mess_check_nullpointer(eqn);
    if ( e->ApEINV.to_clear) {
        MSG_INFO("A previous instance of ApEINV exists.\n" );
        return 0;
    }
    MESS_CALL_IF_EXIST(e->child->ApEINV.generate,e->child, parameters);
    e->ApEINV.to_clear = 1;
    return 0;
}


static int ApEINV_apply(mess_equation e, mess_operation_t op, mess_double_cpx_t p, mess_int_t idx_p, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    mess_matrix tmp, tmp2, U;
    __stable_lyap *eqn = (__stable_lyap *) e->aux;
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);

    switch(op){
        case MESS_OP_NONE:
            if ( e->B && e->K){
                //(ApE-B*K^T)^{-1} = ApE^{-1} + ApE^{-1}*B*(I-K^T*ApE^{-1}*B)^{-1}*K^{T}*ApE^{-1}//
                ret = e->child->ApEINV.apply(e->child, MESS_OP_NONE, p, idx_p, in, out);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->child->ApEINV.apply);
                MESS_INIT_MATRICES(&tmp,&tmp2,&U);
                ret = e->child->ApEINV.apply(e->child, MESS_OP_NONE, p, idx_p, e->B, tmp);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->child->ApEINV.apply);
                ret = mess_matrix_multiply(MESS_OP_NONE, e->K, MESS_OP_NONE, tmp, tmp2);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

                ret = mess_matrix_eye(U, e->B->cols, e->B->cols, MESS_DENSE);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
                ret = mess_matrix_add(-1, tmp2, 1, U);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                ret = mess_matrix_multiply(MESS_OP_NONE, e->K, MESS_OP_NONE, out, eqn->z2);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

                ret = mess_matrix_backslashm(MESS_OP_NONE,U,eqn->z2,eqn->x2);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_backslashm);
                ret = mess_matrix_multiply(MESS_OP_NONE, tmp, MESS_OP_NONE, eqn->x2, tmp2);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
                ret = mess_matrix_add(1, tmp2, 1, out);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                MESS_CLEAR_MATRICES(&tmp,&tmp2,&U);
            } else {
                ret = e->child->ApEINV.apply(e->child, MESS_OP_NONE, p, idx_p, in, out);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),e->child->ApEINV.apply);
            }
            break;
        case MESS_OP_TRANSPOSE:
        case MESS_OP_HERMITIAN:
            if ( e->B && e->K){
                //(ApE-B*K^T)^{-T} = ApE^{-T} + ApE^{-T}*K*(I-B^T*ApE^{-T}*K)^{-1}*B^{T}*ApE^{-T}//
                ret = e->child->ApEINV.apply(e->child, MESS_OP_HERMITIAN, p, idx_p, in, out);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), e->child->ApEINV.apply);
                mess_matrix KT;
                MESS_INIT_MATRICES(&tmp,&tmp2,&U,&KT);
                mess_matrix_ctranspose(e->K,KT);
                ret = e->child->ApEINV.apply(e->child, MESS_OP_HERMITIAN, p, idx_p, KT, tmp);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), e->child->ApEINV.apply);
                ret = mess_matrix_multiply(MESS_OP_HERMITIAN, e->B, MESS_OP_NONE, tmp, tmp2);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mul);

                ret = mess_matrix_eye(U, e->B->cols, e->B->cols, MESS_DENSE);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
                ret = mess_matrix_add(-1, tmp2, 1, U);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                ret = mess_matrix_multiply(MESS_OP_HERMITIAN, e->B, MESS_OP_NONE, out, eqn->z2);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

                ret = mess_matrix_backslashm(MESS_OP_NONE,U,eqn->z2,eqn->x2);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_backslashm);
                ret = mess_matrix_multiply(MESS_OP_NONE, tmp, MESS_OP_NONE, eqn->x2, tmp2);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);
                ret = mess_matrix_add(1, tmp2, 1, out);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                MESS_CLEAR_MATRICES(&tmp,&tmp2,&U,&KT);
            } else {
                ret = e->child->ApEINV.apply(e->child, MESS_OP_HERMITIAN, p, idx_p, in, out);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), e->child->ApEINV.apply);
            }

            break;
        default:
            MSG_ERROR("Unknown Operation Type\n");
            return MESS_ERROR_ARGUMENTS;
    }

    return 0;
}


static int ApEINV_clear(mess_equation e){
    MSG_FNAME(__func__);
    int ret = 0 ;
    mess_check_nullpointer(e);
    if ( e->ApEINV.to_clear == 0) return 0;
    __stable_lyap *eqn = (__stable_lyap *) e->aux;
    mess_check_nullpointer(eqn);
    if ( e->child->ApEINV.to_clear)
        MESS_CALL_IF_EXIST(e->child->ApEINV.clear,e->child);

    e->ApEINV.to_clear = 0;
    return 0;
}

static int AINV_generate(mess_equation e){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    __stable_lyap *eqn = (__stable_lyap *) e->aux;
    mess_check_nullpointer(eqn);
    if ( e->AINV.to_clear) {
        MSG_INFO("A previous instance of AINV exists.\n" );
        return 0;
    }
    MESS_CALL_IF_EXIST(e->child->AINV.generate,e->child);
    e->AINV.to_clear = 1;
    return 0;
}


static int AINV_apply(mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    mess_matrix tmp, tmp2, U;
    __stable_lyap *eqn = (__stable_lyap *) e->aux;
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);

    switch(op){
        case MESS_OP_NONE:
            if ( e->B && e->K){
                //(ApE-B*K^T)^{-1} = ApE^{-1} + ApE^{-1}*B*(I-K^T*ApE^{-1}*B)^{-1}*K^{T}*ApE^{-1}//
                ret = e->child->AINV.apply(e->child, MESS_OP_NONE, in, out);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), e->child->AINV.apply);
                MESS_INIT_MATRICES(&tmp,&tmp2,&U);
                ret = e->child->AINV.apply(e->child, MESS_OP_NONE,  e->B, tmp);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), e->child->AINV.apply);
                ret = mess_matrix_multiply(MESS_OP_NONE, e->K, MESS_OP_NONE, tmp, tmp2);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

                ret = mess_matrix_eye(U, e->B->cols, e->B->cols, MESS_DENSE);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
                ret = mess_matrix_add(-1, tmp2, 1, U);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                ret = mess_matrix_multiply(MESS_OP_NONE, e->K, MESS_OP_NONE, out, eqn->z2);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

                ret = mess_matrix_backslashm(MESS_OP_NONE,U,eqn->z2,eqn->x2);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_backslashm);
                ret = mess_matrix_multiply(MESS_OP_NONE, tmp, MESS_OP_NONE, eqn->x2, tmp2);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
                ret = mess_matrix_add(1, tmp2, 1, out);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                MESS_CLEAR_MATRICES(&tmp,&tmp2,&U);
            } else {
                ret = e->child->AINV.apply(e->child, MESS_OP_NONE,  in, out);                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),e->child->AINV.apply);
            }
            break;
        case MESS_OP_TRANSPOSE:
        case MESS_OP_HERMITIAN:
            if ( e->B && e->K){
                mess_matrix KT; MESS_INIT_MATRICES(&KT);
                mess_matrix_ctranspose(e->K,KT);
                //(ApE-B*K^T)^{-T} = ApE^{-T} + ApE^{-T}*K*(I-B^T*ApE^{-T}*K)^{-1}*B^{T}*ApE^{-T}//
                ret = e->child->AINV.apply(e->child, MESS_OP_HERMITIAN, in, out);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), e->child->AINV.apply);
                MESS_INIT_MATRICES(&tmp,&tmp2,&U);
                ret = e->child->AINV.apply(e->child, MESS_OP_HERMITIAN,  KT, tmp);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), e->child->AINV.apply);
                ret = mess_matrix_multiply(MESS_OP_HERMITIAN, e->B, MESS_OP_NONE, tmp, tmp2);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mul);

                ret = mess_matrix_eye(U, e->B->cols, e->B->cols, MESS_DENSE);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
                ret = mess_matrix_add(-1, tmp2, 1, U);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                ret = mess_matrix_multiply(MESS_OP_HERMITIAN, e->B, MESS_OP_NONE, out, eqn->z2);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

                ret = mess_matrix_backslashm(MESS_OP_NONE,U,eqn->z2,eqn->x2);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_backslashm);
                ret = mess_matrix_multiply(MESS_OP_NONE, tmp, MESS_OP_NONE, eqn->x2, tmp2);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);
                ret = mess_matrix_add(1, tmp2, 1, out);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
                MESS_CLEAR_MATRICES(&tmp,&tmp2,&U,&KT);
            } else {
                ret = e->child->AINV.apply(e->child, MESS_OP_HERMITIAN,  in, out);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), e->child->AINV.apply);
            }

            break;
        default:
            MSG_ERROR("Unknown Operation Type\n");
            return MESS_ERROR_ARGUMENTS;
    }

    return 0;
}


static int AINV_clear(mess_equation e){
    MSG_FNAME(__func__);
    int ret = 0 ;
    mess_check_nullpointer(e);
    if ( e->AINV.to_clear == 0) return 0;
    __stable_lyap *eqn = (__stable_lyap *) e->aux;
    mess_check_nullpointer(eqn);
    if ( e->child->AINV.to_clear)
        MESS_CALL_IF_EXIST(e->child->AINV.clear,e->child);

    e->AINV.to_clear = 0;
    return 0;
}


static int init_rhs(mess_equation e, mess_options opt){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(e);
    __stable_lyap *eqn = (__stable_lyap *) e->aux;
    mess_check_nullpointer(eqn);
    MESS_CALL_IF_EXIST(e->child->init_rhs,e->child,opt);
    return 0;
}


static int parameter(mess_equation e, mess_options opt, mess_status stat){
    /**
     * WE EXPLICITLY CALL THE PARAMETER FUNCTION FROM CHILD,
     * TO ALLOW THE CHILD MODIFICIATIONS:
     *  - INSIDE THEM SELF (SEE FOR EXAMPLE DAE2 FUNCTION HANDLES)
     *  - INSIDE THE PARENT STRUCTURE (equation_stable itself) (SEE FOR EXAMPLE DAE2 FUNCTION HANDLES)
     *
     * THESE OPERATIONS ARE FOR EXAMPLE:
     *  - LIFTING OPERATIONS TO B and K  (equation_stable itself) (SEE FOR EXAMPLE DAE2 FUNCTION HANDLES)
     *  - SHIFT PARAMETER FILTERING (SEE FOR EXAMPLE SO1 and SO2 FUNCTION HANDLES)
     *
     * IN OTHER WORDS, WE LET THE CHILD EQUATION KNOW, THAT WE WANT TO
     * PERFORM A SHIFT PARAMTER COMPUTATION NOW AND ALL NECESSARY
     * MODIFICATION TO FUNCTION HANDLES/DATA ect,
     * SHOULD BE PERFORMED BY THE CHILD STRUCTURE.
     */
    return e->child->parameter(e->child,opt,stat);
}


static int clear(mess_equation e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __stable_lyap *data = (__stable_lyap *) e->aux;
    mess_check_nullpointer(data);
    MESS_CLEAR_MATRICES(&data->ApEX_tmp,&data->AX, &data->KX,&data->x2,&data->z2,&data->Khat,&data->Bhat);
    mess_free(data);
    data = NULL;
    e->aux = NULL;
    return 0;
}


/**
 * @brief Generate a stabilized Lyapunov Equation object.
 * @param[out] e        generated object
 * @param[in]  opt      input options for Lyapunov Equation
 * @param[in]  child    input child equation to use
 * @param[in]  U        input matrix
 * @param[in]  V        input matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_equation_stable function generates a stabilized operator \f$ (A-UV) \f$ to solve.
 * The resulting operator uses the function handles defined in @p child and the Sherman-Morrison-Woodbury
 * formula to perform solving with the stabilized operator.
 * Please not that @p U and @p V must have the same size and @p U and @p V are referenced by the stabilized
 * equation @p e.
 *
 * \sa mess_equation_lyap
 * \sa mess_equation_glyap
 */
int mess_equation_stable(mess_equation e, mess_options opt , mess_equation child, mess_matrix U, mess_matrix V ){
    MSG_FNAME(__func__);
    __stable_lyap *data ;
    int have_mass_matrix = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(e);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(child);

    if ( child->AX.apply == NULL ) {
        MSG_ERROR("child does not provide a AX.apply function\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( child->AINV.apply == NULL ) {
        MSG_ERROR("child does not provide a AINV.apply function\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( child->ApEINV.apply == NULL ) {
        MSG_ERROR("child does not provide a ApEINV.apply function\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( child->parameter == NULL ) {
        MSG_ERROR("child does not provide a parameter function\n");
        return MESS_ERROR_ARGUMENTS;
    }


    if ( child->EX.apply || child->EINV.apply) {
        have_mass_matrix = 1;
    }

    /*-----------------------------------------------------------------------------
     * Set DEFAULT
     *-----------------------------------------------------------------------------*/
    if ( e->clear != NULL ) {
        e->clear(e->aux);
    }

    e->eqn_type             = MESS_EQN_NONE;
    e->dim                  = 0;
    e->AX.generate          = NULL;
    e->AX.apply             = NULL;
    e->AX.clear             = NULL;
    e->EX.generate          = NULL;
    e->EX.apply             = NULL;
    e->EX.clear             = NULL;
    e->AINV.generate        = NULL;
    e->AINV.apply           = NULL;
    e->AINV.clear           = NULL;
    e->EINV.generate        = NULL;
    e->EINV.apply           = NULL;
    e->EINV.clear           = NULL;
    e->ApEX.generate        = NULL;
    e->ApEX.apply           = NULL;
    e->ApEX.clear           = NULL;
    e->ApEINV.generate      = NULL;
    e->ApEINV.apply         = NULL;
    e->ApEINV.clear         = NULL;
    e->parameter            = NULL;
    e->RHS                  = NULL;
    e->B                    = NULL;
    e->C                    = NULL;
    e->K                    = NULL;
    e->clear                = NULL;
    e->aux                  = NULL;

    if(U && V){
        mess_check_same_colsrows(U,V);
    }

    /*-----------------------------------------------------------------------------
     *  Setup the equation
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data,__stable_lyap *, sizeof(__stable_lyap));
    MESS_INIT_MATRICES(&data->ApEX_tmp,&data->AX, &data->KX,&data->x2,&data->z2,&data->Khat,&data->Bhat);
    //data->ApEX_tmp        = NULL;
    //data->AX            = NULL;
    //data->KX            = NULL;
    //data->x2            = NULL;
    //data->z2            = NULL;
    //data->Khat          = NULL;
    //data->Bhat          = NULL;

    e->aux              = data;
    e->eqn_type         = child->eqn_type;
    e->dim              = child->dim;
    e->RHS              = child->RHS;
    e->B                = U;
    e->C                = NULL;
    e->K                = V;
    e->clear            = clear;
    e->AX.generate      = AX_generate;
    e->AX.apply         = AX_apply;
    e->AX.clear         = AX_clear;
    e->AINV.generate    = AINV_generate;
    e->AINV.apply       = AINV_apply;
    e->AINV.clear       = AINV_clear;
    e->ApEINV.generate  = ApEINV_generate;
    e->ApEINV.apply     = ApEINV_apply;
    e->ApEINV.clear     = ApEINV_clear;
    e->parameter        = parameter;
    e->init_rhs         = init_rhs;

    if ( have_mass_matrix ) {
        e->EX.generate   = EX_generate;
        e->EX.apply      = EX_apply;
        e->EX.clear      = EX_clear;
        e->EINV.generate = EINV_generate;
        e->EINV.apply    = EINV_apply;
        e->EINV.clear    = EINV_clear;
    }

    //connect equations together
    e->parent = NULL;
    child->child = NULL;
    e->child = child;
    child->parent = e;

    return 0;
}




