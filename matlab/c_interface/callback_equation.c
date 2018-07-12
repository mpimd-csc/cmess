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
 * @file matlab/c_interface/callback_equation.c
 * @brief
 * @author @mbehr
 */


#include "interface_matlab.h"

typedef struct {
    mxArray * data;
    mxArray * options;

    mxArray * function_AX_generate;
    mxArray * function_AX_apply;
    mxArray * function_AX_clear;
    mxArray * function_EX_generate;
    mxArray * function_EX_apply;
    mxArray * function_EX_clear;
    mxArray * function_AINV_generate;
    mxArray * function_AINV_apply;
    mxArray * function_AINV_clear;
    mxArray * function_EINV_generate;
    mxArray * function_EINV_apply;
    mxArray * function_EINV_clear;
    mxArray * function_ApEX_generate;
    mxArray * function_ApEX_apply;
    mxArray * function_ApEX_clear;
    mxArray * function_ApEINV_generate;
    mxArray * function_ApEINV_apply;
    mxArray * function_ApEINV_clear;

    mxArray * function_parameter;

} __matlab_equation;

static int clear(mess_equation e) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(e);
    __matlab_equation *eqn = (__matlab_equation *) e->aux;
    mess_check_nullpointer(eqn);
    if ( eqn->data)     mxDestroyArray(eqn->data);
    if ( eqn->options)  mxDestroyArray(eqn->options);
    mess_free(eqn);
    eqn = NULL;
    e->aux = NULL;
    return 0;
}

#define GENERATE_FUNCTION(name,component_extern,component_intern) static int name ( mess_equation eqn) { \
    MSG_FNAME(__func__); \
    mess_check_nullpointer(eqn); \
    mxArray *prhs[3];  \
    __matlab_equation *data = (__matlab_equation *) eqn->aux; \
    int ret = 0; \
    mess_check_nullpointer(data); \
    if ( eqn->component_extern.to_clear == 1  ) { return 0; } \
    eqn->component_extern.to_clear = 1; \
    prhs[0] = mxCreateString(#name); \
    prhs[1] = data->data; \
    prhs[2] = data->options;\
    ret = mexCallMATLAB(0, NULL, 3, prhs, "feval"); \
    if ( ret != 0 ) {\
        csc_error_message("Exection of " #name "failed.\n"); \
        return -1; \
    }\
    return 0; \
}

#define CLEAR_FUNCTION(name,component_extern,component_intern) static int name ( mess_equation eqn) { \
    MSG_FNAME(__func__); \
    mess_check_nullpointer(eqn); \
    mxArray *prhs[2];  \
    __matlab_equation *data = (__matlab_equation *) eqn->aux; \
    int ret = 0; \
    mess_check_nullpointer(data); \
    if ( eqn->component_extern.to_clear == 1  ) { return 0; } \
    eqn->component_extern.to_clear = 1; \
    prhs[0] = mxCreateString(#name); \
    prhs[1] = data->data; \
    ret = mexCallMATLAB(0, NULL, 2, prhs, "feval"); \
    if ( ret != 0 ) {\
        csc_error_message("Exection of " #name "failed.\n"); \
        return -1; \
    }\
    return 0; \
}

#define APPLY_FUNCTION(name,component_extern,component_intern) static int name ( mess_equation eqn, mess_operation_t op, mess_matrix in, mess_matrix out ) {\
    MSG_FNAME(__func__);\
    mess_check_nullpointer(eqn); \
    mxArray *prhs[5], *plhs[1]; \
    __matlab_equation *data = (__matlab_equation *) eqn->aux;\
    int ret = 0; \
    mess_check_nullpointer(data);\
    \
    prhs[0] = mxCreateString(#name);\
    prhs[1] = data->data;\
    prhs[2] = data->options;\
    prhs[3] = mess_operation_t_to_mexmess(op);\
    prhs[4] = mess_matrix_to_mexmess(in);\
    ret = mexCallMATLAB(1, plhs, 5, prhs, "feval");\
    if ( ret == 0 ) {\
        mess_matrix tmp; \
        tmp = mess_matrix_from_mexmess(plhs[0]); \
        if (tmp == NULL ) {\
            mxDestroyArray(prhs[3]);\
            mxDestroyArray(prhs[4]);\
            csc_error_message("Failed to pick up the return value.\n");\
            return -1; \
        }\
        mess_matrix_copy(tmp, out);\
        mess_matrix_clear(&tmp);\
    }\
    mxDestroyArray(prhs[3]);\
    mxDestroyArray(prhs[4]);\
    if ( ret != 0 ) {\
        csc_error_message("Application of A_apply failed.\n");\
        return -1;\
    }\
    return 0;\
}



/*-----------------------------------------------------------------------------
 *  ATTENTION MAKE SURE YOU INSERT NO ADDITIONAL SPACES, DUE TO STRINGIFICATION
 *-----------------------------------------------------------------------------*/
GENERATE_FUNCTION(AX_generate,AX,function_A_generate);
GENERATE_FUNCTION(AINV_generate,AINV,function_As_generate);
GENERATE_FUNCTION(EX_generate,EX,function_E_generate);
GENERATE_FUNCTION(EINV_generate,EINV,function_Es_generate);

APPLY_FUNCTION(AX_apply,AX,function_A_appply);
APPLY_FUNCTION(EX_apply,EX,function_E_apply);
APPLY_FUNCTION(AINV_apply,AINV,function_AINV_apply);
APPLY_FUNCTION(EINV_apply,EINV,function_EINV_apply);

CLEAR_FUNCTION(AX_clear,AX,function_A_clear);
CLEAR_FUNCTION(AINV_clear,AINV,function_AINV_clear);
CLEAR_FUNCTION(EX_clear,EX,function_E_clear);
CLEAR_FUNCTION(EINV_clear,EINV,function_EINV_clear);


/*-----------------------------------------------------------------------------
 *  Complex Operators
 *-----------------------------------------------------------------------------*/
#define GENERATE_FUNCTION2(name,component_extern,component_intern) static int name ( mess_equation eqn, mess_vector p) { \
    MSG_FNAME(__func__); \
    mess_check_nullpointer(eqn); \
    mxArray *prhs[4];  \
    __matlab_equation *data = (__matlab_equation *) eqn->aux; \
    int ret = 0; \
    mess_check_nullpointer(data); \
    if ( eqn->component_extern.to_clear == 1  ) { return 0; } \
    eqn->component_extern.to_clear = 1; \
    prhs[0] = mxCreateString(#name); \
    prhs[1] = data->data; \
    prhs[2] = data->options;\
    prhs[3] = mess_vector_to_mexmess(p);\
    ret = mexCallMATLAB(0, NULL, 4, prhs, "feval"); \
    mxDestroyArray(prhs[3]); \
    if ( ret != 0 ) {\
        csc_error_message("Exection of " #name "failed.\n"); \
        return -1; \
    }\
    return 0; \
}

#define APPLY_FUNCTION2(name,component_extern,component_intern) static int name ( mess_equation eqn, mess_operation_t op, mess_double_cpx_t p, mess_int_t idx_p, mess_matrix in, mess_matrix out ) {\
    MSG_FNAME(__func__);\
    mess_check_nullpointer(eqn); \
    mxArray *prhs[7], *plhs[1]; \
    __matlab_equation *data = (__matlab_equation *) eqn->aux;\
    int ret = 0; \
    mess_check_nullpointer(data);\
    prhs[0] = mxCreateString(#name);\
    prhs[1] = data->data;\
    prhs[2] = data->options;\
    prhs[3] = mess_operation_t_to_mexmess(op);\
    prhs[4] = mxCreateComplexScalar(p); \
    prhs[5] = mxCreateDoubleScalar(idx_p+1); \
    prhs[6] = mess_matrix_to_mexmess(in);\
    ret = mexCallMATLAB(1, plhs, 7, prhs, "feval");\
    if ( ret == 0 ) {\
        mess_matrix tmp; \
        tmp = mess_matrix_from_mexmess(plhs[0]); \
        if (tmp == NULL ) {\
            mxDestroyArray(prhs[3]);\
            mxDestroyArray(prhs[4]);\
            mxDestroyArray(prhs[5]);\
            mxDestroyArray(prhs[6]);\
            csc_error_message("Failed to pick up the return value.\n");\
            return -1; \
        }\
        mess_matrix_copy(tmp, out);\
        mess_matrix_clear(&tmp);\
    }\
    mxDestroyArray(prhs[3]);\
    mxDestroyArray(prhs[4]);\
    mxDestroyArray(prhs[5]);\
    mxDestroyArray(prhs[6]);\
    if ( ret != 0 ) {\
        csc_error_message("Application of A_apply failed.\n");\
        return -1;\
    }\
    return 0;\
}


/*-----------------------------------------------------------------------------
 *  ATTENTION MAKE SURE YOU INSERT NO ADDITIONAL SPACES, DUE TO STRINGIFICATION
 *-----------------------------------------------------------------------------*/
GENERATE_FUNCTION2(ApEX_generate,ApEX,function_ApEX_generate);
GENERATE_FUNCTION2(ApEINV_generate,ApEINV,function_ApEINV_generate);

APPLY_FUNCTION2(ApEX_apply,ApEX,function_ApEX_apply);
APPLY_FUNCTION2(ApEINV_apply,ApEINV,function_ApEINV_apply);

CLEAR_FUNCTION(ApEX_clear,ApEX,function_ApEX_clear);
CLEAR_FUNCTION(ApEINV_clear,ApEINV,function_ApEINV_clear);


static int parameter (mess_equation eqn , mess_options opt, mess_status stat )
{
    MSG_FNAME(__func__);
    mxArray *prhs[7], *plhs[1];
    mess_vector rv = NULL ;
    int ret = 0;

    mess_check_nullpointer(eqn);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(stat);
    __matlab_equation *data = (__matlab_equation *) eqn->aux;
    mess_check_nullpointer(data);

    prhs[0] = mxCreateString("parameter");
    prhs[1] = data->data;
    prhs[2] = data->options;
    prhs[3] = mxCreateDoubleScalar(opt->adi_shifts_arp_p);
    prhs[4] = mxCreateDoubleScalar(opt->adi_shifts_arp_m);

    if(eqn->parent && eqn->parent->B && eqn->parent->K){
        prhs[5] = mess_matrix_to_mexmess(eqn->parent->B);
        prhs[6] = mess_matrix_to_mexmess(eqn->parent->K);
    }else{
        prhs[5] = mxCreateDoubleMatrix(0,0, mxREAL);
        prhs[6] = mxCreateDoubleMatrix(0,0, mxREAL);
    }


    ret = mexCallMATLAB(1, plhs, 5, prhs, "feval");
    mxDestroyArray(prhs[3]);
    mxDestroyArray(prhs[4]);
    if ( ret == 0 ) {
        rv = mess_vector_from_mexmess(plhs[0]);
        if ( rv == NULL ) {
            MSG_PRINT("Shift parameter routine returned []. Shift Parameters will be computed by C.-M.E.S.S.\n");
            if(eqn->parent){
                //equation was stabilised by equation_stable
                //in this case we use the stabilized function handles to compute shift parameters
                return mess_lrcfadi_parameter(eqn->parent,opt,stat);
            }
            //equation was not stabilised by equation_stable
            return mess_lrcfadi_parameter(eqn,opt,stat);
        }

        if(opt->adi_shifts_p==NULL){
            ret = mess_vector_init(&(opt->adi_shifts_p));
            ret = mess_vector_alloc((opt->adi_shifts_p),0,MESS_COMPLEX);                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
        }
        ret = mess_vector_copy(rv, opt->adi_shifts_p);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
        mess_vector_clear(&rv);
    } else {
        return MESS_ERROR_GENERAL;
    }

    return 0;
}




int mess_callback_equation_mexmess(mess_equation eqn, mess_options opt , const  mxArray *meqn)
{
    MSG_FNAME(__func__);
    __matlab_equation *data ;
    int ret = 0;
    mess_operation_t op;

    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(meqn);

    /*-----------------------------------------------------------------------------
     *  check equation
     *-----------------------------------------------------------------------------*/
    if ( eqn->clear != NULL ) {
        eqn->clear(eqn);
    }

    if ( eqn->clearRHS) {
        mess_matrix_clear(&eqn->RHS);
        eqn->RHS=NULL;
        eqn->clearRHS = 0 ;
    }

    if ( eqn->clearB) {
        mess_matrix_clear(&eqn->B);
        eqn->B=NULL;
        eqn->clearB = 0;
    }

    if ( eqn->clearC) {
        mess_matrix_clear(&eqn->C);
        eqn->C=NULL;
        eqn->clearC = 0;
    }

    memset(eqn, 0, sizeof(struct mess_equation_st));



    /*-----------------------------------------------------------------------------
     *  get equation type
     *-----------------------------------------------------------------------------*/
    eqn->eqn_type = mess_equation_t_from_mexmess(mxGetProperty(meqn,0,"eqn_type"));

    MSG_PRINT("CALLBACK equation type = %s.\n",mess_equation_t_str(eqn->eqn_type));
    /*-----------------------------------------------------------------------------
     *  Consistency Check
     *-----------------------------------------------------------------------------*/
    op = opt->type;
    if (op == MESS_OP_HERMITIAN) op= MESS_OP_TRANSPOSE;


    /*-----------------------------------------------------------------------------
     *  get right hand side for lyapunov equation
     *-----------------------------------------------------------------------------*/

    if ( eqn->eqn_type == MESS_EQN_LYAP || eqn->eqn_type == MESS_EQN_GLYAP){
        if ( op == MESS_OP_TRANSPOSE) {
            mess_matrix tmp;
            ret = mess_matrix_clear(&eqn->RHS);
            ret = mess_matrix_init(&eqn->RHS);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
            tmp = mess_matrix_from_mexmess(mxGetProperty(meqn, 0, "RHS"));
            if ( tmp == NULL ) {
                csc_error_message("failed to get the right hand side.\n");
            }
            eqn->clearRHS = 1;
            ret = mess_matrix_ctranspose(tmp, eqn->RHS); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
            mess_matrix_clear(&tmp);
        }else{
            eqn->RHS = mess_matrix_from_mexmess(mxGetProperty(meqn, 0, "RHS"));
            if (eqn->RHS == NULL ) {
                csc_error_message("Failed to get the right hand side from the equation.\n");
            }
            eqn->clearRHS = 1;

        }
    }

    /*-----------------------------------------------------------------------------
     *  get B and C for riccati equation
     *-----------------------------------------------------------------------------*/
    if ( eqn->eqn_type == MESS_EQN_RICCATI || eqn->eqn_type == MESS_EQN_GRICCATI){
        eqn->B = mess_matrix_from_mexmess(mxGetProperty(meqn, 0, "B"));
        if ( eqn->B == NULL){
            csc_error_message("Failed to get the matrix B from the equation.\n");
        }
        eqn->clearB = 1;

        eqn-> C = mess_matrix_from_mexmess(mxGetProperty(meqn,0,"C"));
        if ( eqn->C == NULL){
            csc_error_message("Failed to get the matrix C from the equation.\n");
        }
        eqn->clearC = 1;
    }

    /*-----------------------------------------------------------------------------
     *  get dimension
     *-----------------------------------------------------------------------------*/
    eqn->dim = (mess_int_t) mxGetScalar(mxGetProperty(meqn, 0, "dim"));


    /*-----------------------------------------------------------------------------
     *  Setup the equation
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data,__matlab_equation *, sizeof(__matlab_equation));
    memset(data, 0, sizeof(__matlab_equation));

    eqn->aux = data;
    opt->type = op;

    /*-----------------------------------------------------------------------------
     *  set equation type and clear function
     *-----------------------------------------------------------------------------*/
    eqn->clear     = clear;

    /*-----------------------------------------------------------------------------
     *  copy equation class from matlab
     *-----------------------------------------------------------------------------*/
    /*  Pick up data from matlab  */
    data -> data = mxDuplicateArray(meqn);
    if (data->data == NULL ) {
        csc_error_message("Failed to copy equation object\n");
    }

    /*-----------------------------------------------------------------------------
     *  copy options from matlab
     *-----------------------------------------------------------------------------*/
    data->options = mess_options_to_mexmess(opt);
    if ( data->options == NULL ) {
        csc_error_message("Can not create MATLAB options object");
    }

    /*-----------------------------------------------------------------------------
     *  Build Operations
     *-----------------------------------------------------------------------------*/
    eqn->AX.generate = AX_generate;
    eqn->AX.clear = AX_clear;
    eqn->AX.apply = AX_apply;

    eqn->EX.generate = EX_generate;
    eqn->EX.clear = EX_clear;
    eqn->EX.apply = EX_apply;

    eqn->AINV.generate = AINV_generate;
    eqn->AINV.clear = AINV_clear;
    eqn->AINV.apply = AINV_apply;

    eqn->EINV.generate = EINV_generate;
    eqn->EINV.clear = EINV_clear;
    eqn->EINV.apply = EINV_apply;

    eqn->ApEX.generate = ApEX_generate;
    eqn->ApEX.clear = ApEX_clear;
    eqn->ApEX.apply = ApEX_apply;

    eqn->ApEINV.generate = ApEINV_generate;
    eqn->ApEINV.clear = ApEINV_clear;
    eqn->ApEINV.apply = ApEINV_apply;
    eqn->parameter = parameter;

    eqn->clear = clear;

    return 0;
}



