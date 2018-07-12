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
 * @addtogroup interfaces_python
 * @file python/c_interface/equations.c
 * @brief Handle different equation types for @python interface.
 *
 * @author @koehlerm
 * @author @nitin
 * @author @mbehr
 *
 */


#include <Python.h>
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "pymess.h"


static int op_to_int(mess_operation_t op){
    switch(op) {
        case MESS_OP_NONE:
            return 0;
        case MESS_OP_TRANSPOSE:
            return 1;
        case MESS_OP_HERMITIAN:
            return 2;
        default:
            return 0;
    }
    return 0;
}

#define PRE_POST_FUNCTION(cname,pyname) static int cname (mess_equation e) {\
    MSG_FNAME(__func__); \
    PyObject *pyOper;\
    PyObject *pyResult;\
    PyObject *pyError = NULL;\
    /* Check input  */\
    mess_check_nullpointer(e); \
    pyOper = (PyObject*) e->aux;\
    mess_check_nullpointer(pyOper);\
    /*  Call Function  */\
    Py_INCREF(pyOper);\
    PyErr_Clear(); \
    pyResult =  PyObject_CallMethod(pyOper, #pyname, NULL);\
    Py_XDECREF(pyResult);\
    pyError = PyErr_Occurred(); \
    if ( pyError != NULL ) { \
        PyErr_Print(); \
        MSG_ERROR("The Python callback for %s failed.\n", #pyname);\
        Py_DECREF(pyOper);\
        return MESS_ERROR_PYTHON; \
    }\
    Py_DECREF(pyOper);\
    return 0;\
}

#define PRE_POST_FUNCTION2(cname,pyname) static int cname (mess_equation e, mess_vector p) {\
    MSG_FNAME(__func__); \
    PyObject *pyOper;\
    PyObject *pyResult;\
    PyObject *pyP; \
    PyObject *pyError = NULL;\
    /* Check input  */\
    mess_check_nullpointer(e); \
    pyOper = (PyObject*) e->aux;\
    mess_check_nullpointer(pyOper);\
    /*  Call Function  */\
    pyP = vector_to_python(p,1);\
    if ( pyP == NULL) { \
        MSG_ERROR("Cannot convert the parameter vector to Python.\n");\
    }\
    Py_INCREF(pyOper);\
    PyErr_Clear(); \
    pyResult =  PyObject_CallMethod(pyOper, #pyname, "O", pyP);\
    Py_XDECREF(pyResult);\
    Py_XDECREF(pyP); \
    pyError = PyErr_Occurred(); \
    if ( pyError != NULL ) { \
        PyErr_Print(); \
        Py_DECREF(pyOper);\
        return MESS_ERROR_PYTHON; \
    }\
    Py_DECREF(pyOper);\
    return 0;\
}

#define APPLY_FUNCTION(cname,pyname) static int cname (mess_equation e, mess_operation_t op, mess_matrix in, mess_matrix out) {\
    MSG_FNAME(__func__); \
    PyObject *pyOper;\
    PyObject *pyResult;\
    PyObject *pyInput = NULL;\
    PyObject *pyError = NULL;\
    mess_matrix tmp; \
    /* Check input  */\
    mess_check_nullpointer(e); \
    pyOper = (PyObject*) e->aux;\
    mess_check_nullpointer(pyOper);\
    /*  Call Function  */\
    Py_INCREF(pyOper);\
    pyInput = matrix_to_python(in); \
    if ( pyInput == NULL ) { \
        Py_DECREF(pyOper);\
        MSG_ERROR("Cannot convert input matrix into python.\n");\
        return MESS_ERROR_NULLPOINTER; \
    }\
    PyErr_Clear();\
    pyResult =  PyObject_CallMethod(pyOper, #pyname, "iO", op_to_int(op), pyInput);\
    Py_DECREF(pyInput); \
    pyError = PyErr_Occurred(); \
    if ( pyError != NULL ) { \
        PyErr_Print(); \
        Py_XDECREF(pyResult);\
        Py_DECREF(pyOper);\
        return MESS_ERROR_PYTHON; \
    }\
    /* PyObject_Print(pyResult, stdout, 0);*/\
    MESS_MATRIX_RESET(out); \
    tmp = matrix_to_c(pyResult);\
    if ( tmp == NULL ) { \
        Py_XDECREF(pyOper); Py_XDECREF(pyResult);\
        MSG_ERROR("Cannot transfer the result back to C\n"); \
        return MESS_ERROR_NULLPOINTER;\
    }\
    memcpy(out, tmp, sizeof(struct mess_matrix_st)); \
    memset(tmp, 0, sizeof(struct mess_matrix_st)); \
    mess_matrix_clear(&tmp);\
    Py_XDECREF(pyResult);\
    Py_DECREF(pyOper);\
    return 0;\
}

#define APPLY_FUNCTION2(cname,pyname) static int cname (mess_equation e, mess_operation_t op, mess_double_cpx_t p, mess_int_t idx_p,  mess_matrix in, mess_matrix out) {\
    MSG_FNAME(__func__); \
    PyObject *pyOper;\
    PyObject *pyResult;\
    PyObject *pyInput = NULL;\
    PyObject *pyError = NULL;\
    Py_complex p_tmp; \
    mess_matrix tmp; \
    /* Check input  */\
    mess_check_nullpointer(e); \
    pyOper = (PyObject*) e->aux;\
    mess_check_nullpointer(pyOper);\
    /*  Call Function  */\
    p_tmp.real = creal(p);\
    p_tmp.imag = cimag(p);\
    Py_INCREF(pyOper);\
    pyInput = matrix_to_python(in); \
    if ( pyInput == NULL ) { \
        Py_DECREF(pyOper);\
        MSG_ERROR("Cannot convert input matrix into python.");\
        return MESS_ERROR_NULLPOINTER; \
    }\
    PyErr_Clear();\
    pyResult =  PyObject_CallMethod(pyOper, #pyname, "iDiO", op_to_int(op), &p_tmp, idx_p, pyInput);\
    Py_DECREF(pyInput); \
    pyError = PyErr_Occurred(); \
    if ( pyError != NULL ) { \
        PyErr_Print(); \
        Py_XDECREF(pyResult);\
        Py_DECREF(pyOper);\
        return MESS_ERROR_PYTHON; \
    }\
    if ( pyResult == NULL ) { \
        MSG_ERROR("The python routine does not return anything. Expect a matrix.\n"); \
        Py_XDECREF(pyResult); Py_XDECREF(pyOper); \
        return MESS_ERROR_PYTHON;\
    }\
    MESS_MATRIX_RESET(out); \
    tmp = matrix_to_c(pyResult);\
    if ( tmp == NULL ) { \
        Py_XDECREF(pyOper); Py_XDECREF(pyResult); \
        MSG_ERROR("Cannot transfer the result back to C\n"); \
        return MESS_ERROR_NULLPOINTER;\
    }\
    memcpy(out, tmp, sizeof(struct mess_matrix_st)); \
    memset(tmp, 0, sizeof(struct mess_matrix_st)); \
    mess_matrix_clear(&tmp);\
    Py_XDECREF(pyResult);\
    Py_DECREF(pyOper);\
    return 0;\
}

#define SET_IF_IMPLEMENTED(EQN_NAME,C_NAME,PY_NAME) do {\
    if ( pymess_is_reimplemented(pymess_eqn_base, obj, #PY_NAME) == 1) { \
        /*  MSG_PRINT("Set %s to %s\n", #C_NAME, #PY_NAME);*/\
        eqn->EQN_NAME = C_NAME;\
    } else { \
        eqn->EQN_NAME = NULL; \
    }\
} while (0)


PRE_POST_FUNCTION(callback_AX_generate,ax_generate);
PRE_POST_FUNCTION(callback_AX_clear,ax_clear);
APPLY_FUNCTION(callback_AX_apply,ax_apply);


PRE_POST_FUNCTION(callback_EX_generate,ex_generate);
PRE_POST_FUNCTION(callback_EX_clear,ex_clear);
APPLY_FUNCTION(callback_EX_apply,ex_apply);

PRE_POST_FUNCTION(callback_AINV_generate,ainv_generate);
PRE_POST_FUNCTION(callback_AINV_clear,ainv_clear);
APPLY_FUNCTION(callback_AINV_apply,ainv_apply);

PRE_POST_FUNCTION(callback_EINV_generate,einv_generate);
PRE_POST_FUNCTION(callback_EINV_clear,einv_clear);
APPLY_FUNCTION(callback_EINV_apply,einv_apply);

PRE_POST_FUNCTION2(callback_ApEX_generate,apex_generate);
PRE_POST_FUNCTION(callback_ApEX_clear,apex_clear);
APPLY_FUNCTION2(callback_ApEX_apply,apex_apply);

PRE_POST_FUNCTION2(callback_ApEINV_generate,apeinv_generate);
PRE_POST_FUNCTION(callback_ApEINV_clear,apeinv_clear);
APPLY_FUNCTION2(callback_ApEINV_apply,apeinv_apply);




static int parameter(mess_equation eqn , mess_options opt, mess_status stat){
    MSG_FNAME(__func__);
    PyObject *pyOper;
    PyObject *pyResult;
    PyObject *pyError = NULL;
    mess_vector tmp = NULL;
    PyObject *pyB=NULL;
    PyObject *pyK=NULL;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    pyOper = (PyObject*) eqn->aux;
    mess_check_nullpointer(pyOper);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(stat);

    /*-----------------------------------------------------------------------------
     *  check if B and K is available
     *-----------------------------------------------------------------------------*/
    Py_INCREF(pyOper);

    //check if equation was stabilized
    if ( eqn->parent && eqn->parent->B != NULL && eqn->parent->K != NULL ) {
        MSG_INFO("Get B from parent equation.");
        pyB = matrix_to_python(eqn->parent->B);
        if ( pyB == NULL ) {
            Py_DECREF(pyOper);
            MSG_ERROR("Cannot convert input matrix into python.");
            return MESS_ERROR_NULLPOINTER;
        }
        MSG_INFO("Get K from parent equation.");
        pyK = matrix_to_python(eqn->parent->K);
        if ( pyK == NULL ) {
            Py_DECREF(pyOper);
            Py_XDECREF(pyB);
            MSG_ERROR("Cannot convert input matrix into python.");
            return MESS_ERROR_NULLPOINTER;
        }
    }else{
        pyB = Py_None;
        pyK = Py_None;
    }


    /*-----------------------------------------------------------------------------
     *  call parameter routine
     *-----------------------------------------------------------------------------*/
    PyErr_Clear();
    pyResult =  PyObject_CallMethod(pyOper, "parameter", "iiOO",(int) opt->adi_shifts_arp_p, (int) opt->adi_shifts_arp_m,  pyB, pyK);

    Py_XDECREF(pyB);
    Py_XDECREF(pyK);

    pyError = PyErr_Occurred();
    if ( pyError != NULL ) {
        PyErr_Print();
        Py_XDECREF(pyResult);
        Py_DECREF(pyOper);
        return MESS_ERROR_PYTHON;
    }

    if ( pyResult == NULL ) {
        MSG_ERROR("The python routine does not return anything. Expect a matrix/vector or Py_None.\n");
        Py_XDECREF(pyResult); Py_XDECREF(pyOper);
        return MESS_ERROR_PYTHON;
    }else if ( pyResult==Py_None ){
        MSG_INFO("Shift parameter routine returned Py_None. Shift Parameters will be computed by C.-M.E.S.S.\n");
        if(eqn->parent){
            //equation was stabilised by equation_stable
            //in this case we use the stabilized function handles to compute shift parameters
            return mess_lrcfadi_parameter(eqn->parent,opt,stat);
        }
        //equation was not stabilised by equation_stable
        return mess_lrcfadi_parameter(eqn,opt,stat);
    }else{

        tmp = vector_to_c(pyResult);
        if ( tmp == NULL ) {
            Py_XDECREF(pyOper); Py_XDECREF(pyResult);
            MSG_ERROR("Cannot transfer the result back to C\n");
            return MESS_ERROR_NULLPOINTER;
        }
        if(opt->adi_shifts_p==NULL){
            mess_vector_init(&opt->adi_shifts_p);
            mess_vector_alloc(opt->adi_shifts_p, tmp->dim,tmp->data_type);
        }
        mess_vector_copy(tmp, opt->adi_shifts_p);
        mess_vector_clear(&tmp);

        Py_XDECREF(pyResult);
    }
    Py_DECREF(pyOper);
    return 0;

}


static int clear (mess_equation e) {
    if ( e == NULL) return 0;
    PyObject *pyOper = e->aux;
    Py_DECREF(pyOper);
    if (e->clearB){
        mess_matrix_clear(&(e->B));
    }
    e->clearB=0;
    if(e->clearC){
        mess_matrix_clear(&(e->C));
    }
    e->clearC=0;
    if(e->clearRHS){
        mess_matrix_clear(&(e->RHS));
    }
    e->clearRHS=0;
    return 0;
}


/* Typedef for the convert functions  */
typedef mess_equation (*eqn_converter)(PyObject *, mess_freelist);

/**
 * @internal
 * @brief Structure for all @pymess to @mess converters.
 * Structure for all @pymess to @mess converters.
 *
 * @attention Internal use only.
 *
 **/
typedef struct _equations_t {
    char * name;
    eqn_converter converter;
} equations_t;


/* Register all converts  */
static equations_t eqns[] = {
    {"EquationGLyap",               eqn_conv_lyap},
    {"EquationGRiccati",            eqn_conv_riccati},
    {"EquationGLyapSO1",            eqn_conv_lyap_so1},
    {"EquationGRiccatiSO1",         eqn_conv_riccati_so1},
    {"EquationGLyapSO2",            eqn_conv_lyap_so2},
    {"EquationGRiccatiSO2",         eqn_conv_riccati_so2},
    {"EquationGLyapDAE1",           eqn_conv_lyap_dae1},
    {"EquationGRiccatiDAE1",        eqn_conv_riccati_dae1},
    {"EquationGLyapDAE2",           eqn_conv_lyap_dae2},
    {"EquationGRiccatiDAE2",        eqn_conv_riccati_dae2},
    {NULL, NULL}
};

mess_equation mess_equation_from_python ( PyObject * obj , mess_freelist  mem, mess_equation_t eqn_type, int *need_callback)
{
    MSG_FNAME(__func__);
    PyObject *pymess_module = NULL;
    PyObject *pymess_eqn_base = NULL;
    PyObject *pymess_eqn_special = NULL;
    mess_equation eqn;
    int ret = 0;
    int cnt = 0;

    if (obj == NULL )
        return NULL;

    if ( eqn_type == MESS_EQN_LYAP ) eqn_type = MESS_EQN_GLYAP;
    if ( eqn_type == MESS_EQN_RICCATI ) eqn_type = MESS_EQN_GRICCATI;
    Py_INCREF(obj);


    /*-----------------------------------------------------------------------------
     *  Import Modules
     *-----------------------------------------------------------------------------*/
    pymess_module   = PyImport_ImportModule((char *) "pymess");
    pymess_eqn_base = PyObject_GetAttrString(pymess_module,(char*)"Equation");

    if ( pymess_module == NULL || pymess_eqn_base == NULL) {
        Py_XDECREF(pymess_eqn_base);
        Py_XDECREF(pymess_module);
        PyErr_SetString(PyExc_Exception, "Cannot import pymess module or the equation class.");
        return NULL;
    }
    ret = PyObject_IsInstance(obj, pymess_eqn_base);
    if ( ret <= 0) {
        MESS_Py_XDECREF(pymess_eqn_base, pymess_module);
        PyErr_SetString(PyExc_Exception, "Cannot import pymess module or the equation class.");
        return NULL;
    }
    *need_callback = 0;
    /*-----------------------------------------------------------------------------
     * Check for special equations like Lyapunov, riccati, ... defined in C
     *-----------------------------------------------------------------------------*/
    cnt = 0;
    MSG_PRINT("\n");
    while ( eqns[cnt].name != NULL ) {
        MSG_PRINT("Check for.... %s ", eqns[cnt].name);
        pymess_eqn_special = PyObject_GetAttrString(pymess_module, (char *) eqns[cnt].name);
        if ( pymess_eqn_special == NULL ) {
            MSG_ERROR("This should not be NULL, skip entry.");
            Py_XDECREF(pymess_eqn_special);
            return NULL;
            //continue;
        }

        if ( PyObject_IsInstance(obj, pymess_eqn_special) == 1 ){
            eqn = eqns[cnt].converter(obj, mem);
            if ( eqn == NULL ) {
                MESS_Py_XDECREF(pymess_eqn_special, pymess_eqn_base, pymess_module);
                PyErr_SetString(PyExc_Exception, "Cannot convert the options object.");
                return NULL;
            }
            MESS_Py_XDECREF(pymess_eqn_special, pymess_eqn_base, pymess_module);
            MSG_PRINT("OK\n");
            return eqn;
        }
        MSG_PRINT("NO\n");
        Py_DECREF(pymess_eqn_special);
        cnt++;
    }


    /* CLose the Module  */
    *need_callback = 1;
    MSG_PRINT("Callback Equation\n");
    /*-----------------------------------------------------------------------------
     *  User defined Object
     *-----------------------------------------------------------------------------*/
    mess_equation_init(&eqn);
    //memset(eqn, 0, sizeof(struct mess_equation_st));

    eqn->aux = (void *) obj;
    eqn->eqn_type = eqn_type;
    PY_GET_LONG(eqn->dim, obj, "dim");

    /* set clear function   */
    eqn->clear = clear;


    if ( eqn->eqn_type == MESS_EQN_GLYAP ){
        PyObject *pyRHS = PyObject_GetAttrString(obj, "rhs");

        eqn->RHS = matrix_to_c(pyRHS);
        Py_DECREF(pyRHS);
        if ( eqn -> RHS == NULL ) {
            MSG_ERROR("Failed to copy RHS to C-MESS\n");
            mess_equation_clear(&eqn);
            MESS_Py_XDECREF(pymess_eqn_base, pymess_module);
            return NULL;
        }
        eqn->clearRHS=1;
    } else if ( eqn->eqn_type == MESS_EQN_GRICCATI ){
        PyObject *pyB   = PyObject_GetAttrString(obj, "b");
        PyObject *pyC   = PyObject_GetAttrString(obj, "c");

        eqn->B = matrix_to_c(pyB);
        Py_DECREF(pyB);
        if ( !eqn -> B ) {
            Py_DECREF(pyC);
            MSG_ERROR("Failed to copy B to C-MESS\n");
            mess_equation_clear(&eqn);
            MESS_Py_XDECREF(pymess_eqn_base, pymess_module);
            return NULL;
        }
        eqn->clearB=1;

        eqn->C = matrix_to_c(pyC);
        Py_DECREF(pyC);
        if ( !eqn -> C ) {
            MSG_ERROR("Failed to copy C to C-MESS\n");
            mess_equation_clear(&eqn);
            MESS_Py_XDECREF(pymess_eqn_base, pymess_module);
            return NULL;
        }
        eqn->clearC=1;
    } else {
        MSG_ERROR("Unknown equation type\n");
        mess_equation_clear(&eqn);
        MESS_Py_XDECREF(pymess_eqn_base, pymess_module);
        return NULL;
    }

    SET_IF_IMPLEMENTED(AX.generate,callback_AX_generate, ax_generate);
    eqn->AX.apply = callback_AX_apply;
    SET_IF_IMPLEMENTED(AX.clear,callback_AX_clear, ax_clear);

    SET_IF_IMPLEMENTED(EX.generate,callback_EX_generate, ex_generate);
    eqn->EX.apply = callback_EX_apply;
    SET_IF_IMPLEMENTED(EX.clear,callback_EX_clear, ex_clear);

    SET_IF_IMPLEMENTED(AINV.generate,callback_AINV_generate, ainv_generate);
    eqn->AINV.apply = callback_AINV_apply;
    SET_IF_IMPLEMENTED(AINV.clear,callback_AINV_clear, ainv_clear);

    SET_IF_IMPLEMENTED(EINV.generate,callback_EINV_generate, einv_generate);
    eqn->EINV.apply= callback_EINV_apply;
    SET_IF_IMPLEMENTED(EINV.clear,callback_EINV_clear, einv_clear);

    SET_IF_IMPLEMENTED(ApEX.generate,callback_ApEX_generate, apex_generate);
    eqn->ApEX.apply = callback_ApEX_apply;
    SET_IF_IMPLEMENTED(ApEX.clear,callback_ApEX_clear, apex_clear);

    SET_IF_IMPLEMENTED(ApEINV.generate,callback_ApEINV_generate, apeinv_generate);
    eqn->ApEINV.apply = callback_ApEINV_apply;
    SET_IF_IMPLEMENTED(ApEINV.clear,callback_ApEINV_clear,apeinv_apply);

    //parameter is mandatory
    eqn->parameter = parameter;
    eqn->init_rhs = NULL;

    mess_freelist_add_mess_equation(mem,eqn);

    MESS_Py_XDECREF(pymess_eqn_base, pymess_module);

    return eqn;
}

