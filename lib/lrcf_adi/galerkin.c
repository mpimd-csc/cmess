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
 * @file lib/lrcf_adi/galerkin.c
 * @brief Galerkin projection for the @ref mess_lrnm method.
 * @author @koehlerm
 */


#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "mess/mess.h"
#include "mess/error_macro.h"


#define IS_REAL(X) (cimag((X))==0.0)


/**
 * @brief Refine a solution of a matrix equation using Galerkin projection.
 * @param[in] eqn input      equation object of the underlying matrix equation
 * @param[in] opt input  operations object containing additional options
 * @param[in,out] Z  factor \f$ Z \f$ of the solution \f$ X \approx ZZ^T \f$, updated with the refinement
 * @param[in] eqn_type input    type of equation (Lyapunov, Riccati)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_galerkin function projects an equation to the subspace spanned by the
 * matrix \f$ Z \f$. After the projection the small equation is solved with a classical dense
 * matrix equation solver and the correct solution of the small equation is lifted
 * back to the original space.\n
 * It supports standard and generalized Riccati Equations.
 *
 * The projected Riccati Equation will always be solved using a dense Newton-method.\n
 * The function \ref mess_dense_nm_gmare is called.
 *
 *
 */
int mess_lrcfadi_galerkin ( mess_equation eqn, mess_options opt, mess_equation_t eqn_type, mess_matrix Z )
{
    MSG_FNAME(__func__);
    int ret = 0;
    int have_mass_matrix = 0;
    mess_matrix Zorth;
    mess_matrix A,B,E, TMP,XC,C,U,V;
    mess_int_t i = 0, j = 0;
    mess_vector S;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(Z);
    mess_check_real(Z);
    if ( eqn->AX.apply == NULL){
        MSG_ERROR("Need a AX.apply function in the Equation");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( Z->rows != eqn->dim) {
        MSG_ERROR("The factor has an other dimension than the equation.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( eqn->EX.apply ) {
        have_mass_matrix = 1;
    }

    MESS_CALL_IF_EXIST(eqn->AX.generate, eqn);
    if ( have_mass_matrix ){
        MESS_CALL_IF_EXIST(eqn->EX.generate, eqn);
        if ( eqn_type == MESS_EQN_LYAP)
            eqn_type = MESS_EQN_GLYAP;
        if ( eqn_type == MESS_EQN_RICCATI)
            eqn_type = MESS_EQN_GRICCATI;
    }
    switch(eqn_type){
        case MESS_EQN_RICCATI:
        case MESS_EQN_GRICCATI:
            break;
        default:
            MSG_ERROR("Equation type not supported.\n");
            return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     * get the basis and project
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&Zorth);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_orth(Z, Zorth);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_orth);

    switch(eqn_type){
        case MESS_EQN_RICCATI:
            MESS_INIT_MATRICES(&TMP,&A,&B,&C,&U,&V,&XC);
            MESS_INIT_VECTORS(&S);
            ret = mess_vector_alloc(S, Zorth->cols, MESS_REAL);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);

            if ( opt->type == MESS_OP_NONE) {
                // Ak = Zorth' * A * Zorth
                ret = eqn->AX.apply(eqn,MESS_OP_TRANSPOSE,Zorth,TMP);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->AX.apply);
            } else {
                // Ak = Zorth' * A * Zorth
                ret = eqn->AX.apply(eqn,MESS_OP_NONE,Zorth,TMP);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->AX.apply);
            }
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, Zorth, MESS_OP_NONE, TMP, A);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            // B=Zorth'*B;
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, Zorth, MESS_OP_NONE, eqn->B, TMP);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_NONE, TMP, MESS_OP_HERMITIAN, TMP, B);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            // C= CZ
            ret = mess_matrix_multiply(MESS_OP_NONE, eqn->C, MESS_OP_NONE, Zorth, TMP);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, TMP, MESS_OP_NONE, TMP, C);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            // solve
            if ( opt->type == MESS_OP_NONE) {
                // ret = mess_direct_care(A,B,C, XC);           FUNCTION_FAILURE_HANDLE (ret, (ret!=0), mess_direct_care);
                ret =   mess_dense_nm_gmare(A,NULL,B,C,XC);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dense_nm_mare);
            } else {
                // ret = mess_direct_care(A,C,B,XC);            FUNCTION_FAILURE_HANDLE (ret, (ret!=0), mess_direct_care);
                ret =   mess_dense_nm_gmare(A,NULL,C,B,XC);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dense_nm_mare);
            }

            //lift up
            ret = mess_eigen_svd(XC,S,U,V);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_svd);
            // U = U*sqrt(S)
            for (i=0; i<S->dim; i++){
                for (j=0; j < U->rows; j++){
                    U->values[i*U->ld+j] *= sqrt(S->values[i]);
                }
            }
            // Z = Zorth * U
            ret = mess_matrix_multiply(MESS_OP_NONE, Zorth, MESS_OP_NONE, U, Z);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            MESS_CLEAR_MATRICES(&TMP,&A,&B,&C,&U,&V,&XC);
            MESS_CLEAR_VECTORS(&S);
            break;

        case MESS_EQN_GRICCATI:
            MESS_INIT_MATRICES(&TMP,&A,&E,&B,&C,&U,&V,&XC);
            MESS_INIT_VECTORS(&S);
            ret = mess_vector_alloc(S, Zorth->cols, MESS_REAL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);

            // Ak = Zorth' * A * Zorth
            if ( opt->type == MESS_OP_NONE ) {
                ret = eqn->AX.apply(eqn,MESS_OP_HERMITIAN, Zorth,TMP);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->AX.apply);
            } else {
                ret = eqn->AX.apply(eqn,MESS_OP_NONE,Zorth,TMP);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->AX.apply);
            }
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, Zorth, MESS_OP_NONE, TMP, A);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

            // Ek = Zorth' * E * Zorth
            if ( opt->type == MESS_OP_NONE) {
                ret = eqn->EX.apply(eqn,MESS_OP_HERMITIAN,Zorth,TMP);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->EX.apply);
            } else {
                ret = eqn->EX.apply(eqn,MESS_OP_NONE,Zorth,TMP);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->EX.apply);
            }
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, Zorth, MESS_OP_NONE, TMP, E);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

            // B=Zorth'*B;
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, Zorth, MESS_OP_NONE, eqn->B, TMP);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_NONE, TMP, MESS_OP_HERMITIAN, TMP, B);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            // C= CZ
            ret = mess_matrix_multiply(MESS_OP_NONE, eqn->C, MESS_OP_NONE, Zorth, TMP);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, TMP, MESS_OP_NONE, TMP, C);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            // solve
            if ( opt->type == MESS_OP_NONE ) {
                ret = mess_dense_nm_gmare(A,E,B,C,XC);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dense_nm_gmare);
            } else {
                ret = mess_dense_nm_gmare(A,E,C,B,XC);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dense_nm_gmare);
            }

            //lift up
            ret = mess_eigen_svd(XC,S,U,V);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_svd);
            // U = U*sqrt(S)
            for (i=0; i<S->dim; i++){
                for (j=0; j < U->rows; j++){
                    U->values[i*U->ld+j] *= sqrt(S->values[i]);
                }
            }
            // Z = Zorth * U
            ret = mess_matrix_multiply(MESS_OP_NONE, Zorth, MESS_OP_NONE, U, Z);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            mess_matrix_clear(&TMP);
            mess_matrix_clear(&A);
            mess_matrix_clear(&E);
            mess_matrix_clear(&B);
            mess_matrix_clear(&C);
            mess_matrix_clear(&U);
            mess_matrix_clear(&V);
            mess_matrix_clear(&XC);
            mess_vector_clear(&S);
            break;
        default:
            break;
    }

    mess_matrix_clear(&Zorth);
    return 0;
}       /* -----  end of function mess_lrcfadi_galerkin  ----- */
