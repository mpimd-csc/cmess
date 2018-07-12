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
 * @file lib/dynsys/dynsysinit.c
 * @brief Initialize dynamical systems.
 * @author @koehlerm
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Return name of a dynamical system type.
 * @param[in] type  input id of the type
 * @return name of the dynamical system type, otherwise "UNKNOWN"
 *
 * The @ref mess_dynsysstr function returns the name of a dynamical system type.
 *
 */
char *  mess_dynsysstr ( unsigned short type )
{
    switch (type){
        case MESS_DYNSYS_LTI:
            return "standard state space";
        case MESS_DYNSYS_GLTI:
            return "generalized state space";
        case MESS_DYNSYS_2ND:
            return "2nd order";
        default:
            return "UNKNOWN";
    }
    return "UNKNOWN";
}       /* -----  end of function mess_dynsysstr  ----- */

/**
 * @brief Initialize a mess_dynsys object.
 * @param[in,out] sys pointer to the object
 * @return always zero
 *
 * The @ref mess_dynsys_init function initializes a mess_dynsys object.
 *
 */
int mess_dynsys_init(mess_dynsys *sys){
    MSG_FNAME(__func__);

    mess_try_alloc(*sys, mess_dynsys, sizeof(struct mess_dynsys_st));
    (*sys)->dim = 0;
    (*sys)->inputs =0;
    (*sys)->outputs =0;
    (*sys)->A = NULL;
    (*sys)->B = NULL;
    (*sys)->C = NULL;
    (*sys)->E = NULL;
    (*sys)->M = NULL;
    (*sys)->G = NULL;
    (*sys)->K = NULL;
    (*sys)->Cp = NULL;
    (*sys)->Cv = NULL;
    return 0;
}


/**
 * @brief Clean up a mess_dynsys object.
 * @param[in] sys  input pointer to the object
 * @return always zero
 *
 * The @ref mess_dynsys_clear function cleans up a mess_dynsys object.
 *
 */
int mess_dynsys_clear(mess_dynsys *sys){
    // MSG_FNAME(__func__);
    if (sys == NULL) return 0;
    if (*sys == NULL ) return 0;

    if ((*sys)->A != NULL){
        mess_matrix_clear(&(*sys)->A);
    }
    if ((*sys)->B != NULL){
        mess_matrix_clear(&(*sys)->B);
    }
    if ((*sys)->C != NULL){
        mess_matrix_clear(&(*sys)->C);
    }
    if ((*sys)->E != NULL){
        mess_matrix_clear(&(*sys)->E);
    }
    if ((*sys)->M != NULL){
        mess_matrix_clear(&(*sys)->M);
    }
    if ((*sys)->K != NULL){
        mess_matrix_clear(&(*sys)->K);
    }
    if ((*sys)->G != NULL){
        mess_matrix_clear(&(*sys)->G);
    }
    if ((*sys)->Cp != NULL){
        mess_matrix_clear(&(*sys)->Cp);
    }
    if ((*sys)->Cv != NULL){
        mess_matrix_clear(&(*sys)->Cv);
    }

    mess_free(*sys);
    *sys=NULL;
    return 0;
}


/**
 * @brief Generate a LTI system from matrices \f$ A, B\f$ and \f$ C \f$.
 * @param[out] lti  mess_dynsys object
 * @param[in] A      input system matrix
 * @param[in] B     input matrix
 * @param[in] C      input output matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_dynsys_lti function builds a mess_dynsys object for an LTI
 * system defined by \f$ A, B \f$ and \f$ C \f$:
 * \f[ \begin{array}{ccc}
 * \dot{x} &= A x & + B u \\
 * y &= C x&
 *  \end{array} .
 * \f]
 * The matrices are not copied to the object, only a pointer is used. \n
 * When the object is cleared, the matrices will be cleared too.
 *
 */
int mess_dynsys_lti(mess_dynsys lti, mess_matrix A, mess_matrix B, mess_matrix C){
    MSG_FNAME(__func__);


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(lti);
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);

    if ( B->rows != A->rows) {
        MSG_ERROR("B and A must have the same number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( C->cols != A->cols) {
        MSG_ERROR("C and A must have the same number of columns.\n");
        return MESS_ERROR_DIMENSION;
    }


    /*-----------------------------------------------------------------------------
     *  fill the struct
     *-----------------------------------------------------------------------------*/
    lti->A = A;
    lti->B = B;
    lti->C = C;
    lti->dim = A->rows;
    lti->inputs = B->cols;
    lti->outputs = C->rows;
    lti->type = MESS_DYNSYS_LTI;

    return 0;
}


/**
 * @brief Generate a LTI system from \f$ A, B\f$ and \f$ C \f$ (copy version).
 * @param[out] lti LTI system
 * @param[in] A      input system matrix
 * @param[in] B     input matrix
 * @param[in] C      input output matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_dynsys_lti_copy function builds a mess_dynsys object as LTI
 * system defined by the matrices \f$ A, B \f$ and \f$ C \f$:
 * \f[
 * \begin{array}{ccc}
 * \dot{x} &= A x & + B u \\
 * y &= C x&
 *  \end{array} .
 * \f]
 * The matrices will be copied to the object. That means \f$ lti \to A \f$ is independent from \f$ A \f$.
 *
 */
int mess_dynsys_lti_copy(mess_dynsys lti, mess_matrix A, mess_matrix B, mess_matrix C){
    MSG_FNAME(__func__);
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(lti);
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);

    if ( B->rows != A->rows) {
        MSG_ERROR("B and A must have the same number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( C->cols != A->cols) {
        MSG_ERROR("C and A must have the same number of columns.\n");
        return MESS_ERROR_DIMENSION;
    }

    ret= mess_matrix_init(&(lti->A));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret= mess_matrix_init(&(lti->B));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret= mess_matrix_init(&(lti->C));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);


    /*-----------------------------------------------------------------------------
     *  fill the struct
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_copy(A, lti->A);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(B, lti->B);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(C, lti->C);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);

    lti->dim = A->rows;
    lti->inputs = B->cols;
    lti->outputs = C->rows;
    lti->type = MESS_DYNSYS_LTI;

    return 0;
}

/**
 * @brief Generate a generalized LTI system from matrices \f$ A,B,C \f$ and \f$ E \f$.
 * @param[out] lti  mess_dynsys object
 * @param[in] A      input system matrix
 * @param[in] B     input matrix
 * @param[in] C      input output matrix
 * @param[in] E      input mass matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_dynsys_glti function builds a mess_dynsys object for a generalized LTI
 * system defined by \f$ A, B, C \f$ and \f$ E \f$:
 * \f[
 * \begin{array}{ccc}
 * E \dot{x} &= A x & + B u \\
 * y &= C x&
 *  \end{array} .
 * \f]
 * The matrices are not copied to the object, only a  pointer is used. \n
 * When the object is cleared, the matrices will be cleared too.
 *
 */
int mess_dynsys_glti(mess_dynsys lti, mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix E){
    MSG_FNAME(__func__);


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(lti);
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);
    mess_check_square(A);

    if ( B->rows != A->rows) {
        MSG_ERROR("B and A must have the same number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( C->cols != A->cols) {
        MSG_ERROR("C and A must have the same number of columns.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( E == NULL) {
        MSG_WARN("using mess_dynsys_lti to create the object\n");
        return mess_dynsys_lti(lti, A, B, C);
    } else {
        mess_check_square(E);
        if ( E->rows!=A->rows){
            MSG_ERROR("E and A must have the same dimension.\n");
            return MESS_ERROR_DIMENSION;
        }
    }

    /*-----------------------------------------------------------------------------
     *  fill the struct
     *-----------------------------------------------------------------------------*/
    lti->A = A;
    lti->B = B;
    lti->C = C;
    lti->E = E;
    lti->dim = A->rows;
    lti->inputs = B->cols;
    lti->outputs = C->rows;
    lti->type = MESS_DYNSYS_GLTI;

    return 0;
}


/**
 * @brief Generate a LTI system from A, B, C and E (copy version).
 * @param[out] lti LTI system
 * @param[in] A      input system matrix
 * @param[in] B     input matrix
 * @param[in] C      input output matrix
 * @param[in] E      input mass matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_dynsys_lti_copy function builds a mess_dynsys object as generalized LTI
 * system defined by the matrices \f$ A, B, C \f$ and \f$ E \f$:
 * \f[
 * \begin{array}{ccc}
 * E \dot{x} &= A x & + B u \\
 * y &= C x&
 *  \end{array} .
 * \f]
 * The matrices will be copied to the object.
 * That means \f$ lti \to A \f$ is independent from \f$ A \f$ and so on.
 *
 */
int mess_dynsys_glti_copy(mess_dynsys lti, mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix E){
    MSG_FNAME(__func__);
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(lti);
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);
    mess_check_square(A);

    if ( B->rows != A->rows) {
        MSG_ERROR("B and A must have the same number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( C->cols != A->cols) {
        MSG_ERROR("C and A must have the same number of columns.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( E == NULL) {
        MSG_WARN("using mess_dynsys_lti to create the object\n");
        return mess_dynsys_lti_copy(lti, A, B, C);
    } else {
        mess_check_square(E);
        if ( E->rows!=A->rows){
            MSG_ERROR("E and A must have the same dimension.\n");
            return MESS_ERROR_DIMENSION;
        }
    }

    ret= mess_matrix_init(&(lti->A));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret= mess_matrix_init(&(lti->B));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret= mess_matrix_init(&(lti->C));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret= mess_matrix_init(&(lti->E));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);


    /*-----------------------------------------------------------------------------
     *  fill the struct
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_copy(A, lti->A);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(B, lti->B);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(C, lti->C);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(E, lti->E);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);

    lti->dim = A->rows;
    lti->inputs = B->cols;
    lti->outputs = C->rows;
    lti->type = MESS_DYNSYS_GLTI;

    return 0;
}

/**
 * @brief Generate a second order system from matrices \f$ M,G,K,B,Cp \f$ and \f$ Cv \f$.
 * @param[out] sys  mess_dynsys object
 * @param[in] M      input 2nd derivative matrix
 * @param[in] G      input 1st derivative matrix
 * @param[in] K      input state space matrix
 * @param[in] B     input matrix
 * @param[in] Cp     input positon output matrix
 * @param[in] Cv     input velocity output matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_dynsys_2nd function builds a mess_dynsys object for a second order
 * system defined by \f$ M, G, K, B, C_p \f$ and \f$ C_v \f$:
 * \f[
 * \begin{array}{ccccc}
 * M \ddot{x}(t) + G \dot{x}(t) + & K x(t)  &=& B u (t) & \\
 * &y &=& C_p x(t)&  + C_v \dot{x}(t)
 *  \end{array} .
 * \f]
 * The matrices are not copied to the object, only a pointer is used. \n
 * When the object is cleared, the matrices will be cleared too.
 *
 */
int mess_dynsys_2nd(mess_dynsys sys, mess_matrix M, mess_matrix G, mess_matrix K,
        mess_matrix B, mess_matrix Cp, mess_matrix Cv){
    MSG_FNAME(__func__);


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sys);
    mess_check_nullpointer(M);
    mess_check_nullpointer(G);
    mess_check_nullpointer(K);
    mess_check_nullpointer(B);
    mess_check_nullpointer(Cp);
    mess_check_nullpointer(Cv);
    mess_check_square(M);
    mess_check_square(G);
    mess_check_square(K);

    if ( B->rows != M->rows) {
        MSG_ERROR("B and M must have the same number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( Cp->cols != M->cols) {
        MSG_ERROR("Cp and M must have the same number of columns.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( Cv->cols != M->cols) {
        MSG_ERROR("Cp and M must have the same number of columns.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( M->rows != G->rows) {
        MSG_ERROR("M and G must have the same number of rows\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( M->rows != K->rows) {
        MSG_ERROR("M and K must have the same number of rows\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  fill the struct
     *-----------------------------------------------------------------------------*/
    sys->M = M;
    sys->G = G;
    sys->K = K;
    sys->B = B;
    sys->Cp = Cp;
    sys->Cv = Cv;

    sys->dim = M->rows;
    sys->inputs = B->cols;
    sys->outputs = MESS_MAX(Cp->rows, Cv->rows);
    sys->type = MESS_DYNSYS_2ND;

    return 0;
}


/**
 * @brief Generate a second order system from matrices \f$ M,G,K,B,Cp \f$ and \f$ Cv \f$ (copy version).
 * @param[out] sys  mess_dynsys object
 * @param[in] M      input second derivative matrix
 * @param[in] G      input first derivative matrix
 * @param[in] K      input state space matrix
 * @param[in] B     input matrix
 * @param[in] Cp     input position output matrix
 * @param[in] Cv     input velocity output matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_dynsys_2nd_copy function builds a mess_dynsys object for a second order
 * system defined by \f$ M, G, K, B, C_p \f$ and \f$ C_v \f$:
 * \f[
 * \begin{array}{ccccc}
 * M \ddot{x}(t) + G \dot{x}(t) + & K x(t)  &=& B u (t) & \\
 * &y(t) &=& C_p x(t)&  + C_v \dot{x}(t)
 *  \end{array} .
 * \f]
 * The matrices are copied to the object. \n
 * It works like \ref mess_dynsys_lti_copy or \ref mess_dynsys_glti_copy.
 *
 */
int mess_dynsys_2nd_copy(mess_dynsys sys, mess_matrix M, mess_matrix G, mess_matrix K,
        mess_matrix B, mess_matrix Cp, mess_matrix Cv){
    MSG_FNAME(__func__);
    int ret = 0;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sys);
    mess_check_nullpointer(M);
    mess_check_nullpointer(G);
    mess_check_nullpointer(K);
    mess_check_nullpointer(B);
    mess_check_nullpointer(Cp);
    mess_check_nullpointer(Cv);
    mess_check_square(M);
    mess_check_square(G);
    mess_check_square(K);

    if ( B->rows != M->rows) {
        MSG_ERROR("B and M must have the same number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( Cp->cols != M->cols) {
        MSG_ERROR("Cp and M must have the same number of columns.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( Cv->cols != M->cols) {
        MSG_ERROR("Cp and M must have the same number of columns.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( M->rows != G->rows) {
        MSG_ERROR("M and G must have the same number of rows\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( M->rows != K->rows) {
        MSG_ERROR("M and K must have the same number of rows\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  Init matrices
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&(sys->M)); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&(sys->G)); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&(sys->K)); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&(sys->B)); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&(sys->Cp)); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&(sys->Cv)); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);


    /*-----------------------------------------------------------------------------
     *  fill the struct
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_copy(M, sys->M); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(K, sys->K); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(G, sys->G); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(B, sys->B); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(Cp, sys->Cp); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(Cv, sys->Cv); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);

    sys->dim = M->rows;
    sys->inputs = B->cols;
    sys->outputs = MESS_MAX(Cp->rows, Cv->rows);
    sys->type = MESS_DYNSYS_2ND;

    return 0;
}


/**
 * @brief Convert a second order system to a first order system.
 * @param[in] sys2nd    second order system ( input )
 * @param[out] sys1st   first order system (output)
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_dynsys_2nd_to_1st function converts a second order system
 * \f[
 * \begin{array}{ccccc}
 * M \ddot{x}(t) + G \dot{x}(t) + & K x(t)  &=& B u (t) & \\
 * &y &=& C_p x(t)&  + C_v \dot{x}(t)
 *  \end{array} .
 * \f]
 * to a generalized state space first order system
 * \f[ \begin{array}{cccccc}
 * \left[
 * \begin{array}{cc}
 * I & 0 \\ 0 & M
 * \end{array}
 * \right]
 * & \left[
 * \begin{array}{c}
 * \dot{x}(t) \\ \ddot{x}(t)
 * \end{array}
 * \right]
 * &=& \left[
 * \begin{array}{cc}
 * 0 & I \\ -K & -G
 * \end{array}
 * \right]
 * & \left[
 * \begin{array}{c}
 * x(t) \\ \dot{x}(t)
 * \end{array}
 * \right] & +
 * \left[
 * \begin{array}{c}
 * 0 \\ B
 * \end{array}
 * \right] u(t) \\
 * & y(t) &=&
 * \left[
 * \begin{array}{cc}
 * C_p & C_v
 * \end{array}
 * \right]
 * & \left[
 * \begin{array}{c}
 * x(t) \\ \dot{x}(t)
 * \end{array}
 * \right] &
 * \end{array} .
 * \f]
 */
int mess_dynsys_2nd_to_1st ( mess_dynsys sys2nd, mess_dynsys sys1st )
{
    MSG_FNAME(__func__);
    mess_int_t n;
    mess_matrix eye, A, B,C, E;
    int ret = 0 ;

    mess_check_nullpointer(sys1st);
    mess_check_nullpointer(sys2nd);
    if ( !MESS_IS_DYNSYS_2ND(sys2nd)) {
        MSG_ERROR("The input must be a second order system.\n");
        return MESS_ERROR_DYNSYS;
    }
    n = sys2nd->dim;

    ret = mess_matrix_init(&eye); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_eye(eye, n,n, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
    ret = mess_matrix_init(&E); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    // E = [ I 0;0 M];
    ret = mess_matrix_cat(eye, NULL, NULL, sys2nd->M, MESS_CSR, E);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);
    // A= [ 0 I; -K -G]
    ret = mess_matrix_init(&A); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_scale(-1, sys2nd->K); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
    ret = mess_matrix_scale(-1, sys2nd->G); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
    ret = mess_matrix_cat(NULL, eye, sys2nd->K, sys2nd->G, MESS_CSR, A);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_cat);
    ret = mess_matrix_scale(-1, sys2nd->K); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
    ret = mess_matrix_scale(-1, sys2nd->G); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);

    // B = [ 0; Bi];
    ret = mess_matrix_init(&B); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    {
        mess_matrix T;
        ret = mess_matrix_init(&T); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_alloc(T,sys2nd->B->rows, sys2nd->B->cols, 0, MESS_DENSE, MESS_REAL);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_cat(T,NULL,sys2nd->B, NULL, MESS_DENSE, B);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);
        mess_matrix_clear(&T);
    }
    // C=[Cp Cv]
    ret = mess_matrix_init(&C); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_cat(sys2nd->Cp,sys2nd->Cv,NULL,NULL, MESS_DENSE, C);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);
    mess_matrix_clear(&eye);

    ret = mess_dynsys_glti(sys1st, A,B,C,E);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dynsys_glti);
    return 0;
}       /* -----  end of function mess_dynsys_2nd_to_1st  ----- */

/**
 * @brief Print information about a dynamical system.
 * @param[in] sys  input dynamical system
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_dynsys_printinfo function prints basic information about
 * a dynamical system.
 *
 */
int  mess_dynsys_printinfo ( mess_dynsys sys )
{
    MSG_FNAME(__func__);
    mess_check_nullpointer(sys);

    MSG_PRINT("mess_dynsys info:\n");
    MSG_PRINT("-> state space dimension: " MESS_PRINTF_INT "\n", sys->dim);
    MSG_PRINT("-> number of inputs:      " MESS_PRINTF_INT "\n", sys->inputs);
    MSG_PRINT("-> number of outputs:     " MESS_PRINTF_INT "\n", sys->outputs);
    MSG_PRINT("-> type:                  %s\n", mess_dynsysstr(sys->type));
    return 0;
}       /* -----  end of function mess_dynsys_printinfo  ----- */

