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
 * @file lib/dynsys/project.c
 * @brief Project a system ((G)LTI, second order) to another system ((G)LTI, second order).
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
#include "h2/irka_common.c"


/**
 * @brief Project a dynamical system (LTI or GLTI) to a LTI system.
 * @param[in] sys   input system
 * @param[in] V      input right projection matrix
 * @param[in] W      input left projection matrix
 * @param[out] red  reduced system
 * @return zero on succes or a non-zero error value
 *
 * The @ref mess_dynsys_project_to_lti function projects a (G)LTI system
 * \f[
 * \begin{array}{ccc}
 * E \dot{x} &= A x & + B u \\
 * y &= C x&
 *  \end{array}
 * \f]
 * to a LTI system
 * \f[
 * \begin{array}{ccc}
 *  \dot{x} &= A_r x & + B_r u \\
 * y &= C_r x&
 *  \end{array}
 * \f]
 * with the aid  matrices \f$ V \f$ and \f$ W \f$:
 *  \f[ \begin{array}{ccc}
 *      A_r = W^T A V ,\\
 *      B_r = W^T B, \\
 *      C_r=C V.
 *   \end{array}
 * \f]
 *
 */
int  mess_dynsys_project_to_lti ( mess_dynsys sys, mess_matrix V, mess_matrix W, mess_dynsys red )
{
    MSG_FNAME(__func__);
    mess_matrix Ar,Br,Cr;
    mess_matrix A, B,C;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sys);
    mess_check_nullpointer(V);
    mess_check_nullpointer(W);
    mess_check_nullpointer(red);

    if (!( MESS_IS_DYNSYS_LTI(sys) || MESS_IS_DYNSYS_GLTI(sys))){
        MSG_ERROR("input must be a LTI or a gLTI system.\n");
    }
    A=sys->A;
    B=sys->B;
    C=sys->C;

    ret = mess_matrix_init(&Ar); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Br); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Cr); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

    ret = __project_A(W,A,V,Ar); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_A);
    ret = __project_Bmat(W, B,Br); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_Bmat);
    ret = __project_Cmat(V,C,Cr); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_Cmat);

    ret =mess_dynsys_lti(red, Ar,Br,Cr);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dynsys_lti);


    return 0;
}       /* -----  end of function mess_dynsys_project  ----- */

/**
 * @brief Project a GLTI system to a GLTI system.
 * @param[in] sys   input system
 * @param[in] V      input right projection matrix
 * @param[in] W  input left projection matrix
 * @param[out] red  reduced system
 * @return zero on succes or a non-zero error value
 *
 * The @ref mess_dynsys_project_to_lti function projects a GLTI system
 * \f[
 * \begin{array}{ccc}
 * E \dot{x} &= A x & + B u \\
 * y &= C x&
 *  \end{array}
 * \f]
 * to another GLTI system
 * \f[
 * \begin{array}{ccc}
 * E_r \dot{x} &= A_r x & + B_r u \\
 * y &= C_r x&
 *  \end{array}
 * \f]
 * with the aid of matrices
 * \f$ V \f$ and \f$ W \f$ :
 *  \f[ \begin{array}{ccc}
 *      E_r &=& W^T E V, \\
 *      A_r &=& W^T A V , \\
 *      B_r &=& W^TB, \\
 *      C_r &=& C V.
 *     \end{array}
 * \f]
 *
 */
int  mess_dynsys_project_to_glti ( mess_dynsys sys, mess_matrix V, mess_matrix W, mess_dynsys red )
{
    MSG_FNAME(__func__);
    mess_matrix Ar,Br,Cr,Er;
    mess_matrix A, B,C,E;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sys);
    mess_check_nullpointer(V);
    mess_check_nullpointer(W);
    mess_check_nullpointer(red);

    if (!(MESS_IS_DYNSYS_GLTI(sys))){
        MSG_ERROR("input must be a gLTI system.\n");
    }
    A=sys->A;
    B=sys->B;
    C=sys->C;
    E=sys->E;

    ret = mess_matrix_init(&Ar); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Br); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Cr); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Er); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

    ret = __project_A(W,A,V,Ar); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_A);
    ret = __project_A(W,E,V,Er); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_A);
    ret = __project_Bmat(W, B,Br); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_Bmat);
    ret = __project_Cmat(V,C,Cr); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_Cmat);

    ret =mess_dynsys_glti(red, Ar,Br,Cr,Er);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dynsys_glti);


    return 0;
}       /* -----  end of function mess_dynsys_project  ----- */


/**
 * @brief Project a second order system to a first order system.
 * @param[in] sys2nd     input second order system
 * @param[in] V      input right projection matrix
 * @param[in] W      input left projection matrix
 * @param[out] sys1st   first order system
 * @return zero on succes or a non-zero error value
 *
 * The @ref mess_dynsys_project_2nd_to_1st function projects a second order system
 * \f[
 * \begin{array}{ccccc}
 * M \ddot{x}(t) + G \dot{x}(t) + & K x(t)  &=& B u (t) & \\
 * &y &=& C_p x(t)&  + C_v \dot{x}(t)
 *  \end{array}
 * \f]
 * to a first order system
 * \f[
 * \begin{array}{ccc}
 * \dot{z} &= A_r z & + B_r u \\
 * y &= C_r x&
 *  \end{array}
 * \f]
 * using matrices \f$ V \f$ and \f$ W \f$:
 * \f[\begin{array}{ccc}
 *     A_r &= & \left[
 *               \begin{array}{cc}
 *                 W_1^T & W_2^T
 *               \end{array}
 *              \right]
 *              \left[
 *               \begin{array}{cc}
 *                0 & I \\-M^{-1} K & -M^{-1}G
 *               \end{array}
 *              \right]
 *              \left[
 *               \begin{array}{c}
 *                V_1 \\ V_2
 *               \end{array}
 *              \right] \\
 *    B_r & = & \left[
 *               \begin{array}{cc}
 *                W_1^T & W_2^T
 *               \end{array}
 *              \right]
 *              \left[
 *               \begin{array}{c}
 *                0 \\ B
 *               \end{array}
 *              \right] \\
 *   C_r & = & \left[
 *              \begin{array}{cc}
 *                  C_p & C_v
 *              \end{array}
 *             \right]
 *             \left[
 *              \begin{array}{c}
 *               V_1 \\ V_2
 *              \end{array}
 *             \right] \\
 *   z & = & \left[
 *             \begin{array}{c}
 *              x \\ \dot{x}
 *             \end{array}
 *           \right]
 *  \end{array}.
 * \f]
 *
 */
int mess_dynsys_project_2nd_to_1st ( mess_dynsys sys2nd, mess_matrix V, mess_matrix W, mess_dynsys sys1st)
{
    MSG_FNAME(__func__);
    mess_direct Msolver;
    mess_matrix SB1, SB2;
    mess_matrix SC1, SC2;
    mess_matrix T1, T2;
    mess_matrix Ar, Br, Cr;

    int ret = 0;

    mess_check_nullpointer(sys1st);
    mess_check_nullpointer(V);
    mess_check_nullpointer(W);
    mess_check_nullpointer(sys2nd);

    if (!MESS_IS_DYNSYS_2ND(sys2nd)) {
        MSG_ERROR("The input system must be a second order system.");
        return MESS_ERROR_DYNSYS;
    }

    if ( sys2nd->dim*2 != V->rows) {
        MSG_ERROR("V doesn't have the right number of rows\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( sys2nd->dim*2 != W->rows) {
        MSG_ERROR("W doesn't have the right number of rows\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  extract SB1, SB2 from V and SC1, SC2 from W
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&SB1);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&SB2);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&SC1);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&SC2);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&T1);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&T2);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Ar);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Br);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Cr);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

    ret = mess_matrix_rowsub(V,0,sys2nd->dim-1,SB1);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
    ret = mess_matrix_rowsub(V,sys2nd->dim,sys2nd->dim*2-1,SB2);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
    ret = mess_matrix_rowsub(W,0,sys2nd->dim-1,SC1);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
    ret = mess_matrix_rowsub(W,sys2nd->dim,sys2nd->dim*2-1,SC2);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);

    ret = mess_direct_init(&Msolver);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
    ret = mess_direct_lu(sys2nd->M, Msolver);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);

    /*-----------------------------------------------------------------------------
     *  calc AR
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, SC1, MESS_OP_NONE, SB2, Ar);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, sys2nd->G, MESS_OP_NONE, SB2, T1);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_NONE, sys2nd->K, MESS_OP_NONE, SB1, T2);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_add(1,T2,1,T1);                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    ret = mess_direct_solvem(MESS_OP_NONE,Msolver, T1, T2);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, SC2, MESS_OP_NONE, T2, T1);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_add(-1, T1, 1, Ar);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);


    /*-----------------------------------------------------------------------------
     *  Br
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_solvem(MESS_OP_NONE,Msolver, sys2nd->B, T1);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, SC2, MESS_OP_NONE, T1, Br);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

    /*-----------------------------------------------------------------------------
     *  Cr
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_multiply(MESS_OP_NONE, sys2nd->Cp, MESS_OP_NONE, SB1, Cr);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mul);
    ret = mess_matrix_multiply(MESS_OP_NONE, sys2nd->Cv, MESS_OP_NONE, SB2, T1);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mul);
    ret = mess_matrix_add(1,T1,1,Cr);                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);

    ret = mess_dynsys_lti(sys1st, Ar,Br,Cr);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dynsys_lti);

    mess_matrix_clear(&SB1);
    mess_matrix_clear(&SB2);
    mess_matrix_clear(&SC1);
    mess_matrix_clear(&SC2);
    mess_matrix_clear(&T1);
    mess_matrix_clear(&T2);
    mess_direct_clear(&Msolver);
    return 0;
}       /* -----  end of function mess_dynsys_project_2nd_to_1st  ----- */

/**
 * @brief Project a second order system to a second order system.
 * @param[in] sys   input system
 * @param[in] V      input right projection matrix
 * @param[in] W  input left projection matrix
 * @param[out] red  reduced system
 * @return zero on succes or a non-zero error value
 *
 * The @ref mess_dynsys_project_2nd_to_2nd function projects a second order system
 * \f[
 *  \begin{array}{ccccc}
 *  M \ddot{x}(t) + G \dot{x}(t) + & K x(t)  &=& B u (t) & \\
 *  &y &=& C_p x(t)&  + C_v \dot{x}(t)
 *  \end{array}
 * \f]
 * to another second order system
 * \f[ \begin{array}{ccccc}
 *  M_r \ddot{x}(t) + G_r \dot{x}(t) + & K_r x(t)  &=& B_r u (t) & \\
 *  &y &=& Cp_r x(t)&  + Cv_r \dot{x}(t)
 *  \end{array}
 * \f]
 * using matrices \f$ V \f$ and \f$ W \f$:
 *  \f[ \begin{array}{ccc}
 *      M_r &=& W^TMV, \\
 *      G_r &=& W^TGV, \\
 *      K_r &=& W^TKV, \\
 *      B_r &=& W^TB,  \\
 *      Cp_r &=& CpV   \\
 *      Cv_r &=&Cv V
 *   \end{array}\f]
 *
 */
int  mess_dynsys_project_2nd_to_2nd ( mess_dynsys sys, mess_matrix V, mess_matrix W, mess_dynsys red )
{
    MSG_FNAME(__func__);
    mess_matrix Mr,Gr,Kr,Br,Cpr,Cvr;
    mess_matrix M,G,K,B,Cp,Cv;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sys);
    mess_check_nullpointer(V);
    mess_check_nullpointer(W);
    mess_check_nullpointer(red);

    if (!MESS_IS_DYNSYS_2ND(sys)) {
        MSG_ERROR("The input system must be a second order system.");
        return MESS_ERROR_DYNSYS;
    }

    if ( sys->dim != V->rows) {
        MSG_ERROR("V doesn't have the right number of rows\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( sys->dim != W->rows) {
        MSG_ERROR("W doesn't have the right number of rows\n");
        return MESS_ERROR_DIMENSION;
    }

    M=sys->M;
    G=sys->G;
    K=sys->K;
    B=sys->B;
    Cp= sys->Cp;
    Cv= sys->Cv;

    /*-----------------------------------------------------------------------------
     *  init data
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&Mr);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Gr);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Kr);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Br);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Cpr);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&Cvr);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);


    /*-----------------------------------------------------------------------------
     *  project matrices
     *-----------------------------------------------------------------------------*/
    ret = __project_A(W,M,V,Mr);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_A);
    ret = __project_A(W,G,V,Gr);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_A);
    ret = __project_A(W,K,V,Kr);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_A);

    ret = __project_Bmat(W, B,Br);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_Bmat);
    ret = __project_Cmat(V,Cp,Cpr); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_Cmat);
    ret = __project_Cmat(V,Cv,Cvr); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_Cmat);

    ret = mess_dynsys_2nd(red, Mr,Gr,Kr,Br,Cpr, Cvr);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dynsys_2nd);

    return 0;
}       /* -----  end of function mess_dynsys_project  ----- */


