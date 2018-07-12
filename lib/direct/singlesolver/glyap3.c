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
 * @file lib/direct/singlesolver/glyap3.c
 * @brief Generate a solver for dense Lyapunov and Stein equations based on @glyap3
 * @author @koehlerm
 * @author @mbehr
 *
 * This code implements a dense Lyapunov / Stein equation solver.
 *
 * Lyapunov Equations:
 * <ul>
 * <li> Standard:  \f$AX+XA^T+Y=0\f$, \f$A^T X + X A + Y = 0\f$
 * <li> Generalized:  \f$AXE^T+EXA^T+Y=0\f$, \f$A^T XE + E^TX A + Y = 0\f$
 * </ul>
 *
 * Stein Equations:
 * <ul>
 * <li> Standard:  \f$A X A^T - X + Y = 0\f$, \f$A^T X A - X + Y = 0\f$
 * <li> Generalized:  \f$A X A^T - E X E^T + Y = 0\f$, \f$A^T X A - E^T X E + Y = 0\f$
 * </ul>
 *
 * with \f$ A, E, Y \f$ dense real matrices.
 *
 * If \f$(A,E)\f$ the eigenvalues lie in the left open halfplane the solution of the Lyapunov equation
 * is symmetric positive semidefinite provided that the right hand side is symmetric negative semidefinite.
 *
 * Solvers for the following Lyapunov equations are implemented too:
 * <ul>
 * <li> Standard:  \f$AX+XA^T+BB^T=0\f$, \f$A^T X + X A + B^TB = 0\f$
 * <li> Generalized:  \f$AXE^T+EXA^T+BB^T=0\f$, \f$A^T XE + E^TX A + B^TB = 0\f$
 * </ul>
 *
 *
 * @attention Standard Stein Equations are not supported by @glyap3, we use the solver for generalized Stein equations
 * of @glyap3.
 *
 * For references see @cite Pen97 and @cite KoeS14.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"


/*-----------------------------------------------------------------------------
 *  parameters for @glyap3
 *-----------------------------------------------------------------------------*/
#define GLYAP3_NB       64
#define GLYAP3_ISOLVE   1


/**
 *
 * @internal
 * @brief Enumeration to represent the type of supported Lyapunov / Stein equation of @glyap3.
 *
 * The @ref _glyap3_eqn_t  enumeration used to represent the type of supported Lyapunov / Stein equations of @glyap3.
 *
 * @attention Internal use only.
 **/
typedef enum {
    /** Represents a standard Lyapunov Equation \f$ AX+XA^T+C=0 \f$ or \f$ A^T X + X A + C = 0 \f$                      */
    GLYAP3_STD_LYAP = 0,
    /** Represents a generalized Lyapunov Equation \f$ AXE^T+EXA^T+C=0 \f$ or \f$ A^T XE + E^TX A + C = 0 \f$           */
    GLYAP3_GEN_LYAP = 1,
    /** Represents a standard Stein Equation \f$ A X A^T - X + C = 0 \f$ or \f$ A^T X A - X + C = 0 \f$                 */
    GLYAP3_STD_STEIN = 2,
    /** Represents a generalized Stein Equation \f$ A X A^T - E X E^T + C = 0 \f$ or \f$ A^T X A - E^T X E + C = 0 \f$  */
    GLYAP3_GEN_STEIN = 3,
} _glyap3_eqn_t;


/**
 * @internal
 * @brief Return a human readable name of a data type value.
 * @param[in] eqn_t input eqn_type value
 * @return constant pointer to a human readable string containing the name of the data type
 *
 * The @ref glyap3_eqn_t_str function translates a data type constant to the
 * corresponding string and returns a constant pointer to this string.
 *
 * @attention Internal use only.
 *
 */
const char *glyap3_eqn_t_str ( const _glyap3_eqn_t eqn_t ){
    switch (eqn_t){
        case GLYAP3_STD_LYAP:
            return "GLYAP3_STD_LYAP";
        case GLYAP3_GEN_LYAP:
            return "GLYAP3_GEN_LYAP";
        case GLYAP3_STD_STEIN:
            return "GLYAP3_STD_STEIN";
        case GLYAP3_GEN_STEIN:
            return "GLYAP3_GEN_STEIN";
        default:
            return "Unknown _glyap3_eqn_t";
    }
}


/**
 *
 * @internal
 *
 * @attention Internal use only.
 */
typedef struct _glyap3_st {
    mess_matrix Ahat;                   /** Schur Form of Matrix \f$A\f$ of Sylvester Equation. */
    mess_matrix QA;                     /** Schur Transformation Matrix of Matrix \f$A\f$. */
    mess_matrix Ehat;                   /** Schur Form of Matrix \f$E\f$ of Sylvester Equation. */
    mess_matrix QE;                     /** Schur Transformation Matrix of Matrix \f$E\f$. */
    _glyap3_eqn_t eqn_type;             /** Represents the Type of Equation. */
    int semidefinite;                   /** Flag if Eigenvalues of  Matrix Pencil  \f$(A,E)\f$ guarantee a positive semidefinite solution. */
} _glyap3;


/**
 * @internal
 * @brief Clean up handler for @ref _glyap3.
 * @param[in] solver  input pointer to the solver data
 * @return always zero
 *
 * The @ref glyap3_clear function is the clean up function.
 *
 * @attention Internal use only.
 */
static int glyap3_clear(void *solver){
    _glyap3 * sol = (_glyap3*) solver;
    if ( sol != NULL) {
        if(sol->Ahat)   mess_matrix_clear(&(sol->Ahat));
        if(sol->QA)     mess_matrix_clear(&(sol->QA));
        if(sol->Ehat)   mess_matrix_clear(&(sol->Ehat));
        if(sol->QE)     mess_matrix_clear(&(sol->QE));
        mess_free(sol);
    }
    return 0;
}


/**
 * @internal
 * @brief Solve the Lyapunov or Stein Equation using @ref _glyap3 and @glyap3
 * @param[in] op        input operation type
 * @param[in] data      input pointer to the internal data structure
 * @param[in] Y         input real symmetric right hand side
 * @param[in,out] X     solution @p X
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref glyap3_solvemx function solves a dense Lyapunov or Stein Equation and returns
 * the solution @p X.
 * Depending on the value of @ref _glyap3.eqn_type the following equations and @p op
 *
 * <ul>
 *  <li> @ref GLYAP3_STD_LYAP and @ref MESS_OP_NONE \f$ A X + XA^ T + Y = 0\f$.
 *  <li> @ref GLYAP3_STD_LYAP and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X + XA + Y = 0\f$.
 *  <li> @ref GLYAP3_GEN_LYAP and @ref MESS_OP_NONE \f$ A X E^T + E X A^ T + Y = 0\f$.
 *  <li> @ref GLYAP3_GEN_LYAP and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X E + E^T X A + Y = 0\f$.
 * </ul>
 *
 * <ul>
 *  <li> @ref GLYAP3_STD_STEIN and @ref MESS_OP_NONE \f$ A X A^T - X + Y = 0\f$.
 *  <li> @ref GLYAP3_STD_STEIN and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A X A^T - X + Y = 0 \f$.
 *  <li> @ref GLYAP3_GEN_STEIN and @ref MESS_OP_NONE \f$ A X A^T - E X E^ T + Y = 0\f$.
 *  <li> @ref GLYAP3_GEN_STEIN and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X A - E^T X E + Y = 0\f$.
 * </ul>
 *
 * @attention @c dgelyp does not support standard stein equation therefore the case @ref GLYAP3_STD_STEIN is not supported.
 *
 * @attention Internal use only.
 *
 */
static int glyap3_solvemx(mess_operation_t op, void* data, mess_matrix Y, mess_matrix X){
    MSG_FNAME(__func__);
    int ret = 0;
    _glyap3 * sol = (_glyap3*) data;

    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sol);
    mess_check_nullpointer(X);
    mess_check_nullpointer(Y);
    mess_check_real(Y);
    mess_check_same_size(sol->Ahat,Y);

    if(sol->eqn_type == GLYAP3_STD_STEIN){
        return MESS_ERROR_NOSUPPORT;
    }

    /*-----------------------------------------------------------------------------
     *  prepare dgelyp / dgglyp call
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_convert(Y,X,MESS_DENSE);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);

    char dico []    = "C";

    if(sol->eqn_type == GLYAP3_GEN_LYAP || sol->eqn_type == GLYAP3_STD_LYAP){
        strncpy(dico, "C", 2);
    }else{
        strncpy(dico, "D", 2);
    }

    char job []     = "X";
    char fact []    = "F";
    char trans []   = "N";

    if(op == MESS_OP_NONE){
        strncpy(trans, "T", 2);
    }else{
        strncpy(trans, "N", 2);
    }

    char ul []              = "U";

    mess_int_t ISOLVE       = GLYAP3_ISOLVE;
    mess_int_t NB           = GLYAP3_NB;
    double scale            = 1;
    double sep              = 0;
    double ferr             = 0;
    mess_int_t *iwork       = NULL;
    mess_int_t tmp_iwork    = 0;
    double *dwork           = NULL;
    double tmp_dwork        = 0;
    mess_int_t ldwork       = -1;
    mess_int_t info         = 0;


    /*-----------------------------------------------------------------------------
     *  workspace query and solve
     *-----------------------------------------------------------------------------*/
    if(sol->eqn_type == GLYAP3_STD_LYAP || sol->eqn_type == GLYAP3_STD_STEIN){
        // standard equation
        F77_GLOBAL(dgelyp,DGELYP)(dico, job, fact, trans, &(sol->Ahat->rows), sol->Ahat->values, &(sol->Ahat->ld),
                sol->QA->values, &(sol->QA->ld), X->values, &X->ld,
                &scale, &sep, &ferr, NULL, NULL, NULL, &tmp_dwork, &ldwork, &info);

        // allocate workspace
        ldwork = nearbyint(tmp_dwork);
        mess_try_alloc(dwork, double *, sizeof(double)*(ldwork));

        // solve
        F77_GLOBAL(dgelyp,DGELYP)(dico, job, fact, trans, &(sol->Ahat->rows), sol->Ahat->values, &(sol->Ahat->ld),
                sol->QA->values, &(sol->QA->ld), X->values, &X->ld,
                &scale, &sep, &ferr, NULL, NULL, NULL, dwork, &ldwork, &info);

    }else{
        // generalized equation
        F77_GLOBAL(dgglyp,DGGLYP)(dico, job, fact, trans, ul, &ISOLVE, &NB, &(sol->Ahat->rows), sol->Ahat->values, &(sol->Ahat->ld),
                sol->Ehat->values, &(sol->Ehat->ld), sol->QA->values, &(sol->QA->ld),
                sol->QE->values, &(sol->QE->ld), X->values, &(X->ld),
                &scale, &sep, &ferr, NULL, NULL, NULL, &tmp_iwork, &tmp_dwork, &ldwork, &info);

        // allocate workspace
        mess_try_alloc(iwork, mess_int_t *, sizeof(mess_int_t)*(tmp_iwork));
        ldwork = nearbyint(tmp_dwork);
        mess_try_alloc(dwork, double*, sizeof(double)*(ldwork));

        // solve
        F77_GLOBAL(dgglyp,DGGLYP)(dico, job, fact, trans, ul, &ISOLVE, &NB, &(sol->Ahat->rows), sol->Ahat->values, &(sol->Ahat->ld),
                sol->Ehat->values, &(sol->Ehat->ld), sol->QA->values, &(sol->QA->ld),
                sol->QE->values, &(sol->QE->ld), X->values, &(X->ld),
                &scale, &sep, &ferr, NULL, NULL, NULL, iwork, dwork, &ldwork, &info);

    }

    if ( info !=0 ) {
        MSG_ERROR("DGELYP/DGGLYP returned with error " MESS_PRINTF_INT "\n", info);
        return MESS_ERROR_LAPACK;
    }

    if (scale != 1){
        MSG_WARN("DGELYP/DGGLYP returned scale = %e to avoid overflow.\n", scale);
    }

    ret = mess_matrix_scale(-1.0/scale,X);                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_scale);

    /*-----------------------------------------------------------------------------
     *  clean memory
     *-----------------------------------------------------------------------------*/
    if(dwork) mess_free(dwork);
    if(iwork) mess_free(iwork);

    return 0;
}



/**
 * @internal
 * @brief Solve function the Lyapunov or Stein Equation using @ref _glyap3 and @glyap3
 * @param[in] data      input pointer to the internal data structure
 * @param[in] Y         input real symmetric right hand side
 * @param[in,out] X     solution X
 * @return zero on success or a non-zero error value otherwise
 *
 * @sa glyap3_solvemx
 *
 * @attention Internal use only.
 *
 */
static int glyap3_solvem(void* data, mess_matrix Y, mess_matrix X){
    return glyap3_solvemx(MESS_OP_NONE, data, Y, X);
}


/**
 * @internal
 * @brief Solve function the Lyapunov or Stein Equation using @ref _glyap3 and @glyap3
 * @param[in] data      input pointer to the internal data structure
 * @param[in] Y         input real symmetric right hand side
 * @param[in,out] X     solution X
 * @return zero on success or a non-zero error value otherwise
 *
 * @sa glyap3_solvemx
 *
 * @attention Internal use only.
 *
 */
static int glyap3_solvemt(void* data, mess_matrix Y, mess_matrix X){
    return glyap3_solvemx(MESS_OP_TRANSPOSE, data, Y, X);
}


/**
 * @internal
 * @brief Create a solver the Lyapunov or Stein Equation using @ref _glyap3 and @glyap3
 * @param[in] eqn_t     input @ref _glyap3_eqn_t type of equationyyyyyy
 * @param[in] A         input real matrix \f$A\f$
 * @param[in] E         input real matrix \f$ E\f$, (@c NULL if not wanted)
 * @param[out] solver   output solver
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_direct_create_glyap3 function generates a solver based on @glyap3 for standard / generalized Lyapunov and Stein equations.
 * The solver only wokrs correctly with symmetric right hand sides.
 * The type of equation is determined by @p eqn_t. If @p eqn_t is @ref GLYAP3_GEN_STEIN and @ref GLYAP3_GEN_LYAP the matrix
 * @p E must not points to @c NULL. All matrices have to be real and dense.
 *
 * @sa mess_direct_create_generalized_lyapunov
 * @sa mess_direct_create_generalized_stein
 *
 * @attention @c dgelyp does not support standard stein equation therefore the case @ref GLYAP3_STD_STEIN is not supported.
 *
 * @attention Internal use only.
 *
 */
int mess_direct_create_glyap3(_glyap3_eqn_t eqn_t, mess_matrix A, mess_matrix E, mess_direct solver){
    MSG_FNAME(__func__);
    int ret  = 0;
    mess_int_t i, j;
    double singular_tol  = mess_eps();
    mess_vector eva, evb;
    mess_double_cpx_t eval1, eval2;
    _glyap3 * data = NULL;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_square(A);
    mess_check_real(A);
    mess_check_nullpointer(solver);

    if(E){
        mess_check_square(E);
        mess_check_real(E);
        mess_check_same_size(A,E);
        if(eqn_t == GLYAP3_STD_LYAP || eqn_t == GLYAP3_STD_STEIN) {
            MSG_ERROR("E is given but no generalized Lyapuov / Stein equation is wanted!.\n");
            return MESS_ERROR_ARGUMENTS;
        }
    }

    /*-----------------------------------------------------------------------------
     *  alloc data
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data, _glyap3 *, sizeof(_glyap3));
    data->semidefinite = 1;
    data->eqn_type = eqn_t;
    MESS_INIT_MATRICES(&(data->Ahat), &(data->QA), &(data->Ehat), &(data->QE));
    MESS_INIT_VECTORS(&eva,&evb);

    /*-----------------------------------------------------------------------------
     * compute schur form and check eigenvalues.
     *-----------------------------------------------------------------------------*/
    if(!E){
        ret = mess_eigen_schur(A, data->Ahat, data->QA, eva);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_schur);
    }else{
        ret = mess_eigen_gschur(A, E, data->Ahat, data->Ehat, data->QA, data->QE, eva, evb, NULL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), ness_eigen_gschur);
    }

    /*-----------------------------------------------------------------------------
     *  check eigenvalues for unique solvability
     *-----------------------------------------------------------------------------*/
    switch (eqn_t){
        case GLYAP3_STD_LYAP:
            for (i=0; i<eva->dim; ++i){
                ret = mess_vector_get(eva, i, &eval1);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_get);
                for(j=0; j<eva->dim; ++j){
                    mess_vector_get(eva, j, &eval2);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_get);
                    if(cabs(eval1 + eval2) < singular_tol){
                        MSG_ERROR("Unique solvability can not be guaranted, there are two eigenvalues which sums to zero.\n"
                                "lambda_1 = %e + I%e and lambda_2 = %e + I %e\n",eval1, eval2);
                        return MESS_ERROR_ARGUMENTS;
                    }
                }
                data->semidefinite = data->semidefinite && creal(eval1) < 0;
            }
            break;
        case GLYAP3_GEN_LYAP:
            for (i=0; i<eva->dim; ++i){
                ret = mess_vector_get(eva, i, &eval1);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_get);
                for(j=0; j<evb->dim; ++j){
                    mess_vector_get(evb, j, &eval2);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_get);
                    if(cabs(eval1 + eval2) < singular_tol){
                        MSG_ERROR("Unique solvability can not be guaranted, there are two eigenvalues which sums to zero.\n"
                                "lambda_1 = %e + I%e and lambda_2 = %e + I %e\n",eval1, eval2);
                        return MESS_ERROR_ARGUMENTS;
                    }
                }
                data->semidefinite = data->semidefinite && creal(eval1/eval2) < 0;
            }
            break;

        case GLYAP3_STD_STEIN:
            for (i=0; i<eva->dim; ++i){
                ret = mess_vector_get(eva, i, &eval1);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_get);
                for(j=0; j<eva->dim; ++j){
                    mess_vector_get(eva, j, &eval2);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_get);
                    if(cabs(eval1*eval2-1) < singular_tol){
                        MSG_ERROR("Unique solvability can not be guaranted, there are two eigenvalues which product is one.\n"
                                "lambda_1 = %e + I%e and lambda_2 = %e + I %e\n",eval1, eval2);
                        return MESS_ERROR_ARGUMENTS;
                    }
                }
                data->semidefinite = data->semidefinite && cabs(eval1) < 1;
            }
            break;

        case GLYAP3_GEN_STEIN:
            for (i=0; i<eva->dim; ++i){
                ret = mess_vector_get(eva, i, &eval1);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_get);
                for(j=0; j<evb->dim; ++j){
                    mess_vector_get(evb, j, &eval2);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_get);
                    if(cabs(eval1*eval2-1) < singular_tol){
                        MSG_ERROR("Unique solvability can not be guaranted, there are two eigenvalues which product is one.\n"
                                "lambda_1 = %e + I%e and lambda_2 = %e + I %e\n",eval1, eval2);
                        return MESS_ERROR_ARGUMENTS;
                    }
                }
                data->semidefinite = data->semidefinite && cabs(eval1/eval2) < 1;
            }
            break;

        default:
            MSG_ERROR("Unkown _glyap3_eqn_t.\n");
            return MESS_ERROR_NOSUPPORT;
    }


    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_VECTORS(&eva,&evb);

    /*-----------------------------------------------------------------------------
     *  create solver
     *-----------------------------------------------------------------------------*/
    solver->rows = A->rows;
    solver->cols = A->cols;
    SET_SOLVERNAME(solver->name, glyap3_eqn_t_str(eqn_t));
    solver->data    = (void *)data;
    solver->solvem  = glyap3_solvem;
    solver->solvemt = glyap3_solvemt;
    solver->solvemh = glyap3_solvemt;
    solver->clear   = glyap3_clear;
    solver->data_type = MESS_REAL;
    return 0;

}


/**
 * @brief Generate a dense standard / generalized Lyapunov equation solver.
 * @param[in] A         input real matrix \f$A\f$
 * @param[in] E         input real matrix \f$ E\f$, (@c NULL if not wanted)
 * @param[out] solver   output solver
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_direct_create_generalized_lyapunov function generates a solver for dense standard / generalized Lyapunov equations
 * based on the @glyap3 solver. \n
 * The created solver can solve the generalized Lyapunov equations: \f$AXE^T+EXA^T+Y=0\f$ and \f$A^T XE + E^TX A + Y = 0\f$.
 * If @p E points to @c NULL  the created solver can solve standard Lyapunov equations: \f$AX+XA^T+Y=0\f$ and \f$A^T X + X A + Y = 0\f$.
 *
 * @attention In both cases the right hand side \f$Y\f$ must be real and symmetric.
 *
 */
int mess_direct_create_generalized_lyapunov(mess_matrix A, mess_matrix E, mess_direct solver){
    return mess_direct_create_glyap3(E?GLYAP3_GEN_LYAP:GLYAP3_STD_LYAP, A, E, solver);
}


/**
 * @brief Generate a dense standard / generalized Stein equation solver.
 * @param[in] A         input real matrix \f$A\f$
 * @param[in] E         input real matrix \f$E\f$, (@c NULL if not wanted)
 * @param[out] solver   output solver
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_direct_create_generalized_stein function generates a solver for dense standard / generalized Stein equations
 * based on the @glyap3 solver. \n
 * The created solver can solve the generalized Stein equations:  \f$A X A^T - E X E^T + Y = 0\f$ and \f$A^T X A - E^T X E + Y = 0\f$.
 * If @p E points to @c NULL  the created solver can solve standard Stein equations: \f$A X A^T -  X + Y = 0\f$ and \f$A^T X A - X + Y = 0\f$.
 *
 * @note @glyap3 currently does not support standard Stein equations, therefore @ref mess_direct_create_generalized_stein uses the generalized Stein
 * solver @glyap3 with \f$E=I_n\f$.
 *
 */
int mess_direct_create_generalized_stein(mess_matrix A, mess_matrix E, mess_direct solver){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_matrix eye = NULL;

    if(!E){
        mess_check_nullpointer(A);
        mess_check_real(A);
        mess_check_square(A);

        ret = mess_matrix_init(&eye);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_eye(eye, A->rows, A->cols, MESS_DENSE);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);

        ret = mess_direct_create_glyap3(GLYAP3_GEN_STEIN, A, eye, solver);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_glyap3);

        if(eye) mess_matrix_clear(&eye);
    }else{
        ret = mess_direct_create_glyap3(GLYAP3_GEN_STEIN, A, E, solver);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_glyap3);
    }

    return ret;
}



/**
 * @internal
 * @brief Solve the generalized Lyapunov Equation for Cholesky Factor.
 * @param[in] op     input determine type of equation
 * @param[in] data   input solver data
 * @param[in] B      input real low rank factor of the right hand side
 * @param[in,out] Z  input/output low rank factor of the solution \f$X\f$.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref glyapchol_solvemx function solves a dense Lyapunov Equation and returns the
 * Cholesky factor of the solution \f$X\f$. Depending on @p op the following Lyapunov equations are solved:
 *
 * Depending on the value of @ref _glyap3.eqn_type the and @p op following equations are solved
 *
 * <ul>
 *  <li> @ref GLYAP3_STD_LYAP and @ref MESS_OP_NONE \f$ A X + XA^ T + BB^T = 0\f$.
 *  <li> @ref GLYAP3_STD_LYAP and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X + XA + B^TB = 0\f$.
 *  <li> @ref GLYAP3_GEN_LYAP and @ref MESS_OP_NONE \f$ A X E^T + E X A^ T + BB^T = 0\f$.
 *  <li> @ref GLYAP3_GEN_LYAP and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X E + E^T X A + B^TB = 0\f$.
 * </ul>
 *
 * @attention Internal use only.
 *
 */
static int glyapchol_solvemx(mess_operation_t op, void* data, mess_matrix B, mess_matrix Z){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_matrix BB, ZZ;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(B);
    mess_check_real(B);
    mess_check_nullpointer(Z);

    /*-----------------------------------------------------------------------------
     *  create right hand side
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&BB,&ZZ);
    if(op==MESS_OP_NONE){
        ret = mess_matrix_multiply(MESS_OP_NONE, B, MESS_OP_HERMITIAN, B, BB);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    }else{
        ret = mess_matrix_multiply(MESS_OP_HERMITIAN, B, MESS_OP_NONE, B, BB);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    }

    /*-----------------------------------------------------------------------------
     *  solve Lyapunov equation
     *-----------------------------------------------------------------------------*/
    ret = glyap3_solvemx(op, data, BB, ZZ);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), glyap3_solvemx);

    /*-----------------------------------------------------------------------------
     *  compute cholesky factorization
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_cholfactor(ZZ,Z);                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_cholfactor);
    MESS_CLEAR_MATRICES(&BB, &ZZ);

    return 0;
}


/**
 * @internal
 * @brief Solve the generalized Lyapunov Equation for Cholesky Factor.
 * @param[in] data   input solver data
 * @param[in] B      input real low rank factor of the right hand side
 * @param[in,out] Z  input/output low rank factor of the solution \f$X\f$.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref glyapchol_solvem function solves a dense Lyapunov Equation and returns the
 * Cholesky factor of the solution \f$X\f$. Depending on @p op the following Lyapunov equations are solved:
 *
 * Depending on the value of @ref _glyap3.eqn_type the following equations are solved
 *
 * <ul>
 *  <li> @ref GLYAP3_STD_LYAP  \f$ A X + XA^ T + BB^T = 0\f$.
 *  <li> @ref GLYAP3_GEN_LYAP  \f$ A X E^T + E X A^ T + BB^T = 0\f$.
 * </ul>
 *
 *
 * @attention Internal use only.
 *
 */
static int glyapchol_solvem(void* data, mess_matrix B, mess_matrix Z){
    return glyapchol_solvemx(MESS_OP_NONE, data, B, Z);
}


/**
 * @internal
 * @brief Solve the generalized Lyapunov Equation for Cholesky Factor.
 * @param[in] data   input solver data
 * @param[in] B      input real low rank factor of the right hand side
 * @param[in,out] Z  input/output low rank factor of the solution \f$X\f$.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref glyapchol_solvemt function solves a dense Lyapunov Equation and returns the
 * Cholesky factor of the solution \f$X\f$. Depending on @p op the following Lyapunov equations are solved:
 *
 * Depending on the value of @ref _glyap3.eqn_type the following equations are solved:
 *
 ** <ul>
 *  <li> @ref GLYAP3_STD_LYAP \f$ A^T X + XA + B^TB = 0\f$.
 *  <li> @ref GLYAP3_GEN_LYAP \f$ A^T X E + E^T X A + B^TB = 0\f$.
 * </ul>
 *
 *
 * @attention Internal use only.
 *
 */
static int glyapchol_solvemt(void* data, mess_matrix B, mess_matrix Z){
    return glyapchol_solvemx(MESS_OP_TRANSPOSE, data, B, Z);
}


/**
 * @brief Generate  a dense standard / generalized Lyapunov equation solver for the cholesky factor.
 * @param[in] A         input real matrix \f$A\f$
 * @param[in] E         input real matrix \f$ E\f$, (@c NULL if not wanted)
 * @param[out] solver   output solver
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_direct_create_generalized_lyapunovchol function generates a solver for dense standard / generalized Lyapunov equations
 * based on the @glyap3 solver. \n
 * The created solver can solve the generalized Lyapunov equations: \f$AXE^T+EXA^T+Y=0\f$ and \f$A^T X E + E^TX A + B B^T = 0\f$.
 * If @p E points to @c NULL  the created solver can solve standard Lyapunov equations: \f$AX+XA^T+Y=0\f$ and \f$A^T X + X A + B B^T = 0\f$.
 * The solution is always returned as a cholesky factorization \f$ X \approx ZZ^T\f$.
 *
 * @attention In both cases the right hand side \f$Y\f$ must be real and symmetric.
 *
 */
int mess_direct_create_generalized_lyapunovchol(mess_matrix A, mess_matrix E,  mess_direct solver){
    MSG_FNAME(__func__);
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  create solver
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_create_glyap3(E?GLYAP3_GEN_LYAP:GLYAP3_STD_LYAP, A, E, solver);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_glyap3);

    /*-----------------------------------------------------------------------------
     *  Check if stable
     *-----------------------------------------------------------------------------*/
    _glyap3 * data  = (_glyap3 *)  solver->data;

    if ( !data -> semidefinite ) {
        MSG_ERROR("The matrix pencil (A,E) is unstable. A stable pencil is requiered.");
        return MESS_ERROR_EIGENVALUES;
    }

    /*-----------------------------------------------------------------------------
     *  modifiy fields that cholesky factored solvers are used
     *-----------------------------------------------------------------------------*/
    solver->rows = A->rows;
    solver->cols = A->cols;
    mess_free(solver->name);
    SET_SOLVERNAME(solver->name, __func__);
    solver->solvem = glyapchol_solvem;
    solver->solvemt = glyapchol_solvemt;
    solver->solvemh = glyapchol_solvemt;
    solver->data_type = MESS_REAL;

    return 0;
}




/**
 * @internal
 * @brief Solve a standard / generalized Lyapunov / Stein equation.
 * @param[in] eqn_t     input @ref _glyap3_eqn_t type of equation
 * @param[in] op        input determine which equation has to be solved
 * @param[in] A         input real dense matrix \f$A\f$
 * @param[in] E         input real dense matrix \f$E\f$
 * @param[in] Y         input real dense matrix symmetric \f$Y\f$
 * @param[out] Ahat     output real dense matrix Schur Form of \f$A\f$, (optional @c NULL if not wanted)
 * @param[out] QA       output real dense matrix Schur transformation matrix of \f$A\f$, (optional @c NULL if not wanted)
 * @param[out] Ehat     output real dense matrix Schur Form of \f$E\f$, (optional @c NULL if not wanted)
 * @param[out] QE       output real dense matrix Schur transformation matrix of \f$E\f$, (optional @c NULL if not wanted)
 * @param[out] X        output real dense symmetrix matrix, solution of equation \f$ X\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_glyap3x function solves the standard / generalized Lyapunov or Stein equation.
 *
 * Depending on @p eqn_t and @p op the following equation is solved:
 *
 * <ul>
 *  <li> @ref GLYAP3_STD_LYAP and @ref MESS_OP_NONE \f$ A X + XA^ T + Y = 0\f$.
 *  <li> @ref GLYAP3_STD_LYAP and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X + XA + Y = 0\f$.
 *  <li> @ref GLYAP3_GEN_LYAP and @ref MESS_OP_NONE \f$ A X E^T + E X A^ T + Y = 0\f$.
 *  <li> @ref GLYAP3_GEN_LYAP and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X E + E^T X A + Y = 0\f$.
 * </ul>
 *
 * <ul>
 *  <li> @ref GLYAP3_STD_STEIN and @ref MESS_OP_NONE \f$ A X A^T - X + Y = 0\f$.
 *  <li> @ref GLYAP3_STD_STEIN and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A X A^T - X + Y = 0 \f$.
 *  <li> @ref GLYAP3_GEN_STEIN and @ref MESS_OP_NONE \f$ A X A^T - E X E^ T + Y = 0\f$.
 *  <li> @ref GLYAP3_GEN_STEIN and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X A - E^T X E + Y = 0\f$.
 * </ul>
 * The right hand side matrix @p Y must be real and symmetric.
 *
 * In order to solve the equations a (generalized) Schur decomposition of \f$ A \f$ / \f$ (A,E) \f$ is computed internally.
 *
 * If a standard Lyapunov or Stein equation is solved then the real Schur decomposition fullfills: \f$ A = Q_A \hat{A} {Q_A}^H \f$.
 * If a generalized Lyapunov or Stein equation is solved then the real Schur decomposition fullfills: \f$ A = Q_A \hat{A} {Q_E}^H \f$ and  \f$ E = Q_A \hat{E} {Q_E}^H \f$.
 *
 * If @p Ahat, @p QA, @p Ehat or @p QE is not @c NULL then the corresponding matrix of the Schur decomposition is written to that argument.
 *
 * @attention Internal use only.
 *
 */
int mess_glyap3x(_glyap3_eqn_t eqn_t, mess_operation_t op, mess_matrix A, mess_matrix E, mess_matrix Y, mess_matrix Ahat, mess_matrix QA, mess_matrix Ehat, mess_matrix QE, mess_matrix X){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_direct solver;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_operation_type(op);

    mess_check_nullpointer(A);
    mess_check_dense(A);
    mess_check_real(A);
    mess_check_square(A);

    mess_check_nullpointer(Y);
    mess_check_dense(Y);
    mess_check_real(Y);
    mess_check_square(Y);

    mess_check_same_size(A,Y);
    mess_check_nullpointer(X);

    if(eqn_t == GLYAP3_GEN_LYAP || eqn_t == GLYAP3_GEN_STEIN){
        mess_check_nullpointer(E);
        mess_check_dense(E);
        mess_check_real(E);
        mess_check_square(E);
        mess_check_same_size(A,E);
    }

    /*-----------------------------------------------------------------------------
     *  create solver
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_init(&solver);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);


    if(eqn_t == GLYAP3_GEN_LYAP || eqn_t == GLYAP3_STD_LYAP){
        ret = mess_direct_create_generalized_lyapunov(A, E, solver);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_generalized_lyapunov);
    }else{
        ret = mess_direct_create_generalized_stein(A, E, solver);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_generalized_lyapunov);
    }


    /*-----------------------------------------------------------------------------
     *  solve equation
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_solvem(op, solver, Y, X);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);

    /*-----------------------------------------------------------------------------
     *  extract schur decomposition and copy to arguments
     *-----------------------------------------------------------------------------*/
    _glyap3* data = (_glyap3*) solver->data;
    if(Ahat){
        ret = mess_matrix_copy(data->Ahat, Ahat);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    }

    if(QA){
        ret = mess_matrix_copy(data->QA, QA);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    }

    if(eqn_t == GLYAP3_GEN_LYAP || eqn_t == GLYAP3_GEN_STEIN){
    if(Ehat){
        ret = mess_matrix_copy(data->Ehat, Ehat);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    }
    if(QE){
        ret = mess_matrix_copy(data->QE, QE);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    }
    }

    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_clear(&solver);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_clear);
    return 0;

}



/**
 * @internal
 * @brief Solve a standard / generalized Lyapunov / Stein equation from given Schur decomposition.
 * @param[in] eqn_t     input @ref _glyap3_eqn_t type of equation
 * @param[in] op        input determine which equation has to be solved
 * @param[in] Ahat      input real dense matrix Schur Form of \f$A\f$
 * @param[in] QA        input real dense matrix Schur transformation matrix of \f$A\f$
 * @param[in] Ehat      input real dense matrix Schur Form of \f$E\f$, (only necessary for generalized equation)
 * @param[in] QE        input real dense matrix Schur transformation matrix of \f$E\f$, (only necessary for generalized equation)
 * @param[in] Y         input real dense matrix symmetric \f$Y\f$
 * @param[out] X        output real dense symmetrix matrix, solution of equation \f$ X\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_tglyap3x function solves the standard / generalized Lyapunov or Stein equation.
 *
 * Depending on @p eqn_t and @p op the following equation is solved:
 *
 * <ul>
 *  <li> @ref GLYAP3_STD_LYAP and @ref MESS_OP_NONE \f$ A X + XA^ T + Y = 0\f$.
 *  <li> @ref GLYAP3_STD_LYAP and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X + XA + Y = 0\f$.
 *  <li> @ref GLYAP3_GEN_LYAP and @ref MESS_OP_NONE \f$ A X E^T + E X A^ T + Y = 0\f$.
 *  <li> @ref GLYAP3_GEN_LYAP and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X E + E^T X A + Y = 0\f$.
 * </ul>
 *
 * <ul>
 *  <li> @ref GLYAP3_STD_STEIN and @ref MESS_OP_NONE \f$ A X A^T - X + Y = 0\f$.
 *  <li> @ref GLYAP3_STD_STEIN and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A X A^T - X + Y = 0 \f$.
 *  <li> @ref GLYAP3_GEN_STEIN and @ref MESS_OP_NONE \f$ A X A^T - E X E^ T + Y = 0\f$.
 *  <li> @ref GLYAP3_GEN_STEIN and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X A - E^T X E + Y = 0\f$.
 * </ul>
 * The right hand side matrix @p Y must be real and symmetric.
 *
 * In order to solve the equations a (generalized) Schur decomposition of \f$ A \f$ / \f$ (A,E) \f$ is requiered.
 *
 * If a standard Lyapunov or Stein equation is solved then the real Schur decomposition  must fullfill: \f$ A = Q_A \hat{A} {Q_A}^H \f$.
 * If a generalized Lyapunov or Stein equation is solved then the real Schur decomposition must fullfill: \f$ A = Q_A \hat{A} {Q_E}^H \f$ and  \f$ E = Q_A \hat{E} {Q_E}^H \f$.
 *
 *
 * @attention Internal use only.
 *
 */
int mess_tglyap3x(_glyap3_eqn_t eqn_t, mess_operation_t op, mess_matrix Ahat, mess_matrix QA, mess_matrix Ehat, mess_matrix QE, mess_matrix Y,  mess_matrix X){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_matrix eye;
    mess_direct solver;
    _glyap3 * data = NULL;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_operation_type(op);
    mess_check_nullpointer(Ahat);
    mess_check_dense(Ahat);
    mess_check_real(Ahat);
    mess_check_square(Ahat);

    mess_check_nullpointer(QA);
    mess_check_dense(QA);
    mess_check_real(QA);
    mess_check_square(QA);

    mess_check_nullpointer(Y);
    mess_check_dense(Y);
    mess_check_real(Y);
    mess_check_square(Y);

    mess_check_nullpointer(X);

    mess_check_same_size(Ahat,QA);
    mess_check_same_size(Ahat,Y);


    if(eqn_t == GLYAP3_GEN_LYAP || eqn_t == GLYAP3_GEN_STEIN){

        mess_check_dense(Ehat);
        mess_check_real(Ehat);
        mess_check_square(Ehat);

        mess_check_nullpointer(QE);
        mess_check_dense(QE);
        mess_check_real(QE);
        mess_check_square(QE);

        mess_check_same_size(Ahat,Ehat);
        mess_check_same_size(Ahat,QE);
    }


    /*-----------------------------------------------------------------------------
     *  create solver from schur decomposition
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&eye);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_direct_init(&solver);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
    mess_try_alloc(data, _glyap3 *, sizeof(_glyap3));
    data->eqn_type      = eqn_t;

    data->Ahat          = Ahat;
    data->QA            = QA;

    if(eqn_t == GLYAP3_GEN_LYAP || eqn_t == GLYAP3_GEN_STEIN){
        data->Ehat          = Ehat;
        data->QE            = QE;
    }

    //please not that dgelyp cannot handle standard Stein equations, therefore we have to generate an identity matrix
    if(eqn_t == GLYAP3_STD_STEIN){
        ret = mess_matrix_eye(eye,Ahat->rows,Ahat->cols,MESS_DENSE);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
        data->Ehat          = eye;
        data->QE            = QA;
        data->eqn_type      = GLYAP3_GEN_STEIN;
    }

    data->semidefinite  = 0;

    solver->rows = Ahat->rows;
    solver->cols = Ahat->cols;
    SET_SOLVERNAME(solver->name, glyap3_eqn_t_str(eqn_t));
    solver->data    = (void *)data;
    solver->solvem  = glyap3_solvem;
    solver->solvemt = glyap3_solvemt;
    solver->solvemh = glyap3_solvemt;
    solver->clear   = glyap3_clear;
    solver->data_type = MESS_REAL;

    /*-----------------------------------------------------------------------------
     *  solve lyapunov equation
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_solvem(op, solver, Y, X);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);

    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    data->Ahat  = NULL;
    data->QA    = NULL;
    data->Ehat  = NULL;
    data->QE    = NULL;

    ret = mess_direct_clear(&solver);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_clear);
    ret = mess_matrix_clear(&eye);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_clear);
    return 0;
}




/**
 * @brief Solve a standard / generalized Lyapunov  equation.
 * @param[in] op        input determine which equation has to be solved
 * @param[in] A         input real dense matrix \f$A\f$
 * @param[in] E         input real dense matrix \f$E\f$, (optional @c NULL if no generalized equation)
 * @param[in] Y         input real dense matrix symmetric \f$Y\f$
 * @param[out] Ahat     output real dense matrix Schur Form of \f$A\f$, (optional @c NULL if not wanted)
 * @param[out] QA       output real dense matrix Schur transformation matrix of \f$A\f$, (optional @c NULL if not wanted)
 * @param[out] Ehat     output real dense matrix Schur Form of \f$E\f$, (optional @c NULL if not wanted)
 * @param[out] QE       output real dense matrix Schur transformation matrix of \f$E\f$, (optional @c NULL if not wanted)
 * @param[out] X        output real dense symmetrix matrix, solution of equation \f$ X\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_glyap function solves the standard / generalized Lyapunov equation.
 *
 * Depending on @p E and @p op the following equation is solved:
 *
 * <ul>
 *  <li> @p E is @c NULL and @ref MESS_OP_NONE \f$ A X + XA^ T + Y = 0\f$.
 *  <li> @p E is not @c NULL and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X + XA + Y = 0\f$.
 *  <li> @p E is @c NULL  and @ref MESS_OP_NONE \f$ A X E^T + E X A^ T + Y = 0\f$.
 *  <li> @p E is not @c NULL and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X E + E^T X A + Y = 0\f$.
 * </ul>
 *
 * The right hand side matrix @p Y must be real and symmetric.
 *
 * In order to solve the equations a (generalized) Schur decomposition of \f$ A \f$ / \f$ (A,E) \f$ is computed internally.
 *
 * If a standard equation  is solved then the real Schur decomposition fullfills: \f$ A = Q_A \hat{A} {Q_A}^H \f$.
 * If a generalized equation is solved then the real Schur decomposition fullfills: \f$ A = Q_A \hat{A} {Q_E}^H \f$ and  \f$ E = Q_A \hat{E} {Q_E}^H \f$.
 *
 * If @p Ahat, @p QA, @p Ehat or @p QE is not @c NULL then the corresponding matrix of the Schur decomposition is written to that argument.
 *
 */
int mess_glyap(mess_operation_t op, mess_matrix A, mess_matrix E, mess_matrix Y, mess_matrix Ahat, mess_matrix QA, mess_matrix Ehat, mess_matrix QE, mess_matrix X){
    return mess_glyap3x(E? GLYAP3_GEN_LYAP : GLYAP3_STD_LYAP, op, A, E, Y, Ahat, QA, Ehat, QE, X);
}


/**
 * @brief Solve a standard / generalized Stein  equation.
 * @param[in] op        input determine which equation has to be solved
 * @param[in] A         input real dense matrix \f$A\f$
 * @param[in] E         input real dense matrix \f$E\f$, (optional @c NULL if no generalized equation)
 * @param[in] Y         input real dense matrix symmetric \f$Y\f$
 * @param[out] Ahat     output real dense matrix Schur Form of \f$A\f$, (optional @c NULL if not wanted)
 * @param[out] QA       output real dense matrix Schur transformation matrix of \f$A\f$, (optional @c NULL if not wanted)
 * @param[out] Ehat     output real dense matrix Schur Form of \f$E\f$, (optional @c NULL if not wanted)
 * @param[out] QE       output real dense matrix Schur transformation matrix of \f$E\f$, (optional @c NULL if not wanted)
 * @param[out] X        output real dense symmetrix matrix, solution of equation \f$ X\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_gstein function solves the standard / generalized Stein equation.
 *
 * Depending on @p E and @p op the following equation is solved:
 *
 * <ul>
 *  <li> @p E is @c NULL and @ref MESS_OP_NONE \f$ A X A^T - X + Y = 0\f$.
 *  <li> @p E is @c NULL and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A X A^T - X + Y = 0 \f$.
 *  <li> @p E is not @c NULL and @ref MESS_OP_NONE \f$ A X A^T - E X E^ T + Y = 0\f$.
 *  <li> @p E is not @c NULL and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X A - E^T X E + Y = 0\f$.
 * </ul>
 *
 * The right hand side matrix @p Y must be real and symmetric.
 *
 * In order to solve the equations a (generalized) Schur decomposition of \f$ A \f$ / \f$ (A,E) \f$ is computed internally.
 *
 * If a standard equation is solved then the real Schur decomposition fullfills: \f$ A = Q_A \hat{A} {Q_A}^H \f$.
 * If a generalized equation is solved then the real Schur decomposition fullfills: \f$ A = Q_A \hat{A} {Q_E}^H \f$ and  \f$ E = Q_A \hat{E} {Q_E}^H \f$.
 *
 * If @p Ahat, @p QA, @p Ehat or @p QE is not @c NULL then the corresponding matrix of the Schur decomposition is written to that argument.
 *
 */
int mess_gstein(mess_operation_t op, mess_matrix A, mess_matrix E, mess_matrix Y, mess_matrix Ahat, mess_matrix QA, mess_matrix Ehat, mess_matrix QE, mess_matrix X){
    return mess_glyap3x(E? GLYAP3_GEN_STEIN : GLYAP3_STD_STEIN, op, A, E, Y, Ahat, QA, Ehat, QE, X);
}


/**
 * @brief Solve a standard / generalized Lyapunov  equation from given Schur decomposition.
 * @param[in] op        input determine which equation has to be solved
 * @param[in] Ahat      input real dense matrix Schur Form of \f$A\f$
 * @param[in] QA        input real dense matrix Schur transformation matrix of \f$A\f$
 * @param[in] Ehat      input real dense matrix Schur Form of \f$E\f$, (only necessary for generalized equation)
 * @param[in] QE        input real dense matrix Schur transformation matrix of \f$E\f$, (only necessary for generalized equation)
 * @param[in] Y         input real dense matrix symmetric \f$Y\f$
 * @param[out] X        output real dense symmetrix matrix, solution of equation \f$ X\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_tglyap function solves the standard / generalized Lyapunov  equation.
 *
 * Depending on @p Ehat and @p op the following equation is solved:
 *
 * <ul>
 *  <li> @p Ehat is @c NULL and @ref MESS_OP_NONE \f$ A X + XA^ T + Y = 0\f$.
 *  <li> @p Ehat is @c NULL and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X + XA + Y = 0\f$.
 *  <li> @p Ehat is not @c NULL and @ref MESS_OP_NONE \f$ A X E^T + E X A^ T + Y = 0\f$.
 *  <li> @p Ehat is not @c NULL and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X E + E^T X A + Y = 0\f$.
 * </ul>
 *
 * The right hand side matrix @p Y must be real and symmetric.
 *
 * In order to solve the equations a (generalized) Schur decomposition of \f$ A \f$ / \f$ (A,E) \f$ is requiered.
 *
 * If a standard equation is solved then the real Schur decomposition  must fullfill: \f$ A = Q_A \hat{A} {Q_A}^H \f$.
 * If a generalized equation is solved then the real Schur decomposition must fullfill: \f$ A = Q_A \hat{A} {Q_E}^H \f$ and  \f$ E = Q_A \hat{E} {Q_E}^H \f$.
 *
 */
int mess_tglyap(mess_operation_t op, mess_matrix Ahat, mess_matrix QA, mess_matrix Ehat, mess_matrix QE, mess_matrix Y,  mess_matrix X){
    return  mess_tglyap3x(Ehat? GLYAP3_GEN_LYAP : GLYAP3_STD_LYAP,op, Ahat, QA, Ehat, QE, Y, X);
}

/**
 * @brief Solve a standard / generalized Stein  equation from given Schur decomposition.
 * @param[in] op        input determine which equation has to be solved
 * @param[in] Ahat      input real dense matrix Schur Form of \f$A\f$
 * @param[in] QA        input real dense matrix Schur transformation matrix of \f$A\f$
 * @param[in] Ehat      input real dense matrix Schur Form of \f$E\f$, (only necessary for generalized equation)
 * @param[in] QE        input real dense matrix Schur transformation matrix of \f$E\f$, (only necessary for generalized equation)
 * @param[in] Y         input real dense matrix symmetric \f$Y\f$
 * @param[out] X        output real dense symmetrix matrix, solution of equation \f$ X\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_tglyap function solves the standard / generalized Stein equation.
 *
 * Depending on @p Ehat and @p op the following equation is solved:
 *
 * <ul>
 *  <li> @p Ehat is @c NULL and @ref MESS_OP_NONE \f$ A X + XA^ T + Y = 0\f$.
 *  <li> @p Ehat is @c NULL and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X + XA + Y = 0\f$.
 *  <li> @p Ehat is not @c NULL and @ref MESS_OP_NONE \f$ A X E^T + E X A^ T + Y = 0\f$.
 *  <li> @p Ehat is not @c NULL and @ref MESS_OP_TRANSPOSE or @ref MESS_OP_HERMITIAN \f$ A^T X E + E^T X A + Y = 0\f$.
 * </ul>
 *
 * The right hand side matrix @p Y must be real and symmetric.
 *
 * In order to solve the equations a (generalized) Schur decomposition of \f$ A \f$ / \f$ (A,E) \f$ is requiered.
 *
 * If a standard  equation is solved then the real Schur decomposition  must fullfill: \f$ A = Q_A \hat{A} {Q_A}^H \f$.
 * If a generalized  equation is solved then the real Schur decomposition must fullfill: \f$ A = Q_A \hat{A} {Q_E}^H \f$ and  \f$ E = Q_A \hat{E} {Q_E}^H \f$.
 *
 */
int mess_tgstein(mess_operation_t op, mess_matrix Ahat, mess_matrix QA, mess_matrix Ehat, mess_matrix QE, mess_matrix Y,  mess_matrix X){
    return  mess_tglyap3x(Ehat? GLYAP3_GEN_STEIN : GLYAP3_STD_STEIN,op, Ahat, QA, Ehat, QE, Y, X);
}












