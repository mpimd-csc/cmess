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
 * @file lib/direct/interface/direct.c
 * @brief   The ultimate direct solver interface.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <string.h>
#include <pthread.h>
#include "mess/mess.h"
#include "mess/error_macro.h"

/**
 * @brief Initialize a mess_direct object.
 * @param[in,out] direct pointer to the mess_direct object
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_init function initializes a @ref mess_direct object.
 *
 */
int mess_direct_init( mess_direct *direct){
    MSG_FNAME(__func__);
    mess_try_alloc((*direct), struct mess_direct_st *, sizeof(struct mess_direct_st));
    memset(*direct, 0, sizeof(struct mess_direct_st));
    (*direct)->rows = 0;
    (*direct)->cols = 0;
    (*direct)->data = NULL;
    (*direct)->name = NULL;
    (*direct)->solve = NULL;
    (*direct)->solvem = NULL;
    (*direct)->solvemt = NULL;
    (*direct)->solvet = NULL;
    (*direct)->solvemh = NULL;
    (*direct)->solveh = NULL;
    (*direct)->det = NULL;
    (*direct)->detc = NULL;
    (*direct)->getU = NULL;
    (*direct)->getL = NULL;
    (*direct)->getpermp = NULL;
    (*direct)->getpermq = NULL;
    (*direct)->getscalerow = NULL;
    (*direct)->getscalecol = NULL;
    (*direct)->inverse = NULL;
    (*direct)->clear = NULL;
    return 0;
}

/**
 * @brief Clear a mess_direct object.
 * @param[in,out] direct pointer to the mess_direct object
 * @return always zero.
 *
 * The @ref mess_direct_clear function clears a @ref mess_direct object and calls
 * the clear function of the solver.
 *
 */
int mess_direct_clear( mess_direct *direct){
    MSG_FNAME(__func__);
    if ( *direct == NULL) {
        MSG_WARN(" already cleared\n");
        return 0;
    }
    if ( (*direct)->name != NULL) mess_free( (*direct)->name);
    if  ((*direct)->clear != NULL){
        (*direct)->clear((*direct)->data);
    }
    (*direct)->rows = 0;
    (*direct)->cols = 0;
    (*direct)->data = NULL;
    (*direct)->name = NULL;
    (*direct)->solve = NULL;
    (*direct)->solvem = NULL;
    (*direct)->solvemt = NULL;
    (*direct)->solvet = NULL;
    (*direct)->solvemh = NULL;
    (*direct)->solveh = NULL;
    (*direct)->det = NULL;
    (*direct)->detc = NULL;
    (*direct)->getU = NULL;
    (*direct)->getL = NULL;
    (*direct)->getpermp = NULL;
    (*direct)->getpermq = NULL;
    (*direct)->getscalerow = NULL;
    (*direct)->getscalecol = NULL;
    (*direct)->inverse = NULL;
    (*direct)->clear = NULL;

    mess_free((*direct));
    (*direct) = NULL;
    return 0;
}

/**
 * @brief Solve \f$ op(A)X=B \f$, where \f$ B\f$ and \f$ X\f$ are matrices.
 * @param[in] op      input operation applied to A
 * @param[in] solver  input solver to use
 * @param[in] B       input right-hand-side matrix of the equation
 * @param[out] X     solution matrix of the equation
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_solvem function solves \f$ op(A)X=B \f$. The operation @c op can be
 * <ul>
 *  <li> \ref MESS_OP_NONE (\f$ op(A)=A \f$)
 *  <li> \ref MESS_OP_TRANSPOSE (\f$ op(A)=A^T \f$)
 *  <li> \ref MESS_OP_HERMITIAN (\f$ op(A)=A^H \f$).
 * </ul>
 * The function calls the corresponding function of internal solver.
 * If the solver is real then \ref MESS_OP_TRANSPOSE and \ref MESS_OP_HERMITIAN are
 * the same. If a solver does not provide a corresponding internal solver function
 * for multiple right hand sides or the desired operation an error is returned.
 */
int mess_direct_solvem(mess_operation_t op, mess_direct solver, mess_matrix B, mess_matrix X){
    MSG_FNAME(__func__);
    int ret = 0 ;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(B);
    mess_check_nullpointer(X);

    if ( MESS_IS_REAL(solver) && op == MESS_OP_HERMITIAN) {
        op = MESS_OP_TRANSPOSE;
    }

    switch(op) {
        case MESS_OP_NONE:
            if ( solver->solvem == NULL) {
                MSG_ERROR("Solver %s don't provide a solve function for AX=B\n", solver->name);
                return MESS_ERROR_MISSING;
            }
            ret =  solver->solvem(solver->data, B, X);
            break;
        case MESS_OP_TRANSPOSE:
            if ( solver->solvemt == NULL) {
                MSG_ERROR("solver %s don't provide a solve function for A^TX=B\n", solver->name);
                return MESS_ERROR_MISSING;
            }
            ret =  solver->solvemt(solver->data, B,X);

            break;
        case MESS_OP_HERMITIAN:
            if ( solver->solvemh == NULL) {
                MSG_ERROR("Solver %s don't provide a solve function for A^Hx=b\n", solver->name);
                return MESS_ERROR_MISSING;
            }
            ret =  solver->solvemh(solver->data, B, X);
            break;
        default:
            MSG_ERROR("Operation Type not supported\n");
            return MESS_ERROR_ARGUMENTS;
    }
    return ret;
}

/**
 * @brief Solve \f$ op(A)x=b \f$, where \f$ b \f$ and \f$x\f$ are vectors.
 * @param[in] op      input operation applied to A
 * @param[in] solver  input solver to use
 * @param[in] b       input right-hand-side matrix of the equation
 * @param[out] x     solution matrix of the equation
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_solve function solves \f$ op(A)x=b \f$ . The operation @c op can be
 * <ul>
 * <li> \ref MESS_OP_NONE (\f$ op(A)=A \f$)
 * <li> \ref MESS_OP_TRANSPOSE (\f$ op(A)=A^T \f$)
 * <li> \ref MESS_OP_HERMITIAN (\f$ op(A)=A^H \f$).
 * </ul>
 * The function calls the corresponding function of internal solver.
 * If the solver is real then \ref MESS_OP_TRANSPOSE and \ref MESS_OP_HERMITIAN are
 * the same. If a solver does not provide a corresponding internal solver function
 * for vector right-hand-sides or the desired operation an error is returned.
 */
int mess_direct_solve(mess_operation_t op, mess_direct solver, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    int ret = 0 ;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( MESS_IS_REAL(solver) && op == MESS_OP_HERMITIAN) {
        op = MESS_OP_TRANSPOSE;
    }

    switch(op) {
        case MESS_OP_NONE:
            if ( solver->solve == NULL) {
                MSG_ERROR("Solver %s don't provide a solve function for Ax=b\n", solver->name);
                return MESS_ERROR_MISSING;
            }
            ret =  solver->solve(solver->data, b, x);
            break;
        case MESS_OP_TRANSPOSE:
            if ( solver->solvet == NULL) {
                MSG_ERROR("solver %s don't provide a solve function for a^tx=b\n", solver->name);
                return MESS_ERROR_MISSING;
            }
            ret =  solver->solvet(solver->data, b, x);

            break;
        case MESS_OP_HERMITIAN:
            if ( solver->solveh == NULL) {
                MSG_ERROR("Solver %s don't provide a solve function for A^Hx=b\n", solver->name);
                return MESS_ERROR_MISSING;
            }
            ret =  solver->solveh(solver->data, b, x);
            break;
        default:
            MSG_ERROR("Operation Type not supported\n");
            return MESS_ERROR_ARGUMENTS;
    }
    return ret;
}


/**
 * @brief Calculate the determinant of a matrix.
 * @param[in]  solver   input solver build up from the matrix
 * @param[out] m        output mantissa of determinant
 * @param[out] e        output exponent of determinant
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_determinant function computes the determinant of
 * a matrix. It returns an error if the underlying solver does not
 * provide a det function. Unfortunately most of the solvers does not
 * have such a function. This function should be used for real matrices.
 * To avoid overflow problems the determinant is returned as
 * \$det=m*2^e\$.
 *
 */
int mess_direct_determinant(mess_direct solver, double *m, double *e){
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(m);
    mess_check_nullpointer(e);

    if ( solver->det == NULL) {
        MSG_ERROR("Solver %s don't provide a determinant function\n", solver->name);
        return MESS_ERROR_MISSING;
    }
    ret = solver->det ( solver->data, m, e);

    return ret;
}

/**
 * @brief Calculate the determinant of a complex matrix.
 * @param[in]  solver   input solver build up from the matrix
 * @param[out] mr       output mantissa for real part of determinant
 * @param[out] mi       output mantissa for imaginary part of determinant
 * @param[out] e        output exponent for real part and imaginary part of determinant
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_determinantc function computes the determinant of
 * a matrix. It returns an error if the underlying solver does not
 * provide a det function. Unfortunately most of the solvers does not
 * have such a function. This function works for complex matrices.
 * To avoid overflow problems the determinant is returned as
 * \$det=m_r*2^e+m_i*2^e*I\$.
 *
 */
int mess_direct_determinantc(mess_direct solver, double *mr, double *mi, double *e){
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(mr);
    mess_check_nullpointer(mi);
    mess_check_nullpointer(e);

    if ( solver->detc == NULL) {
        MSG_ERROR("Solver %s don't provide a derterminat function\n", solver->name);
        return MESS_ERROR_MISSING;
    }
    ret = solver->detc ( solver->data, mr,mi,e);

    return ret;
}


/**
 * @brief Retrieve the left matrix of a decomposition.
 * @param[in]  solver    input solver
 * @param[out] L    output left matrix matrix L
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_getL function gets the left matrix
 * of a decomposition. It returns an error if the corresponding
 * function is not available for the solver.
 *
 */
int mess_direct_getL(mess_direct solver, mess_matrix L){
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(L);

    if ( solver->getL == NULL) {
        MSG_ERROR("Solver %s don't provide a getL function\n", solver->name);
        return MESS_ERROR_MISSING;
    }
    ret = solver->getL ( solver->data, L);

    return ret;
}

/**
 * @brief Retrieve the right matrix of a decomposition.
 * @param[in]  solver    input solver
 * @param[out] U    output right matrix matrix L
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_getU function gets the right matrix
 * of a decomposition. It returns an error if the corresponding
 * function is not available for the solver.
 *
 */
int mess_direct_getU(mess_direct solver, mess_matrix U){
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(U);

    if ( solver->getU == NULL) {
        MSG_ERROR("Solver %s don't provide a getU function\n", solver->name);
        return MESS_ERROR_MISSING;
    }
    ret = solver->getU ( solver->data, U);

    return ret;
}

/**
 * @brief Get the row permutation of a decomposition.
 * @param[in]  solver  input solver
 * @param[out] p     output row permutation
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_getpermp function gets the row permutation
 * of a decomposition. The array p must be preallocated.
 *
 */
int mess_direct_getpermp(mess_direct solver, mess_int_t *p){
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(p);


    if ( solver->getpermp== NULL) {
        MSG_ERROR("Solver %s don't provide a getpermp function\n", solver->name);
        return MESS_ERROR_MISSING;
    }
    ret = solver->getpermp ( solver->data, p);

    return ret;
}

/**
 * @brief Get the column permutation of a decomposition.
 * @param[in] solver     input solver
 * @param[out]  q   output column permutation
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_getpermq function gets the column permutation
 * of a decomposition. The array q must be preallocated.
 *
 */
int mess_direct_getpermq(mess_direct solver, mess_int_t *q){
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(q);

    if ( solver->getpermq== NULL) {
        MSG_ERROR("Solver %s don't provide a getpermq function\n", solver->name);
        return MESS_ERROR_MISSING;
    }
    ret = solver->getpermq ( solver->data, q);

    return ret;
}

/**
 * @brief Get the row scaling factors of a decomposition.
 * @param[in] solver     input solver
 * @param[out]  r   output row scaling factors
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_getscalerow function gets the row scaling factors
 * of a decomposition.
 *
 */
int mess_direct_getscalerow(mess_direct solver, mess_vector r){
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(r);

    if ( solver->getscalerow== NULL) {
        MSG_ERROR("Solver %s don't provide a getscalerow function\n", solver->name);
        return MESS_ERROR_MISSING;
    }
    ret = solver->getscalerow ( solver->data, r);

    return ret;
}

/**
 * @brief Get the col scaling factors of a decomposition.
 * @param[in] solver    input solver
 * @param[out]  c       output col scaling factors
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_getscalecol function gets the col scaling factors
 * of a decomposition.
 *
 */
int mess_direct_getscalecol(mess_direct solver, mess_vector c){
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(c);

    if ( solver->getscalerow== NULL) {
        MSG_ERROR("Solver %s don't provide a getscalerow function\n", solver->name);
        return MESS_ERROR_MISSING;
    }
    ret = solver->getscalecol ( solver->data, c);

    return ret;
}

/**
 * @brief Compute the inverse of a matrix.
 * @param[in]  solver   input solver for the matrix
 * @param[out] inv     inverse of the matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_inverse function computes the inverse of a matrix.
 * Depending on the solver this is a very time consuming operation.
 *
 */
int mess_direct_inverse ( mess_direct solver, mess_matrix inv )
{
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(inv);

    if ( solver->inverse== NULL) {
        MSG_ERROR("Solver %s don't provide an inverse function\n", solver->name);
        return MESS_ERROR_MISSING;
    }
    ret = solver->inverse ( solver->data, inv);

    return ret ;
}       /* -----  end of function mess_direct_inverse  ----- */


static mess_direct_lupackage_t lu_type = MESS_DIRECT_DEFAULT_LU;
static pthread_mutex_t lu_type_mutex = PTHREAD_MUTEX_INITIALIZER;

static mess_direct_cholpackage_t chol_type = MESS_DIRECT_DEFAULT_CHOLESKY;
static pthread_mutex_t chol_type_mutex = PTHREAD_MUTEX_INITIALIZER;


/**
 * @brief Compute a LU decomposition of a matrix.
 * @param[in] matrix    input matrix which should be decomposed
 * @param[out] solver   solver which contains a LU decomposition of the input matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_lu function computes the LU decomposition of a given matrix. In the default
 * case it selects the solver which works best for the matrix. This means in case of a sparse
 * matrix it tries to use @umfpack or in case of a dense matrix it uses @lapack. For the sparse
 * case are fallback situations implemented, such that it works even if @umfpack is not available.
 * This default selecting behaviour can be changed using \ref mess_direct_lu_select.
 *
 */
int mess_direct_lu ( mess_matrix matrix, mess_direct solver )
{
    MSG_FNAME(__func__);
    mess_direct_lupackage_t lu;
    int ret ;

    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_real_or_complex(matrix);
    mess_check_square(matrix);

    pthread_mutex_lock(&lu_type_mutex);
    lu = lu_type;
    pthread_mutex_unlock(&lu_type_mutex);

    if ( lu == MESS_DIRECT_DEFAULT_LU ){
        if ( MESS_IS_DENSE(matrix)) {
            ret = mess_direct_create_lapack_lu(matrix, solver);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_lapack_lu);
        } else  {
#if defined(MESS_HAVE_UMFPACK)
            ret = mess_direct_create_umfpack(matrix, solver);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_umfpack);
#elif defined(MESS_HAVE_CSPARSE)
            ret = mess_direct_create_csparse_lu(matrix, solver);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_csparse_lu);
#else
            ret = mess_direct_create_sparse_lu(matrix, solver);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sparse_lu);
#endif
        }
    } else {
        switch(lu) {
            case MESS_DIRECT_SPARSE_LU:
                ret = mess_direct_create_sparse_lu(matrix, solver);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sparse_lu);
                break;
            case MESS_DIRECT_DEFAULT_LU:
                /*is this meaning correct*/
            case MESS_DIRECT_LAPACK_LU:
                ret = mess_direct_create_lapack_lu(matrix, solver);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_lapack_lu);
                break;
            case MESS_DIRECT_UMFPACK_LU:
#ifdef MESS_HAVE_UMFPACK
                ret = mess_direct_create_umfpack(matrix, solver);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_umfpack);
                break;
#else
                return MESS_ERROR_NOSUPPORT;
#endif
            case MESS_DIRECT_CSPARSE_LU:
#ifdef  MESS_HAVE_CSPARSE
                ret = mess_direct_create_csparse_lu(matrix, solver);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_csparse_lu);
                break;
#else
                return MESS_ERROR_NOSUPPORT;
#endif
            case MESS_DIRECT_BANDED_LU:
                ret = mess_direct_create_banded(matrix,solver);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_create_banded);
                break;
#ifdef MESS_HAVE_SUPERLU
            case MESS_DIRECT_SUPERLU_LU:
                ret = mess_direct_create_superlu(matrix,solver);
                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_create_superlu);
                break;
#endif
#ifdef MESS_HAVE_MKLPARDISO
            case MESS_DIRECT_MKLPARDISO_LU:
                ret = mess_direct_create_mklpardiso(matrix,solver);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_create_mklpardiso);
                break;
#endif

            default:
                MSG_ERROR("No suitable solver available.");
                return MESS_ERROR_NOSUPPORT;

        }
    }
    return 0;
}       /* -----  end of function mess_direct_lu  ----- */


/**
 * @brief Change the default behaviour of \ref mess_direct_lu
 * @param[in] lu     input type of the LU decomposition used by \ref mess_direct_lu
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_lu_select function selects the solver package which should be used
 * for the \ref mess_direct_lu function. It overrides the automatic selection of the solver.
 * If the automatic selection should be enabled again, use \ref MESS_DIRECT_DEFAULT_LU as @ref mess_direct_lupackage_t.
 *
 */
int mess_direct_lu_select ( mess_direct_lupackage_t lu )
{
    /*-----------------------------------------------------------------------------
     *  check if solver is available
     *-----------------------------------------------------------------------------*/
    MSG_FNAME(__func__);
    switch(lu){
        case MESS_DIRECT_DEFAULT_LU:
            break;
        case MESS_DIRECT_SPARSE_LU:
            break;
        case MESS_DIRECT_LAPACK_LU:
            break;
        case MESS_DIRECT_UMFPACK_LU:
        #ifndef MESS_HAVE_UMFPACK
            MSG_ERROR("Solver not available");
            return MESS_ERROR_NOSUPPORT;
        #endif
            break;
        case MESS_DIRECT_SUPERLU_LU:
        #ifndef MESS_HAVE_SUPERLU
            MSG_ERROR("Solver not available");
            return MESS_ERROR_NOSUPPORT;
        #endif
            break;
        case MESS_DIRECT_CSPARSE_LU:
        #ifndef MESS_HAVE_CSPARSE
            MSG_ERROR("Solver not available");
            return MESS_ERROR_NOSUPPORT;
        #endif
            break;
        case MESS_DIRECT_BANDED_LU:
            break;
        case MESS_DIRECT_MKLPARDISO_LU:
        #ifndef MESS_HAVE_MKLPARDISO
            MSG_ERROR("Solver not available");
            return MESS_ERROR_NOSUPPORT;
        #endif
            break;
        default:
            MSG_ERROR("Solver not available");
            return MESS_ERROR_NOSUPPORT;
    }

    /*-----------------------------------------------------------------------------
     *  set lu_type
     *-----------------------------------------------------------------------------*/
    pthread_mutex_lock (&lu_type_mutex);
    lu_type = lu;
    pthread_mutex_unlock(&lu_type_mutex);
    return 0;
}       /* -----  end of function mess_direct_lu_select  ----- */

/**
 * @brief Compute an Chol decomposition of a matrix.
 * @param[in] matrix    input matrix which should be decomposed
 * @param[out] solver   solver which contains a Cholesky decomposition of the input matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_chol function computes the Cholesky decomposition of a given matrix. In the default
 * case it selects the solver which works best for the matrix. This means in case of a sparse
 * matrix it tries to use @csparse or in case of a dense matrix it uses @lapack. For the sparse
 * case are fallback situations implemented, such that it works even if @csparse is not available.
 * This default selecting behaviour can be changed using \ref mess_direct_chol_select.
 *
 */
int mess_direct_chol ( mess_matrix matrix, mess_direct solver )
{
    MSG_FNAME(__func__);
    mess_direct_cholpackage_t chol;
    int ret ;

    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_real_or_complex(matrix);
    mess_check_square(matrix);

    pthread_mutex_lock(&chol_type_mutex);
    chol  = chol_type;
    pthread_mutex_unlock(&chol_type_mutex);

    if ( chol == MESS_DIRECT_DEFAULT_CHOLESKY ){
        if ( MESS_IS_DENSE(matrix)) {
            ret = mess_direct_create_cholesky(matrix, solver);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_cholesky);
        } else  {
#if defined(MESS_HAVE_CHOLMOD)
            ret = mess_direct_create_cholmod_cholesky(matrix,solver);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_cholmod_cholesky);
#elif defined(MESS_HAVE_CSPARSE)
            ret = mess_direct_create_csparse_cholesky(matrix, solver);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_csparse_cholesky);
#else
            ret = mess_direct_create_cholesky(matrix, solver);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_cholesky);
#endif
        }
    } else {
        switch(chol) {
            case MESS_DIRECT_DEFAULT_CHOLESKY:
            case MESS_DIRECT_LAPACK_CHOLESKY:
                ret = mess_direct_create_cholesky(matrix, solver);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_cholesky);
                break;
            case MESS_DIRECT_CSPARSE_CHOLESKY:
#ifdef  MESS_HAVE_CSPARSE
                ret = mess_direct_create_csparse_cholesky(matrix, solver);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_csparse_cholesky);
                break;
#else
                return MESS_ERROR_NOSUPPORT;
#endif
            case MESS_DIRECT_CHOLMOD_CHOLESKY:
#ifdef MESS_HAVE_CHOLMOD
                ret = mess_direct_create_cholmod_cholesky(matrix,solver);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_cholmod_cholesky);
#else
                return MESS_ERROR_NOSUPPORT;
#endif
#ifdef MESS_HAVE_MKLPARDISO
                ret = mess_direct_create_mklpardiso(matrix,solver);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_create_mklpardiso);
#else
                return MESS_ERROR_NOSUPPORT;
#endif


        }
    }
    return 0;
}       /* -----  end of function mess_direct_chol ----- */



/**
 * @brief Change the default behaviour of \ref mess_direct_chol.
 * @param[in] chol   input type of the Cholesky decomposition used by \ref mess_direct_chol
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_direct_chol_select function selects the solver package which should be used
 * for the \ref mess_direct_chol function. It overrides the automatic selection of the solver.
 * If the automatic selection should be enabled again, use \ref MESS_DIRECT_DEFAULT_CHOLESKY as @ref mess_direct_cholpackage_t.
 *
 */
int mess_direct_chol_select ( mess_direct_cholpackage_t chol )
{
    /*-----------------------------------------------------------------------------
     *  check if solver is available
     *-----------------------------------------------------------------------------*/
    MSG_FNAME(__func__);
    switch(chol){
        case MESS_DIRECT_DEFAULT_CHOLESKY:
            break;
        case MESS_DIRECT_LAPACK_CHOLESKY:
            break;
        case MESS_DIRECT_CSPARSE_CHOLESKY:
        #ifndef MESS_HAVE_CSPARSE
            MSG_ERROR("Solver not available");
            return MESS_ERROR_NOSUPPORT;
        #endif
            break;
        case MESS_DIRECT_CHOLMOD_CHOLESKY:
        #ifndef MESS_HAVE_CHOLMOD
            MSG_ERROR("Solver not available");
            return MESS_ERROR_NOSUPPORT;
        #endif
            break;
        default:
            MSG_ERROR("Solver not available");
            return MESS_ERROR_NOSUPPORT;
    }
    /*-----------------------------------------------------------------------------
     *  set chol_type
     *-----------------------------------------------------------------------------*/
    pthread_mutex_lock (&chol_type_mutex);
    chol_type = chol;
    pthread_mutex_unlock(&chol_type_mutex);
    return 0;
}       /* -----  end of function mess_direct_chol_select  ----- */

