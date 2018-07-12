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
 * @file lib/direct/interface/multidirect.c
 * @brief Generic interface for the multi-LU / multisolvers.
 * @author @koehlerm
 *
 * This file implements the interface to the multidirect solvers.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <assert.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"

#ifdef MESS_DEBUG
#define DTINIT(X) unsigned short X##_dt = (X)->data_type;
#define DTCHECK(X) assert(X##_dt == (X)->data_type);
#else
#define DTINIT(X)
#define DTCHECK(X)
#endif

static mess_multidirect_t selected_solver = MESS_MULTIDIRECT_SPARSE_LU;
static pthread_mutex_t selected_solver_mutex = PTHREAD_MUTEX_INITIALIZER;



/**
 * @brief   Initialize a mess_multidirect object.
 * @param[in,out]   direct  pointer to the mess_multidirect object
 * @return always zero
 *
 * The @ref mess_multidirect_init function initializes a mess_multidirect object. \n
 * If no memory is left, @ref MESS_ERROR_MEMORY is returned.
 *
 */
int mess_multidirect_init( mess_multidirect *direct){
    MSG_FNAME(__func__);
    mess_try_alloc((*direct), struct mess_multidirect_st *, sizeof(struct mess_multidirect_st));
    (*direct)->rows = 0;
    (*direct)->cols = 0;
    (*direct)->indx = 0;
    (*direct)->data = NULL;
    (*direct)->name = NULL;
    (*direct)->solve = NULL;
    (*direct)->solvem = NULL;
    (*direct)->solvemt = NULL;
    (*direct)->solvet = NULL;
    (*direct)->solveh = NULL;
    (*direct)->solvemh = NULL;
    (*direct)->clear = NULL;
    (*direct)->memsize = NULL;
    (*direct)->getdatatype = NULL;
    (*direct)->getL = NULL ;
    (*direct)->getU = NULL ;
    (*direct)->getp = NULL ;
    (*direct)->getq = NULL ;
    return 0;
}

/**
 * @brief Clear a mess_multidirect object.
 * @param[in] direct  input pointer to the object
 * @return always zero
 *
 * The @ref mess_multidirect_clear function clears a @ref mess_multidirect object and
 * call the clear function of the solver.
 *
 */
int mess_multidirect_clear( mess_multidirect *direct){
    MSG_FNAME(__func__);
    if ( *direct == NULL) {
        MSG_INFO(" already cleared\n");
        return 0;
    }
    if ( (*direct)->name != NULL) mess_free( (*direct)->name);
    if  ((*direct)->clear != NULL){
        (*direct)->clear((*direct)->data);
    }
    (*direct)->data = NULL;
    (*direct)->solvem = NULL;
    (*direct)->solvet = NULL;
    (*direct)->solve = NULL;
    (*direct)->solvemt = NULL;
    (*direct)->solveh = NULL;
    (*direct)->solvemh = NULL;
    mess_free((*direct));
    (*direct) = NULL;
    return 0;
}

/**
 * @brief Solve \f$ op(A(i))x = b \f$.
 * @param[in] op      input operation applied to matrix
 * @param[in] solver  input multisolver
 * @param[in] ind     input index of the desired matrix
 * @param[in] b       input right hand side vector
 * @param[out] x     solution vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_multidirect_solve function solves the system
 * \f[ op(A(ind))x=b \f]
 * with \f$ x \f$ and \f$ b \f$ vectors, where \f$ op \f$ is a \ref mess_operation_t, i.e. \f$ op \f$ can be
 * <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 *  |@ref MESS_OP_HERMITIAN      |   \f$ op(A)=A^H\f$         |
 * </center>
 *
 */
int mess_multidirect_solve(mess_operation_t op, mess_multidirect solver, mess_int_t ind, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    if ( ind < 0 || ind >= solver->indx) {
        MSG_ERROR("index " MESS_PRINTF_INT " is out of range\n.", ind);
        return MESS_ERROR_ARGUMENTS;
    }

    DTINIT(b);
    switch( op) {
        case MESS_OP_NONE:
            if ( solver->solve == NULL) {
                MSG_ERROR("Solver %s don't provide a solve function for Ax=b\n", solver->name);
                return MESS_ERROR_MISSING;
            }
            ret = solver->solve(solver->data, ind, b, x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), solver->solve);
            break;
        case MESS_OP_TRANSPOSE:
            if ( solver->solvet == NULL) {
                MSG_ERROR("Solver %s don't provide a solve function for A^Tx=b\n", solver->name);
                return MESS_ERROR_MISSING;
            }
            ret = solver->solvet(solver->data, ind, b, x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), solver->solvet);
            break;
        case MESS_OP_HERMITIAN:
            if ( solver->solveh == NULL) {
                MSG_ERROR("Solver %s don't provide a solve function for A^Hx=b\n", solver->name);
                return MESS_ERROR_MISSING;
            }
            ret = solver->solveh(solver->data, ind, b, x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), solver->solveh);
            break;
    }
    DTCHECK(b);
    return 0;
}

/**
 * @brief Solve \f$ op(A(i))X = B \f$.
 * @param[in] op      input operation applied to matrix
 * @param[in] solver  input multisolver
 * @param[in] ind     input index of the desired matrix
 * @param[in] b       input right hand side matrix
 * @param[out] x     solution matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_multidirect_solvem function solves the system
 * \f[ op(A(ind))X=B \f]
 * with \f$ X \f$ and  \f$ B \f$ matrices, where
 * \f$ op \f$ is a \ref mess_operation_t, i.e. \f$ op \f$ can be
 * <ul>
 * <li> \ref MESS_OP_NONE (\f$ op(A)=A \f$),
 * <li> \ref MESS_OP_TRANSPOSE (\f$ op(A)=A^T \f$),
 * <li> \ref MESS_OP_HERMITIAN (\f$ op(A)=A^H \f$).
 * </ul>
 *
 */
int mess_multidirect_solvem(mess_operation_t op, mess_multidirect solver, mess_int_t ind, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(solver);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( ind < 0 || ind >= solver->indx) {
        MSG_ERROR("index " MESS_PRINTF_INT " is out of range\n.", ind);
        return MESS_ERROR_ARGUMENTS;
    }
    DTINIT(b);

    switch( op) {
        case MESS_OP_NONE:
            if ( solver->solvem == NULL) {
                MSG_ERROR("Solver %s don't provide a solve function for AX=B\n", solver->name);
                return MESS_ERROR_MISSING;
            }
            ret = solver->solvem(solver->data, ind, b, x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), solver->solvem);
            break;
        case MESS_OP_TRANSPOSE:
            if ( solver->solvemt == NULL) {
                MSG_ERROR("Solver %s don't provide a solve function for A^TX=B\n", solver->name);
                return MESS_ERROR_MISSING;
            }
            ret = solver->solvemt(solver->data, ind, b, x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), solver->solvemt);
            break;
        case MESS_OP_HERMITIAN:
            if ( solver->solvemh == NULL) {
                MSG_ERROR("Solver %s don't provide a solve function for A^Hx=B\n", solver->name);
                return MESS_ERROR_MISSING;
            }
            ret = solver->solvemh(solver->data, ind, b, x);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), solver->solvemh);
            break;
    }
    DTCHECK(b);
    return 0;
}



/**
 * @brief Create a multisolver for shifted linear systems \f$ (sA+pI) \f$ or \f$ (sA+pE) \f$.
 * @param[in] matrix  input system matrix
 * @param[in] shiftsl  input left shifts \f$ s \f$, in front of \f$ A \f$, @c NULL is not given
 * @param[in] shiftsr  input right shifts \f$ p \f$, in front of \f$ E \f$
 * @param[out] mlu  multidirect solver to create
 * @param[in] base   input existing solver to get the pattern if it is possible
 * @param[in] shiftmatrix  input matrix to shift if it points to @c NULL, the identity matrix is used
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_multidirect_create function creates a solver for shifted linear systems. Depending on the solver choosen by \ref mess_multidirect_select it
 * call
 * \li \ref mess_multidirect_create_sparse_lu,
 * \li \ref mess_multidirect_create_umfpack,
 *
 * If you specify a base solver you have to ensure that the base solver is created from a matrix with the same pattern
 * as matrix. \n
 * If you do not want to use a base solver set it to @c NULL. \n
 * If you want to shift with the identity matrix set @p shiftmatrix to @c NULL. \n
 * If left or right parameters are not needed set them to @c NULL as well.
 *
 * \sa mess_multidirect_select
 */

int mess_multidirect_create(mess_matrix matrix, mess_vector shiftsl, mess_vector shiftsr , mess_multidirect mlu, mess_direct base, mess_matrix shiftmatrix)
{

    MSG_FNAME(__func__);
    switch(selected_solver) {
#ifdef MESS_HAVE_UMFPACK
        case MESS_MULTIDIRECT_UMFPACK_LU:
            return  mess_multidirect_create_umfpack(matrix, shiftsl, shiftsr, mlu, base, shiftmatrix);
#else
        case MESS_MULTIDIRECT_UMFPACK_LU:
#endif
        case MESS_MULTIDIRECT_SPARSE_LU:
            return mess_multidirect_create_sparse_lu(matrix, shiftsl, shiftsr, mlu, base, shiftmatrix);
        default:
            MSG_ERROR("Solver type unknown.\n");
            return MESS_ERROR_ARGUMENTS;
    }

}


/**
 * @brief Select the multisolver which is used by mess_multidirect_create.
 * @param[in]  solver    input Solver to select
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_multidirect_select function selects one of the following solver to be used as backend in \ref mess_multidirect_create :
 * \li \ref mess_multidirect_create_sparse_lu,
 * \li \ref mess_multidirect_create_umfpack,
 *
 *
 * \sa mess_multidirect_create
 */

int mess_multidirect_select (mess_multidirect_t solver) {
    MSG_FNAME(__func__);
    switch(solver){
        case MESS_MULTIDIRECT_SPARSE_LU:
            break;
        case MESS_MULTIDIRECT_UMFPACK_LU:
        #ifndef MESS_HAVE_UMFPACK
            MSG_ERROR("Solver not available");
            return MESS_ERROR_NOSUPPORT;
        #endif
            break;
        default:
            MSG_ERROR("Solver not available");
            return MESS_ERROR_NOSUPPORT;
    }
    pthread_mutex_lock (&selected_solver_mutex);
    selected_solver = solver;
    pthread_mutex_unlock(&selected_solver_mutex);

    return 0;
}



