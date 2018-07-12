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
 * @file lib/direct/singlesolver/cholmod_chol.c
 * @brief Interface for Tim Davis @cholmod (Cholesky Decomposition)
 * @author @mbehr
 *
 * This file provide an interface for @cholmod. You need a working
 * @suitesparse [2] to use it.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#ifdef MESS_HAVE_CHOLMOD
#include "mess/interface_cholmod.h"
#ifndef SuiteSparse_long
#define SuiteSparse_long long
#endif


#ifndef MESS_MATLAB
/**
 * @brief Converts the return value of cholmod_symmetry to human readable string.
 * @param[in] symmetry  input return value to be converted to a string
 * @return A constant pointer to a human readable string containing meaning of the return value of cholmod_symmetry.
 *
 */
static char* symmetry_char(int symmetry){
    switch(symmetry){
    case -1:
        return "Out of Memory, stype not zero, Matrix is not sorted";
    case 0:
        return "Matrix is rectangular";
    case 2:
        return "Matrix is unsymmetric";
    case 3:
        return "Matrix is symmetric, but with non-pos. diagonal";
    case 4:
        return "Matrix is Hermitian, but with non-pos. diagonal";
    case 5:
        return "Matrix is skew symmetric";
    case 6:
        return "Matrix is symmetric with positive diagonal";
    case 7:
        return "Matrix is Hermitian with positive diagonal";
    default:
        return "Unknown return type";
    }
}
#endif

/** Internal structure to save the factorization. */
struct mess_cholmod_solver {
    mess_int_t dim;
    cholmod_common c;
    cholmod_sparse *A;
    cholmod_factor *L;
};

/**
 * @brief Solve Ax=b.
 * @param[in] data  input solver data
 * @param[in] b  input right hand side
 * @param[out] x solution
 *
 * The @ref mess_cholmod_solve function solves the system Ax=b with x and b vectors.
 *
 */
static int mess_cholmod_solve(void*data, mess_vector b, mess_vector x) {
    MSG_FNAME(__func__);
    struct mess_cholmod_solver * sol = (struct mess_cholmod_solver *)data;
    mess_int_t ret=0;
    cholmod_dense *x_chol=NULL,*b_chol;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);

    /*-----------------------------------------------------------------------------
     *  convert input to cholmod
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_convert_dense_to_cholmod(b,&b_chol,&(sol->c));    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_convert_dense_to_cholmod);

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/

    x_chol=cholmod_l_solve(CHOLMOD_A,sol->L,b_chol,&(sol->c));
    if(!x_chol){
        MSG_ERROR("CHOLMOD: Error during cholmod_l_solve.\n");
        return MESS_ERROR_CHOLMOD;
    }

    /*-----------------------------------------------------------------------------
     *  convert solution to  mess
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_convert_cholmod_to_dense(x_chol,x, &(sol->c));    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_convert_cholmod_to_dense);

    /*-----------------------------------------------------------------------------
     *  clear cholmod dense
     *-----------------------------------------------------------------------------*/

    if(!cholmod_l_free_dense(&b_chol,&(sol->c))){
        MSG_ERROR("CHOLMOD: Error during cholmod_l_free_dense.\n");
        return MESS_ERROR_CHOLMOD;
    }

    if(!cholmod_l_free_dense(&x_chol,&(sol->c))){
        MSG_ERROR("CHOLMOD: Error during cholmod_l_free_dense.\n");
        return MESS_ERROR_CHOLMOD;
    }


    return 0;
}

/**
 * @brief Solve Ax=b.
 * @param[in] data  input solver data
 * @param[in] b  input right hand sides
 * @param[out] x solution
 *
 * The @ref mess_cholmod_solvem function solves the system AX=B with X and B matrices.
 *
 */
static int mess_cholmod_solvem(void*data, mess_matrix b, mess_matrix x) {
    MSG_FNAME(__func__);
    struct mess_cholmod_solver * sol = (struct mess_cholmod_solver *)data;
    mess_int_t ret=0;
    cholmod_dense *x_chol=NULL,*b_chol;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    mess_check_dense(b);
    /*-----------------------------------------------------------------------------
     *  convert input to cholmod
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_convert_dense_to_cholmod(b,&b_chol,&(sol->c));    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert_dense_to_cholmod);

    /*-----------------------------------------------------------------------------
     *  solve system
     *-----------------------------------------------------------------------------*/
    x_chol=cholmod_l_solve(CHOLMOD_A,sol->L,b_chol,&(sol->c));
    if(!x_chol){
        MSG_ERROR("CHOLMOD: Error during cholmod_l_solve.\n");
        return MESS_ERROR_CHOLMOD;
    }

    /*-----------------------------------------------------------------------------
     *  convert solution to @ref mess
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_convert_cholmod_to_dense(x_chol,x, &(sol->c));    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_convert_cholmod_to_dense);

    /*-----------------------------------------------------------------------------
     *  clear cholmod dense
    *-----------------------------------------------------------------------------*/
    if(!cholmod_l_free_dense(&b_chol,&(sol->c))){
        MSG_ERROR("CHOLMOD: Error during cholmod_l_free_dense.\n");
        return MESS_ERROR_CHOLMOD;
    }
    if(!cholmod_l_free_dense(&x_chol,&(sol->c))){
        MSG_ERROR("CHOLMOD: Error during cholmod_l_free_dense.\n");
        return MESS_ERROR_CHOLMOD;
    }

    return 0;

}



/**
 * @brief Gets the row permutation.
 * @param[in] data  input solver data
 * @param[out] p output permutation
 *
 * The @ref mess_cholmod_getpermp function gets the row permutation of the decomposition.
 *
 */
static int mess_cholmod_getpermp(void *data, mess_int_t *p){
    MSG_FNAME(__func__);
    struct mess_cholmod_solver *sol = (struct mess_cholmod_solver *) data;
    mess_int_t i=0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(p);

    /*-----------------------------------------------------------------------------
     *  get permutation
     *-----------------------------------------------------------------------------*/
    SuiteSparse_long *perm = (SuiteSparse_long*) sol->L->Perm;
    if(!perm){
        MSG_ERROR("CHOLMOD attribute Perm in cholmod_factor points to NULL.\n");
        return MESS_ERROR_CHOLMOD;
    }

    for(i=0;i<sol->L->n;++i){
        //printf("Perm[%d]=%d\n",j,perm[j]);
        p[i]=perm[i];
    }

    return 0;
}


/**
 * @brief Compute the inverse.
 * @param[in] data   input pointer to the internal data structure
 * @param[out] inv  output inverse
 *
 * The @ref mess_cholmod_inverse function computes \f$A^{-1}\f$.
 *
 */
static int mess_cholmod_inverse ( void *data, mess_matrix inv )
{
    MSG_FNAME(__func__);
    struct mess_cholmod_solver * sol = (struct mess_cholmod_solver*) data;
    mess_matrix eye;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(inv);
    MESS_MATRIX_RESET(inv);

    /*-----------------------------------------------------------------------------
     *  compute dense inverse
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&eye);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_eye(eye, sol->dim, sol->dim, MESS_DENSE);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
    ret = mess_cholmod_solvem(data, eye, inv);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_cholmod_solvem);
    mess_matrix_clear(&eye);
    return 0;
}       /* -----  end of function mess_cholmod_inverse  ----- */


/**
 * @brief Gets the L (the first factor) of a solver.
 * @param[in] data  input solver data
 * @param[out] L    output L factor
 *
 * The @ref mess_cholmod_getL function gets the L factor from the @ref mess_cholmod_solver object. Output is in @ref MESS_CSR
 * format.
 *
 */
static int mess_cholmod_getL(void *data, mess_matrix L ){
    MSG_FNAME(__func__);
    struct mess_cholmod_solver *sol = (struct mess_cholmod_solver *) data;
    mess_matrix tmp =NULL ;
    cholmod_sparse *L_sparse_chol=NULL;
    cholmod_factor *L_factor=NULL;
    MSG_INFO("get L\n");
    int ret;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(L);

    /*-----------------------------------------------------------------------------
     *  get L
     *-----------------------------------------------------------------------------*/

    MESS_MATRIX_RESET(L);
    ret = mess_matrix_init(&tmp);                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

    L_factor = cholmod_l_copy_factor(sol->L, &(sol->c));
    if(!L_factor){
        MSG_ERROR("CHOLMOD:Error during cholmod_l_copy_factor");
        return MESS_ERROR_CHOLMOD;
    }

    L_sparse_chol =cholmod_l_factor_to_sparse(L_factor,&(sol->c));
    if(!L_sparse_chol){
        MSG_ERROR("CHOLMOD:Error during cholmod_l_factor_to_sparse.\n");
        return MESS_ERROR_CHOLMOD;
    }

    ret = mess_matrix_convert_cholmod_to_csc(L_sparse_chol,tmp,&(sol->c));      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert_cholmod_to_csc);
    ret = mess_matrix_convert_csc_csr(tmp,L);                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert_csc_csr);

    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_clear(&tmp);                                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    if(!cholmod_l_free_sparse(&L_sparse_chol,&(sol->c))){
        MSG_ERROR("CHOLMOD:Error during cholmod_l_free_sparse.\n");
        return MESS_ERROR_CHOLMOD;
    }
    if(!cholmod_l_free_factor(&L_factor,&(sol->c))){
            MSG_ERROR("CHOLMOD:Error during cholmod_l_free_factor.\n");
            return MESS_ERROR_CHOLMOD;
    }


    return 0;
}


/**
 * @brief Gets the U (the second factor) of a solver.
 * @param[in] data  input solver data
 * @param[out] U    output U factor
 *
 * The @ref mess_cholmod_getU function gets the U factor from the @ref mess_cholmod_solver object. Output is in @ref MESS_CSR
 * format.
 *
 */
static int mess_cholmod_getU(void *data, mess_matrix U ){
    MSG_FNAME(__func__);
    struct mess_cholmod_solver *sol = (struct mess_cholmod_solver *) data;
    mess_matrix tmp =NULL ;
    cholmod_sparse *L_sparse_chol=NULL;
    cholmod_factor *L_factor=NULL;
    MSG_INFO("get U\n");
    int ret;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(U);

    /*-----------------------------------------------------------------------------
     *  get U
     *-----------------------------------------------------------------------------*/

    MESS_MATRIX_RESET(U);
    ret = mess_matrix_init(&tmp);                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

    L_factor = cholmod_l_copy_factor(sol->L, &(sol->c));
    if(!L_factor){
        MSG_ERROR("CHOLMOD:Error during cholmod-l_copy_factor");
        return MESS_ERROR_CHOLMOD;
    }

    L_sparse_chol = cholmod_l_factor_to_sparse(L_factor,&(sol->c));
    if(!L_sparse_chol){
        MSG_ERROR("CHOLMOD:Error during cholmod_l_factor_to_sparse.\n");
        return MESS_ERROR_CHOLMOD;
    }
    ret = mess_matrix_convert_cholmod_to_csc(L_sparse_chol,tmp,&(sol->c));      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert_cholmod_to_csc);
    ret = mess_matrix_convert_csc_csr(tmp,U);                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert_csc_csr);
    ret = mess_matrix_ctranspose(U,tmp);                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    ret = mess_matrix_copy(tmp,U);                                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_clear(&tmp);                                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    if(!cholmod_l_free_sparse(&L_sparse_chol,&(sol->c))){
        MSG_ERROR("CHOLMOD:Error during cholmod_l_free_sparse.\n");
        return MESS_ERROR_CHOLMOD;
    }

    if(!cholmod_l_free_factor(&L_factor,&(sol->c))){
            MSG_ERROR("CHOLMOD:Error during cholmod_l_free_factor.\n");
            return MESS_ERROR_CHOLMOD;
    }



    return 0;
}

/**
 * @brief Clears a cholmod solver.
 * @param[in,out] solver solver data
 *
 * The @ref mess_cholmod_clear function is the clean up function.
 *
 */
static int mess_cholmod_clear(void *solver){
    MSG_FNAME(__func__);
    struct mess_cholmod_solver * sol = (struct mess_cholmod_solver*) solver;
    if ( sol != NULL) {
        if(!cholmod_l_free_factor(&(sol->L), &(sol->c))){
            MSG_ERROR("CHOLMOD: Error during cholmod_l_free_factor.\n");
            return MESS_ERROR_CHOLMOD;
        }
        if(!cholmod_l_free_sparse(&(sol->A), &(sol->c))){
            MSG_ERROR("CHOLMOD: Error during cholmod_l_free_sparse.\n");
            return MESS_ERROR_CHOLMOD;
        }
        if(!cholmod_l_finish(&(sol->c))){
            MSG_ERROR("CHOLMOD: Error during cholmod_l_finish.\n");
            return MESS_ERROR_CHOLMOD;
        }
        mess_free( sol );
    }
    return 0;
}

/**
 * @brief Generate a @cholmod based solver.
 * @param[in] matrix  input system matrix
 * @param[out] solver output solver
 *
 * The @ref mess_direct_create_cholmod_cholesky function generate a  Cholesky factorization based solver. It uses
 * @cholmod to reorder and factorize the system matrix.
 */
int mess_direct_create_cholmod_cholesky(mess_matrix matrix, mess_direct solver){
    MSG_FNAME(__func__);
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input and convert matrix
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_square(matrix);
    mess_check_real_or_complex(matrix);

    mess_matrix cscmatrix;
    ret = mess_matrix_init(&cscmatrix); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_convert(matrix,cscmatrix,MESS_CSC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert);

    /*-----------------------------------------------------------------------------
     *  alloc @ref mess_cholmod_solver structure
     *-----------------------------------------------------------------------------*/
    struct mess_cholmod_solver * data;

    mess_try_alloc(data, struct mess_cholmod_solver*, sizeof(struct mess_cholmod_solver));
    cholmod_l_start(&(data->c));
    (data->c).itype         = CHOLMOD_LONG;
    (data->c).dtype         = CHOLMOD_DOUBLE;
    (data->c).supernodal    = CHOLMOD_SUPERNODAL;


    /*-----------------------------------------------------------------------------
     *  convert matrix to cholmod and determine symmetry
     *-----------------------------------------------------------------------------*/
    //convert to cholmod
    ret = mess_matrix_convert_csc_to_cholmod(cscmatrix,&(data->A), &(data->c));             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert_csc_to_cholmod);

    //sort matrix
    if(!cholmod_l_sort(data->A,&(data->c))){
        MSG_ERROR("CHOLMOD: Error during cholmod_l_sort.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    //determine type of symmetry but not if we use the matlab interface
#ifndef MESS_MATLAB
    int symmetry = 0;
    SuiteSparse_long xmatched=0,pmatched=0, nzoffdiag=0, nzdiag=0;
    symmetry = cholmod_l_symmetry(data->A,1,&xmatched,&pmatched,&nzoffdiag,&nzdiag,&(data->c));
    if (symmetry <= CHOLMOD_MM_UNSYMMETRIC ){
        //MSG_ERROR("CHOLMOD: Error during cholmod_l_symmetry: %s\nxmatched: %d\npmatched:%d\nnzoff:%d\nnzdiag:%d\n",symmetry_char(symmetry),xmatched,pmatched,nzoffdiag,nzdiag);
        MSG_ERROR("CHOLMOD: Error during cholmod_l_symmetry: %s\n",symmetry_char(symmetry));
        if(!cholmod_l_free_sparse (&(data->A), &(data->c))){
            MSG_ERROR("CHOLMOD:Error during cholmod_l_free_sparse.\n");
            return MESS_ERROR_CHOLMOD;
        }
        if(!cholmod_l_finish (&(data->c))){
            MSG_ERROR("CHOLMOD:ERROR during cholmod_l_finish.\n");
            return MESS_ERROR_CHOLMOD;
        }
        mess_free(data);
        return MESS_ERROR_SYMMETRIC;
    }
    MSG_INFO("symmetry: %s \n",symmetry_char(symmetry));
    data->A->stype=symmetry;
#endif
    /*-----------------------------------------------------------------------------
     *  analyze and factorize matrix A = L*L'
     *-----------------------------------------------------------------------------*/
    data->L = cholmod_l_analyze(data->A,&(data->c));
    if(!cholmod_l_factorize(data->A,data->L,&(data->c))){
        MSG_ERROR("CHOLMOD: Error during cholmod_l_factorize.\n")
                return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  construct the solver structure
     *-----------------------------------------------------------------------------*/
    data->dim           = matrix->rows;
    solver->rows        = matrix->rows;
    solver->cols        = matrix->cols;
    solver->data        = (void *)data;
    solver->solve       = mess_cholmod_solve;
    solver->solvet      = mess_cholmod_solve;
    solver->solveh      = mess_cholmod_solve;
    solver->solvem      = mess_cholmod_solvem;
    solver->solvemt     = mess_cholmod_solvem;
    solver->solvemh     = mess_cholmod_solvem;
    solver->clear       = mess_cholmod_clear;
    solver->getL        = mess_cholmod_getL;
    solver->getU        = mess_cholmod_getU;
    solver->getpermp    = mess_cholmod_getpermp;
    solver->getpermq    = mess_cholmod_getpermp;
    solver->inverse     = mess_cholmod_inverse;
    solver->det         = NULL;
    solver->detc        = NULL;
    solver->data_type   = matrix->data_type;

    //print
    cholmod_l_print_factor(data->L,"L",&(data->c));
    cholmod_l_print_sparse(data->A, "A", &(data->c));
    cholmod_l_print_common("Common",&(data->c));


    SET_SOLVERNAME(solver->name, __func__);

    /*-----------------------------------------------------------------------------
     *  clear additional data
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_clear(&cscmatrix);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);

    return 0;
}

#endif











