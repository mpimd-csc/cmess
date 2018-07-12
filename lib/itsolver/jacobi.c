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
 * @file lib/itsolver/jacobi.c
 * @brief Solve a linear equation with Jacobi iteration.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mess/mess.h"
#include "mess/error_macro.h"

/**
 * @brief Solve a system of linear equations with the Jacobi method.
 * @param[in] matrix     input matrix \f$ A \f$ of the equation
 * @param[in] pre    input preconditioner (ignored for this method)
 * @param[in] b               input right hand side
 * @param[in,out] x     on input: intial solution,\n on output: computed solution
 * @param[in] opt    input options for the solver \n (details see \ref mess_solver_options)
 * @param[out] stat statistics about the solution process \n (details see \ref mess_solver_status)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_solver_jacobi function implements the Jacobi method to solve a system of linear equations
 * \f[ Ax=b. \f]
 * Details about the algorithm can be found in \cite Mei11 . \n
 * Additionally the convergence can be checked before using the \ref mess_solver_jacobi_convergence
 * function.
 *
 */
int mess_solver_jacobi( mess_matrix matrix, mess_precond pre, mess_vector b, mess_vector x,
        mess_solver_options opt, mess_solver_status stat)
{

    MSG_FNAME(__func__);
    mess_int_t i, j, it, dim, maxit;
    int ret=0;
    double tol, tolb, normb, normr = 0.0,resr;
    double *r;
    double diag;
    int stepdebug = 0;
    mess_matrix work;
    int conv = -1 ;
    // int ret  = 0;


    /*-----------------------------------------------------------------------------
     *  check inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(stat);
    mess_check_square(matrix);
    mess_check_real(matrix);
    mess_check_real(b);
    mess_check_real(x);

    if ( pre ) {
        MSG_INFO("Preconditioner is ignored.\n");
    }

    MESS_MATRIX_CHECKFORMAT(matrix, work, conv, MESS_CSR);
    if ( opt->stepdebug != NULL ) stepdebug = 1;

    if ( x->dim != matrix->rows ) {
        MSG_ERROR("Initial solution vector has the wrong dimension. x->dim = " MESS_PRINTF_INT " A->rows =" MESS_PRINTF_INT " \n", x->dim, matrix->rows);
        return MESS_ERROR_DIMENSION;
    }
    if ( b->dim != matrix->rows ) {
        MSG_ERROR("Right Hand Side vector has the wrong dimension. b->dim = " MESS_PRINTF_INT " A->rows =" MESS_PRINTF_INT " \n", b->dim, matrix->rows);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  main loop
     *-----------------------------------------------------------------------------*/
    dim = matrix->rows;
    tol = opt->tol;
    maxit = opt->maxit;
    it = 0;

    ret = mess_vector_norm2(b,&normb);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    tolb = tol* normb;

    mess_direct_res2(MESS_OP_NONE, work, x, b, &normr, &resr);
    mess_try_alloc(r, double *, sizeof(double) *dim);
    for (i = 0; i <dim; i++) {
        r[i] = 0;
    }
    while ( it < maxit && normr > tolb ){
        normr = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(j,i,diag) reduction(+:normr) default(shared)
#endif
        for ( i = 0; i < dim; i++){
            r[i] = 0.0;
            diag = 0.0;
            for ( j = work->rowptr[i]; j < work->rowptr[i+1]; j++){
                if ( matrix->colptr[j] == i ) {
                    // fMSG_PRINT( stderr, "Found diag in line %lu\n", i);
                    diag =work->values[j];
                }
                r[i] += work->values[j] * x->values [work->colptr[j]];

            }
            r[i] =  b->values[i] - r[i];
            normr += r[i]*r[i];

            if( diag != 0.0) r[i] = r[i]/diag;
        }

        for ( i = 0; i < dim; i++){
            x->values [i] += r[i];
        }
        normr = sqrt(normr);

        if ( stepdebug ){
            opt->stepdebug(normr, normr/normb, it, opt->aux_stepdebug);
        }
        it++;
    }

    stat->it = it;
    stat->res = normr;
    stat->relres = normr/normb;
    stat->need_restart = ( normr < tolb) ? 0 :1;
    stat->converged = ( stat->relres <= tol) ? 1:0;


    //mess_free( dx);
    mess_free( r);
    if ( stat == 0) {
        mess_matrix_clear(&work);
    }
    return 0;
}

/**
 * @brief Check if the Jacobi iteration is guaranteed to converge.
 * @param[in] matrix   input matrix \f$ A \f$ of the equation
 * @param[out] convergence output parameter \n (true if the Jacobi iteration is guaranteed to converge)
 *  @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_solver_jacobi_convergence function checks if the matrix is strictly diagonal dominant. If the
 * matrix is strictly diagonal dominant the Jacobi iteration is guaranteed to converge. In this case convergence
 * is set to true otherwise to false. \n
 * A proof of this convergence criteria can be found in \cite Mei11 Theorem 4.11.
 */
int mess_solver_jacobi_convergence(mess_matrix matrix, int *convergence){
    MSG_FNAME(__func__);
    mess_int_t i,j;
    double sum;
    double globmax = 0.0;
    double diag;
    int conv = -1;
    mess_matrix work;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(convergence);
    mess_check_square(matrix);
    mess_check_real(matrix);

    MESS_MATRIX_CHECKFORMAT(matrix, work, conv, MESS_CSR);

    for ( i = 0; i < matrix->rows; i++){
        sum = 0.0;
        diag = 0.0;
        for ( j = work->rowptr[i]; j < work->rowptr[i+1]; j++){

            if ( work->colptr[j] != i ){
                sum += fabs(work->values[j]);
            }else{
                diag = fabs(work->values[j]);
            }
        }
        if (diag == 0.0) return 0;
        sum = sum / diag;

        if ( sum > globmax) {
            globmax = sum;
        }
    }

    if ( conv == 0 ){
        mess_matrix_clear(&work);
    }

    if ( globmax < 1.0 ) *convergence = -1; else *convergence = 0;
    return 0;
}

