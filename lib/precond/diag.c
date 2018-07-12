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
 * @file lib/precond/diag.c
 * @brief Diagonal/Jacobi preconditioner.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>

#ifdef _OPENMP_H
#include <omp.h>
#endif

/**
 * @internal
 * @brief Solves a linear system arises from a diagonal preconditioner (real case).
 * @param[in] myself    input @ref mess_precond diagonal preconditioner
 * @param[in] opt       input @ref mess_solver_options structure
 * @param[in] in        input the right hand side?
 * @param[out] out      output @ref mess_vector structure contains the solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __mess_precond_diag_solve function solves a linear System arises from a diagonal preconditioner (real case).
 *
 * @attention Internal use only.
 *
 * */
static int __mess_precond_diag_solve(mess_precond myself, mess_solver_options opt, mess_vector in, mess_vector out){
    MSG_FNAME(__func__);
    mess_int_t dim;
    mess_int_t i;
    double *data;
    int ret = 0;

    data = (double *)myself->data;
    dim  = in->dim;
    if ( MESS_IS_REAL(in)){
        ret = mess_vector_toreal_nowarn(out); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        for (i = 0; i < dim; i++){
            out->values[i] = in->values[i]*(data[i]);
        }
    } else {
        ret = mess_vector_tocomplex(out); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        for (i = 0; i < dim; i++){
            out->values_cpx[i] = in->values_cpx[i]*(data[i]);
        }

    }
    return 0;
}

/**
 * @internal
 * @brief Solves a linear system arises from a diagonal preconditioner (complex case).
 * @param[in] myself    input @ref mess_precond diagonal preconditioner
 * @param[in] opt       input @ref mess_solver_options structure
 * @param[in] in        input the right hand side?
 * @param[out] out      output @ref mess_vector structure contains the solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __mess_precond_diag_solvecomplex function solves a linear System arises from a diagonal preconditioner (complex case).
 *
 * @attention Internal use only.
 *
 * */
static int __mess_precond_diag_solvecomplex(mess_precond myself, mess_solver_options opt, mess_vector in, mess_vector out){
    MSG_FNAME(__func__);
    mess_int_t dim;
    mess_int_t i;
    int ret = 0 ;
    mess_double_cpx_t *data;


    ret = mess_vector_tocomplex(out);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
    data = (mess_double_cpx_t *)myself->data;
    dim  = in->dim;
    if ( MESS_IS_REAL(in)) {
        for (i = 0; i < dim; i++){
            out->values_cpx[i] = in->values[i]*(data[i]);
        }
    } else {
        for (i = 0; i < dim; i++){
            out->values_cpx[i] = in->values_cpx[i]*(data[i]);
        }
    }
    return 0;
}

/**
 * @internal
 * @brief Clear a @ref mess_precond structure (real case).
 * @param[in,out] myself input/output @ref mess_precond structure to clear
 * @return zero always
 *
 * The @ref __mess_precond_diag_clear function clears a @ref mess_precond structure (real case).
 *
 * @attention Internal use only.
 *
 **/
static int __mess_precond_diag_clear(mess_precond myself){
    double  *data;
    if ( myself->data != NULL ){
        data = (double *)myself->data;
        mess_free(data);
        myself->data = NULL;
    }
    return 0;
}

/**
 * @internal
 * @brief Clear a @ref mess_precond structure (complex case).
 * @param[in,out] myself input/output @ref mess_precond structure to clear
 * @return zero always
 *
 * The @ref __mess_precond_diag_complexclear function clears a @ref mess_precond structure (complex case).
 *
 * @attention Internal use only.
 *
 **/
static int __mess_precond_diag_complexclear(mess_precond myself){
    mess_double_cpx_t *data;
    if ( myself->data != NULL ){
        data = (mess_double_cpx_t *)myself->data;
        mess_free(data);
        myself->data = NULL;
    }
    return 0;
}

/**
 * @brief Setup a diagonal preconditioner.
 * @param[out] pre   generated preconditioner
 * @param[in] mat     input matrix to compute the preconditioner for
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_precond_diag function generates a diagonal preconditioner from
 * a matrix \f$ mat \f$. \n
 * If the matrix is not in @ref MESS_CSR storage format it is converted to
 * it internally. The original matrix is kept as it is in this case.
 *
 */
int mess_precond_diag(mess_precond pre,  mess_matrix mat){
    MSG_FNAME(__func__);
    mess_matrix work;
    int conv = 0;
    mess_int_t i, j, dim;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(pre);
    mess_check_nullpointer(mat);
    mess_check_square(mat);


    MESS_MATRIX_CHECKFORMAT(mat, work, conv, MESS_CSR);

    dim = work->rows;
    if (MESS_IS_REAL(mat)) {
        double *data;
        mess_try_alloc (data, double *, sizeof(double)*dim);
        for ( i = 0; i < dim; i++) data[i] = 1.0;
        for (i = 0 ; i < work->rows; i++){
            for (j=work->rowptr[i]; j < work->rowptr[i+1]; j++){
                if ( work->colptr [j] == i ) {
                    data[i] = 1.0/ work->values[j];
                }
            }
        }
        pre->solve = __mess_precond_diag_solve;
        pre->data  = (void*) data;
        pre->clear = __mess_precond_diag_clear;

    } else {
        mess_double_cpx_t *data;
        mess_try_alloc (data, mess_double_cpx_t *, sizeof(mess_double_cpx_t )*dim);
        for ( i = 0; i < dim; i++) data[i] = 1.0;
        for (i = 0 ; i < work->rows; i++){
            for (j=work->rowptr[i]; j < work->rowptr[i+1]; j++){
                if ( work->colptr [j] == i ) {
                    data[i] = 1.0/ work->values_cpx[j];
                }
            }
        }
        pre->solve = __mess_precond_diag_solvecomplex;
        pre->data  = (void*) data;
        pre->clear = __mess_precond_diag_complexclear;
    }
    pre->type = MESS_PRECOND_DIAG;

    if (conv == 0) {
        mess_matrix_clear(&work);
    }
    return 0;
}

