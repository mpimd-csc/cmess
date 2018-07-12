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
 * @file lib/precond/ilu0.c
 * @brief Generate an \f$ ILU(0) \f$ preconditioner.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


/**
 * @internal
 * @brief Internal structure for \f$ILU(0)\f$ preconditioner.
 * Internal structure for \f$ILU(0)\f$ preconditioner.
 * @attention Interal use only.
 */
typedef struct __ilu0_data {
    mess_matrix mat;
    mess_int_t * uptr;
} ilu0_data;

/**
 * @internal
 * @brief Solve with \f$ ILU(0) \f$ preconditioner.
 * @param[in] ilu  input ilu factorized matrix
 * @param[in] opt  input option from the iterative solver (unused)
 * @param[out] out output solution of the equation
 * @param[in] in   input right hand side
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref __mess_precond_ilu0_solve function solves
 * \f[ LUx=b \f]
 * with forward and backward substitution.
 *
 * @attention Interal use only.
 */
static int __mess_precond_ilu0_solve(mess_precond ilu, mess_solver_options opt, mess_vector in, mess_vector out){
    MSG_FNAME(__func__);
    mess_int_t i, ii, j,dim;
    mess_int_t *uptr;
    ilu0_data *data = (ilu0_data *)ilu->data;
    double *y;
    double t;


    dim = data->mat->rows;
    uptr = data->uptr;
    mess_try_alloc(y, double *, sizeof(double) * dim);

    // Forward
    for ( i = 0; i < dim; i++){
        t = 0.0;
        for ( j = data->mat->rowptr[i]; j < uptr[i]; j++){
            t = t + data->mat->values[j] * y[data->mat->colptr[j]];
        }
        y[i] = in->values[i] - t;
    }

    // Backward

    for ( i = 0; i < dim; i++){
        ii = dim - i - 1;
        out->values[ii] = y[ii];
        for ( j = uptr[ii]+1; j < data->mat->rowptr[ii+1]; j++){
            out->values[ii]  = out->values[ii] - data->mat->values[j]*out->values[data->mat->colptr[j]];
        }
        // the main diagonal is inverted.
        out->values[ii] = out->values[ii] * data->mat->values[uptr[ii]];
    }

    mess_free(y);
    return 0;
}

/**
 * @brief Clear function for ILU(0).
 * @param[in,out] ilu mess_precond object
 * @return always zero
 *
 * The @ref __mess_precond_ilu0_clear function clears the \f$ ILU(0) \f$ preconditioner.
 *
 */
static int __mess_precond_ilu0_clear(mess_precond ilu){
    ilu0_data *data = (ilu0_data *) ilu->data;
    mess_matrix_clear(&(data->mat));
    mess_free(data->uptr);
    mess_free(data);
    return 0;
}



/**
 * @brief Compute an \f$ ILU(0) \f$ preconditioner.
 * @param[in] matrix    input matrix to compute the preconditioner for
 * @param[out] pre      output generated  ILU(0) preconditioner
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_precond_ilu0 function computes an \f$ ILU(0) \f$ preconditioner.
 *
 */
int mess_precond_ilu0(mess_precond pre, mess_matrix matrix ){
    MSG_FNAME(__func__);
    mess_int_t dim;
    mess_int_t *workidx;
    double *luval;
    mess_int_t *uptr;
    mess_matrix work;
    int conv = -1;
    mess_int_t i, k, j1, j2, j, jrow = 0, jj, jw;
    double tl;
    mess_int_t breakdown = 0;
    ilu0_data *data;
    int ret = 0;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(pre);
    mess_check_real(matrix);
    mess_check_square(matrix);


    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    dim = matrix->rows;
    MESS_MATRIX_CHECKFORMAT(matrix, work, conv, MESS_CSR);

    mess_try_alloc(workidx, mess_int_t *, sizeof(mess_int_t)*dim);
    mess_try_alloc( luval, double *, sizeof(double)* matrix->nnz);
    mess_try_alloc(uptr, mess_int_t * , sizeof(mess_int_t) *dim);
    mess_try_alloc(data, ilu0_data *, sizeof(ilu0_data));

    ret = mess_matrix_init(&(data->mat));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

    data->mat->rows = work->rows;
    data->mat->cols = work->cols;
    data->mat->nnz = work->nnz;
    mess_try_alloc(data->mat->rowptr, mess_int_t*,  sizeof(mess_int_t) * (dim +1));
    memcpy(data->mat->rowptr, work->rowptr, (dim+1)*sizeof(mess_int_t));
    mess_try_alloc(data->mat->colptr, mess_int_t*, sizeof(mess_int_t) * work->nnz);
    memcpy(data->mat->colptr, work->colptr, sizeof(mess_int_t) *   work->nnz);
    data->mat->store_type = MESS_CSR;
    data->mat->data_type = MESS_REAL;
    data->mat->symmetry = MESS_GENERAL;

    for ( i = 0; i < work->nnz; i++){
        luval[i] = work->values[i];
    }
    for ( i = 0; i<dim; i++){
        workidx[i] = -1;
    }

    for ( k = 0; k < dim; k++){
        j1 = work->rowptr[k];
        j2 = work->rowptr[k+1]-1;

        for ( j = j1; j <= j2; j++){
            workidx[matrix->colptr[j]]= j;
        }
        j=j1;
        do {
            jrow=work->colptr[j];
            if ( jrow >= k )    break;

            tl = luval[j] * luval[uptr[jrow]];
            luval[j] = tl;
            for ( jj = uptr[jrow]+1; jj < work->rowptr[jrow+1]; jj++){
                jw = workidx[matrix->colptr[jj]];
                if ( jw != -1 ) luval[jw] = luval[jw] - tl * luval[jj];
            }
            j++;
        } while ( j <= j2 );

        uptr[k] = j;

        if ( (jrow != k) || ( luval[j] == 0.0)){
            breakdown = -k;
            break;
        }

        luval[j] = 1.0/luval[j];

        for ( j = j1; j <= j2; j++){
            workidx[matrix->colptr[j]]= -1;
        }
    }

    data->mat->values = luval;
    data->uptr = uptr;
    pre->clear = __mess_precond_ilu0_clear;
    pre->solve = __mess_precond_ilu0_solve;
    pre->data = (void* )data;
    pre->type = MESS_PRECOND_ILU0;

    mess_free(workidx);

    if (conv == 0)
        mess_matrix_clear(&work);
    if ( breakdown ) {
        MSG_ERROR("ILU-0 breakdown in line " MESS_PRINTF_INT "\n", -breakdown);
    }

    return breakdown;

}/** \}@ */

