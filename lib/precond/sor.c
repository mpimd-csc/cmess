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
 * @file lib/precond/sor.c
 * @brief Generate a SOR based preconditioner.
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

struct __presordata {
    mess_solver_options opt;
    mess_solver_status stat;
    mess_matrix matrix;
};

static int __mess_precond_sor_solve(mess_precond myself, mess_solver_options opt, mess_vector b, mess_vector x){
    struct __presordata *data = (struct __presordata *) myself->data;
    return mess_solver_sor(data->matrix, NULL,  b, x, data->opt, data->stat);
}

static int __mess_precond_ssor_solve(mess_precond myself, mess_solver_options opt, mess_vector b, mess_vector x){
    struct __presordata *data = (struct __presordata *) myself->data;
    return mess_solver_ssor(data->matrix, NULL,  b, x, data->opt, data->stat);
}


static int __mess_precond_sor_clear(mess_precond myself) {
    struct __presordata *data = (struct __presordata *) myself->data;
    mess_solver_options_clear(&(data->opt));
    mess_matrix_clear(&(data->matrix));
    mess_solver_status_clear(&(data->stat));
    mess_free(data);
    return 0;
}


/**
 * @brief Compute a SOR preconditioner.
 * @param[out] pre generated preconditioner
 * @param[in] matrix   input matrix to compute the preconditioner for
 * @param[in] it     input number of steps
 * @param[in] omega  input SOR parameter
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_precond_sor function generates a SOR preconditioner on top of the
 * \ref mess_solver_sor function. \n
 * If a symmetric preconditioner is wanted use \ref mess_solver_ssor instead.
 *
 */
int mess_precond_sor ( mess_precond pre, mess_matrix matrix, mess_int_t it, double omega )
{
    MSG_FNAME(__func__);
    struct __presordata * data;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(pre);
    mess_check_nullpointer(matrix);
    mess_check_square(matrix);
    mess_check_real(matrix);
    if ( omega <= 0.0 || omega >=2.0) {
        MSG_ERROR("omega have to be between 0 and 2\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( it < 0) {
        MSG_ERROR("it have to be positive\n");
        return MESS_ERROR_ARGUMENTS;
    }

    mess_try_alloc(data, struct  __presordata*, sizeof(struct __presordata));
    ret = mess_matrix_init(&(data->matrix));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_solver_options_init(&(data->opt));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_options_init);
    ret = mess_solver_status_init(&(data->stat));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_status_init);
    if ( MESS_IS_CSR(matrix)){
        ret  = mess_matrix_copy(matrix, data->matrix);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    } else {
        ret = mess_matrix_convert(matrix, data->matrix, MESS_CSR);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    }

    data->opt->maxit = it;
    data->opt->omega = omega;
    data->opt->tol = mess_eps();

    pre->data = (void * )data;
    pre->solve = __mess_precond_sor_solve;
    pre->clear = __mess_precond_sor_clear;

    return 0;
}

/**
 * @brief Compute a SSOR preconditioner.
 * @param[out] pre generated preconditioner
 * @param[in] matrix   input matrix to compute the preconditioner for
 * @param[in] it     input number of steps
 * @param[in] omega  input SSOR parameter
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_precond_ssor function generates a SSOR preconditioner on top of the
 * \ref mess_solver_ssor function.\n
 * The non symmetric case is implemented in \ref mess_precond_sor.
 *
 */
int mess_precond_ssor ( mess_precond pre, mess_matrix matrix, mess_int_t it, double omega )
{
    MSG_FNAME(__func__);
    struct __presordata * data;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(pre);
    mess_check_nullpointer(matrix);
    mess_check_square(matrix);
    mess_check_real(matrix);
    if ( omega <= 0.0 || omega >=2.0) {
        MSG_ERROR("omega have to be between 0 and 2\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( it < 0) {
        MSG_ERROR("it have to be positive\n");
        return MESS_ERROR_ARGUMENTS;
    }

    mess_try_alloc(data, struct  __presordata*, sizeof(struct __presordata));
    ret = mess_matrix_init(&(data->matrix));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_solver_options_init(&(data->opt));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_options_init);
    ret = mess_solver_status_init(&(data->stat));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_solver_status_init);
    if ( MESS_IS_CSR(matrix)){
        ret  = mess_matrix_copy(matrix, data->matrix);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    } else {
        ret = mess_matrix_convert(matrix, data->matrix, MESS_CSR);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    }

    data->opt->maxit = it;
    data->opt->omega = omega;
    data->opt->tol = mess_eps();

    pre->data = (void * )data;
    pre->solve = __mess_precond_ssor_solve;
    pre->clear = __mess_precond_sor_clear;

    return 0;
}

