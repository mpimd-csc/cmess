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
 * @file lib/itsolver/status.c
 * @brief  Status / options object for iterative solvers.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


/**
 * @brief Initialize mess_solver_options.
 * @param[in,out] opt pointer to the mess_solver_options object
 * @return always zero
 *
 * The @ref mess_solver_options_init function initializes a mess_solver_options object.
 *
 */
int mess_solver_options_init(mess_solver_options *opt){
    MSG_FNAME(__func__);
    mess_try_alloc(*opt, mess_solver_options, sizeof(struct __mess_solver_options));
    (*opt)->maxit  = 50;
    (*opt)->restarts = 5;
    (*opt)->tol = 1e-7;
    (*opt)->omega = 1.0;
    (*opt)->stepdebug = NULL;
    (*opt)->aux_stepdebug=NULL;
    return 0;
}


/**
 * @brief Clear a mess_solver_options object.
 * @param[in] opt input pointer to the mess_solver_options object
 * @return always zero
 *
 * The @ref mess_solver_options_clear function clears a mess_solver_options object.
 *
 */
int mess_solver_options_clear(mess_solver_options *opt){
    MSG_FNAME(__func__);
    mess_check_nullpointer(*opt);
    mess_free(*opt);
    *opt = NULL;
    return 0;
}


/**
 * @brief Print a mess_solver_options object.
 * @param[in] opt input input mess_solver_options objects
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_solver_options_print function print a mess_solver_options to stdout.
 *
 */
int mess_solver_options_print ( mess_solver_options opt)
{
    MSG_FNAME(__func__);
    mess_check_nullpointer(opt);
    MSG_PRINT("maxit:        " MESS_PRINTF_INT "\n", opt->maxit);
    MSG_PRINT("restarts:     " MESS_PRINTF_INT "\n", opt->restarts);
    MSG_PRINT("relative tol: %.8e\n", opt->tol);
    MSG_PRINT("omega:        %.8e\n", opt->omega);
    return 0;
}       /* -----  end of function mess_solver_options_print  ----- */

/**
 * @brief Initialize a mess_solver_status object.
 * @param[in,out] opt pointer to the mess_solver_status object
 * @return always zero
 *
 * The @ref mess_solver_status_init function initializes a mess_solver_status object.
 *
 */
int mess_solver_status_init(mess_solver_status *opt){
    MSG_FNAME(__func__);
    mess_try_alloc(*opt, mess_solver_status, sizeof(struct __mess_solver_status));
    (*opt)->it = 0 ;
    (*opt)->res = 0.0;
    (*opt)->relres = 0.0;
    (*opt)->need_restart  = 0;
    (*opt)->converged  = 0;
    (*opt)->restarts = 0;
    (*opt)->num_mvp = 0;
    return 0;
}


/**
 * @brief Clear a mess_solver_status object.
 * @param[in] opt input pointer to the mess_solver_status object
 * @return always zero
 *
 * The @ref mess_solver_status_clear function clears a mess_solver_status object.
 *
 */
int mess_solver_status_clear(mess_solver_status *opt){
    MSG_FNAME(__func__);
    mess_check_nullpointer(*opt);
    mess_free(*opt);
    *opt = NULL;
    return 0;
}


/**
 * @brief Print a mess_solver_status object to stdout.
 * @param[in] stat input mess_solver_status object
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_solver_status_print function prints a mess_solver_status object
 * to stdout.
 *
 */
int mess_solver_status_print ( mess_solver_status stat )
{
    MSG_FNAME(__func__);
    mess_check_nullpointer(stat);

    MSG_PRINT("number of iterations: " MESS_PRINTF_INT "\n", stat->it);
    MSG_PRINT("absolute residual:    %.8e\n", stat->res);
    MSG_PRINT("relative residual:    %.8e\n", stat->relres);
    MSG_PRINT("converged:            %d\n", stat->converged);
    MSG_PRINT("need_restart:          %d\n", stat->need_restart);
    MSG_PRINT("restarts:             " MESS_PRINTF_INT "\n", stat->restarts);
    MSG_PRINT("number of mvp's:      %ld\n", (long) stat->num_mvp);
    return 0;
}


/**
 * @brief Reset values in mess_solver_status.
 * @param[in] stat  input mess_solver_status object
 * @return always zero
 *
 * The @ref mess_solver_status_clean function resets all values in mess_solver_status.
 *
 */
int mess_solver_status_clean ( mess_solver_status stat)
{
    MSG_FNAME(__func__);
    mess_check_nullpointer(stat);
    stat->it = 0;
    stat->res = 0;
    stat->relres = 0;
    stat->converged = 0;
    stat->need_restart = 0;
    stat->restarts = 0;
    return 0;
}

