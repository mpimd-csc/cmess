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
 * @file lib/direct/singlesolver/umfpack.c
 * @brief Interface to @umfpack.
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#ifdef MESS_HAVE_UMFPACK
#include <umfpack.h>

#ifdef MESS64
#define umfpack_dl_solve    umfpack_dl_solve
#define umfpack_zl_solve    umfpack_zl_solve
#define umfpack_zl_load_numeric umfpack_zl_load_numeric
#define umfpack_dl_load_numeric umfpack_dl_load_numeric
#define umfpack_zl_defaults     umfpack_zl_defaults
#define umfpack_dl_defaults     umfpack_dl_defaults
#define umfpack_dl_save_numeric umfpack_dl_save_numeric
#define umfpack_zl_save_numeric umfpack_zl_save_numeric
#define umfpack_zl_free_numeric umfpack_zl_free_numeric
#define umfpack_dl_free_numeric umfpack_dl_free_numeric
#define umfpack_dl_symbolic umfpack_dl_symbolic
#define umfpack_dl_numeric  umfpack_dl_numeric
#define umfpack_zl_symbolic     umfpack_zl_symbolic
#define umfpack_zl_numeric  umfpack_zl_numeric
#define umfpack_dl_report_info  umfpack_dl_report_info
#define umfpack_dl_free_symbolic umfpack_dl_free_symbolic
#define umfpack_zl_report_info   umfpack_zl_report_info
#define umfpack_zl_free_symbolic umfpack_zl_free_symbolic
#define umfpack_zl_get_lunz  umfpack_zl_get_lunz
#define umfpack_zl_get_numeric   umfpack_zl_get_numeric
#define umfpack_dl_get_lunz      umfpack_dl_get_lunz
#define umfpack_dl_get_numeric   umfpack_dl_get_numeric


#else
#define umfpack_dl_solve    umfpack_di_solve
#define umfpack_zl_solve    umfpack_zi_solve
#define umfpack_zl_load_numeric umfpack_zi_load_numeric
#define umfpack_dl_load_numeric umfpack_di_load_numeric
#define umfpack_zl_defaults     umfpack_zi_defaults
#define umfpack_dl_defaults     umfpack_di_defaults
#define umfpack_dl_save_numeric umfpack_di_save_numeric
#define umfpack_zl_save_numeric umfpack_zi_save_numeric
#define umfpack_zl_free_numeric umfpack_zi_free_numeric
#define umfpack_dl_free_numeric umfpack_di_free_numeric
#define umfpack_dl_symbolic umfpack_di_symbolic
#define umfpack_dl_numeric  umfpack_di_numeric
#define umfpack_zl_symbolic     umfpack_zi_symbolic
#define umfpack_zl_numeric  umfpack_zi_numeric
#define umfpack_dl_report_info  umfpack_di_report_info
#define umfpack_dl_free_symbolic umfpack_di_free_symbolic
#define umfpack_zl_report_info   umfpack_zi_report_info
#define umfpack_zl_free_symbolic umfpack_zi_free_symbolic
#define umfpack_zl_get_lunz  umfpack_zi_get_lunz
#define umfpack_zl_get_numeric   umfpack_zi_get_numeric
#define umfpack_dl_get_lunz      umfpack_di_get_lunz
#define umfpack_dl_get_numeric   umfpack_di_get_numeric

#endif

/**
 * @brief Generate a direct linear solver for standard linear systems \f$ Ax=b \f$ with @umfpack
 * @param[in] matrix matrix to decompose
 * @param[out] solver output solver
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref  mess_direct_create_umfpack function factorizes a matrix with @umfpack with iterative refinement and provides
 * a solver for linear systems.
 *
 */
int mess_direct_create_umfpack(mess_matrix matrix, mess_direct solver){
    MSG_FNAME(__func__);
    double Control[UMFPACK_CONTROL];
    mess_int_t ret;
    /*-----------------------------------------------------------------------------
     *  check input pointers
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_real_or_complex(matrix);
    mess_check_square(matrix);

    /*-----------------------------------------------------------------------------
     *  create umfpack control
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_COMPLEX(matrix)) {
        umfpack_zl_defaults(Control);
    } else {
        umfpack_dl_defaults(Control);
    }

    /*-----------------------------------------------------------------------------
     *  create umfpack based solver
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_create_umfpack_control(matrix,solver,Control);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_create_umfpack_control);
    //overwrite solver name
    if(solver->name){
        MESS_CLEAR_POINTERS(solver->name);
        SET_SOLVERNAME(solver->name, __func__);
    }
    return ret;
}

#endif
