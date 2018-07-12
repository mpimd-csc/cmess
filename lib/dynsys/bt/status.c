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
 * @file lib/dynsys/bt/status.c
 * @brief Work with the mess_bt_status objects.
 * @author @koehlerm
 *
 *    Description:
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"

#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Initialize a mess_bt_status object.
 * @param[in,out] status pointer to the mess_bt_status object
 * @return always zero
 *
 * The @ref mess_bt_status_init function initializes a mess_bt_status object.
 *
 */
int mess_bt_status_init ( mess_bt_status * status)
{
    MSG_FNAME(__func__);
    int ret = 0 ;
    mess_try_alloc(*status, mess_bt_status, sizeof(struct mess_bt_status_st));
    ret  = mess_status_init(&((*status)->statB));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_status_init);
    ret  = mess_status_init(&((*status)->statC));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_status_init);
    (*status)->time = 0;
    (*status)->time_lyap = 0;
    (*status)->time_VW  = 0;
    (*status)->rdim  =0;
    (*status)->esterror = 0 ;

    return 0;
}       /* -----  end of function mess_bt_status_init  ----- */

/**
 * @brief Clear a mess_bt_status object.
 * @param[in] status     input pointer to the mess_bt_status object
 * @return always zero
 *
 * The @ref mess_bt_status_clear function clears a mess_bt_status object.
 *
 */
int  mess_bt_status_clear ( mess_bt_status *status )
{
    if ( status == NULL) return 0;
    if ( *status == NULL) return 0;
    mess_status_clear(&((*status)->statB));
    mess_status_clear(&((*status)->statC));
    mess_free(*status);
    return 0;
}       /* -----  end of function mess_bt_status_clear  ----- */


/**
 * @brief Print a mess_bt_status object to stdout.
 * @param[in] status     input mess_bt_status object to print
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_bt_status_print function prints a mess_bt_status object to stdout.
 *
 */
int mess_bt_status_print (mess_bt_status status )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(status);
    MSG_PRINT("ROM dimension:    " MESS_PRINTF_INT "\n", status->rdim);
    MSG_PRINT("est. error:       %.10e\n", status->esterror);
    MSG_PRINT("time (overall):   %lgs\n", status->time);
    MSG_PRINT("time (lyap part): %lgs\n", status->time_lyap);
    MSG_PRINT("time (comp. V/W): %lgs\n", status->time_VW);
    MSG_PRINT("LRCF-ADI status for AX+XA^T=-BB^T\n");
    if (status->statB != NULL ) {
        ret = mess_status_print(status->statB);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_status_print);
    }
    if ( status->statC != NULL ) {
        MSG_PRINT("LRCF-ADI status for A^TX+XA=-C^TC\n");
        ret = mess_status_print(status->statC);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_status_print);
    }
    return 0;
}       /* -----  end of function mess_bt_status_print  ----- */

