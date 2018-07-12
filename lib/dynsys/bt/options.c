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
 * @file lib/dynsys/bt/options.c
 * @brief Options for balance truncation.
 * @author @koehlerm
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
 * @brief Initialize a mess_bt_options object.
 * @param[in,out] opt   pointer to the object
 * @return always zero
 *
 * The @ref mess_bt_options_init function initializes a mess_bt_options object and
 * set default values.
 *
 */
int  mess_bt_options_init ( mess_bt_options * opt )
{
    MSG_FNAME(__func__);
    mess_try_alloc(*opt, struct mess_bt_options_st*, sizeof (struct mess_bt_options_st ));
    (*opt)->chooseorder = NULL;
    (*opt)->chooseorder_aux = NULL;
    (*opt)->rdim = 20;
    (*opt)->tol = 1e-7;
    (*opt)->so_type = MESS_BT_POSITION;
    return 0;
}       /* -----  end of function mess_bt_options_init  ----- */


/**
 * @brief Clear a mess_bt_options object.
 * @param[in] opt  input pointer to the object
 * @return always zero
 *
 * The @ref mess_bt_options_clear function clears a mess_bt_options object.
 *
 */
int  mess_bt_options_clear ( mess_bt_options *opt )
{
    if ( opt == NULL) return 0;
    if ( *opt == NULL) return 0;
    mess_free(*opt);
    *opt = NULL;
    return 0;
}       /* -----  end of function mess_bt_options_clear  ----- */


/**
 * @brief Print a mess_bt_options object to stdout.
 * @param[in] opt    input mess_bt_options object
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_bt_options_print function prints a mess_bt_options object to stdout.
 *
 */
int  mess_bt_options_print ( mess_bt_options opt )
{
    MSG_FNAME(__func__);
    mess_check_nullpointer(opt);
    MSG_PRINT("ROM dimension:  " MESS_PRINTF_INT "\n", opt->rdim);
    MSG_PRINT("tolerance:      %.10e\n", opt->tol);
    if ( opt->chooseorder !=NULL){
        MSG_PRINT("chooseorder:    userdefined\n");
    }else {
        MSG_PRINT("chooseorder:    default\n");
    }
    MSG_PRINT("second order type: %s\n", mess_bt_getsotypestr(opt->so_type));
    return 0;
}       /* -----  end of function mess_bt_options_print  ----- */


/**
 * @brief Get the string of the so_type value.
 * @param[in] so_type input     number of the type
 * @return zero on success and a string of the so_type
 *
 * The @ref mess_bt_getsotypestr function gets the string of
 * a so_type value.
 *
 */
const char * mess_bt_getsotypestr( unsigned short so_type  )
{
    switch( so_type) {
        case MESS_BT_POSITION:
            return "position";
        case MESS_BT_VELOCITY:
            return "velocity";
        case MESS_BT_POSITION_VELOCITY:
            return "position-velocity";
        case MESS_BT_VELOCITY_POSITION:
            return "velocity-position";
        default:
            return "unknown";
    }
    return 0;
}       /* -----  end of function mess_bt_getsotypestr(  ----- */

