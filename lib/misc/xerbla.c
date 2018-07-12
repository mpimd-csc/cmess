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
 * @file lib/misc/xerbla.c
 * @brief New error handler for @lapack.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/eigenvalue.h"
#include "mess/error_macro.h"
#include <complex.h>
#include <math.h>

/**
 * @brief Override @lapack xerbla function.
 * @param[in] fname input function name
 * @param[in] info input    info code
 *
 * The @ref xerbla_ function overrides @lapack xerbla function.
 *
 */
void xerbla_( char* fname, mess_int_t *info )
{
    MSG_FNAME(__func__);
    MSG_ERROR("On entry to %s argument " MESS_PRINTF_INT " has the wrong value.\n",  fname, *info);
#ifdef MESS_DEBUG
    abort();
#endif
}       /* -----  end of function xerbla  ----- */

