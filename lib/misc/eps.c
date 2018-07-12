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
 * @file lib/misc/eps.c
 * @brief Compute machine precision.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>

double F77_GLOBAL(dlamch,DLAMCH)(char *what);
/**
 * @brief Calcute the machine epsilon.
 * @return machine epsilon
 *
 * The @ref mess_eps function calculates the machine precision.
 */
double mess_eps()
{
    /* double machEps = 1.0;
       do {
       machEps /= 2.0;
    // If next epsilon yields 1, then break, because current
    // epsilon is the machine epsilon.
    }
    while ((double)(1.0 + (machEps/2.0)) != 1.0);
    // CSS_Mach_Eps = machEps;
    return (double) machEps; */
    return F77_GLOBAL(dlamch,DLAMCH)("E")*2;
}
