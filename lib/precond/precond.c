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
 * @file lib/precond/precond.c
 * @brief Preconditioner interface.
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
 * @brief Initialize a mess_precond object.
 * @param[in,out] precond pointer to a mess_precond object
 * @return always zero
 *
 * The @ref mess_precond_init function sets all values in mess_precond to zero and all pointers to @c NULL.
 *
 */
int mess_precond_init(mess_precond *precond)
{
    MSG_FNAME(__func__);
    mess_try_alloc(*precond, mess_precond, sizeof(struct mess_precond_st));
    (*precond)->clear = NULL;
    (*precond)->data = NULL;
    (*precond)->solve = NULL;
    (*precond)->type = 0;
    return 0;
}

/**
 * @brief Clean up the mess_precond structure to use it again.
 * @param[in,out] precond pointer to the mess_precond
 * @return always zero
 *
 * The @ref mess_precond_clear function removes the mess_precond object from memory and cleans up.
 *
 */
int mess_precond_clear(mess_precond *precond)
{
    // MSG_FNAME(__func__);
    if (*precond == NULL) {
        return 0;
    }
    if  ((*precond)->clear != NULL){
        (*precond)->clear(*precond);
    }
    mess_free(*precond);
    *precond = NULL;
    return 0;
}
