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
 * @file lib/vector/init.c
 * @brief Init function  of @ref mess_vector structures.
 * @author @koehlerm
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>
#include "blas_defs.h"
#ifdef _OPENMP_H
#include <omp.h>
#endif

/**
 * On MacOS X 10.6 the zdotc implementation
 * is defect, this is the work around
 */
#ifdef MESS_USE_APPLE_BLAS
#include <Accelerate/Accelerate.h>
#endif



/**
 * @brief Initialize the @ref mess_vector  structure.
 * @param[in,out] v a pointer to a @ref mess_vector object
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_vector_init function creates a new @ref mess_vector object. \n
 * It allocates the basic data structures and sets all values to 0. \n
 *
 */
int mess_vector_init(mess_vector *v){
    MSG_FNAME(__func__);
    mess_try_alloc(*v , struct mess_vector_st*,sizeof(struct mess_vector_st));
    memset(*v, 0, sizeof(struct mess_vector_st));
    (*v)->data_type = MESS_REAL;
    return 0;
}





