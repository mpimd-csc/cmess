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
 * @file lib/misc/openmp.c
 * @brief OpenMP utilities.
 * @author @koehlerm
 */

#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include <pthread.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "mess/error_macro.h"
#include "mess/mess.h"

/**
 * @brief Display information about OpenMP.
 *
 * The @ref mess_openmp_info function displays some information
 * about the OpenMP environment including
 * <ul>
 * <li> current threads,
 * <li> number of CPUs,
 * <li> maximal threads,
 * <li> thread ID.
 * </ul>
 *
 */
void mess_openmp_info(){
#ifdef _OPENMP
    int curid = omp_get_thread_num();
    MSG_PRINT("<%2d> number of cpus:  %d\n", curid, omp_get_num_procs());
    MSG_PRINT("<%2d> max threads:     %d\n", curid, omp_get_max_threads());
    MSG_PRINT("<%2d> current threads: %d\n", curid, omp_get_num_threads());
    MSG_PRINT("<%2d> thread ID:       %d\n", curid, omp_get_thread_num());

#else
    MSG_PRINT("no OpenMP available\n");
#endif
}



