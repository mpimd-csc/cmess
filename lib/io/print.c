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
 * @file lib/io/matrix_print.c
 * @brief Print formatting.
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>
#include <pthread.h>
#include <math.h>



static mess_print_format_t print_format = MESS_PRINT_FORMAT_SHORT;
static pthread_mutex_t print_format_mutex = PTHREAD_MUTEX_INITIALIZER;


/**
 * @brief Change the default behaviour of printing.
 * @param[in] p     input @ref mess_print_format_t  type to change printing behaviour
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_print_format_select function selects a printing format used for printing @p double and @ref mess_double_cpx_t values of  @ref mess_matrix
 * and @ref mess_vector structures.
 *
 */
int mess_print_format_select(mess_print_format_t p){
    pthread_mutex_lock(&print_format_mutex);
    print_format = p;
    pthread_mutex_unlock(&print_format_mutex);
    return 0;
}


/**
 * @brief Print a double value using @ref MSG_PRINT.
 * @param[in] d     input value to print
 * @return lenght of printed characters.
 *
 */
int mess_print_format_double(double d){
    switch(print_format){
        case MESS_PRINT_FORMAT_SHORT:
            MSG_PRINT("% .6e",d);
            return 3+6+4;
        case MESS_PRINT_FORMAT_LONG:
        default:
            MSG_PRINT("% .15e",d);
            return 3+15+4;
    }
}


/**
 * @brief Print a @ref mess_double_cpx_t value using @ref MSG_PRINT.
 * @param[in] cpx     input value to print
 * @return lenght of printed characters.
 *
 */
int mess_print_format_double_cpx(mess_double_cpx_t cpx){
    switch(print_format){
        case MESS_PRINT_FORMAT_SHORT:
            MSG_PRINT("% .6e %c %.6e i",creal(cpx),cimag(cpx)<0?'-':'+',fabs(cimag(cpx)));
            return (3+6+4)+3+(2+6+4)+2;
        case MESS_PRINT_FORMAT_LONG:
        default:
            MSG_PRINT("% .15e %c %.15e i",creal(cpx),cimag(cpx)<0?'-':'+',fabs(cimag(cpx)));
            return (3+15+4)+3+(2+15+4)+2;
    }
}

