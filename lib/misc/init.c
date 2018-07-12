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
 * @file lib/misc/init.c
 * @brief Initialize @mess.
 * @author @koehlerm
 */



#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "cscutils/error_message.h"
#include <complex.h>

static int init_done = 0;
static int exit_done = 0;

int mess_error_level = LOGLEVEL;

static void print_error(const char * str)
{
    fprintf(stderr, "\033[1;2;31m%s\033[0m",str);
    fflush(stderr);
}

static void print_warn(const char * str ) {
    fprintf(stderr, "\033[36m%s\033[0m", str);
    fflush(stderr);
}

static void print_info(const char * str ) {
    fprintf(stderr, "\033[32m%s\033[0m", str);
    fflush(stderr);
}

#pragma GCC diagnostic ignored "-Wformat-security"
static void print_print(const char * str ) {
    fprintf(stdout,str);
    fflush(stdout);
}


/**
 * @brief Initialize @mess.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_init function intializes @mess. \n
 * You have only to call it, if you want to use the dynamic solver registration in a static library.
 *
 */
__attribute__((constructor)) int mess_init(){
    // MSG_FNAME(__func__);
    if ( init_done ) return 0;
    init_done = 1;
    csc_error_message_handle(print_error);
    csc_warn_message_handle(print_warn);
    csc_info_message_handle(print_info);
    csc_print_message_handle(print_print);
    atexit((void *) mess_exit);
    return 0;
}


/**
 * @brief Clean up some initial @mess structures.
 *
 * The @ref mess_exit function clears the configuration.\n
 * Normaly you do not have to run it.
 *
 */
int mess_exit(){
    if ( exit_done ) return 0;
    exit_done = 1;

    return 0;
}

