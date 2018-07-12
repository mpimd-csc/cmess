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
 * @file include/mess/error_macro.h
 * @author @koehlerm
 * @brief Interface to error macros.
 */
#ifndef  _ERROR_H
#define  _ERROR_H

// #define MESS_HAVE_BACKTRACE
#define BTRSIZE 50

#include "cscutils/error_message.h"
/** @addtogroup error_macro
 * @{ */

#include "mess/errors.h"
/** Print error messages on stderr if LOGLEVEL is bigger than zero. */
#define     MSG_ERROR( text, args... )  if ( mess_error_level > 0){ \
    csc_error_message("%s: %s(%5d) - error: \t" text ,__FILE__,__funcname,__LINE__, ## args); \
    csc_show_backtrace(); \
    fflush(stderr);\
}

/** Print warning messages on stderr if LOGLEVEL is bigger than one. */
#define     MSG_WARN(text, args...)     if ( mess_error_level > 1) { \
    csc_warn_message("%s: %s(%5d) - warning: \t" text ,__FILE__,__funcname,__LINE__, ## args); \
    csc_show_backtrace(); \
    fflush(stderr); \
}

/** Print error messages on stderr if LOGLEVEL is bigger than three.*/
#define     MSG_INFO(text, args...)     if ( mess_error_level > 2) { \
    csc_info_message("%s: %s(%5d) - info: \t" text ,__FILE__,__funcname,__LINE__, ## args); \
    fflush(stderr); \
}


/** Print a warning without the backtrace */
#define     MSG_WARN_NOBT(text, args...)        if ( mess_error_level > 0) { \
    csc_warn_message("%s: %s(%5d) - warning: \t" text ,__FILE__,__funcname,__LINE__, ## args); \
    fflush(stderr); \
}

/** Print usual messages on stdout.*/
#define     MSG_PRINT(text, args...) {               \
    csc_print_message(text, ## args); \
    fflush(stdout); \
}


/** Name of the current function */
#define     MSG_FNAME( nn )         static const char*  __funcname = nn

/** Function failure handler */

#ifdef PYMESS
#undef _POSIX_C_SOURCE
#include <Python.h>
#define FUNCTION_FAILURE_HANDLE( err, cond , fun )                                                  \
    if(PyErr_CheckSignals()!=0){                                                                    \
        MSG_ERROR(" %s returned with got Python Ctrl-C Signal - %s\n", #fun);                       \
        return ((int)MESS_ERROR_PYTHON_SIG);                                                        \
    }                                                                                               \
if ( (cond)) {                                                                                  \
    MSG_ERROR(" %s returned with %d - %s\n", #fun, (int)err, mess_get_error((int)err));         \
    return ((int)err);                                                                          \
}

#define FUNCTION_FAILURE_HANDLE_OMP( err, cond , fun )                                              \
    if ( (cond)) {                                                                                  \
        MSG_ERROR(" %s returned with %d - %s\n", #fun, (int) err, mess_get_error((int)err));        \
    }

#else

#define FUNCTION_FAILURE_HANDLE( err, cond , fun )                                                  \
    if ( (cond)) {                                                                                  \
        MSG_ERROR(" %s returned with %d - %s\n", #fun, (int)err, mess_get_error((int)err));         \
        return ((int)err);                                                                          \
    }

#define FUNCTION_FAILURE_HANDLE_OMP( err, cond , fun )                                              \
    if ( (cond)) {                                                                                  \
        MSG_ERROR(" %s returned with %d - %s\n", #fun, (int) err, mess_get_error((int)err));        \
    }

#endif


/** @} */

#endif

