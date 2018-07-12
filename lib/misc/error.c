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
 * @file lib/misc/error.c
 * @brief Error level handling.
 * @author @koehlerm
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "cscutils/error_message.h"

/**
 * @brief Set error level.
 * @param[in] level input   number of error level to set
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_set_errorlevel function sets the error level of @mess \n
 * Possible error values are:
 * \li \f$ 0 \f$ - display no errors,
 * \li \f$ 1 \f$ - display only errors (default for release),
 * \li \f$ 2 \f$ - display errors and warnings,
 * \li \f$ 3 \f$ - display errors, warnings and infos (default for debug).
 *
 */
int mess_set_errorlevel(int level){
    MSG_FNAME(__func__);
    if ( level < 0 || level > 4){
        MSG_ERROR("level must be between 0 and 4\n");
        return MESS_ERROR_ARGUMENTS;
    }
    mess_error_level = level;
    return 0;
}


/**
 * @brief Return an error message.
 * @param[in] code input number of error
 * @return string with error message described in \f$ code \f$
 *
 * The @ref mess_get_error function returns a string which describes the error named in \f$ code \f$.
 *
 */
char *mess_get_error(int code) {
    /** Successful return */
    switch (code){
        case MESS_ERROR_PYTHON:
            return "Error in the Python interface. Please check your code.";
        case MESS_ERROR_DYNSYS:
            return "Something went wrong with a dynamical system.";
        case MESS_ERROR_CONVERGE:
            return "an iterative algorithm don't converge.";
        case MESS_ERROR_GENERAL:
            return "a general error occured.";
        case MESS_ERROR_MEMLOC:
            return "data is not stored on the required memory location.";
        case MESS_ERROR_LAPACK:
            return "A LAPACK subroutine returned an error.";
        case MESS_ERROR_EIGENVALUES:
            return "At least one of the passed matrices has eigenvalues, which must not appear.";
        case MESS_ERROR_NOSUPPORT:
            return "The called function isn't suported by the compile time configuration.";
        case MESS_ERROR_MISSING:
            return "The called function isn't implemented ( -> direct/multidirect solver).";
        case MESS_ERROR_UMFPACK:
            return "UMFPACK returned an error.";
        case MESS_ERROR_MEMORY:
            return "memory allocation was not possible.";
        case MESS_ERROR_NOTLOAD:
            return "The requested index in a multisolver isn't loaded.";
        case MESS_ERROR_NOPATTERN:
            return "The called function doesn't work on pattern matrices.";
        case MESS_ERROR_NOCOMPLEX:
            return "The called function doesn't work on complex matrices.";
        case MESS_ERROR_SINGULAR:
            return "The given matrix is singular.";
        case MESS_ERROR_DIMENSION:
            return "The dimension of the given objects don't fit.";
        case MESS_ERROR_CONVERT:
            return "Convertation between the matrix formats is impossible or caused an error.";
        case MESS_ERROR_SYMMETRIC:
            return "problems with the symmetric matrices.";
        case MESS_ERROR_ARGUMENTS:
            return "the passed arguments are incorrect.";
        case MESS_ERROR_STORAGETYPE:
            return "the matrix doesn't have the right storage format.";
        case MESS_ERROR_DATATYPE:
            return "Some data don't have have the rigth data type";
        case MESS_ERROR_DATA:
            return "Wrong (input) data (while file i/o).";
        case MESS_ERROR_WRONG_HEADER:
            return "The input file doesn't have a correct header.";
        case MESS_ERROR_FILEIO:
            return "error while reading/writing files.";
        case MESS_ERROR_NULLPOINTER:
            return "oops.... a null pointer. no one is perfect :-).";
        case MESS_ERROR_NOTNULLPOINTER:
            return "Pointer should points to NULL.";
        case MESS_ERROR_CHOLMOD:
            return "Error with CHOLMOD occured.";
        case MESS_ERROR_MKLPARDISO:
            return "Error with MKL-PARDISO occured.";
        case MESS_SUCCESS:
            return "returns with no error. all fine.";
        case MESS_ERROR_MISC:
            return "A generic/misc error occured.";
        case MESS_ERROR_PYTHON_SIG:
            return "Received Ctrl-C from Python.";
        default:
            return "an unspecified error occured.";

    }
    return "";
}

