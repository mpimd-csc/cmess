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
 * @file include/mess/errors.h
 * @brief Interface to error codes.
 * @author @koehlerm
 */



#ifndef ERRORS_H_
#define ERRORS_H_

/** @addtogroup error_codes
 * @{ */
#ifdef __cplusplus
extern "C" {
#endif
    /** Global variable to set the error level/verbosity. */
    extern int mess_error_level;
    char*  mess_get_error(int code);
    int mess_set_errorlevel(int level);

#ifdef __cplusplus
}
#endif

/** Successful return */
#define MESS_SUCCESS                0

/** Error in the Python interface or Executing @python.  */
#define MESS_ERROR_PYTHON           225

/** Dynamical system error */
#define MESS_ERROR_DYNSYS           229

/** Convergence error */
#define MESS_ERROR_CONVERGE         230

/** General error */
#define MESS_ERROR_GENERAL          231

/** Memory location error */
#define MESS_ERROR_MEMLOC           233

/** @lapack returned an error */
#define MESS_ERROR_LAPACK           234

/** The called function has problems with the eigenvalues of the matrix. */
#define MESS_ERROR_EIGENVALUES      235

/** The called function is not supported by the compile time configuration. */
#define MESS_ERROR_NOSUPPORT        236

/** The called function is not implemented. */
#define MESS_ERROR_MISSING          237

/** @umfpack returned an error. */
#define MESS_ERROR_UMFPACK          238

/** Memory allocation was not possible. */
#define MESS_ERROR_MEMORY           239

/** The requested index in a multisolver is not loaded. */
#define MESS_ERROR_NOTLOAD          240

/** The called function does not work on pattern matrices. */
#define MESS_ERROR_NOPATTERN        242

/** The called function does not work on complex matrices. */
#define MESS_ERROR_NOCOMPLEX        243

/** The given matrix is singular. */
#define MESS_ERROR_SINGULAR         244

/** The dimension of the given objects does not fit. */
#define MESS_ERROR_DIMENSION        245

/** Conversion between the matrix formats is impossible. */
#define MESS_ERROR_CONVERT          246

/** Problems with symmetric matrices */
#define MESS_ERROR_SYMMETRIC        247

/** Passed arguments are wrong. */
#define MESS_ERROR_ARGUMENTS        248

/** Matrix does not have the right storage format. */
#define MESS_ERROR_STORAGETYPE      249

/** Matrix does not have the right data type. */
#define MESS_ERROR_DATATYPE         250

/** While reading a matrix there is wrong data in the input stream. */
#define MESS_ERROR_DATA             251

/** Matrix file has the wrong file header. */
#define MESS_ERROR_WRONG_HEADER     252

/** Error while file I/O */
#define MESS_ERROR_FILEIO           253

/** A nullpointer is passed to a function. */
#define MESS_ERROR_NULLPOINTER      254

/** A not nullpointer is passed to a function, but it is assumed to point to null. */
#define MESS_ERROR_NOTNULLPOINTER   255

/** Miscellaneous error */
#define MESS_ERROR_MISC             256

/** @cholmod returned an error */
#define MESS_ERROR_CHOLMOD          257

/** @mklpardiso return an error */
#define MESS_ERROR_MKLPARDISO       258

/** Got Ctrl-C from Python */
#define MESS_ERROR_PYTHON_SIG       259



/** @}  */

#endif /*ERRORS_H_*/


