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
 * @file include/mess/mess.h
 * @brief Main @mess include file.
 *
 * This file is the main header file for @mess. \n
 * Please include this file in every programm you use @mess.
 */

#ifndef MESS_H_
#define MESS_H_

#include <stdint.h>
#include "mess/config.h"

/**
 * @addtogroup general
 * @{ */

#define DEPRECATED  __attribute__ ((deprecated))


#ifndef __cplusplus
#include <complex.h>
#endif

/**
 * @brief Default double complex of @mess.
 *
 * The @ref mess_double_cpx_t typedef is used to avoid incompatability problems when including  @mess from C++. \n
 * C++ treats complex is a template.
 */
typedef double _Complex mess_double_cpx_t;


/**
 * @brief Default float complex of @mess.
 *
 * The @ref mess_float_cpx_t typedef is used to avoid incompatability problems when including  @mess from C++. \n
 * C++ treats complex is a template.
 */
typedef float _Complex mess_float_cpx_t;



#if defined(MESS64)
#include <stdint.h>
#include <inttypes.h>
typedef int64_t mess_int_t;
#define MESS_PRINTF_INT "%"PRId64
#define MESS_PRINTF_INT2D "%2"PRId64
#else
/**
 * @brief Default integer of @mess.
 *
 * The @ref mess_int_t typedef is used to make @mess independent from the size of the operating
 * systems integer. \n
 * It allows to change between 32 and 64  integers without rewriting code. Only
 * the return values of function are always 32bit integers. \n
 * Per default the @ref mess_int_t is equal to @c int. \n
 * If @mess is configured with @c MESS64=ON than it is equal to @c int64_t.
 */
typedef int mess_int_t;



/**
 * @brief String to print a @ref mess_int_t inside @c printf.
 *
 * The @ref MESS_PRINTF_INT is used to display a @ref mess_int_t inside @c printf and similar functions.
 * @sa mess_int_t
 **/
#define MESS_PRINTF_INT "%d"

/**
 * @brief String to print a @ref mess_int_t with 2 digits inside @c printf.
 *
 * The @ref MESS_PRINTF_INT2D is used to display a @ref mess_int_t inside @c printf and similar functions with 2 digits.
 * @sa mess_int_t
 **/
#define MESS_PRINTF_INT2D "%2d"

#endif

/**
 * @brief Type definition for the data type in @mess.
 *
 * The @ref MESS_STORETYPE is used to represent the storage type flag inside @ref mess_matrix and @ref mess_vector.
 **/
typedef unsigned short MESS_STORETYPE;

/**
 * @brief Type definition for the data type in @mess.
 *
 * The @ref MESS_SYMMETRY is used to represent the symmetry type flag inside @ref mess_matrix.
 **/
typedef unsigned short MESS_SYMMETRY;

/**
 * @brief Transpose type for applying a @ref mess_matrix object.
 *
 * The @ref mess_operation_t enumeration is used to define how a subroutine
 * acts on a @ref mess_matrix.
 */
typedef enum {
    /** @ref mess_matrix object is applied non-transposed.  */
    MESS_OP_NONE =  0x00,
    /** @ref mess_matrix object is applied transposed.  */
    MESS_OP_TRANSPOSE,
    /** @ref mess_matrix object is applied hermitian transposed. */
    MESS_OP_HERMITIAN
} mess_operation_t;


/**
 * @brief Enumerator for possible data types.
 *
 * The @ref mess_datatype_t enumeration defines the possible data types inside  @mess \n
 * The @ref mess_datatype_t is used inside @ref mess_matrix, @ref mess_vector and @ref mess_direct.
 *
 **/
typedef enum {
    /** The values of the object are real.  */
    MESS_REAL = 0,
    /** The values of the object are complex.  */
    MESS_COMPLEX
} mess_datatype_t;

/**
 * @brief Enumeration to represent type of norm.
 *
 * The @ref mess_norm_t enumeration is used to represent the type of norm.
 */
typedef enum {
    /** The 2-norm. */
    MESS_2_NORM = 0,
    /** The Frobenius-norm. */
    MESS_FROBENIUS_NORM = 1,
    /** The 1-norm, maximum absolute column sum norm. */
    MESS_1_NORM = 2,
    /** The inf-norm, maximum absolute row sum norm. */
    MESS_INF_NORM = 3,
} mess_norm_t;


/**
 * @brief Enumeration to control the output of print functions.
 *
 * The @ref mess_print_format_t enumeration is used to control the output of print functions.
 */
typedef enum {
    /** Short scientific notation with 6 digits after decimal point.*/
    MESS_PRINT_FORMAT_SHORT = 0,
    /** Long scientific notation with 15 digits after decimal point.*/
    MESS_PRINT_FORMAT_LONG = 1,
} mess_print_format_t;


/**
 * @brief Define for compatibility with the @mm file format.
 *
 * The @ref MESS_INTEGER define is only included for compatibility reasons with
 * the @mm file format.\n
 * It should not be used in programs.
 * @sa mess_datatype_t
 * */
#define MESS_INTEGER MESS_REAL

/**
 * @brief Check if a @ref mess_matrix, @ref mess_vector or @ref mess_direct object is real.
 *
 * The @ref MESS_IS_REAL macro is used to check easily for a real
 * @ref mess_matrix, @ref mess_vector or  @ref mess_direct.
 * */
#define MESS_IS_REAL(m)     ((*m).data_type == MESS_REAL)

/**
 * @brief Check if a @ref mess_matrix, @ref mess_vector or @ref mess_direct object is real.
 *
 * The @ref MESS_IS_INTEGER is a wrapper around \ref MESS_IS_REAL. \n
 * It exists for compatibility reasons with the @mm file format.
 **/
#define MESS_IS_INTEGER(m)      ((*m).data_type == MESS_REAL)
/**
 * @brief Check if a @ref mess_matrix, @ref mess_vector or @ref mess_direct object is complex.
 *
 * The @ref MESS_IS_COMPLEX macro is used to check easily for a complex
 * @ref mess_matrix, @ref mess_vector or @ref mess_direct.
 **/
#define MESS_IS_COMPLEX(m)  ((*m).data_type == MESS_COMPLEX)


/** @} */


#include "mess/errors.h"
#include "mess/matrix.h"
#include "mess/vector.h"

#include "mess/solvername.h"
#include "mess/solver.h"
#include "mess/eigenvalue.h"
#include "mess/direct.h"
#include "mess/lrcf_adi.h"
#include "mess/misc.h"
#include "mess/threadpool.h"
#include "mess/plot_x11.h"
#include "mess/plot_scriptExporter.h"
#include "mess/dynsys.h"
#include "mess/h2.h"
#include "mess/freelist.h"
#include "mess/easyfrontend.h"

#endif /*MESS_H_*/
