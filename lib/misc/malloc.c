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
 * @file lib/misc/malloc.c
 * @brief Alloc functions for @cmess.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include "mess/config.h"
#include "mess/error_macro.h"
#include "mess/mess.h"


typedef void * (*malloc_call_t)(size_t size);
typedef void * (*realloc_call_t)(void *ptr, size_t size);
typedef void * (*calloc_call_t)(size_t nmeb, size_t selem);
typedef void (*free_call_t)(void *ptr);

static malloc_call_t   fn_malloc  = (void *) &malloc;
static realloc_call_t  fn_realloc = (void *) &realloc;
static calloc_call_t   fn_calloc  = (void *) &calloc;
static free_call_t     fn_free    = (void *) &free;


/**
 * @brief Set the malloc implementation to be used.
 * @param[in] malloc_ptr input    Pointer to the malloc function
 * @param[in] realloc_ptr input   Pointer to the realloc function
 * @param[in] calloc_ptr input    Pointer to the calloc function
 * @param[in] free_ptr input      Pointer to the free function
 *
 * The @ref mess_set_malloc function sets a new memory management function. That
 * means that the typical C memory management functions are replaced by the
 * given ones. That given function must have the same calling sequence as the
 * ones provide by libc. If the mess_set_malloc function is never called, the standard
 * libc routines are used. If the default needs to be reseted use
 * \code{.c}
 *  @ref mess_set_malloc(NULL, NULL, NULL, NULL);
 * \endcode
 *
 * The primary intention behind this function is to be able to use the @matlab memory
 * management when compiling a MEX file.
 *
 * \remark
 * It only affects memory allocation done by @mess Memory allocations done by subpackages
 * like @umfpack, @csparse, ... are not covered by this mechanism.
 *
 * \attention
 * The function access internal static variables and is not thread-safe. It must be called
 * outside any multithreaded code.
 *
 */
void mess_set_malloc ( void * malloc_ptr, void * realloc_ptr, void * calloc_ptr, void *free_ptr )
{
    if ( malloc_ptr != NULL ) {
        fn_malloc = malloc_ptr;
    } else {
        fn_malloc = &malloc;
    }
    if ( realloc_ptr != NULL ) {
        fn_realloc = realloc_ptr;
    } else {
        fn_realloc = &realloc;
    }
    if ( calloc_ptr != NULL ) {
        fn_calloc = calloc_ptr;
    } else {
        fn_calloc = &calloc;
    }
    if ( free_ptr != NULL ) {
        fn_free = free_ptr;
    } else {
        fn_free = &free;
    }

    return;
}       /* -----  end of function mess_set_malloc  ----- */


/**
 * @brief Wrapper around malloc.
 * @param[in] size input size of the required memory location
 * @return a pointer to the allocated memory
 *
 * The @ref __mess_malloc function is a configurable wrapper around malloc.
 *
 */
void * __mess_malloc ( size_t size )
{

    return fn_malloc(size) ;
}       /* -----  end of function mess_malloc  ----- */


/**
 * @brief Wrapper around realloc.
 * @param[in] ptr input    Old pointer to the memory location
 * @param[in] newsize input      New size of the memory location
 * @return a pointer to the new resized memory location
 *
 * The @ref __mess_realloc function is a configurable wrapper around realloc.
 *
 */
void * __mess_realloc ( void * ptr, size_t newsize )
{

    return fn_realloc(ptr, newsize);
}       /* -----  end of function __mess_realloc  ----- */


/**
 * @brief Wrapper around calloc.
 * @param[in] nmemb input  Number of elements to allocate
 * @param[in] size input  Size of one element
 * @return a pointer to the allocated memory location.
 *
 * The @ref __mess_calloc function is a configurable wrapper around calloc
 *
 */
void * __mess_calloc ( size_t nmemb, size_t size   )
{

    return fn_calloc(nmemb, size) ;
}       /* -----  end of function __mess_calloc  ----- */


/**
 * @brief Wrapper around free.
 * @param[in] ptr input Pointer to the memory which should be freed
 *
 * The @ref __mess_free function is a configurable wrapper around free.
 *
 */
void __mess_free ( void *ptr )
{
    fn_free(ptr);
}       /* -----  end of function __mess_free  ----- */


