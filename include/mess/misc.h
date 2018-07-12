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
 * @file include/mess/misc.h
 * @brief Miscallenous functions.
 * @author @koehlerm
 */

#ifndef MESS_MISC_H_
#define MESS_MISC_H_

#include <stdio.h>
#include "mess/config.h"

#ifdef __cplusplus
extern "C" {
#endif
    /** @addtogroup misc
     * @{ */
    int mess_init(void);
    int mess_exit(void);

    double mess_eps(void);
    int mess_int_t_size(void);

    void mess_version(void);
    void mess_version_verbose(void);
    int mess_version_major(void);
    int mess_version_minor(void);
    int mess_version_patch(void);
    const char* mess_git_id(void);
    const char* mess_git_branch(void);
    void mess_openmp_info(void);

    int mess_is_debug(void);
    int mess_have_zlib(void);
    int mess_have_bzip2(void);
    int mess_have_umfpack(void);
    int mess_have_amd(void);
    int mess_have_colamd(void);
    int mess_have_cholmod(void);
    int mess_have_csparse(void);
    int mess_have_superlu(void);
    int mess_have_mklpardiso(void);
    int mess_have_arpack(void);
    int mess_have_matio(void);
    int mess_have_openmp(void);
    int mess_have_mess64(void);

    const char *mess_datatype_t_str(const mess_datatype_t data_type);
    const char *mess_storage_t_str(const  mess_storage_t store_type);
    const char *mess_symmetry_t_str(const mess_symmetry_t symmetry);
    const char *mess_operation_t_str(const mess_operation_t op);
    const char *mess_equation_t_str(const mess_equation_t eqn);
    const char *mess_parameter_t_str(const mess_parameter_t para);
    const char *mess_memusage_t_str(const mess_memusage_t mem);
    const char *mess_norm_t_str(const mess_norm_t nrm);
    const char* mess_sign_scale_t_str(const mess_sign_scale_t sc_type);
    const char* mess_direct_lupackage_t_str(const mess_direct_lupackage_t lu_type);
    const char* mess_direct_cholpackage_t_str(const mess_direct_cholpackage_t chol_type);
    const char* mess_multidirect_t_str(const mess_multidirect_t ml_type);
    const char* mess_residual_t_str(const mess_residual_t res_type);

    int mess_print_format_select(mess_print_format_t print_type);
    int mess_print_format_double(double d);
    int mess_print_format_double_cpx(mess_double_cpx_t cpx);


    void mess_print_bytes(mess_int_t size);

    void * __mess_malloc(size_t size);
    void * __mess_realloc( void* ptr, size_t newsize);
    void * __mess_calloc( size_t nmemb, size_t size);
    void __mess_free(void* ptr);
    void mess_set_malloc(void* malloc_ptr, void* reallloc_ptr, void* calloc_ptr, void* free_ptr);

#ifndef MESS_HAVE_NEARBYINT
    double nearbyint(double x);
#endif

#ifndef MESS_HAVE_ISFINITE
#define isfinite(x) \
    (sizeof(x) == sizeof(float))?__finitef(x):(sizeof(x)==sizeof(double))?__finite(x):__finitel(x)
#endif

    //offset functions
    int mess_matrix_st_offset(int* size, int* offsets, char** fieldnames);
    int mess_vector_st_offset(int* size, int* offsets, char** fieldnames);
    int mess_status_st_offset(int* size, int* offsets, char** fieldnames);
    int mess_options_st_offset(int* size, int* offsets, char** fieldnames);
    int mess_equation_st_offset(int* size, int* offsets, char** fieldnames);
    int mess_direct_st_offset(int* size, int* offsets, char** fieldnames);

    /** @}  */


    /** @addtogroup misc_timer
     * @{ */
    double mess_wtime (void);
    double mess_ctime (void);

    /**
     * @brief Representation of a timer.
     *
     * The @ref mess_timer_st structure represents a timer.
     */
    struct mess_timer_st {
        struct timespec *timer; /**<  Pointer to timer object */
        double dtimer;          /**<  Determined time */
    };

    /**
     * @brief Type definition of a timer.
     *
     * The @ref mess_timer type definition defines a timer.
     *
     **/
    typedef struct mess_timer_st *mess_timer;

    int mess_timer_init(mess_timer* timer);
    int mess_timer_clear(mess_timer* timer);
    double mess_timer_get(mess_timer timer);
    double mess_timer_getandreset(mess_timer timer);
    double mess_timer_diff(mess_timer timerStart, mess_timer timerEnd);

    /** @} */

    /** @addtogroup misc_macro
     * @{ */

    /**
     * @brief Macro for Maximium of two numbers.
     * The @ref MESS_MAX macro return the maximum of two numbers.
     */
#define MESS_MAX(a,b) (((a)<(b)) ? (b):(a))

    /**
     * @brief Macro for Minimum of two numbers.
     * The @ref MESS_MIN macro return the minimum of two numbers.
     *
     * */
#define MESS_MIN(a,b) (((a)<(b)) ? (a):(b))

    /**
     * @brief Macro to init a list of matrices.
     *
     * The @ref MESS_INIT_MATRICES macro allow to init a list of matrices within one call. It uses
     * a variable argument list macro. The calling sequence is
     * \code{.c}
     * MESS_INIT_MATRICES(&A, &B, &C, ...);
     * \endcode
     */
#define  MESS_INIT_MATRICES(...) { \
    void * _stopper = (void * ) -1; \
    mess_matrix*  _M[] = { __VA_ARGS__, _stopper}; \
    int _i; \
    for ( _i = 0; _M[_i] != _stopper ; _i++) {\
        mess_matrix_init(_M[_i]);\
    }\
}

    /**
     * @brief Macro to reset a list of matrices.
     *
     * The @ref MESS_RESET_MATRICES macro allow to reset  a list of matrices within one call. It uses
     * a variable argument list macro. The calling sequence is
     * \code{.c}
     * MESS_RESET_MATRICES(A, B, C, ...);
     * \endcode
     */
#define  MESS_RESET_MATRICES(...) { \
    void * _stopper = (void * ) -1; \
    mess_matrix  _M[] = { __VA_ARGS__, _stopper}; \
    int _i; \
    for ( _i = 0; _M[_i] != _stopper ; _i++) {\
        MESS_MATRIX_RESET(_M[_i]);\
    }\
}

    /**
     * @brief Macro to clear a list of matrices.
     *
     * The @ref MESS_CLEAR_MATRICES macro allow to clear a list of matrices within one call. It uses
     * a variable argument list macro. The calling sequence is
     * \code{.c}
     * MESS_CLEAR_MATRICES(&A, &B, &C, ...);
     * \endcode
     * Before the matrices are cleared it is check if they are already @c NULL.
     */
#define  MESS_CLEAR_MATRICES(...) { \
    void * _stopper = (void * ) -1; \
    mess_matrix*  _M[] = { __VA_ARGS__, _stopper}; \
    int _i; \
    for ( _i = 0; _M[_i] != _stopper ; _i++) {\
        if ( *_M[_i] != NULL) mess_matrix_clear(_M[_i]);\
    }\
}

    /**
     * @brief Macro to init a list of @ref mess_vector instances.
     *
     * The @ref MESS_INIT_VECTORS macro allow to init a list of vectors within one call. It uses
     * a variable argument list macro. The calling sequence is
     * \code{.c}
     * MESS_INIT_VECTORS(&v1, &v2, &v3, ...);
     * \endcode
     */
#define  MESS_INIT_VECTORS(...) { \
    void * _stopper = (void * ) -1; \
    mess_vector*  _M[] = { __VA_ARGS__, _stopper}; \
    int _i; \
    for ( _i = 0; _M[_i] != _stopper ; _i++) {\
        mess_vector_init(_M[_i]);\
    }\
}

    /**
     * @brief Macro to clear a list of vectors.
     *
     * The @ref MESS_CLEAR_VECTORS macro allow to clear a list of vectors within one call. It uses
     * a variable argument list macro. The calling sequence is
     * \code{.c}
     * MESS_CLEAR_VECTORS(&v, &v2, &v3, ...);
     * \endcode
     * Before the vectors are cleared it is check if they are already @c NULL.
     * */
#define  MESS_CLEAR_VECTORS(...) { \
    void * _stopper = (void * ) -1; \
    mess_vector*  _M[] = { __VA_ARGS__, _stopper}; \
    int _i; \
    for ( _i = 0; _M[_i] != _stopper ; _i++) {\
        if ( *_M[_i] != NULL) mess_vector_clear(_M[_i]);\
    }\
}

    /**
     * @brief Macro to init  a list of direct solvers.
     *
     * The @ref MESS_INIT_DIRECTS macro allow to init  a list of direct solvers within one call. It uses
     * a variable argument list macro. The calling sequence is
     * \code{.c}
     * MESS_INIT_DIRECTS(&direct1, &direct2, &direct3, ...);
     * \endcode
     */
#define  MESS_INIT_DIRECTS(...) { \
    void * _stopper = (void * ) -1; \
    mess_direct*  _M[] = { __VA_ARGS__, _stopper}; \
    int _i; \
    for ( _i = 0; _M[_i] != _stopper ; _i++) {\
        mess_direct_init(_M[_i]);\
    }\
}

    /**
     * @brief Macro to clear a list of direct solvers.
     *
     * The @ref MESS_CLEAR_DIRECTS macro allow to clear a list of direct solvers within one call. It uses
     * a variable argument list macro. The calling sequence is
     * \code{.c}
     * MESS_CLEAR_DIRECTS(&direct1, &direct2, &direct3, ...);
     * \endcode
     * */
#define  MESS_CLEAR_DIRECTS(...) { \
    void * _stopper = (void * ) -1; \
    mess_direct*  _M[] = { __VA_ARGS__, _stopper}; \
    int _i; \
    for ( _i = 0; _M[_i] != _stopper ; _i++) {\
        if ( *_M[_i] != NULL) mess_direct_clear(_M[_i]);\
    }\
}

    /**
     * @brief Macro to init a list of @ref mess_equation instances.
     *
     * The @ref MESS_INIT_EQUATIONS macro allow to init  a list of @ref mess_equation instances within one call. It uses
     * a variable argument list macro. The calling sequence is
     * \code{.c}
     * MESS_INIT_EQUATIONS(&eqn1, &eqn2, &eqn3, ...);
     * \endcode
     */
#define  MESS_INIT_EQUATIONS(...) { \
    void * _stopper = (void * ) -1; \
    mess_equation*  _M[] = { __VA_ARGS__, _stopper}; \
    int _i; \
    for ( _i = 0; _M[_i] != _stopper ; _i++) {\
        mess_equation_init(_M[_i]);\
    }\
}

    /**
     * @brief Macro to clear a list of @ref mess_equation instances.
     *
     * The @ref MESS_CLEAR_EQUATIONS macro allow to clear a list of direct solvers within one call. It uses
     * a variable argument list macro. The calling sequence is
     * \code{.c}
     * MESS_CLEAR_EQUATIONS(&eqn1, &eqn2, &eqn3, ...);
     * \endcode
     * */
#define  MESS_CLEAR_EQUATIONS(...) { \
    void * _stopper = (void * ) -1; \
    mess_equation*  _M[] = { __VA_ARGS__, _stopper}; \
    int _i; \
    for ( _i = 0; _M[_i] != _stopper ; _i++) {\
        if ( *_M[_i] != NULL) mess_equation_clear(_M[_i]);\
    }\
}

    /**
     * @brief Macro to clear a list of pointers.
     *
     * The @ref MESS_CLEAR_POINTERS macro allow to free a list of pointers within one call. It uses
     * a variable argument list macro. The calling sequence is
     * \code{.c}
     * MESS_CLEAR_POINTERS(A, B, C, ...);
     * \endcode
     * Before the pointers are cleared it is check if they are already @c NULL.
     * */
#define  MESS_CLEAR_POINTERS(...) { \
    void * _stopper = (void * ) -1; \
    void*  _M[] = { __VA_ARGS__, _stopper}; \
    int _i; \
    for ( _i = 0; _M[_i] != _stopper ; _i++) {\
        if ( _M[_i] != NULL){\
            free(_M[_i]);\
            _M[_i]=NULL;\
        }\
    }\
}

    /**
     * @brief Macro to free memory.
     *
     * The @ref mess_free macro is a wrapper around free which checks for a @c NULL pointer first
     * and sets the freed pointer to @c NULL afterwards.
     *
     **/
#define mess_free(X) do {\
    if ( (X) != NULL ) { \
        __mess_free((X));\
        (X) = NULL; \
    }\
} while (0)

    /**
     * @brief Macro to allocate memory (malloc).
     *
     * The @ref mess_try_alloc macro allocates new memory of possible with C's malloc.
     * Otherwise it returns an error.
     *
     **/
#define mess_try_alloc(mem, cast,  X) { \
    (mem) = (cast) __mess_malloc((X)); \
    if ( (X)>0 &&  (mem) == NULL ) { \
        MSG_ERROR("no memory left to allocate: %s size: %lu\n", #mem, (unsigned long)(X));\
        return (MESS_ERROR_MEMORY) ;\
    } \
}
    /**
     * @brief Macro to allocate memory (malloc).
     *
     * The @ref mess_try_alloc macro allocates new memory of possible with C's malloc. In contrast
     * to \ref mess_try_alloc it does not return an error.
     *
     **/
#define mess_try_alloc2(mem, cast,  X) { \
    (mem) = (cast) __mess_malloc((X)); \
    if ( (X) > 0 && (mem) == NULL ) { \
        MSG_ERROR("no memory left to allocate: %s size: %lu\n", #mem, (unsigned long)(X));\
    } \
}

    /**
     * @brief Macro to allocate memory (calloc).
     *
     * The @ref mess_try_calloc macro allocates new memory of possible with C's calloc.
     * Otherwise it returns an error.
     *
     **/
#define mess_try_calloc(mem, cast, n, size) { \
    (mem) =(cast) __mess_calloc((n),(size));\
    if ( ((n)*(size))>0 && (mem) == NULL) { \
        MSG_ERROR("no memory left to allocate: %s size: %lu\n", #mem, (unsigned long)(size)*(n));\
        return (MESS_ERROR_MEMORY);\
    }\
}

    /**
     * @brief Macro to allocate memory (calloc).
     *
     * The @ref mess_try_calloc macro allocates new memory of possible with C's calloc.
     * In contrast to \ref mess_try_calloc this function does not return an error.
     *
     **/
#define mess_try_calloc2(mem, cast, n, size) { \
    (mem) =(cast) __mess_calloc((n),(size));\
    if ( ((n)*(size))> 0 && (mem) == NULL) { \
        MSG_ERROR("no memory left to allocate: %s size: %lu\n", #mem, (unsigned long)(size)*(n));\
    }\
}
    /**
     * @brief Macro to allocate memory (realloc).
     *
     * The @ref mess_try_realloc macro allocates new memory of possible with C's realloc.
     * Otherwise it returns an error.
     *
     **/
#define mess_try_realloc(mem, cast, X) { \
    (mem) = (cast) __mess_realloc((mem),(X)); \
    if ((X) > 0  && (mem) == NULL){ \
        MSG_ERROR("no memory left to allocate: %s size: %lu\n", #mem,(unsigned long)(X));\
        return (MESS_ERROR_MEMORY);\
    }\
}

    /**
     * @brief Macro to allocate memory (realloc).
     *
     * The @ref mess_try_realloc macro allocates new memory of possible with C's realloc.
     * In contrast to \ref mess_try_realloc this function does not return an error.
     *
     **/
#define mess_try_realloc2(mem, cast, X) { \
    (mem) = (cast) __mess_realloc((mem),(X)); \
    if ( (X) > 0 && (mem) == NULL){ \
        MSG_ERROR("no memory left to allocate: %s size: %lu\n", #mem,(unsigned long)(X));\
    }\
}

    /**
     * @brief Macro to allocate memory (malloc).
     *
     * The @ref mess_try_alloc_omp macro allocates new memory of possible with C's malloc.
     * Otherwise it returns an error.
     *
     **/
#define mess_try_alloc_omp(mem, cast,  X) do { \
    (mem) = (cast) __mess_malloc((X));\
    if ( (X) > 0 &&  (mem) == NULL ) { \
        MSG_ERROR("no memory left to allocate: %s size: %lu\n", #mem, (unsigned long)(X));\
        err = MESS_ERROR_MEMORY; \
    }\
} while (0)

    /**
     * @internal
     * @brief Macro to call a static function (if available) with arguments.
     *
     * The @ref MESS_CALL_IF_EXIST allow to call a static function @p func (if not points to @c NULL)
     * with given list of arguments. If the function @p func points to @c NULL, the function is not called.
     * Please be sure that a variable @p ret of type @ref mess_int_t is predefined for
     * the @ref FUNCTION_FAILURE_HANDLE.
     * @brief Internal use only.
     *
     */
#define MESS_CALL_IF_EXIST(func, args...) \
        if ( func != NULL){     \
            ret = func(args);   \
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), func); \
        }

    /**
     * @brief Macro to check if an object points to @c NULL.
     *
     * Check if an object points to @c NULL. */
#define mess_check_nullpointer(X) if ((X)==NULL) { MSG_ERROR("%s points to NULL\n", #X); return MESS_ERROR_NULLPOINTER; }

    /**
     * @brief Macro to check if an object does not point to @c NULL.
     *
     * Check if an object does not point to @c NULL. */
#define mess_check_notnullpointer(X) if ((X)!=NULL) { MSG_ERROR("%s points not to NULL\n", #X); return MESS_ERROR_NOTNULLPOINTER; }

    /**
     * @brief Macro to check if a matrix or vector is real.
     *
     * Check if a matrix or vector is real. */
#define mess_check_real(X) if (!MESS_IS_REAL((X))) { MSG_ERROR("%s have to be real, current = %s \n", #X, mess_datatype_t_str((X)->data_type)); return MESS_ERROR_DATATYPE; }

    /**
     * @brief Macro to check if a matrix or vector is complex.
     *
     * Check if a matrix or vector is complex. */
#define mess_check_complex(X) if (!MESS_IS_COMPLEX((X))) { MSG_ERROR("%s have to be complex\n", #X); return MESS_ERROR_DATATYPE; }

    /**
     * @brief Macro to check if a matrix is square.
     *
     * Check if a matrix is square. */
#define mess_check_square(X) if ( (X)->rows!=(X)->cols ) { MSG_ERROR("%s have to be square, rows = %ld, cols = %ld\n", #X, (long) (X)->rows, (long) (X)->cols); return MESS_ERROR_DIMENSION;}

    /**
     * @brief Macro to check if a matrix is stored in the CSC storage.
     *
     * Check if a matrix is stored in the CSC storage.*/
#define mess_check_csc(X) if(!MESS_IS_CSC((X))) { MSG_ERROR("%s have to be a CSC stored matrix\n", #X); return MESS_ERROR_STORAGETYPE;}

    /**
     * @brief Macro to check if a matrix is stored in the CSR storage.
     *
     * Check if a matrix is stored in the CSR storage.*/
#define mess_check_csr(X) if(!MESS_IS_CSR((X))) { MSG_ERROR("%s have to be a CSR stored matrix\n", #X); return MESS_ERROR_STORAGETYPE;}

    /**
     * @brief Macro to check if a matrix is stored in the DENSE storage.
     *
     * Check if a matrix is stored in the DENSE storage.*/
#define mess_check_dense(X) if(!MESS_IS_DENSE((X))) { MSG_ERROR("%s have to be a dense matrix\n", #X); return MESS_ERROR_STORAGETYPE;}

    /**
     * @brief Macro to check if a matrix is stored in the sparse storage format.
     *
     * Check if a matrix is stored in the sparse storage format.*/
#define mess_check_sparse(X) if(!(MESS_IS_CSC((X)) || MESS_IS_CSR((X)) || MESS_IS_COORD((X)))) { MSG_ERROR("%s have to be a sparse matrix\n", #X); return MESS_ERROR_STORAGETYPE;}


    /**
     * @brief Macro to check if a matrix is stored in the COORD storage.
     *
     * Check if a matrix  is stored in the COORD storage.*/
#define mess_check_coord(X) if(!MESS_IS_COORD((X))) { MSG_ERROR("%s have to be a coordinate matrix\n", #X); return MESS_ERROR_STORAGETYPE;}

    /**
     * @brief Macro to check if a matrix is not stored in the CSC storage.
     *
     * Check if a matrix is not stored in the CSC storage.*/
#define mess_check_nocoord(X) if(MESS_IS_COORD((X))) { MSG_ERROR("%s is a coordinate matrix. Coordinate matrices are not supported.\n", #X); return MESS_ERROR_STORAGETYPE;}

    /**
     * @brief Macro to check if  an object is positive.
     *
     * Check if an object is positive.*/
#define mess_check_positive(X) if ( (X) <= 0 ) { MSG_ERROR("%s have to be positive\n", #X); return MESS_ERROR_ARGUMENTS;}

    /**
     * @brief Macro to check if  an object is nonnegative.
     *
     * Check if an object is nonnegative.*/
#define mess_check_nonnegative(X) if ( (X) < 0 ) { MSG_ERROR("%s have to be nonnegative\n", #X); return MESS_ERROR_ARGUMENTS;}

    /**
     * @brief Macro to check if two matrices have the same size.
     *
     * Check if two matrices have the same size.*/
#define mess_check_same_size(X,Y) if (( (X)->rows != (Y)->rows) || ( (X)->cols != (Y)->cols)) { MSG_ERROR("%s and %s must have the same size. (" MESS_PRINTF_INT ", " MESS_PRINTF_INT ") <-> (" MESS_PRINTF_INT ", " MESS_PRINTF_INT "))\n",#X,#Y,X->rows,X->cols, Y->rows, Y->cols); return ( MESS_ERROR_DIMENSION); }

    /**
     * @brief Macro to check if two structs have the same dim.
     *
     * Check if two structs have the same dim.*/
#define mess_check_same_dim(X,Y) if (( (X)->dim != (Y)->dim) ) { MSG_ERROR("%s and %s must have the same dim. " MESS_PRINTF_INT " <-> " MESS_PRINTF_INT ")\n",#X,#Y, X->dim, Y->dim); return ( MESS_ERROR_DIMENSION); }

    /**
     * @brief Macro to check if two matrices or vectors have the same number of rows.
     *
     * Check if two matrices or vectors have the same number of rows.*/
#define mess_check_same_rows(X,Y) if (( (X)->rows != (Y)->rows) ) { MSG_ERROR("%s and %s must have the number of rows. (" MESS_PRINTF_INT ", " MESS_PRINTF_INT ") <-> (" MESS_PRINTF_INT ", " MESS_PRINTF_INT "))\n",#X,#Y,X->rows,X->cols, Y->rows, Y->cols); return ( MESS_ERROR_DIMENSION); }

    /**
     * @brief Macro to check if two matrices or vectors have the same number of columns.
     *
     * Check if two matrices or vectors have the same number of columns.*/
#define mess_check_same_cols(X,Y) if (( (X)->cols != (Y)->cols) ) { MSG_ERROR("%s and %s must have the number of columns. (" MESS_PRINTF_INT ", " MESS_PRINTF_INT ") <-> (" MESS_PRINTF_INT ", " MESS_PRINTF_INT "))\n",#X,#Y,X->rows,X->cols, Y->rows, Y->cols); return ( MESS_ERROR_DIMENSION); }

    /**
     * @brief Macro to check if the number of rows of a matrix is identical to the number of columns of another matrix.
     *
     * Check if the number of rows of a matrix is identical to the number of columns of another matrix.*/
#define mess_check_same_colsrows(X,Y) if (( (X)->cols != (Y)->rows) ) { MSG_ERROR("Number of cols of %s must be the same as the number of rows of %s . (" MESS_PRINTF_INT ", " MESS_PRINTF_INT ") <-> (" MESS_PRINTF_INT ", " MESS_PRINTF_INT "))\n",#X,#Y,X->rows,X->cols, Y->rows, Y->cols); return ( MESS_ERROR_DIMENSION); }

    /**
     * @brief Macro to check if the number of columns of a matrix is identical to the dimension of a vector.
     *
     * Check if the number of columns of a matrix is identical to the dimension of a vector.*/
#define mess_check_same_colsdim(X,Y) if (( (X)->cols != (Y)->dim) ) { MSG_ERROR("Number of cols of %s must be the same as the dimension of %s . (" MESS_PRINTF_INT ", " MESS_PRINTF_INT ") <-> (" MESS_PRINTF_INT "))\n",#X,#Y,X->rows,X->cols, Y->dim); return ( MESS_ERROR_DIMENSION); }

    /**
     * @brief Macro to check if the number of rows of a matrix is identical to the dimension of a vector.
     *
     * Check if the number of rows of a matrix is identical to the dimension of a vector.*/
#define mess_check_same_rowsdim(X,Y) if (( (X)->rows != (Y)->dim) ) { MSG_ERROR("Number of rows of %s must be the same as the dimension of %s . (" MESS_PRINTF_INT ", " MESS_PRINTF_INT ") <-> (" MESS_PRINTF_INT "))\n",#X,#Y,X->rows,X->cols, Y->dim); return ( MESS_ERROR_DIMENSION); }

    /**
     * @brief Macro to check if two objects have the same datatype.
     *
     * Check if two objects have the same datatype. */
#define mess_check_same_datatype(X,Y) if (( (X)->data_type != (Y)->data_type) ) { MSG_ERROR("%s and %s must have the same data type.\n",#X,#Y); return( MESS_ERROR_DATATYPE); }

    /**
     * @brief Macro to check if two objects have the same storage type.
     *
     * Check if two objects have the same storage type. */
#define mess_check_same_storetype(X,Y) if (( (X)->store_type != (Y)->store_type) ) { MSG_ERROR("%s and %s must have the same storage type.\n",#X,#Y); return ( MESS_ERROR_STORAGETYPE); }

    /**
     * @brief Macro to check if a matrix or a vector is real or complex.
     *
     * Check if a matrix or a vector is real or complex.*/
#define mess_check_real_or_complex(X) if ( (!MESS_IS_REAL(X) && ! MESS_IS_COMPLEX(X)) ) { MSG_ERROR("%s must be real or complex.\n",#X); return ( MESS_ERROR_DATATYPE); }

    /**
     * @brief Macro to check if an operation on a matrix is correct.
     *
     * Check if an operation on a matrix is correct.*/
#define mess_check_operation_type(X) if ( ((X) != MESS_OP_NONE) && ( (X)!=MESS_OP_TRANSPOSE) && ( (X) != MESS_OP_HERMITIAN )) { \
    MSG_ERROR("The given mess_operation_t is not supported. Invalid Argument\n");\
    return ( MESS_ERROR_ARGUMENTS);  }

    /**
     * @brief Macro to check if a  given storage type is valid.
     *
     * Check if a given storage type is valid.*/
#define mess_check_storage_type(X) if ( ((X) != MESS_DENSE) && ( (X)!=MESS_CSC) && ( (X) != MESS_CSR )  && ( (X) != MESS_COORD )) { \
    MSG_ERROR("The given mess_storage_t is not supported. Invalid Argument\n");\
    return ( MESS_ERROR_ARGUMENTS);  }

    /**
     * @brief Macro to check if an datatype  is correct.
     *
     * Check if an datatype is correct.*/
#define mess_check_datatype(DT) if ( ((DT) != MESS_REAL) && ( (DT)!=MESS_COMPLEX) ) { \
    MSG_ERROR("The given mess_datatype_t is not supported. Invalid Argument\n");\
    return ( MESS_ERROR_DATATYPE);  }

    /**
     * @brief Check if a value is true.
     *
     * The macro checks if a given value is true. In this way it is similar to assert.
     * */
#define mess_check_true(X) if ( (X) == 0 ) { MSG_ERROR("Assertion  %s == true failed.\n", #X); return MESS_ERROR_ARGUMENTS; }

    /** @} */


#ifdef __cplusplus
    }
#endif
#endif
