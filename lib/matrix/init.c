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
 * @file lib/matrix/init.c
 * @brief Basic management of mess_matrix structures.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <stddef.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>


/**
 * @brief Remove a mess_matrix object from memory.
 * @param[in,out] matrix pointer to the mess_matrix structure
 * @return always zero (if @mess is compiled using the DEBUG option then it returns an error if the matrix
 * is already freed)
 *
 * The @ref mess_matrix_clear function removes the mess_matrix and cleans up. \n
 * If the matrix is a view on to another one only the helper data is freed. The data is freed when the parent matrix
 * is freed. \n
 * In contrast to other functions the @ref mess_matrix is passed as a pointer
 * instead of a direct passing.
 */
int mess_matrix_clear(mess_matrix *matrix)
{
#ifdef MESS_DEBUG
    MSG_FNAME(__func__);
#endif
    if (matrix  == NULL){
#ifdef MESS_DEBUG
        MSG_INFO("matrix already points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
#else
        return 0;
#endif
    }
    if (*matrix == NULL) {
#ifdef MESS_DEBUG
        MSG_INFO("*matrix already points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
#else
        return 0;
#endif
    }
    if ((*matrix)->colptr != NULL )mess_free((*matrix)->colptr);
    if ((*matrix)->rowptr != NULL )mess_free((*matrix)->rowptr);
    if ((*matrix)->values != NULL )mess_free((*matrix)->values);
    if ((*matrix)->values_cpx != NULL )mess_free((*matrix)->values_cpx);


    mess_free(*matrix);
    *matrix = NULL;
    return 0;
}

/**
 * @brief Initialize the mess_matrix structure.
 * @param[in,out] matrix a pointer to a mess_matrix object
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_init function creates a new mess_matrix object. \n
 * It allocates the basic data structures and sets all values to 0. \n
 * The array for the matrix data is allocated. This have to be done
 * separately by \ref mess_matrix_alloc. \n
 * In contrast to other functions the @ref mess_matrix is passed as a pointer
 * instead of a direct passing.
 *
 */
int mess_matrix_init(mess_matrix *matrix)
{
    MSG_FNAME(__func__);
    mess_try_alloc(*matrix , struct mess_matrix_st*,sizeof(struct mess_matrix_st));
    memset(*matrix, 0, sizeof(struct mess_matrix_st));
    (*matrix)->cols  = 0;
    (*matrix)->rows  = 0;
    (*matrix)->ld  = 0;
    (*matrix)->data_type  = 0;
    (*matrix)->nnz  = 0;
    (*matrix)->store_type = MESS_UNKNOWN;
    (*matrix)->colptr = NULL;
    (*matrix)->rowptr = NULL;
    (*matrix)->values = NULL;
    (*matrix)->values_cpx = NULL;

    return 0;
}

/**
 * @brief Copy a matrix to another one.
 * @param[in] in    input source matrix
 * @param[out] out   destination matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_copy function copies a matrix from in to out. \n
 * If a view of a matrix is copied, the leading dimension changes to the
 * default of \ref mess_matrix_alloc.
 *
 *
 */
int mess_matrix_copy(mess_matrix in, mess_matrix out) {
    MSG_FNAME(__func__);

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    MESS_MATRIX_RESET(out);
    memset(out, 0, sizeof(struct mess_matrix_st));
    memcpy(out, in, sizeof ( struct mess_matrix_st));
    /* Set out->{colptr, rowptr, values, values_cpx} to NULL because otherwise the original data is freed by mess_matrix_alloc */
    out->colptr = NULL;
    out->rowptr = NULL;
    out->values = NULL;
    out->values_cpx = NULL;


    if ( MESS_IS_CSR(in)) {
        // Copy CSR
        mess_try_alloc(out->rowptr, mess_int_t*, sizeof(mess_int_t)*(in->rows+1));
        mess_try_alloc(out->colptr, mess_int_t*, sizeof(mess_int_t)*(in->nnz));

        memcpy(out->rowptr, in->rowptr, (in->rows+1)*sizeof(mess_int_t));
        memcpy(out->colptr, in->colptr, (in->nnz)*sizeof(mess_int_t));
        if ( MESS_IS_REAL(in)){
            mess_try_alloc(out->values, double*, (in->nnz) * sizeof(double));
            memcpy(out->values, in->values, (in->nnz)*sizeof(double));
            out->data_type = MESS_REAL;
        }
        if ( MESS_IS_COMPLEX(in)){
            mess_try_alloc(out->values_cpx, mess_double_cpx_t*, (in->nnz) * sizeof(mess_double_cpx_t ));
            memcpy(out->values_cpx, in->values_cpx, (in->nnz)*sizeof(mess_double_cpx_t));
            out->data_type= MESS_COMPLEX;
        }
    } else   if ( MESS_IS_CSC(in)) {
        // Copy CSC
        mess_try_alloc(out->rowptr, mess_int_t*, sizeof(mess_int_t)*(in->nnz));
        mess_try_alloc(out->colptr, mess_int_t*, sizeof(mess_int_t)*(in->cols+1));

        memcpy(out->rowptr, in->rowptr, (in->nnz)*sizeof(mess_int_t));
        memcpy(out->colptr, in->colptr, (in->cols+1)*sizeof(mess_int_t));
        if ( MESS_IS_REAL(in)){
            mess_try_alloc(out->values, double*, (in->nnz) * sizeof(double));
            memcpy(out->values, in->values, (in->nnz)*sizeof(double));
            out->data_type = MESS_REAL;
        }
        if ( MESS_IS_COMPLEX(in)){
            mess_try_alloc(out->values_cpx, mess_double_cpx_t*, (in->nnz) * sizeof(mess_double_cpx_t ));
            memcpy(out->values_cpx, in->values_cpx, (in->nnz)*sizeof(mess_double_cpx_t));
            out->data_type= MESS_COMPLEX;
        }
    } else   if ( MESS_IS_COORD(in)) {
        // COORD
        mess_try_alloc(out->rowptr, mess_int_t*, sizeof(mess_int_t)*(in->nnz));
        mess_try_alloc(out->colptr, mess_int_t*, sizeof(mess_int_t)*(in->nnz));
        memcpy(out->rowptr, in->rowptr, (in->nnz)*sizeof(mess_int_t));
        memcpy(out->colptr, in->colptr, (in->nnz)*sizeof(mess_int_t));
        if ( MESS_IS_REAL(in)){
            mess_try_alloc(out->values, double*, sizeof(double)*(in->nnz));
            memcpy(out->values, in->values, (in->nnz)*sizeof(double));
        }
        if ( MESS_IS_COMPLEX(in)){
            mess_try_alloc(out->values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*(in->nnz));
            memcpy(out->values_cpx, in->values_cpx, (in->nnz)*sizeof(mess_double_cpx_t));
        }
    } else if ( MESS_IS_DENSE(in)) {
        // DENSE
        //ret = mess_matrix_alloc(out, in->rows, in->cols, in->nnz, in->store_type, in->data_type);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        if ( MESS_IS_REAL(in)){
            mess_try_alloc(out->values, double*, sizeof(double)*(out->ld*out->cols));
            F77_GLOBAL(dlacpy,DLACPY)("A",&in->rows,&in->cols,in->values, &in->ld, out->values, &out->ld);
        } else   if ( MESS_IS_COMPLEX(in)){
            mess_try_alloc(out->values_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*(out->ld*out->cols));
            F77_GLOBAL(zlacpy,ZLACPY)("A",&in->rows,&in->cols,in->values_cpx, &in->ld, out->values_cpx, &out->ld);
        }
    } else {
        MSG_ERROR("unknown/unsupported storage type\n");
        return MESS_ERROR_STORAGETYPE;
    }
    return 0;
}

/**
 * @brief Return the size of a matrix in bytes.
 * @param[in] matrix input matrix
 * @return approximate size of the input matrix in bytes.
 *
 * The @ref mess_matrix_memsize function determines the size of a matrix in bytes. \n
 * If the matrix is a view the returned size is 0, because the effectively
 * the matrix needs no additional memory.
 *
 * @sa mess_matrix_memsize_nnz
 *
 */
mess_int_t mess_matrix_memsize(mess_matrix matrix){
    MSG_FNAME(__func__);
    mess_int_t  total = 0;
    mess_int_t sval =0;

    mess_check_nullpointer(matrix);


    if ( MESS_IS_COMPLEX(matrix)) {
        sval = sizeof(mess_double_cpx_t);
    } else if ( MESS_IS_REAL(matrix)){
        sval = sizeof(double);
    } else {
        sval =0;
    }

    if ( MESS_IS_CSR(matrix)) {
        total = matrix->nnz * (sval + sizeof(mess_int_t)) + (matrix->rows+1)*sizeof(mess_int_t);
    } else if ( MESS_IS_CSC(matrix)) {
        total = matrix->nnz * (sval + sizeof(mess_int_t)) + (matrix->cols+1)*sizeof(mess_int_t);
    } else if ( MESS_IS_COORD(matrix)) {
        total = matrix->nnz * (sval + 2*sizeof(mess_int_t));
    } else if ( MESS_IS_DENSE(matrix)) {
        total = matrix->cols*matrix->rows * sval;
    } else {
        MSG_WARN("can't determine size. wrong storage size.");
        total = 0;
    }
    return total;
}

/**
 * @brief Determine memory size of a matrix from parameters.
 * @param[in] rows   input number of rows
 * @param[in] cols    input number of cols
 * @param[in] nnz   input number of nonzero elements
 * @param[in] type   input storage type
 * @param[in] dtype   input data type of the matrix
 * @return  approximate size of the desired matrix
 *
 * The @ref mess_matrix_memsize_nnz function determines the size of a matrix in byte
 * from given parameters.
 *
 * @sa mess_matrix_memsize
 *
 */
mess_int_t mess_matrix_memsize_nnz( mess_int_t rows, mess_int_t cols, mess_int_t nnz,
        mess_int_t type, mess_int_t dtype)
{
    MSG_FNAME(__func__);
    mess_int_t  total = 0;
    mess_int_t  esize;

    if ( dtype == MESS_REAL){
        esize = sizeof(double);
    } else if ( dtype == MESS_COMPLEX) {
        esize = sizeof(mess_double_cpx_t);
    } else {
        esize = 0;
    }

    if (type == MESS_CSR ) {
        total = nnz * (esize + sizeof(mess_int_t)) + (rows+1)*sizeof(mess_int_t);
    } else if (type == MESS_CSC) {
        total = nnz * (esize + sizeof(mess_int_t)) + (cols+1)*sizeof(mess_int_t);
    } else if ( type == MESS_COORD) {
        total = nnz * (esize + 2*sizeof(mess_int_t));
    } else if ( type == MESS_DENSE) {
        total = cols* rows * esize;
    } else {
        MSG_WARN("can't determine size. wrong storage size.");
        total = 0;
    }
    return total;
}

/**
 * @brief Change the data type of a matrix to complex.
 * @param[in,out] m matrix to be changed to complex
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_tocomplex function converts an arbitrary matrix to a complex
 * matrix. \n
 * If the matrix is already complex nothing will change. \n
 * If the matrix is real the whole values array will be copied to a complex array. \n
 * If the matrix is a view on to an other matrix, the parent matrix
 * will be converted to a complex one.
 *
 */
int mess_matrix_tocomplex(mess_matrix m) {
    MSG_FNAME(__func__);
    mess_int_t i;
    mess_check_nullpointer(m);

    if ( MESS_IS_COMPLEX(m)) return 0;

    if (MESS_IS_DENSE(m)){
        mess_try_alloc(m->values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * m->cols*m->ld);
        for ( i = 0; i < m->cols*m->ld; i++) m->values_cpx[i] = m->values[i] + 0*i;
        mess_free( m->values);
        m->values = NULL;
        m->data_type = MESS_COMPLEX;
    } else {
        mess_try_alloc(m->values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * m->nnz);
        for ( i = 0; i < m->nnz; i++) m->values_cpx[i] = m->values[i] + 0*i;
        mess_free( m->values);
        m->values = NULL;
        m->data_type = MESS_COMPLEX;
    }
    return 0;
}

/**
 * @brief Change the data type of a matrix to real.
 * @param[in,out] m matrix to be changed to real
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_toreal function converts an arbitrary matrix to a real
 * matrix. \n
 * If the matrix is already real nothing will change. \n
 * If the matrix is complex all imaginary parts are cut of. \n
 * If the matrix is a view on to an other matrix, the parent matrix
 * will be converted to a complex one.
 *
 */
int mess_matrix_toreal(mess_matrix m) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i;
    mess_check_nullpointer(m);

    if ( MESS_IS_REAL(m)) return 0;

    if (MESS_IS_DENSE(m)) {
        mess_try_alloc(m->values, double *, sizeof(double) * m->cols*m->ld);
        for ( i = 0; i < m->cols*m->ld; i++) m->values[i] = creal(m->values_cpx[i]);
        mess_free( m->values_cpx);
        m->values_cpx = NULL;
        m->data_type = MESS_REAL;
    } else {
        mess_try_alloc(m->values, double *, sizeof(double) * m->nnz);
        for ( i = 0; i < m->nnz; i++) m->values[i] = creal(m->values_cpx[i]);
        mess_free( m->values_cpx);
        m->values_cpx = NULL;
        m->data_type = MESS_REAL;
        ret = mess_matrix_eliminate_zeros( m );        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eliminate_zeros);
    }
    return 0;
}


/**
 * @brief Convert a matrix to a specified data type.
 * @param[in,out]  mat   matrix to be converted
 * @param[in] dt   input desired output type
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_totype function converts a matrix to
 * a specified data type.\n
 * It is a wrapper around \ref mess_matrix_toreal and \ref mess_matrix_tocomplex.
 *
 */
int  mess_matrix_totype ( mess_matrix mat, mess_datatype_t dt )
{
    MSG_FNAME(__func__);

    mess_check_nullpointer(mat);
    mess_check_real_or_complex(mat);

    if ( dt == mat->data_type) return 0;
    switch(dt){
        case MESS_REAL:
            return mess_matrix_toreal(mat);
            break;
        case MESS_COMPLEX:
            return mess_matrix_tocomplex(mat);
            break;
        default:
            MSG_ERROR("Unsupported/unknown data type\n");
            return MESS_ERROR_DATATYPE;
    }
    return 0;
}    /* -----  end of function mess_matrix_totype  ----- */






