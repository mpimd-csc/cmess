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
 * @file lib/matrix/sort.c
 * @brief Sort the data structures inside a matrix.
 * @author @koehlerm
 *
 *
 * This file implements the data structure sorting inside the matrix data structures.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>

/*-----------------------------------------------------------------------------
 *  helper stuff for sorting
 *-----------------------------------------------------------------------------*/

/** structure to sort pairs (@ref mess_int_t,@c double) */
struct pair_mess_int_t_double {
    mess_int_t l;
    double d;
};

/** structure to sort pairs (@ref mess_int_t,@c double @c complex) */
struct pair_mess_int_t_complex {
    mess_int_t l;
    mess_double_cpx_t d;
};

/** structure to sort triples (@ref mess_int_t, @ref mess_int_t,@c double) */
struct pair_mess_int_t2_double {
    mess_int_t l1;
    mess_int_t l2;
    double d;
};

/** structure to sort triples (@ref mess_int_t, @ref mess_int_t,@c complex) */
struct pair_mess_int_t2_complex {
    mess_int_t l1;
    mess_int_t l2;
    mess_double_cpx_t d;
};

/** structure to sort pairs (@ref mess_int_t, @ref mess_int_t) */
struct pair_mess_int_t2 {
    mess_int_t l1;
    mess_int_t l2;
};


/**
 * @brief Compare function for qsort (@ref mess_int_t,@c double).
 * @param[in] left  input left side
 * @param[in] right  input right side
 * @return -1 if left<right, 0 if left==right, 1 if left > right
 *
 * The @ref __cmp_mess_int_t_double function is a compare function to
 * sort pairs (@ref mess_int_t, @c double) using lexicographic order on the first entry.
 *
 */
static int __cmp_mess_int_t_double ( const void * left, const void *right)
{
    struct pair_mess_int_t_double *l = (struct pair_mess_int_t_double*) left;
    struct pair_mess_int_t_double *r = (struct pair_mess_int_t_double*) right;

    if ( l->l < r->l){
        return -1;
    } else if ( l->l == r->l) {
        return 0;
    } else {
        return 1;
    }
}

/**
 * @brief Compare function for qsort (@ref mess_int_t, @c complex).
 * @param[in] left  input left side
 * @param[in] right  input right side
 * @return -1 if left<right, 0 if left==right, 1 if left > right
 *
 * The @ref __cmp_mess_int_t_complex function is a compare function to
 * sort pairs (@ref mess_int_t, @c double @c complex) using lexicographic order on the first entry.
 *
 */
static int __cmp_mess_int_t_complex ( const void * left, const void *right)
{
    struct pair_mess_int_t_complex *l = (struct pair_mess_int_t_complex*) left;
    struct pair_mess_int_t_complex *r = (struct pair_mess_int_t_complex*) right;

    if ( l->l < r->l){
        return -1;
    } else if ( l->l == r->l) {
        return 0;
    } else {
        return 1;
    }
}


/**
 * @brief Compare function for qsort (@ref mess_int_t,@ref mess_int_t, @c double).
 * @param[in] left  input left side
 * @param[in] right  input right side
 * @return -1 if left<right, 0 if left==right, 1 if left > right
 *
 * The @ref __cmp_mess_int_t2_double function is a compare function to
 * sort triples (@ref mess_int_t, @ref mess_int_t,@c double) using lexicographic order on the first and second entry.
 *
 */
static int __cmp_mess_int_t2_double ( const void * left, const void *right)
{
    struct pair_mess_int_t2_double *l = (struct pair_mess_int_t2_double*) left;
    struct pair_mess_int_t2_double *r = (struct pair_mess_int_t2_double*) right;

    if ( l->l1 < r->l1) {
        return -1;
    } else  if ( l->l1 == r->l1) {
        if ( l->l2 < r->l2) {
            return -1;
        } else if ( l->l2 == r->l2) {
            return 0;
        } else {
            return 1;
        }
    } else {
        return 1;
    }
}

/**
 * @brief Compare function for qsort (@ref mess_int_t,@ref mess_int_t,@c double @c complex).
 * @param[in] left  input left side
 * @param[in] right  input right side
 * @return -1 if left<right, 0 if left==right, 1 if left > right
 *
 * The @ref __cmp_mess_int_t2_complex function is a compare function to
 * sort triples (mess_int_t, @ref mess_int_t, @c double @c complex).
 *
 */
static int __cmp_mess_int_t2_complex ( const void * left, const void *right)
{
    struct pair_mess_int_t2_complex *l = (struct pair_mess_int_t2_complex*) left;
    struct pair_mess_int_t2_complex *r = (struct pair_mess_int_t2_complex*) right;

    if ( l->l1 < r->l1) {
        return -1;
    } else  if ( l->l1 == r->l1) {
        if ( l->l2 < r->l2) {
            return -1;
        } else if ( l->l2 == r->l2) {
            return 0;
        } else {
            return 1;
        }
    } else {
        return 1;
    }
}


/**
 * @brief Sort the internal structure of a sparse matrix.
 * @param[in,out] mat matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_sort function sorts the internal structure
 * of a sparse matrix. \n
 * It uses the qsort function provided by libc. \n
 * Sorting means that the row or column indices in a compressed
 * sparse matrix appear in ascending order. This is necessary for many algorithms
 * to work correctly or more efficient. In case of a matrix-vector product
 * it helps to get a better data locality.
 *
 * @sa mess_matrix_dupl
 *
 */
int  mess_matrix_sort ( mess_matrix mat )
{
    MSG_FNAME(__func__);
    struct pair_mess_int_t_double *pld = NULL;
    struct pair_mess_int_t_complex *plc = NULL;
    struct pair_mess_int_t2_double *pl2d = NULL;
    struct pair_mess_int_t2_complex *pl2c = NULL;

    mess_int_t i, j;
    mess_int_t pos = 0;

    mess_check_nullpointer (mat);

    if ( MESS_IS_DENSE(mat)) return 0;

    /*-----------------------------------------------------------------------------
     *  sort CSR matrices
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_CSR(mat)){
        if ( MESS_IS_REAL(mat)) {
            mess_try_alloc(pld, struct pair_mess_int_t_double *, sizeof(struct pair_mess_int_t_double)* mat->cols);
            for (i=0; i < mat->rows; i++){
                pos = 0;
                for (j = mat->rowptr[i]; j < mat->rowptr[i+1]; j++){
                    pld[pos].l = mat->colptr[j];
                    pld[pos].d = mat->values[j];
                    pos++;
                }
                qsort(pld, pos, sizeof(struct pair_mess_int_t_double), __cmp_mess_int_t_double);
                pos = 0;
                for (j = mat->rowptr[i]; j < mat->rowptr[i+1]; j++){
                    mat->colptr[j] = pld[pos].l;
                    mat->values[j] = pld[pos].d;
                    pos++;
                }
            }
            mess_free(pld);
        } else if ( MESS_IS_COMPLEX(mat)){
            mess_try_alloc(plc, struct pair_mess_int_t_complex *, sizeof(struct pair_mess_int_t_complex)* mat->cols);
            for (i=0; i < mat->rows; i++){
                pos = 0;
                for (j = mat->rowptr[i]; j < mat->rowptr[i+1]; j++){
                    plc[pos].l = mat->colptr[j];
                    plc[pos].d = mat->values_cpx[j];
                    pos++;
                }
                qsort(plc, pos, sizeof(struct pair_mess_int_t_complex), __cmp_mess_int_t_complex);
                pos = 0;
                for (j = mat->rowptr[i]; j < mat->rowptr[i+1]; j++){
                    mat->colptr[j]     = plc[pos].l;
                    mat->values_cpx[j] = plc[pos].d;
                    pos++;
                }
            }
            mess_free(plc);
        } else {
            MSG_ERROR("data type unknown\n");
            return MESS_ERROR_DATATYPE;
        }
    }
    /*-----------------------------------------------------------------------------
     *  sort CSC matrices
     *-----------------------------------------------------------------------------*/
    else if ( MESS_IS_CSC(mat)){
        if ( MESS_IS_REAL(mat)) {
            mess_try_alloc(pld, struct pair_mess_int_t_double *, sizeof(struct pair_mess_int_t_double)* mat->rows);
            for (i=0; i < mat->cols; i++){
                pos = 0;
                for (j = mat->colptr[i]; j < mat->colptr[i+1]; j++){
                    pld[pos].l = mat->rowptr[j];
                    pld[pos].d = mat->values[j];
                    pos++;
                }
                qsort(pld, pos, sizeof(struct pair_mess_int_t_double), __cmp_mess_int_t_double);
                pos = 0;
                for (j = mat->colptr[i]; j < mat->colptr[i+1]; j++){
                    mat->rowptr[j] = pld[pos].l;
                    mat->values[j] = pld[pos].d;
                    pos++;
                }
            }
            mess_free(pld);
        } else if ( MESS_IS_COMPLEX(mat)){
            mess_try_alloc(plc, struct pair_mess_int_t_complex *, sizeof(struct pair_mess_int_t_complex)* mat->rows);
            for (i=0; i < mat->cols; i++){
                pos = 0;
                for (j = mat->colptr[i]; j < mat->colptr[i+1]; j++){
                    plc[pos].l = mat->rowptr[j];
                    plc[pos].d = mat->values_cpx[j];
                    pos++;
                }
                qsort(plc, pos, sizeof(struct pair_mess_int_t_complex), __cmp_mess_int_t_complex);
                pos = 0;
                for (j = mat->colptr[i]; j < mat->colptr[i+1]; j++){
                    mat->rowptr[j]     = plc[pos].l;
                    mat->values_cpx[j] = plc[pos].d;
                    pos++;
                }
            }
            mess_free(plc);
        } else {
            MSG_ERROR("data type unknown\n");
            return MESS_ERROR_DATATYPE;
        }
    }

    /*-----------------------------------------------------------------------------
     *  sort COORD
     *-----------------------------------------------------------------------------*/
    else if ( MESS_IS_COORD(mat)) {
        if ( MESS_IS_REAL(mat)) {
            mess_try_alloc(pl2d, struct pair_mess_int_t2_double *, sizeof(struct pair_mess_int_t2_double)* mat->nnz);
            pos  = 0;
            for ( i = 0; i< mat->nnz ; i++)  {
                pl2d[pos].l1 = mat->rowptr [i];
                pl2d[pos].l2 = mat->colptr [i];
                pl2d[pos].d  = mat->values [i];
                pos++;
            }
            qsort(pl2d, pos, sizeof(struct pair_mess_int_t2_double), __cmp_mess_int_t2_double);
            pos = 0;
            for ( i = 0; i< mat->nnz ; i++)  {
                mat->rowptr [i] = pl2d[pos].l1;
                mat->colptr [i] = pl2d[pos].l2;
                mat->values [i] = pl2d[pos].d;
                pos++;
            }
            mess_free(pl2d);
        } else if ( MESS_IS_COMPLEX(mat)) {
            mess_try_alloc(pl2c, struct pair_mess_int_t2_complex *, sizeof(struct pair_mess_int_t2_complex)* mat->nnz);
            pos  = 0;
            for ( i = 0; i< mat->nnz ; i++)  {
                pl2c[pos].l1 = mat->rowptr [i];
                pl2c[pos].l2 = mat->colptr [i];
                pl2c[pos].d  = mat->values_cpx [i];
                pos++;
            }
            qsort(pl2c, pos, sizeof(struct pair_mess_int_t2_complex), __cmp_mess_int_t2_complex);
            pos = 0;
            for ( i = 0; i< mat->nnz ; i++)  {
                mat->rowptr [i] = pl2c[pos].l1;
                mat->colptr [i] = pl2c[pos].l2;
                mat->values_cpx [i] = pl2c[pos].d;
                pos++;
            }
            mess_free(pl2c);
        } else {
            MSG_ERROR("data type unknown\n");
            return MESS_ERROR_DATATYPE;
        }
    } else {
        MSG_ERROR("Storage type isn't supported at the moment.\n");
        return MESS_ERROR_STORAGETYPE;
    }

    return 0;
}       /* -----  end of function mess_matrix_sort  ----- */

