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
 * @file lib/matrix/setelement.c
 * @brief Set an element in a matrix.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>

/** structure to sort pairs (mess_int_t, double) */
struct pair_mess_int_t_double {
    mess_int_t l;
    double d;
};

/** structure to sort pairs (mess_int_t, complex) */
struct pair_mess_int_t_complex {
    mess_int_t l;
    mess_double_cpx_t d;
};

/**
 * @brief Compare function for qsort (mess_int_t,double).
 * @param[in] left    input left hand side
 * @param[in] right   input right hand side
 * @return -1, 0 or 1 to indicate if left<right, left == right or left > right.
 *
 * The @ref __cmp_mess_int_t_double function is a compare function to
 * sort pairs (mess_int_t, double) using qsort.
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
 * @brief Compare function for qsort (mess_int_t,complex).
 * @param[in] left    input left hand side
 * @param[in] right   input right hand side
 * @return -1, 0 or 1 to indicate if left<right, left == right or left > right.
 *
 * The @ref __cmp_mess_int_t_complex function is a compare function to
 * sort pairs (mess_int_t, complex) using qsort.
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
 * @brief Set a complex entry in a matrix.
 * @param[in,out] matrix    matrix to set a value
 * @param[in] row    input row index
 * @param[in] col    input column index
 * @param[in] value  input value to set at \f$ (row, col) \f$
 * @return zero on success or a non zero error code otherwise

 *
 * The @ref mess_matrix_setelement_complex function sets a value in a matrix.\n
 * This function should not be used to  assemble a sparse matrix because of performance problems. \n
 * The input need to be a complex matrix. Otherwise it has to be converted to a complex one before.
 *
 * @sa mess_matrix_setelement
 * @sa mess_matrix_getelement
 *
 */
int mess_matrix_setelement_complex ( mess_matrix matrix, mess_int_t row, mess_int_t col, mess_double_cpx_t value )
{
    MSG_FNAME(__func__);
    mess_int_t pos;
    mess_check_nullpointer(matrix);
    mess_check_complex(matrix);


    if ( row >= matrix->rows || row < 0 ) {
        MSG_ERROR("row index out of range row = " MESS_PRINTF_INT " , matrix->rows=" MESS_PRINTF_INT "\n", row, matrix->rows);
        return(MESS_ERROR_DIMENSION);
    }
    if ( col >= matrix->cols || col < 0 ) {
        MSG_ERROR("column index out of range col = " MESS_PRINTF_INT " , matrix->cols=" MESS_PRINTF_INT "\n", col, matrix->cols);
        return(MESS_ERROR_DIMENSION);
    }

    if ( MESS_IS_DENSE(matrix)) {
        matrix->values_cpx[row+col*matrix->ld]=value;
    } else if ( MESS_IS_CSR(matrix)) {
        if ( value != 0.0 ) {
            mess_int_t i;
            int found = 0;
            for ( i= matrix->rowptr[row] ; i < matrix->rowptr[row+1]; i++){
                if ( matrix->colptr[i] == col ) {
                    matrix->values_cpx[i] = value;
                    found = 1;
                }
            }
            if (!found) {

                mess_int_t elementstocopy;
                struct pair_mess_int_t_complex *pld=NULL;
                mess_int_t ii,j;
                mess_try_realloc(matrix->values_cpx, mess_double_cpx_t*, (matrix->nnz+1) * sizeof(mess_double_cpx_t));
                mess_try_realloc(matrix->colptr, mess_int_t *, (matrix->nnz+1) * sizeof(mess_int_t));
                matrix->values_cpx[matrix->nnz] = 0;
                matrix->colptr[matrix->nnz] = 0;
                elementstocopy = matrix->nnz-matrix->rowptr[row];
                if ( matrix -> nnz > 0 ) {
                    memmove(&(matrix->values_cpx[matrix->rowptr[row]+1]), &(matrix->values_cpx[matrix->rowptr[row]]), elementstocopy * sizeof(mess_double_cpx_t));
                    memmove(&(matrix->colptr[matrix->rowptr[row]+1]), &(matrix->colptr[matrix->rowptr[row]]), elementstocopy * sizeof(mess_int_t));
                }
                for (ii = row+1; ii <= matrix->rows; ii++) matrix->rowptr[ii]++;
                matrix->values_cpx[matrix->rowptr[row]]=value;
                matrix->colptr[matrix->rowptr[row]]=col;
                matrix->nnz++;

                mess_try_alloc(pld, struct pair_mess_int_t_complex *, sizeof(struct pair_mess_int_t_complex)* matrix->cols);
                pos = 0;
                for (j = matrix->rowptr[row]; j < matrix->rowptr[row+1]; j++){
                    pld[pos].l = matrix->colptr[j];
                    pld[pos].d = matrix->values_cpx[j];
                    pos++;
                }
                qsort(pld, pos, sizeof(struct pair_mess_int_t_complex), __cmp_mess_int_t_complex);
                pos = 0;
                for (j = matrix->rowptr[row]; j < matrix->rowptr[row+1]; j++){
                    matrix->colptr[j] = pld[pos].l;
                    matrix->values_cpx[j] = pld[pos].d;
                    pos++;
                }
                mess_free(pld);
            }
        } else {
            // remove the element
            if ( matrix->nnz < 1 ) return(0);
            mess_int_t i, fp=-1;
            mess_int_t elementstocopy;
            for (i=matrix->rowptr[row]; i<matrix->rowptr[row+1]; i++) {
                if ( matrix->colptr[i]==col) {fp = i; break; }
            }
            if ( fp < 0 ) return(0);
            elementstocopy = matrix->nnz - fp - 1;
            memmove(&(matrix->values_cpx[fp]),&(matrix->values_cpx[fp+1]), sizeof(mess_double_cpx_t)*elementstocopy) ;
            memmove(&(matrix->colptr[fp]),&(matrix->colptr[fp+1]), sizeof(mess_int_t)*elementstocopy);
            matrix->nnz--;
            mess_try_realloc(matrix->values_cpx, mess_double_cpx_t* , sizeof(mess_double_cpx_t) * matrix->nnz);
            mess_try_realloc(matrix->colptr, mess_int_t* , sizeof(mess_int_t) * matrix->nnz);

            for ( i = row +1; i <= matrix->rows; i++) matrix->rowptr[i]--;
        }
    }  else if ( MESS_IS_CSC(matrix)) {
        if ( value != 0.0 ) {
            mess_int_t i;
            int found = 0;
            for ( i= matrix->colptr[col] ; i < matrix->colptr[col+1]; i++){
                if ( matrix->rowptr[i] == row ) {
                    matrix->values_cpx[i] = value;
                    found = 1;
                }
            }
            if (!found) {
                mess_int_t elementstocopy;
                struct pair_mess_int_t_complex *pld=NULL;
                mess_int_t j;
                mess_try_realloc(matrix->values_cpx, mess_double_cpx_t *, (matrix->nnz+1) * sizeof(mess_double_cpx_t));
                mess_try_realloc(matrix->rowptr, mess_int_t *, (matrix->nnz+1) * sizeof(mess_int_t));
                elementstocopy = matrix->nnz-matrix->colptr[col];
                if ( matrix->nnz > 0 ) {
                    memmove(&(matrix->values_cpx[matrix->colptr[col]+1]), &(matrix->values_cpx[matrix->colptr[col]]), elementstocopy * sizeof(mess_double_cpx_t));
                    memmove(&(matrix->rowptr[matrix->colptr[col]+1]), &(matrix->rowptr[matrix->colptr[col]]), elementstocopy * sizeof(mess_int_t));
                }
                for (i = col+1; i <= matrix->cols; i++) matrix->colptr[i]++;
                matrix->values_cpx[matrix->colptr[col]]=value;
                matrix->rowptr[matrix->colptr[col]]=row;
                matrix->nnz++;

                mess_try_alloc(pld, struct pair_mess_int_t_complex *, sizeof(struct pair_mess_int_t_complex)* matrix->rows);
                pos = 0;
                for (j = matrix->colptr[col]; j < matrix->colptr[col+1]; j++){
                    pld[pos].l = matrix->rowptr[j];
                    pld[pos].d = matrix->values_cpx[j];
                    pos++;
                }
                qsort(pld, pos, sizeof(struct pair_mess_int_t_complex), __cmp_mess_int_t_complex);
                pos = 0;
                for (j = matrix->colptr[col]; j < matrix->colptr[col+1]; j++){
                    matrix->rowptr[j] = pld[pos].l;
                    matrix->values_cpx[j] = pld[pos].d;
                    pos++;
                }
                mess_free(pld);
            }
        } else {
            if ( matrix->nnz < 1 ) return(0);
            // remove the element
            mess_int_t i, fp=-1;
            mess_int_t elementstocopy;
            for (i=matrix->colptr[col]; i<matrix->colptr[col+1]; i++) {
                if ( matrix->rowptr[i]==row) {fp = i; break; }
            }
            if ( fp < 0 ) return(0);
            elementstocopy = matrix->nnz - fp - 1;
            memmove(&(matrix->values_cpx[fp]),&(matrix->values_cpx[fp+1]), sizeof(mess_double_cpx_t)*elementstocopy) ;
            memmove(&(matrix->rowptr[fp]),&(matrix->rowptr[fp+1]), sizeof(mess_int_t)*elementstocopy);
            matrix->nnz--;
            mess_try_realloc(matrix->values_cpx, mess_double_cpx_t* , sizeof(mess_double_cpx_t) * matrix->nnz);
            mess_try_realloc(matrix->rowptr, mess_int_t* , sizeof(mess_int_t) * matrix->nnz);
            for ( i = col +1; i <= matrix->cols; i++) matrix->colptr[i]--;
        }
    } else if ( MESS_IS_COORD( matrix  )) {
        if ( value !=0.0 ) {
            int found = 0;
            mess_int_t i ;
            for ( i = 0 ; i< matrix->nnz; i++) {
                if ( matrix->rowptr[i] == row && matrix->colptr[i] == col) { matrix->values_cpx[i] = value; found = 1;  break; }
            }
            if ( !found ) {
                matrix->nnz++;
                mess_try_realloc(matrix->rowptr, mess_int_t *, sizeof(mess_int_t) * matrix->nnz);
                mess_try_realloc(matrix->colptr, mess_int_t *, sizeof(mess_int_t) * matrix->nnz);
                mess_try_realloc(matrix->values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * matrix->nnz);
                matrix->rowptr[matrix->nnz-1] = row;
                matrix->colptr[matrix->nnz-1] = col;
                matrix->values_cpx[matrix->nnz-1] = value;
            }
        } else {
            if ( matrix->nnz < 1) return (0);
            mess_int_t i , fp = -1;
            mess_int_t elementstocopy;
            for ( i = 0 ; i< matrix->nnz; i++) {
                if ( matrix->rowptr[i] == row && matrix->colptr[i] == col) { fp = i ;  break; }
            }
            if ( fp < 0 ) return(0);
            elementstocopy = matrix->nnz - fp - 1;
            memmove(&(matrix->rowptr[fp]), &(matrix->rowptr[fp+1]), elementstocopy*sizeof(mess_int_t));
            memmove(&(matrix->colptr[fp]), &(matrix->colptr[fp+1]), elementstocopy*sizeof(mess_int_t));
            memmove(&(matrix->values_cpx[fp]), &(matrix->values_cpx[fp+1]), elementstocopy*sizeof(mess_double_cpx_t));
            matrix->nnz--;
        }
    }

    return(0);
}       /* -----  end of function mess_matrix_setelement  ----- */


/**
 * @brief Set a real entry in a matrix.
 * @param[in,out] matrix    matrix to set a value
 * @param[in] row    input row index
 * @param[in] col    input column index
 * @param[in] value  input value to set at \f$ (row, col) \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_setelement function sets a value in a matrix. \n
 * This function should not be used to assemble a sparse matrix because of performance problems. \n
 *
 * @sa mess_matrix_setelement_complex
 * @sa mess_matrix_getelement
 *
 */
int mess_matrix_setelement ( mess_matrix matrix, mess_int_t row, mess_int_t col, double value )
{
    MSG_FNAME(__func__);
    mess_int_t pos;
    mess_check_nullpointer(matrix);
    mess_check_real_or_complex(matrix);

    if(MESS_IS_COMPLEX(matrix)){
        return mess_matrix_setelement_complex(matrix, row, col, value);
    }

    if ( row >= matrix->rows || row < 0 ) {
        MSG_ERROR("row index out of range row = " MESS_PRINTF_INT " , matrix->rows=" MESS_PRINTF_INT "\n", row, matrix->rows);
        return(MESS_ERROR_DIMENSION);
    }
    if ( col >= matrix->cols || col < 0 ) {
        MSG_ERROR("column index out of range col = " MESS_PRINTF_INT " , matrix->cols=" MESS_PRINTF_INT "\n", col, matrix->cols);
        return(MESS_ERROR_DIMENSION);
    }
    if ( MESS_IS_DENSE(matrix)) {
        matrix->values[row+col*matrix->ld]=value;
    } else if ( MESS_IS_CSR(matrix)) {
        if ( value != 0.0 ) {
            mess_int_t i;
            int found = 0;
            for ( i= matrix->rowptr[row] ; i < matrix->rowptr[row+1]; i++){
                if ( matrix->colptr[i] == col ) {
                    matrix->values[i] = value;
                    found = 1;
                }
            }
            if (!found){

                mess_int_t elementstocopy;
                struct pair_mess_int_t_double *pld=NULL;
                mess_int_t ii,j;
                mess_try_realloc(matrix->values, double *, (matrix->nnz+1) * sizeof(double));
                mess_try_realloc(matrix->colptr, mess_int_t *, (matrix->nnz+1) * sizeof(mess_int_t));
                matrix->values[matrix->nnz] = 0;
                matrix->colptr[matrix->nnz] = 0;
                elementstocopy = matrix->nnz-matrix->rowptr[row];
                if ( matrix -> nnz > 0 ) {
                    memmove(&(matrix->values[matrix->rowptr[row]+1]), &(matrix->values[matrix->rowptr[row]]), elementstocopy * sizeof(double));
                    memmove(&(matrix->colptr[matrix->rowptr[row]+1]), &(matrix->colptr[matrix->rowptr[row]]), elementstocopy * sizeof(mess_int_t));
                }
                for (ii = row+1; ii <= matrix->rows; ii++) matrix->rowptr[ii]++;
                matrix->values[matrix->rowptr[row]]=value;
                matrix->colptr[matrix->rowptr[row]]=col;
                matrix->nnz++;

                mess_try_alloc(pld, struct pair_mess_int_t_double *, sizeof(struct pair_mess_int_t_double)* matrix->cols);
                pos = 0;
                for (j = matrix->rowptr[row]; j < matrix->rowptr[row+1]; j++){
                    pld[pos].l = matrix->colptr[j];
                    pld[pos].d = matrix->values[j];
                    pos++;
                }
                qsort(pld, pos, sizeof(struct pair_mess_int_t_double), __cmp_mess_int_t_double);
                pos = 0;
                for (j = matrix->rowptr[row]; j < matrix->rowptr[row+1]; j++){
                    matrix->colptr[j] = pld[pos].l;
                    matrix->values[j] = pld[pos].d;
                    pos++;
                }
                mess_free(pld);
            }
        } else {
            // remove the element
            if ( matrix->nnz < 1 ) return(0);
            mess_int_t i, fp=-1;
            mess_int_t elementstocopy;
            for (i=matrix->rowptr[row]; i<matrix->rowptr[row+1]; i++) {
                if ( matrix->colptr[i]==col) {fp = i; break; }
            }
            if ( fp < 0 ) return(0);
            elementstocopy = matrix->nnz - fp - 1;
            memmove(&(matrix->values[fp]),&(matrix->values[fp+1]), sizeof(double)*elementstocopy) ;
            memmove(&(matrix->colptr[fp]),&(matrix->colptr[fp+1]), sizeof(mess_int_t)*elementstocopy);
            matrix->nnz--;
            mess_try_realloc(matrix->values, double* , sizeof(double) * matrix->nnz);
            mess_try_realloc(matrix->colptr, mess_int_t* , sizeof(mess_int_t) * matrix->nnz);

            for ( i = row +1; i <= matrix->rows; i++) matrix->rowptr[i]--;
        }
    }  else if ( MESS_IS_CSC(matrix)) {
        if ( value != 0.0 ) {
            mess_int_t i;
            int found = 0;
            for ( i= matrix->colptr[col] ; i < matrix->colptr[col+1]; i++){
                if ( matrix->rowptr[i] == row ) {
                    matrix->values[i] = value;
                    found = 1;
                }
            }
            if (!found) {
                mess_int_t elementstocopy;
                struct pair_mess_int_t_double *pld=NULL;
                mess_int_t j;
                mess_try_realloc(matrix->values, double *, (matrix->nnz+1) * sizeof(double));
                mess_try_realloc(matrix->rowptr, mess_int_t *, (matrix->nnz+1) * sizeof(mess_int_t));
                elementstocopy = matrix->nnz-matrix->colptr[col];
                if ( matrix->nnz > 0 ) {
                    memmove(&(matrix->values[matrix->colptr[col]+1]), &(matrix->values[matrix->colptr[col]]), elementstocopy * sizeof(double));
                    memmove(&(matrix->rowptr[matrix->colptr[col]+1]), &(matrix->rowptr[matrix->colptr[col]]), elementstocopy * sizeof(mess_int_t));
                }
                for (i = col+1; i <= matrix->cols; i++) matrix->colptr[i]++;
                matrix->values[matrix->colptr[col]]=value;
                matrix->rowptr[matrix->colptr[col]]=row;
                matrix->nnz++;

                mess_try_alloc(pld, struct pair_mess_int_t_double *, sizeof(struct pair_mess_int_t_double)* matrix->rows);
                pos = 0;
                for (j = matrix->colptr[col]; j < matrix->colptr[col+1]; j++){
                    pld[pos].l = matrix->rowptr[j];
                    pld[pos].d = matrix->values[j];
                    pos++;
                }
                qsort(pld, pos, sizeof(struct pair_mess_int_t_double), __cmp_mess_int_t_double);
                pos = 0;
                for (j = matrix->colptr[col]; j < matrix->colptr[col+1]; j++){
                    matrix->rowptr[j] = pld[pos].l;
                    matrix->values[j] = pld[pos].d;
                    pos++;
                }
                mess_free(pld);
            }
        } else {
            if ( matrix->nnz < 1 ) return(0);
            // remove the element
            mess_int_t i, fp=-1;
            mess_int_t elementstocopy;
            for (i=matrix->colptr[col]; i<matrix->colptr[col+1]; i++) {
                if ( matrix->rowptr[i]==row) {fp = i; break; }
            }
            if ( fp < 0 ) return(0);
            elementstocopy = matrix->nnz - fp - 1;
            memmove(&(matrix->values[fp]),&(matrix->values[fp+1]), sizeof(double)*elementstocopy) ;
            memmove(&(matrix->rowptr[fp]),&(matrix->rowptr[fp+1]), sizeof(mess_int_t)*elementstocopy);
            matrix->nnz--;
            mess_try_realloc(matrix->values, double* , sizeof(double) * matrix->nnz);
            mess_try_realloc(matrix->rowptr, mess_int_t* , sizeof(mess_int_t) * matrix->nnz);
            for ( i = col +1; i <= matrix->cols; i++) matrix->colptr[i]--;
        }
    } else if ( MESS_IS_COORD( matrix  )) {
        if ( value !=0.0 ) {
            int found = 0;
            mess_int_t i ;
            for ( i = 0 ; i< matrix->nnz; i++) {
                if ( matrix->rowptr[i] == row && matrix->colptr[i] == col) { matrix->values[i] = value; found = 1;  break; }
            }
            if ( !found ) {
                matrix->nnz++;
                mess_try_realloc(matrix->rowptr, mess_int_t *, sizeof(mess_int_t) * matrix->nnz);
                mess_try_realloc(matrix->colptr, mess_int_t *, sizeof(mess_int_t) * matrix->nnz);
                mess_try_realloc(matrix->values, double *, sizeof(double) * matrix->nnz);
                matrix->rowptr[matrix->nnz-1] = row;
                matrix->colptr[matrix->nnz-1] = col;
                matrix->values[matrix->nnz-1] = value;
            }
        } else {
            if ( matrix->nnz < 1) return (0);
            mess_int_t i , fp = -1;
            mess_int_t elementstocopy;
            for ( i = 0 ; i< matrix->nnz; i++) {
                if ( matrix->rowptr[i] == row && matrix->colptr[i] == col) { fp = i ;  break; }
            }
            if ( fp < 0 ) return(0);
            elementstocopy = matrix->nnz - fp - 1;
            memmove(&(matrix->rowptr[fp]), &(matrix->rowptr[fp+1]), elementstocopy*sizeof(mess_int_t));
            memmove(&(matrix->colptr[fp]), &(matrix->colptr[fp+1]), elementstocopy*sizeof(mess_int_t));
            memmove(&(matrix->values[fp]), &(matrix->values[fp+1]), elementstocopy*sizeof(double));
            matrix->nnz--;
        }
    }

    return(0);
}       /* -----  end of function mess_matrix_setelement  ----- */

