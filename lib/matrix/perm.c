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
 * @file lib/matrix/perm.c
 * @brief Different matrix permutation operations.
 * @author @koehlerm
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "shellsort.c"


/**
 * @internal
 *
 * @brief Local function to invert a permutation.
 * @param[in] p      input integer array to be permuted
 * @param[in] n      input number of elements in the permutation
 * @return the permuted array, allocated by malloc
 *
 * The @ref perm_inv function inverts a permutation.
 *
 * @attention Interal use only.
 */
static mess_int_t *perm_inv (const mess_int_t *p,mess_int_t n)
{
    MSG_FNAME(__func__);
    mess_int_t k, *pinv ;
    if (!p) return (NULL);
    mess_try_alloc2( pinv, mess_int_t*, n*sizeof (mess_int_t)) ;
    if (!pinv) return (NULL) ;
    for (k = 0 ; k < n ; k++) {
        pinv [p [k]] = k ;
    }
    return (pinv) ;
}


/**
 * @internal
 *
 * @brief Permute a Compressed Sparse Row matrix.
 * @param[in,out] matrix matrix to permute
 * @param[in] p      input row permutation
 * @param[in] qinv   input inverse column permutation
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __mess_matrix_perm_csr function permutes a Compressed
 * Sparse Row matrix. The algorithm works not in-place.
 *
 * @attention Internal use only.
 */
static int  __mess_matrix_perm_csr ( mess_matrix matrix, const mess_int_t *p, const mess_int_t *qinv )
{
    MSG_FNAME(__func__);
    mess_int_t *colptr=NULL;
    mess_int_t *rowptr=NULL;
    double *values=NULL;
    mess_double_cpx_t *values_cpx=NULL;
    mess_int_t i= 0, j = 0;
    mess_int_t pos = 0;
    mess_int_t row = 0, col = 0;
    int ret = 0;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_csr(matrix);

    if (MESS_IS_REAL(matrix)) {
        mess_try_alloc(values, double *, sizeof(double) * matrix->nnz);
    } else if ( MESS_IS_COMPLEX( matrix)) {
        mess_try_alloc(values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * matrix->nnz);
    } else{
         MSG_ERROR("unknown datatype\n");
         return MESS_ERROR_DATATYPE;
    }
    mess_try_alloc(colptr, mess_int_t*, sizeof (mess_int_t) * matrix->nnz);
    mess_try_alloc(rowptr, mess_int_t *,sizeof (mess_int_t) * (matrix->rows+1));
    rowptr[0] = 0;

    /*-----------------------------------------------------------------------------
     *  permute
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_REAL(matrix)) {
        for ( i = 0; i<matrix->rows; i++) {
            row = ( p == NULL) ? i : p[i];
            for (j = matrix->rowptr[row]; j<matrix->rowptr[row+1]; j++){
                col =  (qinv == NULL) ? matrix->colptr[j] : qinv[matrix->colptr[j]];
                colptr[pos] = col;
                values[pos] = matrix->values[j];
                pos++;
            }
            rowptr[i+1] = pos;
        }

    } else {
        for ( i = 0; i<matrix->rows; i++) {
            row = ( p == NULL) ? i : p[i];
            for (j = matrix->rowptr[row]; j<matrix->rowptr[row+1]; j++){
                col =  (qinv == NULL) ? matrix->colptr[j] : qinv[matrix->colptr[j]];
                colptr[pos] = col;
                values_cpx[pos] = matrix->values_cpx[j];
                pos++;
            }
            rowptr[i+1] = pos;
        }
    }

    mess_free(matrix->colptr);
    mess_free(matrix->rowptr);
    if ( matrix->values != NULL) mess_free(matrix->values);
    if ( matrix->values_cpx != NULL) mess_free(matrix->values_cpx);

    matrix->colptr = colptr;
    matrix->rowptr = rowptr;
    matrix->values = values;
    matrix->values_cpx = values_cpx;

    ret  = mess_matrix_sort(matrix);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_sort);


    return 0;
}       /* -----  end of function __mess_matrix_perm_csr  ----- */

/**
 * @internal
 *
 * @brief Permute a Compressed Sparse Column matrix (local function).
 * @param[in,out] matrix matrix to permute
 * @param[in] pinv   input row permutation
 * @param[in] q   input inverse column permutation
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __mess_matrix_perm_csc function permutes a Compressed
 * Sparse Column matrix. The algorithm works not in-place.
 *
 * @attention Internal use only.
 */
static int  __mess_matrix_perm_csc ( mess_matrix matrix, const mess_int_t *pinv, const mess_int_t *q )
{
    MSG_FNAME(__func__);
    mess_int_t *colptr=NULL;
    mess_int_t *rowptr=NULL;
    double *values=NULL;
    mess_double_cpx_t *values_cpx=NULL;
    mess_int_t i= 0, j = 0;
    mess_int_t pos = 0;
    mess_int_t row = 0, col = 0;
    int ret = 0;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_csc(matrix);

    if (MESS_IS_REAL(matrix)) {
        mess_try_alloc(values, double *, sizeof(double) * matrix->nnz);
    } else if ( MESS_IS_COMPLEX( matrix)) {
        mess_try_alloc(values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * matrix->nnz);
    } else{
         MSG_ERROR("unknown datatype\n");
         return MESS_ERROR_DATATYPE;
    }


    mess_try_alloc(rowptr, mess_int_t*, sizeof (mess_int_t) * matrix->nnz);
    mess_try_alloc(colptr, mess_int_t *,sizeof (mess_int_t) * (matrix->cols+1));
    colptr[0] = 0;

    /*-----------------------------------------------------------------------------
     *  permute
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_REAL(matrix)) {
        for ( i = 0; i<matrix->cols; i++) {
            col = ( q == NULL) ? i : q[i];
            for (j = matrix->colptr[col]; j<matrix->colptr[col+1]; j++){
                row =  (pinv == NULL) ? matrix->rowptr[j] : pinv[matrix->rowptr[j]];
                rowptr[pos] = row;
                values[pos] = matrix->values[j];
                pos++;
            }
            colptr[i+1] = pos;
        }

    } else {
        for ( i = 0; i<matrix->cols; i++) {
            col = ( q == NULL) ? i : q[i];
            for (j = matrix->colptr[col]; j<matrix->colptr[col+1]; j++){
                row =  (pinv == NULL) ? matrix->rowptr[j] : pinv[matrix->rowptr[j]];
                rowptr[pos] = row;
                values_cpx[pos] = matrix->values_cpx[j];
                pos++;
            }
            colptr[i+1] = pos;
        }

    }

    mess_free(matrix->colptr);
    mess_free(matrix->rowptr);
    if ( matrix->values != NULL)mess_free(matrix->values);
    if ( matrix->values_cpx != NULL)mess_free(matrix->values_cpx);

    matrix->colptr = colptr;
    matrix->rowptr = rowptr;
    matrix->values = values;
    matrix->values_cpx = values_cpx;

    ret  = mess_matrix_sort(matrix);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_sort);


    return 0;
}       /* -----  end of function __mess_matrix_perm_csr  ----- */


/**
 * @internal
 * @brief Permute a dense matrix.
 * @param[in,out]  matrix matrix to permute
 * @param[in] p  input row permutation
 * @param[in] qinv   input column permutation
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __mess_matrix_perm_dense function permutes a dense matrix. The algorithm
 * does not work in-place but using @openmp
 *
 * @attention Internal use only.
 */
static int  __mess_matrix_perm_dense ( mess_matrix matrix, mess_int_t *p, mess_int_t *qinv )
{
    MSG_FNAME(__func__);
    double *values = NULL;
    mess_double_cpx_t * values_cpx = NULL;
    mess_int_t i, j;
    mess_int_t col, row;

    mess_check_nullpointer(matrix);
    mess_check_real_or_complex(matrix);

    if ( MESS_IS_REAL(matrix)){
            mess_try_alloc(values, double *, sizeof(double)*matrix->ld*matrix->cols);
#ifdef _OPENMP
#pragma omp parallel for private(j,i,col,row) default(shared)
#endif
            for ( j = 0; j < matrix->cols; j++){
                col = (qinv==NULL) ? j : qinv[j];
                for (i=0; i < matrix->rows; i++){
                    row = (p==NULL) ? i : p[i];
                    values[col*matrix->ld+row] = matrix->values[j*matrix->ld+i];
                }
            }
            mess_free(matrix->values);
            matrix->values = values;
    } else if ( MESS_IS_COMPLEX(matrix)){
            mess_try_alloc(values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*matrix->ld*matrix->cols);
#ifdef _OPENMP
#pragma omp parallel for private(j,i,col,row) default(shared)
#endif
            for ( j = 0; j < matrix->cols; j++){
                col = (qinv==NULL) ? j : qinv[j];
                for (i=0; i < matrix->rows; i++){
                    row = (p==NULL) ? i : p[i];
                    values_cpx[col*matrix->ld+row] = matrix->values_cpx[j*matrix->ld+i];
                }
            }
            mess_free(matrix->values_cpx);
            matrix->values_cpx = values_cpx;
        }

    return 0;
}       /* -----  end of function __mess_perm_dense  ----- */



/**
 * @brief Permute a general matrix.
 * @param[in,out] matrix matrix to permute
 * @param[in] p          input row permutation of length \f$ matrix \to rows \f$
 * @param[in] q              input column permutation of length \f$ matrix \to cols\f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_perm function permutes a matrix in rows and column space. \n
 * If \f$ p \f$ or \f$ q \f$  is @c NULL, it is assumed that it would be the identity permutation.\n
 * The underlying algorithms do not work in-place at the moment.
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_permcopy
 */
int mess_matrix_perm(mess_matrix matrix, mess_int_t *p, mess_int_t *q) {
    MSG_FNAME(__func__);
    mess_int_t *pinv = NULL;
    mess_int_t *qinv = NULL;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    if ( p == NULL && q == NULL) return 0;

    /*-----------------------------------------------------------------------------
     *  perpare
     *-----------------------------------------------------------------------------*/
    pinv = perm_inv(p, matrix->rows);
    qinv = perm_inv(q, matrix->cols);

    if ( MESS_IS_CSR(matrix)){
        ret = __mess_matrix_perm_csr(matrix, p, qinv);
        mess_free(qinv);
        mess_free(pinv);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_perm_csr);
    } else if (MESS_IS_CSC(matrix)) {
        ret = __mess_matrix_perm_csc(matrix, pinv, q);
        mess_free(pinv);
        mess_free(qinv);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_perm_csc);
    } else if (MESS_IS_DENSE(matrix)){
        ret = __mess_matrix_perm_dense(matrix, pinv, qinv);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_perm_dense);
        mess_free(pinv);
        mess_free(qinv);
    } else {
        mess_free(qinv);
        mess_free(pinv);
        MSG_ERROR("Storage type: %s  unsupported/unknown.\n", mess_storage_t_str(matrix->store_type));
        return MESS_ERROR_STORAGETYPE;
    }
    return 0;
}




/**
 * @brief Copy and permute a matrix.
 * @param[in]  matrix   input matrix
 * @param[in]  p     input row permutation of length \f$ matrix \to rows \f$
 * @param[in]  q     input column permutation of length \f$ matrix \to cols \f$
 * @param[out] out  permuted matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_permcopy function copies and permutes a matrix.  \n
 * It is similar to \ref mess_matrix_perm. \n
 * If a permutation vector is set @c NULL, the identity is assumed.
 *
 * \attention
 * This subroutine currently only works for Compressed Sparse Rows storage.
 *
 * \sa  mess_matrix_perm
 *
 */
int mess_matrix_permcopy(mess_matrix matrix, mess_int_t *p, mess_int_t *q, mess_matrix out){
    MSG_FNAME(__func__);
    double *values = NULL;
    mess_double_cpx_t *values_cpx = NULL;
    mess_int_t *colptr = NULL;
    mess_int_t *rowptr = NULL;
    mess_int_t *qinv;
    mess_int_t i, j ;
    mess_int_t pos = 0;
    mess_int_t row, col;
    int ret = 0 ;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(out);
    mess_check_csr(matrix);

    if ( p == NULL && q == NULL) {
        ret =  mess_matrix_copy(matrix, out);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        return 0;
    }

    MESS_MATRIX_RESET(out);
    qinv = perm_inv(q, matrix->cols);

    if (MESS_IS_REAL(matrix)) {
        mess_try_alloc(values, double *, sizeof(double) * matrix->nnz);
    } else if ( MESS_IS_COMPLEX( matrix)) {
        mess_try_alloc(values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * matrix->nnz);
    }
    mess_try_alloc(colptr, mess_int_t *, sizeof(mess_int_t) * matrix->nnz);
    mess_try_alloc(rowptr, mess_int_t *, sizeof(mess_int_t) * (matrix->rows + 1));

    rowptr[0] = 0;
    for ( i = 0; i<matrix->rows; i++) {
        row = ( p == NULL) ? i : p[i];
        for (j = matrix->rowptr[row]; j<matrix->rowptr[row+1]; j++){
            col =  (qinv == NULL) ? matrix->colptr[j] : qinv[matrix->colptr[j]];
            colptr[pos] = col;
            if ( MESS_IS_REAL(matrix)) {
                values[pos] = matrix->values[j];
            }
            if ( MESS_IS_COMPLEX(matrix)){
                values_cpx[pos] = matrix->values_cpx[j];
            }
            pos++;
        }
        rowptr[i+1] = pos;
        __shellsort(colptr, values, values_cpx, rowptr[i], rowptr[i+1]-1, matrix->data_type);
    }

    out->colptr = colptr;
    out->rowptr = rowptr;
    out->values = values;
    out->values_cpx = values_cpx;
    out->cols = matrix->cols;
    out->data_type = matrix->data_type;
    out->nnz = matrix->nnz;
    out->rows = matrix->rows;
    out->store_type = MESS_CSR;
    out->symmetry = MESS_GENERAL;

    mess_free(qinv);
    return 0;

}


/**
 * @internal
 *
 * @brief Exchange two vectors/columns (real valued).
 * @param[in,out] col1 pointer to the first  column
 * @param[in,out] col2 pointer to the second column
 * @param[in] N  input length of col1 and col2
 * @return Nothing.
 *
 * The @ref ex_col_real function exchanges two real valued columns or arrays of
 * length N.
 *
 * @attention Internal use only.
 */
static void ex_col_real(double *col1, double *col2, mess_int_t N){
    mess_int_t i;
    double t;
    for ( i = 0; i< N; i++){
        t = col1[i];
        col1[i]=col2[i];
        col2[i]=t;
    }
}

/**
 * @internal
 *
 * @brief Exchange two vectors/columns (complex valued).
 * @param[in,out] col1 pointer to the first column
 * @param[in,out] col2 pointer to the second column
 * @param[in] N  input length of col1 and col2
 * @return Nothing.
 *
 * The @ref ex_col_real function exchanges two complex valued columns or arrays
 * of length N.
 *
 * @attention Internal use only.
 */
static void ex_col_cpx(mess_double_cpx_t *col1, mess_double_cpx_t *col2, mess_int_t N){
    mess_int_t i;
    mess_double_cpx_t t;
    for ( i = 0; i< N; i++){
        t = col1[i];
        col1[i]=col2[i];
        col2[i]=t;
    }
}


/**
 * @brief Permute the columns of a matrix inplace.
 * @param[in,out] A matrix to permute
 * @param[in] perm  input column permutation
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_colperm function permutes the columns of a matrix \f$ A \f$ by the given
 * permutation \f$ perm \f$, where \f$ perm[i]=j \f$ means the \f$ i-th \f$ column of \f$ A \f$ will be moved to the
 * \f$ j \f$-th column in the permuted \f$ A \f$. \n
 * The algorithm works inplace. \n
 * If OpenMP is available the permuation process will be done in parallel if the permutation allows it.
 *
 * @sa mess_matrix_colpermcopy
 */
int mess_matrix_colperm( mess_matrix A, mess_int_t *perm)
{
    MSG_FNAME(__func__);
    mess_int_t jmin, start, next;
    mess_int_t ok = 0;
    double *tmp = NULL;
    mess_double_cpx_t *tmp_cpx = NULL ;
    mess_int_t *iperm;
    mess_int_t *iperm2;
    mess_int_t *scycle;
    mess_int_t Lscy;
    mess_int_t i,j;
    mess_int_t sc = 0;
    mess_int_t oldnext;
    mess_int_t M, N;
    int cpx = 0;

    mess_check_nullpointer(A);
    mess_check_nullpointer(perm);
    mess_check_dense(A);
    mess_check_real_or_complex(A);

    if ( MESS_IS_COMPLEX(A)) cpx = 1;
    M=A->rows;
    N=A->cols;

    iperm = perm_inv(perm, A->cols);
    mess_try_alloc(iperm2, mess_int_t *, sizeof(mess_int_t)*A->cols);
    memcpy(iperm2, iperm, sizeof(mess_int_t)*A->cols);
    Lscy = A->cols/2 +2 ;
    mess_try_alloc(scycle, mess_int_t *, sizeof(mess_int_t)*Lscy);

    for ( i =0;i<Lscy ;i++ ) { scycle[i] = -1;}

    /*-----------------------------------------------------------------------------
     *  search the start of the cycles
     *-----------------------------------------------------------------------------*/
    while ( !ok ) {
        jmin = 0;
        while( jmin < N && iperm2[jmin]==jmin ) jmin++;
        if ( jmin == N) { ok = 1; continue; }
        start = jmin;
        next = jmin;
        scycle[sc++]=jmin;
        do {
            oldnext = next;
            next = iperm2[next];
            iperm2[oldnext] = oldnext;
        } while ( start != next);
    }
    ok = 0;

    /*-----------------------------------------------------------------------------
     *  permute the matrix
     *-----------------------------------------------------------------------------*/
    // if ( cpx ) {mess_try_alloc(tmp_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*M);}
    // else {mess_try_alloc(tmp, double*, sizeof(double)*M);}
#ifdef _OPENMP
#pragma omp parallel for private(j, i, tmp, tmp_cpx, jmin, start, next, oldnext) default (shared)
#endif
    for ( j = 0; j < sc; j++){
        tmp=NULL;tmp_cpx=NULL;
        if ( cpx ) {
            mess_try_alloc2( tmp_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*M);
        } else {
            mess_try_alloc2( tmp,  double*, sizeof(double)*M);
        }
        jmin = scycle[j];
        if (cpx ){
            for ( i = 0; i < M ; i++ ) { tmp_cpx [ i ] = A->values_cpx[perm[jmin]*A->ld+i];}
        } else {
            for ( i = 0; i < M ; i++ ) { tmp [ i ] = A->values[perm[jmin]*A->ld+i];}
        }
        start = jmin;
        next = jmin;
        do {
            if ( MESS_IS_REAL(A))   ex_col_real(&(A->values[next*A->ld]),tmp, M);
            else ex_col_cpx(&(A->values_cpx[next*A->ld]),tmp_cpx, M);
            oldnext = next;
            next = iperm[next];
            iperm[oldnext] = oldnext;
        } while ( start != next);
        if ( cpx  )mess_free(tmp_cpx); else
            mess_free(tmp);
    }
    mess_free(iperm);
    mess_free(iperm2);
    mess_free(scycle);
    return 0;
}

/**
 * @brief Permute the columns of a matrix while copying.
 * @param[in] A  input matrix to permute
 * @param[in] perm  input column permutation
 * @param[out] B permuted matrix
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_colpermcopy function permutes the columns of a matrix \f$ A \f$  by a given
 * permutation \f$ perm \f$ and stores the permuted matrix in \f$ B \f$,
 * where \f$ perm[i]=j  \f$ means the \f$ i \f$-th column of \f$ A \f$ will be moved to the
 * \f$ j \f$-th column in the permuted matrix \f$ B \f$.
 *
 * @sa mess_matrix_colperm
 *
 */
int mess_matrix_colpermcopy( mess_matrix A, mess_int_t *perm, mess_matrix B)
{
    MSG_FNAME(__func__);
    mess_int_t i,j;
    mess_int_t p;
    int ret;

    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_dense(A);
    mess_check_real_or_complex(A);


    MESS_MATRIX_RESET(B);
    ret=mess_matrix_alloc(B, A->rows, A->cols, A->nnz, MESS_DENSE, A->data_type);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    if ( MESS_IS_REAL(A)) {
#ifdef _OPENMP
#pragma omp parallel for private (j,i,p) default (shared)
#endif
        for(j = 0; j < A->cols; j++){
            p = (perm == NULL) ? j : perm[j];
            for (i = 0; i < A->rows; i++){
                B->values[j*B->ld+i] = A->values[p*A->ld+i];
            }
        }
    }else {
#ifdef _OPENMP
#pragma omp parallel for private (j,i,p) default (shared)
#endif
        for(j = 0; j < A->cols; j++){
            p = (perm == NULL) ? j : perm [j];
            for (i = 0; i < A->rows; i++){
                B->values_cpx[j*B->ld+i] = A->values_cpx[p*A->ld+i];
            }
        }
    }
    return 0;
}

