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
 * @file lib/formats/cholmod.c
 * @brief Import and export vectors and dense/csc matrices from @mess to @cholmod.
 * @author @mbehr
 */


#include "mess/mess.h"
#include "mess/interface_cholmod.h"
#include "mess/config.h"
#include "mess/error_macro.h"

#ifdef MESS_HAVE_CHOLMOD
#ifdef MESS_USE_SUITESPARSE3
#define LONG long
#else
#define LONG SuiteSparse_long
#endif

/**
 * @brief Convert mess_vector to cholmod_dense.
 * @param[in] v     input vector
 * @param[out] v_cholmod cholmod_dense vector \f$ v \f$
 * @param[in] c      input cholmod_common structure
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_convert_dense_to_cholmod converts a @ref mess_vector data structure to @c cholmod_dense data structure.
 */
int mess_vector_convert_dense_to_cholmod(mess_vector v, cholmod_dense **v_cholmod, cholmod_common *c ){
    MSG_FNAME(__func__);

    mess_int_t i;
    mess_check_nullpointer(c);
    mess_check_nullpointer(v);
    mess_check_real_or_complex(v);

    if(c->itype!=CHOLMOD_LONG || c->dtype!=CHOLMOD_DOUBLE){
        MSG_ERROR("Argument error in cholmod_common structure. c.itype has to be CHOLMOD_LONG and c.dtype has to be CHOLMOD_DOUBLE\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if(MESS_IS_REAL(v)){
        *v_cholmod = cholmod_l_allocate_dense(v->dim,1,v->dim,CHOLMOD_REAL,c);
        mess_check_nullpointer(*v_cholmod);
        double *x = (*v_cholmod)->x;
        for ( i = 0; i < v->dim; i++) {
            x[i]=v->values[i];
        }
    }
    else{
        *v_cholmod = cholmod_l_allocate_dense(v->dim,1,v->dim,CHOLMOD_COMPLEX,c);
        mess_check_nullpointer(*v_cholmod);
        mess_double_cpx_t *x = (*v_cholmod)->x;
        for ( i = 0; i < v->dim; i++) {
            x[i]=v->values_cpx[i];
        }

    }
    return 0;
}

/**
 * @brief Convert cholmod_dense structure to mess_vector.
 * @param[in] v_cholmod input cholmod_dense structure
 * @param[out] v        converted vector
 * @param[in] c      input cholmod_common structure
 * @return zero on success or a non-zero error value otherwise
 *
 * The mess_vector_convert_cholmod_to_dense converts a cholmod_dense data structure to a @ref mess_vector data structure. \n
 * If cholmod_dense structure has not exactly one column an error will be returned.
 */
int mess_vector_convert_cholmod_to_dense(cholmod_dense *v_cholmod, mess_vector v, cholmod_common *c){
    MSG_FNAME(__func__);

    int ret=0;
    mess_int_t i;
    mess_check_nullpointer(c);
    mess_check_nullpointer(v_cholmod);
    mess_check_nullpointer(v);
    mess_check_real_or_complex(v);

    if(v_cholmod->ncol!=1){
        MSG_ERROR("Cholmod_dense structure has not exactly one columns\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if(c->itype!=CHOLMOD_LONG || c->dtype!=CHOLMOD_DOUBLE){
        MSG_ERROR("Argument error in cholmod_common structure. c.itype has to be CHOLMOD_LONG and c.dtype has to be CHOLMOD_DOUBLE\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if(v_cholmod->xtype==CHOLMOD_REAL){
        ret = mess_vector_toreal(v);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal);
        ret = mess_vector_resize(v,v_cholmod->nrow);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);
        double *x = v_cholmod->x;
        for ( i = 0; i < v_cholmod->nrow; i++) {
            v->values[i]=x[i];
        }
    }else if (v_cholmod->xtype==CHOLMOD_COMPLEX){
        ret = mess_vector_tocomplex(v);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
        ret = mess_vector_resize(v,v_cholmod->nrow);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);
        mess_double_cpx_t *x = v_cholmod->x;
        for ( i = 0; i < v_cholmod->nrow; i++) {
            v->values_cpx[i]=x[i];
        }
    }else{
        MSG_ERROR("unsupported cholmod datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return ret;
}


/**
 * @brief Convert dense mess_matrix to cholmod_dense.
 * @param[in] A     input dense matrix
 * @param[out] A_cholmod    cholmod_dense matrix \f$ A \f$
 * @param[in] c      input cholmod_common structure
 * @return zero on success or a non-zero error value otherwise
 *
 * The mess_matrix_convert_dense_to_cholmod converts a dense @ref mess_matrix data structure to a cholmod_dense data
 * structure.
 *
 */
int mess_matrix_convert_dense_to_cholmod(mess_matrix A, cholmod_dense **A_cholmod, cholmod_common *c ){
    MSG_FNAME(__func__);

    mess_int_t i,j;
    mess_check_nullpointer(c);
    mess_check_nullpointer(A);
    mess_check_real_or_complex(A);
    mess_check_dense(A);

    if(c->itype!=CHOLMOD_LONG || c->dtype!=CHOLMOD_DOUBLE){
        MSG_ERROR("Argument error in cholmod_common structure. c.itype has to be CHOLMOD_LONG and c.dtype has to be CHOLMOD_DOUBLE\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if(MESS_IS_REAL(A)){
        *A_cholmod = cholmod_l_allocate_dense(A->rows,A->cols,A->ld,CHOLMOD_REAL,c);
        mess_check_nullpointer(*A_cholmod);
        double *x = (*A_cholmod)->x;
        for ( i = 0; i < A->cols; i++) {
            for (j=0; j < A->rows; j++){
                x[i*(*A_cholmod)->d+j]=A->values[i*A->ld+j];
            }
        }
    }else{
        *A_cholmod = cholmod_l_allocate_dense(A->rows,A->cols,A->ld,CHOLMOD_COMPLEX,c);
        mess_check_nullpointer(*A_cholmod);
        mess_double_cpx_t * x = (*A_cholmod)->x;
        for ( i = 0; i < A->cols; i++) {
            for (j=0; j < A->rows; j++){
                x[i*(*A_cholmod)->d+j]=A->values_cpx[i*A->ld+j];
            }
        }
    }
    return 0;
}


/**
 * @brief Convert cholmod_dense to dense mess_matrix.
 * @param[in] A_cholmod     cholmod_dense input matrix
 * @param[out] A         dense matrix
 * @param[in] c      input cholmod_common structure
 *
 * @return zero on success or a non-zero error value otherwise
 *
 * The mess_matrix_convert_dense_to_cholmod converts a cholmod_dense data structure to a @ref mess_matrix data structure.
 *
 */
int mess_matrix_convert_cholmod_to_dense(cholmod_dense *A_cholmod, mess_matrix A, cholmod_common *c){
    MSG_FNAME(__func__);

    int ret=0;
    mess_int_t i,j;
    mess_check_nullpointer(c);
    mess_check_nullpointer(A_cholmod);
    mess_check_nullpointer(A);

    if(c->itype!=CHOLMOD_LONG || c->dtype!=CHOLMOD_DOUBLE){
        MSG_ERROR("Argument error in cholmod_common structure. c.itype has to be CHOLMOD_LONG and c.dtype has to be CHOLMOD_DOUBLE\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if(A_cholmod->xtype==CHOLMOD_REAL){
        ret=mess_matrix_alloc(A, A_cholmod->nrow, A_cholmod->ncol, A_cholmod->nrow*A_cholmod->ncol,MESS_DENSE,MESS_REAL); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
        double *x = A_cholmod->x;
        for ( i = 0; i < A->cols; i++) {
            for (j=0; j < A->rows; j++){
                A->values[i*A->ld+j]=x[i*A_cholmod->d+j];
            }
        }
    }else if (A_cholmod->xtype==CHOLMOD_COMPLEX){
        ret=mess_matrix_alloc(A, A_cholmod->nrow, A_cholmod->ncol, A_cholmod->nrow*A_cholmod->ncol,MESS_DENSE,MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
        mess_double_cpx_t * x = A_cholmod->x;
        for ( i = 0; i < A->cols; i++) {
            for (j=0; j < A->rows; j++){
                A->values_cpx[i*A->ld+j]=x[i*A_cholmod->d+j];
            }
        }
    }else{
        MSG_ERROR("unsupported cholmod datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return ret;
}

/**
 * @brief Convert CSC mess_matrix to cholmod_sparse.
 * @param[in] A     input CSC matrix
 * @param[out] A_cholmod    cholmod_sparse matrix \f$ A \f$
 * @param[in] c      input cholmod_common structure
 * @return zero on success or a non-zero error value otherwise
 *
 * The mess_matrix_convert_csc_to_cholmod converts a CSC @ref mess_matrix to cholmod_sparse data structure.
 *
 *
 */
int mess_matrix_convert_csc_to_cholmod(mess_matrix A, cholmod_sparse **A_cholmod, cholmod_common *c){
    MSG_FNAME(__func__);

    mess_int_t i,j;
    mess_check_nullpointer(c);
    mess_check_nullpointer(A);
    mess_check_real_or_complex(A);
    mess_check_csc(A);

    if(c->itype!=CHOLMOD_LONG || c->dtype!=CHOLMOD_DOUBLE){
        MSG_ERROR("Argument error in cholmod_common structure. c.itype has to be CHOLMOD_LONG and c.dtype has to be CHOLMOD_DOUBLE\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if(MESS_IS_REAL(A)){
        //0 not  sorted, 1 packed, 0 unsymmetric
        *A_cholmod = cholmod_l_allocate_sparse(A->rows,A->cols,A->nnz,0,1,0,CHOLMOD_REAL,c);
        mess_check_nullpointer(*A_cholmod);
        //colptr
        LONG *colptr = (*A_cholmod)->p;
        for (i = 0 ; i <= A->cols; i++){
            colptr[i]=A->colptr[i];
        }
        //rowptr
        LONG    *rowptr = (*A_cholmod)->i;
        //dataptr
        double  *x = (*A_cholmod)->x;
        for (i = 0 ; i < A->cols; i++){
            for (j = A->colptr[i]; j<A->colptr[i+1]; j++){
                rowptr[j]=A->rowptr[j];
                x[j]=A->values[j];
            }
        }
    }else{
        //0 not sorted, 1 packed, 0 unsymmetric
        (*A_cholmod) = cholmod_l_allocate_sparse(A->rows,A->cols,A->nnz,0,1,0,CHOLMOD_COMPLEX,c);
        mess_check_nullpointer(*A_cholmod);
        //colptr
        LONG *colptr = (*A_cholmod)->p;
        for (i = 0 ; i <= A->cols; i++){
            colptr[i]=A->colptr[i];
        }
        //rowptr
        LONG *rowptr  = (*A_cholmod)->i;
        //dataptr
        mess_double_cpx_t *x = (*A_cholmod)->x;
        for (i = 0 ; i < A->cols; i++){
            for (j = A->colptr[i]; j<A->colptr[i+1]; j++){
                rowptr[j]=A->rowptr[j];
                x[j]=A->values_cpx[j];
            }
        }
    }
    return 0;
}

/**
 * @brief Convert cholmod_sparse to CSC mess_matrix.
 * @param[in] A_cholmod input cholmod_sparse matrix
 * @param[out] A        converted CSC matrix
 * @param[in] c      input cholmod_common structure
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_convert_cholmod_to_csc converts a cholmod_sparse data structure to a CSC @ref mess_matrix data structure.
 */
int mess_matrix_convert_cholmod_to_csc(cholmod_sparse *A_cholmod, mess_matrix A, cholmod_common *c){

    MSG_FNAME(__func__);

    int ret;
    mess_int_t i,j;
    mess_check_nullpointer(c);
    mess_check_nullpointer(A_cholmod);

    if(c->itype!=CHOLMOD_LONG || c->dtype!=CHOLMOD_DOUBLE){
        MSG_ERROR("Argument error in cholmod_common structure. c.itype has to be CHOLMOD_LONG and c.dtype has to be CHOLMOD_DOUBLE\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if(A_cholmod->xtype==CHOLMOD_REAL){
        MESS_MATRIX_RESET(A);
        ret=mess_matrix_alloc(A, A_cholmod->nrow, A_cholmod->ncol, A_cholmod->nzmax,MESS_CSC,MESS_REAL);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
        //mess_try_alloc(A->colptr, mess_int_t*, sizeof(mess_int_t)*(A_cholmod->ncol+1));
        //mess_try_alloc(A->rowptr, mess_int_t*, sizeof(mess_int_t)*(A_cholmod->nzmax));
        //mess_try_alloc(A->values, double*,       sizeof(double)*(A_cholmod->nzmax));

        //column pointer
        LONG *colptr = A_cholmod->p;
        for (i = 0 ; i <= A_cholmod->ncol; i++){
            A->colptr[i]=colptr[i];
        }

        //rowptr
        LONG *rowptr = A_cholmod->i;

        //dataptr
        double  *x  = A_cholmod->x;
        if(A_cholmod->packed){
            for (i = 0; i < A_cholmod->ncol; i++){
                for (j = colptr[i]; j<colptr[i+1]; j++){
                    A->rowptr[j]=rowptr[j];
                    A->values[j]=x[j];
                }
            }
        }else{
            LONG *nz = A_cholmod->nz;
            for (i = 0; i < A_cholmod->ncol; i++){
                for (j = colptr[i]; j<nz[i]; j++){
                    A->rowptr[j]=rowptr[j];
                    A->values[j]=x[j];
                }
            }
        }
    }else if(A_cholmod->xtype==CHOLMOD_COMPLEX){
        ret=mess_matrix_alloc(A, A_cholmod->nrow, A_cholmod->ncol, A_cholmod->nzmax,MESS_CSC,MESS_COMPLEX);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
        //mess_try_alloc(A->colptr,             mess_int_t*, sizeof(mess_int_t)*(A_cholmod->ncol+1));
        //mess_try_alloc(A->rowptr,             mess_int_t*, sizeof(mess_int_t)*(A_cholmod->nzmax));
        //mess_try_alloc(A->values_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*(A_cholmod->nzmax));

        //column pointer
        LONG *colptr = A_cholmod->p;

        for (i = 0 ; i <= A_cholmod->ncol; i++){
            A->colptr[i]=colptr[i];
        }

        //rowptr
        LONG *rowptr = A_cholmod->i;

        //dataptr
        mess_double_cpx_t  *x  = A_cholmod->x;
        if(A_cholmod->packed){
            for (i = 0; i < A_cholmod->ncol; i++){
                for (j = colptr[i]; j<colptr[i+1]; j++){
                    A->rowptr[j]=rowptr[j];
                    A->values_cpx[j]=x[j];
                }
            }
        }else{
            LONG *nz = A_cholmod->nz;
            for (i = 0; i < A_cholmod->ncol; i++){
                for (j = colptr[i]; j<nz[i]; j++){
                    A->rowptr[j]=rowptr[j];
                    A->values_cpx[j]=x[j];
                }
            }
        }
    }else{
        MSG_ERROR("unsupported cholmod datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}
#endif

