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
 * @file lib/matrix/getcol.c
 * @brief Extract or set a column in a matrix.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


/**
 * @brief Get a column from a matrix.
 * @param[in] matrix    input matrix
 * @param[in] col       input number of the column to get
 * @param[out] c        output vector containing col
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_getcol function gets a column from a matrix. The @matlab equivalent is
 * @verbatim c=A(:,col) @endverbatim
 * If you want to extract a column from a complex matrix the output vector is converted to complex,
 * In the case of a @ref MESS_CSR matrix this operation might
 * be inefficient due to access limitations.
 * The indexing is zero based.
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_setcol
 * @sa mess_matrix_getrow
 *
 */
int mess_matrix_getcol(mess_matrix matrix, mess_int_t col, mess_vector c ){
    MSG_FNAME(__func__);
    mess_int_t i,j;
    int ret  = 0;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(c);
    mess_check_real_or_complex(matrix);
    mess_check_real_or_complex(c);

    if ( col < 0 || col >= matrix->cols) {
        MSG_ERROR("col is out of range. Should be between 0 and " MESS_PRINTF_INT ".\n", matrix->cols);
        return (MESS_ERROR_DIMENSION);
    }

    if ( c->dim != matrix->rows) {
        ret = mess_vector_resize(c, matrix->rows); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_resize);
    }

    /*-----------------------------------------------------------------------------
     *  dense matrices
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_DENSE (matrix)){
        if ( MESS_IS_REAL(matrix)){
            if ( MESS_IS_REAL(c) ) {
                for ( i = 0; i < matrix->rows; i++) {   c->values[i] = matrix->values[matrix->ld*col+i]; }
            } else if ( MESS_IS_COMPLEX(c)){
                for ( i = 0; i < matrix->rows; i++) {   c->values_cpx[i] = matrix->values[matrix->ld*col+i]; }
            }
        } else if ( MESS_IS_COMPLEX(matrix)){
            ret = mess_vector_tocomplex(c); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            for ( i = 0; i < matrix->rows; i++) {
                c->values_cpx[i] = matrix->values_cpx[matrix->ld*col+i];
            }
        }
    }


    /*-----------------------------------------------------------------------------
     *  CSC can used directly
     *-----------------------------------------------------------------------------*/
    else if ( MESS_IS_CSC(matrix)) {
        if ( MESS_IS_REAL(c)){
            for ( i = 0; i < matrix->rows; i++) c->values[i] = 0;
        } else if ( MESS_IS_COMPLEX(c)){
            for ( i = 0; i < matrix->rows; i++) c->values_cpx[i] = 0;
        }
        if ( MESS_IS_REAL(matrix)){
            if ( MESS_IS_REAL(c)) {
                for ( i = matrix->colptr[col]; i < matrix->colptr[col+1]; i++){ c->values[matrix->rowptr[i]]=matrix->values[i]; }
            } else if ( MESS_IS_COMPLEX(c)) {
                for ( i = matrix->colptr[col]; i < matrix->colptr[col+1]; i++){ c->values_cpx[matrix->rowptr[i]]=matrix->values[i];}
            }
        } else if ( MESS_IS_COMPLEX(matrix)){
            ret = mess_vector_tocomplex(c); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            for ( i = matrix->colptr[col]; i < matrix->colptr[col+1]; i++){ c->values_cpx[matrix->rowptr[i]]=matrix->values_cpx[i];}
        }

        /*-----------------------------------------------------------------------------
         *  CSR
         *-----------------------------------------------------------------------------*/
    } else if (MESS_IS_CSR(matrix)) {
        if ( MESS_IS_REAL(c)){
            for ( i = 0; i < matrix->rows; i++) c->values[i] = 0;
        } else if ( MESS_IS_COMPLEX(c)){
            for ( i = 0; i < matrix->rows; i++) c->values_cpx[i] = 0;
        }
        if ( MESS_IS_REAL(matrix)){
            if ( MESS_IS_REAL(c)) {
                for ( i = 0; i < matrix->rows; i++) {
                    for ( j = matrix->rowptr[i]; j < matrix->rowptr[i+1]; j++){
                        if ( matrix->colptr[j] == col ) {
                            c->values[i] = matrix->values[j];
                        }
                    }
                }
            } else if ( MESS_IS_COMPLEX(c)) {
                for ( i = 0; i < matrix->rows; i++) {
                    for ( j = matrix->rowptr[i]; j < matrix->rowptr[i+1]; j++){
                        if ( matrix->colptr[j] == col ) {
                            c->values_cpx[i] = matrix->values[j];
                        }
                    }
                }
            }
        } else if ( MESS_IS_COMPLEX(matrix)){
            ret = mess_vector_tocomplex(c); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            for ( i = 0; i < matrix->rows; i++) {
                for ( j = matrix->rowptr[i]; j < matrix->rowptr[i+1]; j++){
                    if ( matrix->colptr[j] == col ) {
                        c->values_cpx[i] = matrix->values_cpx[j];
                    }
                }
            }
        }
    } else {
        MSG_ERROR("Unsupported/Unknown Storage type : %s\n", mess_storage_t_str(matrix->store_type));
        return(MESS_ERROR_STORAGETYPE);
    }

    return 0;
}


/**
 * @brief Set a column in a matrix.
 * @param[in,out] matrix    input/output matrix where to set the column
 * @param[in] col           input number of the column to set
 * @param[in] colv          input vector to fill the column
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_setcol function overwrites a column of a matrix with a given vector.
 * In the case that you want to  copy a complex column in a real matrix, the complex information is lost and only
 * the real part is copied.
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_getcol
 * @sa mess_matrix_setrow
 */
int mess_matrix_setcol ( mess_matrix matrix, mess_int_t col, mess_vector colv )
{
    MSG_FNAME(__func__);
    mess_int_t i=0;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(colv);
    mess_check_real_or_complex(matrix);
    mess_check_real_or_complex(colv);


    if ( col < 0) {
        MSG_ERROR("col has an illegal values: " MESS_PRINTF_INT "\n", col);
        return ( MESS_ERROR_DIMENSION);
    }
    if ( colv->dim < matrix->rows) {
        MSG_ERROR("colv is too short to fit in a column (colv->dim=" MESS_PRINTF_INT ", matrix->rows= " MESS_PRINTF_INT ")\n", colv->dim, matrix->rows);
        return(MESS_ERROR_DIMENSION);

    }
    if ( colv->dim > matrix->rows ){
        MSG_ERROR("colv is too long to fit in a column (colv->dim=" MESS_PRINTF_INT ", matrix->rows= " MESS_PRINTF_INT ")\n", colv->dim, matrix->rows);
        return(MESS_ERROR_DIMENSION);
    }
    if ( col >= matrix->cols){
        MSG_WARN("resize matrix\n");
        ret = mess_matrix_resize(matrix, matrix->rows, col+1);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_resize);
    }

    if(MESS_IS_DENSE(matrix)){

        /*-----------------------------------------------------------------------------
         *  complex case
         *-----------------------------------------------------------------------------*/
        if ( MESS_IS_COMPLEX(matrix)){
            if ( MESS_IS_COMPLEX(colv)){
                for (i=0; i < colv->dim;i++){ matrix->values_cpx[col*matrix->ld+i] = colv->values_cpx[i];}
            } else if ( MESS_IS_REAL(colv)) {
                for (i=0; i < colv->dim;i++){ matrix->values_cpx[col*matrix->ld+i] = colv->values[i]; }
            }

            /*-----------------------------------------------------------------------------
             *  real case
             *-----------------------------------------------------------------------------*/
        } else if (MESS_IS_REAL(matrix)) {
            if ( MESS_IS_COMPLEX(colv)){
                MSG_WARN("Copy a complex column in a real matrix, the imaginary information will be lost.\n");
                for (i=0; i<colv->dim;i++){ matrix->values[col*matrix->ld+i] = creal(colv->values_cpx[i]); }
            } else if ( MESS_IS_REAL(colv)) {
                for (i=0; i<colv->dim;i++){ matrix->values[col*matrix->ld+i] = colv->values[i]; }
            }
        }
    }
    /*-----------------------------------------------------------------------------
     *  csr case
     *-----------------------------------------------------------------------------*/
    else if(MESS_IS_CSR(matrix)){
        if (MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(colv)){
            MSG_WARN("A complex vector is set to a col of a real matrix. The imaginary part will be lost.\n");
        }
        double * newvalptr = NULL;
        mess_double_cpx_t * newvalptrc = NULL;
        mess_int_t* newcolptr, *newrowptr;
        mess_int_t  j, ind=0, nonzero;

        if (MESS_IS_REAL(matrix)) {
            mess_try_alloc(newvalptr,double *,sizeof(double )*(matrix->nnz+colv->dim));
        } else if ( MESS_IS_COMPLEX(matrix)) {
            mess_try_alloc(newvalptrc,mess_double_cpx_t *,sizeof(mess_double_cpx_t)*(matrix->nnz+colv->dim));
        }
        mess_try_alloc(newcolptr, mess_int_t* , sizeof(mess_int_t)*(matrix->nnz+colv->dim));
        mess_try_alloc(newrowptr, mess_int_t* , sizeof(mess_int_t)*(matrix->rows+1));

        //set initial value
        newrowptr[0]=0;
        ind = 0;
        for(j=0;j<matrix->rows;++j){
            nonzero=0;
            if ( MESS_IS_REAL(colv)) {
                if(colv->values[j] != 0.0 ){
                    nonzero=1;
                }
            } else if ( MESS_IS_COMPLEX(colv)) {
                if (colv->values_cpx[j]!=0.0) {
                    nonzero=1;
                }
            }

            for(i=matrix->rowptr[j];i<matrix->rowptr[j+1];++i){
                if(matrix->colptr[i]<col){
                    //copy elemnts before row to set
                    if ( MESS_IS_COMPLEX(matrix)) {
                        newvalptrc[ind]=matrix->values_cpx[i];
                    } else if ( MESS_IS_REAL(matrix)) {
                        newvalptr[ind]=matrix->values[i];
                    }
                    newcolptr[ind]=matrix->colptr[i];
                    ++ind;
                }
                else if(col==matrix->colptr[i]){
                    //if matrix(k,l)~=0 and col(k)~=0 set matrix(k,l)=row(k)
                    if ( nonzero ) {
                        if ( MESS_IS_COMPLEX(matrix) && MESS_IS_COMPLEX (colv)) {
                            newvalptrc[ind]=colv->values_cpx[j];
                        } else if ( MESS_IS_COMPLEX(matrix) && MESS_IS_REAL(colv)) {
                            newvalptrc[ind]=colv->values[j];
                        } else if ( MESS_IS_REAL(matrix) && MESS_IS_REAL(colv)) {
                            newvalptr[ind]=colv->values[j];
                        } else if ( MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(colv)){
                            newvalptr[ind]=creal(colv->values_cpx[j]);
                        }
                        newcolptr[ind]=col;
                        ++ind;
                        nonzero=0;
                    }
                }
                else if (matrix->colptr[i]>col){
                    if(nonzero){
                        //elements after row to set
                        if ( MESS_IS_COMPLEX(matrix) && MESS_IS_COMPLEX (colv)) {
                            newvalptrc[ind]=colv->values_cpx[j];
                        } else if ( MESS_IS_COMPLEX(matrix) && MESS_IS_REAL(colv)) {
                            newvalptrc[ind]=colv->values[j];
                        } else if ( MESS_IS_REAL(matrix) && MESS_IS_REAL(colv)) {
                            newvalptr[ind]=colv->values[j];
                        } else if ( MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(colv)){
                            newvalptr[ind]=creal(colv->values_cpx[j]);
                        }
                        newcolptr[ind]=col;
                        ++ind;
                        nonzero=0;
                    }
                    if ( MESS_IS_COMPLEX(matrix)) {
                        newvalptrc[ind]=matrix->values_cpx[i];
                    } else if ( MESS_IS_REAL(matrix)) {
                        newvalptr[ind]=matrix->values[i];
                    }
                    newcolptr[ind]=matrix->colptr[i];
                    ++ind;
                }
            }

            if(nonzero){
                if ( MESS_IS_COMPLEX(matrix) && MESS_IS_COMPLEX (colv)) {
                    newvalptrc[ind]=colv->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(matrix) && MESS_IS_REAL(colv)) {
                    newvalptrc[ind]=colv->values[j];
                } else if ( MESS_IS_REAL(matrix) && MESS_IS_REAL(colv)) {
                    newvalptr[ind]=colv->values[j];
                } else if ( MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(colv)){
                    newvalptr[ind]=creal(colv->values_cpx[j]);
                }
                newcolptr[ind]=col;
                ++ind;
            }
            //number of elements in row
            newrowptr[j+1]=ind;
        }

        //set nnz
        matrix->nnz=ind;

        //resize arrays
        if ( MESS_IS_REAL(matrix)) {
            mess_try_realloc(newvalptr,double* ,sizeof(double)*matrix->nnz);
        } else if (MESS_IS_COMPLEX(matrix)) {
            mess_try_realloc(newvalptrc,mess_double_cpx_t* ,sizeof(mess_double_cpx_t)*matrix->nnz);
        }
        mess_try_realloc(newcolptr,mess_int_t*,sizeof(mess_int_t)*matrix->nnz);

        //set new pointer
        mess_free(matrix->colptr);
        mess_free(matrix->rowptr);
        if ( matrix->values) mess_free(matrix->values);
        if ( matrix->values_cpx) mess_free(matrix->values_cpx);

        matrix->rowptr=newrowptr;
        matrix->colptr=newcolptr;
        matrix->values=newvalptr;
        matrix->values_cpx=newvalptrc;
    }else if(MESS_IS_CSC(matrix)){
        /*-----------------------------------------------------------------------------
         *  complex case
         *-----------------------------------------------------------------------------*/
        if ( MESS_IS_COMPLEX(matrix)){
            if ( MESS_IS_COMPLEX(colv)){

                mess_double_cpx_t* newvalptr;
                mess_int_t * newrowptr;
                mess_int_t j, colnnz=0, diffcolnnz=0;
                mess_try_alloc(newvalptr, mess_double_cpx_t* , sizeof(mess_double_cpx_t)*(matrix->nnz+colv->dim));
                mess_try_alloc(newrowptr, mess_int_t* , sizeof(mess_int_t)*(matrix->nnz+colv->dim));

                //copy values before row to set
                for(i=0; i<matrix->colptr[col]; ++i){
                    newvalptr[i]=matrix->values_cpx[i];
                    newrowptr[i]=matrix->rowptr[i];
                }

                //set new row
                for(j=0; j<colv->dim;++j){
                    if(colv->values_cpx[j]){
                        newvalptr[i]=colv->values_cpx[j];
                        newrowptr[i]=j;
                        ++i;
                        ++colnnz;
                    }
                }

                //set rows after the row to set
                for(j=matrix->colptr[col+1];j<matrix->colptr[matrix->cols];++j){
                    newvalptr[i]=matrix->values_cpx[j];
                    newrowptr[i]=matrix->rowptr[j];
                    ++i;
                }

                //set nnz elements
                diffcolnnz=colnnz-(matrix->colptr[col+1]-matrix->colptr[col]);
                matrix->nnz+=diffcolnnz;

                //set new rowptr
                for(j=col+1;j<=matrix->cols;++j){
                    matrix->colptr[j]+=diffcolnnz;
                }

                //resize arrays
                mess_try_realloc(newvalptr, mess_double_cpx_t*,sizeof(mess_double_cpx_t)*matrix->nnz);
                mess_try_realloc(newrowptr, mess_int_t* , sizeof(mess_int_t)*matrix->nnz);

                //delete old pointers
                mess_free(matrix->values_cpx);
                mess_free(matrix->rowptr);

                //set new pointers
                matrix->values_cpx=newvalptr;
                matrix->rowptr = newrowptr;

            } else if ( MESS_IS_REAL(colv)) {

                mess_double_cpx_t* newvalptr;
                mess_int_t * newrowptr;
                mess_int_t j, colnnz=0, diffcolnnz=0;
                mess_try_alloc(newvalptr, mess_double_cpx_t* , sizeof(mess_double_cpx_t)*(matrix->nnz+colv->dim));
                mess_try_alloc(newrowptr, mess_int_t* , sizeof(mess_int_t)*(matrix->nnz+colv->dim));

                //copy values before row to set
                for(i=0; i<matrix->colptr[col]; ++i){
                    newvalptr[i]=matrix->values_cpx[i];
                    newrowptr[i]=matrix->rowptr[i];
                }

                //set new row
                for(j=0; j<colv->dim;++j){
                    if(colv->values[j]){
                        newvalptr[i]=colv->values[j];
                        newrowptr[i]=j;
                        ++i;
                        ++colnnz;
                    }
                }

                //set rows after the row to set
                for(j=matrix->colptr[col+1];j<matrix->colptr[matrix->cols];++j){
                    newvalptr[i]=matrix->values_cpx[j];
                    newrowptr[i]=matrix->rowptr[j];
                    ++i;
                }

                //set nnz elements
                diffcolnnz=colnnz-(matrix->colptr[col+1]-matrix->colptr[col]);
                matrix->nnz+=diffcolnnz;

                //set new rowptr
                for(j=col+1;j<=matrix->cols;++j){
                    matrix->colptr[j]+=diffcolnnz;
                }

                //resize arrays
                mess_try_realloc(newvalptr, mess_double_cpx_t*,sizeof(mess_double_cpx_t)*matrix->nnz);
                mess_try_realloc(newrowptr, mess_int_t* , sizeof(mess_int_t)*matrix->nnz);

                //delete old pointers
                mess_free(matrix->values_cpx);
                mess_free(matrix->rowptr);

                //set new pointers
                matrix->values_cpx=newvalptr;
                matrix->rowptr = newrowptr;


            }
            /*-----------------------------------------------------------------------------
             *  real case
             *-----------------------------------------------------------------------------*/
        } else if (MESS_IS_REAL(matrix)) {
            if ( MESS_IS_COMPLEX(colv)){
                MSG_WARN("A complex vector is set to a col of a real matrix. The imaginary part will be lost.\n");

                double * newvalptr;
                mess_int_t * newrowptr;
                mess_int_t j, colnnz=0, diffcolnnz=0;
                mess_try_alloc(newvalptr, double * , sizeof(double )*(matrix->nnz+colv->dim));
                mess_try_alloc(newrowptr, mess_int_t* , sizeof(mess_int_t)*(matrix->nnz+colv->dim));

                //copy values before row to set
                for(i=0; i<matrix->colptr[col]; ++i){
                    newvalptr[i]=matrix->values[i];
                    newrowptr[i]=matrix->rowptr[i];
                }

                //set new row
                for(j=0; j<colv->dim;++j){
                    double val = creal(colv->values_cpx[j]);
                    if(val){
                        newvalptr[i]=val;
                        newrowptr[i]=j;
                        ++i;
                        ++colnnz;
                    }
                }

                //set rows after the row to set
                for(j=matrix->colptr[col+1];j<matrix->colptr[matrix->cols];++j){
                    newvalptr[i]=matrix->values[j];
                    newrowptr[i]=matrix->rowptr[j];
                    ++i;
                }

                //set nnz elements
                diffcolnnz=colnnz-(matrix->colptr[col+1]-matrix->colptr[col]);
                matrix->nnz+=diffcolnnz;

                //set new rowptr
                for(j=col+1;j<=matrix->cols;++j){
                    matrix->colptr[j]+=diffcolnnz;
                }

                //resize arrays
                mess_try_realloc(newvalptr, double *,sizeof(double )*matrix->nnz);
                mess_try_realloc(newrowptr, mess_int_t* , sizeof(mess_int_t)*matrix->nnz);

                //delete old pointers
                mess_free(matrix->values);
                mess_free(matrix->rowptr);

                //set new pointers
                matrix->values=newvalptr;
                matrix->rowptr = newrowptr;

            } else if ( MESS_IS_REAL(colv)) {

                double * newvalptr;
                mess_int_t * newrowptr;
                mess_int_t j, colnnz=0, diffcolnnz=0;
                mess_try_alloc(newvalptr, double * , sizeof(double )*(matrix->nnz+colv->dim));
                mess_try_alloc(newrowptr, mess_int_t* , sizeof(mess_int_t)*(matrix->nnz+colv->dim));

                //copy values before row to set
                for(i=0; i<matrix->colptr[col]; ++i){
                    newvalptr[i]=matrix->values[i];
                    newrowptr[i]=matrix->rowptr[i];
                }

                //set new row
                for(j=0; j<colv->dim;++j){
                    if(colv->values[j]){
                        newvalptr[i]=colv->values[j];
                        newrowptr[i]=j;
                        ++i;
                        ++colnnz;
                    }
                }

                //set rows after the row to set
                for(j=matrix->colptr[col+1];j<matrix->colptr[matrix->cols];++j){
                    newvalptr[i]=matrix->values[j];
                    newrowptr[i]=matrix->rowptr[j];
                    ++i;
                }

                //set nnz elements
                diffcolnnz=colnnz-(matrix->colptr[col+1]-matrix->colptr[col]);
                matrix->nnz+=diffcolnnz;

                //set new rowptr
                for(j=col+1;j<=matrix->cols;++j){
                    matrix->colptr[j]+=diffcolnnz;
                }

                //resize arrays
                mess_try_realloc(newvalptr, double *,sizeof(double )*matrix->nnz);
                mess_try_realloc(newrowptr, mess_int_t* , sizeof(mess_int_t)*matrix->nnz);

                //delete old pointers
                mess_free(matrix->values);
                mess_free(matrix->rowptr);

                //set new pointers
                matrix->values=newvalptr;
                matrix->rowptr = newrowptr;

            }
        }
    }else{
        MSG_ERROR("Unsupported Storagetype: %s \n", mess_storage_t_str(matrix->store_type));
        return(MESS_ERROR_STORAGETYPE);
    }
    return (0);
}

