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
 * @file lib/matrix/getrow.c
 * @brief Extract or set a row out of a matrix.
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
 * @brief Get a row from a matrix.
 * @param[in] matrix    input matrix
 * @param[in] row       input number of the row to get
 * @param[out] r        output vector containing row
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_getrow function gets a the @c row from a matrix. The @matlab equivalent is @verbatim r=A(row,:)@endverbatim.
 * 
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_setrow
 * @sa mess_matrix_getcol
 *
 */
int mess_matrix_getrow(mess_matrix matrix, mess_int_t row, mess_vector r ){
    MSG_FNAME(__func__);
    mess_int_t i,j;
    int ret = 0 ;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(r);

    if ( row < 0 || row >= matrix->rows) {
        MSG_ERROR("row is out of range. Should be between 0 and " MESS_PRINTF_INT ".\n", matrix->rows);
        return(MESS_ERROR_DIMENSION);
    }
    if ( r->dim != matrix->cols) {
        ret = mess_vector_resize(r, matrix->cols); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }

    if (MESS_IS_REAL(matrix)) {
        if ( MESS_IS_DENSE (matrix)){
            if ( MESS_IS_REAL(r)) {
                for ( i = 0; i < matrix->cols; i++) {
                    r->values[i] = matrix->values[i*matrix->ld+row];
                }
            } else {
                for ( i = 0; i < matrix->cols; i++) {
                    r->values_cpx[i] = matrix->values[i*matrix->ld+row];
                }
            }
        }  else if ( MESS_IS_CSC(matrix) ) {
            if (MESS_IS_REAL(r)) {
                for ( i = 0 ; i < matrix->cols; i++ ) {
                    r->values[i] = 0.0 ;
                    for ( j = matrix->colptr[i]; j < matrix->colptr[i+1]; j++ ) {
                        if ( matrix->rowptr[j] == row ) {
                            r->values[i] = matrix->values[j];
                        }
                    }
                }

            } else {
                for ( i = 0 ; i < matrix->cols; i++ ) {
                    r->values_cpx[i] = 0.0;
                    for ( j = matrix->colptr[i]; j < matrix->colptr[i+1]; j++ ) {
                        if ( matrix->rowptr[j] == row ) {
                            r->values_cpx[i] = matrix->values[j];
                        }
                    }
                }
            }
        }  else if ( MESS_IS_CSR(matrix)) {
            if (MESS_IS_REAL(r)) {
                for ( i = 0; i < matrix->cols; i++) r->values[i] = 0;
                for ( i = matrix->rowptr[row]; i < matrix->rowptr[row+1]; i++) r->values[matrix->colptr[i]]=matrix->values[i];
            } else {
                for ( i = 0; i < matrix->cols; i++) r->values_cpx[i] = 0;
                for ( i = matrix->rowptr[row]; i < matrix->rowptr[row+1]; i++) r->values_cpx[matrix->colptr[i]]=matrix->values[i];

            }
        } else {
            MSG_ERROR("Unsupported Storage type: %s \n", mess_storage_t_str(matrix->store_type));
            return(MESS_ERROR_STORAGETYPE);
        }
    } else if ( MESS_IS_COMPLEX(matrix)){
        ret = mess_vector_tocomplex(r);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
        if ( MESS_IS_DENSE (matrix)){
            for ( i = 0; i < matrix->cols; i++) {
                r->values_cpx[i] = matrix->values_cpx[i*matrix->ld+row];
            }
        }  else if ( MESS_IS_CSC(matrix) ) {
            for ( i = 0 ; i < matrix->cols; i++ ) {
                r->values_cpx[i] = 0.0;
                for ( j = matrix->colptr[i]; j < matrix->colptr[i+1]; j++ ) {
                    if ( matrix->rowptr[j] == row ) {
                        r->values_cpx[i] = matrix->values_cpx[j];
                    }
                }
            }
        }  else if ( MESS_IS_CSR(matrix)) {
            for ( i = 0; i < matrix->cols; i++) r->values_cpx[i] = 0;
            for ( i = matrix->rowptr[row]; i < matrix->rowptr[row+1]; i++) r->values_cpx[matrix->colptr[i]]=matrix->values_cpx[i];
        } else {
            MSG_ERROR("Unsupported Storage type: %s \n", mess_storage_t_str(matrix->store_type));
            return(MESS_ERROR_STORAGETYPE);
        }

    }
    return(0);
}


/**
 * @brief Set a row in a matrix.
 * @param[in,out] matrix    input/output matrix where to set the column
 * @param[in] row           input number of the column to set
 * @param[in] rowv          input vector to fill the column
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_setrow function overwrites a row of a matrix with a given vector.
 * 
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_setcol
 * @sa mess_matrix_getrow
 */
int mess_matrix_setrow ( mess_matrix matrix, mess_int_t row, mess_vector rowv )
{
    MSG_FNAME(__func__);
    mess_int_t i=0;
    int ret = 0 ;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(rowv);
    mess_check_real_or_complex(matrix);
    mess_check_real_or_complex(rowv);

    if ( row < 0) {
        MSG_ERROR("row has an illegal values: " MESS_PRINTF_INT "\n", row);
        return(MESS_ERROR_DIMENSION);
    }
    if ( rowv->dim < matrix->cols) {
        MSG_ERROR("rowv is too short to fit in a column (rowv->dim=" MESS_PRINTF_INT ", matrix->cols= " MESS_PRINTF_INT ")\n", rowv->dim, matrix->rows);
        return(MESS_ERROR_DIMENSION);
    }
    if ( rowv->dim > matrix->cols ){
        MSG_ERROR("rowv is too mess_int_t to fit in a column (rowv->dim=" MESS_PRINTF_INT ", matrix->rows= " MESS_PRINTF_INT ")\n", rowv->dim, matrix->rows);
        return(MESS_ERROR_DIMENSION);
    }
    if ( row >= matrix->rows){
        MSG_WARN("resize matrix");
        ret = mess_matrix_resize(matrix, row+1, matrix->cols);  FUNCTION_FAILURE_HANDLE(ret, (ret != 0 ), mess_matrix_resize);
    }



    /*-----------------------------------------------------------------------------
     *  dense case
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_DENSE(matrix)){
        /*-----------------------------------------------------------------------------
         *  complex case
         *-----------------------------------------------------------------------------*/
        if ( MESS_IS_COMPLEX(matrix)){
            if ( MESS_IS_COMPLEX(rowv)){
                for (i=0; i < rowv->dim;i++){
                    matrix->values_cpx[i*matrix->ld+row] = rowv->values_cpx[i];
                }
            } else if ( MESS_IS_REAL(rowv)) {
                for (i=0; i < rowv->dim;i++){
                    matrix->values_cpx[i*matrix->ld+row] = rowv->values[i];
                }
            }
            /*-----------------------------------------------------------------------------
             *  real case
             *-----------------------------------------------------------------------------*/
        } else if (MESS_IS_REAL(matrix)) {
            if ( MESS_IS_COMPLEX(rowv)){
                MSG_WARN("A complex vector is set to a row of a real matrix. The imaginary part will be lost.\n");
                for (i=0; i<rowv->dim;i++){
                    matrix->values[i*matrix->ld+row] = creal(rowv->values_cpx[i]);
                }
            } else if ( MESS_IS_REAL(rowv)) {
                for (i=0; i<rowv->dim;i++){
                    matrix->values[i*matrix->ld+row] = rowv->values[i];
                }
            }
        }
    }
    /*-----------------------------------------------------------------------------
     *  csr case
     *-----------------------------------------------------------------------------*/
    else if(MESS_IS_CSR(matrix)){
        /*-----------------------------------------------------------------------------
         *  complex case
         *-----------------------------------------------------------------------------*/
        if ( MESS_IS_COMPLEX(matrix)){
            if ( MESS_IS_COMPLEX(rowv)){

                mess_double_cpx_t* newvalptr;
                mess_int_t * newcolptr;
                mess_int_t j, rownnz=0, diffrownnz=0;
                mess_try_alloc(newvalptr, mess_double_cpx_t* , sizeof(mess_double_cpx_t)*(matrix->nnz+rowv->dim));
                mess_try_alloc(newcolptr, mess_int_t* , sizeof(mess_int_t)*(matrix->nnz+rowv->dim));

                //copy values before row to set
                for(i=0; i<matrix->rowptr[row]; ++i){
                    newvalptr[i]=matrix->values_cpx[i];
                    newcolptr[i]=matrix->colptr[i];
                }

                //set new row
                for(j=0; j<rowv->dim;++j){
                    if(rowv->values_cpx[j]){
                        newvalptr[i]=rowv->values_cpx[j];
                        newcolptr[i]=j;
                        ++i;
                        ++rownnz;
                    }
                }

                //set rows after the row to set
                for(j=matrix->rowptr[row+1];j<matrix->rowptr[matrix->rows];++j){
                    newvalptr[i]=matrix->values_cpx[j];
                    newcolptr[i]=matrix->colptr[j];
                    ++i;
                }

                //set nnz elements
                diffrownnz=rownnz-(matrix->rowptr[row+1]-matrix->rowptr[row]);
                matrix->nnz+=diffrownnz;

                //set new rowptr
                for(j=row+1;j<=matrix->rows;++j){
                    matrix->rowptr[j]+=diffrownnz;
                }

                //resize arrays
                mess_try_realloc(newvalptr, mess_double_cpx_t*,sizeof(mess_double_cpx_t)*matrix->nnz);
                mess_try_realloc(newcolptr, mess_int_t* , sizeof(mess_int_t)*matrix->nnz);

                //delete old pointers
                mess_free(matrix->values_cpx);
                mess_free(matrix->colptr);

                //set new pointers
                matrix->values_cpx=newvalptr;
                matrix->colptr = newcolptr;

            } else if ( MESS_IS_REAL(rowv)) {

                mess_double_cpx_t* newvalptr;
                mess_int_t * newcolptr;
                mess_int_t j, rownnz=0, diffrownnz=0;
                mess_try_alloc(newvalptr, mess_double_cpx_t* , sizeof(mess_double_cpx_t)*(matrix->nnz+rowv->dim));
                mess_try_alloc(newcolptr, mess_int_t* , sizeof(mess_int_t)*(matrix->nnz+rowv->dim));

                //copy values before row to set
                for(i=0; i<matrix->rowptr[row]; ++i){
                    newvalptr[i]=matrix->values_cpx[i];
                    newcolptr[i]=matrix->colptr[i];
                }

                //set new row
                for(j=0; j<rowv->dim;++j){
                    if(rowv->values[j]){
                        newvalptr[i]=rowv->values[j];
                        newcolptr[i]=j;
                        ++i;
                        ++rownnz;
                    }
                }

                //set rows after the row to set
                for(j=matrix->rowptr[row+1];j<matrix->rowptr[matrix->rows];++j){
                    newvalptr[i]=matrix->values_cpx[j];
                    newcolptr[i]=matrix->colptr[j];
                    ++i;
                }

                //set nnz elements
                diffrownnz=rownnz-(matrix->rowptr[row+1]-matrix->rowptr[row]);
                matrix->nnz+=diffrownnz;

                //set new rowptr
                for(j=row+1;j<=matrix->rows;++j){
                    matrix->rowptr[j]+=diffrownnz;
                }

                //resize arrays
                mess_try_realloc(newvalptr, mess_double_cpx_t*,sizeof(mess_double_cpx_t)*matrix->nnz);
                mess_try_realloc(newcolptr, mess_int_t* , sizeof(mess_int_t)*matrix->nnz);

                //delete old pointers
                mess_free(matrix->values_cpx);
                mess_free(matrix->colptr);

                //set new pointers
                matrix->values_cpx=newvalptr;
                matrix->colptr = newcolptr;



            }
            /*-----------------------------------------------------------------------------
             *  real case
             *-----------------------------------------------------------------------------*/
        } else if (MESS_IS_REAL(matrix)) {
            if ( MESS_IS_COMPLEX(rowv)){
                MSG_WARN("A complex vector is set to a row of a real matrix. The imaginary part will be lost.\n");
                double * newvalptr;
                mess_int_t * newcolptr;
                mess_int_t j, rownnz=0, diffrownnz=0;
                mess_try_alloc(newvalptr, double* , sizeof(double)*(matrix->nnz+rowv->dim));
                mess_try_alloc(newcolptr, mess_int_t* , sizeof(mess_int_t)*(matrix->nnz+rowv->dim));

                //copy values before row to set
                for(i=0; i<matrix->rowptr[row]; ++i){
                    newvalptr[i]=matrix->values[i];
                    newcolptr[i]=matrix->colptr[i];
                }

                //set new row
                double val;
                for(j=0; j<rowv->dim;++j){
                    val = creal(rowv->values_cpx[j]);
                    if(val){
                        newvalptr[i]=val;
                        newcolptr[i]=j;
                        ++i;
                        ++rownnz;
                    }
                }

                //set rows after the row to set
                for(j=matrix->rowptr[row+1];j<matrix->rowptr[matrix->rows];++j){
                    newvalptr[i]=matrix->values[j];
                    newcolptr[i]=matrix->colptr[j];
                    ++i;
                }

                //set nnz elements
                diffrownnz=rownnz-(matrix->rowptr[row+1]-matrix->rowptr[row]);
                matrix->nnz+=diffrownnz;

                //set new rowptr
                for(j=row+1;j<=matrix->rows;++j){
                    matrix->rowptr[j]+=diffrownnz;
                }

                //resize arrays
                mess_try_realloc(newvalptr, double* ,sizeof(double)*matrix->nnz);
                mess_try_realloc(newcolptr, mess_int_t* , sizeof(mess_int_t)*matrix->nnz);

                //delete old pointers
                mess_free(matrix->values);
                mess_free(matrix->colptr);

                //set new pointers
                matrix->values=newvalptr;
                matrix->colptr = newcolptr;
            }
            else if ( MESS_IS_REAL(rowv)) {
                double * newvalptr;
                mess_int_t * newcolptr;
                mess_int_t j, rownnz=0, diffrownnz=0;
                mess_try_alloc(newvalptr, double * , sizeof(double)*(matrix->nnz+rowv->dim));
                mess_try_alloc(newcolptr, mess_int_t* , sizeof(mess_int_t)*(matrix->nnz+rowv->dim));

                //copy values before row to set
                for(i=0; i<matrix->rowptr[row]; ++i){
                    newvalptr[i]=matrix->values[i];
                    newcolptr[i]=matrix->colptr[i];
                }

                //set new row
                for(j=0; j<rowv->dim;++j){
                    if(rowv->values[j]){
                        newvalptr[i]=rowv->values[j];
                        newcolptr[i]=j;
                        ++i;
                        ++rownnz;
                    }
                }

                //set rows after the row to set
                for(j=matrix->rowptr[row+1];j<matrix->rowptr[matrix->rows];++j){
                    newvalptr[i]=matrix->values[j];
                    newcolptr[i]=matrix->colptr[j];
                    ++i;
                }

                //set nnz elements
                diffrownnz=rownnz-(matrix->rowptr[row+1]-matrix->rowptr[row]);
                matrix->nnz+=diffrownnz;

                //set new rowptr
                for(j=row+1;j<=matrix->rows;++j){
                    matrix->rowptr[j]+=diffrownnz;
                }

                //resize arrays
                mess_try_realloc(newvalptr, double *,sizeof(double)*matrix->nnz);
                mess_try_realloc(newcolptr, mess_int_t* , sizeof(mess_int_t)*matrix->nnz);

                //delete old pointers
                mess_free(matrix->values);
                mess_free(matrix->colptr);

                //set new pointers
                matrix->values=newvalptr;
                matrix->colptr = newcolptr;
            }
        }


    }
    /*-----------------------------------------------------------------------------
     *  csc case
     *-----------------------------------------------------------------------------*/
    else if(MESS_IS_CSC(matrix)){
        if (MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(rowv)){
            MSG_WARN("A complex vector is set to a col of a real matrix. The imaginary part will be lost.\n");
        }
        double * newvalptr = NULL;
        mess_double_cpx_t * newvalptrc = NULL;
        mess_int_t* newcolptr, *newrowptr;
        mess_int_t j, ind=0, nonzero;

        if (MESS_IS_REAL(matrix)) {
            mess_try_alloc(newvalptr,double *,sizeof(double )*(matrix->nnz+rowv->dim));
        } else if ( MESS_IS_COMPLEX(matrix)) {
            mess_try_alloc(newvalptrc,mess_double_cpx_t *,sizeof(mess_double_cpx_t)*(matrix->nnz+rowv->dim));
        }
        mess_try_alloc(newrowptr, mess_int_t* , sizeof(mess_int_t)*(matrix->nnz+rowv->dim));
        mess_try_alloc(newcolptr, mess_int_t* , sizeof(mess_int_t)*(matrix->cols+1));

        //set initial value
        newcolptr[0]=0;
        ind = 0;
        for(j=0;j<matrix->cols;++j){
            nonzero=0;
            if ( MESS_IS_REAL(rowv)) {
                if(rowv->values[j] != 0.0 ){
                    nonzero=1;
                }
            } else if ( MESS_IS_COMPLEX(rowv)) {
                if (rowv->values_cpx[j]!=0.0) {
                    nonzero=1;
                }
            }

            for(i=matrix->colptr[j];i<matrix->colptr[j+1];++i){
                if(matrix->rowptr[i]<row){
                    //copy elemnts before row to set
                    if ( MESS_IS_COMPLEX(matrix)) {
                        newvalptrc[ind]=matrix->values_cpx[i];
                    } else if ( MESS_IS_REAL(matrix)) {
                        newvalptr[ind]=matrix->values[i];
                    }
                    newrowptr[ind]=matrix->rowptr[i];
                    ++ind;
                }
                else if(row==matrix->rowptr[i]){
                    //if matrix(k,l)~=0 and col(k)~=0 set matrix(k,l)=row(k)
                    if ( nonzero ) {
                        if ( MESS_IS_COMPLEX(matrix) && MESS_IS_COMPLEX (rowv)) {
                            newvalptrc[ind]=rowv->values_cpx[j];
                        } else if ( MESS_IS_COMPLEX(matrix) && MESS_IS_REAL(rowv)) {
                            newvalptrc[ind]=rowv->values[j];
                        } else if ( MESS_IS_REAL(matrix) && MESS_IS_REAL(rowv)) {
                            newvalptr[ind]=rowv->values[j];
                        } else if ( MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(rowv)){
                            newvalptr[ind]=creal(rowv->values_cpx[j]);
                        }
                        newrowptr[ind]=row;
                        ++ind;
                        nonzero=0;
                    }
                }
                else if (matrix->rowptr[i]>row){
                    if(nonzero){
                        //elements after row to set
                        if ( MESS_IS_COMPLEX(matrix) && MESS_IS_COMPLEX (rowv)) {
                            newvalptrc[ind]=rowv->values_cpx[j];
                        } else if ( MESS_IS_COMPLEX(matrix) && MESS_IS_REAL(rowv)) {
                            newvalptrc[ind]=rowv->values[j];
                        } else if ( MESS_IS_REAL(matrix) && MESS_IS_REAL(rowv)) {
                            newvalptr[ind]=rowv->values[j];
                        } else if ( MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(rowv)){
                            newvalptr[ind]=creal(rowv->values_cpx[j]);
                        }
                        newrowptr[ind]=row;
                        ++ind;
                        nonzero=0;
                    }
                    if ( MESS_IS_COMPLEX(matrix)) {
                        newvalptrc[ind]=matrix->values_cpx[i];
                    } else if ( MESS_IS_REAL(matrix)) {
                        newvalptr[ind]=matrix->values[i];
                    }
                    newrowptr[ind]=matrix->rowptr[i];
                    ++ind;
                }
            }

            if(nonzero){
                if ( MESS_IS_COMPLEX(matrix) && MESS_IS_COMPLEX (rowv)) {
                    newvalptrc[ind]=rowv->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(matrix) && MESS_IS_REAL(rowv)) {
                    newvalptrc[ind]=rowv->values[j];
                } else if ( MESS_IS_REAL(matrix) && MESS_IS_REAL(rowv)) {
                    newvalptr[ind]=rowv->values[j];
                } else if ( MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(rowv)){
                    newvalptr[ind]=creal(rowv->values_cpx[j]);
                }
                newrowptr[ind]=row;
                ++ind;
            }
            //number of elements in row
            newcolptr[j+1]=ind;
        }

        //set nnz
        matrix->nnz=ind;

        //resize arrays
        if ( MESS_IS_REAL(matrix)) {
            mess_try_realloc(newvalptr,double* ,sizeof(double)*matrix->nnz);
        } else if (MESS_IS_COMPLEX(matrix)) {
            mess_try_realloc(newvalptrc,mess_double_cpx_t* ,sizeof(mess_double_cpx_t)*matrix->nnz);
        }
        mess_try_realloc(newrowptr,mess_int_t*,sizeof(mess_int_t)*matrix->nnz);

        //set new pointer
        mess_free(matrix->colptr);
        mess_free(matrix->rowptr);
        if ( matrix->values) mess_free(matrix->values);
        if ( matrix->values_cpx) mess_free(matrix->values_cpx);

        matrix->rowptr=newrowptr;
        matrix->colptr=newcolptr;
        matrix->values=newvalptr;
        matrix->values_cpx=newvalptrc;

    }else{
        MSG_ERROR("Unsupported Storagetype: %s \n", mess_storage_t_str(matrix->store_type));
        return(MESS_ERROR_STORAGETYPE);
    }
    return(0);
}

