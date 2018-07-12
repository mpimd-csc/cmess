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
 * @file lib/matrix/colops.c
 * @brief Perform operations directly on the columns of a matrix.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>
#include <math.h>






/**
 * @brief Compute the euclidean norm of a matrix column.
 * @param[in] Q     input matrix \f$Q\f$
 * @param[in] col   input index \f$col\f$ of the column
 * @param[out] norm output euclidean norm the desired column of \f$Q\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_colnorm function computes
 * \f[norm =\Vert Q(:,col)\Vert_2 . \f]
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 */
int mess_matrix_colnorm(mess_matrix Q, mess_int_t col, double *norm ) {
    MSG_FNAME(__func__);
    mess_int_t one = 1;
    mess_check_nullpointer(Q);
    mess_check_real_or_complex(Q);

    if ( col <0 || col >= Q->cols) {
        MSG_ERROR("The column index ( " MESS_PRINTF_INT " ) is out of range. \n", col);
        return (MESS_ERROR_ARGUMENTS);
    }
    if ( MESS_IS_DENSE(Q)){
        if ( MESS_IS_REAL(Q)){
            *norm = (F77_GLOBAL(dnrm2,DNRM2)(&(Q->rows), &(Q->values[Q->ld*col]), &one));
        } else if ( MESS_IS_COMPLEX(Q)) {
            *norm = (F77_GLOBAL(dznrm2,DZNRM2)(&(Q->rows), &(Q->values_cpx[Q->ld*col]), &one));
        }
    } else if ( MESS_IS_CSC(Q)){
        if (MESS_IS_REAL(Q)){
            mess_int_t i;
            double nrm = 0;
            for ( i = Q->colptr[col]; i < Q->colptr[col+1]; i++){
                nrm += Q->values[i]*Q->values[i];
            }
            *norm = sqrt(nrm);
        } else {
            mess_int_t i;
            double nrm = 0;
            for ( i = Q->colptr[col]; i < Q->colptr[col+1]; i++){
                nrm += pow(cabs(Q->values_cpx[i]),2);
            }
            *norm = sqrt(nrm);
        }
    } else if  (MESS_IS_CSR(Q)){
        if (MESS_IS_REAL(Q)){
            mess_int_t i,j;
            double nrm = 0;

            for ( i = 0 ; i < Q->rows;i++ ) {
                for ( j = Q->rowptr[i]; j < Q->rowptr[i+1]; j++  ) {
                    if ( Q->colptr[j] == col ) {
                        nrm += Q->values[j]*Q->values[j];
                    }
                }
            }
            *norm = sqrt(nrm);
        } else {
            mess_int_t i,j;
            double nrm = 0;
            for ( i = 0 ; i < Q->rows;i++ ) {
                for ( j = Q->rowptr[i]; j < Q->rowptr[i+1]; j++  ) {
                    if ( Q->colptr[j] == col ) {
                        nrm += pow(cabs(Q->values_cpx[j]),2);
                    }
                }
            }
            *norm = sqrt(nrm);
        }
    } else {
        MSG_ERROR("Unsupported Storage scheme: %s\n", mess_storage_t_str(Q->store_type));
        return MESS_ERROR_STORAGETYPE;
    }
    return (0);
}

/**
 * @brief Scale a column of a matrix.
 * @param[in,out] Q     input/output matrix \f$Q\f$
 * @param[in] col       input index \f$col\f$ of the column
 * @param[in] scale     input factor \f$scale\f$ to scale
 * @return zero on success or a non-zero error value otherwise
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * The @ref mess_matrix_colscale function scales a column of a matrix:
 * \f[ Q(:,col) \leftarrow scale\cdot Q(:,col). \f]
 */
int mess_matrix_colscale(mess_matrix Q, mess_int_t col, mess_double_cpx_t scale) {
    MSG_FNAME(__func__);
    mess_int_t one = 1;
    mess_check_nullpointer(Q);
    mess_check_real_or_complex(Q);
    int ret = 0;

    if ( col <0 || col >= Q->cols) {
        MSG_ERROR("The column index ( " MESS_PRINTF_INT " ) is out of range. \n", col);
        return ( MESS_ERROR_ARGUMENTS);
    }

    if ( MESS_IS_DENSE(Q)){
        if ( MESS_IS_REAL(Q) && cimag(scale) == 0.0 ){
            F77_GLOBAL(dscal,DSCAL)(&(Q->rows), (double *) &scale, &(Q->values[col*Q->ld]), &one);
        } else {
            ret = mess_matrix_tocomplex(Q);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
            F77_GLOBAL(zscal,ZSCAL)(&(Q->rows), (mess_double_cpx_t *) &scale, &(Q->values_cpx[col*Q->ld]), &one);
        }
    } else if ( MESS_IS_CSC(Q)){
        if (MESS_IS_REAL(Q) && cimag(scale) == 0.0 ){
            mess_int_t i;
            for ( i = Q->colptr[col]; i < Q->colptr[col+1]; i++){
                Q->values[i]*=creal(scale);
            }
        } else {
            mess_int_t i;
            ret = mess_matrix_tocomplex(Q);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
            for ( i = Q->colptr[col]; i < Q->colptr[col+1]; i++){
                Q->values_cpx[i]*=scale;
            }
        }
    } else if  (MESS_IS_CSR(Q)){
        if (MESS_IS_REAL(Q) && cimag(scale) == 0.0 ){
            mess_int_t i,j;
            for ( i = 0 ; i < Q->rows;i++ ) {
                for ( j = Q->rowptr[i]; j < Q->rowptr[i+1]; j++  ) {
                    if ( Q->colptr[j] == col ) {
                        Q->values[j]*=creal(scale);
                    }
                }
            }
        } else {
            mess_int_t i,j;
            ret = mess_matrix_tocomplex(Q);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
            for ( i = 0 ; i < Q->rows;i++ ) {
                for ( j = Q->rowptr[i]; j < Q->rowptr[i+1]; j++  ) {
                    if ( Q->colptr[j] == col ) {
                        Q->values_cpx[j]*=scale;
                    }
                }
            }
        }
    } else {
        MSG_ERROR("Unsupported Storage scheme: %s\n", mess_storage_t_str(Q->store_type));
        return MESS_ERROR_STORAGETYPE;
    }
    return (0);
}

/**
 * @brief Compute the dot product between two columns of a real matrix.
 * @param[in] Q     input matrix \f$Q\f$
 * @param[in] col1  input index \f$col_1\f$ of first column
 * @param[in] col2  input index \f$col_2\f$ of second column
 * @param[out] dot  output dot product
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_coldot function computes the dot product between two columns of a real matrix:
 * \f[dot= Q(:,col_1)^TQ(:,col_2).\f]
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_coldotc
 * @sa mess_matrix_coldotE
 * @sa mess_matrix_coldotcE
 * @sa mess_matrix_colvecdot
 * @sa mess_matrix_colvecdotc
 */
int mess_matrix_coldot(mess_matrix Q, mess_int_t col1, mess_int_t col2, double *dot){
    MSG_FNAME(__func__);
    mess_int_t one =1;
    mess_check_nullpointer(Q);
    mess_check_real(Q);

    if ( col1 <0 || col1 >= Q->cols) {
        MSG_ERROR("The column index 1 ( " MESS_PRINTF_INT " ) is out of range. \n", col1);
        return(MESS_ERROR_ARGUMENTS);
    }
    if ( col2 <0 || col2 >= Q->cols) {
        MSG_ERROR("The column index 2 ( " MESS_PRINTF_INT " ) is out of range. \n", col2);
        return(MESS_ERROR_ARGUMENTS);
    }

    if (MESS_IS_DENSE(Q)) {
        *dot = F77_GLOBAL(ddot,DDOT)(&(Q->rows), &(Q->values[col1*Q->ld]), &one, &(Q->values[col2*Q->ld]), &one);
    } else {
        mess_int_t i,j;
        if ( MESS_IS_CSR(Q)){
            mess_int_t setboth;
            double val1, val2;
            //initial value
            *dot=0;
            for ( i = 0 ; i < Q->rows;i++ ) {
                val1=0; val2=0;
                setboth=-1;
                for ( j = Q->rowptr[i]; j < Q->rowptr[i+1]; j++  ) {
                    if ( Q->colptr[j] == col1){
                        val1 = Q->values[j];
                        ++setboth;
                        if(setboth){
                            *dot+=val1*val2;
                            break;
                        }
                    }
                    if( Q->colptr[j] == col2){
                        val2 = Q->values[j];
                        ++setboth;
                        if(setboth){
                            *dot+=val1*val2;
                            break;
                        }
                    }
                }
            }

        } else if ( MESS_IS_CSC(Q)){

            //inital value
            *dot=0;
            //j indices for col2
            j=Q->colptr[col2];
            //i indices for col1
            for(i=Q->colptr[col1];i<Q->colptr[col1+1];++i){
                while(Q->rowptr[j]<Q->rowptr[i] && j < Q->colptr[col2+1]){
                    ++j;
                }
                //j out of range of col2
                if(j>=Q->colptr[col2+1]){
                    break;
                }

                if(Q->rowptr[j]==Q->rowptr[i]){
                    *dot+=Q->values[i]*Q->values[j];
                }
            }


        } else {
            MSG_ERROR("Unsupported Storage scheme: %s\n", mess_storage_t_str(Q->store_type));
            return MESS_ERROR_STORAGETYPE;
        }
    }
    return(0);
}

/**
 * @brief Compute the dot product between two columns of a complex matrix.
 * @param[in] Q     input matrix \f$Q\f$
 * @param[in] col1  input index \f$col_1\f$ of first column
 * @param[in] col2  input index \f$col_2\f$ of second column
 * @param[out] dot  output dot product
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_coldotc function computes the dot product between two columns of a complex matrix:
 * \f[dot = Q(:,col_1)^HQ(:,col_2).\f]
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_coldot
 * @sa mess_matrix_coldotE
 * @sa mess_matrix_coldotcE
 * @sa mess_matrix_colvecdot
 * @sa mess_matrix_colvecdotc
 */
int mess_matrix_coldotc(mess_matrix Q, mess_int_t col1, mess_int_t col2, mess_double_cpx_t *dot){
    MSG_FNAME(__func__);
    mess_int_t one =1;
    mess_check_nullpointer(Q);
    mess_check_complex(Q);

    if ( col1 <0 || col1 >= Q->cols) {
        MSG_ERROR("The column index 1 ( " MESS_PRINTF_INT " ) is out of range. \n", col1);
        return(MESS_ERROR_ARGUMENTS);
    }
    if ( col2 <0 || col2 >= Q->cols) {
        MSG_ERROR("The column index 2 ( " MESS_PRINTF_INT " ) is out of range. \n", col2);
        return(MESS_ERROR_ARGUMENTS);
    }

    if (MESS_IS_DENSE(Q)) {
#ifdef ZDOTC_MKL
        F77_GLOBAL(zdotc,ZDOTC)(dot, &(Q->rows), &(Q->values_cpx[col1*Q->ld]), &one, &(Q->values_cpx[col2*Q->ld]), &one);
#else
        *dot = F77_GLOBAL(zdotc,ZDOTC)(&(Q->rows), &(Q->values_cpx[col1*Q->ld]), &one, &(Q->values_cpx[col2*Q->ld]), &one);
#endif
    } else {
        mess_int_t i,j;
        if ( MESS_IS_CSR(Q)){
            mess_int_t setboth;
            mess_double_cpx_t val1, val2;
            //initial value
            *dot=0;
            for ( i = 0 ; i < Q->rows;i++ ) {
                val1=0; val2=0;
                setboth=-1;
                for ( j = Q->rowptr[i]; j < Q->rowptr[i+1]; j++  ) {
                    if ( Q->colptr[j] == col1){
                        val1 = Q->values_cpx[j];
                        ++setboth;
                        if(setboth){
                            *dot+=conj(val1)*val2;
                            break;
                        }
                    }
                    if( Q->colptr[j] == col2){
                        val2 = Q->values_cpx[j];
                        ++setboth;
                        if(setboth){
                            *dot+=conj(val1)*val2;
                            break;
                        }
                    }
                }
            }

        } else if ( MESS_IS_CSC(Q)){

            //inital value
            *dot=0;
            //j indices for col2
            j=Q->colptr[col2];
            //i indices for col1
            for(i=Q->colptr[col1];i<Q->colptr[col1+1];++i){
                while(Q->rowptr[j]<Q->rowptr[i] && j < Q->colptr[col2+1]){
                    ++j;
                }
                //j out of range of col2
                if(j>=Q->colptr[col2+1]){
                    break;
                }

                if(Q->rowptr[j]==Q->rowptr[i]){
                    *dot+=conj(Q->values_cpx[i])*Q->values_cpx[j];
                }
            }

        } else {
            MSG_ERROR("Unsupported Storage scheme: %s\n", mess_storage_t_str(Q->store_type));
            return MESS_ERROR_STORAGETYPE;
        }
    }
    return(0);
}

/**
 * @brief Update a column of a matrix by \f$\alpha\f$ times another column.
 * @param[in,out] Q     input/output matrix \f$Q\f$
 * @param[in] alpha     input coefficient \f$ \alpha \f$
 * @param[in] colc      input index \f$col_c\f$ of the column to add
 * @param[in] col1      input index \f$col_1\f$ of the column to update
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_colaxpy2 function updates a column \f$ col_1 \f$ of \f$ Q \f$ with an axpy-operation with a second
 * column \f$ col_c \f$.
 * It computes \f[ Q(:,col_1)\leftarrow Q(:,col_1)+\alpha\cdot Q(:,col_c).\f]
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_colaxpy
 * @sa mess_matrix_colvecaxpy
 *
 */
int mess_matrix_colaxpy2(mess_matrix Q, mess_double_cpx_t alpha, mess_int_t colc, mess_int_t col1){
    MSG_FNAME(__func__);
    mess_int_t one = 1;
    mess_check_nullpointer(Q);
    mess_check_real_or_complex(Q);
    if ( col1 <0 || col1 >= Q->cols) {
        MSG_ERROR("The column index 1 ( " MESS_PRINTF_INT " ) is out of range. \n", col1);
        return (MESS_ERROR_ARGUMENTS);
    }
    if ( colc <0 || colc >= Q->cols) {
        MSG_ERROR("The column index c ( " MESS_PRINTF_INT " ) is out of range. \n", colc);
        return (MESS_ERROR_ARGUMENTS);
    }

    if(MESS_IS_DENSE(Q)){
        if (MESS_IS_REAL(Q)){
            mess_int_t i;
            double  realalpha= creal(alpha);
            if(cimag(alpha)){
                MSG_WARN("Complex scalar for a real matrix, the imaginary information will be lost.\n");
            }
            if ( colc == col1 ) {
                for (i = 0; i < Q->rows; i++) {
                    Q->values[col1*Q->ld+i]  += realalpha * Q->values[colc*Q->ld+i];
                }
            } else {
                F77_GLOBAL(daxpy,DAXPY)(&(Q->rows),  &realalpha, &(Q->values[colc*Q->ld]), &one, &(Q->values[col1*Q->ld]), &one);
            }
        }
        else if (MESS_IS_COMPLEX(Q)){
            mess_int_t i;
            if (colc == col1) {
                for (i = 0; i < Q->rows; i++) {
                    Q->values_cpx[col1*Q->ld+i]  += alpha * Q->values_cpx[colc*Q->ld+i];
                }
            } else {
                F77_GLOBAL(zaxpy,ZAXPY)(&(Q->rows), &alpha, &(Q->values_cpx[colc*Q->ld]), &one, &(Q->values_cpx[col1*Q->ld]), &one);
            }
        }
    }else if(MESS_IS_CSC(Q)){
        if(MESS_IS_REAL(Q)){
            double  realalpha= creal(alpha);
            if(cimag(alpha)){
                MSG_WARN("Complex scalar for a real matrix, the imaginary information will be lost.\n");
            }

            double val, * newvalptr;
            mess_int_t i,j,max, k=Q->colptr[col1], endcolc=Q->colptr[colc+1], endcol1=Q->colptr[col1+1], colnnz=0, *newrowptr, *newcolptr;

            //if realalpha==0 nothing will be changed
            if(realalpha==0){
                return 0;
            }

            //alloc memory for the worst case
            //nnz elements before col1 + full column + elements after col1
            max= (Q->colptr[col1]-Q->colptr[0]) + Q->cols + (Q->colptr[Q->cols]-Q->colptr[col1+1]);
            mess_try_alloc(newvalptr,double* , sizeof(double)*max);
            mess_try_alloc(newrowptr,mess_int_t*,sizeof(mess_int_t)*max);
            mess_try_alloc(newcolptr,mess_int_t*,sizeof(mess_int_t)*(Q->cols+1));

            //copy values before row to set
            for(i=0;i<=col1;++i){newcolptr[i]=Q->colptr[i];}
            for(i=0;i<Q->colptr[col1];++i){newrowptr[i]=Q->rowptr[i]; newvalptr[i]=Q->values[i];}

            //compute new col
            for(i=Q->colptr[col1], j=Q->colptr[colc];i<endcol1 && j<endcolc;++i){

                if(Q->rowptr[i]==Q->rowptr[j]){
                    val = Q->values[i]+realalpha*Q->values[j];
                    if(val){
                        newvalptr[k]=val;
                        newrowptr[k]=Q->rowptr[i];
                        ++colnnz;
                        ++k;
                    }
                    ++j;
                }else{
                    if (Q->rowptr[i]<Q->rowptr[j]){
                        newvalptr[k]=Q->values[i];
                        newrowptr[k]=Q->rowptr[i];
                        ++k;
                        ++colnnz;
                    }else{
                        if(Q->rowptr[j]<Q->rowptr[i]){
                            val = realalpha*Q->values[j];
                            if(val){
                                newvalptr[k]=val;
                                newrowptr[k]=Q->rowptr[j];
                                ++colnnz;
                                ++k;
                            }
                            ++j;
                            --i; //dont change position of i in the next loop
                        }
                    }
                }
            }

            //compute new elements after endcol1 or endcolc indices
            for(;j<endcolc;++j){
                val = realalpha*Q->values[j];
                if(val){
                    newvalptr[k]=val;
                    newrowptr[k]=Q->rowptr[j];
                    ++k;
                    ++colnnz;
                }
            }

            for(;i<endcol1;++i){
                newvalptr[k]=Q->values[i];
                newrowptr[k]=Q->rowptr[i];
                ++k;
                ++colnnz;
            }

            //set new colptr, colnnz describes the change of nnz elements
            colnnz=colnnz-(Q->colptr[col1+1]-Q->colptr[col1]);

            //set new colptrs after col1
            for(i=col1+1;i<=Q->cols;++i){
                newcolptr[i]=Q->colptr[i]+colnnz;
            }

            //set new nnz
            Q->nnz+=colnnz;

            //copy elemnts after col1
            for(i=Q->colptr[col1+1];i<Q->colptr[Q->cols];++i){
                newvalptr[k]=Q->values[i];
                newrowptr[k]=Q->rowptr[i];
                ++k;
            }

            //resize arrays
            mess_try_realloc(newvalptr,double* ,sizeof(double)*Q->nnz);
            mess_try_realloc(newrowptr,mess_int_t*,sizeof(mess_int_t)*Q->nnz);

            //set new pointer
            mess_free(Q->colptr);
            mess_free(Q->rowptr);
            mess_free(Q->values);
            Q->rowptr=newrowptr;
            Q->colptr=newcolptr;
            Q->values=newvalptr;

        }else if(MESS_IS_COMPLEX(Q)){
            mess_double_cpx_t val, * newvalptr;
            mess_int_t i,j,max, k=Q->colptr[col1], endcolc=Q->colptr[colc+1], endcol1=Q->colptr[col1+1], colnnz=0, *newrowptr, *newcolptr;

            //if alpha==0 nothing will be changed
            if(alpha==0){
                return 0;
            }

            //alloc memory for the worst case
            //nnz elements before col1 + full column + elements after col1
            max= (Q->colptr[col1]-Q->colptr[0]) + Q->cols + (Q->colptr[Q->cols]-Q->colptr[col1+1]);
            mess_try_alloc(newvalptr,mess_double_cpx_t* , sizeof(mess_double_cpx_t)*max);
            mess_try_alloc(newrowptr,mess_int_t*,sizeof(mess_int_t)*max);
            mess_try_alloc(newcolptr,mess_int_t*,sizeof(mess_int_t)*(Q->cols+1));

            //copy values before row to set
            for(i=0;i<=col1;++i){newcolptr[i]=Q->colptr[i];}
            for(i=0;i<Q->colptr[col1];++i){newrowptr[i]=Q->rowptr[i]; newvalptr[i]=Q->values_cpx[i];}

            //compute new col
            for(i=Q->colptr[col1], j=Q->colptr[colc];i<endcol1 && j<endcolc;++i){

                if(Q->rowptr[i]==Q->rowptr[j]){
                    val = Q->values_cpx[i]+alpha*Q->values_cpx[j];
                    if(val){
                        newvalptr[k]=val;
                        newrowptr[k]=Q->rowptr[i];
                        ++colnnz;
                        ++k;
                    }
                    ++j;
                }else{
                    if (Q->rowptr[i]<Q->rowptr[j]){
                        newvalptr[k]=Q->values_cpx[i];
                        newrowptr[k]=Q->rowptr[i];
                        ++k;
                        ++colnnz;
                    }else{
                        if(Q->rowptr[j]<Q->rowptr[i]){
                            val = alpha*Q->values_cpx[j];
                            if(val){
                                newvalptr[k]=val;
                                newrowptr[k]=Q->rowptr[j];
                                ++colnnz;
                                ++k;
                            }
                            ++j;
                            --i; //dont change position of i in the next loop
                        }
                    }
                }
            }

            //compute new elements after endcol1 or endcolc indices
            for(;j<endcolc;++j){
                val = alpha*Q->values_cpx[j];
                if(val){
                    newvalptr[k]=val;
                    newrowptr[k]=Q->rowptr[j];
                    ++k;
                    ++colnnz;
                }
            }

            for(;i<endcol1;++i){
                newvalptr[k]=Q->values_cpx[i];
                newrowptr[k]=Q->rowptr[i];
                ++k;
                ++colnnz;
            }

            //set new colptr, colnnz describes the change of nnz elements
            colnnz=colnnz-(Q->colptr[col1+1]-Q->colptr[col1]);

            //set new colptrs after col1
            for(i=col1+1;i<=Q->cols;++i){
                newcolptr[i]=Q->colptr[i]+colnnz;
            }

            //set new nnz
            Q->nnz+=colnnz;

            //copy elemnts after col1
            for(i=Q->colptr[col1+1];i<Q->colptr[Q->cols];++i){
                newvalptr[k]=Q->values_cpx[i];
                newrowptr[k]=Q->rowptr[i];
                ++k;
            }

            //resize arrays
            mess_try_realloc(newvalptr,mess_double_cpx_t* ,sizeof(mess_double_cpx_t)*Q->nnz);
            mess_try_realloc(newrowptr,mess_int_t*,sizeof(mess_int_t)*Q->nnz);

            //set new pointer
            mess_free(Q->colptr);
            mess_free(Q->rowptr);
            mess_free(Q->values_cpx);
            Q->rowptr=newrowptr;
            Q->colptr=newcolptr;
            Q->values_cpx=newvalptr;


        }
    }else if(MESS_IS_CSR(Q)){
        if(MESS_IS_REAL(Q)){

            double  realalpha= creal(alpha);
            if(cimag(alpha)){
                MSG_WARN("Complex scalar for a real matrix, the imaginary information will be lost.\n");
            }

            double val, * newvalptr;
            mess_int_t i,j, k=0, l, max,*newrowptr, *newcolptr;

            //if realalpha==0 nothing will be changed
            if(realalpha==0){
                return 0;
            }

            //alloc memory for the worst case
            //nnz elements before col1 + full column + elements after col1
            max= (Q->rowptr[col1]-Q->rowptr[0]) + 2*Q->cols + (Q->rowptr[Q->cols]-Q->rowptr[col1+1]);
            mess_try_alloc(newvalptr,double* , sizeof(double)*max);
            mess_try_alloc(newcolptr,mess_int_t*,sizeof(mess_int_t)*max);
            mess_try_alloc(newrowptr,mess_int_t*,sizeof(mess_int_t)*(Q->rows+1));

            mess_int_t indcolc=-1;      //default indices of colc not found value
            mess_int_t indcol1=-1;      //default indices of col1 not found value
            mess_int_t nnzrow;          //number of nonzero elements in a row
            newrowptr[0]=0;

            //rowwise
            for(i=0;i<Q->rows;++i){
                val=0;
                nnzrow=0;
                indcol1=-1; indcolc=-1; //default for entry not set in that row
                //look for indices of colc and indices of col1 in that row
                for(l=Q->rowptr[i];l<Q->rowptr[i+1];++l){
                    if(Q->colptr[l]==colc){
                        indcolc=l;
                        if(indcol1!=-1){
                            break;
                        }
                    }else if(Q->colptr[l]==col1){
                        indcol1=l;
                        if(indcolc!=-1){
                            break;
                        }
                    }
                }

                for(j=Q->rowptr[i];j<Q->rowptr[i+1];++j){
                    if(Q->colptr[j]<col1){
                        newvalptr[k]=Q->values[j];
                        newcolptr[k]=Q->colptr[j];
                        ++k; ++nnzrow;
                    }else if(Q->colptr[j]==col1){
                        val=Q->values[j];
                        if(indcolc!=-1){
                            val+=realalpha*Q->values[indcolc]; //compute val<-Q(i,col1)+realalpha*Q(i,colc) if Q(i,colc) exists otherwise take val<-Q(i,col1)
                        }
                        if(val){ //update of col1
                            newvalptr[k]=val;
                            newcolptr[k]=col1;
                            ++k;val=0; ++nnzrow;
                        }
                        indcol1=-1; indcolc=-1;
                    }else{
                        //Q->colptr[j]>col1 and indcol1=-1 must hold, because Q(i,col1) does not exist
                        if(indcolc!=-1){
                            val = realalpha*Q->values[indcolc]; //compute val<-realalpha*Q(i,colc) if en
                            indcolc=-1;
                            if(val){ // add realalpha*Q(i,colc) entry in Q(i,col1), the entry was not there before
                                newvalptr[k]=val;
                                newcolptr[k]=col1;
                                ++k;val=0; ++nnzrow;
                            }
                        }
                        //copy values after col1
                        newvalptr[k]=Q->values[j];
                        newcolptr[k]=Q->colptr[j];
                        ++k; ++nnzrow;
                    }
                }
                if(indcolc!=-1){ //thats the case if colc<col1 and Q(i,col1) does not exists
                    val=realalpha*Q->values[indcolc];
                    indcolc=-1;
                    if(val){
                        newcolptr[k]=col1;
                        newvalptr[k]=val;
                        ++k; ++nnzrow; val=0;
                    }
                }
                newrowptr[i+1]=nnzrow+newrowptr[i];
            }

            Q->nnz=k;

            //resize arrays
            mess_try_realloc(newvalptr,double* ,sizeof(double)*Q->nnz);
            mess_try_realloc(newcolptr,mess_int_t*,sizeof(mess_int_t)*Q->nnz);

            //set new pointer
            mess_free(Q->colptr);
            mess_free(Q->rowptr);
            mess_free(Q->values);
            Q->rowptr=newrowptr;
            Q->colptr=newcolptr;
            Q->values=newvalptr;

        }else if (MESS_IS_COMPLEX(Q)){


            mess_double_cpx_t val, * newvalptr;
            mess_int_t i,j, k=0, l, max, *newrowptr, *newcolptr;

            //if alpha==0 nothing will be changed
            if(alpha==0){
                return 0;
            }

            //alloc memory for the worst case
            //nnz elements before col1 + full column + elements after col1
            max= (Q->rowptr[col1]-Q->rowptr[0]) + 2*Q->cols + (Q->rowptr[Q->cols]-Q->rowptr[col1+1]);
            mess_try_alloc(newvalptr,mess_double_cpx_t * , sizeof(mess_double_cpx_t)*max);
            mess_try_alloc(newcolptr,mess_int_t*,sizeof(mess_int_t)*max);
            mess_try_alloc(newrowptr,mess_int_t*,sizeof(mess_int_t)*(Q->rows+1));

            mess_int_t indcolc=-1;      //default indices of colc not found value
            mess_int_t indcol1=-1;      //default indices of col1 not found value
            mess_int_t nnzrow;          //number of nonzero elements in a row
            newrowptr[0]=0;

            //rowwise
            for(i=0;i<Q->rows;++i){
                val=0;
                nnzrow=0;
                indcol1=-1; indcolc=-1; //default for entry not set in that row
                //look for indices of colc and indices of col1 in that row
                for(l=Q->rowptr[i];l<Q->rowptr[i+1];++l){
                    if(Q->colptr[l]==colc){
                        indcolc=l;
                        if(indcol1!=-1){
                            break;
                        }
                    }else if(Q->colptr[l]==col1){
                        indcol1=l;
                        if(indcolc!=-1){
                            break;
                        }
                    }
                }

                for(j=Q->rowptr[i];j<Q->rowptr[i+1];++j){
                    if(Q->colptr[j]<col1){
                        newvalptr[k]=Q->values_cpx[j];
                        newcolptr[k]=Q->colptr[j];
                        ++k; ++nnzrow;
                    }else if(Q->colptr[j]==col1){
                        val=Q->values_cpx[j];
                        if(indcolc!=-1){
                            val+=alpha*Q->values_cpx[indcolc]; //compute val<-Q(i,col1)+realalpha*Q(i,colc) if Q(i,colc) exists otherwise take val<-Q(i,col1)
                        }
                        if(val){ //update of col1
                            newvalptr[k]=val;
                            newcolptr[k]=col1;
                            ++k;val=0; ++nnzrow;
                        }
                        indcol1=-1; indcolc=-1;
                    }else{
                        //Q->colptr[j]>col1 and indcol1=-1 must hold, because Q(i,col1) does not exist
                        if(indcolc!=-1){
                            val = alpha*Q->values_cpx[indcolc]; //compute val<-realalpha*Q(i,colc) if en
                            indcolc=-1;
                            if(val){ // add realalpha*Q(i,colc) entry in Q(i,col1), the entry was not there before
                                newvalptr[k]=val;
                                newcolptr[k]=col1;
                                ++k;val=0; ++nnzrow;
                            }
                        }
                        //copy values after col1
                        newvalptr[k]=Q->values_cpx[j];
                        newcolptr[k]=Q->colptr[j];
                        ++k; ++nnzrow;
                    }
                }
                if(indcolc!=-1){ //thats the case if colc<col1 and Q(i,col1) does not exists
                    val=alpha*Q->values_cpx[indcolc];
                    indcolc=-1;
                    if(val){
                        newcolptr[k]=col1;
                        newvalptr[k]=val;
                        ++k; ++nnzrow; val=0;
                    }
                }
                newrowptr[i+1]=nnzrow+newrowptr[i];
            }

            Q->nnz=k;

            //resize arrays
            mess_try_realloc(newvalptr,mess_double_cpx_t* ,sizeof(mess_double_cpx_t)*Q->nnz);
            mess_try_realloc(newcolptr,mess_int_t*,sizeof(mess_int_t)*Q->nnz);

            //set new pointer
            mess_free(Q->colptr);
            mess_free(Q->rowptr);
            mess_free(Q->values_cpx);
            Q->rowptr=newrowptr;
            Q->colptr=newcolptr;
            Q->values_cpx=newvalptr;

        }
    }else{
        MSG_ERROR("Unsupported Storage scheme: %s\n", mess_storage_t_str(Q->store_type));
        return MESS_ERROR_STORAGETYPE;
    }
    return (0);
}


/**
 * @brief Compute the dot product between a matrix column and a vector (real).
 * @param[in] Q     input the matrix \f$Q\f$
 * @param[in] col   input index \f$col\f$ of the column
 * @param[in] v     input vector \f$v\f$
 * @param[out] dot  output dot product
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_colvecdot function computes
 * \f[dot=Q(:,col)^Tv\f]
 * for real data.
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_coldot
 * @sa mess_matrix_coldotc
 * @sa mess_matrix_coldotE
 * @sa mess_matrix_coldotcE
 * @sa mess_matrix_colvecdotc
 */
int mess_matrix_colvecdot(mess_matrix Q, mess_int_t col, mess_vector v, double *dot ){
    MSG_FNAME(__func__);
    mess_check_nullpointer(Q);
    mess_check_nullpointer(v);
    mess_check_real(Q);
    mess_check_real(v);

    if ( col <0 || col >= Q->cols) {
        MSG_ERROR("The column index ( " MESS_PRINTF_INT " ) is out of range. \n", col);
        return(MESS_ERROR_ARGUMENTS);
    }
    if ( v->dim != Q->rows) {
        MSG_ERROR("The dimension of Q and v doesn't fit.\n");
        return (MESS_ERROR_DIMENSION);
    }
    if ( v->data_type != Q->data_type){
        MSG_ERROR("Q and v must have the same data type\n");
        return MESS_ERROR_DATATYPE;
    }
    if (MESS_IS_DENSE(Q)) {
#if defined(MESS_USE_OPENBLAS) && defined(__x86_64__)
        mess_int_t i;
        double dotx = 0 ;
        for ( i = 0; i<Q->rows; i++){
            dotx += Q->values[col*Q->ld+i]*v->values[i];
        }
        *dot = dotx;
#else
        mess_int_t one =1;
        *dot = ( F77_GLOBAL(ddot,DDOT)(&(Q->rows), &(Q->values[col*Q->ld]), &one, v->values, &one));
#endif
    } else if (MESS_IS_CSR(Q)){
        mess_int_t i,j;
        double d = 0;
        for ( i = 0 ; i < Q->rows;i++ ) {
            for ( j = Q->rowptr[i]; j < Q->rowptr[i+1]; j++  ) {
                if ( Q->colptr[j] == col ) {
                    d += Q->values[j] * v->values[i];
                }
            }
        }
        *dot = d;
    } else if (MESS_IS_CSC(Q)){
        mess_int_t i;
        double d = 0;
        for ( i = Q->colptr[col]; i < Q->colptr[col+1]; i++) {
            d+=Q->values[i]*v->values[Q->rowptr[i]];
        }
        *dot = d;
    } else {
        MSG_ERROR("Unsupported Storage scheme: %s\n", mess_storage_t_str(Q->store_type));
        return MESS_ERROR_STORAGETYPE;
    }
    return (0);
}

/**
 * @brief Computes the dot product between a matrix column and a vector (complex).
 * @param[in] Q     input the matrix \f$Q\f$
 * @param[in] col   input index \f$col\f$ of the column
 * @param[in] v     input vector \f$v\f$
 * @param[out] dot  output dot product
 *
 * The @ref mess_matrix_colvecdotc function computes
 * \f[dot=Q(:,col)^Hv\f]
 * for complex valued data.
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_coldot
 * @sa mess_matrix_coldotc
 * @sa mess_matrix_coldotE
 * @sa mess_matrix_coldotcE
 * @sa mess_matrix_colvecdot
 *
 */
int  mess_matrix_colvecdotc(mess_matrix Q, mess_int_t col, mess_vector v, mess_double_cpx_t *dot )
{
    MSG_FNAME(__func__);
    mess_check_nullpointer(Q);
    mess_check_nullpointer(v);

    if ( col <0 || col >= Q->cols) {
        MSG_ERROR("The column index ( " MESS_PRINTF_INT " ) is out of range. \n", col);
        return(MESS_ERROR_ARGUMENTS);
    }
    if ( v->dim != Q->rows) {
        MSG_ERROR("The dimension of Q and v doesn't fit.\n");
        return (MESS_ERROR_DIMENSION);
    }

    if (MESS_IS_DENSE(Q)) {
        if(MESS_IS_COMPLEX((Q))){
            if(MESS_IS_COMPLEX(v)){
#if defined(MESS_USE_OPENBLAS) && defined(__x86_64__)
                mess_int_t i;
                mess_double_cpx_t dotx = 0 ;
                for ( i = 0; i<Q->rows; i++){
                    dotx += conj(Q->values_cpx[col*Q->ld+i])*v->values_cpx[i];
                }
                *dot = dotx;
#else
                mess_int_t one =1;
#ifdef ZDOTC_MKL
                F77_GLOBAL(zdotc,ZDOTC)(dot, &(Q->rows), &(Q->values_cpx[col*Q->ld]), &one, v->values_cpx, &one);
#else
                *dot = F77_GLOBAL(zdotc,ZDOTC)(&(Q->rows), &(Q->values_cpx[col*Q->ld]), &one, v->values_cpx, &one);
#endif
#endif
            }else if(MESS_IS_REAL(v)){
                mess_int_t i;
                mess_double_cpx_t dotx = 0 ;
                for ( i = 0; i<Q->rows; i++){
                    dotx += Q->values_cpx[col*Q->ld+i]*v->values[i];
                }
                *dot = dotx;
            }else{
                MSG_ERROR("Unknown datatype for v!\n");
                return MESS_ERROR_DATATYPE;
            }
        }else if(MESS_IS_REAL(Q)){
            if(MESS_IS_COMPLEX(v)){
                mess_int_t i;
                mess_double_cpx_t dotx = 0 ;
                for ( i = 0; i<Q->rows; i++){
                    dotx += Q->values[col*Q->ld+i]*v->values_cpx[i];
                }
                *dot = dotx;
            }else if(MESS_IS_REAL(v)){//case real Q real v
                MSG_ERROR("Q and v are real. Please use mess_matrix_colvecdot.");
                return MESS_ERROR_DATA;
            }else{
                MSG_ERROR("Unknown dataype for v!\n");
                return MESS_ERROR_DATATYPE;
            }
        }else{
            MSG_ERROR("Unknown datatype for Q!\n");
            return MESS_ERROR_DATATYPE;

        }
    } else if (MESS_IS_CSR(Q)){
        if(MESS_IS_COMPLEX(Q)){
            if(MESS_IS_COMPLEX(v)){
                mess_int_t i,j;
                mess_double_cpx_t d = 0;
                for ( i = 0 ; i < Q->rows;i++ ) {
                    for ( j = Q->rowptr[i]; j < Q->rowptr[i+1]; j++  ) {
                        if ( Q->colptr[j] == col ) {
                            d += conj( Q->values_cpx[j]) * v->values_cpx[i];
                        }
                    }
                }
                *dot = d;
            }else if(MESS_IS_REAL(v)){
                mess_int_t i,j;
                mess_double_cpx_t d = 0;
                for ( i = 0 ; i < Q->rows;i++ ) {
                    for ( j = Q->rowptr[i]; j < Q->rowptr[i+1]; j++  ) {
                        if ( Q->colptr[j] == col ) {
                            d += conj( Q->values_cpx[j]) * v->values[i];
                        }
                    }
                }
                *dot = d;
            }else{
                MSG_ERROR("Unknown datatype for v!\n");
                return MESS_ERROR_DATATYPE;
            }
        }else if(MESS_IS_REAL(Q)){
            if(MESS_IS_COMPLEX(v)){
                mess_int_t i,j;
                mess_double_cpx_t d = 0;
                for ( i = 0 ; i < Q->rows;i++ ) {
                    for ( j = Q->rowptr[i]; j < Q->rowptr[i+1]; j++  ) {
                        if ( Q->colptr[j] == col ) {
                            d += Q->values[j] * v->values_cpx[i];
                        }
                    }
                }
                *dot = d;
            }else if(MESS_IS_REAL(v)){
                //case real Q real v
                MSG_ERROR("Q and v are real. Please use mess_matrix_colvecdot.");
                return MESS_ERROR_DATA;
            }else{
                MSG_ERROR("Unknown datatype for v!\n");
                return MESS_ERROR_DATATYPE;
            }
        }else{
            MSG_ERROR("Unknown datatype for Q!\n");
            return MESS_ERROR_DATATYPE;
        }
    } else if (MESS_IS_CSC(Q)){
        if(MESS_IS_COMPLEX(Q)){
            if(MESS_IS_COMPLEX(v)){
                mess_int_t i;
                mess_double_cpx_t d = 0;
                for ( i = Q->colptr[col]; i < Q->colptr[col+1]; i++) {
                    d+=conj(Q->values_cpx[i])*v->values_cpx[Q->rowptr[i]];
                }
                *dot = d;
            }else if (MESS_IS_REAL(v)){
                mess_int_t i;
                mess_double_cpx_t d = 0;
                for ( i = Q->colptr[col]; i < Q->colptr[col+1]; i++) {
                    d+=conj(Q->values_cpx[i])*v->values[Q->rowptr[i]];
                }
                *dot = d;
            }else{
                MSG_ERROR("Unknown datatype for v!");
                return MESS_ERROR_DATATYPE;
            }
        }else if(MESS_IS_REAL(Q)){
            if(MESS_IS_COMPLEX(v)){
                mess_int_t i;
                mess_double_cpx_t d = 0;
                for ( i = Q->colptr[col]; i < Q->colptr[col+1]; i++) {
                    d+=Q->values[i]*v->values_cpx[Q->rowptr[i]];
                }
                *dot = d;
            }else if(MESS_IS_REAL(v)){
                MSG_ERROR("Q and v are real. Please use mess_matrix_colvecdot.");
                return MESS_ERROR_DATA;
            }else{
                MSG_ERROR("Unknown datatype for Q!\n");
                return MESS_ERROR_DATATYPE;
            }
        }else{

        }
    } else {
        MSG_ERROR("Unsupported Storage scheme: %s\n", mess_storage_t_str(Q->store_type));
        return MESS_ERROR_STORAGETYPE;
    }
    return (0);
}


/**
 * @brief Update a vector using a column of a matrix.
 * @param[in] alpha     input coefficient \f$\alpha\f$
 * @param[in] col       input index \f$col\f$ of the column
 * @param[in] Q         input matrix \f$Q\f$
 * @param[in,out] v     input/output vector to update \f$v\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_colvecaxpy function computes
 * \f[v \leftarrow v+\alpha\cdot Q(:,col).\f]
 * The operation with exchanged roles is implemented in \ref mess_matrix_colaxpy. \n
 * In case of real \f$ Q  \f$ and real \f$ v \f$ the imaginary part of \f$ \alpha \f$ will be ignored.
 * If you want to avoid this effect convert \f$ v \f$ to complex vector before computation. \n
 * In case of complex \f$ Q \f$ and real \f$ v \f$, \f$ v \f$ will be a complex vector after computation.
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_colaxpy
 * @sa mess_matrix_colaxpy2
 *
 */
int mess_matrix_colvecaxpy(mess_double_cpx_t alpha, mess_int_t col,mess_matrix Q, mess_vector v){
    MSG_FNAME(__func__);
    mess_int_t one = 1;
    int ret = 0 ;
    double alpha_real;

    mess_check_nullpointer(Q);
    mess_check_nullpointer(v);
    mess_check_real_or_complex(Q);
    mess_check_real_or_complex(v);

    if ( col <0 || col >= Q->cols) {
        MSG_ERROR("The column index  ( " MESS_PRINTF_INT " ) is out of range. \n", col);
        return (MESS_ERROR_ARGUMENTS);
    }
    if ( v->dim != Q->rows) {
        MSG_ERROR("v and Q must have the same dimension.\n");
        return (MESS_ERROR_DIMENSION);
    }

    if(MESS_IS_DENSE(Q)){
        if (MESS_IS_REAL(Q) && MESS_IS_REAL(v)) {
            if ( cimag(alpha) != 0.0 ) {
                MSG_WARN("Real Matrix Q and real vector v. The imaginary part of alpha will be ignored!");
            }
            alpha_real = creal(alpha);
#if defined(MESS_USE_OPENBLAS) && defined(__x86_64__)
            mess_int_t i;
            for ( i = 0; i < Q->rows; i++){
                v->values[i] += alpha_real * Q->values[col*Q->ld + i];
            }
#else
            F77_GLOBAL(daxpy,DAXPY)(&(Q->rows), (double *) &alpha_real, &(Q->values[col*Q->ld]), &one, v->values, &one);
#endif
        }
        else if (MESS_IS_COMPLEX(Q) && MESS_IS_COMPLEX(v)){
            F77_GLOBAL(zaxpy,ZAXPY)(&(Q->rows), &alpha, &(Q->values_cpx[col*Q->ld]), &one, v->values_cpx, &one);
        } else if ( MESS_IS_COMPLEX(Q) && MESS_IS_REAL(v)){
            ret = mess_vector_tocomplex(v);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            for ( one = 0 ; one < v->dim; one ++) {
                v->values_cpx[one] += alpha*Q->values_cpx[Q->ld*col+one];
            }
        } else if (MESS_IS_REAL(Q) && MESS_IS_COMPLEX(v)) {
            for ( one = 0 ; one < v->dim; one ++) {
                v->values_cpx[one] += alpha*Q->values[Q->ld*col+one];
            }
        }
    }else if(MESS_IS_CSR(Q)){
        if(MESS_IS_REAL(Q) && MESS_IS_REAL(v)){
            mess_int_t i,j;
            for(i=0;i<Q->rows;++i){
                for(j=Q->rowptr[i];j<Q->rowptr[i+1];++j){
                    if(Q->colptr[j]==col){
                        v->values[i]+=creal(alpha)*Q->values[j];
                        break;
                    }
                }
            }
        }else if(MESS_IS_REAL(Q) && MESS_IS_COMPLEX(v)){
            mess_int_t i,j;
            for(i=0;i<Q->rows;++i){
                for(j=Q->rowptr[i];j<Q->rowptr[i+1];++j){
                    if(Q->colptr[j]==col){
                        v->values_cpx[i]+=alpha*Q->values[j];
                        break;
                    }
                }
            }
        }else if(MESS_IS_COMPLEX(Q) && MESS_IS_REAL(v)){
            mess_vector_tocomplex(v);
            mess_int_t i,j;
            for(i=0;i<Q->rows;++i){
                for(j=Q->rowptr[i];j<Q->rowptr[i+1];++j){
                    if(Q->colptr[j]==col){
                        v->values_cpx[i]+=alpha*Q->values_cpx[j];
                        break;
                    }
                }
            }
        }else{ //MESS_IS_COMPLEX(Q) && MESS_IS_COMPLEX(v) must hold
            mess_int_t i,j;
            for(i=0;i<Q->rows;++i){
                for(j=Q->rowptr[i];j<Q->rowptr[i+1];++j){
                    if(Q->colptr[j]==col){
                        v->values_cpx[i]+=alpha*Q->values_cpx[j];
                        break;
                    }
                }
            }
        }

    }else if(MESS_IS_CSC(Q)){
        if(MESS_IS_REAL(Q) && MESS_IS_REAL(v)){
            mess_int_t i;
            for(i=Q->colptr[col];i<Q->colptr[col+1];++i){
                v->values[Q->rowptr[i]]+=creal(alpha)*Q->values[i];
            }

        }else if(MESS_IS_REAL(Q) && MESS_IS_COMPLEX(v)){
            mess_int_t i;
            for(i=Q->colptr[col];i<Q->colptr[col+1];++i){
                v->values_cpx[Q->rowptr[i]]+=alpha*Q->values[i];
            }

        }else if(MESS_IS_COMPLEX(Q) && MESS_IS_REAL(v)){
            mess_vector_tocomplex(v);
            mess_int_t i;
            for(i=Q->colptr[col];i<Q->colptr[col+1];++i){
                v->values_cpx[Q->rowptr[i]]+=alpha*Q->values_cpx[i];
            }
        }else{ //MESS_IS_COMPLEX(Q) && MESS_IS_COMPLEX(v) must hold
            mess_int_t i;
            for(i=Q->colptr[col];i<Q->colptr[col+1];++i){
                v->values_cpx[Q->rowptr[i]]+=alpha*Q->values_cpx[i];
            }
        }

    }else {
        MSG_ERROR("The input matrix has an unknown datatype: Q=%s, v=%s.\n", mess_datatype_t_str(Q->data_type), mess_datatype_t_str(v->data_type));
        return(MESS_ERROR_DATATYPE);
    }
    return(0);
}

/**
 * @brief Perform an axpy-update on a column of a matrix.
 * @param[in] alpha     input coefficient \f$\alpha\f$
 * @param[in] v         input vector \f$v\f$ to add
 * @param[in] col       input index \f$col\f$ of the column
 * @param[in,out] Q     input/output matrix \f$Q\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_colaxpy function updates a column of matrix using an axpy-operation
 * \f[Q(:,col) \leftarrow \alpha\cdot v+Q(:,col)\f].
 *
 * @attention This function expects a dense matrix.
 *
 * @sa mess_matrix_colvecaxpy
 * @sa mess_matrix_colaxpy2
 *
 */
int  mess_matrix_colaxpy ( mess_double_cpx_t alpha, mess_vector v, mess_int_t col, mess_matrix Q )
{
    MSG_FNAME(__func__);
    mess_int_t i;
    mess_int_t one= 1 ;
    int ret = 0 ;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(Q);
    mess_check_nullpointer(v);
    mess_check_dense(Q);
    if (col < 0 || col >= Q->cols){
        MSG_ERROR("The column number is out of range.\n");
        return ( MESS_ERROR_DIMENSION);
    }
    if ( v->dim != Q->rows){
        MSG_ERROR("The input vector has the wrong dimension. Is:" MESS_PRINTF_INT " must: " MESS_PRINTF_INT "\n", v->dim, Q->rows);
        return (MESS_ERROR_DIMENSION);
    }


    /*-----------------------------------------------------------------------------
     *  do it
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_REAL(Q) && MESS_IS_REAL(v) ) {
        F77_GLOBAL(daxpy,DAXPY)(&(Q->rows), (double * ) &alpha, v->values, &one, &(Q->values[Q->ld*col]), &one);
    } else if ( MESS_IS_REAL(Q) && MESS_IS_COMPLEX(v)) {
        ret = mess_matrix_tocomplex(Q);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
#ifdef _OPENMP
#pragma omp parallel for private(i) default(shared)
#endif
        for (i=0; i < Q->rows; i++){
            Q->values_cpx[col*Q->ld+i] = Q->values_cpx[col*Q->ld+i] + alpha * v->values_cpx[i];
        }
    } else if ( MESS_IS_COMPLEX(Q) && MESS_IS_REAL(v)) {
#ifdef _OPENMP
#pragma omp parallel for private(i) default(shared)
#endif
        for (i=0; i < Q->rows; i++){
            Q->values_cpx[col*Q->ld+i] = Q->values_cpx[col*Q->ld+i] + alpha * v->values[i];
        }
    } else if ( MESS_IS_COMPLEX(Q) && MESS_IS_COMPLEX(v) ) {
        F77_GLOBAL(zaxpy,ZAXPY)(&(Q->rows), &alpha, v->values_cpx, &one, &(Q->values_cpx[Q->ld*col]), &one);
    } else {
        MSG_ERROR("Unknown / Unsupported data type combination\n");
        return(MESS_ERROR_DATATYPE);
    }


    return (0);
}   /* -----  end of function mess_matrix_colaxpy  ----- */

/**
 * @brief Compute the weighted dot product of two matrix columns (real).
 * @param[in] Q         input real matrix \f$Q\f$
 * @param[in] E         input symmetric positive definite weight matrix \f$E\f$
 * @param[in] col1      input index \f$col_1\f$ of the first column
 * @param[in] col2      input index \f$col_2\f$ of the second column
 * @param[out] dot      output dot product
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_coldotE function computes
 * \f[dot = Q(:,col_1)^T*E*Q(:,col_2)\f]
 * for real matrices.
 *
 * @attention \f$ E \f$ need to be symmetric positive definite.
 *
 * @sa mess_matrix_coldot
 * @sa mess_matrix_coldotc
 * @sa mess_matrix_coldotcE
 *
 */
int mess_matrix_coldotE(mess_matrix Q, mess_matrix E, mess_int_t col1, mess_int_t col2, double *dot ) {
    MSG_FNAME(__func__);
    mess_int_t one =1;
    mess_vector tx, ty ;
    int ret = 0;

    mess_check_nullpointer(Q);
    mess_check_nullpointer(E);
    mess_check_real(Q);
    mess_check_real(E);
    mess_check_square(E);
    mess_check_same_rows(Q,E); //check if E->cols==Q->rows <-> E->rows==Q->rows (E is square)

    if ( col1 <0 || col1 >= Q->cols) {
        MSG_ERROR("The column index 1 ( " MESS_PRINTF_INT " ) is out of range. \n", col1);
        return (MESS_ERROR_ARGUMENTS);
    }
    if ( col2 <0 || col2 >= Q->cols) {
        MSG_ERROR("The column index 2 ( " MESS_PRINTF_INT " ) is out of range. \n", col2);
        return (MESS_ERROR_ARGUMENTS);
    }

    if(MESS_IS_DENSE(Q) && MESS_IS_DENSE(E)){
        MESS_INIT_VECTORS(&ty,&tx);
        ret = mess_vector_alloc(ty, E->cols, MESS_REAL);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_alloc(tx, 1, MESS_REAL);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        mess_free(tx->values);
        tx->values = &(Q->values[col2*Q->ld]);
        tx->dim = Q ->rows;
        ret = mess_matrix_mvp(MESS_OP_NONE,E, tx, ty);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
        *dot = F77_GLOBAL(ddot,DDOT)(&(Q->rows), &(Q->values[col1*Q->ld]), &one, ty->values /* &(Q->values[col2*Q->rows])*/, &one);
        tx->values=NULL;
    }else{
        mess_vector temp;
        MESS_INIT_VECTORS(&temp,&ty,&tx);
        ret = mess_vector_alloc(temp,Q->rows,MESS_REAL);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_alloc(ty,Q->rows,MESS_REAL);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_alloc(tx,Q->rows,MESS_REAL);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
        ret = mess_matrix_getcol(Q,col1,temp);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_getcol);
        ret = mess_matrix_getcol(Q,col2,ty);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_getcol);
        ret = mess_matrix_mvp(MESS_OP_NONE,E,temp,tx);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mvp);
        ret = mess_vector_clear(&temp);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_clear);
        ret = mess_vector_dot(tx,ty,dot);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_dot);
    }
    mess_vector_clear(&tx);
    mess_vector_clear(&ty);
    return (0);
}


/**
 * @brief Compute the weighted dot product of two  matrix columns (complex).
 * @param[in] Q         input complex  matrix \f$Q\f$
 * @param[in] E         input symmetric positive definite weight matrix \f$E\f$
 * @param[in] col1      input index \f$col_1\f$ of the first column
 * @param[in] col2      input index \f$col_2\f$ of the second column
 * @param[out] dot      output dot product
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_coldotcE function computes
 * \f[dot = Q(:,col_1)^H*E*Q(:,col_2)\f]
 * for complex matrices.
 *
 * @attention \f$ E \f$ need to be symmetric positive definite.
 *
 * @sa mess_matrix_coldot
 * @sa mess_matrix_coldotc
 * @sa mess_matrix_coldotE
 *
 */
int mess_matrix_coldotcE(mess_matrix Q, mess_matrix E, mess_int_t col1, mess_int_t col2, mess_double_cpx_t *dot ) {
    MSG_FNAME(__func__);
    mess_int_t one =1;
    mess_vector tx, ty ;
    int ret = 0;

    mess_check_nullpointer(Q);
    mess_check_nullpointer(E);
    mess_check_complex(Q);
    mess_check_real(E);
    mess_check_square(E);
    mess_check_same_rows(Q,E); //check if E->cols==Q->rows <-> E->rows==Q->rows (E is square)

    if ( col1 <0 || col1 >= Q->cols) {
        MSG_ERROR("The column index 1 ( " MESS_PRINTF_INT " ) is out of range. \n", col1);
        return (MESS_ERROR_ARGUMENTS);
    }
    if ( col2 <0 || col2 >= Q->cols) {
        MSG_ERROR("The column index 2 ( " MESS_PRINTF_INT " ) is out of range. \n", col2);
        return (MESS_ERROR_ARGUMENTS);
    }

    if(MESS_IS_DENSE(Q) && MESS_IS_DENSE(E)){
        MESS_INIT_VECTORS(&ty,&tx);
        ret = mess_vector_alloc(ty, E->cols, MESS_COMPLEX);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(tx, Q->rows, Q->data_type);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_matrix_getcol(Q,col2,tx);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
        ret = mess_matrix_mvp(MESS_OP_NONE,E, tx, ty);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
#ifdef ZDOTC_MKL
        F77_GLOBAL(zdotc,ZDOTC)(dot, &(Q->rows),&(Q->values_cpx[col1*Q->ld]), &one, ty->values_cpx /* &(Q->values_cpx[col2*Q->rows])*/, &one);
#else
        *dot= F77_GLOBAL(zdotc,ZDOTC)(&(Q->rows),&(Q->values_cpx[col1*Q->ld]), &one, ty->values_cpx /* &(Q->values_cpx[col2*Q->rows])*/, &one);
#endif
    }else{
        mess_vector temp;
        MESS_INIT_VECTORS(&temp,&ty,&tx);
        ret = mess_vector_alloc(temp,Q->rows,MESS_COMPLEX);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_alloc(ty,Q->rows,MESS_COMPLEX);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_alloc(tx,Q->rows,MESS_COMPLEX);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
        ret = mess_matrix_getcol(Q,col2,temp);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_getcol);
        ret = mess_matrix_getcol(Q,col1,ty);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_getcol);
        ret = mess_matrix_mvp(MESS_OP_NONE,E,temp,tx);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mvp);
        ret = mess_vector_clear(&temp);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_clear);
        ret = mess_vector_dotc(ty,tx,dot);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_dotc);
    }
    MESS_CLEAR_VECTORS(&tx,&ty);
    return (0);
}


/**
 * @brief Compute the \f$E\f$ norm of a matrix column.
 * @param[in] Q         input matrix \f$Q\f$
 * @param[in] E         input symmetric positive definite weight matrix \f$E\f$
 * @param[in] col       input index \f$col\f$ of the column
 * @param[out] norm     output norm
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_colnormE function computes
 * \f[norm = \Vert Q(:,col)\Vert_E = \sqrt{ Q(:,col)^TEQ(:,col)}\f]
 * or
 * \f[norm = \Vert Q(:,col)\Vert_E = \sqrt{ Q(:,col)^HEQ(:,col)}\f]
 * for real or complex matrix \f$ Q \f$ and real matrix \f$ E \f$.
 *
 * @attention \f$ E \f$ needs to be symmetric positive definite.
 *
 * @sa mess_matrix_colnorm
 *
 */
int mess_matrix_colnormE(mess_matrix Q, mess_matrix E,mess_int_t col, double *norm){
    MSG_FNAME(__func__);
    int ret = 0 ;
    mess_check_nullpointer(Q);
    mess_check_dense(Q);
    mess_check_dense(E);
    mess_check_real_or_complex(Q);
    mess_check_real(E);

    if ( col <0 || col >= Q->cols) {
        MSG_ERROR("The column index 1 ( " MESS_PRINTF_INT " ) is out of range. \n", col);
        return (MESS_ERROR_ARGUMENTS);
    }

    if(MESS_IS_REAL(Q)){
        ret =  mess_matrix_coldotE(Q,E,col,col, norm);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldotE);
        *norm = sqrt(*norm);
    }else{
        mess_double_cpx_t nrmcpx;
        ret =  mess_matrix_coldotcE(Q,E,col,col, &nrmcpx);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldotcE);
        *norm = sqrt(creal(nrmcpx));
    }
    return(0);
}


