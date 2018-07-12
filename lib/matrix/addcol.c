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
 * @file lib/matrix/addcol.c
 * @brief Add new columns to a matrix.
 * @author @koehlerm
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"


/**
 * @brief Add columns to an existing matrix.
 * @param[in,out] matrix    input/output matrix where the new columns should be added
 * @param[in] toadd         input matrix with the columns to be added
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_addcols function adds columns to an existing matrix. \n
 * It works like the @verbatim Z = [Z V]; @endverbatim statement in @matlab. \n
 * @attention The function currently only supports dense matrices.
 *
 * @sa mess_matrix_mgs_add
 * @sa mess_matrix_addcols2p
 * @sa mess_matrix_addcols1
 */
int mess_matrix_addcols(mess_matrix matrix, mess_matrix toadd){
    MSG_FNAME(__func__);
    mess_int_t start, n;
    int conv = -1;
    int ret = 0;
    mess_matrix toadd_work = NULL;

    mess_check_nullpointer(matrix);
    mess_check_real_or_complex(matrix);
    mess_check_dense(matrix);

    if ( toadd == NULL ) { return (0);}
    if ( toadd->rows == 0 || toadd->cols == 0){  return 0; }
    if ( matrix->rows==0 && matrix->cols ==0) {
        ret = mess_matrix_convert(toadd, matrix, MESS_DENSE);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
        return (0);
    }
    if ( matrix->rows != toadd->rows) {
        MSG_ERROR("The row dimension of the matrices have to be the same.\n");
        return (MESS_ERROR_DIMENSION);
    }
    mess_check_dense(matrix);
    MESS_MATRIX_CHECKFORMAT(toadd, toadd_work, conv, MESS_DENSE);

    if ( MESS_IS_COMPLEX(matrix) || MESS_IS_COMPLEX(toadd)) {
        ret = mess_matrix_tocomplex(matrix);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        ret = mess_matrix_tocomplex(toadd_work);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
    }


    start = matrix->cols*matrix->ld;
    if ( MESS_IS_COMPLEX(matrix)) {
        n = matrix->cols+toadd_work->cols;
        mess_try_realloc(matrix->values_cpx, mess_double_cpx_t * ,sizeof ( mess_double_cpx_t) * n * matrix->ld);
        F77_GLOBAL(zlacpy,ZLACPY)("A",&toadd_work->rows, &toadd_work->cols, toadd_work->values_cpx, &toadd_work->ld,matrix->values_cpx+start,&matrix->ld);
        matrix->cols = n;
        matrix->nnz = n * matrix->rows;
    } else {
        n = matrix->cols+toadd_work->cols;
        mess_try_realloc(matrix->values, double * ,sizeof (double) * n * matrix->ld);
        F77_GLOBAL(dlacpy,DLACPY)("A",&toadd_work->rows, &toadd_work->cols, toadd_work->values, &toadd_work->ld,matrix->values+start,&matrix->ld);
        matrix->cols = n;
        matrix->nnz = n * matrix->rows;
    }
    if ( conv == 0) {
        mess_matrix_clear(&toadd_work);
    }


    return (0);
}

/**
 * @brief Add the sum of two column matrices to a real matrix.
 * @param[in,out] Z     input/output matrix where to add the sum
 * @param[in]    s1     input scaling factor for V1
 * @param[in]    V1     input matrix V1
 * @param[in]    s2     input scaling factor for V2
 * @param[in]    V2     input matrix V2
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_addcols2p function forms
 * @verbatim Z = [Z , s1*V1+s2*V2]  @endverbatim
 * @attention The function currently only supports real dense matrices.
 * @sa mess_matrix_addcols2p
 * @sa mess_matrix_addcols1
 */
int mess_matrix_addcols2p ( mess_matrix Z , double s1, mess_matrix V1, double s2, mess_matrix V2 )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t n,j,i, start;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(Z);
    mess_check_nullpointer(V1);
    mess_check_nullpointer(V2);
    mess_check_dense(Z);
    mess_check_dense(V1);
    mess_check_dense(V2);
    mess_check_real(Z);
    mess_check_real(V1);
    mess_check_real(V2);
    mess_check_same_size(V1,V2);

    if ( Z->rows == 0 || Z->cols == 0 ) {
        ret = mess_matrix_copy(V1, Z);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        ret = mess_matrix_add(s2,V2,s1,Z);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        return (0);
    }

    if ( Z ->rows != V1->rows){
        MSG_ERROR("The columns to add needs the same number of rows.\n");
        return (MESS_ERROR_DIMENSION);
    }

    start = Z->cols * Z->ld;
    n = Z->cols+V1->cols;
    mess_try_realloc(Z->values,double *,sizeof ( double) * n * Z->ld);
    Z->cols = n;
    Z->nnz = n * Z->rows;


#ifdef _OPENMP
#pragma omp parallel for private(j,i) default(shared)
#endif
    for ( j = 0; j < V1->cols; j++ ){
        for ( i = 0; i < V1->rows; i++){
            Z->values[start+j*Z->ld+i] = s1*V1->values[j*V1->ld+i] + s2 *V2->values[j*V2->ld+i];
        }
    }
    return (0);
}       /* -----  end of function mess_matrix_addcols2p  ----- */

/**
 * @brief Add a scaled column to a real matrix.
 * @param[in,out] Z     input/output matrix where to add the column
 * @param[in]    s1     input scaling factor  V1
 * @param[in]    V1     input matrix V1
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_addcols1 function forms
 * @verbatim  Z = [Z , s1*V1] @endverbatim
 * @attention The function currently only supports real dense matrices.
 * @sa mess_matrix_addcols2p
 * @sa mess_matrix_addcols
 *
 */
int mess_matrix_addcols1 ( mess_matrix Z , double s1, mess_matrix V1)
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t n,j,i, start;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(Z);
    mess_check_nullpointer(V1);
    mess_check_dense(Z);
    mess_check_dense(V1);
    mess_check_real(Z);
    mess_check_real(V1);

    if ( Z->rows == 0 || Z->cols == 0 ) {
        ret = mess_matrix_copy(V1, Z);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        ret = mess_matrix_scale(s1,Z);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
        return (0);
    }

    if ( Z ->rows != V1->rows){
        MSG_ERROR("The columns to add needs the same number of rows.\n");
        return (MESS_ERROR_DIMENSION);
    }

    start = Z->cols*Z->ld;
    n = Z->cols+V1->cols;
    mess_try_realloc(Z->values,double *,sizeof ( double) * n * Z->ld);
    Z->cols = n;
    Z->nnz = n * Z->rows;


#ifdef _OPENMP
#pragma omp parallel for private(j,i) default(shared)
#endif
    for ( j = 0; j < V1->cols; j++ ){
        for ( i = 0 ; i < V1->rows; i ++) {
            Z->values[start+j*Z->ld+i] = s1*V1->values[j*V1->ld+i] ;
        }
    }
    return (0);
}       /* -----  end of function mess_matrix_addcols1  ----- */

