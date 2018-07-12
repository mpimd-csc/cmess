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
 * @file lib/matrix/triul.c
 * @brief Upper triangular part of @ref matrix.
 * @author @mbehr
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 *
 * @brief Upper triangular part of @ref mess_matrix.
 * @param[in,out] mat   input/output @ref MESS_DENSE @ref mess_matrix.
 * @param[in] k         input scalar specifies the k-th diagonal of @p mat.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_triu function sets everything below the @p k-th
 * diagonal to zero.
 * Only @ref MESS_DENSE matrices are supported.
 * If @p k is zero everything below the main diagonal is set to zero.
 * if @p k is positive/negative everything below the @p k -th upper/lower main diagonal is set to zero.
 *
 */
int mess_matrix_triu(mess_matrix mat, mess_int_t k){
    mess_int_t i,j=0;
    int ret = 0;
    MSG_FNAME(__func__);


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(mat);
    mess_check_real_or_complex(mat);
    mess_check_dense(mat);

    /*-----------------------------------------------------------------------------
     *  check for early return
     *-----------------------------------------------------------------------------*/
    if (k <= -(mat->rows-1)){
        MSG_INFO("k="MESS_PRINTF_INT",  rows = "MESS_PRINTF_INT", %s will return original matrix\n", k, mat->rows, __func__)
        return 0;
    }

    if (mat->cols <= k){
        MSG_INFO("k="MESS_PRINTF_INT",  cols = "MESS_PRINTF_INT", %s will set all entries to 0\n", k, mat->cols, __func__)
        ret = mess_matrix_zeros(mat);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);
        return 0;
    }


    /*-----------------------------------------------------------------------------
     *  apply function to a matrix
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(mat)){
        for(j=0;j<mat->cols;++j){
            for(i=0;i<mat->rows;++i){
                if(j-i<k){
                    mat->values[i+j*mat->ld]=0;
                }
            }
        }
    }else{
        for(j=0;j<mat->cols;++j){
            for(i=0;i<mat->rows;++i){
                if(j-i<k){
                    mat->values_cpx[i+j*mat->ld]=0;
                }
            }
        }
    }

    return 0;
} /* -----  end of function mess_matrix_triu  ----- */


/**
 *
 * @brief Lower triangular part of @ref mess_matrix.
 * @param[in,out] mat   input/output @ref MESS_DENSE @ref mess_matrix.
 * @param[in] k         input scalar specifies the k-th diagonal of @p mat.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_tril function sets everything above the @p k-th
 * diagonal to zero.
 * Only @ref MESS_DENSE matrices are supported.
 * If @p k is zero everything above the main diagonal is set to zero.
 * if @p k is positive/negative everything above the @p k -th upper/lower main diagonal is set to zero.
 *
 */
int mess_matrix_tril(mess_matrix mat, mess_int_t k){
    mess_int_t i,j=0;
    int ret = 0;
    MSG_FNAME(__func__);


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(mat);
    mess_check_real_or_complex(mat);
    mess_check_dense(mat);

    /*-----------------------------------------------------------------------------
     *  check for early return
     *-----------------------------------------------------------------------------*/
    if ( mat->cols <= k +1 ){
        MSG_INFO("k="MESS_PRINTF_INT",  cols = "MESS_PRINTF_INT", %s will return original matrix\n", k, mat->cols, __func__);
        return 0;
    }

    if (k <= -(mat->rows)){
        MSG_INFO("k="MESS_PRINTF_INT",  rows = "MESS_PRINTF_INT", %s will set all entries to 0\n", k, mat->rows, __func__);
        ret = mess_matrix_zeros(mat);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);
        return 0;
    }


    /*-----------------------------------------------------------------------------
     *  apply function to a matrix
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(mat)){
        for(j=0;j<mat->cols;++j){
            for(i=0;i<mat->rows;++i){
                if(j-i>k){
                    mat->values[i+j*mat->ld]=0;
                }
            }
        }
    }else{
        for(j=0;j<mat->cols;++j){
            for(i=0;i<mat->rows;++i){
                if(j-i>k){
                    mat->values_cpx[i+j*mat->ld]=0;
                }
            }
        }
    }

    return 0;
} /* -----  end of function mess_matrix_tril  ----- */

