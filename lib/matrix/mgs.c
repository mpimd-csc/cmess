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
 * @file lib/matrix/mgs.c
 * @brief The modified Gram-Schmidt algorithm.
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
#include <math.h>



/**
 * @brief Perform a modified Gram-Schmidt orthogonalization.
 * @param[in] A     input matrix \f$A\f$
 * @param[out] Q    output orthogonal matrix \f$Q\f$
 * @param[out] R    output upper triangular coefficient matrix \f$R\f$ (set to @c NULL if not wanted)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_mgs function computes an orthonormal basis for \f$ A \f$ using
 * the modified Gram-Schmidt process. \n
 * If \f$ R \f$ is not equal to @c NULL the coefficients
 * are computed, too. \n
 * Sparse matrices are converted to dense ones immediately.
 *
 * @sa mess_matrix_mgs_inplace
 * @sa mess_matrix_mgs_add
 * @sa mess_matrix_orth
 *
 */
int  mess_matrix_mgs ( mess_matrix A, mess_matrix Q, mess_matrix R )
{
    MSG_FNAME(__func__);
    int wantR;
    int ret = 0;
    mess_int_t k,j;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(Q);

    wantR = 0;
    if ( R != NULL){
        wantR = 1;
        MESS_MATRIX_RESET(R);
        ret = mess_matrix_alloc(R,A->cols,A->cols,A->cols*A->cols, MESS_DENSE, A->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    }

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_convert(A,Q,MESS_DENSE); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    if(MESS_IS_REAL(A)){
        double nrmR;
        if ( wantR == 0){
            for (k = 0; k < Q->cols; k++){
                ret =  mess_matrix_colnorm(Q,k, &nrmR);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);
                ret =  mess_matrix_colscale(Q,k,1.0/nrmR);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colscale);
                for (j=k+1;j<Q->cols;j++){
                    ret = mess_matrix_coldot(Q,k,j, &nrmR);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldot);
                    ret = mess_matrix_colaxpy2(Q,-nrmR,k,j);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
                }
            }
        } else {
            for (k = 0; k < Q->cols; k++){
                ret  = mess_matrix_colnorm(Q,k, &nrmR );            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);
                R->values[k+k*R->ld] = nrmR;
                ret = mess_matrix_colscale(Q,k,1.0/nrmR);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colscale);
                for (j=k+1;j<Q->cols;j++){
                    ret = mess_matrix_coldot(Q,k,j, &nrmR);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldot);
                    R->values[k+j*R->ld] = nrmR;
                    ret = mess_matrix_colaxpy2(Q,-nrmR,k,j);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
                }
            }
        }
    }else if(MESS_IS_COMPLEX(A)){
        mess_double_cpx_t nrmR;
        if ( wantR == 0){
            for (k = 0; k < Q->cols; k++){
                ret =  mess_matrix_coldotc(Q,k,k, &nrmR);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);
                nrmR = sqrt(creal(nrmR));
                ret =  mess_matrix_colscale(Q,k,1.0/nrmR);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colscale);
                for (j=k+1;j<Q->cols;j++){
                    ret = mess_matrix_coldotc(Q,k,j, &nrmR);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldot);
                    ret = mess_matrix_colaxpy2(Q,-nrmR,k,j);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
                }
            }
        } else {
            for (k = 0; k < Q->cols; k++){
                ret  = mess_matrix_coldotc(Q,k,k, &nrmR );          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);
                nrmR = sqrt(creal(nrmR));
                R->values_cpx[k+k*R->ld] = nrmR;
                ret = mess_matrix_colscale(Q,k,1.0/nrmR);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colscale);
                for (j=k+1;j<Q->cols;j++){
                    ret = mess_matrix_coldotc(Q,k,j, &nrmR);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldot);
                    R->values_cpx[k+j*R->ld] = nrmR;
                    ret = mess_matrix_colaxpy2(Q,-nrmR,k,j);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
                }
            }
        }
    }else{
        MSG_ERROR("Unsupported datatype for A!");
        return MESS_ERROR_DATATYPE;
    }
    return (0);
}       /* -----  end of function mess_matrix_mgs  ----- */


/**
 * @brief Perform a modified Gram-Schmidt orthogonalization inplace on a matrix.
 * @param[in,out] Q on input: matrix to orthonormalize \n
 *                  on output: orthonormal matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_mgs_inplace function computes an orthonormal basis for \f$ Q \f$ using
 * the modified Gram-Schmidt process. All operations are inplace on \f$ Q \f$ and need no extra
 * memory.
 * @sa mess_matrix_mgs
 * @sa mess_matrix_mgs_add
 * @sa mess_matrix_orth
 *
 */
int  mess_matrix_mgs_inplace ( mess_matrix Q )
{
    MSG_FNAME(__func__);
    mess_int_t k,j;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(Q);
    mess_check_dense(Q);

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/

    if(MESS_IS_REAL(Q)){
        double nrmR=0;
        for (k = 0; k < Q->cols; k++){
            ret = mess_matrix_colnorm(Q,k, &nrmR);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);
            ret = mess_matrix_colscale(Q,k,1.0/nrmR);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colscale);
            for (j=k+1;j<Q->cols;j++){
                ret = mess_matrix_coldot(Q,k,j, &nrmR);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldot);
                ret = mess_matrix_colaxpy2(Q,-nrmR,k,j);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
            }
        }
    }else if(MESS_IS_COMPLEX(Q)){
        mess_double_cpx_t nrmR=0;
        for (k = 0; k < Q->cols; k++){
            ret =  mess_matrix_coldotc(Q,k,k, &nrmR);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);
            nrmR = sqrt(creal(nrmR));
            ret = mess_matrix_colscale(Q,k,1.0/nrmR);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colscale);
            for (j=k+1;j<Q->cols;j++){
                ret = mess_matrix_coldotc(Q,k,j, &nrmR);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldot);
                ret = mess_matrix_colaxpy2(Q,-nrmR,k,j);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
            }
        }
    }else{
        MSG_ERROR("Unsupported datatype for Q!");
        return MESS_ERROR_DATATYPE;
    }

    return(0);
}       /* -----  end of function mess_matrix_mgs  ----- */

