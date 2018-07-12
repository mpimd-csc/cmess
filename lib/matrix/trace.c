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
 * @file lib/matrix/trace.c
 * @brief Compute the trace of a matrix.
 * @author @koehlerm
 * @author @dykstra
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/eigenvalue.h"
#include "mess/error_macro.h"
#include <complex.h>

/**
 * @brief Compute the trace of a real matrix.
 * @param[in] A input matrix \f$A\f$
 * @param[out] tr   output trace of \f$A\f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_trace function computes the trace of a real matrix.
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_tracec
 *
 */
int  mess_matrix_trace ( mess_matrix A, double *tr )
{
    MSG_FNAME(__func__);
    double t = 0;
    mess_int_t i=0,j=0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(tr);
    mess_check_square(A);
    mess_check_real(A);

    if ( MESS_IS_DENSE(A)){
        for (i=0; i < A->rows; i++){
            t += A->values[A->ld*i+i];
        }
    } else if ( MESS_IS_CSR(A)){
        for (i=0; i<A->rows;i++){
            for ( j=A->rowptr[i];j<A->rowptr[i+1];j++){
                if (A->colptr[j] == i){
                    t += A->values[j];
                    break;
                }
            }
        }
    }  else if ( MESS_IS_CSC(A)){
        for (i=0; i<A->cols;i++){
            for ( j=A->colptr[i];j<A->colptr[i+1];j++){
                if (A->rowptr[j] == i){
                    t += A->values[j];
                    break;
                }
            }
        }
    } else {
        MSG_ERROR("storage type isn't supportedi: %s\n", mess_storage_t_str(A->store_type));
        return (MESS_ERROR_STORAGETYPE);
    }

    *tr = t;
    return(0);
}       /* -----  end of function mess_matrix_trace  ----- */

/**
 * @brief Compute the trace of a complex matrix.
 * @param[in] A input matrix \f$A\f$
 * @param[out] tr   output trace \f$A\f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_tracec function computes the trace of a complex matrix.
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_trace
 */
int  mess_matrix_tracec ( mess_matrix A, mess_double_cpx_t *tr )
{
    MSG_FNAME(__func__);
    mess_double_cpx_t t = 0;
    mess_int_t i=0,j=0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(tr);
    mess_check_complex(A);
    mess_check_square(A);

    if ( MESS_IS_DENSE(A)){
        for (i=0; i < A->rows; i++){
            t += A->values_cpx[A->ld*i+i];
        }
    } else if ( MESS_IS_CSR(A)){
        for (i=0; i<A->rows;i++){
            for ( j=A->rowptr[i];j<A->rowptr[i+1];j++){
                if (A->colptr[j] == i){
                    t += A->values_cpx[j];
                    break;
                }
            }
        }
    }  else if ( MESS_IS_CSC(A)){
        for (i=0; i<A->cols;i++){
            for ( j=A->colptr[i];j<A->colptr[i+1];j++){
                if (A->rowptr[j] == i){
                    t += A->values_cpx[j];
                    break;
                }
            }
        }
    } else {
        MSG_ERROR("storage type isn't supported: %s\n", mess_storage_t_str(A->store_type));
        return(MESS_ERROR_STORAGETYPE);
    }

    *tr = t;
    return(0);
}       /* -----  end of function mess_matrix_trace  ----- */



/**
 * @brief Compute frobenius inner product of two dense matrices.
 * @param[in] op input operation for A (transposed / hermitian / none)
 * @param[in] A dense input matrix
 * @param[in] B dense input matrix
 * @param[out] fro output pointer to resulting frobenius inner product
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_fro_inn function computes the frobenius inner product of two dense matrices:
 * \f[ fro = \langle op(A),B \rangle_F = \sum\limits_{i,j} \overline{op(A)_{i,j}} B_{i,j} \f]
 * The output argument is a void pointer. If both input matrices were real, it has to be casted to (double *),
 * otherwise to (@ref mess_double_cpx_t *).
 *
 * \attention This function completely relies on dense computations. It returns an error if it is called with
 * sparse matrices.
 * 
 */
int mess_matrix_fro_inn(mess_operation_t op, mess_matrix A, mess_matrix B, void *fro) {
    MSG_FNAME(__func__);
    mess_matrix Help=NULL;
    int ret,i,j;
    double * result=NULL;
    mess_double_cpx_t * resultc=NULL;

    mess_check_dense(A);
    mess_check_dense(B);
    if(op != MESS_OP_NONE){
        if (( (A)->rows != (B)->cols) || ( (A)->cols != (B)->rows)) { 
            MSG_ERROR("A^T and B must have the same size. (" MESS_PRINTF_INT ", " MESS_PRINTF_INT ")^T <-> (" MESS_PRINTF_INT ", " MESS_PRINTF_INT "))\n", A->rows,A->cols, B->rows, B->cols); 
            return (MESS_ERROR_DIMENSION); 
        }
    } else {
        mess_check_same_size(A,B);
    }
    if(MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
        ret = mess_matrix_init(&Help);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
        ret = mess_matrix_copy(B,Help);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);
        ret = mess_matrix_tocomplex(Help);                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_tocomplex);
        ret = mess_matrix_fro_inn(op, A, Help, fro);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_fro_inn);
        mess_matrix_clear(&Help);
        return 0;
    }
    if(MESS_IS_COMPLEX(B) && MESS_IS_REAL(A)){
        ret = mess_matrix_init(&Help);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
        ret = mess_matrix_copy(A,Help);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);
        ret = mess_matrix_tocomplex(Help);                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_tocomplex);
        ret = mess_matrix_fro_inn(op, Help, B, fro);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_fro_inn);
        mess_matrix_clear(&Help);
        return 0;
    }
    
    if(MESS_IS_REAL(A)){
        result = (double *) fro;
        *result = 0.0;
    } else {
        resultc = (mess_double_cpx_t *) fro;
        *resultc = 0.0;
    }
    
    
    if(op == MESS_OP_NONE){
        if(MESS_IS_REAL(A)){
            for(j=0;j<A->cols;j++){
                for(i=0;i<A->rows;i++){
                    *result = *result + (A->values[i+j*A->ld] * B->values[i+j*B->ld]);
                }
            }
        } else {
            for(j=0;j<A->cols;j++){
                for(i=0;i<A->rows;i++){
                    *resultc = *resultc + (conj(A->values_cpx[i+j*A->ld]) * B->values_cpx[i+j*B->ld]);
                }
            }
        }
    } else {
        if(MESS_IS_REAL(A)){
            for(i=0;i<A->rows;i++){
                for(j=0;j<A->cols;j++){
                    *result = *result + (A->values[j+i*A->ld] * B->values[i+j*B->ld]);
                }
            }
        } else {
            if(op == MESS_OP_HERMITIAN){
                for(i=0;i<A->rows;i++){
                    for(j=0;j<A->cols;j++){
                        *resultc = *resultc + (A->values_cpx[j+i*A->ld] * B->values_cpx[i+j*B->ld]);
                    }
                }
            } else {
                for(i=0;i<A->rows;i++){
                    for(j=0;j<A->cols;j++){
                        *resultc = *resultc + (conj(A->values_cpx[j+i*A->ld]) * B->values_cpx[i+j*B->ld]);
                    }
                }
            }
        }
    }

    return 0;
}