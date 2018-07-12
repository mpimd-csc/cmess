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
 * @file lib/matrix/gaxpy.c
 * @brief Compute the matrix-vector product \f$ y \leftarrow y + op(A)x \f$
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>

/*-----------------------------------------------------------------------------
 *  include kernels
 *-----------------------------------------------------------------------------*/
#include "gaxpy_kernels/csr.c"
#include "gaxpy_kernels/csc.c"
#include "gaxpy_kernels/dense.c"


/**
 * @brief Compute a matrix-vector product with update of the output.
 * @param[in] op    input operation on matrix
 * @param[in] A     input matrix  \f$A\f$
 * @param[in] x     input right hand side vector  \f$x\f$
 * @param[out] y    output vector \f$y\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_gaxpy function computes the matrix-vector product with an update of the output, i.e.
 * \f[ y \leftarrow y + op(A)x \f]
 * where \f$ op(A) \f$ can be
 * <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 *  |@ref MESS_OP_HERMITIAN      |   \f$ op(A)=A^H\f$         |
 * </center>
 *
 * @attention This function does not yet support the MESS_COORD storage type.
 *
 * @sa mess_matrix_mvp
 *
 */
int mess_matrix_gaxpy ( mess_operation_t op, mess_matrix A, mess_vector x, mess_vector y )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(x);
    mess_check_nullpointer(y);
    mess_check_real_or_complex(A);
    mess_check_real_or_complex(x);
    mess_check_real_or_complex(y);
    mess_check_operation_type(op);

    /*-----------------------------------------------------------------------------
     * select and run the kernel
     *-----------------------------------------------------------------------------*/
    switch( op) {

        /*-----------------------------------------------------------------------------
         *  y = Ax
         *-----------------------------------------------------------------------------*/
        case MESS_OP_NONE:
            {
                if ( x->dim != A->cols) {
                    MSG_ERROR("dimension error: x->dim = " MESS_PRINTF_INT ", matrix->cols = " MESS_PRINTF_INT "\n", x->dim, A->cols);
                    return(MESS_ERROR_DIMENSION);
                }
                if ( y->dim != A->rows) {
                    MSG_ERROR("dimension error: y->dim = " MESS_PRINTF_INT ", matrix->rows = " MESS_PRINTF_INT "\n", x->dim, A->rows);
                    return(MESS_ERROR_DIMENSION);
                }
                switch ( A->store_type ){
                    case MESS_CSR:
                        ret = __gaxpy_csr_ge(A, x, y);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpy_csr_ge);
                        break;
                    case MESS_CSC:
                        ret = __gaxpy_csc_ge(A, x, y);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpy_csc_ge);
                        break;
                    case MESS_DENSE:
                        ret = __gaxpy_dense(A, x, y);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpy_dense);
                        break;
                    default:
                        MSG_ERROR("unknown/unsupported storage type: %s\n", mess_storage_t_str(A->store_type));
                        return(MESS_ERROR_STORAGETYPE);
                }
            }
            break;

            /*-----------------------------------------------------------------------------
             *  y = A'x with ' hermitian transpose
             *-----------------------------------------------------------------------------*/
        case MESS_OP_HERMITIAN:
            {
                if ( x->dim != A->rows) {
                    MSG_ERROR("dimension error: x->dim = " MESS_PRINTF_INT ", matrix->rows = " MESS_PRINTF_INT "\n", x->dim, A->rows);
                    return(MESS_ERROR_DIMENSION);
                }
                if ( y->dim != A->cols) {
                    MSG_ERROR("dimension error: y->dim = " MESS_PRINTF_INT ", matrix->cols = " MESS_PRINTF_INT "\n", x->dim, A->cols);
                    return(MESS_ERROR_DIMENSION);
                }

                switch ( A->store_type ){
                    case MESS_CSR:
                        ret = __gaxpyh_csr_ge(A, x, y);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpyh_csr_ge);
                        break;
                    case MESS_CSC:
                        ret = __gaxpyh_csc_ge(A, x, y);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpyh_csc_ge);
                        break;
                    case MESS_DENSE:
                        ret = __gaxpyh_dense(A, x, y);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpyh_dense);
                        break;
                    default:
                        MSG_ERROR("unknown/unsupported storage type: %s\n", mess_storage_t_str(A->store_type));
                        return (MESS_ERROR_STORAGETYPE);
                }
            }
            break;

            /*-----------------------------------------------------------------------------
             *  y = A'x with ' tranpose with out conj.
             *-----------------------------------------------------------------------------*/
        case MESS_OP_TRANSPOSE:
            {
                if ( x->dim != A->rows) {
                    MSG_ERROR("dimension error: x->dim = " MESS_PRINTF_INT ", matrix->rows = " MESS_PRINTF_INT "\n", x->dim, A->rows);
                    return(MESS_ERROR_DIMENSION);
                }
                if ( y->dim != A->cols) {
                    MSG_ERROR("dimension error: y->dim = " MESS_PRINTF_INT ", matrix->cols = " MESS_PRINTF_INT "\n", x->dim, A->cols);
                    return(MESS_ERROR_DIMENSION);
                }

                switch ( A->store_type ){
                    case MESS_CSR:
                        ret = __gaxpyt_csr_ge(A, x, y);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpyt_csr_ge);
                        break;
                    case MESS_CSC:
                        ret = __gaxpyt_csc_ge(A, x, y);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpyt_csc_ge);
                        break;
                    case MESS_DENSE:
                        ret = __gaxpyt_dense(A, x, y);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpyt_dense);
                        break;
                    default:
                        MSG_ERROR("unknown/unsupported storage type: %s\n", mess_storage_t_str(A->store_type));
                        return (MESS_ERROR_STORAGETYPE);
                }

            }
            break;
        default:
            // never called because of mess_check_operation_type macro.
            ;

    }
    return(ret);
}       /* ----- end of function mess_matrix_gaxpy  ----- */


