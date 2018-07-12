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
 * @file lib/matrix/mvp.c
 * @brief Compute the matrix-vector product \f$ y=Ax \f$.
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
#include <assert.h>


/*-----------------------------------------------------------------------------
 *  include the kernels
 *-----------------------------------------------------------------------------*/
#include "gaxpy_kernels/csr.c"
#include "gaxpy_kernels/csc.c"
#include "gaxpy_kernels/dense.c"


/**
 * @brief Compute a matrix-vector product.
 * @param[in] op     input operation applied to \f$A\f$
 * @param[in] A      input matrix \f$A\f$
 * @param[in] x      input right hand side vector \f$x\f$
 * @param[out] y    output vector \f$y\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_mvp function computes the matrix-vector product, i.e.
 * \f$ y = op(A)x \f$
 * where \f$ op(A) \f$ can be
 * <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 *  |@ref MESS_OP_HERMITIAN      |   \f$ op(A)=A^H\f$         |
 * </center>
 *
 * The function chooses automatically the right kernel for the multiplication depending
 * on the storage type of the matrix. Currently only general matrices without any symmetric
 * structure property are supported.
 * @sa mess_matrix_gaxpy
 *
 */
int mess_matrix_mvp(mess_operation_t op,mess_matrix A, mess_vector x,mess_vector y){
    MSG_FNAME(__func__);
    int ret = 0;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(y);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(x);
    mess_check_real_or_complex(y);
    mess_check_real_or_complex(A);
    mess_check_operation_type(op);


    switch( op) {

        /*-----------------------------------------------------------------------------
         *  y = Ax
         *-----------------------------------------------------------------------------*/
        case MESS_OP_NONE:
            {
                if (x->dim != A->cols) {
                    MSG_ERROR("The dimension of x doesn't match.\n");
                    MSG_ERROR(" A->cols = " MESS_PRINTF_INT " \t x->dim = " MESS_PRINTF_INT "\n", A->cols, x->dim);
                    return (MESS_ERROR_DIMENSION);
                }
                if (y->dim != A->rows) {
                    MSG_WARN("The dimension of y doesn't match. Resized.\n");
                    MSG_WARN(" current dim = " MESS_PRINTF_INT " need = " MESS_PRINTF_INT "\n", y->dim, A->rows);
                    ret = mess_vector_resize(y, A->rows);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
                }
                ret = mess_vector_zeros(y);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_zeros);

                if ( MESS_IS_COMPLEX(x) || MESS_IS_COMPLEX(A)) {
                    ret =  mess_vector_tocomplex(y);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
                } else {
                    ret =  mess_vector_toreal_nowarn(y);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
                }
                switch (A->store_type){
                    case MESS_CSR:{
                                      switch (A->symmetry) {
                                          case MESS_GENERAL:
                                              ret = __gaxpy_csr_ge(A, x, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpy_csr_ge);
                                              break;
                                          case MESS_SYMMETRIC:
                                              ret = __gaxpy_csr_sym(A, x, y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpy_csr_sym);
                                              break;
                                          case MESS_SKEWSYMMETRIC:
                                              ret = __gaxpy_csr_sym(A, x, y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpy_csr_sym);
                                              break;
                                          default:
                                              MSG_ERROR("Unknown symmetry type: %s\n", mess_symmetry_t_str(A->symmetry));
                                              return (MESS_ERROR_SYMMETRIC);
                                      }
                                      break;
                                  }
                    case MESS_CSC: {
                                       switch (A->symmetry) {
                                           case MESS_GENERAL:
                                               ret = __gaxpy_csc_ge(A, x, y);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpy_csc_ge);
                                               break;
                                           case MESS_SYMMETRIC:
                                           case MESS_SKEWSYMMETRIC:
                                               MSG_ERROR("Matrix-Vector product for Symmertric CSC matrices is not implemented yet.\n");
                                               return MESS_ERROR_MISSING;
                                           default:
                                               MSG_ERROR("Unknown symmetry type: %s\n", mess_symmetry_t_str(A->symmetry));
                                               return (MESS_ERROR_SYMMETRIC);
                                       }
                                       break;
                                   }
                    case MESS_COORD: {
                                         mess_matrix tmp;
                                         MESS_MATRIX_CHECKFORMAT(A, tmp, ret, MESS_CSR);
                                         if ( ret == 0 ) {
                                             ret = mess_matrix_mvp(MESS_OP_NONE,tmp,x, y);
                                             mess_matrix_clear(&tmp);
                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                                         }else {
                                             MSG_ERROR ("error converting MESS_COORD to MESS_CSR\n");
                                             return(MESS_ERROR_CONVERT);
                                         }
                                         break;
                                     }
                    case MESS_DENSE:{
                                        ret = __gaxpy_dense(A, x, y);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpy_dense);
                                        break;
                                    }

                    default:
                                    MSG_ERROR("unkown/unsupported storage type: %s\n", mess_storage_t_str(A->store_type));
                                    return (MESS_ERROR_STORAGETYPE);
                }
            }
            break;

            /*-----------------------------------------------------------------------------
             *  y = A'x with ' hermitian transpose
             *-----------------------------------------------------------------------------*/
        case MESS_OP_HERMITIAN:
            {
                if (x->dim != A->rows) {
                    MSG_ERROR("The dimension of x doesn't match. x->dim = " MESS_PRINTF_INT " matrix->rows = " MESS_PRINTF_INT "\n", x->dim, A->rows);
                    return (MESS_ERROR_DIMENSION);
                }
                if (y->dim != A->cols) {
                    MSG_WARN("The dimension of y doesn't match. Resized.\n");
                    MSG_WARN(" current dim = " MESS_PRINTF_INT " need = " MESS_PRINTF_INT "\n", y->dim, A->cols);
                    ret = mess_vector_resize(y, A->cols); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
                }

                if ( MESS_IS_REAL(A) && MESS_IS_REAL(x)){
                    ret =  mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
                }  else {
                    ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
                }

                ret = mess_vector_zeros(y);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_zeros);
                switch ( A->store_type){
                    case MESS_CSR: {
                                       switch (A->symmetry) {
                                           case MESS_GENERAL:
                                               ret = __gaxpyh_csr_ge(A, x, y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpyh_csr_ge);
                                               break;
                                           case MESS_SYMMETRIC:
                                               ret = __gaxpy_csr_sym(A, x, y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpy_csr_sym);
                                               break;
                                           case MESS_SKEWSYMMETRIC:
                                               ret = __gaxpy_csr_sym(A, x, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpy_csr_sym);
                                               break;
                                           default:
                                               MSG_ERROR("Unknown symmetry type: %s\n", mess_symmetry_t_str(A->symmetry));
                                               return (MESS_ERROR_SYMMETRIC);
                                       }
                                       break;
                                   }
                    case MESS_CSC:{
                                      switch (A->symmetry) {
                                          case MESS_GENERAL:
                                              ret = __gaxpyh_csc_ge(A, x, y);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpyh_csc_ge);
                                              break;
                                          case MESS_SYMMETRIC:
                                          case MESS_SKEWSYMMETRIC:
                                              MSG_ERROR("Matrix-Vector product for Symmertric CSC matrices is not implemented yet.\n");
                                              return MESS_ERROR_MISSING;
                                          default:
                                              MSG_ERROR("Unknown symmetry type: %s\n", mess_symmetry_t_str(A->symmetry));
                                              return (MESS_ERROR_SYMMETRIC);
                                      }
                                      break;
                                  }

                    case MESS_COORD:{
                                        mess_matrix tmp;
                                        MESS_MATRIX_CHECKFORMAT(A, tmp, ret, MESS_CSR);
                                        if ( ret == 0 ) {
                                            ret = mess_matrix_mvp(MESS_OP_HERMITIAN,tmp,x, y);
                                            mess_matrix_clear(&tmp);
                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                                        }else {
                                            MSG_ERROR ("error converting MESS_COORD to MESS_CSR\n");
                                            return(MESS_ERROR_CONVERT);
                                        }
                                        break;
                                    }
                    case MESS_DENSE:{
                                        ret = __gaxpyh_dense(A, x, y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpyh_dense);
                                        break;
                                    }
                    default:
                                    MSG_ERROR("unknown/unsupported storage type: %s\n", mess_storage_t_str(A->store_type));
                                    return(MESS_ERROR_STORAGETYPE);
                }

            }
            break;

            /*-----------------------------------------------------------------------------
             *  y = A'x with ' tranpose with out conj.
             *-----------------------------------------------------------------------------*/
        case MESS_OP_TRANSPOSE:
            {
                if (x->dim != A->rows) {
                    MSG_ERROR("The dimension of x doesn't match. x->dim = " MESS_PRINTF_INT " matrix->rows = " MESS_PRINTF_INT "\n", x->dim, A->rows);
                    return( MESS_ERROR_DIMENSION);
                }
                if (y->dim != A->cols) {
                    MSG_WARN("The dimension of y doesn't match. Resized.\n");
                    MSG_WARN(" current dim = " MESS_PRINTF_INT " need = " MESS_PRINTF_INT "\n", y->dim, A->cols);
                    ret = mess_vector_resize(y, A->cols); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
                }

                if ( MESS_IS_REAL(A) && MESS_IS_REAL(x)){
                    ret =  mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
                }  else {
                    ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
                }

                ret = mess_vector_zeros(y);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_zeros);
                switch ( A->store_type){
                    case MESS_CSR: {
                                       switch (A->symmetry) {
                                           case MESS_GENERAL:
                                               ret = __gaxpyt_csr_ge(A, x, y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpyt_csr_ge);
                                               break;
                                           default:
                                               MSG_ERROR("Unknown/Unsupported symmetry type: %s\n", mess_symmetry_t_str(A->symmetry));
                                               return (MESS_ERROR_SYMMETRIC);
                                       }
                                       break;
                                   }
                    case MESS_CSC:{
                                      switch (A->symmetry) {
                                          case MESS_GENERAL:
                                              ret = __gaxpyt_csc_ge(A, x, y);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpyt_csc_ge);
                                              break;
                                          default:
                                              MSG_ERROR("Unknown/Unsupported symmetry type: %s\n", mess_symmetry_t_str(A->symmetry));
                                              return (MESS_ERROR_SYMMETRIC);

                                      }
                                      break;
                                  }

                    case MESS_COORD:{
                                        mess_matrix tmp;
                                        MESS_MATRIX_CHECKFORMAT(A, tmp, ret, MESS_CSR);
                                        if ( ret == 0 ) {
                                            ret = mess_matrix_mvp(MESS_OP_TRANSPOSE,tmp,x, y);
                                            mess_matrix_clear(&tmp);
                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                                        }else {
                                            MSG_ERROR ("error converting MESS_COORD to MESS_CSR\n");
                                            return (MESS_ERROR_CONVERT);
                                        }
                                        break;
                                    }
                    case MESS_DENSE:{
                                        ret = __gaxpyt_dense(A, x, y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __gaxpyt_dense);
                                        break;
                                    }
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

    return (ret);
}


