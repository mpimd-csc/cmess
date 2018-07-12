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
 * @file lib/matrix/gaxpy_kernels/dense.c
 * @brief Kernels for \f$ y\leftarrow Ax+y \f$ in @ref MESS_DENSE storage.
 * @author @koehlerm
 */


/**
 * @internal
 * @brief Dense matrix-vector product for general matrices in @ref MESS_DENSE fomart.
 * @param[in] matrix    input matrix in DENSE storage
 * @param[in] x         input right hand side vector
 * @param[in,out] y     input/output vector
 * @return zero on success or a non zero error code
 *
 * The @ref __gaxpy_dense function computes
 * \f[ y \leftarrow Ax +y \f]
 * for a general matrix in @ref MESS_DENSE format. Because this is only the
 * kernel function it does not check any input or output arguments.
 *
 * @attention Internal use only.
 */
static int __gaxpy_dense(mess_matrix matrix,  mess_vector x, mess_vector y) {
    MSG_FNAME(__func__);

    mess_int_t one = 1;
    mess_int_t two = 2;
    int ret = 0;

    if (MESS_IS_REAL(matrix) && MESS_IS_REAL(x) && MESS_IS_REAL(y)){
        double alpha = 1.0;
        double beta = 1.0;
        F77_GLOBAL(dgemv,DGEMV)("N",&matrix->rows, &matrix->cols, &alpha, matrix->values, &matrix->ld, x->values, &one, &beta, y->values, &one );
    } else if ( MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(x) ) {
        double alpha = 1.0;
        double beta = 1;
        ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(dgemv,DGEMV)("N",&matrix->rows, &matrix->cols, &alpha, matrix->values, &matrix->ld, (double *) x->values_cpx, &two, &beta, (double *) y->values_cpx, &two );
        F77_GLOBAL(dgemv,DGEMV)("N",&matrix->rows, &matrix->cols, &alpha, matrix->values, &matrix->ld, (double *) x->values_cpx+1, &two, &beta, (double *) y->values_cpx+1, &two );
    } else if ( MESS_IS_COMPLEX(matrix)&& MESS_IS_REAL(x)) {
        mess_double_cpx_t *x2;
        mess_double_cpx_t alpha = 1.0;
        mess_double_cpx_t beta = 1.0;
        mess_int_t i;
        mess_try_alloc(x2, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*x->dim);
        for (i=0; i < x->dim ; i++){ x2[i] = x->values[i];  }
        ret =   mess_vector_tocomplex(y);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(zgemv,ZGEMV)("N",&matrix->rows, &matrix->cols, &alpha, matrix->values_cpx, &matrix->ld, x2, &one, &beta, y->values_cpx, &one );
        mess_free(x2);
    }  else if ( MESS_IS_COMPLEX(matrix)&& MESS_IS_COMPLEX(x)) {
        mess_double_cpx_t alpha = 1.0;
        mess_double_cpx_t beta = 1.0;
        ret =   mess_vector_tocomplex(y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(zgemv,ZGEMV)("N",&matrix->rows, &matrix->cols, &alpha, matrix->values_cpx, &matrix->ld, x->values_cpx, &one, &beta, y->values_cpx, &one );
    }else if(MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(y)){
        double alpha = 1.0;
        double beta = 1.0;
        F77_GLOBAL(dgemv,DGEMV)("N",&matrix->rows, &matrix->cols, &alpha, matrix->values, &matrix->ld, x->values, &one, &beta, (double *) y->values_cpx, &two );
    }

    return (0);
}


/**
 * @internal
 * @brief Hermitian transposed dense matrix-vector product for general matrices in @ref MESS_DENSE fomart.
 * @param[in] matrix    input matrix in DENSE storage
 * @param[in] x         input right hand side vector
 * @param[in,out] y     input/output vector
 * @return zero on success or a non zero error code
 *
 * The @ref __gaxpyh_dense function computes
 * \f[ y \leftarrow A^H x +y \f]
 * for a general matrix in @ref MESS_DENSE format. Because this is only the
 * kernel function it does not check any input or output arguments.
 *
 *
 * @attention Internal use only.
 */
static int __gaxpyh_dense(mess_matrix matrix, mess_vector x, mess_vector y) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t one = 1;
    mess_int_t two = 2;

    if ( MESS_IS_REAL(matrix) && MESS_IS_REAL(x) && MESS_IS_REAL(y)) {
        double alpha = 1.0;
        double beta = 1.0;
        F77_GLOBAL(dgemv,DGEMV)("T",&matrix->rows, &matrix->cols, &alpha, matrix->values, &matrix->ld, x->values, &one, &beta, y->values, &one );
    } else if ( MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(x) ) {
        double alpha = 1.0;
        double beta = 1;
        ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(dgemv,DGEMV)("T",&matrix->rows, &matrix->cols, &alpha, matrix->values, &matrix->ld, (double *) x->values_cpx, &two, &beta, (double *) y->values_cpx, &two );
        F77_GLOBAL(dgemv,DGEMV)("T",&matrix->rows, &matrix->cols, &alpha, matrix->values, &matrix->ld, (double *) x->values_cpx+1, &two, &beta, (double *) y->values_cpx+1, &two );
    }  else if ( MESS_IS_COMPLEX(matrix)&& MESS_IS_REAL(x)) {
        mess_double_cpx_t *x2;
        mess_double_cpx_t alpha = 1.0;
        mess_double_cpx_t beta = 1.0;
        mess_int_t i;
        mess_try_alloc(x2, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*x->dim);
        for (i=0; i < x->dim ; i++){ x2[i] = x->values[i];  }
        ret =   mess_vector_tocomplex(y);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(zgemv,ZGEMV)("C",&matrix->rows, &matrix->cols, &alpha, matrix->values_cpx, &matrix->ld, x2, &one, &beta, y->values_cpx, &one );
        mess_free(x2);
    } else if (MESS_IS_COMPLEX(matrix) && MESS_IS_COMPLEX(x)) {
        mess_double_cpx_t alpha = 1.0;
        mess_double_cpx_t beta = 1.0;
        ret =  mess_vector_tocomplex(y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(zgemv,ZGEMV)("C",&matrix->rows, &matrix->cols, &alpha, matrix->values_cpx, &matrix->ld, x->values_cpx, &one, &beta, y->values_cpx, &one );
    } else if (MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(y)){
        double alpha = 1.0;
        double beta = 1.0;
        F77_GLOBAL(dgemv,DGEMV)("T",&matrix->rows, &matrix->cols, &alpha, matrix->values, &matrix->ld, x->values, &one, &beta, (double *) y->values_cpx, &two );
    }

    return (0);
}

/**
 * @internal
 * @brief Transpose dense matrix-vector product for general matrices in DENSE fomart.
 * @param[in] matrix     input matrix in DENSE storage
 * @param[in] x      input right hand side vector
 * @param[in,out] y     output vector
 * @return zero on success or a non zero error code
 *
 * The @ref __gaxpyt_dense function computes
 * \f[ y \leftarrow A^Tx +y \f]
 * for a general matrix in @ref MESS_DENSE format. If the matrix is real this function
 * is the same as @ref __gaxpyh_dense but in the case of a complex matrix it only transposes
 * it without using the complex conjugate values.
 *
 * Because this is only the
 * kernel function it does not check any input or output arguments.
 *
 * @attention Internal use only.
 */
static int __gaxpyt_dense(mess_matrix matrix, mess_vector x, mess_vector y) {
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t one = 1;
    mess_int_t two = 2;
    if ( MESS_IS_REAL(matrix) && MESS_IS_REAL(x) && MESS_IS_REAL(y)) {
        double alpha = 1.0;
        double beta = 1.0;
        F77_GLOBAL(dgemv,DGEMV)("T",&matrix->rows, &matrix->cols, &alpha, matrix->values, &matrix->ld, x->values, &one, &beta, y->values, &one );
    } else if ( MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(x) ) {
        double alpha = 1.0;
        double beta = 1;
        ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(dgemv,DGEMV)("T",&matrix->rows, &matrix->cols, &alpha, matrix->values, &matrix->ld, (double *) x->values_cpx, &two, &beta, (double *) y->values_cpx, &two );
        F77_GLOBAL(dgemv,DGEMV)("T",&matrix->rows, &matrix->cols, &alpha, matrix->values, &matrix->ld, (double *) x->values_cpx+1, &two, &beta, (double *) y->values_cpx+1, &two );
    }  else if ( MESS_IS_COMPLEX(matrix)&& MESS_IS_REAL(x)) {
        mess_double_cpx_t *x2;
        mess_double_cpx_t alpha = 1.0;
        mess_double_cpx_t beta = 1.0;
        mess_int_t i;
        mess_try_alloc(x2, mess_double_cpx_t *, sizeof(mess_double_cpx_t)*x->dim);
        for (i=0; i < x->dim ; i++){ x2[i] = x->values[i];  }
        ret =   mess_vector_tocomplex(y);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(zgemv,ZGEMV)("T",&matrix->rows, &matrix->cols, &alpha, matrix->values_cpx, &matrix->ld, x2, &one, &beta, y->values_cpx, &one );
        mess_free(x2);
    } else if (MESS_IS_COMPLEX(matrix) && MESS_IS_COMPLEX(x)) {
        mess_double_cpx_t alpha = 1.0;
        mess_double_cpx_t beta = 1.0;
        ret =   mess_vector_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret =  mess_vector_tocomplex(y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(zgemv,ZGEMV)("T",&matrix->rows, &matrix->cols, &alpha, matrix->values_cpx, &matrix->ld, x->values_cpx, &one, &beta, y->values_cpx, &one );
    }  else if (MESS_IS_REAL(matrix) && MESS_IS_COMPLEX(y)){
        double alpha = 1.0;
        double beta = 1.0;
        F77_GLOBAL(dgemv,DGEMV)("T",&matrix->rows, &matrix->cols, &alpha, matrix->values, &matrix->ld, x->values, &one, &beta, (double *) y->values_cpx, &two );
    }

    return (0);
}


