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
 * @file lib/matrix/gaxpy_kernels/csc.c
 * @brief Kernels for \f$ y\leftarrow Ax+y \f$ in @ref MESS_CSC storage.
 * @author @koehlerm
 */

/**
 * @internal
 * @brief Sparse matrix-vector product for general matrices in @ref MESS_CSC format.
 * @param[in] A         input matrix in @ref MESS_CSC storage
 * @param[in] x         input right hand side vector \f$x\f$
 * @param[in,out] y     input/output vector \f$y\f$
 * @return zero on success or a non zero error code
 *
 * The @ref __gaxpy_csc_ge function computes
 * \f[ y \leftarrow Ax +y \f]
 * for a general matrix in @ref MESS_CSC format. Because this is only the
 * kernel function it does not check any input or output arguments.
 *
 * @attention Internal use only.
 */
static int __gaxpy_csc_ge(mess_matrix A , mess_vector x, mess_vector y) {
    mess_int_t i,j;
    int ret;
    MSG_FNAME(__func__);

    /*-----------------------------------------------------------------------------
     *  All data is real.
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_REAL(A) && MESS_IS_REAL(x) && MESS_IS_REAL(y)) {
        double t;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,t)
#endif
        for ( j = 0; j < A->cols ; j++){
            t= x->values[j];
            for (i = A->colptr[j]; i < A->colptr[j+1]; i++){
#ifdef _OPENMP
#pragma omp atomic
#endif
                y->values[A->rowptr[i]] += t * A->values[i];
            }
        }

    }

    /*-----------------------------------------------------------------------------
     *  Real matrix and complex right hand side
     *-----------------------------------------------------------------------------*/
    else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(x)) {
        ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        mess_double_cpx_t t;
#ifdef _OPENMP
#pragma omp parallel for private(i,j,t)
#endif
        for ( j = 0; j < A->cols ; j++){
            t= x->values_cpx[j];
            for (i = A->colptr[j]; i < A->colptr[j+1]; i++){
#ifdef _OPENMP
#pragma omp critical
#endif
                y->values_cpx[A->rowptr[i]] += t * A->values[i];
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  Complex Matrix and real right hand side
     *-----------------------------------------------------------------------------*/
    else if ( MESS_IS_COMPLEX(A) && MESS_IS_REAL(x)){
        ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
#ifdef _OPENMP
        mess_double_cpx_t t;
#pragma omp parallel for private(i,j,t)
#else
        register mess_double_cpx_t t;
#endif
        for ( j = 0; j < A->cols ; j++){
            t= x->values[j];
            for (i = A->colptr[j]; i < A->colptr[j+1]; i++){
#ifdef _OPENMP
#pragma omp critical
#endif
                y->values_cpx[A->rowptr[i]] += t * A->values_cpx[i];
            }
        }
    } else {
        ret = mess_vector_tocomplex(x);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_tocomplex(y);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        if(MESS_IS_COMPLEX(A)){
#ifdef _OPENMP
            mess_double_cpx_t t;
#pragma omp parallel for private(i,j,t)
#else
            register mess_double_cpx_t t;
#endif
            for ( j = 0; j < A->cols ; j++){
                t= x->values_cpx[j];
                for (i = A->colptr[j]; i < A->colptr[j+1]; i++){
#ifdef _OPENMP
#pragma omp critical
#endif
                    y->values_cpx[A->rowptr[i]] += t * A->values_cpx[i];
                }
            }
        }else{
#ifdef _OPENMP
            mess_double_cpx_t t;
#pragma omp parallel for private(i,j,t)
#else
            register mess_double_cpx_t t;
#endif
            for ( j = 0; j < A->cols ; j++){
                t= x->values_cpx[j];
                for (i = A->colptr[j]; i < A->colptr[j+1]; i++){
#ifdef _OPENMP
#pragma omp critical
#endif
                    y->values_cpx[A->rowptr[i]] += t * A->values[i];
                }
            }

        }

    }
    return (0);
}

/**
 * @internal
 * @brief Hermitian transposed sparse matrix-vector product for general matrices in @ref MESS_CSC format.
 * @param[in] A    input matrix in Compressed Sparse Column storage
 * @param[in] x         input right hand side vector \f$x\f$
 * @param[in,out] y     input/output vector \f$y\f$
 * @return zero on success or a non zero error code
 *
 * The @ref __gaxpyh_csc_ge function computes
 * \f[ y \leftarrow A^H x +y \f]
 * for a general matrix in @ref MESS_CSC format. Because this is only the
 * kernel function it does not check any input or output arguments.
 *
 * @attention Internal use only.
 *
 */
static int __gaxpyh_csc_ge( mess_matrix  A, mess_vector x, mess_vector y){
    mess_int_t i,j;
    int ret = 0;
    MSG_FNAME(__func__);


    /*-----------------------------------------------------------------------------
     *  all data is real
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_REAL(A) && MESS_IS_REAL(x) && MESS_IS_REAL(y)){
#ifdef _OPENMP
        double t;
#pragma omp parallel for private ( i, j, t)
#else
        register double t;
#endif
        for ( i = 0 ; i < A->cols; i++){
            t = 0.0;
            for ( j = A->colptr[i]; j < A->colptr[i+1]; j++){
                t += A->values[j]*x->values[A->rowptr[j]];
            }
            y->values[i]=t + y->values[i];
        }
    }

    /*-----------------------------------------------------------------------------
     *  Real matrix and complex right hand side
     *-----------------------------------------------------------------------------*/
    else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(x)) {
        ret = mess_vector_tocomplex(y);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
#ifdef _OPENMP
        mess_double_cpx_t t;
#pragma omp parallel for private ( i, j, t)
#else
        register mess_double_cpx_t t;
#endif
        for ( i = 0 ; i < A->cols; i++){
            t = 0.0;
            for ( j = A->colptr[i]; j < A->colptr[i+1]; j++){
                t += A->values[j]*x->values_cpx[A->rowptr[j]];
            }
            y->values_cpx[i]=t+y->values_cpx[i];
        }
    }

    /*-----------------------------------------------------------------------------
     *  Complex matrix with real right hand side
     *-----------------------------------------------------------------------------*/
    else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(x)) {

        ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE( ret,(ret!=0), mess_vector_tocomplex);
#ifdef _OPENMP
        mess_double_cpx_t t;
#pragma omp parallel for private ( i, j, t)
#else
        register mess_double_cpx_t t;
#endif
        for ( i = 0 ; i < A->cols; i++){
            t = 0.0;
            for ( j = A->colptr[i]; j < A->colptr[i+1]; j++){
                t += conj(A->values_cpx[j])*x->values[A->rowptr[j]];
            }
            y->values_cpx[i]=t+y->values_cpx[i];
        }
    } else {
        ret = mess_vector_tocomplex(x);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0) , mess_vector_tocomplex);
        ret = mess_vector_tocomplex(y);     FUNCTION_FAILURE_HANDLE( ret,(ret!=0), mess_vector_tocomplex);
        if(MESS_IS_COMPLEX(A)){
#ifdef _OPENMP
            mess_double_cpx_t t;
#pragma omp parallel for private ( i, j, t)
#else
            register mess_double_cpx_t t;
#endif
            for ( i = 0 ; i < A->cols; i++){
                t = 0.0;
                for ( j = A->colptr[i]; j < A->colptr[i+1]; j++){
                    t += conj(A->values_cpx[j])*x->values_cpx[A->rowptr[j]];
                }
                y->values_cpx[i]=t+y->values_cpx[i];
            }
        }else{
#ifdef _OPENMP
            mess_double_cpx_t t;
#pragma omp parallel for private ( i, j, t)
#else
            register mess_double_cpx_t t;
#endif
            for ( i = 0 ; i < A->cols; i++){
                t = 0.0;
                for ( j = A->colptr[i]; j < A->colptr[i+1]; j++){
                    t += conj(A->values[j])*x->values_cpx[A->rowptr[j]];
                }
                y->values_cpx[i]=t+y->values_cpx[i];
            }
        }
    }
    return (0);
}

/**
 * @internal
 * @brief Transposed sparse matrix-vector product for general matrices in @ref MESS_CSC format.
 * @param[in] A         input matrix in Compressed Sparse Column storage
 * @param[in] x         input right hand side vector \f$x\f$
 * @param[in,out] y     input/output vector \f$y\f$
 * @return zero on success or a non zero error code
 *
 * The @ref __gaxpyt_csc_ge function computes
 * \f[ y \leftarrow A^Tx +y \f]
 * for a general matrix in @ref MESS_CSC. If the matrix is real this function
 * is the same as @ref __gaxpyh_csc_ge but in the case of a complex matrix it only transposes
 * it without using the complex conjugate values.
 *
 * Because this is only the
 * kernel function it does not check any input or output arguments.
 *
 * @attention Internal use only.
 */
static int __gaxpyt_csc_ge( mess_matrix  A, mess_vector x, mess_vector y){
    mess_int_t i,j;
    int ret = 0;
    MSG_FNAME(__func__);

    if ( MESS_IS_REAL(A) && MESS_IS_REAL(x) && MESS_IS_REAL(y)){
        return(__gaxpyh_csc_ge(A, x, y));
    } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(x)) {
        return(__gaxpyh_csc_ge(A,x,y));
    } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(x)) {
        ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE( ret,(ret!=0), mess_vector_tocomplex);
#ifdef _OPENMP
        mess_double_cpx_t t;
#pragma omp parallel for private ( i, j, t)
#else
        register mess_double_cpx_t t;
#endif
        for ( i = 0 ; i < A->cols; i++){
            t = 0.0;
            for ( j = A->colptr[i]; j < A->colptr[i+1]; j++){
                t += A->values_cpx[j]*x->values[A->rowptr[j]];
            }
            y->values_cpx[i]=t+y->values_cpx[i];
        }
    } else {
        ret = mess_vector_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0) , mess_vector_tocomplex);
        ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE( ret,(ret!=0), mess_vector_tocomplex);
        if(MESS_IS_COMPLEX(A)){
#ifdef _OPENMP
            mess_double_cpx_t t;
#pragma omp parallel for private ( i, j, t)
#else
            register mess_double_cpx_t t;
#endif
            for ( i = 0 ; i < A->cols; i++){
                t = 0.0;
                for ( j = A->colptr[i]; j < A->colptr[i+1]; j++){
                    t += A->values_cpx[j]*x->values_cpx[A->rowptr[j]];
                }
                y->values_cpx[i]=t+y->values_cpx[i];
            }
        }else{
#ifdef _OPENMP
            mess_double_cpx_t t;
#pragma omp parallel for private ( i, j, t)
#else
            register mess_double_cpx_t t;
#endif
            for ( i = 0 ; i < A->cols; i++){
                t = 0.0;
                for ( j = A->colptr[i]; j < A->colptr[i+1]; j++){
                    t += A->values[j]*x->values_cpx[A->rowptr[j]];
                }
                y->values_cpx[i]=t+y->values_cpx[i];
            }
        }
    }
    return (0);
}

