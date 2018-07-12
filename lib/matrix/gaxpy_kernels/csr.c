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
 * @file lib/matrix/gaxpy_kernels/csr.c
 * @brief Kernels for \f$ y\leftarrow Ax+y \f$ in @ref MESS_CSR storage.
 * @author @koehlerm
 */

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * @internal
 * @brief Sparse matrix-vector product for general matrices in @ref MESS_CSR format.
 * @param[in] A         input matrix in Compressed Sparse Row storage
 * @param[in] x         input right hand side vector \f$x\f$
 * @param[in,out] y     input/output vector \f$y\f$
 * @return zero on success or a non zero error code
 *
 * The @ref __gaxpy_csr_ge function computes
 * \f[ y \leftarrow Ax +y \f]
 * for a general matrix in @ref MESS_CSR format. Because this is only the
 * kernel function it does not check any input or output arguments.
 *
 * @attention Internal use only.
 */
static int __gaxpy_csr_ge( mess_matrix A, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    mess_int_t i,j;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  ALl data is real
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_REAL(A) && MESS_IS_REAL(x) && MESS_IS_REAL(y)) {
#ifdef _OPENMP
        double t;
#pragma omp parallel for private ( i, j, t)
#else
        register double t;
#endif
        for ( i = 0 ; i < A->rows; i++){
            t = 0.0;
            for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                t += A->values[j]*x->values[A->colptr[j]];
            }
            y->values[i]=t+y->values[i];
        }
    }

    /*-----------------------------------------------------------------------------
     *  Real matrix and complex right hand side
     *-----------------------------------------------------------------------------*/
    else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(x)) {
        ret = mess_vector_tocomplex(y);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
#ifdef _OPENMP
        mess_double_cpx_t t;
#pragma omp parallel for private ( i, j, t)
#else
        register mess_double_cpx_t t;
#endif
        for ( i = 0 ; i < A->rows; i++){
            t = 0.0;
            for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                t += A->values[j]*x->values_cpx[A->colptr[j]];
            }
            y->values_cpx[i]=t+y->values_cpx[i];
        }
    }

    /*-----------------------------------------------------------------------------
     *  Complex matrix with real right hand side
     *-----------------------------------------------------------------------------*/
    else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(x)){
        ret = mess_vector_tocomplex(y);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
#ifdef _OPENMP
        mess_double_cpx_t t;
#pragma omp parallel for private ( i, j, t)
#else
        register mess_double_cpx_t t;
#endif
        for ( i = 0 ; i < A->rows; i++){
            t = 0.0;
            for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                t += A->values_cpx[j]*x->values[A->colptr[j]];
            }
            y->values_cpx[i]=t+y->values_cpx[i];
        }
    } else {
        ret = mess_vector_tocomplex(x);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_tocomplex(y);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);

        if(MESS_IS_COMPLEX(A)){
#ifdef _OPENMP
            mess_double_cpx_t t;
#pragma omp parallel for private ( i, j, t)
#else
            register mess_double_cpx_t t;
#endif
            for ( i = 0 ; i < A->rows; i++){
                t = 0.0;
                for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                    t += A->values_cpx[j]*x->values_cpx[A->colptr[j]];
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
            for ( i = 0 ; i < A->rows; i++){
                t = 0.0;
                for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                    t += A->values[j]*x->values_cpx[A->colptr[j]];
                }
                y->values_cpx[i]=t+y->values_cpx[i];
            }

        }
    }
    return (0);
}

/**
 * @internal
 * @brief Hermitian transposed sparse matrix-vector product for general matrices in @ref MESS_CSR format.
 * @param[in] A         input matrix in Compressed Sparse Row storage
 * @param[in] x         input right hand side vector \f$x\f$
 * @param[in,out] y     input/output vector \f$y\f$
 * @return zero on success or a non zero error code
 *
 * The @ref __gaxpyh_csr_ge function computes
 * \f[ y \leftarrow A^H x +y \f]
 * for a general matrix in @ref MESS_CSR format. Because this is only the
 * kernel function it does not check any input or output arguments.
 *
 * @attention Internal use only.
 * @warning untested parallel code
 */
static int __gaxpyh_csr_ge( mess_matrix A, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    mess_int_t i,j;
    int ret = 0;


    /*-----------------------------------------------------------------------------
     *  All data is real
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_REAL(A) && MESS_IS_REAL(x)) {

        /*-----------------------------------------------------------------------------
         *  Real Output vector
         *-----------------------------------------------------------------------------*/
        if (MESS_IS_REAL(y)) {
#ifdef _OPENMP
            double r,t;
#pragma omp parallel for private(i,j,r,t)
#else
            register double r, t;
#endif
            for ( i = 0 ; i < A->rows; i++){
                t = x->values[i];
                for (j = A->rowptr[i]; j<A->rowptr[i+1]; j++){
                    r = A->values[j]*t;
#ifdef _OPENMP
#pragma omp atomic
#endif
                    y->values[A->colptr[j]] += r;
                }
            }
        }

        /*-----------------------------------------------------------------------------
         *  Complex Output vector
         *-----------------------------------------------------------------------------*/
        else {
#ifdef _OPENMP
            double r,t;
#pragma omp parallel for private(i,j,r,t)
#else
            register double r,t;
#endif
            for ( i = 0 ; i < A->rows; i++){
                t = x->values[i];
                for (j = A->rowptr[i]; j<A->rowptr[i+1]; j++){
                    r = A->values[j]*t;
#ifdef _OPENMP
#pragma omp critical
#endif
                    y->values_cpx[A->colptr[j]] += r;
                }
            }
        }
        /*-----------------------------------------------------------------------------
         *
         *-----------------------------------------------------------------------------*/
    } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(x)) {
        ret = mess_vector_tocomplex(y);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);

#ifdef _OPENMP
        mess_double_cpx_t t,r;
#pragma omp parallel for private(i,j,r,t)
#else
        register mess_double_cpx_t t,r;
#endif
        for ( i = 0 ; i < A->rows; i++){
            t = x->values_cpx[i];
            for (j = A->rowptr[i]; j<A->rowptr[i+1]; j++){
                r = A->values[j]*t;
#ifdef _OPENMP
#pragma omp critical
#endif
                y->values_cpx[A->colptr[j]] += r;
            }
        }

    } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL (x)) {
        ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
#ifdef _OPENMP
        mess_double_cpx_t t,r;
#pragma omp parallel for private(i,j,r,t)
#else
        register mess_double_cpx_t t,r;
#endif
        for ( i = 0 ; i < A->rows; i++){
            t = x->values[i];
            for (j = A->rowptr[i]; j<A->rowptr[i+1]; j++){
                r = conj(A->values_cpx[j])*t;
#ifdef _OPENMP
#pragma omp critical
#endif
                y->values_cpx[A->colptr[j]] += r;
            }
        }
    } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(x)){ // COMPLEX matrix COMPLEX x
        ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
#ifdef _OPENMP
        mess_double_cpx_t t,r;
#pragma omp parallel for private(i,j,r,t)
#else
        register mess_double_cpx_t t,r;
#endif
        for ( i = 0 ; i < A->rows; i++){
            t = x->values_cpx[i];
            for (j = A->rowptr[i]; j<A->rowptr[i+1]; j++){
                r = conj(A->values_cpx[j])*t;
#ifdef _OPENMP
#pragma omp critical
#endif
                y->values_cpx[A->colptr[j]] += r;
            }
        }
    } else {
        MSG_ERROR("unknown data types matrix=%s, x=%s\n", mess_datatype_t_str(A->data_type), mess_datatype_t_str(x->data_type));
        return (MESS_ERROR_DATATYPE);
    }
    return (0);
}


/**
 * @internal
 * @brief Hermitian transposed sparse matrix-vector product for general matrices in @ref MESS_CSR format.
 * @param[in] A         input matrix in Compressed Sparse Row storage
 * @param[in] x         input right hand side vector \f$x\f$
 * @param[in,out] y     input/output vector \f$y\f$
 * @return zero on success or a non zero error code
 *
 * The @ref __gaxpyt_csr_ge function computes
 * \f[ y \leftarrow A^T x +y \f]
 * for a general matrix in @ref MESS_CSR format. Because this is only the
 * kernel function it does not check any input or output arguments. In the
 * case of a real matrix this is the same as  \ref __gaxpyh_csr_ge.
 *
 * @attention Internal use only.
 */
static int __gaxpyt_csr_ge( mess_matrix A, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    mess_int_t i,j;
    int ret = 0;

    if ( MESS_IS_REAL(A) && MESS_IS_REAL(x) ) {
        return (__gaxpyh_csr_ge(A, x, y));
    } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(x)) {
        return (__gaxpyh_csr_ge(A, x, y));

    } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL (x)) {
        ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);

#ifdef _OPENMP
        mess_double_cpx_t t;
#pragma omp parallel for private(t,i,j)
#else
        register mess_double_cpx_t t;
#endif
        for ( i = 0 ; i < A->rows; i++){
            t = x->values[i];
            for (j = A->rowptr[i]; j<A->rowptr[i+1]; j++){
#ifdef _OPENMP
#pragma omp critical
#endif
                y->values_cpx[A->colptr[j]] += A->values_cpx[j]*t;;
            }
        }
    } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(x)){ // COMPLEX matrix COMPLEX x

        ret = mess_vector_tocomplex(y);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
#ifdef _OPENMP
        mess_double_cpx_t t;
#pragma omp parallel for private(t,i,j)
#else
        register mess_double_cpx_t t;
#endif
        for ( i = 0 ; i < A->rows; i++){
            t = x->values_cpx[i];
            for (j = A->rowptr[i]; j<A->rowptr[i+1]; j++){
#ifdef _OPENMP
#pragma omp critical
#endif
                y->values_cpx[A->colptr[j]] += A->values_cpx[j]*t;
            }
        }

    } else {
        MSG_ERROR("unknown data types matrix=%s, x=%s\n", mess_datatype_t_str(A->data_type), mess_datatype_t_str(x->data_type));
        return(MESS_ERROR_DATATYPE);
    }
    return(0);
}


/**
 * @internal
 * @brief Sparse matrix-vector product for (skew-)symmertric matrices in @ref MESS_CSR format.
 * @param[in] A         input matrix in Compressed Sparse Row storage
 * @param[in] xv        input right hand side vector
 * @param[in,out] yv    input/output vector
 * @return zero on success or a non zero error code
 *
 * The @ref __gaxpy_csr_sym function computes
 * \f[ y \leftarrow Ax +y \f]
 * for a (skew-)symmetric matrix in @ref MESS_CSR format. Because this is only the
 * kernel function it does not check any input or output arguments.
 *
 * @attention Internal use only.
 */
__attribute__((unused)) static int __gaxpy_csr_sym( mess_matrix A, mess_vector xv, mess_vector yv){
    MSG_FNAME(__func__);
    mess_int_t i,j;
    double t, r;
    double skew = 1.0;
    int num_th;
    double *x = xv->values;
    double *y = yv->values;
    mess_check_real(A);
    if ( MESS_IS_SKEWSYMMETRIC ( A) ) skew = -1.0;

#ifdef _OPENMP
    num_th = omp_get_max_threads();
#else
    num_th = 1;
#endif
    if ( num_th == 1){
        for ( i = 0 ; i < A->rows; i++){
            t = 0.0;
            for (j = A->rowptr[i]; j<A->rowptr[i+1]; j++){
                t += (A->values[j]*x[A->colptr[j]]);
                if ( i != A->colptr[j] ){
                    y[A->colptr[j]] += skew * A->values[j] * x[i];
                }
            }
            y[i] += t;
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for private (i, j, t)
#endif
        for ( i = 0 ; i < A->rows; i++){
            t = 0.0 ;
            for (j = A->rowptr[i]; j<A->rowptr[i+1]; j++){
                t += (A->values[j]*x[A->colptr[j]]);
            }
            y[i] += t;
        }

#ifdef _OPENMP
#pragma omp parallel for private (i, j, t, r)
#endif
        for ( i = 0; i < A->rows; i++){
            r = x[i];
            for ( j = A->rowptr[i]; j < A->rowptr[i+1]; j++){
                if ( i != A->colptr[j] ){
                    t = skew * A->values[j] * r;
#ifdef _OPENMP
#pragma omp atomic
#endif
                    y[A->colptr[j]] += t;
                }
            }
        }
    }
    return (0);
}



