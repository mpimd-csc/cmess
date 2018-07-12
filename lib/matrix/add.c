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
 * @file lib/matrix/add.c
 * @brief Add two matrices.
 * @author @koehlerm
 *
 * This file contains various functions to perform matrix-matrix additions and
 * similar problems.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>



/**
 * @internal
 * @brief Add a real matrix to a real one (Kernel for @ref MESS_CSC and @ref MESS_CSR matrices).
 * @param[in] n         input dimension of the first index pointer (row or column pointer)
 * @param[in] alpha     input factor \f$alpha\f$
 * @param[in] Asptr     input array with the start indices of the rows/columns in \f$A\f$
 * @param[in] Aiptr     input array with the column/rows indices of the values of  \f$A\f$
 * @param[in] Avalues   input array containing the values of \f$A\f$
 * @param[in] beta      input factor \f$\beta\f$
 * @param[in] Bsptr     input array with the start indices of the rows/columns in \f$B\f$
 * @param[in] Biptr     input array with the column/rows indices of the values of \f$B\f$
 * @param[in] Bvalues   input array containing the values of \f$B\f$
 * @param[out] sptr     output array with the start indices of the rows/columns in the result
 * @param[out] iptr     output array with the column/row indices of the values in the result
 * @param[out] values   output array with the values of the result
 * @return Nothing

 * The @ref add_kernel_real function adds two real matrices stored whether in @ref MESS_CSR or @ref MESS_CSC storage.
 *
 * @attention Internal use only.
 *
 */
static void add_kernel_real(mess_int_t n, double alpha, mess_int_t *Asptr, mess_int_t *Aiptr, double *Avalues,
        double beta,  mess_int_t *Bsptr, mess_int_t *Biptr, double *Bvalues,
        mess_int_t *sptr, mess_int_t *iptr, double *values)
{
    mess_int_t i, rcA,reA,rcB,reB, wpos;
    sptr[0] = 0;
    for (i = 0; i < n; i++) {
        rcA = Asptr[i];
        reA = Asptr[i+1];
        rcB = Bsptr[i];
        reB = Bsptr[i+1];
        wpos = 0;
        // union of two lines
        while ( rcA < reA && rcB < reB){
            if ( Aiptr[rcA] < Biptr[rcB]){
                values[sptr[i]+wpos] = alpha*Avalues[rcA];
                iptr  [sptr[i]+wpos] = Aiptr[rcA];
                wpos++;
                rcA++;
            } else if ( Aiptr[rcA] == Biptr[rcB]){
                values[sptr[i]+wpos] = alpha*Avalues[rcA] + beta*Bvalues[rcB];
                iptr  [sptr[i]+wpos] = Aiptr[rcA];
                wpos++;
                rcA++;
                rcB++;
            } else {
                values[sptr[i]+wpos] = beta*Bvalues[rcB];
                iptr  [sptr[i]+wpos] = Biptr[rcB];
                wpos++;
                rcB++;
            }
        }
        while (rcA<reA) {
            values[sptr[i]+wpos] = alpha*Avalues[rcA];
            iptr  [sptr[i]+wpos] = Aiptr[rcA];
            wpos++;
            rcA++;
        }
        while (rcB<reB) {
            values[sptr[i]+wpos] = beta*Bvalues[rcB];
            iptr  [sptr[i]+wpos] = Biptr[rcB];
            wpos++;
            rcB++;
        }
        sptr[i+1] = sptr[i]+wpos;
    }
}

/**
 * @internal
 * @brief Add a complex matrix to a complex one (Kernel for @ref MESS_CSC and @ref MESS_CSR matrices).
 * @param[in] n         input dimension of the first indexpointer (row or column pointer)
 * @param[in] alpha     input factor \f$\alpha\f$
 * @param[in] Asptr     input array with the start indices of the rows/columns in \f$A\f$
 * @param[in] Aiptr     input array with the column/rows indices of the values of  \f$A\f$
 * @param[in] Avalues   input array containing the values of \f$A\f$
 * @param[in] beta      input factor \f$\beta\f$
 * @param[in] Bsptr     input array with the start indices of the rows/columns in \f$B\f$
 * @param[in] Biptr     input array with the column/rows indices of the values of \f$B\f$
 * @param[in] Bvalues   input array containing the values of \f$B\f$
 * @param[out] sptr     output array with the start indices of the rows/columns in the result
 * @param[out] iptr     output array with the column/row indices of the values in the result
 * @param[out] values   output array with the values of the result
 * @return Nothing
 *
 * The @ref add_kernel_complex function adds two complex matrices stored whether in @ref MESS_CSR or @ref MESS_CSC storage.
 * @attention Internal use only.
 *
 */
static void add_kernel_complex(mess_int_t n, mess_double_cpx_t alpha, mess_int_t *Asptr, mess_int_t *Aiptr, mess_double_cpx_t *Avalues,
        mess_double_cpx_t beta,  mess_int_t *Bsptr, mess_int_t *Biptr, mess_double_cpx_t *Bvalues,
        mess_int_t *sptr, mess_int_t *iptr, mess_double_cpx_t *values)
{
    mess_int_t i, rcA,reA,rcB,reB, wpos;
    sptr[0] = 0;
    for (i = 0; i < n; i++) {
        rcA = Asptr[i];
        reA = Asptr[i+1];
        rcB = Bsptr[i];
        reB = Bsptr[i+1];
        wpos = 0;
        // union of two lines
        while ( rcA < reA && rcB < reB){
            if ( Aiptr[rcA] < Biptr[rcB]){
                values[sptr[i]+wpos] = alpha*Avalues[rcA];
                iptr  [sptr[i]+wpos] = Aiptr[rcA];
                wpos++;
                rcA++;
            } else if ( Aiptr[rcA] == Biptr[rcB]){
                values[sptr[i]+wpos] = alpha*Avalues[rcA] + beta*Bvalues[rcB];
                iptr  [sptr[i]+wpos] = Aiptr[rcA];
                wpos++;
                rcA++;
                rcB++;
            } else {
                values[sptr[i]+wpos] = beta*Bvalues[rcB];
                iptr  [sptr[i]+wpos] = Biptr[rcB];
                wpos++;
                rcB++;
            }
        }
        while (rcA<reA) {
            values[sptr[i]+wpos] = alpha*Avalues[rcA];
            iptr  [sptr[i]+wpos] = Aiptr[rcA];
            wpos++;
            rcA++;
        }
        while (rcB<reB) {
            values[sptr[i]+wpos] = beta*Bvalues[rcB];
            iptr  [sptr[i]+wpos] = Biptr[rcB];
            wpos++;
            rcB++;
        }
        sptr[i+1] = sptr[i]+wpos;
    }
}

/**
 * @internal
 * @brief Add a complex matrix to a real one (Kernel for @ref MESS_CSC and @ref MESS_CSR matrices).
 * @param[in] n         input dimension of the first index pointer (row or column pointer)
 * @param[in] alpha     input factor \f$\alpha\f$
 * @param[in] Asptr     input array with the start indices of the rows/columns in \f$A\f$
 * @param[in] Aiptr     input array with the column/rows indices of the values of  \f$A\f$
 * @param[in] Avalues   input array containing the values of \f$A\f$
 * @param[in] beta      input factor \f$\beta\f$
 * @param[in] Bsptr     input array with the start indices of the rows/columns in \f$B\f$
 * @param[in] Biptr     input array with the column/rows indices of the values of \f$B\f$
 * @param[in] Bvalues   input array containing the values of \f$B\f$
 * @param[out] sptr     output array with the start indices of the rows/columns in the result
 * @param[out] iptr     output array with the column/row indices of the values in the result
 * @param[out] values   output array with the values of the result
 * @return Nothing

 * The @ref add_kernel_complex_real function adds a complex matrix to a real one. The result will be complex again.
 * The matrices stored whether in @ref MESS_CSR or @ref MESS_CSC storage.
 *
 * @attention Internal use only.
 *
 */
static void add_kernel_complex_real(mess_int_t n, mess_double_cpx_t alpha, mess_int_t *Asptr, mess_int_t *Aiptr, mess_double_cpx_t *Avalues,
        mess_double_cpx_t beta,  mess_int_t *Bsptr, mess_int_t *Biptr, double *Bvalues,
        mess_int_t *sptr, mess_int_t *iptr, mess_double_cpx_t *values)
{
    mess_int_t i, rcA,reA,rcB,reB, wpos;
    sptr[0] = 0;
    for (i = 0; i < n; i++) {
        rcA = Asptr[i];
        reA = Asptr[i+1];
        rcB = Bsptr[i];
        reB = Bsptr[i+1];
        wpos = 0;
        // union of two lines
        while ( rcA < reA && rcB < reB){
            if ( Aiptr[rcA] < Biptr[rcB]){
                values[sptr[i]+wpos] = alpha*Avalues[rcA];
                iptr  [sptr[i]+wpos] = Aiptr[rcA];
                wpos++;
                rcA++;
            } else if ( Aiptr[rcA] == Biptr[rcB]){
                values[sptr[i]+wpos] = alpha*Avalues[rcA] + beta*Bvalues[rcB];
                iptr  [sptr[i]+wpos] = Aiptr[rcA];
                wpos++;
                rcA++;
                rcB++;
            } else {
                values[sptr[i]+wpos] = beta*Bvalues[rcB];
                iptr  [sptr[i]+wpos] = Biptr[rcB];
                wpos++;
                rcB++;
            }
        }
        while (rcA<reA) {
            values[sptr[i]+wpos] = alpha*Avalues[rcA];
            iptr  [sptr[i]+wpos] = Aiptr[rcA];
            wpos++;
            rcA++;
        }
        while (rcB<reB) {
            values[sptr[i]+wpos] = beta*Bvalues[rcB];
            iptr  [sptr[i]+wpos] = Biptr[rcB];
            wpos++;
            rcB++;
        }
        sptr[i+1] = sptr[i]+wpos;
    }
}

/**
 * @internal
 * @brief Add a real matrix to a complex one (Kernel for @ref MESS_CSC and @ref MESS_CSR matrices).
 * @param[in] n         input dimension of the first index pointer (row or column pointer)
 * @param[in] alpha     input factor \f$\alpha\f$
 * @param[in] Asptr     input array with the start indices of the rows/columns in \f$A\f$
 * @param[in] Aiptr     input array with the column/rows indices of the values of  \f$A\f$
 * @param[in] Avalues   input array containing the values of \f$A\f$
 * @param[in] beta      input factor \f$\beta\f$
 * @param[in] Bsptr     input array with the start indices of the rows/columns in \f$B\f$
 * @param[in] Biptr     input array with the column/rows indices of the values of \f$B\f$
 * @param[in] Bvalues   input array containing the values of \f$B\f$
 * @param[out] sptr     output array with the start indices of the rows/columns in the result
 * @param[out] iptr     output array with the column/row indices of the values in the result
 * @param[out] values   output array with the values of the result
 * @return Nothing

 * The @ref add_kernel_real_complex function adds a real matrix to a complex one. The result will be complex again.
 * The matrices stored whether in @ref MESS_CSR or @ref MESS_CSC storage.
 * @attention Internal use only.
 *
 */
static void add_kernel_real_complex(mess_int_t n, mess_double_cpx_t alpha, mess_int_t *Asptr, mess_int_t *Aiptr, double *Avalues,
        mess_double_cpx_t beta,  mess_int_t *Bsptr, mess_int_t *Biptr, mess_double_cpx_t *Bvalues,
        mess_int_t *sptr, mess_int_t *iptr, mess_double_cpx_t *values)
{
    mess_int_t i, rcA,reA,rcB,reB, wpos;
    sptr[0] = 0;
    for (i = 0; i < n; i++) {
        rcA = Asptr[i];
        reA = Asptr[i+1];
        rcB = Bsptr[i];
        reB = Bsptr[i+1];
        wpos = 0;
        // union of two lines
        while ( rcA < reA && rcB < reB){
            if ( Aiptr[rcA] < Biptr[rcB]){
                values[sptr[i]+wpos] = alpha*Avalues[rcA];
                iptr  [sptr[i]+wpos] = Aiptr[rcA];
                wpos++;
                rcA++;
            } else if ( Aiptr[rcA] == Biptr[rcB]){
                values[sptr[i]+wpos] = alpha*Avalues[rcA] + beta*Bvalues[rcB];
                iptr  [sptr[i]+wpos] = Aiptr[rcA];
                wpos++;
                rcA++;
                rcB++;
            } else {
                values[sptr[i]+wpos] = beta*Bvalues[rcB];
                iptr  [sptr[i]+wpos] = Biptr[rcB];
                wpos++;
                rcB++;
            }
        }
        while (rcA<reA) {
            values[sptr[i]+wpos] = alpha*Avalues[rcA];
            iptr  [sptr[i]+wpos] = Aiptr[rcA];
            wpos++;
            rcA++;
        }
        while (rcB<reB) {
            values[sptr[i]+wpos] = beta*Bvalues[rcB];
            iptr  [sptr[i]+wpos] = Biptr[rcB];
            wpos++;
            rcB++;
        }
        sptr[i+1] = sptr[i]+wpos;
    }
}


/**
 * @internal
 * @brief Add two real matrices (kernel for @ref MESS_CSR matrices).
 * @param[in] alpha     input coefficient in front of \f$A\f$
 * @param[in] A         input matrix \f$A\f$
 * @param[in] beta      input coefficient in front of \f$B\f$
 * @param[in,out] B     input/output matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __mess_matrix_add_bothsparse_csr function add \f$ \alpha \f$ times a matrix \f$A\f$ to \f$ \beta \f$ times a matrix
 * \f$B\f$ and store it in \f$B\f$. This is only the kernel for @ref MESS_CSR matrices which is wrapped by \ref mess_matrix_addc.
 * \f[
 * B \leftarrow \alpha A + \beta B
 * \f]
 * @attention Internal use only.
 */
static int __mess_matrix_add_bothsparse_csr(double alpha, mess_matrix A, double beta, mess_matrix B){
    MSG_FNAME(__func__);
    mess_matrix wA, wB;
    int convA = -1, convB = -1;
    mess_int_t nnz = 0;
    mess_int_t *rowptr;
    mess_int_t *colptr;
    double *values;

    MESS_MATRIX_CHECKFORMAT(A, wA, convA, MESS_CSR);
    MESS_MATRIX_CHECKFORMAT(B, wB, convB, MESS_CSR);

    mess_try_alloc(rowptr, mess_int_t *,sizeof(mess_int_t)*(wA->rows+1));
    mess_try_alloc(colptr, mess_int_t *,sizeof(mess_int_t)*(wA->nnz+wB->nnz));
    mess_try_alloc(values, double *,sizeof(double)*(wA->nnz+wB->nnz));

    add_kernel_real(wA->rows, alpha, wA->rowptr, wA->colptr, wA->values, beta,  wB->rowptr, wB->colptr, wB->values, rowptr, colptr, values);
    nnz = rowptr[wA->rows];
    mess_try_realloc(values,double *, nnz * sizeof ( double ));
    mess_try_realloc(colptr,mess_int_t *, nnz * sizeof ( mess_int_t ));

    if ( B->colptr != NULL) mess_free(B->colptr);
    if ( B->rowptr != NULL) mess_free(B->rowptr);
    if ( B->values != NULL ) mess_free(B->values);
    if ( B->values_cpx != NULL) mess_free(B->values_cpx);
    B->colptr = colptr;
    B->rowptr = rowptr;
    B->values = values;
    B->values_cpx = NULL;
    B->nnz = nnz;
    B->data_type = MESS_REAL;
    B->store_type = MESS_CSR;
    B->symmetry = MESS_GENERAL;

    if ( convA == 0 ){ mess_matrix_clear(&wA);  }
    if ( convB == 0 ){ mess_matrix_clear(&wB);  }
    return(0);
}

/**
 * @internal
 * @brief Add two real matrices (kernel for @ref MESS_CSC matrices).
 * @param[in] alpha     input coefficient in front of \f$A\f$
 * @param[in] A         input matrix \f$A\f$
 * @param[in] beta      input coefficient in front of \f$B\f$
 * @param[in,out] B     input/output matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __mess_matrix_add_bothsparse_csc function add \f$\alpha\f$ times a matrix \f$A\f$ to \f$\beta\f$ times a matrix
 * \f$B\f$ and store it in \f$B\f$. This is only the kernel for @ref MESS_CSC matrices which is wrapped by \ref mess_matrix_addc.
 * \f[
 *  B \leftarrow \alpha A + \beta B
 * \f]
 * @attention Internal use only.
 */
static int __mess_matrix_add_bothsparse_csc(double alpha, mess_matrix A, double beta, mess_matrix B){
    MSG_FNAME(__func__);
    mess_matrix wA, wB;
    int convA = -1, convB = -1;
    mess_int_t nnz = 0;
    mess_int_t *rowptr;
    mess_int_t *colptr;
    double *values;

    MESS_MATRIX_CHECKFORMAT(A, wA, convA, MESS_CSC);
    MESS_MATRIX_CHECKFORMAT(B, wB, convB, MESS_CSC);

    mess_try_alloc(colptr, mess_int_t *,sizeof(mess_int_t)*(wA->cols+1));
    mess_try_alloc(rowptr, mess_int_t *,sizeof(mess_int_t)*(wA->nnz+wB->nnz));
    mess_try_alloc(values, double *,sizeof(double)*(wA->nnz+wB->nnz));

    add_kernel_real(wA->cols , alpha, wA->colptr, wA->rowptr, wA->values,
            beta,  wB->colptr, wB->rowptr, wB->values,
            colptr, rowptr, values);
    nnz = colptr[wA->cols];
    mess_try_realloc(values,double *, nnz * sizeof ( double ));
    mess_try_realloc(rowptr,mess_int_t *, nnz * sizeof ( mess_int_t ));

    if ( B->colptr != NULL) mess_free(B->colptr);
    if ( B->rowptr != NULL) mess_free(B->rowptr);
    if ( B->values != NULL ) mess_free(B->values);
    if ( B->values_cpx != NULL) mess_free(B->values_cpx);
    B->colptr = colptr;
    B->rowptr = rowptr;
    B->values = values;
    B->values_cpx = NULL;
    B->nnz = nnz;
    B->data_type = MESS_REAL;
    B->store_type = MESS_CSC;
    B->symmetry = MESS_GENERAL;

    if ( convA == 0 ){ mess_matrix_clear(&wA);  }
    if ( convB == 0 ){ mess_matrix_clear(&wB);  }
    return(0);
}

/**
 * @internal
 * @brief Add two matrices (kernel for @ref MESS_CSR matrices), both complex version.
 * @param[in] alpha     input coefficient in front of \f$A\f$
 * @param[in] A         input matrix \f$A\f$
 * @param[in] beta      input coefficient in front of \f$B\f$
 * @param[in,out] B     input/output matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref __mess_matrix_add_bothsparse_complex_csr function add \f$\alpha\f$ times a matrix \f$A\f$ to \f$\beta\f$ times a matrix
 * \f$B\f$ and store it in \f$B\f$. This is only the kernel for @ref MESS_CSR matrices which is wrapped by @ref mess_matrix_add.
 * \f[
 *  B \leftarrow \alpha A + \beta B
 * \f]
 * This function assumes that \f$A\f$ and \f$B\f$ are complex matrices.
 *
 * @attention Internal use only.
 */
static int __mess_matrix_add_bothsparse_complex_csr(mess_double_cpx_t alpha, mess_matrix A, mess_double_cpx_t beta, mess_matrix B){
    MSG_FNAME(__func__);
    mess_matrix wA, wB;
    int convA = -1, convB = -1;
    mess_int_t nnz = 0;
    mess_int_t *rowptr;
    mess_int_t *colptr;
    mess_double_cpx_t *values;

    MESS_MATRIX_CHECKFORMAT(A, wA, convA, MESS_CSR);
    MESS_MATRIX_CHECKFORMAT(B, wB, convB, MESS_CSR);

    mess_try_alloc(rowptr, mess_int_t *,sizeof(mess_int_t)*(wA->rows+1))
        mess_try_alloc(colptr, mess_int_t *,sizeof(mess_int_t)*(wA->nnz+wB->nnz));
    mess_try_alloc(values, mess_double_cpx_t *,sizeof(mess_double_cpx_t)*(wA->nnz+wB->nnz));

    add_kernel_complex(wA->rows, alpha, wA->rowptr, wA->colptr, wA->values_cpx,
            beta,  wB->rowptr, wB->colptr, wB->values_cpx,
            rowptr, colptr, values);

    nnz = rowptr[wA->rows];
    mess_try_realloc(values, mess_double_cpx_t *,nnz * sizeof ( mess_double_cpx_t ));
    mess_try_realloc(colptr, mess_int_t* ,nnz * sizeof ( mess_int_t ));

    if ( B->colptr != NULL) mess_free(B->colptr);
    if ( B->rowptr != NULL) mess_free(B->rowptr);
    if ( B->values != NULL ) mess_free(B->values);
    if ( B->values_cpx != NULL) mess_free(B->values_cpx);
    B->rowptr = rowptr;
    B->colptr = colptr;
    B->values = NULL;
    B->values_cpx = values;
    B->nnz = nnz;
    B->data_type = MESS_COMPLEX;
    B->store_type = MESS_CSR;
    B->symmetry = MESS_GENERAL;

    if ( convA == 0 ){ mess_matrix_clear(&wA);  }
    if ( convB == 0 ){ mess_matrix_clear(&wB);  }
    return(0);
}

/**
 * @internal
 * @brief Add two matrices (kernel for @ref MESS_CSC matrices), both complex version.
 * @param[in] alpha     input coefficient in front of \f$A\f$
 * @param[in] A         input matrix \f$A\f$
 * @param[in]  beta     input coefficient in front of \f$B\f$
 * @param[in,out] B     input/output matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref __mess_matrix_add_bothsparse_complex_csc function add \f$alpha\f$ times a matrix \f$A\f$ to \f$\beta\f$ times a matrix
 * \f$B\f$ and store it in \f$B\f$. This is only the kernel for @ref MESS_CSC matrices which is wrapped by @ref mess_matrix_add.
 * \f[
 *  B \leftarrow \alpha A + \beta B
 * \f]
 * This function assumes that \f$A\f$ and \f$B\f$ are complex matrices.
 *
 * @attention Internal use only.
 */
static int __mess_matrix_add_bothsparse_complex_csc(mess_double_cpx_t alpha, mess_matrix A, mess_double_cpx_t beta, mess_matrix B){
    MSG_FNAME(__func__);
    mess_matrix wA, wB;
    int convA = -1, convB = -1;
    mess_int_t nnz = 0;
    mess_int_t *rowptr;
    mess_int_t *colptr;
    mess_double_cpx_t *values;

    MESS_MATRIX_CHECKFORMAT(A, wA, convA, MESS_CSC);
    MESS_MATRIX_CHECKFORMAT(B, wB, convB, MESS_CSC);

    mess_try_alloc(colptr, mess_int_t *,sizeof(mess_int_t)*(wA->cols+1))
    mess_try_alloc(rowptr, mess_int_t *,sizeof(mess_int_t)*(wA->nnz+wB->nnz));
    mess_try_alloc(values, mess_double_cpx_t *,sizeof(mess_double_cpx_t)*(wA->nnz+wB->nnz));

    add_kernel_complex(wA->cols, alpha, wA->colptr, wA->rowptr, wA->values_cpx,
            beta,  wB->colptr, wB->rowptr, wB->values_cpx,
            colptr, rowptr, values);

    nnz = colptr[wA->cols];
    mess_try_realloc(values, mess_double_cpx_t *,nnz * sizeof ( mess_double_cpx_t ));
    mess_try_realloc(rowptr, mess_int_t* ,nnz * sizeof ( mess_int_t ));

    if ( B->colptr != NULL) mess_free(B->colptr);
    if ( B->rowptr != NULL) mess_free(B->rowptr);
    if ( B->values != NULL ) mess_free(B->values);
    if ( B->values_cpx != NULL) mess_free(B->values_cpx);
    B->rowptr = rowptr;
    B->colptr = colptr;
    B->values = NULL;
    B->values_cpx = values;
    B->nnz = nnz;
    B->data_type = MESS_COMPLEX;
    B->store_type = MESS_CSC;
    B->symmetry = MESS_GENERAL;

    if ( convA == 0 ){ mess_matrix_clear(&wA);  }
    if ( convB == 0 ){ mess_matrix_clear(&wB);  }
    return(0);
}

/**
 * @internal
 * @brief Add two matrices (kernel for @ref MESS_CSR matrices), real - complex version.
 * @param[in] alpha     input coefficient in front of \f$A\f$
 * @param[in] A         input matrix \f$A\f$
 * @param[in] beta      input coefficient in front of \f$B\f$
 * @param[in,out] B     input/output matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref __mess_matrix_add_bothsparse_real_complex_csr function add \f$alpha\f$ times a matrix \f$A\f$ to \f$\beta\f$ times a matrix
 * \f$B\f$ and store it in \f$B\f$. This is only the kernel for @ref MESS_CSR matrices which is wrapped by @ref mess_matrix_add.
 * \f[
 *  B \leftarrow \alpha A + \beta B
 * \f]
 * This function assumes that \f$A\f$ is real and \f$B\f$ is complex.
 *
 * @attention Internal use only.
 */
static int __mess_matrix_add_bothsparse_real_complex_csr(mess_double_cpx_t alpha, mess_matrix A, mess_double_cpx_t beta, mess_matrix B){
    MSG_FNAME(__func__);
    mess_matrix wA, wB;
    int convA = -1, convB = -1;
    mess_int_t nnz = 0;
    mess_int_t *rowptr;
    mess_int_t *colptr;
    mess_double_cpx_t *values;

    MESS_MATRIX_CHECKFORMAT(A, wA, convA, MESS_CSR);
    MESS_MATRIX_CHECKFORMAT(B, wB, convB, MESS_CSR);

    mess_try_alloc(rowptr, mess_int_t *,sizeof(mess_int_t)*(wA->rows+1))
    mess_try_alloc(colptr, mess_int_t *,sizeof(mess_int_t)*(wA->nnz+wB->nnz));
    mess_try_alloc(values, mess_double_cpx_t *,sizeof(mess_double_cpx_t)*(wA->nnz+wB->nnz));

    add_kernel_real_complex(wA->rows, alpha, wA->rowptr, wA->colptr, wA->values,
            beta,  wB->rowptr, wB->colptr, wB->values_cpx,
            rowptr, colptr, values);

    nnz = rowptr[wA->rows];
    mess_try_realloc(values, mess_double_cpx_t *,nnz * sizeof ( mess_double_cpx_t ));
    mess_try_realloc(colptr, mess_int_t* ,nnz * sizeof ( mess_int_t ));

    if ( B->colptr != NULL) mess_free(B->colptr);
    B->colptr = colptr;
    if ( B->rowptr != NULL) mess_free(B->rowptr);
    B->rowptr = rowptr;
    if ( B->values != NULL ) mess_free(B->values);
    if ( B->values_cpx != NULL) mess_free(B->values_cpx);
    B->values = NULL;
    B->values_cpx = values;
    B->nnz = nnz;
    B->data_type = MESS_COMPLEX;
    B->store_type = MESS_CSR;
    B->symmetry = MESS_GENERAL;

    if ( convA == 0 ){ mess_matrix_clear(&wA); }
    if ( convB == 0 ){ mess_matrix_clear(&wB); }
    return(0);
}

/**
 * @internal
 * @brief Add two matrices (kernel for @ref MESS_CSC matrices), real - complex version.
 * @param[in] alpha     input coefficient in front of \f$A\f$
 * @param[in] A         input matrix \f$A\f$
 * @param[in] beta      input coefficient in front of \f$B\f$
 * @param[in,out] B     input/output matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref __mess_matrix_add_bothsparse_real_complex_csc function add \f$alpha\f$ times a matrix \f$A\f$ to \f$\beta\f$ times a matrix
 * \f$B\f$ and store it in \f$B\f$. This is only the kernel for @ref MESS_CSC matrices which is wrapped by @ref mess_matrix_add.
 * \f[
 * B \leftarrow \alpha A + \beta B
 * \f]
 * This function assumes that \f$A\f$ is real and \f$B\f$ is complex.
 *
 * @attention Internal use only.
 */
static int __mess_matrix_add_bothsparse_real_complex_csc(mess_double_cpx_t alpha, mess_matrix A, mess_double_cpx_t beta, mess_matrix B){
    MSG_FNAME(__func__);
    mess_matrix wA, wB;
    int convA = -1, convB = -1;
    mess_int_t nnz = 0;
    mess_int_t *rowptr;
    mess_int_t *colptr;
    mess_double_cpx_t *values;

    MESS_MATRIX_CHECKFORMAT(A, wA, convA, MESS_CSC);
    MESS_MATRIX_CHECKFORMAT(B, wB, convB, MESS_CSC);

    mess_try_alloc(colptr, mess_int_t *,sizeof(mess_int_t)*(wA->cols+1))
    mess_try_alloc(rowptr, mess_int_t *,sizeof(mess_int_t)*(wA->nnz+wB->nnz));
    mess_try_alloc(values, mess_double_cpx_t *,sizeof(mess_double_cpx_t)*(wA->nnz+wB->nnz));

    add_kernel_real_complex(wA->cols, alpha, wA->colptr, wA->rowptr, wA->values,
            beta,  wB->colptr, wB->rowptr, wB->values_cpx,
            colptr,rowptr, values);

    nnz = colptr[wA->cols];
    mess_try_realloc(values, mess_double_cpx_t *,nnz * sizeof ( mess_double_cpx_t ));
    mess_try_realloc(rowptr, mess_int_t* ,nnz * sizeof ( mess_int_t ));

    if ( B->colptr != NULL) mess_free(B->colptr);
    B->colptr = colptr;
    if ( B->rowptr != NULL) mess_free(B->rowptr);
    B->rowptr = rowptr;
    if ( B->values != NULL ) mess_free(B->values);
    if ( B->values_cpx != NULL) mess_free(B->values_cpx);
    B->values = NULL;
    B->values_cpx = values;
    B->nnz = nnz;
    B->data_type = MESS_COMPLEX;
    B->store_type = MESS_CSC;
    B->symmetry = MESS_GENERAL;

    if ( convA == 0 ){ mess_matrix_clear(&wA); }
    if ( convB == 0 ){ mess_matrix_clear(&wB); }
    return(0);
}


/**
 * @internal
 * @brief Add two matrices (kernel for @ref MESS_CSR matrices), complex-real version.
 * @param[in] alpha     input coefficient in front of \f$A\f$
 * @param[in] A         input matrix \f$A\f$
 * @param[in] beta      input coefficient in front of \f$B\f$
 * @param[in,out] B     input/output matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref __mess_matrix_add_bothsparse_complex_real_csr function add \f$alpha\f$ times a matrix \f$A\f$ to \f$\beta\f$ times a matrix
 * \f$B\f$ and store it in \f$B\f$. This is only the kernel for @ref MESS_CSR matrices which is wrapped by @ref mess_matrix_add.
 * \f[
 * B \leftarrow \alpha A + \beta B
 * \f]
 * This function assumes that \f$A\f$ is complex and \f$B\f$ is real on input.
 *
 * @attention Internal use only.
 */
static int __mess_matrix_add_bothsparse_complex_real_csr(mess_double_cpx_t alpha, mess_matrix A, mess_double_cpx_t beta, mess_matrix B){
    MSG_FNAME(__func__);
    mess_matrix wA, wB;
    int convA = -1, convB = -1;
    mess_int_t nnz = 0;
    mess_int_t *rowptr;
    mess_int_t *colptr;
    mess_double_cpx_t *values;

    MESS_MATRIX_CHECKFORMAT(A, wA, convA, MESS_CSR);
    MESS_MATRIX_CHECKFORMAT(B, wB, convB, MESS_CSR);

    mess_try_alloc(rowptr, mess_int_t *,sizeof(mess_int_t)*(wA->rows+1))
        mess_try_alloc(colptr, mess_int_t *,sizeof(mess_int_t)*(wA->nnz+wB->nnz));
    mess_try_alloc(values, mess_double_cpx_t *,sizeof(mess_double_cpx_t)*(wA->nnz+wB->nnz));

    add_kernel_complex_real(wA->rows, alpha, wA->rowptr, wA->colptr, wA->values_cpx,
            beta,  wB->rowptr, wB->colptr, wB->values,
            rowptr, colptr, values);

    nnz = rowptr[wA->rows];
    mess_try_realloc(values, mess_double_cpx_t *,nnz * sizeof ( mess_double_cpx_t ));
    mess_try_realloc(colptr, mess_int_t* ,nnz * sizeof ( mess_int_t ));

    if ( B->colptr != NULL) mess_free(B->colptr);
    B->colptr = colptr;
    if ( B->rowptr != NULL) mess_free(B->rowptr);
    B->rowptr = rowptr;
    if ( B->values != NULL ) mess_free(B->values);
    if ( B->values_cpx != NULL) mess_free(B->values_cpx);
    B->values = NULL;
    B->values_cpx = values;
    B->nnz = nnz;
    B->data_type = MESS_COMPLEX;
    B->store_type = MESS_CSR;
    B->symmetry = MESS_GENERAL;

    if ( convA == 0 ){ mess_matrix_clear(&wA); }
    if ( convB == 0 ){ mess_matrix_clear(&wB); }
    return(0);
}


/**
 * @internal
 * @brief Add two matrices (kernel for @ref MESS_CSC matrices), complex-real version.
 * @param[in] alpha     input coefficient in front of \f$A\f$
 * @param[in] A         input matrix \f$A\f$
 * @param[in] beta      input coefficient in front of \f$B\f$
 * @param[in,out] B     input/output matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref __mess_matrix_add_bothsparse_complex_real_csc function add \f$alpha\f$ times a matrix \f$A\f$ to \f$\beta\f$ times a matrix
 * \f$B\f$ and store it in \f$B\f$. This is only the kernel for @ref MESS_CSC matrices which is wrapped by @ref mess_matrix_add.
 * \f[
 *  B \leftarrow \alpha A + \beta B
 * \f]
 * This function assumes that \f$A\f$ is complex and \f$B\f$ is real on input.
 *
 * @attention Internal use only.
 */
static int __mess_matrix_add_bothsparse_complex_real_csc(mess_double_cpx_t alpha, mess_matrix A, mess_double_cpx_t beta, mess_matrix B){
    MSG_FNAME(__func__);
    mess_matrix wA, wB;
    int convA = -1, convB = -1;
    mess_int_t nnz = 0;
    mess_int_t *rowptr;
    mess_int_t *colptr;
    mess_double_cpx_t *values;

    MESS_MATRIX_CHECKFORMAT(A, wA, convA, MESS_CSC);
    MESS_MATRIX_CHECKFORMAT(B, wB, convB, MESS_CSC);

    mess_try_alloc(colptr, mess_int_t *,sizeof(mess_int_t)*(wA->cols+1))
    mess_try_alloc(rowptr, mess_int_t *,sizeof(mess_int_t)*(wA->nnz+wB->nnz));
    mess_try_alloc(values, mess_double_cpx_t *,sizeof(mess_double_cpx_t)*(wA->nnz+wB->nnz));

    add_kernel_complex_real(wA->cols, alpha, wA->colptr, wA->rowptr, wA->values_cpx, beta,  wB->colptr, wB->rowptr, wB->values, colptr, rowptr, values);

    nnz = colptr[wA->cols];
    mess_try_realloc(values, mess_double_cpx_t *,nnz * sizeof ( mess_double_cpx_t ));
    mess_try_realloc(rowptr, mess_int_t* ,nnz * sizeof ( mess_int_t ));

    if ( B->colptr != NULL) mess_free(B->colptr);
    B->colptr = colptr;
    if ( B->rowptr != NULL) mess_free(B->rowptr);
    B->rowptr = rowptr;
    if ( B->values != NULL ) mess_free(B->values);
    if ( B->values_cpx != NULL) mess_free(B->values_cpx);
    B->values = NULL;
    B->values_cpx = values;
    B->nnz = nnz;
    B->data_type = MESS_COMPLEX;
    B->store_type = MESS_CSC;
    B->symmetry = MESS_GENERAL;

    if ( convA == 0 ){ mess_matrix_clear(&wA); }
    if ( convB == 0 ){ mess_matrix_clear(&wB); }
    return(0);
}



/**
 *
 * @brief Add two real matrices.
 * @param[in] alpha     input coefficient in front of \f$A\f$
 * @param[in] A         input matrix \f$A\f$
 * @param[in] beta      input coefficient in front of \f$B\f$
 * @param[in,out] B     input/output matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_add function add  \f$\alpha\f$ times a  mess_matrix  \f$A\f$ to  \f$\beta\f$ times a  matrix
 *  \f$B\f$ and store it in \f$B\f$.
 * \f[
 * B \leftarrow \alpha A + \beta B
 * \f]
 * If there are new zero elements they have to stay in the pattern for compatibility.
 * The matrix must be sorted internally. This can be achieved using \ref mess_matrix_sort.
 * This function is only a wrapper around \ref mess_matrix_addc, which provides an interface
 * for real and complex matrices.
 *
 */
int mess_matrix_add(double alpha, mess_matrix A, double beta, mess_matrix B){
    return mess_matrix_addc(alpha,A,beta,B);
}

/**
 * @brief Add two matrices.
 * @param[in] alpha     input coefficient in front of \f$A\f$
 * @param[in] A         input matrix \f$A\f$
 * @param[in] beta      input coefficient in front of \f$B\f$
 * @param[in,out] B     input/output matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_addc function add \f$\alpha\f$ times a matrix \f$A\f$ to \f$\beta\f$ times a matrix
 * \f$B\f$ and store it in \f$B\f$.
 * \f[
 * B \leftarrow \alpha A + \beta B
 * \f]
 * The matrices can be real or complex. Depending on the data, the right internal algorithm is
 * chosen. If there are new zero elements in the matrix, they will be kept due to produces the correct
 * of the sum. If the matrices are sparse, they must have sorted internal data structures. If
 * this is not the case, the result will be wrong. The internal data structure can be sorted using
 * \ref mess_matrix_sort before.
 *
 */
int mess_matrix_addc(mess_double_cpx_t alpha, mess_matrix A, mess_double_cpx_t beta, mess_matrix B){
    MSG_FNAME(__func__);
    int ret =  0;

    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_real_or_complex(A);
    mess_check_real_or_complex(B);
    mess_check_same_size(A,B);


    /*-----------------------------------------------------------------------------
     *  Detect the all real case
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_REAL(A) && MESS_IS_REAL(B)){
        double a,b;
        a = creal(alpha);
        b = creal(beta);
        if ( cimag(alpha) != 0.0 || cimag(beta) !=0.0){
            ret = mess_matrix_tocomplex(B);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_tocomplex);
            return mess_matrix_addc(alpha,A,beta,B);
            // MSG_WARN("Try to multipy a real matrix by an complex value. Ignoring imaginary part.\n");
        }
        if (MESS_IS_DENSE(A) && MESS_IS_DENSE(B)){
            F77_GLOBAL(dgeadd,DGEADD)(&A->rows,&A->cols, &a,A->values, &A->ld, &b, B->values,&B->ld);
        } else if (MESS_IS_CSR(A) && MESS_IS_CSR(B)){
            ret = __mess_matrix_add_bothsparse_csr(a, A, b, B); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_add_bothsparse_csr);
        } else if (MESS_IS_CSC(A) && MESS_IS_CSC(B)){
            ret = __mess_matrix_add_bothsparse_csc(a, A, b, B); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_add_bothsparse_csc);
        } else {
            mess_storage_t st = B->store_type;
            mess_matrix workA,workB;
            int statA,statB;
            MESS_MATRIX_CHECKFORMAT(A,workA,statA,MESS_CSR);
            MESS_MATRIX_CHECKFORMAT(B,workB,statB,MESS_CSR);
            ret  = __mess_matrix_add_bothsparse_csr (alpha, workA, beta, workB);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_add_bothsparse_csr);
            if ( B->store_type != st) {
                MSG_WARN("converting back to %s\n", mess_storage_t_str(st));
                mess_matrix tmp;
                ret = mess_matrix_init(&tmp);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(workB, tmp, st);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
                ret = mess_matrix_copy (tmp, B);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
                mess_matrix_clear(&tmp);
            }
            if ( workB != B) {
                ret = mess_matrix_copy(workB, B);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            }
            if(statA==0) mess_matrix_clear(&workA);
            if(statB==0) mess_matrix_clear(&workB);

        }
    }
    /*-----------------------------------------------------------------------------
     *  At least one value is complex
     *-----------------------------------------------------------------------------*/
    else {
        if (MESS_IS_DENSE(A) && MESS_IS_DENSE(B)){
            if ( MESS_IS_REAL(A)){
                F77_GLOBAL(dzgeadd,DZGEADD)(&A->rows,&A->cols, &alpha,A->values, &A->ld, &beta, B->values_cpx,&B->ld);
            } else {
                ret = mess_matrix_tocomplex(B);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
                F77_GLOBAL(zgeadd,ZGEADD)(&A->rows,&A->cols, &alpha,A->values_cpx, &A->ld, &beta, B->values_cpx,&B->ld);
            }
        } else if ( MESS_IS_CSR(A) && MESS_IS_CSR(B)){
            if ( MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                ret = __mess_matrix_add_bothsparse_real_complex_csr(alpha,A,beta,B);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_add_bothsparse_real_complex_csr);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                ret = __mess_matrix_add_bothsparse_complex_real_csr(alpha,A,beta,B);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_add_bothsparse_complex_real_csr);
            } else {
                ret = __mess_matrix_add_bothsparse_complex_csr(alpha,A,beta,B);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_add_bothsparse_complex_csr);
            }
        } else if ( MESS_IS_CSC(A) && MESS_IS_CSC(B)){
            if ( MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                ret = __mess_matrix_add_bothsparse_real_complex_csc(alpha,A,beta,B);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_add_bothsparse_real_complex_csc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                ret = __mess_matrix_add_bothsparse_complex_real_csc(alpha,A,beta,B);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_add_bothsparse_complex_real_csc);
            } else {
                ret = __mess_matrix_add_bothsparse_complex_csc(alpha,A,beta,B);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_add_bothsparse_complex_csc);
            }

        } else {
            mess_storage_t st = B->store_type;
            if ( MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                ret = __mess_matrix_add_bothsparse_real_complex_csr(alpha,A,beta,B);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_add_bothsparse_real_complex_csr);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                ret = __mess_matrix_add_bothsparse_complex_real_csr(alpha,A,beta,B);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_add_bothsparse_complex_real_csr);
            } else {
                ret = __mess_matrix_add_bothsparse_complex_csr(alpha,A,beta,B);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_add_bothsparse_complex_csr);
            }

            if ( B->store_type != st) {
                MSG_WARN("converting back to %s\n", mess_storage_t_str(st));
                mess_matrix tmp;
                ret = mess_matrix_init(&tmp);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B, tmp, st);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
                ret = mess_matrix_copy (tmp, B);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
                mess_matrix_clear(&tmp);
            }

        }
    }
    return(0);
}

