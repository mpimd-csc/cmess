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
 * @file lib/matrix/equal.c
 * @brief Operations on the pattern of matrices and norm of the difference of matrices.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>

/**
 * @brief Check if two matrices are equal on data structure level.
 * @param[in] mat1   input first matrix
 * @param[in] mat2   input second matrix
 * @returns zero if the matrices are not equal or an error occured, a non zero value otherwise
 *
 * The @ref mess_matrix_equal function checks if two matrices are equal on
 * data structure level.\n
 * If an error occurs this is treated as non equal matrices. \n
 * If you compare sparse matrices please ensure that they are sorted
 * using \ref mess_matrix_sort before. \n
 * If matrices have different storage type they are non equal, too. The leading dimension of dense matrices is ignored.
 *
 */
int mess_matrix_equal ( mess_matrix mat1, mess_matrix mat2 )
{
    mess_int_t i,j;

    if ( mat1 == NULL) return 0;
    if ( mat2 == NULL) return 0;
    if ( mat1->rows != mat2->rows) return 0;
    if ( mat1->cols != mat2->cols) return 0;
    if ( mat1->store_type != mat2->store_type) return 0;
    if ( mat1->symmetry   != mat2->symmetry) return 0;
    if ( mat1->data_type  != mat2->data_type) return 0;

    if ( MESS_IS_DENSE(mat1)) {
        if ( MESS_IS_REAL(mat1)) {
            for ( j = 0 ;  j < mat1->cols; j++ ) {
                for ( i =0 ; i < mat1->rows; i++ ) {
                    if ( mat1->values[i+j*mat1->ld] != mat2->values[i+j*mat2->ld]){
                        return 0;
                    }
                }
            }
        } else {
            for ( j = 0 ;  j < mat1->cols; j++ ) {
                for ( i =0 ; i < mat1->rows; i++ ) {
                    if ( mat1->values_cpx[i+j*mat1->ld] != mat2->values_cpx[i+j*mat2->ld]){
                        return 0;
                    }
                }
            }

        }
        return 1;
    } else if ( MESS_IS_CSR(mat1)) {
        if ( mat1->nnz != mat2->nnz) return 0;
        for ( i = 0 ; i <= mat1->rows ; i++) {
            if ( mat1->rowptr[i] != mat2->rowptr[i]) {
                return 0;
            }
        }
        for ( i = 0 ; i < mat1->nnz ; i++) {
            if ( mat1->colptr[i] != mat2->colptr[i]) {
                return 0;
            }
        }
        if ( MESS_IS_REAL(mat1)){
            for ( i = 0 ; i < mat1->nnz ; i++) {
                if ( mat1->values[i] != mat2->values[i]) {
                    return 0;
                }
            }
        } else {
            for ( i = 0 ; i < mat1->nnz ; i++) {
                if ( mat1->values_cpx[i] != mat2->values_cpx[i]) {
                    return 0;
                }
            }
        }
        return 1;
    } else if ( MESS_IS_CSC(mat1)) {
        if ( mat1->nnz != mat2->nnz) return 0;
        for ( i = 0 ; i < mat1->nnz ; i++) {
            if ( mat1->rowptr[i] != mat2->rowptr[i]) {
                return 0;
            }
        }
        for ( i = 0 ; i <= mat1->cols ; i++) {
            if ( mat1->colptr[i] != mat2->colptr[i]) {
                return 0;
            }
        }
        if ( MESS_IS_REAL(mat1)){
            for ( i = 0 ; i < mat1->nnz ; i++) {
                if ( mat1->values[i] != mat2->values[i]) {
                    return 0;
                }
            }
        } else {
            for ( i = 0 ; i < mat1->nnz ; i++) {
                if ( mat1->values_cpx[i] != mat2->values_cpx[i]) {
                    return 0;
                }
            }
        }
        return 1;
    } else if ( MESS_IS_COORD(mat1)) {
        if ( mat1->nnz != mat2->nnz) return 0;
        for ( i = 0 ; i < mat1->nnz ; i++) {
            if ( mat1->rowptr[i] != mat2->rowptr[i]) {
                return 0;
            }
        }
        for ( i = 0 ; i < mat1->nnz ; i++) {
            if ( mat1->colptr[i] != mat2->colptr[i]) {
                return 0;
            }
        }
        if ( MESS_IS_REAL(mat1)){
            for ( i = 0 ; i < mat1->nnz ; i++) {
                if ( mat1->values[i] != mat2->values[i]) {
                    return 0;
                }
            }
        } else {
            for ( i = 0 ; i < mat1->nnz ; i++) {
                if ( mat1->values_cpx[i] != mat2->values_cpx[i]) {
                    return 0;
                }
            }
        }
        return 1;
    } else {
        return 0;
    }
    return 1;
}



/**
 * @brief Check if two matrices are equal on data structure level (with verbose output).
 * @param[in] mat1   input first matrix
 * @param[in] mat2   input second matrix
 * @returns zero if the matrices are not equal or an error occured, a non zero value otherwise
 *
 * The @ref mess_matrix_equal_verbose function checks if two matrices are equal on
 * data structure level.\n
 * If an error occurs this is treated as non equal matrices. \n
 * If you compare sparse matrices please ensure that they are sorted
 * using \ref mess_matrix_sort before. If matrices have different storage type they
 * are non equal, too. The leading dimension of dense matrices is ignored. \n
 * In contrast to \ref mess_matrix_equal this function produces some additional output to identify where
 * the matrices differ.
 *
 */
int mess_matrix_equal_verbose ( mess_matrix mat1, mess_matrix mat2 )
{
    mess_int_t i,j;

    if ( mat1 == NULL) { printf("mess_matrix_equal: mat1 is NULL\n"); return 0;}
    if ( mat2 == NULL) { printf("mess_matrix_equal: mat2 is NULL\n"); return 0;}
    if ( mat1->rows != mat2->rows) { printf("mess_matrix_equal: number of rows are different. %d <-> %d\n", (int)mat1->rows, (int)mat2->rows); return 0;}
    if ( mat1->cols != mat2->cols) { printf("mess_matrix_equal: number of cols are different. %d <-> %d\n", (int)mat1->cols, (int)mat2->cols); return 0;}
    if ( mat1->store_type != mat2->store_type) { printf("mess_matrix_equal: storage types not equal:  %s <-> %s\n", mess_storage_t_str(mat1->store_type), mess_storage_t_str(mat2->store_type)); return 0;}
    if ( mat1->symmetry   != mat2->symmetry) { printf("mess_matrix_equal: matrices have a different symmetry.\n"); return 0;}
    if ( mat1->data_type  != mat2->data_type) { printf("mess_matrix_equal: data type not equal. %s <-> %s\n", mess_datatype_t_str(mat1->data_type), mess_datatype_t_str(mat2->data_type)); return 0;}


    if ( MESS_IS_DENSE(mat1)) {
        if ( MESS_IS_REAL(mat1)) {
            for ( j = 0 ;  j < mat1->cols; j++ ) {
                for ( i =0 ; i < mat1->rows; i++ ) {
                    if ( mat1->values[i+j*mat1->ld] != mat2->values[i+j*mat2->ld]){
                        printf("mess_matrix_equal: value(%d,%d) wrong: %lg <-> %lg\n", (int) i, (int) j, mat1->values[i+j*mat1->ld],mat2->values[i+j*mat2->ld]);
                        return 0;
                    }
                }
            }
        } else {
            for ( j = 0 ;  j < mat1->cols; j++ ) {
                for ( i =0 ; i < mat1->rows; i++ ) {
                    if ( mat1->values_cpx[i+j*mat1->ld] != mat2->values_cpx[i+j*mat2->ld]){
                        printf("mess_matrix_equal: value(%d,%d) wrong: %lg+%lg*i <-> %lg+%lg*i\n", (int)i, (int) j, creal(mat1->values[i+j*mat1->ld]),cimag(mat1->values[i+j*mat1->ld]),creal(mat2->values[i+j*mat2->ld]), cimag(mat2->values[i+j*mat2->ld]));
                        return 0;
                    }
                }
            }
        }
        return 1;
    } else if ( MESS_IS_CSR(mat1)) {
        if ( mat1->nnz != mat2->nnz) return 0;
        for ( i = 0 ; i <= mat1->rows ; i++) {
            if ( mat1->rowptr[i] != mat2->rowptr[i]) {
                printf("mess_matrix_equal: rowptr[%d] wrong: %d <-> %d\n",(int) i, (int) mat1->rowptr[i],(int) mat2->rowptr[i]);
                return 0;
            }
        }
        for ( i = 0 ; i < mat1->nnz ; i++) {
            if ( mat1->colptr[i] != mat2->colptr[i]) {
                printf("mess_matrix_equal: colptr[%d] wrong: %d <-> %d\n",(int) i, (int) mat1->rowptr[i],(int) mat2->rowptr[i]);
                return 0;
            }
        }
        if ( MESS_IS_REAL(mat1)){
            for ( i = 0 ; i < mat1->nnz ; i++) {
                if ( mat1->values[i] != mat2->values[i]) {
                    printf("mess_matrix_equal: values[%d] wrong: %lg <-> %lg\n",(int) i, mat1->values[i],mat2->values[i]);
                    return 0;
                }
            }
        } else {
            for ( i = 0 ; i < mat1->nnz ; i++) {
                if ( mat1->values_cpx[i] != mat2->values_cpx[i]) {
                    printf("mess_matrix_equal: values[%d] wrong: %lg+%lg*i <-> %lg+%lg*i\n",(int) i,  creal(mat1->values_cpx[i]),cimag(mat1->values_cpx[i]), creal(mat2->values_cpx[i]),cimag(mat2->values_cpx[i]));
                    return 0;
                }
            }
        }
        return 1;
    } else if ( MESS_IS_CSC(mat1)) {
        if ( mat1->nnz != mat2->nnz) return 0;
        for ( i = 0 ; i < mat1->nnz ; i++) {
            if ( mat1->rowptr[i] != mat2->rowptr[i]) {
                printf("mess_matrix_equal: rowptr[%d] wrong: %d <-> %d\n",(int) i, (int) mat1->rowptr[i],(int) mat2->rowptr[i]);
                return 0;
            }
        }
        for ( i = 0 ; i <= mat1->cols ; i++) {
            if ( mat1->colptr[i] != mat2->colptr[i]) {
                printf("mess_matrix_equal: colptr[%d] wrong: %d <-> %d\n",(int) i, (int) mat1->rowptr[i],(int) mat2->rowptr[i]);
                return 0;
            }
        }
        if ( MESS_IS_REAL(mat1)){
            for ( i = 0 ; i < mat1->nnz ; i++) {
                if ( mat1->values[i] != mat2->values[i]) {
                    printf("mess_matrix_equal: values[%d] wrong: %lg <-> %lg\n",(int) i, mat1->values[i],mat2->values[i]);
                    return 0;
                }
            }
        } else {
            for ( i = 0 ; i < mat1->nnz ; i++) {
                if ( mat1->values_cpx[i] != mat2->values_cpx[i]) {
                    printf("mess_matrix_equal: values[%d] wrong: %lg+%lg*i <-> %lg+%lg*i\n",(int) i,  creal(mat1->values_cpx[i]),cimag(mat1->values_cpx[i]), creal(mat2->values_cpx[i]),cimag(mat2->values_cpx[i]));
                    return 0;
                }
            }
        }
        return 1;
    } else if ( MESS_IS_COORD(mat1)) {
        if ( mat1->nnz != mat2->nnz) return 0;
        for ( i = 0 ; i < mat1->nnz ; i++) {
            if ( mat1->rowptr[i] != mat2->rowptr[i]) {
                printf("mess_matrix_equal: rowptr[%d] wrong: %d <-> %d\n",(int) i, (int) mat1->rowptr[i],(int) mat2->rowptr[i]);
                return 0;
            }
        }
        for ( i = 0 ; i < mat1->nnz ; i++) {
            if ( mat1->colptr[i] != mat2->colptr[i]) {
                printf("mess_matrix_equal: colptr[%d] wrong: %d <-> %d\n",(int) i, (int) mat1->rowptr[i],(int) mat2->rowptr[i]);
                return 0;
            }
        }
        if ( MESS_IS_REAL(mat1)){
            for ( i = 0 ; i < mat1->nnz ; i++) {
                if ( mat1->values[i] != mat2->values[i]) {
                    printf("mess_matrix_equal: values[%d] wrong: %lg <-> %lg\n",(int) i, mat1->values[i],mat2->values[i]);
                    return 0;
                }
            }
        } else {
            for ( i = 0 ; i < mat1->nnz ; i++) {
                if ( mat1->values_cpx[i] != mat2->values_cpx[i]) {
                    printf("mess_matrix_equal: values[%d] wrong: %lg+%lg*i <-> %lg+%lg*i\n",(int) i,  creal(mat1->values_cpx[i]),cimag(mat1->values_cpx[i]), creal(mat2->values_cpx[i]),cimag(mat2->values_cpx[i]));
                    return 0;
                }
            }
        }
        return 1;
    } else {
        return 0;
    }
    return 1;
}



/**
 * @brief Check if two matrices have the same pattern.
 * @param[in] A input matrix \f$A\f$
 * @param[in] B input matrix \f$B\f$
 * @return 1 if the pattern are equal, 0 if the pattern differs and all values larger than
 * one are error codes.
 *
 * The @ref mess_matrix_equalpattern function checks if the matrices \f$ A \f$ and \f$ B \f$
 * have the same pattern. \n
 * If they have the same pattern @c 1 will be returned. \n
 * If the pattern is not the same @c 0 will be returned. \n
 * All values bigger than @c 1 are error codes.
 *
 * The input matrices \f$A\f$ and \f$B\f$  have to be sorted internally. That
 * means the colptr of a @ref MESS_CSR or the rowptr of a @ref MESS_CSC \ref mess_matrix  must be ordered in the same way.
 * If the matrices \f$A\f$ or \f$B\f$  are sparse they need to be sorted using \ref mess_matrix_sort
 * before. If both are sparse the need to have the same storage format.
 *
 */
int mess_matrix_equalpattern(mess_matrix A, mess_matrix B) {
    MSG_FNAME(__func__);
    mess_int_t i = 0;
    mess_int_t eqn = 0;

    // if both are the matrix, return 1 even if they are NULL
    if ( A == B ) return(1);
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);

    if ( A->rows != B->rows ) return(0);
    if ( A->cols != B->cols ) return(0);
    if ( A->nnz != B->nnz) return(0);

    if ( MESS_IS_CSR(A) && MESS_IS_CSR(B)) {
#ifdef _OPENMP
#pragma omp parallel for private (i) reduction(+:eqn)
#endif
        for ( i = 0;  i < A->rows+1; i++ ) {
            if ( A->rowptr[i] != B->rowptr[i]) eqn++;
        }
        if ( eqn != 0) return(0);

#ifdef _OPENMP
#pragma omp parallel for private (i) reduction(+:eqn)
#endif
        for ( i = 0;  i < A->nnz; i++ ) {
            if ( A->colptr[i] != B->colptr[i]) eqn++;
        }
        if ( eqn != 0 ) return(0);
        return(1);
    } else if (MESS_IS_CSC(A) && MESS_IS_CSC(B)) {
#ifdef _OPENMP
#pragma omp parallel for private (i) reduction(+:eqn)
#endif
        for ( i = 0;  i < A->cols+1; i++ ) {
            if ( A->colptr[i] != B->colptr[i]) eqn++;
        }
        if ( eqn != 0 ) return(0);
#ifdef _OPENMP
#pragma omp parallel for private (i) reduction(+:eqn)
#endif
        for ( i = 0;  i < A->nnz; i++ ) {
            if ( A->rowptr[i] != B->rowptr[i]) eqn++;
        }
        if ( eqn != 0) return(0);
        return(1);
    } else if (MESS_IS_DENSE(A) && MESS_IS_DENSE(B)){
        // Dense matrices have the same pattern by definition.
        return(1);
    } else  {
        MSG_ERROR("storage format not supported: A=%s, B=%s\n", mess_storage_t_str(A->store_type), mess_storage_t_str(B->store_type));
        return( MESS_ERROR_NOSUPPORT);
    }
    return(0);
}


/**
 * @brief Join the pattern of two sparse matrices.
 * @param[in,out] A input/output matrix \f$A\f$
 * @param[in,out] B input/output matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_joinpatterns function joins the pattern of two sparse matrices \f$ A \f$ and \f$ B \f$. At exit
 * \f$ \mathcal{P}(A) = \mathcal{P}(B) \f$. This is done by computing
 * \f[ A \leftarrow A + 0\cdot B \f]
 * and
 * \f[ B \leftarrow B + 0\cdot A.\f]
 * using \ref mess_matrix_add.
 *
 */
int mess_matrix_joinpatterns(mess_matrix A, mess_matrix B)
{
    MSG_FNAME(__func__);
    int ret = 0 ;

    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_same_size(A,B);
    mess_check_same_datatype(A,B);
    mess_check_same_storetype(A,B);

    if ( MESS_IS_COMPLEX(A)) {
        ret = mess_matrix_addc(0.0,A, 1, B);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
        ret = mess_matrix_addc(0,B, 1, A);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
    } else if ( MESS_IS_REAL(A)) {
        ret = mess_matrix_add(0.0,A, 1, B);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_add(0,B, 1, A);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    } else {
        MSG_ERROR("Unknown / unsupported data type: %s \n",mess_datatype_t_str(A->data_type));
        return (MESS_ERROR_DATATYPE);
    }
    ret = mess_matrix_sort(A);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
    ret = mess_matrix_sort(B);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);

    return (0);
}



/**
 * @internal
 * @brief Helper structure for @ref mat_diff_mvp.
 *
 * The @ref mat_diff is a helper structure for holding
 * two \@ref mess_matrix structures \f$A\f$ and \f$B\f$.

 * @attention Internal use only.
 */
struct mat_diff {
    mess_matrix A;  /**< First matrix \f$A\f$ */
    mess_matrix B;  /**< Second matrix \f$B\f$ */
};


/**
 * @internal
 * @brief Implement the \f$ y= (A-B)x \f$  matrix-vector product for the Arnoldi process.
 * @param[in] data  input pointer to the matrix data
 * @param[in] op    input (not used)
 * @param[in] x     input right hand side vector \f$x\f$
 * @param[out] y   output vector  \f$y\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mat_diff_mvp function computes
 * \f[y=(A-B)x.\f]
 *
 *
 * @attention Internal use only.
 *
 */
static int mat_diff_mvp(void *data, mess_operation_t op,  mess_vector x, mess_vector y) {
    MSG_FNAME(__func__);
    struct mat_diff * d = (struct mat_diff*) data;
    int ret = 0;
    ret = mess_vector_zeros(y);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_zeros);
    ret = mess_matrix_gaxpy(MESS_OP_NONE,d->A, x, y);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_vector_scale(-1, y);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_scale);
    ret = mess_matrix_gaxpy(MESS_OP_NONE,d->B, x, y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    return (0);
}

/**
 * @brief Compute the \f$ 2 \f$-norm of the difference of two matrices.
 * @param[in] A   input matrix \f$A\f$
 * @param[in] B     input matrix \f$B\f$
 * @param[out] nrm  output \f$ 2 \f$-norm of the difference of matrix \f$A\f$ and matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_diffnorm function computes the \f$ 2 \f$-norm of the difference of two matrices
 * \f[ nrm=\vert\vert  A-B\vert\vert_2. \f]
 * It uses an Arnoldi process to estimate the largest eigenvalue. If both matrices
 * are @ref MESS_DENSE the difference matrix is computed and a normal eigensolver is used.
 *
 */
int mess_matrix_diffnorm ( mess_matrix A, mess_matrix B, double *nrm )
{
    MSG_FNAME(__func__);
    int ret = 0 ;
    struct mat_diff dat;
    mess_mvpcall mvpcall;
    mess_vector sv = NULL;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(nrm);
    mess_check_real_or_complex(A);
    mess_check_real_or_complex(B);
    mess_check_same_size(A,B);
    *nrm = -1;


    /*-----------------------------------------------------------------------------
     *  at least one dense matrix
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_DENSE(A) || MESS_IS_DENSE (B)){
        mess_matrix DA, DB;
        int convA = -1;
        MESS_MATRIX_CHECKFORMAT(A,DA,convA, MESS_DENSE);
        ret = mess_matrix_init(&DB);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_convert(B, DB, MESS_DENSE);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
        ret = mess_matrix_add(-1, DA, 1, DB);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_norm2(DB,nrm);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm2);
        if ( convA == 0 ) mess_matrix_clear(&DA);
        mess_matrix_clear(&DB);
        return(0);
    }

    dat.A = A;
    dat.B = B;
    ret = mess_mvpcall_operator(&mvpcall, A->rows, A->data_type, mat_diff_mvp, &dat);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_mvpcall_operator);
    mvpcall->data_type = A->data_type;

    ret = mess_vector_init(&sv);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sv, A->rows, A->data_type);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(sv);                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);

    ret = mess_eigen_arnoldi_template_nrm(mvpcall, MESS_MIN(50,A->rows-1), sv,  nrm);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);

    mess_mvpcall_clear(&mvpcall);
    mess_vector_clear(&sv);

    return(0);
}    /* -----  end of function mess_matrix_diffnorm  ----- */


/**
 * @brief Compute the Frobenius norm of the difference of two matrices.
 * @param[in] A   input matrix \f$A\f$
 * @param[in] B     input matrix \f$B\f$
 * @param[out] nrm  output Frobenius norm of the difference of matrix \f$A\f$ and matrix \f$B\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_diffnormf function computes the Frobenius norm of the difference of two matrices
 * \f[ nrm=\vert\vert  A-B\vert\vert_F. \f]
 *
 */
int mess_matrix_diffnormf ( mess_matrix A, mess_matrix B, double *nrm )
{
    MSG_FNAME(__func__);
    mess_int_t ret = 0 ;
    mess_matrix diff;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(nrm);
    mess_check_real_or_complex(A);
    mess_check_real_or_complex(B);
    mess_check_same_size(A,B);
    *nrm = -1;


    /*-----------------------------------------------------------------------------
     *  compute difference and frobenius norm
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&diff);
    ret = mess_matrix_copy(B,diff);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    ret = mess_matrix_add(1.0,A,-1.0,diff);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
    ret = mess_matrix_normf(diff,nrm);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);


    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&diff);

    return(0);
}    /* -----  end of function mess_matrix_diffnormf  ----- */



