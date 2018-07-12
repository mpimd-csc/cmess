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
 * @file lib/direct/trisolve.c
 * @brief Solving with triangular systems.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"

/**
 * @internal
 * @brief LU solver for complex CSR stored \f$ L \f$ and \f$ U \f$.
 * @param[in] dim    input dimension
 * @param[in] lval   input values of \f$ L \f$
 * @param[in] lcolptr  input column pointer for \f$ L \f$
 * @param[in] lrowptr  input row pointer for \f$ L \f$
 * @param[in] uval   input values of \f$ U \f$
 * @param[in] ucolptr  input column pointer for \f$ U \f$
 * @param[in] urowptr  input row pointer for \f$ U \f$
 * @param[in] p  input row permuation
 * @param[in] q      input column permutaion
 * @param[in] b      input right hand side
 * @param[out] x    solution
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_lusolvem_kernelcsr_complex function solves
 * \f[ LUx=b \f]
 * for \f$ B \f$ and \f$ X \f$ matrices.
 * \f$ L \f$ and  \f$ U \f$ are given as vectors.
 * @attention Internal use only.
 */
int mess_solver_lusolvem_kernelcsr_complex (mess_int_t dim,
        mess_double_cpx_t *lval, mess_int_t *lcolptr, mess_int_t *lrowptr,
        mess_double_cpx_t *uval, mess_int_t *ucolptr, mess_int_t *urowptr,
        mess_int_t *p, mess_int_t *q,
        mess_matrix b, mess_matrix x)
{
    MSG_FNAME(__func__);
    mess_vector y;
    mess_matrix work;
    int conv = -1;
    mess_int_t i,j;
    int ret=0;
    unsigned short outdt = MESS_REAL;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(lval);
    mess_check_nullpointer(lcolptr);
    mess_check_nullpointer(lrowptr);
    mess_check_nullpointer(uval);
    mess_check_nullpointer(ucolptr);
    mess_check_nullpointer(urowptr);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
    MESS_MATRIX_RESET(x);

    outdt = MESS_COMPLEX;
    ret = mess_matrix_alloc(x, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, outdt );
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    /*-----------------------------------------------------------------------------
     *  real solver and real right hand sides
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_REAL(b)){
        ret = 0 ;
#ifdef _OPENMP
#pragma omp parallel for private (y, j, i)  default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret += mess_vector_init(&y);
            ret += mess_vector_alloc(y, work->rows, MESS_COMPLEX);
            // perm LU->p
            for (j = 0 ; j < work->rows ; j++) { y->values_cpx[j] = work->values[(p ? p [j] : j)+i*work->rows] ;}
            ret += mess_solver_lsolve_kernelcsr_complex(dim, lval, lcolptr, lrowptr,y->values_cpx);
            ret += mess_solver_usolve_kernelcsr_complex(dim, uval, ucolptr, urowptr, y->values_cpx);
            // iperm LU->q
            for (j = 0 ; j < work->rows; j++){ x->values_cpx[(q ? q [j] : j)+i*work->rows] = y->values_cpx[j] ;}
            mess_vector_clear(&y);
        }
    } else {
        ret = 0 ;
#ifdef _OPENMP
#pragma omp parallel for private (y, j, i)  default(shared) reduction(+:ret)
#endif
        for ( i = 0; i < work->cols; i++) {
            ret += mess_vector_init(&y);
            ret += mess_vector_alloc(y, work->rows, MESS_COMPLEX);
            // perm LU->p
            for (j = 0 ; j < work->rows ; j++) { y->values_cpx[j] = work->values_cpx[(p ? p [j] : j)+i*work->rows] ;}
            ret += mess_solver_lsolve_kernelcsr_complex(dim, lval, lcolptr, lrowptr,y->values_cpx);
            ret += mess_solver_usolve_kernelcsr_complex(dim, uval, ucolptr, urowptr, y->values_cpx);
            // iperm LU->q
            for (j = 0 ; j < work->rows; j++){ x->values_cpx[(q ? q [j] : j)+i*work->rows] = y->values_cpx[j] ;}
            mess_vector_clear(&y);
        }

    }

    if ( conv == 0) {
        mess_matrix_clear(&work);
    }

    if ( ret!=0) {
        MSG_ERROR("an error occured while solving\n");
        return MESS_ERROR_GENERAL;
    }

    return 0;

}       /* -----  end of function mess_solver_lusolvem_kernel_real  ----- */

/**
 * @internal
 * @brief Solve \f$ Lx=b \f$ (CSR, kernel only, real).
 * @param[in]     dim  input row dimension of \f$ L \f$
 * @param[in]     values   input values of \f$ L \f$
 * @param[in]     colptr   input column pointer for \f$ L \f$
 * @param[in]     rowptr   input row pointer for \f$ L \f$
 * @param[in,out] y on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_lsolve_kernelcsr_real function solves the system
 * \f[ Lx=b \f]
 * where \f$ L \f$ is a real lower triangular matrix stored in Compressed Row Storage.
 * @attention Internal use only.
 */
int mess_solver_lsolve_kernelcsr_real(mess_int_t dim, double *values, mess_int_t* colptr, mess_int_t *rowptr , double *y)
{
    // L\b --> y
    mess_int_t i, j;
    mess_int_t end;
    for ( i = 0; i < dim; i++){
        end = (rowptr[i+1])-1;
        for ( j = rowptr[i]; j < end; j++){
            y[i] -=  values[j] * y[colptr[j]];
        }
        y[i] /=  values[end];
    }
    return 0;
}


/**
 * @internal
 * @brief Solve \f$ Lx=b \f$ (CSR, kernel only, real \f$ L \f$, complex \f$ b \f$).
 * @param[in]     dim  input row dimension of \f$ L  \f$
 * @param[in]     values   input values of \f$ L  \f$
 * @param[in]     colptr   input column pointer for \f$ L  \f$
 * @param[in]     rowptr  input row pointer for \f$ L  \f$
 * @param[in,out] y on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_lsolve_kernelcsr_real_complex function solves the system \f[ Lx=b \f] where \f$ L \f$ is a
 * real lower triangular matrix stored in Compressed Row Storage.
 * @attention Internal use only.
 */
int mess_solver_lsolve_kernelcsr_real_complex(mess_int_t dim, double *values, mess_int_t* colptr, mess_int_t *rowptr , mess_double_cpx_t *y){
    // L\b --> y
    mess_int_t i, j;
    mess_int_t end;
    for ( i = 0; i < dim; i++){
        end = (rowptr[i+1])-1;
        for ( j = rowptr[i]; j < end; j++){
            y[i] -=  values[j] * y[colptr[j]];
        }
        y[i] /=  values[end];
    }
    return 0;
}


/**
 * @internal
 * @brief Solve \f$ Lx=b \f$ (CSR, kernel only, complex).
 * @param[in]     dim  input row dimension of  \f$ L  \f$
 * @param[in]     values   input values of  \f$ L  \f$
 * @param[in]     colptr   input column pointer for \f$   L  \f$
 * @param[in]     rowptr  input row pointer for \f$ L  \f$
 * @param[in,out] y on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_lsolve_kernelcsr_complex function solves the system \f[ Lx=b \f] where \f$ L \f$ is a
 * complex lower triangular matrix stored in Compressed Row Storage.
 * @attention Internal use only.
 */
int mess_solver_lsolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t* colptr, mess_int_t *rowptr , mess_double_cpx_t *y){
    mess_int_t i, j;
    mess_int_t end;
    for ( i = 0; i < dim; i++){
        end = (rowptr[i+1])-1;
        for ( j = rowptr[i]; j < end; j++){
            y[i] -=  values[j] * y[colptr[j]];
        }
        y[i] /=  values[end];
    }
    return 0;
}

/**
 * @brief Solve \f$ Lx=b \f$ (user interface).
 * @param[in] L  input lower triangular matrix
 * @param[in,out] y on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_lsolve function solves the system \f[ Lx=b. \f] It is the user interface
 * for \ref mess_solver_lsolve_kernelcsr_real and \ref mess_solver_lsolve_kernelcsr_complex.
 *
 */
int mess_solver_lsolve(mess_matrix L,  mess_vector y){
    // L\b --> y
    mess_int_t dim;
    MSG_FNAME(__func__);
    int ret =0 ;
    mess_int_t i_one=1;
    double d_one=1.0;
    mess_double_cpx_t cpx_one=1.0;

    mess_check_nullpointer(L);
    mess_check_nullpointer(y);
    mess_check_real_or_complex(L);
    mess_check_real_or_complex(y);
    mess_check_square(L)

        dim=L->rows;
    if ( y->dim != dim){
        MSG_ERROR("Dimension of y does not fit: " MESS_PRINTF_INT " <-> " MESS_PRINTF_INT "\n", y->dim, dim);
        return MESS_ERROR_DIMENSION;
    }

    if(MESS_IS_DENSE(L)){
        if (MESS_IS_REAL(L) && MESS_IS_REAL(y)) {
            F77_GLOBAL(dtrsm,DTRSM)("L","L","N","N", &(y->dim),&i_one, &(d_one), L->values, &(L->ld), y->values, &(y->dim));
        } else{
            ret = mess_vector_tocomplex(y);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
            ret = mess_matrix_tocomplex(L);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_tocomplex);
            F77_GLOBAL(ztrsm,ZTRSM)("L","L","N","N", &(y->dim),&i_one, &(cpx_one), L->values_cpx, &(L->ld), y->values_cpx, &(y->dim));
        }
    }else if(MESS_IS_CSR(L)){
        if (MESS_IS_REAL(L) && MESS_IS_REAL(y)){
            mess_solver_lsolve_kernelcsr_real(dim, L->values, L->colptr, L->rowptr , y->values);
        }else if (MESS_IS_REAL(L) && MESS_IS_COMPLEX(y)){
            mess_solver_lsolve_kernelcsr_real_complex(dim, L->values, L->colptr,L->rowptr,y->values_cpx);
        }else{
            ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
            mess_solver_lsolve_kernelcsr_complex(dim, L->values_cpx, L->colptr,L->rowptr,y->values_cpx);
        }
    }else if(MESS_IS_CSC(L)){
        if (MESS_IS_REAL(L) && MESS_IS_REAL(y)){
            mess_solver_utsolve_kernelcsr_real(dim, L->values, L->rowptr, L->colptr , y->values);
        }else if (MESS_IS_REAL(L) && MESS_IS_COMPLEX(y)){
            mess_solver_utsolve_kernelcsr_real_complex(dim, L->values, L->rowptr, L->colptr , y->values_cpx);
        }else{
            ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
            mess_solver_utsolve_kernelcsr_complex(dim, L->values_cpx, L->rowptr, L->colptr , y->values_cpx);
        }
    }else{
        MSG_ERROR("unsupported storagetype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}

/**
 * @brief Solve \f$ LX=B \f$ (user interface).
 * @param[in] L  input lower triangular matrix
 * @param[in,out] Y on input: right hand side matrix \f$ B \f$ \n
 *                        on output: solution matrix \f$ X \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_lsolvem function solves the system \f[ LX=B. \f] It is the user interface
 * for \ref mess_solver_lsolve_kernelcsr_real and \ref mess_solver_lsolve_kernelcsr_complex.
 *
 */
int mess_solver_lsolvem(mess_matrix L,  mess_matrix Y){
    // L\B --> Y
    MSG_FNAME(__func__);
    mess_int_t dim, i;
    int ret = 0;

    mess_check_nullpointer(L);
    mess_check_nullpointer(Y);
    mess_check_real_or_complex(L);
    mess_check_real_or_complex(Y);
    mess_check_dense(Y);
    mess_check_same_rows(L,Y);
    mess_check_square(L);

    dim=L->rows;
    if ( MESS_IS_DENSE(L) ){
        if (MESS_IS_REAL(L) && MESS_IS_REAL(Y)) {
            double alpha=1.0;
            F77_GLOBAL(dtrsm,DTRSM)("L","L","N","N", &Y->rows,&Y->cols, &alpha, L->values, &L->ld, Y->values, &Y->ld);
        } else {
            ret = mess_matrix_tocomplex(Y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_tocomplex);
            ret = mess_matrix_tocomplex(L); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_tocomplex);
            mess_double_cpx_t alpha=1.0;
            F77_GLOBAL(ztrsm,ZTRSM)("L","L","N","N", &Y->rows,&Y->cols, &alpha, L->values_cpx, &L->ld, Y->values_cpx, &Y->ld);
        }
    } else if(MESS_IS_CSR(L)){
        if (MESS_IS_REAL(L) && MESS_IS_REAL(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_lsolve_kernelcsr_real(dim, L->values, L->colptr, L->rowptr , Y->values+i*Y->ld);
            }
        }else if (MESS_IS_REAL(L) && MESS_IS_COMPLEX(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_lsolve_kernelcsr_real_complex(dim, L->values, L->colptr, L->rowptr , Y->values_cpx+i*Y->ld);
            }
        }else{
            ret = mess_matrix_tocomplex(Y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
            for(i=0; i<Y->cols; ++i){
                mess_solver_lsolve_kernelcsr_complex(dim, L->values_cpx, L->colptr, L->rowptr , Y->values_cpx+i*Y->ld);
            }
        }
    }else if(MESS_IS_CSC(L)){
        if (MESS_IS_REAL(L) && MESS_IS_REAL(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_utsolve_kernelcsr_real(dim, L->values, L->rowptr, L->colptr , Y->values+i*Y->ld);
            }
        }else if (MESS_IS_REAL(L) && MESS_IS_COMPLEX(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_utsolve_kernelcsr_real_complex(dim, L->values, L->rowptr, L->colptr , Y->values_cpx+i*Y->ld);
            }
        }else{
            ret = mess_matrix_tocomplex(Y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
            for(i=0; i<Y->cols; ++i){
                mess_solver_utsolve_kernelcsr_complex(dim, L->values_cpx, L->rowptr, L->colptr , Y->values_cpx+i*Y->ld);
            }
        }
    }else{
        MSG_ERROR("unsupported storagetype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}

/**
 * @internal
 * @brief Solve \f$ Ux=b \f$ (CSR, kernel only, real).
 * @param[in]      dim  input row dimension of \f$ U  \f$
 * @param[in]      values  input values of  \f$ U  \f$
 * @param[in]      colptr  input column pointer for \f$ U  \f$
 * @param[in]      rowptr  input row pointer for \f$ U  \f$
 * @param[in,out]  y    on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_usolve_kernelcsr_real function solves the system \f[ Ux=b \f] where \f$ U \f$ is a
 * real upper triangular matrix stored in Compressed Row Storage.
 * @attention Internal use only.
 */
int mess_solver_usolve_kernelcsr_real(mess_int_t dim, double *values, mess_int_t *colptr, mess_int_t *rowptr, double *y)
{
    mess_int_t i, j;
    for ( i = dim-1; i >= 0; i--){
        for ( j = rowptr[i]+1; j < rowptr[i+1]; j++){
            y[i]  -= values[j]*y[ colptr[j] ];
        }
        y[i] /=  values[rowptr[i] ];
    }
    return 0;
}

/**
 * @internal
 * @brief Solve \f$ Ux=b \f$ (CSR, kernel only, real \f$ U \f$, complex \f$ y \f$).
 * @param[in]      dim  input row dimension of \f$ U  \f$
 * @param[in]      values  input values of \f$ U  \f$
 * @param[in]      colptr  input column pointer for \f$ U  \f$
 * @param[in]      rowptr  input row pointer for \f$ U  \f$
 * @param[in,out]  y    on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_usolve_kernelcsr_real_complex function solves the system \f[ Ux=b \f] where \f$ U \f$ is a
 * real upper triangular matrix stored in Compressed Row Storage.
 *
 * @attention Internal use only.
 */
int mess_solver_usolve_kernelcsr_real_complex(mess_int_t dim, double *values, mess_int_t *colptr, mess_int_t *rowptr, mess_double_cpx_t *y){
    mess_int_t i, j;
    for ( i = dim-1; i >= 0; i--){
        for ( j = rowptr[i]+1; j < rowptr[i+1]; j++){
            y[i]  -= values[j]*y[ colptr[j] ];
        }
        y[i] /=  values[rowptr[i] ];
    }
    return 0;
}


/**
 * @internal
 * @brief Solve \f$ Ux=b \f$ (CSR, kernel only, complex).
 * @param[in]      dim  input row dimension of \f$ U \f$
 * @param[in]      values   input values of \f$ U \f$
 * @param[in]      colptr   input column pointer for \f$ U \f$
 * @param[in]      rowptr  input row pointer for \f$ U \f$
 * @param[in,out]  y    on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_usolve_kernelcsr_complex function solves the system \f[ Ux=b \f] where \f$ U \f$ is a
 * real lower triangular matrix stored in Compressed Row Storage.
 *
 * @attention Internal use only.
 */
int mess_solver_usolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t *colptr, mess_int_t *rowptr, mess_double_cpx_t *y){
    mess_int_t i, j;
    for ( i = dim-1; i >= 0; i--){
        for ( j = rowptr[i]+1; j < rowptr[i+1]; j++){
            y[i]  -= values[j]*y[ colptr[j] ];
        }
        y[i] /=  values[rowptr[i] ];
    }
    return 0;
}


/**
 * @brief Solve \f$ Ux=b \f$ (user interface).
 * @param[in] U      input upper triangular matrix
 * @param[in,out]   y on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_usolve function solves the system  \f[ Ux=b. \f] It is the user interface
 * for \ref mess_solver_usolve_kernelcsr_real and \ref mess_solver_usolve_kernelcsr_complex.
 *
 */
int mess_solver_usolve (mess_matrix U, mess_vector y){
    // U\y -> x
    mess_int_t dim;
    int ret = 0 ;
    mess_int_t i_one=1;double d_one=1.0; mess_double_cpx_t cpx_one=1.0;
    MSG_FNAME(__func__);

    mess_check_nullpointer(U);
    mess_check_nullpointer(y);

    mess_check_real_or_complex(U);
    mess_check_real_or_complex(y);
    mess_check_square(U)

        dim=U->rows;
    if ( y->dim != U->rows){
        MSG_WARN("Dimension of y does not fit:  " MESS_PRINTF_INT " <-> " MESS_PRINTF_INT "\n", y->dim, U->rows);
        return MESS_ERROR_DIMENSION;
    }

    if(MESS_IS_DENSE(U)){
        if (MESS_IS_REAL(U) && MESS_IS_REAL(y)) {
            F77_GLOBAL(dtrsm,DTRSM)("L","U","N","N", &(y->dim),&i_one, &(d_one), U->values, &(U->ld), y->values, &(y->dim));
        } else{
            ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
            ret = mess_matrix_tocomplex(U); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_tocomplex);
            F77_GLOBAL(ztrsm,ZTRSM)("L","U","N","N", &(y->dim),&i_one, &(cpx_one), U->values_cpx, &(U->ld), y->values_cpx, &(y->dim));
        }
    }else if(MESS_IS_CSR(U)){
        if (MESS_IS_REAL(U) && MESS_IS_REAL(y)){
            mess_solver_usolve_kernelcsr_real(dim, U->values, U->colptr, U->rowptr , y->values);
        }else if (MESS_IS_REAL(U) && MESS_IS_COMPLEX(y)){
            mess_solver_usolve_kernelcsr_real_complex(dim, U->values, U->colptr,U->rowptr,y->values_cpx);
        }else{
            ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
            mess_solver_usolve_kernelcsr_complex(dim, U->values_cpx, U->colptr,U->rowptr,y->values_cpx);
        }
    }else if(MESS_IS_CSC(U)){
        if (MESS_IS_REAL(U) && MESS_IS_REAL(y)){
            mess_solver_ltsolve_kernelcsr_real(dim, U->values, U->rowptr, U->colptr , y->values);
        }else if (MESS_IS_REAL(U) && MESS_IS_COMPLEX(y)){
            mess_solver_ltsolve_kernelcsr_real_complex(dim, U->values, U->rowptr, U->colptr , y->values_cpx);
        }else{
            ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
            mess_solver_ltsolve_kernelcsr_complex(dim, U->values_cpx, U->rowptr, U->colptr , y->values_cpx);
        }
    }else{
        MSG_ERROR("unsupported storagetype\n");
        return MESS_ERROR_DATATYPE;

    }
    return 0;
}

/**
 * @brief Solve \f$ UX=B \f$ (user interface).
 * @param[in] U  input upper triangular matrix
 * @param[in,out] Y on input: right hand side matrix \f$ B \f$ \n
 *                        on output: solution matrix \f$ X \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_usolvem function solves the system \f[ UX=B. \f] It is the user interface
 * for \ref mess_solver_usolve_kernelcsr_real and \ref mess_solver_usolve_kernelcsr_complex.
 *
 */
int mess_solver_usolvem(mess_matrix U,  mess_matrix Y){
    // U\B --> Y
    MSG_FNAME(__func__);
    mess_int_t dim, i;
    int ret = 0;

    mess_check_nullpointer(U);
    mess_check_nullpointer(Y);
    mess_check_real_or_complex(U);
    mess_check_real_or_complex(Y);
    mess_check_dense(Y);
    mess_check_same_rows(U,Y);
    mess_check_square(U);

    dim=U->rows;
    if ( MESS_IS_DENSE(U) ){
        if (MESS_IS_REAL(U) && MESS_IS_REAL(Y)) {
            double alpha=1.0;
            F77_GLOBAL(dtrsm,DTRSM)("L","U","N","N", &Y->rows,&Y->cols, &alpha, U->values, &U->ld, Y->values, &Y->ld);
        } else {
            ret = mess_matrix_tocomplex(Y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_tocomplex);
            ret = mess_matrix_tocomplex(U); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_tocomplex);
            mess_double_cpx_t alpha=1.0;
            F77_GLOBAL(ztrsm,ZTRSM)("L","U","N","N", &Y->rows,&Y->cols, &alpha, U->values_cpx, &U->ld, Y->values_cpx, &Y->ld);
        }
    } else if(MESS_IS_CSR(U)){
        if (MESS_IS_REAL(U) && MESS_IS_REAL(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_usolve_kernelcsr_real(dim, U->values, U->colptr, U->rowptr , Y->values+i*Y->ld);
            }
        }else if (MESS_IS_REAL(U) && MESS_IS_COMPLEX(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_usolve_kernelcsr_real_complex(dim, U->values, U->colptr, U->rowptr , Y->values_cpx+i*Y->ld);
            }
        }else{
            ret = mess_matrix_tocomplex(Y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
            for(i=0; i<Y->cols; ++i){
                mess_solver_usolve_kernelcsr_complex(dim, U->values_cpx, U->colptr, U->rowptr , Y->values_cpx+i*Y->ld);
            }
        }
    }else if(MESS_IS_CSC(U)){
        if (MESS_IS_REAL(U) && MESS_IS_REAL(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_ltsolve_kernelcsr_real(dim, U->values, U->rowptr, U->colptr , Y->values+i*Y->ld);
            }
        }else if (MESS_IS_REAL(U) && MESS_IS_COMPLEX(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_ltsolve_kernelcsr_real_complex(dim, U->values, U->rowptr, U->colptr , Y->values_cpx+i*Y->ld);
            }
        }else{
            ret = mess_matrix_tocomplex(Y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
            for(i=0; i<Y->cols; ++i){
                mess_solver_ltsolve_kernelcsr_complex(dim, U->values_cpx, U->rowptr, U->colptr , Y->values_cpx+i*Y->ld);
            }
        }
    }else{
        MSG_ERROR("unsupported storagetype\n");
        return MESS_ERROR_DATATYPE;
    }

    return 0;
}


/**
 * @internal
 * @brief Solve \f$ L^Tx=b \f$ (CSR, kernel only, real).
 * @param[in]      dim  input row dimension of \f$ L \f$
 * @param[in]      values   input values of \f$ L \f$
 * @param[in]      colptr   input column pointer for \f$ L \f$
 * @param[in]      rowptr   input row pointer for \f$ L \f$
 * @param[in,out]  y    on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_ltsolve_kernelcsr_real function solves the system \f[ L^T x=b \f] where \f$ L \f$ is a
 * real lower triangular matrix stored in Compressed Row Storage. The transposition is done implicitly.
 *
 * @attention Internal use only.
 */
int mess_solver_ltsolve_kernelcsr_real(mess_int_t dim, double *values, mess_int_t* colptr, mess_int_t *rowptr , double *y)
{
    mess_int_t i, j;
    for ( i = dim-1; i >= 0; i--){
        y[i] = y[i]/ values[rowptr[i+1]-1];
        for ( j = rowptr[i]; j < rowptr[i+1]-1; j++){
            y[colptr[j]]  -= y[i] * values[j];
        }
    }
    return 0;
}

/**
 * @internal
 * @brief Solve \f$ L^Tx=b \f$ (CSR, kernel only, real \f$ L \f$, complex \f$ y \f$).
 * @param[in]     dim  input row dimension of \f$ L \f$
 * @param[in]     values   input values of \f$ L \f$
 * @param[in]     colptr   input column pointer for \f$ L \f$
 * @param[in]     rowptr  input row pointer for \f$ L \f$
 * @param[in,out] y on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_ltsolve_kernelcsr_real_complex function solves the system \f[ L^T x=b \f] where \f$ L \f$ is a
 * real lower triangular matrix stored in Compressed Row Storage. The transposition is done implicitly.
 *
 * @attention Internal use only.
 */
int mess_solver_ltsolve_kernelcsr_real_complex(mess_int_t dim, double *values, mess_int_t* colptr, mess_int_t *rowptr , mess_double_cpx_t *y)
{
    mess_int_t i, j;
    for ( i = dim-1; i >= 0; i--){
        y[i] = y[i]/ values[rowptr[i+1]-1];
        for ( j = rowptr[i]; j < rowptr[i+1]-1; j++){
            y[colptr[j]]  -= y[i] * values[j];
        }
    }
    return 0;
}



/**
 * @internal
 * @brief Solve \f$ L^Tx=b \f$ (CSR, kernel only, complex).
 * @param[in]      dim  input row dimension of \f$ L \f$
 * @param[in]      values   input values of \f$ L \f$
 * @param[in]      colptr   input column pointer for \f$ L \f$
 * @param[in]      rowptr  input row pointer for \f$ L \f$
 * @param[in,out]  y    on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_ltsolve_kernelcsr_complex function solves the system \f[ L^Tx=b \f] where \f$ L \f$ is a
 * complex lower triangular matrix stored in Compressed Row Storage. The transposition is done implicitly.
 *
 * @attention Internal use only.
 */
int mess_solver_ltsolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t* colptr, mess_int_t *rowptr , mess_double_cpx_t *y)
{
    mess_int_t i, j;
    for ( i = dim-1; i >= 0; i--){
        y[i] = y[i]/ values[rowptr[i+1]-1];
        for ( j = rowptr[i]; j < rowptr[i+1]-1; j++){
            y[colptr[j]]  -= y[i] * values[j];
        }
    }
    return 0;
}

/**
 * @internal
 * @brief Solve \f$ L^Hx=b \f$ (CSR, kernel only, complex).
 * @param[in]      dim  input row dimension of \f$ L \f$
 * @param[in]      values  input values of \f$ L \f$
 * @param[in]      colptr  input column pointer for \f$ L \f$
 * @param[in]      rowptr   input row pointer for \f$ L \f$
 * @param[in,out]  y    on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_lhsolve_kernelcsr_complex function solves the system \f[ L^H x=b \f] where \f$ L \f$ is a
 * complex lower triangular matrix stored in Compressed Row Storage. The transposition is done implicitly.
 *
 * @attention Internal use only.
 */
int mess_solver_lhsolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t* colptr, mess_int_t *rowptr , mess_double_cpx_t *y)
{
    mess_int_t i, j;
    for ( i = dim-1; i >= 0; i--){
        y[i] = y[i]/ conj(values[rowptr[i+1]-1]);
        for ( j = rowptr[i]; j < rowptr[i+1]-1; j++){
            y[colptr[j]]  -= y[i] * conj(values[j]);
        }
    }
    return 0;
}


/**
 * @internal
 * @brief Solve  \f$ \bar{L}x=b \f$ (CSR, kernel only, complex).
 * @param[in]     dim  input row dimension of \f$ L \f$
 * @param[in]     values   input values of \f$ L \f$
 * @param[in]     colptr   input column pointer for \f$ L \f$
 * @param[in]     rowptr   input row pointer for \f$ L \f$
 * @param[in,out] y on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_lcsolve_kernelcsr_complex function solves the system \f[ \bar{L} x=b \f] where \f$ L \f$ is a
 * complex lower triangular matrix stored in Compressed Row Storage.
 *
 * @attention Internal use only.
 */
int mess_solver_lcsolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t* colptr, mess_int_t *rowptr , mess_double_cpx_t *y)
{
    // conj(L)\b --> y
    mess_int_t i, j;
    mess_int_t end;
    for ( i = 0; i < dim; i++){
        end = (rowptr[i+1])-1;
        for ( j = rowptr[i]; j < end; j++){
            y[i] -=  conj(values[j]) * y[colptr[j]];
        }
        y[i] /=  conj(values[end]);
    }
    return 0;
}

/**
 * @internal
 * @brief Solve \f$ \bar{U} x=b \f$ (CSR, kernel only, complex).
 * @param[in]      dim  input row dimension of \f$ U \f$
 * @param[in]      values   input values of \f$ U \f$
 * @param[in]      colptr  input column pointer for \f$ U \f$
 * @param[in]      rowptr   input row pointer for \f$ U \f$
 * @param[in,out]  y    on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_ucsolve_kernelcsr_complex function solves the system \f[ \bar{U} x=b \f] where \f$ U \f$ is a
 * complex upper triangular matrix stored in Compressed Row Storage.
 *
 * @attention Internal use only.
 */
int mess_solver_ucsolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t *colptr, mess_int_t *rowptr, mess_double_cpx_t *y)
{
    mess_int_t i, j;
    for ( i = dim-1; i >= 0; i--){
        for ( j = rowptr[i]+1; j < rowptr[i+1]; j++){
            y[i]  -= conj(values[j])*y[ colptr[j] ];
        }
        y[i] /=  conj(values[rowptr[i] ]);
    }
    return 0;
}



/**
 * @brief Solve \f$ L^Tx=b \f$ (user interface).
 * @param[in]   L  input lower triangular matrix
 * @param[in,out] y on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_ltsolve function solves the system \f[ L^Tx=b. \f] It is the user interface
 * for \ref mess_solver_ltsolve_kernelcsr_real and \ref mess_solver_ltsolve_kernelcsr_complex.
 *
 */
int mess_solver_ltsolve(mess_matrix L, mess_vector y)
{
    // L'\y --> y
    mess_int_t dim;
    int ret = 0;
    mess_int_t i_one=1;double d_one=1.0; mess_double_cpx_t cpx_one=1.0;
    MSG_FNAME(__func__);

    mess_check_nullpointer(L);
    mess_check_nullpointer(y);
    mess_check_real_or_complex(L);
    mess_check_real_or_complex(y);
    mess_check_square(L)

        dim=L->cols;
    if ( y->dim != dim){
        MSG_ERROR("Dimension of  y does not fit: " MESS_PRINTF_INT " <-> " MESS_PRINTF_INT "\n", y->dim, dim);
        return MESS_ERROR_DIMENSION;
    }


    if(MESS_IS_DENSE(L)){
        if (MESS_IS_REAL(L) && MESS_IS_REAL(y)) {
            F77_GLOBAL(dtrsm,DTRSM)("L","L","T","N", &(y->dim),&i_one, &(d_one), L->values, &(L->ld), y->values, &(y->dim));
        } else{
            ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
            ret = mess_matrix_tocomplex(L); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_tocomplex);
            F77_GLOBAL(ztrsm,ZTRSM)("L","L","T","N", &(y->dim),&i_one, &(cpx_one), L->values_cpx, &(L->ld), y->values_cpx, &(y->dim));
        }
    }else if(MESS_IS_CSR(L)){
        if (MESS_IS_REAL(L) && MESS_IS_REAL(y)){
            mess_solver_ltsolve_kernelcsr_real(dim, L->values, L->colptr, L->rowptr , y->values);
        }else if (MESS_IS_REAL(L) && MESS_IS_COMPLEX(y)){
            mess_solver_ltsolve_kernelcsr_real_complex(dim, L->values, L->colptr,L->rowptr,y->values_cpx);
        }else{
            ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
            mess_solver_ltsolve_kernelcsr_complex(dim, L->values_cpx, L->colptr,L->rowptr,y->values_cpx);
        }
    }else if(MESS_IS_CSC(L)){
        if (MESS_IS_REAL(L) && MESS_IS_REAL(y)){
            mess_solver_usolve_kernelcsr_real(dim, L->values, L->rowptr, L->colptr , y->values);
        }else if (MESS_IS_REAL(L) && MESS_IS_COMPLEX(y)){
            mess_solver_usolve_kernelcsr_real_complex(dim, L->values, L->rowptr, L->colptr , y->values_cpx);
        }else{
            ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
            mess_solver_usolve_kernelcsr_complex(dim, L->values_cpx, L->rowptr, L->colptr , y->values_cpx);
        }
    }else{
        MSG_ERROR("unsupported storagetype\n");
        return MESS_ERROR_DATATYPE;
    }

    return 0;
}



/**
 * @brief Solve \f$ L^Hx=b \f$ (user interface).
 * @param[in]   L  input lower triangular matrix
 * @param[in,out] y on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_lhsolve function solves the system \f[ L^Hx=b. \f] It is the user interface
 * for \ref mess_solver_ltsolve_kernelcsr_real and \ref mess_solver_lhsolve_kernelcsr_complex.
 *
 */
int mess_solver_lhsolve(mess_matrix L, mess_vector y)
{
    // L^H\y --> y
    mess_int_t dim;
    mess_int_t i_one=1; mess_double_cpx_t cpx_one=1.0;
    MSG_FNAME(__func__);

    mess_check_nullpointer(L);
    mess_check_nullpointer(y);

    mess_check_square(L);
    mess_check_real_or_complex(L);
    mess_check_real_or_complex(y);

    dim=L->cols;
    if ( y->dim != dim){
        MSG_ERROR("Dimension of  y does not fit: " MESS_PRINTF_INT " <-> " MESS_PRINTF_INT "\n", y->dim, dim);
        return MESS_ERROR_DIMENSION;
    }

    if (MESS_IS_REAL(L)){
        return mess_solver_ltsolve(L,y);
    } else if ( MESS_IS_COMPLEX(L) ) {
        mess_vector_tocomplex(y);

        if(MESS_IS_DENSE(L)){
            F77_GLOBAL(ztrsm,ZTRSM)("L","L","C","N", &(y->dim),&i_one, &(cpx_one), L->values_cpx, &(L->ld), y->values_cpx, &(y->dim));
        }else if(MESS_IS_CSR(L)){
            mess_solver_lhsolve_kernelcsr_complex(dim, L->values_cpx, L->colptr,L->rowptr,y->values_cpx);
        }else if(MESS_IS_CSC(L)){
            mess_solver_ucsolve_kernelcsr_complex(dim, L->values_cpx, L->rowptr, L->colptr , y->values_cpx);
        }
    } else {
        MSG_ERROR("unkown datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}

/**
 * @brief Solve \f$ L^TX=B  \f$ (user interface).
 * @param[in] L  input lower triangular matrix
 * @param[in,out] Y on input: right hand side matrix \f$ B \f$ \n
 *                        on output: solution matrix \f$ X \f$
 * @return zero on success, a non zero error code otherwise
 *
 * The @ref mess_solver_ltsolvem function solves the system \f[ L^T X=B. \f] It is the user interface
 * for \ref mess_solver_ltsolve_kernelcsr_real and \ref mess_solver_ltsolve_kernelcsr_complex.
 *
 *
 */
int mess_solver_ltsolvem(mess_matrix L,  mess_matrix Y){
    // L'\B --> Y
    MSG_FNAME(__func__);
    mess_int_t dim, i;
    int ret = 0;

    mess_check_nullpointer(L);
    mess_check_nullpointer(Y);
    mess_check_dense(Y);
    mess_check_same_rows(L,Y);
    mess_check_square(L);

    dim=L->cols;
    if ( MESS_IS_DENSE(L) ){
        if (MESS_IS_REAL(L) && MESS_IS_REAL(Y)) {
            double alpha=1.0;
            F77_GLOBAL(dtrsm,DTRSM)("L","L","T","N", &Y->rows,&Y->cols, &alpha, L->values, &L->ld, Y->values, &Y->ld);

        } else {
            ret = mess_matrix_tocomplex(Y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
            ret = mess_matrix_tocomplex(L); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
            mess_double_cpx_t alpha=1.0;
            F77_GLOBAL(ztrsm,ZTRSM)("L","L","T","N", &Y->rows,&Y->cols, &alpha, L->values_cpx, &L->ld, Y->values_cpx, &Y->ld);
        }

    } else if(MESS_IS_CSR(L)){
        if (MESS_IS_REAL(L) && MESS_IS_REAL(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_ltsolve_kernelcsr_real(dim, L->values, L->colptr, L->rowptr , Y->values+i*Y->ld);
            }
        }else if (MESS_IS_REAL(L) && MESS_IS_COMPLEX(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_ltsolve_kernelcsr_real_complex(dim, L->values, L->colptr, L->rowptr , Y->values_cpx+i*Y->ld);
            }
        }else{
            ret = mess_matrix_tocomplex(Y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
            for(i=0; i<Y->cols; ++i){
                mess_solver_ltsolve_kernelcsr_complex(dim, L->values_cpx, L->colptr, L->rowptr , Y->values_cpx+i*Y->ld);
            }
        }
    }else if(MESS_IS_CSC(L)){
        if (MESS_IS_REAL(L) && MESS_IS_REAL(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_usolve_kernelcsr_real(dim, L->values, L->rowptr, L->colptr , Y->values+i*Y->ld);
            }
        }else if (MESS_IS_REAL(L) && MESS_IS_COMPLEX(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_usolve_kernelcsr_real_complex(dim, L->values, L->rowptr, L->colptr , Y->values_cpx+i*Y->ld);
            }
        }else{
            ret = mess_matrix_tocomplex(Y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
            for(i=0; i<Y->cols; ++i){
                mess_solver_usolve_kernelcsr_complex(dim, L->values_cpx, L->rowptr, L->colptr , Y->values_cpx+i*Y->ld);
            }
        }
    }else {
        MSG_ERROR("unkown datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}



/**
 * @brief Solve \f$ L^H X=B \f$ (user interface).
 * @param[in] L  input lower triangular matrix
 * @param[in,out] Y on input: right hand side matrix \f$ B \f$ \n
 *                        on output: solution matrix \f$ X \f$
 * @return zero on success, a non zero error code otherwise
 *
 * The @ref mess_solver_lhsolvem function solves the system \f[ L^H X=B. \f] It is the user interface
 * for \ref mess_solver_ltsolve_kernelcsr_real and \ref mess_solver_lhsolve_kernelcsr_complex.
 *
 *
 */
int mess_solver_lhsolvem(mess_matrix L,  mess_matrix Y){
    // L^H\B --> Y
    MSG_FNAME(__func__);
    mess_int_t dim, i;

    mess_check_nullpointer(L);
    mess_check_nullpointer(Y);
    mess_check_dense(Y);
    mess_check_same_rows(L,Y);

    mess_check_square(L);
    mess_check_real_or_complex(L);
    mess_check_real_or_complex(Y);

    dim=L->cols;

    if (MESS_IS_REAL(L)){
        return mess_solver_ltsolvem(L,Y);
    } else if ( MESS_IS_COMPLEX(L) ) {
        mess_matrix_tocomplex(Y);

        if ( MESS_IS_DENSE(L) ){
            mess_double_cpx_t alpha=1.0;
            F77_GLOBAL(ztrsm,ZTRSM)("L","L","C","N", &Y->rows,&Y->cols, &alpha, L->values_cpx, &L->ld, Y->values_cpx, &Y->ld);
        } else if(MESS_IS_CSR(L)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_lhsolve_kernelcsr_complex(dim, L->values_cpx, L->colptr, L->rowptr , Y->values_cpx+i*Y->ld);
            }
        }else if(MESS_IS_CSC(L)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_ucsolve_kernelcsr_complex(dim, L->values_cpx, L->rowptr, L->colptr , Y->values_cpx+i*Y->ld);
            }
        }else{
            MSG_ERROR("unsupported storagetype\n");
            return MESS_ERROR_DATATYPE;
        }
    }

    return 0;
}

/**
 * @internal
 * @brief Solve \f$ U^Tx=b \f$ (CSR, kernel only, real).
 * @param[in]      dim  input row dimension of \f$ U \f$
 * @param[in]      values  input values of \f$ U \f$
 * @param[in]      colptr  input column pointer for \f$ U \f$
 * @param[in]      rowptr   input row pointer for \f$ U \f$
 * @param[in,out]  x    on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_utsolve_kernelcsr_real function solves the system \f[ U^Tx=b \f] where \f$ U \f$ is a
 * real upper triangular matrix stored in Compressed Row Storage. The transposition is done implicitly.
 *
 * @attention Internal use only.
 */
int mess_solver_utsolve_kernelcsr_real(mess_int_t dim, double *values, mess_int_t *colptr, mess_int_t *rowptr, double *x)
{
    mess_int_t i, j;
    for ( i = 0; i < dim ; i++){
        x[i] /= values[rowptr[i]];
        for ( j = rowptr[i]+1; j < rowptr[i+1]; j++){
            x[colptr[j]]  -= values[j]*x[i];
        }
    }
    return 0;
}

/**
 * @internal
 * @brief Solve \f$ U^Tx=b \f$ (CSR, kernel only, real \f$ U \f$, complex \f$ y \f$).
 * @param[in]      dim  input row dimension of \f$ U \f$
 * @param[in]      values  input values of \f$ U \f$
 * @param[in]      colptr   input column pointer for \f$ U \f$
 * @param[in]      rowptr   input row pointer for \f$ U \f$
 * @param[in,out]  x    on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_utsolve_kernelcsr_real_complex function solves the system \f[ U^Tx=b \f] where \f$ U \f$ is a
 * real upper triangular matrix stored in Compressed Row Storage. The transposition is done implicitly.
 *
 * @attention Internal use only.
 */
int mess_solver_utsolve_kernelcsr_real_complex(mess_int_t dim, double *values, mess_int_t *colptr, mess_int_t *rowptr, mess_double_cpx_t *x)
{
    mess_int_t i, j;
    for ( i = 0; i < dim ; i++){
        x[i] /= values[rowptr[i]];
        for ( j = rowptr[i]+1; j < rowptr[i+1]; j++){
            x[colptr[j]]  -= values[j]*x[i];
        }
    }
    return 0;
}


/**
 * @internal
 * @brief Solve \f$ U^Tx=b \f$ (CSR, kernel only, complex).
 * @param[in]     dim  input row dimension of \f$ U \f$
 * @param[in]     values  input values of \f$ U \f$
 * @param[in]     colptr   input column pointer for \f$ U \f$
 * @param[in]     rowptr   input row pointer for \f$ U \f$
 * @param[in,out] x on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_utsolve_kernelcsr_complex function solves the system \f[ U^T x=b \f] where \f$ U \f$ is a
 * complex upper triangular matrix stored in Compressed Row Storage. The transposition is done implicitly.
 *
 * @attention Internal use only.
 */
int mess_solver_utsolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t *colptr, mess_int_t *rowptr, mess_double_cpx_t *x)
{
    mess_int_t i, j;
    for ( i = 0; i  <dim ; i++){
        x[i] /= values[rowptr[i]];
        for ( j = rowptr[i]+1; j < rowptr[i+1]; j++){
            x[colptr[j]]  -= values[j]*x[i];
        }
    }

    return 0;
}

/**
 * @internal
 * @brief Solve \f$U^H x=b \f$ (CSR, kernel only, complex).
 * @param[in]     dim  input row dimension of \f$ U \f$
 * @param[in]     values   input values of \f$ U \f$
 * @param[in]     colptr   input column pointer for \f$ U \f$
 * @param[in]     rowptr   input row pointer for \f$ U \f$
 * @param[in,out] x on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_uhsolve_kernelcsr_complex function solves the system \f[ U^H x=b \f] where \f$ U \f$ is a
 * complex upper triangular matrix stored in Compressed Row Storage. The transposition is done implicitly.
 *
 * @attention Internal use only.
 */
int mess_solver_uhsolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t *colptr, mess_int_t *rowptr, mess_double_cpx_t *x)
{
    mess_int_t i, j;
    for ( i = 0; i  <dim ; i++){
        x[i] /= conj(values[rowptr[i]]);
        for ( j = rowptr[i]+1; j < rowptr[i+1]; j++){
            x[colptr[j]]  -= conj(values[j])*x[i];
        }
    }

    return 0;
}






/**
 * @brief Solve \f$ U^Tx=b \f$ (user interface).
 * @param[in] U  input upper triangular matrix
 * @param[in,out] y on input: right hand side \f$ b \f$ \n
 *                  on output: solution \f$ x \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_solver_utsolve function solves the system \f[ U^T x=b. \f] It is the user interface
 * for \ref mess_solver_utsolve_kernelcsr_real and \ref mess_solver_utsolve_kernelcsr_complex.
 *
 */

int mess_solver_utsolve (mess_matrix U, mess_vector y)
{
    // U'\y -> y
    int ret=0;
    mess_int_t dim;
    mess_int_t i_one=1;double d_one=1.0; mess_double_cpx_t cpx_one=1.0;
    MSG_FNAME(__func__);

    mess_check_nullpointer(U);
    mess_check_nullpointer(y);

    mess_check_square(U);
    mess_check_real_or_complex(U);
    mess_check_real_or_complex(y);

    dim=U->cols;

    if ( y->dim != U->cols){
        MSG_WARN("resize y from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", y->dim, U->rows);
        mess_vector_resize(y, U->cols);
    }

    if(MESS_IS_DENSE(U)){
        if (MESS_IS_REAL(U) && MESS_IS_REAL(y)) {
            F77_GLOBAL(dtrsm,DTRSM)("L","U","T","N", &(y->dim),&i_one, &(d_one), U->values, &(U->ld), y->values, &(y->dim));
        } else{
            ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
            ret = mess_matrix_tocomplex(U); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_tocomplex);
            F77_GLOBAL(ztrsm,ZTRSM)("L","U","T","N", &(y->dim),&i_one, &(cpx_one), U->values_cpx, &(U->ld), y->values_cpx, &(y->dim));
        }
    }else if(MESS_IS_CSR(U)){
        if (MESS_IS_REAL(U) && MESS_IS_REAL(y)){
            mess_solver_utsolve_kernelcsr_real(dim, U->values, U->colptr, U->rowptr , y->values);
        }else if (MESS_IS_REAL(U) && MESS_IS_COMPLEX(y)){
            mess_solver_utsolve_kernelcsr_real_complex(dim, U->values, U->colptr,U->rowptr,y->values_cpx);
        }else{
            ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
            mess_solver_utsolve_kernelcsr_complex(dim, U->values_cpx, U->colptr,U->rowptr,y->values_cpx);
        }
    }else if(MESS_IS_CSC(U)){
        if (MESS_IS_REAL(U) && MESS_IS_REAL(y)){
            mess_solver_lsolve_kernelcsr_real(dim, U->values, U->rowptr, U->colptr , y->values);
        }else if (MESS_IS_REAL(U) && MESS_IS_COMPLEX(y)){
            mess_solver_lsolve_kernelcsr_real_complex(dim, U->values, U->rowptr, U->colptr , y->values_cpx);
        }else{
            ret = mess_vector_tocomplex(y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
            mess_solver_lsolve_kernelcsr_complex(dim, U->values_cpx, U->rowptr, U->colptr , y->values_cpx);
        }
    } else {
        MSG_ERROR("unkown datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}


/**
 * @brief Solve \f$ U^Hx=b \f$ (user interface).
 * @param[in] U  input upper triangular matrix
 * @param[in,out] y on input: right hand side \f$ b \f$ \n
 *                        on output: solution \f$ x \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_solver_uhsolve function solves the system \f[ U^Hx=b \f] It is the user interface
 * for \ref mess_solver_uhsolve_kernelcsr_complex and \ref mess_solver_utsolve_kernelcsr_real.
 *
 */
int mess_solver_uhsolve (mess_matrix U, mess_vector y)
{
    // U'\y -> y
    mess_int_t dim;
    mess_int_t i_one=1; mess_double_cpx_t cpx_one=1.0;
    MSG_FNAME(__func__);

    mess_check_nullpointer(U);
    mess_check_nullpointer(y);
    mess_check_square(U);
    mess_check_real_or_complex(U);
    mess_check_real_or_complex(y);


    dim=U->cols;
    if ( y->dim != dim){
        MSG_ERROR("Dimension of  y does not fit: " MESS_PRINTF_INT " <-> " MESS_PRINTF_INT "\n", y->dim, dim);
        return MESS_ERROR_DIMENSION;
    }

    if (MESS_IS_REAL(U)){
        return mess_solver_utsolve(U,y);
    } else if ( MESS_IS_COMPLEX(U) ) {
        mess_vector_tocomplex(y);

        if(MESS_IS_DENSE(U)){
            F77_GLOBAL(ztrsm,ZTRSM)("L","U","C","N", &(y->dim),&i_one, &(cpx_one), U->values_cpx, &(U->ld), y->values_cpx, &(y->dim));
        }else if(MESS_IS_CSR(U)){
            mess_solver_uhsolve_kernelcsr_complex(dim, U->values_cpx, U->colptr,U->rowptr,y->values_cpx);
        }else if(MESS_IS_CSC(U)){
            mess_solver_lcsolve_kernelcsr_complex(dim, U->values_cpx, U->rowptr, U->colptr , y->values_cpx);
        }
    } else {
        MSG_ERROR("unkown datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}

/**
 * @brief Solve \f$ U^T X=B \f$ (user interface).
 * @param[in] U  input upper triangular matrix
 * @param[in,out] Y on input: right hand side matrix \f$ B \f$ \n
 *                        on output: solution matrix \f$ X \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_utsolvem function solves the system \f[ U^TX=B \f] It is the user interface
 * for \ref mess_solver_utsolve_kernelcsr_real and \ref mess_solver_utsolve_kernelcsr_complex.
 *
 */
int mess_solver_utsolvem(mess_matrix U,  mess_matrix Y){
    // U'\B --> Y
    MSG_FNAME(__func__);
    mess_int_t i;
    mess_double_cpx_t cpx_one=1.0; double d_one=1;
    int ret = 0;
    mess_int_t dim;

    mess_check_nullpointer(U);
    mess_check_nullpointer(Y);
    mess_check_dense(Y);
    mess_check_same_rows(U,Y);

    mess_check_real_or_complex(U)
        mess_check_real_or_complex(Y);
    mess_check_square(U);

    dim = U->cols;

    if ( MESS_IS_DENSE(U) ){
        if (MESS_IS_REAL(U) && MESS_IS_REAL(Y)) {
            F77_GLOBAL(dtrsm,DTRSM)("L","U","T","N", &Y->rows,&Y->cols, &d_one, U->values, &U->ld, Y->values, &Y->ld);

        } else {
            ret = mess_matrix_tocomplex(Y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
            ret = mess_matrix_tocomplex(U); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
            F77_GLOBAL(ztrsm,ZTRSM)("L","U","T","N", &Y->rows,&Y->cols, &cpx_one, U->values_cpx, &U->ld, Y->values_cpx, &Y->ld);
        }

    } else if(MESS_IS_CSR(U)){
        if (MESS_IS_REAL(U) && MESS_IS_REAL(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_utsolve_kernelcsr_real(dim, U->values, U->colptr, U->rowptr , Y->values+i*Y->ld);
            }
        }else if (MESS_IS_REAL(U) && MESS_IS_COMPLEX(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_utsolve_kernelcsr_real_complex(dim, U->values, U->colptr, U->rowptr , Y->values_cpx+i*Y->ld);
            }
        }else{
            ret = mess_matrix_tocomplex(Y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
            for(i=0; i<Y->cols; ++i){
                mess_solver_utsolve_kernelcsr_complex(dim, U->values_cpx, U->colptr, U->rowptr , Y->values_cpx+i*Y->ld);
            }
        }
    }else if(MESS_IS_CSC(U)){
        if (MESS_IS_REAL(U) && MESS_IS_REAL(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_lsolve_kernelcsr_real(dim, U->values, U->rowptr, U->colptr , Y->values+i*Y->ld);
            }
        }else if (MESS_IS_REAL(U) && MESS_IS_COMPLEX(Y)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_lsolve_kernelcsr_real_complex(dim, U->values, U->rowptr, U->colptr , Y->values_cpx+i*Y->ld);
            }
        }else{
            ret = mess_matrix_tocomplex(Y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
            for(i=0; i<Y->cols; ++i){
                mess_solver_lsolve_kernelcsr_complex(dim, U->values_cpx, U->rowptr, U->colptr , Y->values_cpx+i*Y->ld);
            }
        }
    } else {
        MSG_ERROR("unkown datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}

/**
 * @brief Solve \f$ U^H X=B \f$ (user interface).
 * @param[in] U  input upper triangular matrix
 * @param[in,out] Y on input: right hand side matrix \f$ B \f$ \n
 *                        on output: solution matrix \f$ X \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_solver_uhsolvem function solves the system \f[ U^H X=B. \f] It is the user interface
 * for \ref mess_solver_utsolve_kernelcsr_real and \ref mess_solver_uhsolve_kernelcsr_complex.
 *
 */
int mess_solver_uhsolvem(mess_matrix U,  mess_matrix Y){
    // U'\B --> Y
    MSG_FNAME(__func__);
    mess_int_t i;
    mess_int_t dim;

    mess_check_nullpointer(U);
    mess_check_nullpointer(Y);
    mess_check_dense(Y);
    mess_check_same_rows(U,Y);

    mess_check_square(U);
    mess_check_real_or_complex(U);
    mess_check_real_or_complex(Y);

    dim = U->cols;

    if (MESS_IS_REAL(U)){
        return mess_solver_utsolvem(U,Y);
    } else if ( MESS_IS_COMPLEX(U) ) {
        mess_matrix_tocomplex(Y);

        if ( MESS_IS_DENSE(U) ){
            mess_double_cpx_t alpha=1.0;
            F77_GLOBAL(ztrsm,ZTRSM)("L","U","C","N", &Y->rows,&Y->cols, &alpha, U->values_cpx, &U->ld, Y->values_cpx, &Y->ld);
        } else if(MESS_IS_CSR(U)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_uhsolve_kernelcsr_complex(dim, U->values_cpx, U->colptr, U->rowptr , Y->values_cpx+i*Y->ld);
            }
        }else if(MESS_IS_CSC(U)){
            for(i=0; i<Y->cols; ++i){
                mess_solver_lcsolve_kernelcsr_complex(dim, U->values_cpx, U->rowptr, U->colptr , Y->values_cpx+i*Y->ld);
            }
        } else {
            MSG_ERROR("unkown datatype\n");
            return MESS_ERROR_DATATYPE;
        }
    }
    return 0;
}

