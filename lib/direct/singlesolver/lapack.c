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
 * @file lib/direct/singlesolver/lapack.c
 * @brief Generate a dense @lapack based LU-solver.
 * @author @koehlerm
 *
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
 * @brief Internal structure for the @lapack based solver.
 *
 * @attention Internal use only.
 */
struct lapack_solver {
    double *val;                /**< the real lu decomposed matrix  */
    mess_double_cpx_t  *val_cpx;   /**< the complex lu decomposed matrix  */
    mess_int_t *ipiv;           /**< row permutation */
    mess_int_t n;               /**< dimension of the system */
    unsigned short cpx;         /**< flag real/complex */
};



/**
 * @brief Clear the @lapack solver.
 * @param[in,out] data  pointer to the internal data
 * @return always zero
 *
 * The @ref lapack_clear function clears a @lapack solver.
 *
 */
static int lapack_clear(void *data){
    struct lapack_solver *sol = (struct lapack_solver *) data;
    if ( sol->val != NULL)    mess_free(sol->val);
    if ( sol->val_cpx != NULL)    mess_free(sol->val_cpx);
    mess_free(sol->ipiv);
    mess_free(sol);
    return 0;
}


/**
 * @brief Solve \f$ Ax=b \f$ using @lapack.
 * @param[in] data  input pointer to the internal data structure
 * @param[in] b      input right hand side
 * @param[in,out] x     solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapack_solve function solves \f$ Ax=b \f$ using @lapack.
 *
 */
static int lapack_solve(void * data, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct lapack_solver * sol = (struct lapack_solver*) data;
    mess_int_t info = 0;
    mess_int_t one = 1;
    int ret = 0;

    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);
    ret = mess_vector_copy(b,x);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);

    if (sol->cpx == 0 && MESS_IS_REAL(x)) {
        // solve real system
        F77_GLOBAL(dgetrs,DGETRS)("N",&(sol->n), &one, sol->val, &(sol->n), sol->ipiv, x->values, &(x->dim),&info);
    } else if ( sol->cpx == 0 && MESS_IS_COMPLEX(x)){
        F77_GLOBAL(dzgetrs,DZGETRS)("N",&(sol->n), &one, sol->val, &(sol->n), sol->ipiv, x->values_cpx, &(x->dim),&info);
    } else {
        // solve complex system
        ret = mess_vector_tocomplex(x);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(zgetrs,ZGETRS)("N",&(sol->n), &one, sol->val_cpx , &(sol->n), sol->ipiv, x->values_cpx, &(x->dim), &info);
    }
    if ( info < 0) {
        MSG_ERROR("error calling DGETRS/ZGETRS. Invalid argument: " MESS_PRINTF_INT "\n", -info);
    }

    return 0;
}

/**
 * @brief Solve \f$ A^Tx=b \f$ using @lapack.
 * @param[in] data  input pointer to the internal data structure
 * @param[in] b      input right hand side
 * @param[in,out] x     solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapack_solvet function solves \f$ A^Tx=b \f$ using @lapack.
 *
 */
static int lapack_solvet(void * data, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct lapack_solver * sol = (struct lapack_solver*) data;
    mess_int_t info = 0;
    mess_int_t one = 1;
    int ret = 0;

    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    // mess_check_real_or_complex(x);
    mess_check_real_or_complex(b);

    ret = mess_vector_copy(b,x);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);

    if (sol->cpx == 0 && MESS_IS_REAL(x)) {
        // solve real system
        F77_GLOBAL(dgetrs,DGETRS)("T",&(sol->n), &one, sol->val, &(sol->n), sol->ipiv, x->values, &(x->dim),&info);
    } else if ( sol->cpx == 0 && MESS_IS_COMPLEX(x)){
        F77_GLOBAL(dzgetrs,DZGETRS)("T",&(sol->n), &one, sol->val, &(sol->n), sol->ipiv, x->values_cpx, &(x->dim),&info);
    } else {
        // solve complex system
        ret = mess_vector_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(zgetrs,ZGETRS)("T",&(sol->n), &one, sol->val_cpx, &(sol->n), sol->ipiv, x->values_cpx, &(x->dim), &info);

    }
    if ( info < 0) {
        MSG_ERROR("error calling DGETRS/ZGETRS. Invalid argument: " MESS_PRINTF_INT "\n", -info);
    }

    return 0;
}

/**
 * @brief Solve \f$ A^Hx=b \f$ using @lapack.
 * @param[in] data  input pointer to the internal data structure
 * @param[in] b      input right hand side
 * @param[in,out] x     solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapack_solveh function solves \f$ A^Hx=b \f$ using @lapack.
 *
 */
static int lapack_solveh(void * data, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    struct lapack_solver * sol = (struct lapack_solver*) data;
    mess_int_t info = 0;
    mess_int_t one = 1;
    int ret = 0;

    mess_check_nullpointer(sol);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    // mess_check_real_or_complex(x);
    mess_check_real_or_complex(b);

    ret = mess_vector_copy(b,x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);

    if (sol->cpx == 0 && MESS_IS_REAL(x)) {
        // solve real system
        F77_GLOBAL(dgetrs,DGETRS)("T",&(sol->n), &one, sol->val, &(sol->n), sol->ipiv, x->values, &(x->dim),&info);
    } else if ( sol->cpx == 0 && MESS_IS_COMPLEX(x)){
        F77_GLOBAL(dzgetrs,DZGETRS)("T",&(sol->n), &one, sol->val, &(sol->n), sol->ipiv, x->values_cpx, &(x->dim),&info);
    } else {
        // solve complex system
        ret = mess_vector_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        F77_GLOBAL(zgetrs,ZGETRS)("C",&(sol->n), &one, sol->val_cpx, &(sol->n), sol->ipiv, x->values_cpx, &(x->dim), &info);

    }
    if ( info < 0) {
        MSG_ERROR("error calling DGETRS/ZGETRS. Invalid argument: " MESS_PRINTF_INT "\n", -info);
    }

    return 0;
}

/**
 * @brief Solve \f$ AX=B \f$ using @lapack (matrix version).
 * @param[in] data  input pointer to the internal data structure
 * @param[in] b      input right hand side
 * @param[in,out] x     solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapack_solvem function solves \f$ AX=B \f$ using @lapack where \f$ B \f$ and \f$ X \f$ are matrices.
 *
 */
static int lapack_solvem(void * data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct lapack_solver * sol = (struct lapack_solver*) data;
    mess_int_t info = 0;
    int ret =0;
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);

    ret = mess_matrix_copy(b,x);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);

    if (sol->cpx == 0 && MESS_IS_REAL(x)) {
        // solve real system
        F77_GLOBAL(dgetrs,DGETRS)("N",&(sol->n), &(x->cols), sol->val, &(sol->n), sol->ipiv, x->values, &(x->ld),&info);
    } else if ( sol->cpx == 0 && MESS_IS_COMPLEX(x)){
        F77_GLOBAL(dzgetrs,DZGETRS)("N",&(sol->n), &(x->cols), sol->val, &(sol->n), sol->ipiv, x->values_cpx, &(x->ld),&info);
    } else {
        // solve complex system
        ret = mess_matrix_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        F77_GLOBAL(zgetrs,ZGETRS)("N",&(sol->n), &(x->cols), sol->val_cpx, &(sol->n), sol->ipiv, x->values_cpx, &(x->ld), &info);

    }
    if ( info < 0) {
        MSG_ERROR("error calling DGETRS/ZGETRS. Invalid argument: " MESS_PRINTF_INT "\n", -info);
    }

    return 0;

}

/**
 * @brief Solve \f$ A^TX=B \f$ using @lapack.
 * @param[in] data  input pointer to the internal data structure
 * @param[in] b      input right hand side
 * @param[in,out] x     solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapack_solvemt function solves \f$ A^TX=B \f$ using @lapack where \f$ B \f$ and \f$ X \f$ are matrices .
 *
 */
static int lapack_solvemt(void * data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct lapack_solver * sol = (struct lapack_solver*) data;
    mess_int_t info = 0;
    int ret = 0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    // mess_check_real_or_complex(x);
    mess_check_real_or_complex(b);

    ret = mess_matrix_copy(b,x);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);

    if (sol->cpx == 0 && MESS_IS_REAL(x)) {
        // solve real system
        F77_GLOBAL(dgetrs,DGETRS)("T",&(sol->n), &(x->cols), sol->val, &(sol->n), sol->ipiv, x->values, &(x->ld),&info);
    } else if ( sol->cpx == 0 && MESS_IS_COMPLEX(x)){
        F77_GLOBAL(dzgetrs,DZGETRS)("T",&(sol->n), &(x->cols), sol->val, &(sol->n), sol->ipiv, x->values_cpx, &(x->ld),&info);
    } else {
        // solve complex system
        ret = mess_matrix_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        F77_GLOBAL(zgetrs,ZGETRS)("T",&(sol->n), &(x->cols), sol->val_cpx, &(sol->n), sol->ipiv, x->values_cpx, &(x->ld), &info);

    }
    if ( info < 0) {
        MSG_ERROR("error calling DGETRS/ZGETRS. Invalid argument: " MESS_PRINTF_INT "\n", -info);
    }

    return 0;

}

/**
 * @brief Solve \f$ A^HX=B \f$ using @lapack.
 * @param[in] data  input pointer to the internal data structure
 * @param[in] b      input right hand side
 * @param[in,out] x     solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapack_solvemt function solves \f$ A^HX=B \f$ using @lapack where \f$ B \f$ and \f$ X \f$ are matrices.
 *
 */
static int lapack_solvemh(void * data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct lapack_solver * sol = (struct lapack_solver*) data;
    mess_int_t info = 0;
    int ret =0 ;

    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);
    mess_check_real_or_complex(b);

    ret = mess_matrix_copy(b,x);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);

    if (sol->cpx == 0 && MESS_IS_REAL(x)) {
        // solve real system
        F77_GLOBAL(dgetrs,DGETRS)("T",&(sol->n), &(x->cols), sol->val, &(sol->n), sol->ipiv, x->values, &(x->ld),&info);
    } else if ( sol->cpx == 0 && MESS_IS_COMPLEX(x)){
        F77_GLOBAL(dzgetrs,DZGETRS)("T",&(sol->n), &(x->cols), sol->val, &(sol->n), sol->ipiv, x->values_cpx, &(x->ld),&info);
    } else {
        // solve complex system
        ret = mess_matrix_tocomplex(x); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
        F77_GLOBAL(zgetrs,ZGETRS)("C",&(sol->n), &(x->cols), sol->val_cpx, &(sol->n), sol->ipiv, x->values_cpx, &(x->ld), &info);

    }
    if ( info < 0) {
        MSG_ERROR("error calling DGETRS/ZGETRS. Invalid argument: " MESS_PRINTF_INT "\n", -info);
    }

    return 0;

}



/**
 * @brief Compute the inverse of a matrix using @lapack.
 * @param[in] data  input pointer to the internal data
 * @param[in,out] inv output inverse
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapack_inverse function computes \f$ A^{-1} \f$ using @lapack.
 *
 */
static int lapack_inverse( void *data, mess_matrix inv) {
    MSG_FNAME(__func__);
    struct lapack_solver * sol = (struct lapack_solver*) data;
    mess_int_t info = 0;
    mess_check_nullpointer(sol);
    int ret =  0;
    mess_int_t n ;
    n = sol->n ;

    if (sol->cpx == 0) {
        double *work;
        mess_int_t lwork;
        ret = mess_matrix_alloc(inv, n,n,n*n, MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
        F77_GLOBAL(dlacpy,DLACPY)("A", &n, &n, sol->val, &sol->n, inv->values, &inv->ld);
        lwork = n * 32;
        mess_try_alloc(work, double*, sizeof(double) * lwork);
        F77_GLOBAL(dgetri,DGETRI)(&(sol->n), inv->values, &(inv->ld), sol->ipiv, work,&lwork, &info);
        mess_free(work);
    } else {
        mess_double_cpx_t *work;
        mess_int_t lwork;
        ret = mess_matrix_alloc(inv, n,n,n*n, MESS_DENSE, MESS_COMPLEX);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
        F77_GLOBAL(zlacpy,ZLACPY)("A", &n, &n, sol->val_cpx, &sol->n, inv->values_cpx, &inv->ld);
        lwork = n * 32;
        mess_try_alloc(work, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * lwork);
        F77_GLOBAL(zgetri,ZGETRI)(&(sol->n), inv->values_cpx, &(inv->ld), sol->ipiv, work,&lwork, &info);
        mess_free(work);

    }
    if ( info < 0) {
        MSG_ERROR("error calling DGETRI/ZGETRI. Invalid argument: " MESS_PRINTF_INT "\n", info);
    }

    return 0;
}

/**
 * @brief Compute the determinant of a matrix using @lapack.
 * @param[in] data  input pointer to the internal data
 * @param[out] m output mantissa of determinant
 * @param[out] e output exponent of determinant
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapack_det function computes \f$ \det( A ) \f$ using @lapack.
 *
 */
static int lapack_det( void *data, double * m, double *e) {
    MSG_FNAME(__func__);
    struct lapack_solver * sol = (struct lapack_solver*) data;
    mess_check_nullpointer(sol);
    mess_int_t i;
    mess_int_t n;
    n = sol->n ;
    double s;
    int exp;

    if (sol->cpx) {
        MSG_ERROR("Complex Matrix for real determinant computation.");
        return MESS_ERROR_ARGUMENTS;
    }

    //take produkt of diag(U)
    *m = 1.0;  *e= 0.0;
    for(i=0;i<n;i++){
        //decompose value in val=s*2^exp and update exponent and mantissa
        s = frexp(sol->val[i+i*n],&exp);
        (*e)+=exp;
        (*m)*=s;
        //decompose mantissa again and update
        s = frexp(*m,&exp);
        (*m)=s;
        (*e)+=exp;
    }

    //if det nonzero, we have to find the correct sign, by counting the number transpositions
    if(*m){
        for(i=0;i<n;++i){
            if(i!=(sol->ipiv[i]-1))
                (*m)=-(*m);
        }
    }
    return 0;
}


/**
 * @brief Compute the determinant of a matrix using @lapack.
 * @param[in] data  input pointer to the internal data
 * @param[out] mr   output mantissa for real part of determinant
 * @param[out] mi   output mantissa for imaginary part of determinant
 * @param[out] e    output exponent of determinant
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapack_detc function computes \f$ \det( A ) \f$ using @lapack.
 *
 */
static int lapack_detc( void *data, double *mr, double *mi, double *e) {
    MSG_FNAME(__func__);
    struct lapack_solver * sol = (struct lapack_solver*) data;
    mess_check_nullpointer(sol);
    mess_int_t i;
    mess_int_t n;
    n = sol->n;
    double t1,t2,t3;
    int exp;


    if (!sol->cpx) {
        MSG_ERROR("Real Matrix for complex determinant computation.");
        return MESS_ERROR_ARGUMENTS;
    }

    //take produkt of diag(U)
    *mr=1.0; *mi=0.0; *e=0.0;
    for(i=0;i<n;++i){
        //complex multiplikation with mantissa
        t1  = creal(sol->val_cpx[i+i*sol->n]);
        t2  = cimag(sol->val_cpx[i+i*sol->n]);
        t3  = (*mr)*t1 - (*mi)*t2;
        *mi = (*mr)*t2 + (*mi)*t1;
        *mr = t3;

        //normalize
        t3 = hypot(*mr,*mi); *mr/=t3; *mi/=t3;

        //put absoulute value into exponent
        t2 = frexp(t3,&exp); (*mr)*=t2; (*mi)*=t2; *e+=exp;
    }

    //if det nonzero, we have to find the correct sign, by counting the number transpositions
    if(fabs(*mr)>0 && fabs(*mi)>0){
        for(i=0;i<n;++i){
            if(i!=(sol->ipiv[i]-1)){ *mr=-*mr;*mi=-*mi;}
        }
    }
    return 0;
}




/**
 * @brief Get the L (the first factor) of a solver.
 * @param[in] data  input solver data
 * @param[in,out] L     output L factor
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapack_getL function gets the \f$ L \f$ factor of \f$ A =PLU \f$. Output is in MESS_DENSE
 * format.
 *
 */
static int lapack_getL(void *data, mess_matrix L ){
    MSG_FNAME(__func__);
    struct lapack_solver * sol = (struct lapack_solver*) data;
    mess_check_nullpointer(sol);
    mess_int_t i;
    mess_int_t n;
    n = sol->n;
    int ret=0;

    ret = mess_matrix_reset(L);                                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_reset);
    ret = mess_matrix_alloc(L,n,n,n*n, MESS_DENSE,sol->cpx? MESS_COMPLEX:MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);

    if(MESS_IS_COMPLEX(L)){
        for(i=0;i<n;++i){ L->values_cpx[i+i*L->ld]=1;}
        --n;
        F77_GLOBAL(zlacpy,ZLACPY)("L",&n,&n,sol->val_cpx+1,&(sol->n),L->values_cpx+1,&(L->ld));
    }else{
        for(i=0;i<n;++i){ L->values[i+i*L->ld]=1;}
        --n;
        F77_GLOBAL(dlacpy,DLACPY)("L",&n,&n,sol->val+1,&(sol->n),L->values+1,&(L->ld));
    }
    return ret;
}


/**
 * @brief Get the U (the second factor) of a solver.
 * @param[in] data  input solver data
 * @param[in,out] U     output U factor
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapack_getU function gets the U factor of \f$ A =PLU \f$. Output is in MESS_DENSE
 * format.
 *
 */
static int lapack_getU(void *data, mess_matrix U ){
    MSG_FNAME(__func__);
    struct lapack_solver * sol = (struct lapack_solver*) data;
    mess_check_nullpointer(sol);
    mess_int_t n;
    n = sol->n;
    int ret=0;

    ret = mess_matrix_reset(U);                                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_reset);
    ret = mess_matrix_alloc(U,n,n,n*n,MESS_DENSE,sol->cpx? MESS_COMPLEX:MESS_REAL );    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);

    if(MESS_IS_COMPLEX(U)){
        F77_GLOBAL(zlacpy,ZLACPY)("U",&(sol->n),&(sol->n),sol->val_cpx,&(sol->n),U->values_cpx,&(U->ld));
    }else{
        F77_GLOBAL(dlacpy,DLACPY)("U",&(sol->n),&(sol->n),sol->val,&(sol->n),U->values,&(U->ld));
    }
    return ret;
}


/**
 * @brief Get the row permutation.
 * @param[in] data  input solver data
 * @param[in,out] p output permutation
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref lapack_getpermp function gets the row permutation of the decomposition  \f$ PA =LU \f$.
 *
 */
static int lapack_getpermp(void *data, mess_int_t *p){
    MSG_FNAME(__func__);
    struct lapack_solver *sol = (struct lapack_solver *) data;
    mess_int_t i=0;
    mess_int_t temp;
    mess_check_nullpointer(data);
    mess_check_nullpointer(p);


    //initialize
    for (i = 0 ; i < sol->n; ++i) {
        p[i]=i;
    }

    //compute row permutation of A
    for(i=sol->n-1;i>=0;i--){
        temp = p[i];
        p[i] = p[sol->ipiv[i]-1];
        p[sol->ipiv[i]-1] = temp;
    }


    return 0;
}


/**
 * @brief Generate a dense @lapack based LU-solver.
 * @param[in] matrix     input matrix to decompose
 * @param[out] solver   generated solver
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_direct_create_lapack_lu function generates a LU solver for dense matrices based on @lapack. \n
 * It uses DGETRF/ZGETRF to compute the LU decomposition.
 *
 *
 */
int mess_direct_create_lapack_lu (mess_matrix matrix, mess_direct solver){
    MSG_FNAME(__func__);
    mess_int_t info = 0;
    int ret = 0;
    struct lapack_solver * sol ;
    mess_matrix mwork;
    int conv = 0;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_square(matrix);
    mess_check_real_or_complex(matrix);

    MESS_MATRIX_CHECKFORMAT(matrix, mwork, conv, MESS_DENSE);

    mess_try_alloc(sol, struct lapack_solver *, sizeof(struct lapack_solver));
    sol->n = matrix->rows;
    mess_try_alloc(sol->ipiv, mess_int_t *, sizeof(mess_int_t)*sol->n);

    if ( MESS_IS_COMPLEX ( matrix )) {
        sol->cpx = 1;
        mess_try_alloc(sol->val_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * matrix->rows * matrix->cols);
        sol->val = NULL;
        F77_GLOBAL(zlacpy,ZLACPY)("All", &(matrix->rows), &(matrix->cols), mwork->values_cpx, &(mwork->ld), sol->val_cpx, &(sol->n));
        F77_GLOBAL(zgetrf,ZGETRF)(&(sol->n),&(sol->n), sol->val_cpx, &(sol->n), sol->ipiv, &info);
    } else {
        sol->cpx = 0;
        mess_try_alloc(sol->val, double *, sizeof(double) * matrix->rows * matrix->cols);
        sol->val_cpx=NULL;
        F77_GLOBAL(dlacpy,DLACPY)("All", &(matrix->rows), &(matrix->cols), mwork->values, &(mwork->ld), sol->val, &(sol->n));
        F77_GLOBAL(dgetrf,DGETRF)(&(sol->n),&(sol->n), sol->val, &(sol->n), sol->ipiv, &info);
    }


    solver->data        = (void *) sol;
    solver->clear       = lapack_clear;
    solver->solve       = lapack_solve;
    solver->solvet      = lapack_solvet;
    solver->solveh      = lapack_solveh;
    solver->solvem      = lapack_solvem;
    solver->solvemt     = lapack_solvemt;
    solver->solvemh     = lapack_solvemh;
    solver->inverse     = lapack_inverse;
    solver->det         = lapack_det;
    solver->detc        = lapack_detc;
    solver->getL        = lapack_getL;
    solver->getU        = lapack_getU;
    solver->getpermp    = lapack_getpermp;
    solver->getpermq    = NULL;
    solver->rows        = matrix->rows;
    solver->cols        = matrix->cols;
    solver->data_type   = matrix->data_type;
    SET_SOLVERNAME(solver->name, __func__);

    if ( conv == 0 ) {
        mess_matrix_clear(&mwork);
    }
    if ( info != 0 ) {
        MSG_ERROR("An error occured in DGETRF/ZGETRF: " MESS_PRINTF_INT "\n", info);
        ret = (int)(info);
    }

    return ret;
}

