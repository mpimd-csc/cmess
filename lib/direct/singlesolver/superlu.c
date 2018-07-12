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
 * @file lib/direct/singlesolver/superlu.c
 * @brief Interface to @superlu.
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#ifdef MESS_HAVE_SUPERLU

#if !defined(MESS_HAVE_SUPERLU_MT_30)
//Serial VERSION of Superlu
#include <slu_zdefs.h>

#ifdef MESS_HAVE_SUPERLU_43
    #define SINGLE SLU_SINGLE
    #define DOUBLE SLU_DOUBLE
    #define EXTRA SLU_EXTRA
#endif

//define generators for functionns double versions of routines
#define GENERATOR(name, rettype)                \
    rettype d##name(name##_ARGS);               \
    rettype s##name(name##_ARGS);               \
    rettype c##name(name##_ARGS);

//calling sequence is changed in Superlu VERSION >= 5.0
#ifdef MESS_HAVE_SUPERLU_50
#define gssvx_ARGS                                                                  \
        superlu_options_t *options, SuperMatrix *A, int *perm_c, int *perm_r,       \
        int *etree, char *equed, double *R, double *C,                              \
        SuperMatrix *L, SuperMatrix *U, void *work, int lwork,                      \
        SuperMatrix *B, SuperMatrix *X, double *recip_pivot_growth,                 \
        double *rcond, double *ferr, double *berr,                                  \
        GlobalLU_t *Glu, mem_usage_t *mem_usage, SuperLUStat_t *stat, int *info
#else
#define gssvx_ARGS                                                                  \
        superlu_options_t *options, SuperMatrix *A, int *perm_c, int *perm_r,       \
        int *etree, char *equed, double *R, double *C,                              \
        SuperMatrix *L, SuperMatrix *U, void *work, int lwork,                      \
        SuperMatrix *B, SuperMatrix *X, double *recip_pivot_growth,                 \
        double *rcond, double *ferr, double *berr,                                  \
        mem_usage_t *mem_usage, SuperLUStat_t *stat, int *info
#endif

#define Create_CompCol_Matrix_ARGS                                                  \
        SuperMatrix *A, int m, int n, int nnz,                                      \
        void *nzval, int *rowind, int *colptr,                                      \
        Stype_t stype, Dtype_t dtype, Mtype_t mtype

#define Create_Dense_Matrix_ARGS                                                    \
        SuperMatrix *X, int m, int n, void *x, int ldx,                             \
        Stype_t stype, Dtype_t dtype, Mtype_t mtype

GENERATOR(gssvx, void);
GENERATOR(Create_CompCol_Matrix, void);
GENERATOR(Create_Dense_Matrix, void);

#ifdef MESS_HAVE_SUPERLU_50
#define CALL_DGSSVX(OPTIONS,AMAT,PERMC,PERMR,ETREE,EQUED,RSCAL,CSCAL,LMAT,UMAT,WORK,LWORK,BMAT,XMAT,RECIP,RCOND,FERR,BERR,GLU,MEMUSAGE,STAT,INFO)       \
    dgssvx(OPTIONS,AMAT,PERMC,PERMR,ETREE,EQUED,RSCAL,CSCAL,LMAT,UMAT,WORK,LWORK,BMAT,XMAT,RECIP,RCOND,FERR,BERR,GLU,MEMUSAGE,STAT,INFO);

#define CALL_ZGSSVX(OPTIONS,AMAT,PERMC,PERMR,ETREE,EQUED,RSCAL,CSCAL,LMAT,UMAT,WORK,LWORK,BMAT,XMAT,RECIP,RCOND,FERR,BERR,GLU,MEMUSAGE,STAT,INFO)       \
    zgssvx(OPTIONS,AMAT,PERMC,PERMR,ETREE,EQUED,RSCAL,CSCAL,LMAT,UMAT,WORK,LWORK,BMAT,XMAT,RECIP,RCOND,FERR,BERR,GLU,MEMUSAGE,STAT,INFO);
#else
#define CALL_DGSSVX(OPTIONS,AMAT,PERMC,PERMR,ETREE,EQUED,RSCAL,CSCAL,LMAT,UMAT,WORK,LWORK,BMAT,XMAT,RECIP,RCOND,FERR,BERR,GLU,MEMUSAGE,STAT,INFO)       \
    dgssvx(OPTIONS,AMAT,PERMC,PERMR,ETREE,EQUED,RSCAL,CSCAL,LMAT,UMAT,WORK,LWORK,BMAT,XMAT,RECIP,RCOND,FERR,BERR,MEMUSAGE,STAT,INFO);

#define CALL_ZGSSVX(OPTIONS,AMAT,PERMC,PERMR,ETREE,EQUED,RSCAL,CSCAL,LMAT,UMAT,WORK,LWORK,BMAT,XMAT,RECIP,RCOND,FERR,BERR,GLU,MEMUSAGE,STAT,INFO)       \
    zgssvx(OPTIONS,AMAT,PERMC,PERMR,ETREE,EQUED,RSCAL,CSCAL,LMAT,UMAT,WORK,LWORK,BMAT,XMAT,RECIP,RCOND,FERR,BERR,MEMUSAGE,STAT,INFO);
#endif



struct superlu_solver {
    SuperMatrix A;
    SuperMatrix L;
    SuperMatrix U;
    superlu_options_t options;
    SuperLUStat_t stat;
    mess_int_t* perm_r;
    mess_int_t* perm_c;
    mess_int_t* etree;
    double * R; //scaling factors
    double * C; //scaling factors
    char equed;
    GlobalLU_t Glu;
    mem_usage_t mem_usage;
    mess_int_t info;
    mess_int_t dim;
};

static int superlu_solve_vector(trans_t trans, void*data,mess_vector b,mess_vector x){
    MSG_FNAME(__func__);
    struct superlu_solver * sol = (struct superlu_solver *)data;
    int ret = 0;
    mess_int_t nrhs=1;
    double ferr,berr;
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->dim != sol->dim) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->dim, sol->dim);
        return MESS_ERROR_DIMENSION;
    }
    if(MESS_IS_REAL(b)){
        mess_vector_toreal_nowarn(x);
    }

    if ( x->dim != sol->dim){
        ret = mess_vector_resize(x, sol->dim); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }

    /*-----------------------------------------------------------------------------
     * convert data to SUPERLU (SUPERLU Structures work on same pointers) and solve
     *-----------------------------------------------------------------------------*/

    sol->options.Trans = trans;
    if(MESS_IS_REAL(b)&& sol->A.Dtype==SLU_D){
        // copy right hand side
        mess_vector b2;
        mess_vector_init(&b2);
         mess_vector_copy(b,b2);

        // create superlu matrices
        SuperMatrix X,B;
        dCreate_Dense_Matrix(&B, b2->dim,nrhs,b2->values,b2->dim,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&X, x->dim,nrhs,x->values,x->dim,SLU_DN,SLU_D,SLU_GE);

        // solve system
        CALL_DGSSVX(&(sol->options),&(sol->A),sol->perm_c,sol->perm_r,sol->etree, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
                ,NULL,0,&B,&X,NULL,NULL,&ferr,&berr,&(sol->Glu),&(sol->mem_usage),&(sol->stat),&(sol->info));

        // clean memory
        Destroy_SuperMatrix_Store(&B);
        Destroy_SuperMatrix_Store(&X);
        MESS_CLEAR_VECTORS(&b2);
    }else if(MESS_IS_COMPLEX(b) && sol->A.Dtype==SLU_D){
        // copy right hand side
        mess_vector breal,bcpx, xreal,xcpx;
        MESS_INIT_VECTORS(&breal,&bcpx,&xreal,&xcpx);

        mess_vector_alloc(breal,b->dim,MESS_REAL);
        mess_vector_alloc(bcpx, b->dim,MESS_REAL);
        mess_vector_alloc(xreal,b->dim,MESS_REAL);
        mess_vector_alloc(xcpx, b->dim,MESS_REAL);
        mess_vector_realpart(b,breal);
        mess_vector_imagpart(b,bcpx);

        // create superlu matrices
        SuperMatrix Xreal, Xcpx, Breal, Bcpx;
        dCreate_Dense_Matrix(&Breal,breal->dim,nrhs,breal->values,breal->dim,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&Bcpx, bcpx->dim,nrhs,bcpx->values,bcpx->dim,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&Xreal, xreal->dim,nrhs,xreal->values,xreal->dim,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&Xcpx, xcpx->dim,nrhs,xcpx->values,xcpx->dim,SLU_DN,SLU_D,SLU_GE);

        // solve system
        CALL_DGSSVX(&(sol->options),&(sol->A),sol->perm_c,sol->perm_r,sol->etree, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
                ,NULL,0,&Breal,&Xreal,NULL,NULL,&ferr,&berr,&(sol->Glu),&(sol->mem_usage),&(sol->stat),&(sol->info));
        CALL_DGSSVX(&(sol->options),&(sol->A),sol->perm_c,sol->perm_r,sol->etree, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
                ,NULL,0,&Bcpx,&Xcpx,NULL,NULL,&ferr,&berr,&(sol->Glu),&(sol->mem_usage),&(sol->stat),&(sol->info));

        // get solution
        ret = mess_vector_complex_from_parts(xreal,xcpx,x);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_complex_from_parts);

        // clean memory
        Destroy_SuperMatrix_Store(&Breal);
        Destroy_SuperMatrix_Store(&Bcpx);
        Destroy_SuperMatrix_Store(&Xreal);
        Destroy_SuperMatrix_Store(&Xcpx);
        MESS_CLEAR_VECTORS(&breal,&bcpx,&xreal,&xcpx);

    }else{
        // copy rhs
        SuperMatrix X,B;
        mess_vector b2;
        mess_vector_init(&b2);
        mess_vector_copy(b,b2);
        mess_vector_tocomplex(b2);
        mess_vector_tocomplex(x);

        // create SUPERLU matrices
        zCreate_Dense_Matrix(&B, b2->dim,nrhs, (doublecomplex*) b2->values_cpx, b2->dim,SLU_DN,SLU_Z,SLU_GE);
        zCreate_Dense_Matrix(&X, x->dim, nrhs, (doublecomplex*) x->values_cpx,  x->dim, SLU_DN,SLU_Z,SLU_GE);

        // solve system
        CALL_ZGSSVX(&(sol->options),&(sol->A),sol->perm_c,sol->perm_r,sol->etree, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
               ,NULL,0,&B,&X,NULL,NULL,&ferr,&berr,&(sol->Glu),&(sol->mem_usage),&(sol->stat),&(sol->info));


        // clean memory
        Destroy_SuperMatrix_Store(&B);
        Destroy_SuperMatrix_Store(&X);
        MESS_CLEAR_VECTORS(&b2);
    }
    return 0;

 }


int superlu_solve_matrix(trans_t trans, void*data,mess_matrix b,mess_matrix x){
    MSG_FNAME(__func__);
    struct superlu_solver * sol = (struct superlu_solver *)data;
    int ret = 0;
    mess_int_t nrhs=b->cols;
    double *ferr,*berr;
    int conv;
    mess_matrix work;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->rows != sol->dim) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->rows, sol->dim);
        return MESS_ERROR_DIMENSION;
    }

    MESS_MATRIX_CHECKFORMAT(b, work, conv, MESS_DENSE);
    MESS_MATRIX_RESET(x);
    ret = mess_matrix_alloc(x, work->rows, work->cols,work->rows*work->cols, MESS_DENSE, (sol->A.Dtype==SLU_D)?MESS_COMPLEX:MESS_REAL );

    if(MESS_IS_REAL(work))
        mess_matrix_toreal(x);

    /*-----------------------------------------------------------------------------
     * convert data to SUPERLU (SUPERLU Structures work on same pointers) and solve
     *-----------------------------------------------------------------------------*/

    sol->options.Trans = trans;
    mess_try_alloc(ferr,double*,sizeof(double)*nrhs);
    mess_try_alloc(berr,double*,sizeof(double)*nrhs);

    if(MESS_IS_REAL(work)&& sol->A.Dtype==SLU_D){
        // copy right hand side
        mess_matrix b2; mess_matrix_init(&b2); mess_matrix_copy(work,b2);

        // create superlu matrices
        SuperMatrix X,B;
        dCreate_Dense_Matrix(&B, b2->rows,nrhs,b2->values,b2->ld,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&X, x->rows,nrhs,x->values,x->ld,SLU_DN,SLU_D,SLU_GE);

        // solve system
        CALL_DGSSVX(&(sol->options),&(sol->A),sol->perm_c,sol->perm_r,sol->etree, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
                ,NULL,0,&B,&X,NULL,NULL,ferr,berr,&(sol->Glu),&(sol->mem_usage),&(sol->stat),&(sol->info));

        // clean memory
        Destroy_SuperMatrix_Store(&B);
        Destroy_SuperMatrix_Store(&X);
        MESS_CLEAR_MATRICES(&b2);
        MESS_CLEAR_POINTERS(ferr,berr);
    }else if(MESS_IS_COMPLEX(work) && sol->A.Dtype==SLU_D){
        // copy right hand side
        mess_matrix breal,bcpx, xreal,xcpx;
        MESS_INIT_MATRICES(&breal,&bcpx,&xcpx,&xreal);
        mess_matrix_realpart(work,breal);
        mess_matrix_imagpart(work,bcpx);
        ret = mess_matrix_alloc(xcpx,work->rows,work->cols,work->nnz,MESS_DENSE,MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
        ret = mess_matrix_alloc(xreal,work->rows,work->cols,work->nnz,MESS_DENSE,MESS_REAL); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);

        // create superlu matrices
        SuperMatrix Xreal, Xcpx, Breal, Bcpx;
        dCreate_Dense_Matrix(&Breal,breal->rows,nrhs,breal->values,breal->rows,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&Bcpx, bcpx->rows,nrhs,bcpx->values,bcpx->rows,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&Xreal, xreal->rows,nrhs,xreal->values,xreal->rows,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&Xcpx, xcpx->rows,nrhs,xcpx->values,xcpx->rows,SLU_DN,SLU_D,SLU_GE);


        // solve system
        CALL_DGSSVX(&(sol->options),&(sol->A),sol->perm_c,sol->perm_r,sol->etree, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
                ,NULL,0,&Breal,&Xreal,NULL,NULL,ferr,berr,&(sol->Glu),&(sol->mem_usage),&(sol->stat),&(sol->info));
        CALL_DGSSVX(&(sol->options),&(sol->A),sol->perm_c,sol->perm_r,sol->etree, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
                ,NULL,0,&Bcpx,&Xcpx,NULL,NULL,ferr,berr,&(sol->Glu),&(sol->mem_usage),&(sol->stat),&(sol->info));

        // get solution
        ret = mess_matrix_complex_from_parts(xreal,xcpx,x);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_complex_from_parts);

        // clean memory
        Destroy_SuperMatrix_Store(&Breal);
        Destroy_SuperMatrix_Store(&Bcpx);
        Destroy_SuperMatrix_Store(&Xreal);
        Destroy_SuperMatrix_Store(&Xcpx);
        MESS_CLEAR_MATRICES(&breal,&bcpx,&xreal,&xcpx);
        MESS_CLEAR_POINTERS(ferr,berr);
    }else{
        // copy rhs
        SuperMatrix X,B;
        mess_matrix b2;
        MESS_INIT_MATRICES(&b2);mess_matrix_copy(work,b2);mess_matrix_tocomplex(b2);
        mess_matrix_tocomplex(x);

        // create SUPERLU matrices
        zCreate_Dense_Matrix(&B, b2->rows,nrhs, (doublecomplex*) b2->values_cpx, b2->ld,SLU_DN,SLU_Z,SLU_GE);
        zCreate_Dense_Matrix(&X, x->rows, nrhs, ( doublecomplex*) x->values_cpx,  x->ld, SLU_DN,SLU_Z,SLU_GE);

        // solve system
        CALL_ZGSSVX(&(sol->options),&(sol->A),sol->perm_c,sol->perm_r,sol->etree, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
               ,NULL,0,&B,&X,NULL,NULL,ferr,berr,&(sol->Glu),&(sol->mem_usage),&(sol->stat),&(sol->info));

        // clean memory
        Destroy_SuperMatrix_Store(&B);
        Destroy_SuperMatrix_Store(&X);
        MESS_CLEAR_MATRICES(&b2);
        MESS_CLEAR_POINTERS(ferr,berr);
    }

    if ( conv == 0) {
        mess_matrix_clear(&work);
    }


    return 0;

}


/**
 * @brief Solve \f$ Ax=b \f$.
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref superlu_solve function solves \f$ Ax=b \f$.
 *
 */
static int superlu_solve(void*data, mess_vector b, mess_vector x) {
    return superlu_solve_vector(NOTRANS,data,b,x);
}

/**
 * @brief Solve \f$ A^Tx=b \f$.
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref superlu_solvet function solves \f$ A^Tx=b \f$.
 *
 */
static int superlu_solvet(void*data, mess_vector b, mess_vector x) {
    return superlu_solve_vector(TRANS,data,b,x);
}

/**
 * @brief Solve \f$ A^Hx=b \f$.
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref superlu_solveh function solves \f$ A^Hx=b \f$.
 *
 */
static int superlu_solveh(void*data, mess_vector b, mess_vector x) {
#if !defined(SUPERLU_REFINEMENT_BUG_FIXED)
    MSG_FNAME(__func__);
    int ret = 0;
    ret = mess_vector_conj(b);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_conj);
    ret = superlu_solve_vector(TRANS,data,b,x);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),superlu_solve_vector);
    ret = mess_vector_conj(b);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_conj);
    ret = mess_vector_conj(x);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_conj);
    return ret;
#else
    return superlu_solve_vector(CONJ,data,b,x);
#endif

}

/**
 * @brief Solve \f$ AX=B \f$ (matrix version).
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref superlu_solvem function solves \f$ AX=B \f$.
 *
 */
static int superlu_solvem(void*data, mess_matrix b, mess_matrix x) {
    return superlu_solve_matrix(NOTRANS,data,b,x);
}

/**
 * @brief Solve \f$ A^TX=B \f$ (matrix version).
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref superlu_solvemt function solves \f$ AX=B \f$.
 *
 */
static int superlu_solvemt(void*data, mess_matrix b, mess_matrix x) {
    return superlu_solve_matrix(TRANS,data,b,x);
}

/**
 * @brief Solve \f$ A^HX=B \f$ (matrix version).
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref superlu_solvemh function solves \f$ AX=B \f$.
 *
 */
static int superlu_solvemh(void*data, mess_matrix b, mess_matrix x) {
#if !defined(SUPERLU_REFINEMENT_BUG_FIXED)
    MSG_FNAME(__func__);
    int ret = 0;
    ret = mess_matrix_conj(b);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_conj);
    ret = superlu_solve_matrix(TRANS,data,b,x);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),superlu_solve_matrix);
    ret = mess_matrix_conj(b);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_conj);
    ret = mess_matrix_conj(x);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_conj);
    return ret;
#else
    return superlu_solve_matrix(CONJ,data,b,x);
#endif

}


/**
 * @brief Clear a SUPERLU object.
 * @param [in,out]solver pointer to the data object
 * @return always zero
 *
 * The @ref superlu_clear function clears an superlu object.
 *
 */
static int superlu_clear(void *solver){
    //MSG_FNAME(__func__);
    struct superlu_solver * sol = (struct superlu_solver*) solver;
    if ( sol != NULL) {
        Destroy_CompCol_Matrix(&(sol->A));
        Destroy_SuperNode_Matrix(&(sol->L));
        Destroy_CompCol_Matrix(&(sol->U));
        StatFree(&(sol->stat));
        if(sol->perm_r) mess_free(sol->perm_r);
        if(sol->perm_c) mess_free(sol->perm_c);
        if(sol->etree)  mess_free(sol->etree);
        if(sol->R) mess_free(sol->R);
        if(sol->C) mess_free(sol->C);
        mess_free( sol );
    }
    return 0;
}

/**
 * @brief Generate a direct linear solver for standard linear systems \f$ Ax=b \f$ with @superlu.
 * @param[in] matrix matrix to decompose
 * @param[out] solver output solver
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref mess_direct_create_superlu function factorizes a matrix with @superlu and provides
 * a solver for linear systems.
 *
 */
int mess_direct_create_superlu(mess_matrix matrix, mess_direct solver){
    MSG_FNAME(__func__);
    struct superlu_solver * sol;
    int ret;
    mess_matrix temp;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_square(matrix);
    mess_check_real_or_complex(matrix);

    /*-----------------------------------------------------------------------------
     *  init solver structure
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(sol, struct superlu_solver*, sizeof(struct superlu_solver));
    sol->dim = matrix->rows;
    sol->info = 0;
    mess_try_alloc(sol->perm_r, mess_int_t*, sizeof(mess_int_t)*sol->dim);
    mess_try_alloc(sol->perm_c, mess_int_t*, sizeof(mess_int_t)*sol->dim);
    mess_try_alloc(sol->etree, mess_int_t*, sizeof(mess_int_t)*sol->dim);
    mess_try_alloc(sol->R, double*, sizeof(double)*sol->dim);
    mess_try_alloc(sol->C, double*, sizeof(double)*sol->dim);

    /*-----------------------------------------------------------------------------
     * set SUPERLU options
     *-----------------------------------------------------------------------------*/
    set_default_options(&(sol->options));
    sol->options.Equil = YES; //NO
    sol->options.ColPerm = COLAMD;
    sol->options.DiagPivotThresh = 1.0;
    sol->options.Trans = NOTRANS;
    sol->options.IterRefine = EXTRA;
    sol->options.SymmetricMode = NO;
    sol->options.PivotGrowth = NO; //can provide NULL pointer for  recip_pivot_growth
    sol->options.ConditionNumber = NO; //can provide NULL pointer for rcond
    sol->options.PrintStat = NO;

    /*-----------------------------------------------------------------------------
     *  init Status SUPERLU
     *-----------------------------------------------------------------------------*/
    StatInit(&(sol->stat));

    /*-----------------------------------------------------------------------------
     *  create SUPERLU matrix
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&temp);
    ret = mess_matrix_convert(matrix,temp,MESS_CSC);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert);
    if(MESS_IS_REAL(temp)){
        dCreate_CompCol_Matrix(&(sol->A),sol->dim, sol->dim,temp->nnz,temp->values,temp->rowptr,temp->colptr,SLU_NC,SLU_D,SLU_GE);
    }else{
        zCreate_CompCol_Matrix(&(sol->A),sol->dim, sol->dim,temp->nnz,(doublecomplex*)temp->values_cpx,temp->rowptr,temp->colptr,SLU_NC,SLU_Z,SLU_GE);
    }
    //pointers are hold by SuperMatrix A and cleared by SUPERLU Destroy_CompCol_Matrix function//
    temp->values_cpx=NULL;
    temp->values=NULL;
    temp->rowptr=NULL;
    temp->colptr=NULL;

    /*-----------------------------------------------------------------------------
     *  factorization
     *-----------------------------------------------------------------------------*/
    //factorize
    SuperMatrix X,B;
    if(MESS_IS_REAL(temp)){
        dCreate_Dense_Matrix(&B, temp->rows, 0, NULL, temp->rows,SLU_DN,SLU_D,SLU_GE); //only factorization
        dCreate_Dense_Matrix(&X, temp->rows, 0, NULL, temp->rows,SLU_DN,SLU_D,SLU_GE);
        CALL_DGSSVX(&(sol->options),&(sol->A),sol->perm_c,sol->perm_r,sol->etree, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
            ,NULL,0,&B,&X,NULL,NULL,NULL,NULL,&(sol->Glu),&(sol->mem_usage),&(sol->stat),&(sol->info));
        solver->data_type = MESS_REAL;
    }else{
        zCreate_Dense_Matrix(&B, temp->rows, 0, NULL, temp->rows,SLU_DN,SLU_Z,SLU_GE); //only factorization
        zCreate_Dense_Matrix(&X, temp->rows, 0, NULL, temp->rows,SLU_DN,SLU_Z,SLU_GE);
        CALL_ZGSSVX(&(sol->options),&(sol->A),sol->perm_c,sol->perm_r,sol->etree, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
            ,NULL,0,&B,&X,NULL,NULL,NULL,NULL,&(sol->Glu),&(sol->mem_usage),&(sol->stat),&(sol->info));
        solver->data_type = MESS_COMPLEX;
    }

    sol->options.Fact = FACTORED; //Matrix is now factorized

    /*-----------------------------------------------------------------------------
     *  set solver functions
     *-----------------------------------------------------------------------------*/
    SET_SOLVERNAME(solver->name, __func__);
    solver->data = (void *)sol;
    solver->solve = superlu_solve;
    solver->solvet = superlu_solvet;
    solver->solveh = superlu_solveh;
    solver->solvem = superlu_solvem;
    solver->solvemt = superlu_solvemt;
    solver->solvemh = superlu_solvemh;
    solver->clear = superlu_clear;
    solver->getL = NULL;
    solver->getU = NULL;
    solver->getpermp = NULL;
    solver->getpermq = NULL;
    solver->inverse = NULL;

    /*-----------------------------------------------------------------------------
     * free memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&(temp));
    Destroy_Dense_Matrix(&X);
    Destroy_Dense_Matrix(&B);

    return 0;
}

#else
//Superlu MT Version >= 3.0
#include "slu_mt_zdefs.h"

//define generators for functionns double versions of routines
#define GENERATOR(name, rettype)                \
    rettype d##name(name##_ARGS);               \
    rettype s##name(name##_ARGS);               \
    rettype c##name(name##_ARGS);

#define PGENERATOR(name, rettype)               \
    rettype pd##name(name##_ARGS);              \
    rettype ps##name(name##_ARGS);              \
    rettype pc##name(name##_ARGS);

//create matrices functions
#define Create_CompCol_Matrix_ARGS                                                  \
        SuperMatrix *A, int m, int n, int nnz,                                      \
        void *nzval, int *rowind, int *colptr,                                      \
        Stype_t stype, Dtype_t dtype, Mtype_t mtype

#define Create_Dense_Matrix_ARGS                                                    \
        SuperMatrix *X, int m, int n, void *x, int ldx,                             \
        Stype_t stype, Dtype_t dtype, Mtype_t mtype

#define gssvx_ARGS                                                                                      \
        mess_int_t nprocs, superlumt_options_t* superlumt_options, SuperMatrix *A,                      \
        mess_int_t *perm_c, mess_int_t *perm_r, equed_t *equed, double *R, double *C,                   \
        SuperMatrix *L, SuperMatrix *U, SuperMatrix *B, SuperMatrix *X, double *recip_pivot_growth,     \
        double *rcond, double *ferr, double *berr, superlu_memusage_t *superlu_memusage, mess_int_t *info

GENERATOR(Create_CompCol_Matrix, void);
GENERATOR(Create_Dense_Matrix, void);
PGENERATOR(gssvx,void);

#define ORDERING_NATURAL  0
#define ORDERING_MINIMUM_DEGREE_ATPA 1
#define ORDERING_MINIMUM_DEGREE_ATMA 2
#define ORDERING_MINIMUM_DEGREE_AMD 3


struct superlu_solver {
    SuperMatrix A;
    SuperMatrix L;
    SuperMatrix U;
    superlumt_options_t options;
    double * R; //scaling factors
    double * C; //scaling factors
    mess_int_t * perm_c;
    mess_int_t * perm_r;
    mess_int_t * etree;
    equed_t equed;
    double recip_pivot_growth;
    double rcond;
    superlu_memusage_t mem_usage;
    mess_int_t info;
    mess_int_t dim;
};


static int superlu_solve_vector(trans_t trans, void*data,mess_vector b,mess_vector x){
    MSG_FNAME(__func__);
    struct superlu_solver * sol = (struct superlu_solver *)data;
    int ret = 0;
    mess_int_t nrhs=1;
    double ferr,berr;
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->dim != sol->dim) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->dim, sol->dim);
        return MESS_ERROR_DIMENSION;
    }
    if(MESS_IS_REAL(b))
        mess_vector_toreal_nowarn(x);

    if ( x->dim != sol->dim){
        ret = mess_vector_resize(x, sol->dim); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }

    /*-----------------------------------------------------------------------------
     * convert data to SUPERLU (SUPERLU Structures work on same pointers) and solve
     *-----------------------------------------------------------------------------*/

    sol->options.trans = trans;
    if(MESS_IS_REAL(b)&& sol->A.Dtype==SLU_D){
        // copy right hand side
        mess_vector b2; mess_vector_init(&b2,b->dim,b->data_type); mess_vector_copy(b,b2);

        // create superlu matrices
        SuperMatrix X,B;
        dCreate_Dense_Matrix(&B, b2->dim,nrhs,b2->values,b2->dim,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&X, x->dim,nrhs,x->values,x->dim,SLU_DN,SLU_D,SLU_GE);

        // solve system
        pdgssvx(sol->options.nprocs,&(sol->options),&(sol->A),sol->perm_c,sol->perm_r, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
            ,&B,&X,&(sol->recip_pivot_growth),&(sol->rcond),&ferr,&berr,&(sol->mem_usage),&(sol->info));

        // clean memory
        Destroy_SuperMatrix_Store(&B);
        Destroy_SuperMatrix_Store(&X);
        MESS_CLEAR_VECTORS(&b2);
    }else if(MESS_IS_COMPLEX(b) && sol->A.Dtype==SLU_D){
        // copy right hand side
        mess_vector breal,bcpx, xreal,xcpx;
        mess_vector_init(&breal,b->dim,MESS_REAL);
        mess_vector_init(&bcpx, b->dim,MESS_REAL);
        mess_vector_init(&xreal,b->dim,MESS_REAL);
        mess_vector_init(&xcpx, b->dim,MESS_REAL);
        mess_vector_realpart(b,breal);
        mess_vector_imagpart(b,bcpx);

        // create superlu matrices
        SuperMatrix Xreal, Xcpx, Breal, Bcpx;
        dCreate_Dense_Matrix(&Breal,breal->dim,nrhs,breal->values,breal->dim,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&Bcpx, bcpx->dim,nrhs,bcpx->values,bcpx->dim,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&Xreal, xreal->dim,nrhs,xreal->values,xreal->dim,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&Xcpx, xcpx->dim,nrhs,xcpx->values,xcpx->dim,SLU_DN,SLU_D,SLU_GE);

        // solve system
        pdgssvx(sol->options.nprocs,&(sol->options),&(sol->A),sol->perm_c,sol->perm_r, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
            ,&Breal,&Xreal,&(sol->recip_pivot_growth),&(sol->rcond),&ferr,&berr,&(sol->mem_usage),&(sol->info));
        pdgssvx(sol->options.nprocs,&(sol->options),&(sol->A),sol->perm_c,sol->perm_r, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
            ,&Bcpx,&Xcpx,&(sol->recip_pivot_growth),&(sol->rcond),&ferr,&berr,&(sol->mem_usage),&(sol->info));

        // get solution
        ret = mess_vector_complex_from_parts(xreal,xcpx,x);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_complex_from_parts);

        // clean memory
        Destroy_SuperMatrix_Store(&Breal);
        Destroy_SuperMatrix_Store(&Bcpx);
        Destroy_SuperMatrix_Store(&Xreal);
        Destroy_SuperMatrix_Store(&Xcpx);
        MESS_CLEAR_VECTORS(&breal,&bcpx,&xreal,&xcpx);

    }else{
        // copy rhs
        SuperMatrix X,B;
        mess_vector b2;
        mess_vector_init(&b2,b->dim,MESS_COMPLEX);mess_vector_copy(b,b2);mess_vector_tocomplex(b2);
        mess_vector_tocomplex(x);

        // create SUPERLU matrices
        zCreate_Dense_Matrix(&B, b2->dim,nrhs, (doublecomplex*) b2->values_cpx, b2->dim,SLU_DN,SLU_Z,SLU_GE);
        zCreate_Dense_Matrix(&X, x->dim, nrhs, (doublecomplex*) x->values_cpx,  x->dim, SLU_DN,SLU_Z,SLU_GE);

        // solve system
        pzgssvx(sol->options.nprocs,&(sol->options),&(sol->A),sol->perm_c,sol->perm_r, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
            ,&B,&X,&(sol->recip_pivot_growth),&(sol->rcond),&ferr,&berr,&(sol->mem_usage),&(sol->info));

        // clean memory
        Destroy_SuperMatrix_Store(&B);
        Destroy_SuperMatrix_Store(&X);
        MESS_CLEAR_VECTORS(&b2);
    }
    return 0;

 }

/**
 * @brief Solve \f$ Ax=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The superlu_solve function solves \f$ Ax=b \f$.
 *
 */
static int superlu_solve(void*data, mess_vector b, mess_vector x) {
    return superlu_solve_vector(NOTRANS,data,b,x);
}

/**
 * @brief Solve \f$ A^Tx=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The superlu_solvet function solves \f$ A^Tx=b \f$.
 *
 */
static int superlu_solvet(void*data, mess_vector b, mess_vector x) {
    return superlu_solve_vector(TRANS,data,b,x);
}

/**
 * @brief Solve \f$ A^Hx=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The superlu_solveh function solves \f$ A^Hx=b \f$.
 *
 */
static int superlu_solveh(void*data, mess_vector b, mess_vector x) {
    return superlu_solve_vector(CONJ,data,b,x);
}

int superlu_solve_matrix(trans_t trans, void*data,mess_matrix b,mess_matrix x){
    MSG_FNAME(__func__);
    struct superlu_solver * sol = (struct superlu_solver *)data;
    int ret = 0;
    mess_int_t nrhs=b->cols;
    double *ferr,*berr;
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->rows != sol->dim) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->rows, sol->dim);
        return MESS_ERROR_DIMENSION;
    }
    if(MESS_IS_REAL(b))
        mess_matrix_toreal(x);

    if ( x->rows != sol->dim){
        ret = mess_matrix_resize(x, sol->dim, b->cols);FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_resize);
    }

    /*-----------------------------------------------------------------------------
     * convert data to SUPERLU (SUPERLU Structures work on same pointers) and solve
     *-----------------------------------------------------------------------------*/

    sol->options.trans = trans;
    mess_try_alloc(ferr,double*,sizeof(double)*nrhs);
    mess_try_alloc(berr,double*,sizeof(double)*nrhs);

    if(MESS_IS_REAL(b)&& sol->A.Dtype==SLU_D){
        // copy right hand side
        mess_matrix b2; mess_matrix_init(&b2); mess_matrix_copy(b,b2);

        // create superlu matrices
        SuperMatrix X,B;
        dCreate_Dense_Matrix(&B, b2->rows,nrhs,b2->values,b2->ld,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&X, x->rows,nrhs,x->values,x->ld,SLU_DN,SLU_D,SLU_GE);

        // solve system
        pdgssvx(sol->options.nprocs,&(sol->options),&(sol->A),sol->perm_c,sol->perm_r, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
            ,&B,&X,&(sol->recip_pivot_growth),&(sol->rcond),ferr,berr,&(sol->mem_usage),&(sol->info));

        // clean memory
        Destroy_SuperMatrix_Store(&B);
        Destroy_SuperMatrix_Store(&X);
        MESS_CLEAR_MATRICES(&b2);
        MESS_CLEAR_POINTERS(ferr,berr);
    }else if(MESS_IS_COMPLEX(b) && sol->A.Dtype==SLU_D){
        // copy right hand side
        mess_matrix breal,bcpx, xreal,xcpx;
        MESS_INIT_MATRICES(&breal,&bcpx,&xcpx,&xreal);
        mess_matrix_realpart(b,breal);
        mess_matrix_imagpart(b,bcpx);
        ret = mess_matrix_alloc(xcpx,b->rows,b->cols,b->nnz,MESS_DENSE,MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
        ret = mess_matrix_alloc(xreal,b->rows,b->cols,b->nnz,MESS_DENSE,MESS_REAL); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);

        // create superlu matrices
        SuperMatrix Xreal, Xcpx, Breal, Bcpx;
        dCreate_Dense_Matrix(&Breal,breal->rows,nrhs,breal->values,breal->rows,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&Bcpx, bcpx->rows,nrhs,bcpx->values,bcpx->rows,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&Xreal, xreal->rows,nrhs,xreal->values,xreal->rows,SLU_DN,SLU_D,SLU_GE);
        dCreate_Dense_Matrix(&Xcpx, xcpx->rows,nrhs,xcpx->values,xcpx->rows,SLU_DN,SLU_D,SLU_GE);

        // solve system
        pdgssvx(sol->options.nprocs,&(sol->options),&(sol->A),sol->perm_c,sol->perm_r, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
            ,&Breal,&Xreal,&(sol->recip_pivot_growth),&(sol->rcond),ferr,berr,&(sol->mem_usage),&(sol->info));
        pdgssvx(sol->options.nprocs,&(sol->options),&(sol->A),sol->perm_c,sol->perm_r, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
            ,&Bcpx,&Xcpx,&(sol->recip_pivot_growth),&(sol->rcond),ferr,berr,&(sol->mem_usage),&(sol->info));

        // get solution
        ret = mess_matrix_complex_from_parts(xreal,xcpx,x);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_complex_from_parts);

        // clean memory
        Destroy_SuperMatrix_Store(&Breal);
        Destroy_SuperMatrix_Store(&Bcpx);
        Destroy_SuperMatrix_Store(&Xreal);
        Destroy_SuperMatrix_Store(&Xcpx);
        MESS_CLEAR_MATRICES(&breal,&bcpx,&xreal,&xcpx);
        MESS_CLEAR_POINTERS(ferr,berr);
    }else{
        // copy rhs
        SuperMatrix X,B;
        mess_matrix b2;
        MESS_INIT_MATRICES(&b2);mess_matrix_copy(b,b2);mess_matrix_tocomplex(b2);
        mess_matrix_tocomplex(x);

        // create SUPERLU matrices
        zCreate_Dense_Matrix(&B, b2->rows,nrhs, (doublecomplex*) b2->values_cpx, b2->ld,SLU_DN,SLU_Z,SLU_GE);
        zCreate_Dense_Matrix(&X, x->rows, nrhs, ( doublecomplex*) x->values_cpx,  x->ld, SLU_DN,SLU_Z,SLU_GE);

        // solve system
        pzgssvx(sol->options.nprocs,&(sol->options),&(sol->A),sol->perm_c,sol->perm_r, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
            ,&B,&X,&(sol->recip_pivot_growth),&(sol->rcond),ferr,berr,&(sol->mem_usage),&(sol->info));

        // clean memory
        Destroy_SuperMatrix_Store(&B);
        Destroy_SuperMatrix_Store(&X);
        MESS_CLEAR_MATRICES(&b2);
        MESS_CLEAR_POINTERS(ferr,berr);
    }
    return 0;

}

/**
 * @brief Solve \f$ Ax=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The superlu_solvem function solves \f$ Ax=b \f$.
 *
 */
static int superlu_solvem(void*data, mess_matrix b, mess_matrix x) {
    return superlu_solve_matrix(NOTRANS,data,b,x);
}

/**
 * @brief Solve \f$ A^Tx=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The superlu_solvemt function solves \f$ Ax=b \f$.
 *
 */
static int superlu_solvemt(void*data, mess_matrix b, mess_matrix x) {
    return superlu_solve_matrix(TRANS,data,b,x);
}

/**
 * @brief Solve \f$ A^Hx=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The superlu_solvemh function solves \f$ Ax=b \f$.
 *
 */
static int superlu_solvemh(void*data, mess_matrix b, mess_matrix x) {
    return superlu_solve_matrix(CONJ,data,b,x);
}

/**
 * @brief Clear an SUPERLU object
 * @param [in,out]solver pointer to the data object
 * @return always zero
 *
 * The umfpack_clear function clears an UMFPACK object.
 *
 */
static int superlu_clear(void *solver){
    //MSG_FNAME(__func__);
    struct superlu_solver * sol = (struct superlu_solver*) solver;
    if ( sol != NULL) {
        Destroy_CompCol_Matrix(&(sol->A));
        Destroy_SuperNode_Matrix(&(sol->L));
        Destroy_CompCol_Matrix(&(sol->U));
        if(sol->R) mess_free(sol->R);
        if(sol->C) mess_free(sol->C);
        if(sol->perm_r) mess_free(sol->perm_r);
        if(sol->perm_c) mess_free(sol->perm_c);
        if(sol->options.etree) mess_free(sol->options.etree);
        if(sol->options.colcnt_h) mess_free(sol->options.colcnt_h);
        if(sol->options.part_super_h) mess_free(sol->options.part_super_h);
        mess_free( sol );
    }
    solver = NULL;
    return 0;
}



int mess_direct_create_superlu(mess_matrix matrix, mess_direct solver){
    MSG_FNAME(__func__);
    struct superlu_solver * sol;
    int ret;
    mess_matrix temp;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(solver);
    mess_check_square(matrix);
    mess_check_real_or_complex(matrix);

    /*-----------------------------------------------------------------------------
     *  init solver structure
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(sol, struct superlu_solver*, sizeof(struct superlu_solver));
    sol->dim = matrix->rows;
    sol->info = 0;
    mess_try_alloc(sol->R, double*, sizeof(double)*sol->dim);
    mess_try_alloc(sol->C, double*, sizeof(double)*sol->dim);
    mess_try_alloc(sol->perm_r, mess_int_t*, sizeof(mess_int_t)*sol->dim);
    mess_try_alloc(sol->perm_c, mess_int_t*, sizeof(mess_int_t)*sol->dim);
    mess_try_alloc(sol->options.etree, mess_int_t*, sizeof(mess_int_t)*sol->dim);
    mess_try_alloc(sol->options.colcnt_h, mess_int_t*, sizeof(mess_int_t)*sol->dim);
    mess_try_alloc(sol->options.part_super_h, mess_int_t*, sizeof(mess_int_t)*sol->dim);

    /*-----------------------------------------------------------------------------
     * set SUPERLU options
     *-----------------------------------------------------------------------------*/
    sol->options.nprocs = omp_get_num_procs();
    sol->options.fact = DOFACT;
    sol->options.trans = NOTRANS;
    sol->options.refact = NO;
    sol->options.panel_size = sp_ienv(1);
    sol->options.relax = sp_ienv(2);
    sol->options.diag_pivot_thresh = 1.0;
    sol->options.usepr = NO; //determine perm_r by partial pivoting
    sol->options.drop_tol = 0.0; //dummy
    sol->options.SymmetricMode = NO;
    sol->options.PrintStat = NO;
    sol->options.perm_c = sol->perm_c;
    sol->options.perm_r = sol->perm_r;
    sol->options.work = NULL; //let superlu do allocations
    sol->options.lwork = 0;


    /*-----------------------------------------------------------------------------
     *  create SUPERLU matrix
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&temp);
    ret = mess_matrix_convert(matrix,temp,MESS_CSC);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert);
    if(MESS_IS_REAL(temp)){
        dCreate_CompCol_Matrix(&(sol->A),sol->dim, sol->dim,temp->nnz,temp->values,temp->rowptr,temp->colptr,SLU_NC,SLU_D,SLU_GE);
    }else{
        zCreate_CompCol_Matrix(&(sol->A),sol->dim, sol->dim,temp->nnz,(doublecomplex*)temp->values_cpx,temp->rowptr,temp->colptr,SLU_NC,SLU_Z,SLU_GE);
    }
    //pointers are hold by SuperMatrix A and cleared by SUPERLU Destroy_CompCol_Matrix function//
    temp->values_cpx=NULL;
    temp->values=NULL;
    temp->rowptr=NULL;
    temp->colptr=NULL;

    //get permc
    get_perm_c(ORDERING_MINIMUM_DEGREE_AMD, &(sol->A), sol->options.perm_c);

    /*-----------------------------------------------------------------------------
     *  factorization
     *-----------------------------------------------------------------------------*/
    //factorize
    SuperMatrix X,B;
    if(MESS_IS_REAL(temp)){
        dCreate_Dense_Matrix(&B, temp->rows, 0, NULL, temp->rows,SLU_DN,SLU_D,SLU_GE); //only factorization
        dCreate_Dense_Matrix(&X, temp->rows, 0, NULL, temp->rows,SLU_DN,SLU_D,SLU_GE);
        pdgssvx(sol->options.nprocs,&(sol->options),&(sol->A),sol->perm_c,sol->perm_r, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
            ,&B,&X,&(sol->recip_pivot_growth),&(sol->rcond),NULL,NULL,&(sol->mem_usage),&(sol->info));
        solver->data_type = MESS_REAL;
    }else{
        zCreate_Dense_Matrix(&B, temp->rows, 0, NULL, temp->rows,SLU_DN,SLU_Z,SLU_GE); //only factorization
        zCreate_Dense_Matrix(&X, temp->rows, 0, NULL, temp->rows,SLU_DN,SLU_Z,SLU_GE);
        pzgssvx(sol->options.nprocs,&(sol->options),&(sol->A),sol->perm_c,sol->perm_r, &(sol->equed),sol->R,sol->C, &(sol->L), &(sol->U)
            ,&B,&X,&(sol->recip_pivot_growth),&(sol->rcond),NULL,NULL,&(sol->mem_usage),&(sol->info));
        solver->data_type = MESS_COMPLEX;
    }

    sol->options.fact = FACTORED; //Matrix is now factorized
    sol->options.refact = YES;

    /*-----------------------------------------------------------------------------
     *  set solver functions
     *-----------------------------------------------------------------------------*/
    SET_SOLVERNAME(solver->name, __func__);
    solver->data = (void *)sol;
    solver->solve = superlu_solve;
    solver->solvet = superlu_solvet;
    solver->solveh = superlu_solveh;
    solver->solvem = superlu_solvem;
    solver->solvemt = superlu_solvemt;
    solver->solvemh = superlu_solvemh;
    solver->clear = superlu_clear;
    solver->getL = NULL;
    solver->getU = NULL;
    solver->getpermp = NULL;
    solver->getpermq = NULL;
    solver->inverse = NULL;

    /*-----------------------------------------------------------------------------
     * free memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&(temp));
    Destroy_SuperMatrix_Store(&X);
    Destroy_SuperMatrix_Store(&B);

    return 0;
}




#endif

#endif


