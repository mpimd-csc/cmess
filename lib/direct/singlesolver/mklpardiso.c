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
 * @file lib/direct/singlesolver/mklpardiso.c
 * @brief Interface to @mklpardiso
 * @author @mbehr
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//define for 32 bit pardiso call
#ifdef PARDISO_32

void pardiso(void* pt, mess_int_t* maxfct, mess_int_t* mnum, mess_int_t* mtype, mess_int_t* phase, mess_int_t* n,
        void* a, mess_int_t* ia, mess_int_t* ja, mess_int_t* idum, mess_int_t* nrhs, mess_int_t* iparm,
        mess_int_t* msglevel, void* b, void* x, mess_int_t* error);
#undef MKLPARDISO
#define MKLPARDISO pardiso

#endif

// if available take 64 bit pardiso call, note you need 64 bit numbers to make this correctly
#if defined(PARDISO_64) && defined(MESS64)

void pardiso_64(void* pt, mess_int_t* maxfct, mess_int_t* mnum, mess_int_t* mtype, mess_int_t* phase, mess_int_t* n,
        void* a, mess_int_t* ia, mess_int_t* ja, mess_int_t* idum, mess_int_t* nrhs, mess_int_t* iparm,
        mess_int_t* msglevel, void* b, void* x, mess_int_t* error);

#undef MKLPARDISO
#define MKLPARDISO pardiso_64

#endif


typedef struct mklpardiso_solver {
    mess_matrix work;
    mess_int_t mtype;
    void *pt[64];
    mess_int_t iparm[64];
    mess_int_t solver;
    mess_int_t maxfct;
    mess_int_t mnum;
    mess_int_t phase;
    mess_int_t error;
    mess_int_t message_level;
    mess_int_t num_procs;
} mklpardiso_solver;



static mess_int_t mklpardiso_error(mess_int_t error){
    if(error==0) return 0;
    MSG_FNAME(__func__);
    switch (error){
        case -1:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT ", Input inconsistent.\n",error); break;
        case -2:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT ", Not enought memory.\n",error); break;
        case -3:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT ", Reordering problem.\n",error); break;
        case -4:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT ", Zero pivot, numerical fact. or iterative refinement problem.\n",error); break;
        case -5:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT ", Unclassified (internal) error.\n",error); break;
        case -6:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT ", Reordering failed.\n",error); break;
        case -7:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT ", Diagonal matrix problem.\n",error); break;
        case -8:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT ", 32-bit integer overflow problem.\n",error); break;
        case -9:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT ", Not enough memory for OOC.\n",error); break;
        case -10:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT ", Error opening OOC files.\n", error); break;
        case -11:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT ", Read/Write Error with OOC files.\n",error); break;
        case -12:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT ", pardiso_64 called from 32-bit library.\n",error); break;
        default:
            MSG_ERROR("MKL-PARDISO Error:" MESS_PRINTF_INT "\n",error); break;
    }
    return MESS_ERROR_MKLPARDISO;
}


static int mklpardiso_solve_vector(void*data, mess_int_t op, mess_vector b, mess_vector x){
    MSG_FNAME(__func__);
    mklpardiso_solver * sol = ( mklpardiso_solver *)data;
    int ret = 0;
    mess_int_t idum, one=1;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->dim != sol->work->rows) {
        MSG_ERROR("b has the wrong dimension (b->dim = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->dim, sol->work->rows);
        return MESS_ERROR_DIMENSION;
    }
    if ( x->dim != sol->work->rows){
        ret = mess_vector_resize(x, sol->work->rows); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }

    /*-----------------------------------------------------------------------------
     *  solve
     *-----------------------------------------------------------------------------*/
    //set transposed or nontransposed
    sol->iparm[11]=op;

    if(MESS_IS_REAL(sol->work) && MESS_IS_REAL(b)){
        ret = mess_vector_toreal(x);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal);
        MKLPARDISO(sol->pt,&(sol->maxfct),&(sol->mnum),&(sol->mtype),&(sol->phase),&(sol->work->rows),
                sol->work->values,sol->work->rowptr,sol->work->colptr,&idum,&one,
                sol->iparm,&(sol->message_level),b->values,x->values,&(sol->error));
    }else if(MESS_IS_REAL(sol->work) && MESS_IS_COMPLEX(b)){
        mess_vector xr,xi,br,bi;
        mess_vector_init(&xr,b->dim, MESS_REAL);
        mess_vector_init(&xi,b->dim, MESS_REAL);
        mess_vector_init(&br,b->dim, MESS_REAL);
        mess_vector_init(&bi,b->dim, MESS_REAL);
        mess_vector_realpart(b,br);
        mess_vector_imagpart(b,bi);

        //solve for real part
        MKLPARDISO(sol->pt,&(sol->maxfct),&(sol->mnum),&(sol->mtype),&(sol->phase),&(sol->work->rows),
                sol->work->values,sol->work->rowptr,sol->work->colptr,&idum,&one,
                sol->iparm,&(sol->message_level),br->values,xr->values,&(sol->error));

        //solve for imag part
        MKLPARDISO(sol->pt,&(sol->maxfct),&(sol->mnum),&(sol->mtype),&(sol->phase),&(sol->work->rows),
                sol->work->values,sol->work->rowptr,sol->work->colptr,&idum,&one,
                sol->iparm,&(sol->message_level),bi->values,xi->values,&(sol->error));

        //get solution
        mess_vector_complex_from_parts(xr,xi,x);

        //clear vectors
        MESS_CLEAR_VECTORS(&xr,&xi,&br,&bi);
    }else{
        if(MESS_IS_REAL(x)){mess_vector_tocomplex(x);}

        //solve
        MKLPARDISO(sol->pt,&(sol->maxfct),&(sol->mnum),&(sol->mtype),&(sol->phase),&(sol->work->rows),
                sol->work->values_cpx,sol->work->rowptr,sol->work->colptr,&idum,&one,
                sol->iparm,&(sol->message_level),b->values_cpx,x->values_cpx,&(sol->error));
    }

    ret = mklpardiso_error(sol->error);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),pardiso);

    return ret;
}

static int mklpardiso_solve_matrix(void*data, mess_int_t op, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    mklpardiso_solver * sol = ( mklpardiso_solver *)data;
    int ret = 0;
    mess_int_t idum, conv=0;
    mess_matrix bwork;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(b);
    mess_check_nullpointer(x);

    if ( b->rows != sol->work->rows) {
        MSG_ERROR("b has the wrong dimension (b->rows = " MESS_PRINTF_INT ", solver->dim = " MESS_PRINTF_INT ") \n", b->rows, sol->work->rows);
        return MESS_ERROR_DIMENSION;
    }
    if ( x->rows != sol->work->rows){
        ret = mess_matrix_resize(x, sol->work->rows, sol->work->cols); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_resize);
    }

    MESS_MATRIX_CHECKFORMAT(b, bwork, conv, MESS_DENSE);
    MESS_MATRIX_RESET(x);
    ret = mess_matrix_alloc(x, bwork->rows, bwork->cols,bwork->rows*bwork->cols, MESS_DENSE, sol->work->data_type );


    /*-----------------------------------------------------------------------------
     *  solve
     *-----------------------------------------------------------------------------*/
    //set transposed or nontransposed
    sol->iparm[11]  = op;

    if(MESS_IS_REAL(sol->work) && MESS_IS_REAL(bwork)){
        ret = mess_matrix_toreal(x);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_toreal);

        MKLPARDISO(sol->pt,&(sol->maxfct),&(sol->mnum),&(sol->mtype),&(sol->phase),&(sol->work->rows),
                sol->work->values,sol->work->rowptr,sol->work->colptr,&idum,&(bwork->cols),
                sol->iparm,&(sol->message_level),bwork->values,x->values,&(sol->error));

    }else if(MESS_IS_REAL(sol->work) && MESS_IS_COMPLEX(bwork)){
        mess_matrix xr,xi,br,bi;
        MESS_INIT_MATRICES(&xr,&xi,&br,&bi);
        mess_matrix_realpart(bwork,br); mess_matrix_totype(br,MESS_DENSE);
        mess_matrix_imagpart(bwork,bi); mess_matrix_totype(br,MESS_DENSE);
        mess_matrix_alloc(xr,br->rows,br->cols,br->nnz,MESS_DENSE,MESS_REAL);
        mess_matrix_alloc(xi,bi->rows,bi->cols,bi->nnz,MESS_DENSE,MESS_REAL);

        //solve for real part
        MKLPARDISO(sol->pt,&(sol->maxfct),&(sol->mnum),&(sol->mtype),&(sol->phase),&(sol->work->rows),
                sol->work->values,sol->work->rowptr,sol->work->colptr,&idum,&(br->cols),
                sol->iparm,&(sol->message_level),br->values,xr->values,&(sol->error));

        //solve for imag part:Y
        MKLPARDISO(sol->pt,&(sol->maxfct),&(sol->mnum),&(sol->mtype),&(sol->phase),&(sol->work->rows),
                sol->work->values,sol->work->rowptr,sol->work->colptr,&idum,&(bi->cols),
                sol->iparm,&(sol->message_level),bi->values,xi->values,&(sol->error));

        //get solution
        mess_matrix_complex_from_parts(xr,xi,x);

        //clear matrix
        MESS_CLEAR_MATRICES(&xr,&xi,&br,&bi);

    }else{
        if(MESS_IS_REAL(x)){mess_matrix_tocomplex(x);}

        //solve
        MKLPARDISO(sol->pt,&(sol->maxfct),&(sol->mnum),&(sol->mtype),&(sol->phase),&(sol->work->rows),
                sol->work->values_cpx,sol->work->rowptr,sol->work->colptr,&idum,&(bwork->cols),
                sol->iparm,&(sol->message_level),bwork->values_cpx,x->values_cpx,&(sol->error));

    }
    ret = mklpardiso_error(sol->error);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),pardiso);

    if ( conv == 0) {
        mess_matrix_clear(&bwork);
    }



    return ret;
}


/**
 * @brief Solve \f$ Ax=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The  @ref mklpardiso_solve function solves \f$ Ax=b \f$.
 *
 */
static int mklpardiso_solve(void*data, mess_vector b, mess_vector x){
    return mklpardiso_solve_vector(data, 0, b, x);
}

/**
 * @brief Solve \f$ A^Tx=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref mklpardiso_solvet function solves \f$ A^Tx=b \f$.
 *
 */
static int mklpardiso_solvet(void*data, mess_vector b, mess_vector x){
    return mklpardiso_solve_vector(data, 2, b, x);
}

/**
 * @brief Solve \f$ A^Hx=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref mklpardiso_solveh function solves \f$ A^Hx=b \f$.
 *
 */
static int mklpardiso_solveh(void*data, mess_vector b, mess_vector x){
    return mklpardiso_solve_vector(data, 1, b, x);
}

/**
 * @brief Solve \f$ AX=B \f$ (matrix version)
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref mklpardiso_solvem function solves \f$ Ax=b \f$.
 *
 */
static int mklpardiso_solvem(void*data, mess_matrix b, mess_matrix x){
    return mklpardiso_solve_matrix(data, 0, b, x);
}

/**
 * @brief Solve \f$ A^Tx=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref mklpardiso_solvemt function solves \f$ A^Tx=b \f$.
 *
 */
static int mklpardiso_solvemt(void*data, mess_matrix b, mess_matrix x){
    return mklpardiso_solve_matrix(data, 2, b, x);
}

/**
 * @brief Solve \f$ A^Hx=b \f$
 * @param [in] data pointer to the data object
 * @param [in] b right hand side
 * @param [in,out] x solution
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref mklpardiso_solvemh function solves \f$ A^Hx=b \f$.
 *
 */
static int mklpardiso_solvemh(void*data, mess_matrix b, mess_matrix x){
    return mklpardiso_solve_matrix(data, 1, b, x);
}

/**
 * @brief Clear an @ref mklpardiso_solver object
 * @param [in,out] solver pointer to the data object
 * @return always zero
 *
 * The @ref mklpardiso_clear function clears an @ref mklpardiso_solver object.
 *
 */
static int mklpardiso_clear(void *solver){
    //MSG_FNAME(__func__);
    mess_int_t dummy=0;
    mklpardiso_solver * sol = ( mklpardiso_solver*) solver;
    if ( sol != NULL) {
        sol->phase = -1;
        if(MESS_IS_REAL(sol->work)){
            MKLPARDISO(sol->pt,&(sol->maxfct),&(sol->mnum),&(sol->mtype),&(sol->phase),&(sol->work->rows),
                    sol->work->values,sol->work->rowptr,sol->work->colptr,&dummy,&dummy,
                    sol->iparm,&(sol->message_level),NULL,NULL,&(sol->error));
        }else{
            MKLPARDISO(sol->pt,&(sol->maxfct),&(sol->mnum),&(sol->mtype),&(sol->phase),&(sol->work->rows),
                    sol->work->values,sol->work->rowptr,sol->work->colptr,&dummy,&dummy,
                    sol->iparm,&(sol->message_level),NULL,NULL,&(sol->error));
        }
        mess_matrix_clear(&sol->work);
        mess_free( sol );
    }

    return 0;
}

/**
 * @brief Generate a direct linear solver for standard linear systems \f$ Ax=b \f$ with @mklpardiso.
 * @param[in] matrix            input matrix to decompose
 * @param[in,out] direct        input/output direct solver
 * @param[in] iparm             input array of solver parameters
 * @param[in] mtype             input matrix type
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref  mess_direct_create_mklpardiso_control function creates a @ref mess_direct solver using @mklpardiso
 * with specified solver parameters in @p iparm and matrix in type specified with @p mtype.
 *
 */
int mess_direct_create_mklpardiso_control(mess_matrix matrix, mess_direct direct, mess_int_t*iparm, mess_int_t mtype){
    MSG_FNAME(__func__);
    mklpardiso_solver * data;
    mess_int_t ret, i, nrhs=0, idum;

    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(direct);
    mess_check_nullpointer(iparm);
    mess_check_square(matrix);
    mess_check_real_or_complex(matrix);

    /*-----------------------------------------------------------------------------
     *  init solver structure and  paradiso
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data,  mklpardiso_solver*, sizeof(mklpardiso_solver));

    //convert to csr matrix
    MESS_INIT_MATRICES(&(data->work));
    ret = mess_matrix_convert(matrix,data->work,MESS_CSR);                                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert);

    /*-----------------------------------------------------------------------------
     *  setup pardiso control parameters
     *-----------------------------------------------------------------------------*/
    for(i=0;i<64;++i){
        data->iparm[i] = iparm[i];
    }

    data->maxfct        = 1;    /* Number of numerical factorizations */
    data->mnum          = 1;    /* Which factorization to use. */
    data->error         = 0;    /* Initialize error flag.  */
    data->message_level = 0;    /* Print statistical informations. */
    data->mtype         = mtype;

    /*-----------------------------------------------------------------------------
     *  Initialize the internal solver memory pointer. This is only necessary
     *  for the FIRST call of the PARDISO solver.
     *-----------------------------------------------------------------------------*/
    for(i=0;i<64;i++){
        data->pt[i]=0;
    }

    /*-----------------------------------------------------------------------------
     *  reordering and symbolic factorization
     *-----------------------------------------------------------------------------*/
    data->phase = 12;

    if(MESS_IS_REAL(data->work)){
        MKLPARDISO(data->pt, &(data->maxfct), &(data->mnum), &(data->mtype), &(data->phase),
                &(data->work->rows), data->work->values, data->work->rowptr, data->work->colptr,&idum,&nrhs,
                data->iparm,&(data->message_level),NULL,NULL,&(data->error));
    }else{
        MKLPARDISO(data->pt, &(data->maxfct), &(data->mnum), &(data->mtype), &(data->phase),
                &(data->work->rows), data->work->values_cpx, data->work->rowptr, data->work->colptr,&idum,&nrhs,
                data->iparm,&(data->message_level),NULL,NULL,&(data->error));
    }
    ret = mklpardiso_error(data->error);                                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_mklparadiso);

    /*-----------------------------------------------------------------------------
     *  set solver
     *-----------------------------------------------------------------------------*/
    data->phase         = 33;
    SET_SOLVERNAME(direct->name, __func__);
    direct->data_type   = matrix->data_type;
    direct->data        = (void *)data;
    direct->rows        = matrix->rows;
    direct->cols        = matrix->cols;
    direct->solve       = mklpardiso_solve;
    direct->solvet      = mklpardiso_solvet;
    direct->solveh      = mklpardiso_solveh;
    direct->solvem      = mklpardiso_solvem;
    direct->solvemt     = mklpardiso_solvemt;
    direct->solvemh     = mklpardiso_solvemh;
    direct->clear       = mklpardiso_clear;
    direct->getL        = NULL;
    direct->getU        = NULL;
    direct->getpermp    = NULL;
    direct->getpermq    = NULL;
    direct->inverse     = NULL;
    direct->getscalerow = NULL;
    direct->inverse     = NULL;
    direct->det         = NULL;
    direct->detc        = NULL;



    return 0;
}


/**
 * @brief Generate a direct linear solver for standard linear systems \f$ Ax=b \f$ with @mklpardiso.
 * @param[in] matrix            matrix to decompose
 * @param[in,out] direct        direct solver
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref  mess_direct_create_mklpardiso function creates a @ref mess_direct solver using @mklpardiso.
 *
 */
int mess_direct_create_mklpardiso(mess_matrix matrix, mess_direct direct){
    MSG_FNAME(__func__);
    mess_int_t iparm[64];
    mess_int_t mtype;
    mess_int_t ret =0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(direct);

    /*-----------------------------------------------------------------------------
     *  iparm array and mtype
     *-----------------------------------------------------------------------------*/

    iparm[0]      = 1;    /* No solver default */
    iparm[1]      = 2;    /* Fill-in reordering from METIS */
    iparm[2]      = 0;    /* Reserved set to 0 */
    iparm[3]      = 0;    /* No iterative-direct algorithm */
    iparm[4]      = 0;    /* No user fill-in reducing permutation */
    iparm[5]      = 0;    /* Write solution into x */
    iparm[6]      = 0;    /* Not in use */
    iparm[7]      = 2;    /* Max number of iterative refinement steps */
    iparm[8]      = 0;    /* Not in use */
    iparm[9]      = 13;   /* Perturb the pivot elements with 1E-13 */
    iparm[10]     = 1;    /* Use nonsymmetric permutation and scaling MPS */
    iparm[11]     = 0;    /* Conjugate/transpose solve */
    iparm[12]     = 1;    /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13]     = 0;    /* Output: Number of perturbed pivots */
    iparm[14]     = 0;    /* Not in use */
    iparm[15]     = 0;    /* Not in use */
    iparm[16]     = 0;    /* Not in use */
    iparm[17]     = -1;   /* Output: Number of nonzeros in the factor LU */
    iparm[18]     = -1;   /* Output: Mflops for LU factorization */
    iparm[19]     = 0;    /* Output: Numbers of CG Iterations */
    iparm[20]     = 0;    /* Pivoting for symmetric matrices */
    iparm[21]     = 0;    /* Output number of positive eigenvalues, for symmetric indefinite matrices */
    iparm[22]     = 0;    /* Output number of negativ  eigenvalues, for symmetric indefinite matrices */
    iparm[23]     = 0;    /* Parallel Factorization Control */
    iparm[24]     = 0;    /* Parallel forward backward solve */
    iparm[25]     = 0;    /* Reserved: Set to zero */
    iparm[26]     = 1;    /* Matrix checker */
    iparm[27]     = 0;    /* Single or Double Precision */
    iparm[28]     = 0;    /* Reserved: Set to zero */
    iparm[29]     = 0;    /* Number of zero or negative pivots */
    iparm[30]     = 0;    /* Partial Solve */
    iparm[31]     = 0;    /* Reserved: Set to zero */
    iparm[32]     = 0;    /* Reserved: Set to zero */
    iparm[33]     = 0;    /* Optimal number of OpenMP threads for conditional numerical reproducibility. */
    iparm[34]     = 1;    /* One or zero based indexing */
    iparm[35]     = 0;    /* Schur complement matrix computation */
    iparm[36]     = 0;    /* Format for storage */
    iparm[37]     = 0;    /* Reserved: Set to zero */
    iparm[38]     = 0;    /* Reserved: Set to zero */
    iparm[39]     = 0;    /* Reserved: Set to zero */
    iparm[40]     = 0;    /* Reserved: Set to zero */
    iparm[41]     = 0;    /* Reserved: Set to zero */
    iparm[42]     = 0;    /* Reserved: Set to zero */
    iparm[43]     = 0;    /* Reserved: Set to zero */
    iparm[44]     = 0;    /* Reserved: Set to zero */
    iparm[45]     = 0;    /* Reserved: Set to zero */
    iparm[46]     = 0;    /* Reserved: Set to zero */
    iparm[47]     = 0;    /* Reserved: Set to zero */
    iparm[48]     = 0;    /* Reserved: Set to zero */
    iparm[49]     = 0;    /* Reserved: Set to zero */
    iparm[50]     = 0;    /* Reserved: Set to zero */
    iparm[51]     = 0;    /* Reserved: Set to zero */
    iparm[52]     = 0;    /* Reserved: Set to zero */
    iparm[53]     = 0;    /* Reserved: Set to zero */
    iparm[54]     = 0;    /* Reserved: Set to zero */
    iparm[55]     = 0;    /* Diagonal and pivoting control */
    iparm[56]     = 0;    /* Reserved: Set to zero */
    iparm[57]     = 0;    /* Reserved: Set to zero */
    iparm[58]     = 0;    /* Reserved: Set to zero */
    iparm[59]     = 0;    /* In Core mode */
    iparm[60]     = 0;    /* Reserved: Set to zero */
    iparm[61]     = 0;    /* Reserved: Set to zero */
    iparm[62]     = 0;    /* Output size of the minimum ooc memory for numerical factorization and solution */
    iparm[63]     = 0;    /* Reserved: Set to zero */

    mtype         = MESS_IS_REAL(matrix)? 11: 13;    /*Real unsymmetric or Complex Unsymmetric Matrix */

    /*-----------------------------------------------------------------------------
     *  create solver
     *-----------------------------------------------------------------------------*/

     ret = mess_direct_create_mklpardiso_control(matrix, direct, iparm, mtype);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_create_mklpardiso_control);

     //overwrite solver name
    if(direct->name){
        MESS_CLEAR_POINTERS(direct->name);
        SET_SOLVERNAME(direct->name, __func__);
    }

    return 0;
}


