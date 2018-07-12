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
 * @file lib/eigen/eig_hessenberg.c
 * @brief  Eigenvalue computations for upper hessenberg matrices from @lapack.
 * @author @mbehr
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>


/**
 * @brief Compute Schur decomposition of an upoer Hessenberg matrix.
 * @param[in] hess          input upper Hessenberg Matrix
 * @param[out] ev           output vector containing eigenvalues of @p hessA
 * @param[out] Q            output unitary transformation matrix of Schur decomposition
 * @param[out] T            output, (quasi)-upper triangular matrix of Schur decompsotion
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_hessenberg_schur function computes the Schur decomposition of an upper
 * Hessenberg matrix @p hess.
 * If @p hess is @ref MESS_REAL the real Schur decomposition is computed.
 * This means that @p T is quasi upper triangular at return.
 * If @p hess is @ref MESS_COMPLEX the complex Schur decomposition is computed.
 * It holds that \f$ hess = Q T Q^H \f$ and @p Q is orthogonal / unitary respectively.
 * The function calls @c DHSEQR / @c ZHSEQR from @lapack.
 * If you want to get the eigenvector too, use @ref mess_eigen_schur_to_evd with the
 * returned Schur decomposition.
 * If the Schur decomposition is not wanted, call @ref mess_eigen_hessenberg_schur
 * with  @p Q and @p T equals @c NULL.
 *
 * @sa mess_eigen_hessenberg_abs_largest
 * @sa mess_eigen_schur_to_evd
 *
 * @warning It is not checked that @p hess is upper Hessenberg.
 */
int mess_eigen_hessenberg_schur(mess_matrix hess, mess_vector ev, mess_matrix Q, mess_matrix T){
    MSG_FNAME(__func__);

    int ret = 0;
    int schur_wanted = 0;

    char JOB[] = "E", COMPZ[] = "N";
    mess_int_t ILO, IHI;
    mess_int_t lwork = -1;
    mess_int_t info;

    double *WR=NULL, *WI=NULL, *WORK=NULL, *hess_vals=NULL, workspace;
    mess_double_cpx_t *WORK_CPX=NULL, *hess_vals_cpx=NULL, workspace_cpx;

    /*-----------------------------------------------------------------------------
     *  check input and prepare
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(hess);
    mess_check_square(hess);
    mess_check_real_or_complex(hess);
    mess_check_dense(hess);
    mess_check_nullpointer(ev);

    if((!Q && T) || (Q && !T)){
        MSG_ERROR("Q points to NULL and T not or vice versa. Either both matrices are needed or none of them.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  prepeare dhesqr call
     *-----------------------------------------------------------------------------*/
    if(Q && T){
        schur_wanted = 1;
        JOB[0] = 'S';
        COMPZ[0] = 'I';

        //matrix T contains schur form of hessenberg
        MESS_MATRIX_RESET(T);
        ret = mess_matrix_copy(hess, T);                                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);

        //matrix Q orthogonal matrix for schurtransformation of hessenberg
        MESS_MATRIX_RESET(Q);
        ret = mess_matrix_alloc(Q, hess->rows, hess->cols, hess->nnz, MESS_DENSE, hess->data_type);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    }else{
        schur_wanted = 0;
        JOB[0] = 'E';
        COMPZ[0] = 'N';

        //hessenberg matrix is overwritten by dhseqr
        if(MESS_IS_REAL(hess)){
            mess_try_alloc(hess_vals, double*, sizeof(double)*hess->ld*hess->cols);
            memcpy(hess_vals, hess->values, sizeof(double)*hess->ld*hess->cols);
        }else{
            mess_try_alloc(hess_vals_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*hess->ld*hess->cols);
            memcpy(hess_vals_cpx, hess->values_cpx, sizeof(mess_double_cpx_t)*hess->ld*hess->cols);
        }
    }

    ILO = 1;
    IHI = hess->rows;

    /*-----------------------------------------------------------------------------
     *  workspace query and call
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(hess)){

        //prepare memory for eigenvalues
        mess_try_alloc(WR, double*, sizeof(double)*hess->rows);
        mess_try_alloc(WI, double*, sizeof(double)*hess->rows);

        // workspace query
        F77_GLOBAL(dhseqr,DHSEQR)(JOB, COMPZ, &(hess->rows), &ILO, &IHI, schur_wanted? T->values:hess_vals, schur_wanted? &(T->ld): &(hess->ld), WR, WI, schur_wanted?Q->values:NULL, schur_wanted?&(Q->ld):&(hess->ld)/*dummy*/, &workspace, &lwork, &info);
        if(info){goto clear;}

        //workspace alloc
        lwork = nearbyint(workspace+1);
        mess_try_alloc(WORK, double*, sizeof(double)*lwork);

        // call
        F77_GLOBAL(dhseqr,DHSEQR)(JOB, COMPZ, &(hess->rows), &ILO, &IHI, schur_wanted? T->values:hess_vals, schur_wanted? &(T->ld): &(hess->ld), WR, WI, schur_wanted?Q->values:NULL, schur_wanted?&(Q->ld):&(hess->ld)/*dummy*/, WORK, &lwork, &info);
        if(info){goto clear;}

        // get vector from lapack data
        ret =  mess_vector_from_lapack(ev, hess->rows, WR, WI);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_from_lapack);

    }else{

        // prepare memory for eigenvalues
        ret = mess_vector_reset(ev);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_reset);
        ret = mess_vector_alloc(ev, hess->rows, MESS_COMPLEX);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);

        // workspace query
        F77_GLOBAL(zhseqr,ZHSEQR)(JOB, COMPZ, &(hess->rows), &ILO, &IHI, schur_wanted? T->values_cpx:hess_vals_cpx, schur_wanted? &(T->ld): &(hess->ld), ev->values_cpx, schur_wanted?Q->values_cpx:NULL, schur_wanted?&(Q->ld):&(hess->ld)/*dummy*/, &workspace_cpx, &lwork, &info);
        if(info){goto clear;}

        //workspace alloc
        lwork = nearbyint(creal(workspace_cpx)+1);
        mess_try_alloc(WORK_CPX, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*lwork);

        // call
        F77_GLOBAL(zhseqr,ZHSEQR)(JOB, COMPZ, &(hess->rows), &ILO, &IHI, schur_wanted? T->values_cpx:hess_vals_cpx, schur_wanted? &(T->ld): &(hess->ld), ev->values_cpx, schur_wanted?Q->values_cpx:NULL, schur_wanted?&(Q->ld):&(hess->ld)/*dummy*/, WORK_CPX, &lwork, &info);
        if(info){goto clear;}

    }


    /*-----------------------------------------------------------------------------
     *  clear stuff
     *-----------------------------------------------------------------------------*/
clear:
    if(WR)              mess_free(WR);
    if(WI)              mess_free(WI);
    if(WORK)            mess_free(WORK);
    if(WORK_CPX)        mess_free(WORK_CPX);
    if(hess_vals)       mess_free(hess_vals);
    if(hess_vals_cpx)   mess_free(hess_vals_cpx);


    /*-----------------------------------------------------------------------------
     *  error if info !=0
     *-----------------------------------------------------------------------------*/
    if(info){
        MSG_ERROR("dhseqr/zhseqr returned with error = " MESS_PRINTF_INT "\n", info);
        return MESS_ERROR_LAPACK;
    }

    return 0;
}


/**
 * @brief Compute the largest eigenvalue in magnitude and corresponding eigenvector of a upper Hessenberg matrix,
 * @param[in] hess          input upper Hessenberg Matrix
 * @param[out] evector      output normalized eigenvector corresponding to largest eigenvalue
 * @param[out] largest_ev   output largest eigenvalue in magnitude of @p hess
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_hessenberg_abs_largest function computes the largest
 * eigenvalue in magnitude and a correspnding eigenvector of a
 * upper Hessenberg matrix @p hess.
 *
 * @sa mess_eigen_hessenberg_abs_largest
 * @sa mess_eigen_schur_to_evd
 *
 * @warning It is not checked that @p hess is upper Hessenberg.
 */
int mess_eigen_hessenberg_abs_largest(mess_matrix hess, mess_vector evector, mess_double_cpx_t *largest_ev){
    MSG_FNAME(__func__);

    int ret = 0;
    char SIDE[] = "R", HOWMNY[] = "S";
    MESS_LAPACK_LOGICAL* SELECT=NULL;
    mess_int_t maxind, i, MM=0, M=0, info;
    double dummy, *VALS_EVEC=NULL, *WORK=NULL;
    mess_double_cpx_t *VALS_EVEC_CPX=NULL, *WORK_CPX=NULL;
    mess_vector evals;
    mess_matrix Q, T;

    /*-----------------------------------------------------------------------------
     *  check input and prepare
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(hess);
    mess_check_square(hess);
    mess_check_real_or_complex(hess);
    mess_check_dense(hess);
    mess_check_nullpointer(evector);
    mess_check_real_or_complex(evector);
    mess_check_nullpointer(largest_ev);

    mess_vector_reset(evector);
    mess_vector_alloc(evector, hess->rows, hess->data_type);
    /*-----------------------------------------------------------------------------
     *  compute schur decomposition
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&Q,&T);
    MESS_INIT_VECTORS(&evals);

    ret = mess_eigen_hessenberg_schur(hess, evals, Q, T);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_hessenberg_schur);

    /*-----------------------------------------------------------------------------
     *  compute index of largest eigenvalue in magnitude and get largest eigenvalue
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_max(evals, &dummy, &maxind);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_max);
    ret = mess_vector_get(evals, maxind, largest_ev);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_get);

    //mess_vector_print(evals);
    //printf("index       = %d\n", maxind);
    //printf("largest ev  = %e+%eI\n", creal(*largest_ev), cimag(*largest_ev));

    /*-----------------------------------------------------------------------------
     *  prepare call of dtrevc/ztrevc
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(SELECT, MESS_LAPACK_LOGICAL*, sizeof(MESS_LAPACK_LOGICAL)*hess->rows);
    for(i=0;i<hess->rows;++i){SELECT[i]=0;}
    SELECT[maxind]=1;

    /*-----------------------------------------------------------------------------
     *  compute eigenvector of largest eigenvalue in magnitude
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(T)){
        //alloc workspace for dtrvec
        mess_try_alloc(WORK, double*, 3*sizeof(double)*T->rows);
        MM=2;
        mess_try_alloc(VALS_EVEC, double*, MM*sizeof(double)*T->rows);

        //call
        F77_GLOBAL(dtrevc,DTREVC)(SIDE, HOWMNY, SELECT, &(T->rows), T->values, &(T->ld), NULL, &(T->rows) /*dummy*/, VALS_EVEC, &(T->rows), &MM, &M, WORK, &info);
        if(info){goto clear;}
    }else{
        //alloc workspace for ztrvec
        mess_try_alloc(WORK_CPX, mess_double_cpx_t*, 2*sizeof(mess_double_cpx_t)*T->rows);
        mess_try_alloc(WORK, double*, sizeof(double)*T->rows);
        MM=1;
        mess_try_alloc(VALS_EVEC_CPX, mess_double_cpx_t*, MM*sizeof(mess_double_cpx_t)*T->rows);

        //call
        F77_GLOBAL(ztrevc,ZTREVC)(SIDE, HOWMNY, SELECT, &(T->rows), T->values_cpx, &(T->ld), NULL, &(T->rows) /*dummy*/, VALS_EVEC_CPX, &(T->rows), &MM, &M, WORK_CPX, WORK, &info);
        if(info){goto clear;}
    }

    /*-----------------------------------------------------------------------------
     * get eigenvector of schur matrix, we use the evals vector.
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(T)){
        if(M==1){
            ret = mess_vector_from_lapack(evals, T->rows, VALS_EVEC, NULL);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_from_lapack);
        }else{
            ret = mess_vector_from_lapack(evals, T->rows, VALS_EVEC, VALS_EVEC+T->rows);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_from_lapack);
        }
    }else{
        ret = mess_vector_from_farray(evals, T->rows, NULL, VALS_EVEC_CPX);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_from_farray);
    }

    /*-----------------------------------------------------------------------------
     *  transform eigenvector from T to eigenvector of hessenberg matrix
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_mvp(MESS_OP_NONE, Q, evals, evector);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

    /*-----------------------------------------------------------------------------
     *  normalize vector
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_norm2(evector, &dummy);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_norm2);
    ret = mess_vector_scale(1.0/dummy, evector);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_scale);

    /*-----------------------------------------------------------------------------
     *  clear stuff
     *-----------------------------------------------------------------------------*/
clear:
    if(SELECT)          mess_free(SELECT);
    if(VALS_EVEC)       mess_free(VALS_EVEC);
    if(VALS_EVEC_CPX)   mess_free(VALS_EVEC_CPX);
    if(WORK)            mess_free(WORK);
    if(WORK_CPX)        mess_free(WORK_CPX);

    MESS_CLEAR_MATRICES(&Q,&T);
    MESS_CLEAR_VECTORS(&evals);


    /*-----------------------------------------------------------------------------
     *  error if info !=0
     *-----------------------------------------------------------------------------*/
    if(info){
        MSG_ERROR("dtrevc/ztrevc returned with error = " MESS_PRINTF_INT "\n", info);
        return MESS_ERROR_LAPACK;
    }

    return 0;
}




