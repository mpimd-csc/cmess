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
 * @file lib/lrcf_adi/sgn/ber_sgn_fac.c
 * @brief Solve a Algebraic Bernoulli equation using the sign function method.
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>
#include <math.h>

/**
 * @brief Solve a standard Bernoulli equation \f$ A^{T} X + X A - XB^{T}BX =0 \f$ via sign function iteration.
 * @param[in] A \f$(n \times n)\f$-matrix
 * @param[in] B \f$(p \times n)\f$-matrix
 * @param[in,out] maxit on input: maximum number of iterations <br>
 *                      on output: number of performed iterations
 * @param[in,out] tol   on input: convergence tolerance <br>
 *                      on output: achieved tolerance
 * @param[in] Z numerical full rank factor \f$ Z \f$ of the solution \f$ X \approx ZZ^T \f$
 * @param[in] scale type of scaling for speedup
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref mess_lrcf_signber function solves a standard algebraic bernoulli equation
 * \f[ A^{T} X + X A + XB^{T} BX = 0 \f]
 * via the sign function iteration (\cite BarBQ07) . The solution is computed as a low rank factor \f$ Z \f$ of the
 * solution \f$ X \approx ZZ^T \f$.
 *
 * The algorithm stops if the maximum iteration number \f$ maxit \f$ is reached or relative Error in Frobenius-norm of the
 * is smaller than \f$ tol \f$
 * \f[ \Vert A_{k+1} - A_{k} \Vert_F < tol \Vert A_{k+1} \Vert_F .\f]
 * At exit \f$ maxit \f$ contains the number of performed iterations and \f$ tol \f$ the achieved tolerance.
 *
 * \attention The function works in dense arithmetics and does not support sparse matrices at all.
 *
 * \sa mess_lrcf_gsignlyap
 * \sa mess_lrcf_signlyap
 * \sa mess_lrcf_gsignber
 * \sa mess_sign_scale_t
 */
int mess_lrcf_signber(mess_matrix A, mess_matrix B, mess_int_t *maxit,double *tol, mess_matrix Z, mess_sign_scale_t scale){

    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t n=0, onemore=0, iter=0, convergence=0, j=0, rrank=0,  *perm, *invperm, rdim=0;
    double d=.0, Err=1.0, nrm1=.0, nrm2=.0,  r_tol, r11= .0, eps = mess_eps(),m,e;
    mess_matrix M1, M2, M3, M4, Binf;
    mess_direct sol;

    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(maxit);
    *maxit = *maxit <=0 ? 50: *maxit;
    mess_check_nullpointer(tol);
    mess_check_positive(*tol);
    mess_check_nullpointer(A); mess_check_dense(A); mess_check_real(A); mess_check_square(A);
    mess_check_nullpointer(B); mess_check_dense(B); mess_check_real(B);
    mess_check_same_rows(A,B);
    mess_check_nullpointer(Z);
    n = A->rows;
    r_tol = n*eps*10;

    /*-----------------------------------------------------------------------------
     * prepare iteration
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&M1, &M2, &M3, &M4, &Binf);
    ret = mess_matrix_alloc(M1, n, n, n*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(M2, n, n, n*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(M3, n, n, n*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(M4, n, n, n*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_copy(B,Binf);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(A,M4);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    mess_try_alloc(perm,mess_int_t * ,sizeof(mess_int_t)*(n));
    mess_try_alloc(invperm,mess_int_t*,sizeof(mess_int_t)*(n));
    for(iter=0;iter<n;++iter){ perm[iter]=0; invperm[iter]=0;}

    /*-----------------------------------------------------------------------------
     * begin iteration
     *-----------------------------------------------------------------------------*/
    for(iter=0;iter<*maxit;++iter){
        if((! convergence)|| ( convergence && (onemore<2))){

            /*-----------------------------------------------------------------------------
             * compute A^(-1) ->M1
             *-----------------------------------------------------------------------------*/
            mess_direct_init(&sol);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
            ret = mess_direct_create_lapack_lu(M4,sol);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_lapack_lu);
            ret = mess_matrix_copy(M4,M3);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            ret = mess_direct_inverse(sol,M1);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_inverse);
            if(Err>1e-1){
                switch (scale){
                    case MESS_SIGN_SCALE_NONE:
                        d=1;
                        break;
                    case MESS_SIGN_SCALE_FRO:
                        ret = mess_matrix_normf(M3,&nrm1);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_normf);
                        ret = mess_matrix_normf(M1,&nrm2);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);
                        d=sqrt(nrm1/nrm2);
                        break;
                    case MESS_SIGN_SCALE_DET:
                        ret = mess_direct_determinant(sol, &m, &e);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_determinant);
                        d = pow(fabs(m),1./((double)n))*pow(2,e/((double) n));
                        break;
                    default:
                        MSG_ERROR("Unknown scaling type %d\n",scale);
                        return MESS_ERROR_ARGUMENTS;
                }
            }else{
                d=1;
            }
            MSG_INFO("it = " MESS_PRINTF_INT " \t d = %lg \n", iter, d);

            /*-----------------------------------------------------------------------------
             * compute (A/d +d*A^(-1))/2=(M3/d +d*M1)/2->M3
             *-----------------------------------------------------------------------------*/
            ret = mess_matrix_add(d/2.0,M1,1.0/(2.0*d),M4);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);

            /*-----------------------------------------------------------------------------
             * compute Binf=[Binf; d*A^-(1)*Binf)]/sqrt(2*d) and column compression
             * via qr decomposition on low rank factor
             *-----------------------------------------------------------------------------*/
            ret = mess_matrix_multiply(MESS_OP_NONE, M1, MESS_OP_NONE, Binf, M2);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_scale(d,M2);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
            ret = mess_matrix_cat(Binf,M2,NULL,NULL,MESS_DENSE,M1);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);
            ret = mess_matrix_scale(1.0/sqrt(2.0*d), M1);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
            ret = mess_matrix_ctranspose(M1, M2);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
            ret = mess_direct_dgeqp3 (M2, NULL, Binf, perm);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_dgeqp3);

            //compute rank r
            r11 = fabs(Binf->values[0])*r_tol;
            rdim = MESS_MIN(Binf->rows, Binf->cols);
            rrank=0;
            for(j=0; j < rdim; ++j){
                if(fabs(Binf->values[j+j*Binf->ld])>r11) ++rrank;
            }
            //compute inverse permutation
            for(j=0; j<n; ++j){invperm[perm[j]]=j;}

            /*-----------------------------------------------------------------------------
             *  compute Binf = Binf(1:r,q)'
             *-----------------------------------------------------------------------------*/
            ret = mess_matrix_rowsub(Binf,0,rrank-1,M2);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            ret = mess_matrix_colperm(M2,invperm);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_colperm);
            ret = mess_matrix_ctranspose(M2,Binf);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_ctranspose);

            /*------------------------------------------------------------------------
             * compute  relative error
             *-----------------------------------------------------------------------*/
            ret = mess_matrix_add(1.0,M4,-1.0,M3);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
            ret = mess_matrix_normf(M3,&nrm1);                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);
            ret = mess_matrix_normf(M4,&nrm2);                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);
            Err = nrm1/nrm2;

            MSG_INFO("Step " MESS_PRINTF_INT "  conv. crit. = %20.15e\n", iter, Err);

            ret = mess_direct_clear(&sol);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_clear);
            if(Err<=*tol){ convergence=1; ++onemore;}else{convergence=0;}
        } else {
            break;
        }
    }

    /*------------------------------------------------------------------------
     * finalize computations
     *-----------------------------------------------------------------------*/
    if(iter==*maxit && Err > *tol){
        MSG_INFO("LYAP_SGN_FAC: No convergence in " MESS_PRINTF_INT  " iterations.\n", *maxit);
        MESS_CLEAR_MATRICES(&M1,&M2,&M3,&M4,&Binf);
        MESS_CLEAR_POINTERS(perm,invperm);
        return MESS_ERROR_CONVERGE;
    }
    *maxit = iter;
    *tol = Err;

    /*-----------------------------------------------------------------------------
     * compute basis von ker(E'-A') via (E-A)P=QR <-> P'(E'-A')Q=R'
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_eye(M1,n,n,MESS_DENSE);                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_eye);
    ret = mess_matrix_add(1.0,M1,-1.0,M4);                                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
    ret = mess_direct_dgeqp3(M4,M1,M2,NULL);                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_dgeqp3);

    r11 = fabs(M2->values[0])*r_tol;
    for(j=0; j < M2->cols; ++j){
        if(fabs(M2->values[j+j*M2->ld])<r11) break;
    }
    if(j==M2->cols){
        MSG_INFO("The Matrix A is probably stable. Lowrank solution factor is empty.");
        ret = mess_matrix_resize(Z,n,0);                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_resize);
        MESS_CLEAR_MATRICES(&M1,&M2,&M3,&M4,&Binf);
        MESS_CLEAR_POINTERS(perm,invperm);
        return 0;
    }

    ret = mess_matrix_colsub(M1,j,M2->cols-1,M3);                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_colsub);
    //M3 is a Basis von ker(E'-A')
    /*-----------------------------------------------------------------------------
     *  compute [Q,R,P]=qr(B'*QY) where QY is a basis of ker(E'-A')
     *-----------------------------------------------------------------------------*/
    //check if rank of Binf is at least dimension of ker(E'-A')
    if(Binf->cols < M3->cols){
        MSG_ERROR("Rank of B is smaller than the number of instable Eigenvalues. System is propably not stabilizable.\n")
            return MESS_ERROR_ARGUMENTS;
    }
    ret = mess_matrix_multiply(MESS_OP_TRANSPOSE,Binf,MESS_OP_NONE,M3,M2);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
    mess_free(perm);
    mess_try_alloc(perm,mess_int_t * ,sizeof(mess_int_t)*(M2->cols));
    ret = mess_direct_dgeqp3(M2,M1,M4,perm);                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);

    /*-----------------------------------------------------------------------------
     *  compute low rank factor
     *-----------------------------------------------------------------------------*/

    r11 = r_tol*fabs(M4->values[0]);
    rrank =0;
    for ( j=0;j<M4->cols ;++j ){
        if(fabs(M4->values[j+j*M4->ld])<r11)break;
        rrank++;
    }
    ret = mess_matrix_colperm(M3,perm);                                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_colperm);
    ret = mess_matrix_ctranspose(M3,M2);                                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    ret = mess_matrix_sub(M4,0,rrank-1,0,rrank-1,M1);                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_sub);
    ret = mess_solver_utsolvem(M1,M2);                                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_solver_utsolvem);
    MESS_MATRIX_RESET(Z);
    ret = mess_matrix_alloc(Z,M2->rows,M2->cols,M2->nnz,M2->store_type,M2->data_type);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
    ret = mess_matrix_ctranspose(M2,Z);                                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    ret = mess_matrix_scale(sqrt(2),Z);                                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
    /*------------------------------------------------------------------------
     * clear storage
     *-----------------------------------------------------------------------*/
    MESS_CLEAR_POINTERS(perm,invperm);
    MESS_CLEAR_MATRICES(&M1, &M2, &M3, &M4, &Binf);
    return (0);
}

/**
 * @brief Solve a generalized Bernoulli equation \f$ A^{T} XE + E^{T}X A - E^{T}XB^{T}BXE =0 \f$ via sign function iteration.
 * @param[in] A \f$(n \times n)\f$-matrix
 * @param[in] E \f$(n \times n)\f$-matrix
 * @param[in] B \f$(p \times n)\f$-matrix
 * @param[in,out] maxit on input: maximum number of iterations <br>
 *                      on output: number of performed iterations
 * @param[in,out] tol   on input: convergence tolerance <br>
 *                      on output: achieved tolerance
 * @param[in] Z numerical full rank factor \f$ Z \f$ of the solution \f$ X \approx ZZ^T \f$
 * @param[in] scale type of scaling for speedup
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref mess_lrcf_gsignber function solves a generalized algebraic bernoulli  equation
 * \f[  A^{T} XE + E^{T}X A - E^{T}XB^{T}BXE =0  \f]
 * via the sign function iteration (\cite BarBQ07) . The solution is computed as a low rank factor \f$ Z \f$ of the
 * solution \f$ X \approx ZZ^T \f$.
 *
 * The algorithm stops if the maximum iteration number \f$ maxit \f$ is reached or relative Error in Frobenius-norm of the
 * is smaller than \f$ tol \f$
 * \f[ \Vert A_{k+1} - A_{k} \Vert_F < tol \Vert A_{k+1} \Vert_F .\f]
 * At exit \f$ maxit \f$ contains the number of performed iterations and \f$ tol \f$ the achieved tolerance.
 *
 * \attention The function works in dense arithmetics and does not support sparse matrices at all.
 *
 * \sa mess_lrcf_signlyap
 * \sa mess_lrcf_gsignlyap
 * \sa mess_lrcf_signber
 * \sa mess_sign_scale_t
 */
int mess_lrcf_gsignber(mess_matrix A,mess_matrix E, mess_matrix B, mess_int_t *maxit,double *tol, mess_matrix Z, mess_sign_scale_t scale){

    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t n=0, onemore=0, iter=0, convergence=0, j=0, rrank=0,  *perm, *invperm, rdim=0;
    double d=.0, Err=1.0, nrm1=.0, nrm2=.0,  r_tol, r11= .0, eps = mess_eps(),m,e;
    mess_matrix M1,M2,M3,M4,Binf;
    mess_direct sol;

    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(maxit);
    *maxit = *maxit <=0 ? 50: *maxit;
    mess_check_nullpointer(tol);
    mess_check_positive(*tol);
    mess_check_nullpointer(A); mess_check_dense(A); mess_check_real(A); mess_check_square(A);
    mess_check_nullpointer(E); mess_check_dense(E); mess_check_real(E); mess_check_square(E);
    mess_check_nullpointer(B); mess_check_dense(B); mess_check_real(B);
    mess_check_same_rows(A,B);
    mess_check_nullpointer(Z);
    n = A->rows;
    r_tol = n*eps*10;

    /*-----------------------------------------------------------------------------
     * prepare iteration
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&M1, &M2, &M3, &M4, &Binf);
    ret = mess_matrix_alloc(M1, n, n, n*n, MESS_DENSE, MESS_REAL);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(M2, n, n, n*n, MESS_DENSE, MESS_REAL);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(M3, n, n, n*n, MESS_DENSE, MESS_REAL);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(M4, n, n, n*n, MESS_DENSE, MESS_REAL);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_copy(B,Binf);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(A,M4);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    mess_try_alloc(perm,mess_int_t * ,sizeof(mess_int_t)*(n));
    mess_try_alloc(invperm,mess_int_t*,sizeof(mess_int_t)*(n));
    for(iter=0;iter<n;++iter){ perm[iter]=0; invperm[iter]=0;}

    /*-----------------------------------------------------------------------------
     *  compute scaling factor for generalized Version
     *-----------------------------------------------------------------------------*/
    //determinant scaling
    double de=0;
    if (scale==MESS_SIGN_SCALE_DET){
        ret = mess_direct_init(&sol);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_direct_init);
        ret = mess_direct_create_lapack_lu(E,sol);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_direct_create_lapack_lu);
        ret = mess_direct_getU(sol,Z);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_getU);
        for(j=0;j<n;++j){
            de+=log(fabs(Z->values[j+(Z->ld)*j]));
        }
        de/=(double)n;
        de = exp(de);
        ret = mess_direct_clear(&sol);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_clear);
    }


    /*-----------------------------------------------------------------------------
     * begin iteration
     *-----------------------------------------------------------------------------*/
    for(iter=0;iter<*maxit;++iter){
        if((! convergence)|| ( convergence && (onemore<2))){

            /*-----------------------------------------------------------------------------
             * compute E*A^(-1) ->M2
             *-----------------------------------------------------------------------------*/
            mess_direct_init(&sol);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
            ret = mess_direct_create_lapack_lu(M4,sol);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_lapack_lu);
            ret = mess_matrix_copy(M4,M3);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            ret = mess_direct_inverse(sol,M1);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_inverse);
            ret = mess_matrix_multiply(MESS_OP_NONE,E,MESS_OP_NONE,M1,M2);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_NONE,M2,MESS_OP_NONE,E,M1);              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);

            if(Err>1e-1){
                switch (scale){
                    case MESS_SIGN_SCALE_NONE:
                        d=1;
                        break;
                    case MESS_SIGN_SCALE_FRO:
                        ret = mess_matrix_normf(M3,&nrm1);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_normf);
                        ret = mess_matrix_normf(M1,&nrm2);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);
                        d=sqrt(nrm1/nrm2);
                        break;
                    case MESS_SIGN_SCALE_DET:
                        ret = mess_direct_determinant(sol, &m, &e);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_determinant);
                        d = pow(fabs(m),1./((double)n))*pow(2,e/((double) n));
                        d /=de;
                        break;
                        /*
                           ret = mess_direct_getU(sol,Z);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_getU);
                           d =0;
                           for(j=0;j<n;++j){
                           d+=log(fabs(Z->values[j+(Z->ld)*j]));
                           }
                           d/=(double)n;
                           d = exp(d)/de;
                           break;
                           */
                    default:
                        MSG_ERROR("Unknown scaling type %d\n",scale);
                        return MESS_ERROR_ARGUMENTS;
                }
            }else{
                d=1;
            }
            //MSG_INFO("it = " MESS_PRINTF_INT " \t d = %lg \n", iter, d);

            /*-----------------------------------------------------------------------------
             * compute (A/d +d*E*A^(-1)E)/2->M1
             *-----------------------------------------------------------------------------*/
            ret = mess_matrix_add(d/2.0,M1,1.0/(2.0*d),M4);                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);

            /*-----------------------------------------------------------------------------
             * compute Binf=[Binf; d*E*A^-(1)*Binf)]/sqrt(2*d) and column compression
             * via qr decomposition on low rank factor
             *-----------------------------------------------------------------------------*/
            ret = mess_matrix_multiply(MESS_OP_NONE, M2, MESS_OP_NONE, Binf, M1);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_scale(d,M1);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
            ret = mess_matrix_cat(Binf,M1,NULL,NULL,MESS_DENSE,M2);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);
            ret = mess_matrix_scale(1.0/sqrt(2.0*d), M2);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
            ret = mess_matrix_ctranspose(M2, M1);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
            ret = mess_direct_dgeqp3 (M1, NULL, Binf, perm);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_dgeqp3);

            //compute rank r
            r11 = fabs(Binf->values[0])*r_tol;
            rdim = MESS_MIN(Binf->rows, Binf->cols);
            rrank=0;
            for(j=0; j < rdim; ++j){
                if(fabs(Binf->values[j+j*Binf->ld])>r11) ++rrank;
            }
            //compute inverse permutation
            for(j=0; j<n; ++j){invperm[perm[j]]=j;}

            /*-----------------------------------------------------------------------------
             *  compute Binf = Binf(1:r,q)'
             *-----------------------------------------------------------------------------*/
            ret = mess_matrix_rowsub(Binf,0,rrank-1,M2);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            ret = mess_matrix_colperm(M2,invperm);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_colperm);
            ret = mess_matrix_ctranspose(M2,Binf);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_ctranspose);

            /*------------------------------------------------------------------------
             * compute  relative error
             *-----------------------------------------------------------------------*/
            ret = mess_matrix_add(1.0,M4,-1.0,M3);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
            ret = mess_matrix_normf(M3,&nrm1);                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);
            ret = mess_matrix_normf(M4,&nrm2);                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);
            Err = nrm1/nrm2;

            MSG_INFO("Step " MESS_PRINTF_INT "  conv. crit. = %20.15e\n", iter, Err);

            ret = mess_direct_clear(&sol);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_clear);
            if(Err<=*tol){ convergence=1; ++onemore;}else{convergence=0;}
        } else {
            break;
        }
    }

    /*------------------------------------------------------------------------
     * finalize computations
     *-----------------------------------------------------------------------*/
    if(iter==*maxit && Err > *tol){
        MSG_INFO("LYAP_SGN_FAC: No convergence in " MESS_PRINTF_INT  " iterations.\n", *maxit);
        MESS_CLEAR_MATRICES(&M1,&M2,&M3,&M4,&Binf);
        MESS_CLEAR_POINTERS(perm,invperm);
        return MESS_ERROR_CONVERGE;
    }
    *maxit = iter;
    *tol = Err;

    /*-----------------------------------------------------------------------------
     * compute basis von ker(E'-A') via (E-A)P=QR <-> P'(E'-A')Q=R'
     *-----------------------------------------------------------------------------*/
    //ret = mess_matrix_eye(M1,n,n,MESS_DENSE);                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_eye);
    ret = mess_matrix_add(1.0,E,-1.0,M4);                                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
    ret = mess_direct_dgeqp3(M4,M1,M2,NULL);                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_dgeqp3);

    r11 = fabs(M2->values[0])*r_tol;
    for(j=0; j < M2->cols; ++j){
        if(fabs(M2->values[j+j*M2->ld])<r11) break;
    }
    if(j==M2->cols){
        MSG_INFO("The Matrix A is probably stable. Lowrank solution factor is empty.");
        ret = mess_matrix_resize(Z,n,0);                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_resize);
        MESS_CLEAR_MATRICES(&M1,&M2,&M3,&M4,&Binf);
        MESS_CLEAR_POINTERS(perm,invperm);
        return 0;
    }

    ret = mess_matrix_colsub(M1,j,M2->cols-1,M3);                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_colsub);
    //M3 is a Basis von ker(E'-A')
    /*-----------------------------------------------------------------------------
     *  compute [Q,R,P]=qr(B'*QY) where QY is a basis of ker(E'-A')
     *-----------------------------------------------------------------------------*/
    //check if rank of Binf is at least dimension of ker(E'-A')
    if(Binf->cols < M3->cols){
        MSG_ERROR("Rank of B is smaller than the number of instable Eigenvalues. System is propably not stabilizable.\n")
            return MESS_ERROR_ARGUMENTS;
    }
    ret = mess_matrix_multiply(MESS_OP_TRANSPOSE,Binf,MESS_OP_NONE,M3,M2);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
    mess_free(perm);
    mess_try_alloc(perm,mess_int_t * ,sizeof(mess_int_t)*(M2->cols));
    ret = mess_direct_dgeqp3(M2,M1,M4,perm);                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);

    //double tempnrm;
    //mess_matrix_normf(M2,&tempnrm);MSG_PRINT("QY=%e\n",tempnrm);
    /*-----------------------------------------------------------------------------
     *  compute low rank factor
     *-----------------------------------------------------------------------------*/

    r11 = r_tol*fabs(M4->values[0]);
    rrank =0;
    for ( j=0;j<M4->cols ;++j ){
        if(fabs(M4->values[j+j*M4->ld])<r11)break;
        rrank++;
    }
    ret = mess_matrix_colperm(M3,perm);                                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_colperm);
    ret = mess_matrix_ctranspose(M3,M2);                                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    ret = mess_matrix_sub(M4,0,rrank-1,0,rrank-1,M1);                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_sub);
    ret = mess_solver_utsolvem(M1,M2);                                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_solver_utsolvem);

    //mess_matrix_normf(M1,&tempnrm);MSG_PRINT("nrm=%e\n",tempnrm);
    //mess_matrix_normf(M2,&tempnrm);MSG_PRINT("nrm=%e\n",tempnrm);

    MESS_MATRIX_RESET(Z);
    ret = mess_matrix_alloc(Z,M2->rows,M2->cols,M2->nnz,M2->store_type,M2->data_type);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
    ret = mess_matrix_ctranspose(M2,Z);                                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    ret = mess_matrix_scale(sqrt(2),Z);                                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
    /*------------------------------------------------------------------------
     * clear storage
     *-----------------------------------------------------------------------*/
    MESS_CLEAR_POINTERS(perm,invperm);
    MESS_CLEAR_MATRICES(&M1, &M2, &M3, &M4, &Binf);

    return (0);
}

