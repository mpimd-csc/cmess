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
 * @file lib/lrcf_adi/sgn/lyap_sgn_fac.c
 * @brief Solve a Lyapunov Equation using the sign function method.
 * @author @koehlerm
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
 * @brief Solve a generalized  Lyapunov Equation \f$ A^{T} X E + E^{T} X A + B^{T} B = 0 \f$ via sign function iteration.
 * @param[in] A input \f$ (n \times n) \f$-matrix
 * @param[in] E input \f$ (n \times n) \f$-matrix
 * @param[in] B input \f$ (p \times n) \f$-matrix
 * @param[in,out] maxit on input: maximum number of iterations \n
 *                      on output: number of performed iterations
 * @param[in,out] tol   on input: convergence tolerance \n
 *                      on output: achieved tolerance
 * @param[out] Z numerical full rank factor \f$ Z \f$ of the solution \f$ X \approx ZZ^T \f$
 * @param[out] scale type of scaling for speedup
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcf_gsignfac function solves a generalized Lyapunov Equation
 * \f[A^{T} X E + E^{T} X A + B^{T} B = 0 \f]
 * via the sign function iteration (\cite BenQ99) . The solution is computed as a low rank factor \f$ Z \f$ of the solution
 * \f$ X \approx ZZ^T \f$.
 *
 * The algorithm stops if the maximum iteration number \f$ maxit \f$ is reached or if the Frobenius-norm of the current
 * iterate is smaller than \f$ tol \f$
 * \f[ \Vert A_i + E \Vert_F < tol .\f]
 * At exit \f$ maxit \f$ contains the number of performed iterations and \f$ tol \f$ the achieved tolerance.
 *
 * \attention The function works in dense arithmetics and does not support sparse matrices at all.
 *
 * \sa mess_lrcf_signfac
 * \sa mess_lrcf_adi
 * \sa mess_sign_scale_t
 */
int mess_lrcf_gsignfac(mess_matrix A, mess_matrix E, mess_matrix B,mess_int_t *maxit, double *tol, mess_matrix Z, mess_sign_scale_t scale){

    MSG_FNAME(__func__);
    int ret=0;
    int convergence = 0;
    mess_int_t n=0;
    mess_int_t onemore=0;
    mess_int_t iter=0;
    double d=0;
    double dU=0;
    mess_int_t j=0;
    mess_int_t r=0;
    mess_int_t *perm,*invperm;
    double Enrm=.0;
    double Err=.0;
    double nrm1=.0;
    double nrm2=.0;
    double r_tol = 1e-8;
    double diagz=.0;
    double m11= .0;
    mess_matrix M1,M2, M3, M4;
    mess_direct sol;
    mess_int_t _maxit = 0;
    double _tol = 0.0;

    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(maxit);
    mess_check_nullpointer(tol);
    _maxit = *maxit;
    _tol = *tol;

    if ( _maxit == 0) {
        // MSG_INFO("Using default _maxit = 50\n");
        _maxit = 50;
    }
    mess_check_positive(_maxit);

    mess_check_nullpointer(A);
    mess_check_nullpointer(E);
    mess_check_same_size(A,E);
    mess_check_real(A);
    mess_check_real(E);
    mess_check_dense(A);
    mess_check_dense(E);
    mess_check_dense(B);


    mess_check_positive(_maxit);
    mess_check_positive(_tol);
    mess_check_nullpointer(Z);

    //check sizes of matrices
    if(!(A->rows==A->cols&&A->cols==E->rows&&E->rows==E->cols&&E->cols==B->cols)){
        MSG_ERROR("Input Matrices have wrong size\n");
        return MESS_ERROR_ARGUMENTS;
    }
    n = A->rows;

    /*-----------------------------------------------------------------------------
     * prepare iteration
     *-----------------------------------------------------------------------------*/

    //initialize M1 M2 M3 M4
    MESS_INIT_MATRICES(&M1,&M2,&M3,&M4);
    ret = mess_matrix_alloc(M1,n, n, n*n, MESS_DENSE, MESS_REAL);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(M2,n, n, n*n, MESS_DENSE, MESS_REAL);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(M3,n, n, n*n, MESS_DENSE, MESS_REAL);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(M4,n, n, n*n, MESS_DENSE, MESS_REAL);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    // MSG_INFO("Matrices initialze\n");
    //allocate storage for Z and B->Z
    ret = mess_matrix_copy(B,Z);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);

    //allocate storage for perm invperm
    mess_try_alloc(perm,mess_int_t * ,sizeof(mess_int_t) * (n));
    mess_try_alloc(invperm,mess_int_t*,sizeof(mess_int_t)*(n));
    for(j =0; j<n;++j){ perm[j]=0; invperm[j]=0;}

    //initialize solver for inverse of A
    ret=mess_direct_init(&sol);                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);

    if(scale==MESS_SIGN_SCALE_DET){
        double m, e;
        ret = mess_direct_create_lapack_lu(E,sol);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_create_lapack_lu);
        ret = mess_direct_determinant(sol,&m,&e);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_determinant);
        dU = pow(fabs(m),1./((double)n))*pow(2,e/((double) n));
        ret = mess_direct_clear(&sol);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_clear);
        ret = mess_direct_init(&sol);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
    }

    //compute fnorm of E
    ret=mess_matrix_normf(E, &Enrm);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);
    if ( _tol == 0.0 ) {
        _tol = sqrt(E->rows *  mess_eps()) * Enrm;
        // MSG_INFO("Using default _tolerance heuristic, _tol = %lg\n", _tol);
    }

    //compute A+E->M1
    ret=mess_matrix_copy(A,M3);                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    ret=mess_matrix_copy(E,M1);                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    ret=mess_matrix_add(1,A,1,M1);                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_add);

    //A+E->M1 compute frobenuis norm of M1
    ret=mess_matrix_normf(M1,&Err);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);

    //_tol=sqrt(n*mess_eps())*Enrm;
    if(Err<=_tol){
        convergence=1;
    } else {
        convergence=0;
    }
    // MSG_INFO("Enrm = %lg Err = %lg\n", Enrm, Err);
    // MSG_INFO("Iteration Matrix Sgn Function starts\n");
    /*-----------------------------------------------------------------------------
     * begin iteration
     *-----------------------------------------------------------------------------*/

    for(iter=0;iter<_maxit;++iter){
        if((! convergence)|| ( convergence && (onemore<2))){

            /*-----------------------------------------------------------------------------
             * compute A^(-1) ->M1
             *-----------------------------------------------------------------------------*/
            mess_direct_clear(&sol);
            mess_direct_init(&sol);
            ret=mess_direct_create_lapack_lu(M3,sol);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_lapack_lu);
            /*-----------------------------------------------------------------------------
             * compute E*A^(-1)*E
             *-----------------------------------------------------------------------------*/
            //A^(-1)*E=M1*E->M2
            ret = mess_direct_solvem(MESS_OP_NONE,sol,E, M2);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solvem);
            //E*A^(-1)*E=E*M2->M1
            ret= mess_matrix_multiply(MESS_OP_NONE, E, MESS_OP_NONE, M2, M1);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);

            if(Err>0.1){
                switch (scale){
                    case MESS_SIGN_SCALE_NONE:
                        d=1;
                        break;
                    case MESS_SIGN_SCALE_FRO:
                        ret=mess_matrix_normf(M3,&nrm1); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);
                        ret=mess_matrix_normf(M1,&nrm2); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);
                        d=sqrt(nrm1/nrm2);
                        break;
                    case MESS_SIGN_SCALE_DET:
                        ret = mess_direct_getU(sol,M4); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_getU);
                        double m, e;
                        ret = mess_direct_determinant(sol,&m,&e);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_determinant);
                        d = pow(fabs(m),1./((double)n))*pow(2,e/((double) n));
                        d/=dU;
                        /*
                        d =0;
                        for(j=0;j<n;++j){
                            d+=log(fabs(M4->values[j+(M4->ld)*j]));
                        }
                        d/=(double)n;
                        d = exp(d);
                        d/=dU;
                        */
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
             * compute (A/d +d*E*A^(-1)*E)/2=(M3/d +d*M1)/2->M3
             *-----------------------------------------------------------------------------*/
            ret = mess_matrix_add(d/2.0,M1,1.0/(2.0*d),M3);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);

            /*-----------------------------------------------------------------------------
             * compute Z =[Z; d*Z*Y*E]/sqrt(2*d)= [Z; d*Z*M2]/sqrt(2*d) -> [Z/sqrt(2*d); M1/sqrt(2/d)]->M2
             *-----------------------------------------------------------------------------*/

            ret = mess_matrix_multiply(MESS_OP_NONE, Z, MESS_OP_NONE, M2, M1);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_scale(d,M1);                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
            // ret = mess_matrix_scale(1/sqrt(2*d), Z); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_scalee);
            ret = mess_matrix_cat(Z,NULL,M1,NULL,MESS_DENSE,M2);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);
            ret = mess_matrix_scale(1.0/sqrt(2.0*d), M2);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);

            //
            /*-----------------------------------------------------------------------------
             * compute  [~,Z,p]=qr(full(Z),0);
             *-----------------------------------------------------------------------------*/
            ret=mess_direct_dgeqp3(M2, NULL, M1, perm);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_direct_dgeqp3);

            /*-----------------------------------------------------------------------------
             * r = sum(diag(Z)>r_tol*abs(Z(1,1))) = sum(diag(M1)>r_tol*(abs(M(1,1))))
             *-----------------------------------------------------------------------------*/
            m11 = fabs(M1->values[0])*r_tol;
            r=0;
            for(j=0; j < MESS_MIN(M1->rows,M1->cols);++j){
                diagz = fabs(M1->values[j+j*M1->ld]);
                if(diagz>m11 ) ++r;
            }
            // MSG_INFO("rank = " MESS_PRINTF_INT "\n", r);

            /*-----------------------------------------------------------------------------
             *  for j=1:rc,  q(j) = find(p==j);  end / inverse permutation of p
             *-----------------------------------------------------------------------------*/
            for(j=0; j<n; ++j){
                invperm[perm[j]]=j;
            }
            /*-----------------------------------------------------------------------------
             *  compute Z = Z(1:r,q)
             *-----------------------------------------------------------------------------*/
            ret= mess_matrix_rowsub(M1,0,r-1,Z); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            ret= mess_matrix_colperm(Z,invperm); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_colperm);

            /* {
                double nrmR;
                mess_matrix_norm2(Z,&nrmR);
                // printf("nrmR= %.15e\n",nrmR);
            } */

            /*------------------------------------------------------------------------
             * compute  Err  = norm(A + E,'fro')
             *-----------------------------------------------------------------------*/
            ret = mess_matrix_copy(M3, M2);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            ret = mess_matrix_add(1,E,1,M2); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
            ret=mess_matrix_normf(M2,&Err); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_normf);

            // MSG_INFO("Step " MESS_PRINTF_INT " conv. crit. = %20.15e\n", iter, Err);

            // mess_matrix_printinfo(Z);
            if(Err<=_tol){ convergence=1; ++onemore;}else{convergence=0;}
        } else {
            break;
        }
    }

    /*------------------------------------------------------------------------
     * finalize computations
     *-----------------------------------------------------------------------*/
    // MSG_INFO("Step " MESS_PRINTF_INT " conv. crit. = %20.15e\n", iter, Err);

    //check
    if(iter==_maxit){
        if (nrm1 > _tol){
            MSG_INFO("LYAP_SGN_FAC: No convergence in " MESS_PRINTF_INT " iterations.\n", _maxit);
        }
    }

    *maxit = iter;
    *tol = Err;

    // ret=mess_matrix_scale(sqrt(2),E); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
    mess_direct_clear(&sol);
    mess_direct_init(&sol);
    ret=mess_direct_create_lapack_lu(E,sol);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_lapack_lu);
    ret=mess_direct_inverse(sol,M1);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_inverse);
    ret= mess_matrix_multiply(MESS_OP_NONE, Z, MESS_OP_NONE, M1, M2);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret=mess_matrix_copy(M2,Z);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret=mess_matrix_scale(1.0/sqrt(2),Z);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
    // mess_matrix_printinfo(Z);

    /*------------------------------------------------------------------------
     * clear storage
     *-----------------------------------------------------------------------*/
    mess_free(perm);
    mess_free(invperm);
    ret=mess_matrix_clear(&M1);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret=mess_matrix_clear(&M2);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret=mess_matrix_clear(&M3);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    ret=mess_matrix_clear(&M4);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    mess_direct_clear(&sol);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_clear);
    return(0);
}

/**
 * @brief Solve a standard Lyapunov Equation \f$ A^{T} X + X A + B^{T} B =0 \f$ via sign function iteration.
 * @param[in] A input \f$(n \times n)\f$-matrix
 * @param[in] B input \f$(p \times n)\f$-matrix
 * @param[in,out] maxit on input: maximum number of iterations \n
 *                      on output: number of performed iterations
 * @param[in,out] tol   on input: convergence tolerance \n
 *                      on output: achieved tolerance
 * @param[in] Z input numerical full rank factor \f$ Z \f$ of the solution \f$ X \approx ZZ^T \f$
 * @param[in] scale input type of scaling for speedup
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcf_signfac function solves a standard Lyapunov Equation
 * \f[ A^{T} X + X A + B^{T} B = 0 \f]
 * via the sign function iteration (\cite BenQ99) . The solution is computed as a low rank factor \f$ Z \f$ of the
 * solution \f$ X \approx ZZ^T \f$.
 *
 * The algorithm stops if the maximum iteration number \f$ maxit \f$ is reached or Frobenius-norm of the current
 * iterate is smaller than \f$ tol \f$
 * \f[ \Vert A_i + I \Vert_F < tol .\f]
 * At exit \f$ maxit \f$ contains the number of performed iterations and \f$ tol \f$ the achieved tolerance.
 *
 * \attention The function works in dense arithmetics and does not support sparse matrices at all.
 *
 * \sa mess_lrcf_gsignfac
 * \sa mess_lrcf_adi
 * \sa mess_sign_scale_t
 */
int mess_lrcf_signfac(mess_matrix A, mess_matrix B, mess_int_t *maxit,double *tol, mess_matrix Z, mess_sign_scale_t scale){
    MSG_FNAME(__func__);
    int ret=0;
    int convergence = 0;
    mess_int_t n=0;
    mess_int_t onemore=0;
    mess_int_t iter=0;
    double d=.0;
    mess_int_t j=0;
    mess_int_t r=0;
    mess_int_t *perm,*invperm;
    mess_int_t diagind=0; //for checking the indices during computing frobenis norm
    double Err=.0;
    double nrm1=.0;
    double nrm2=.0;
    double r_tol = 1e-8;
    double diagz=.0;
    double m11= .0;
    double m,e;
    mess_matrix M1,M2, M3;
    mess_direct sol;
    mess_int_t _maxit;
    double _tol;

    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(maxit);
    mess_check_nullpointer(tol);
    _maxit = *maxit;
    _tol = *tol;

    if ( _maxit <= 0) {
        // MSG_INFO("Using default _maxit = 50\n");
        _maxit = 50;
    }

    mess_check_positive(_tol);
    mess_check_nullpointer(Z);
    mess_check_nullpointer(A);
    mess_check_real(A);
    mess_check_dense(A);
    mess_check_dense(B);

    //check sizes of matrices
    if(!(A->rows==A->cols&&A->cols==B->cols)){
        MSG_ERROR("Input Matrices have wrong size\n");
        return MESS_ERROR_ARGUMENTS;
    }
    n = A->rows;
    if ( _tol == 0.0) {
        _tol = n*sqrt(mess_eps());
        // MSG_INFO("Using default heuristic for _tol, _tol = %lg\n", _tol);
    }

    /*-----------------------------------------------------------------------------
     * prepare iteration
     *-----------------------------------------------------------------------------*/

    //initialize M1 M2 M3
    MESS_INIT_MATRICES(&M1,&M2,&M3);
    ret = mess_matrix_alloc(M1,n, n, n*n, MESS_DENSE, MESS_REAL);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(M2,n, n, n*n, MESS_DENSE, MESS_REAL);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(M3,n, n, n*n, MESS_DENSE, MESS_REAL);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    // MSG_INFO("Matrices initialze\n");
    //allocate storage for Z and B->Z
    ret = mess_matrix_copy(B,Z);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);
    ret = mess_matrix_copy(A,M3);                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);

    //allocate storage for perm invperm
    mess_try_alloc(perm,mess_int_t * ,sizeof(mess_int_t) * (n));
    mess_try_alloc(invperm,mess_int_t*,sizeof(mess_int_t)*(n));
    for(j =0; j<n;++j){ perm[j]=0; invperm[j]=0;}

    //mess direct init qr sol
    mess_direct_init(&sol);

    /*-----------------------------------------------------------------------------------
     * compute frobenuis norm of ||A+E||=||A+I||
     *---------------------------------------------------------------------------------*/
    Err=.0;
    for (j = 0; j < n*n; ++j){
        if(j==diagind){
            Err += (M3->values[j]+1) *(M3->values[j]+1);
            diagind+=(n+1);
        }else{
            Err += (M3->values[j] *M3->values[j]);
        }
    }
    Err=sqrt(Err);

    if(Err<=_tol){
        convergence=1;
    } else {
        convergence=0;
    }
    // MSG_INFO("Err = %lg\nerror", Err);

    // MSG_INFO("Iteration Matrix Sgn Function starts\n");
    /*-----------------------------------------------------------------------------
     * begin iteration
     *-----------------------------------------------------------------------------*/
    for(iter=0;iter<_maxit;++iter){
        if((! convergence)|| ( convergence && (onemore<2))){

            /*-----------------------------------------------------------------------------
             * compute A^(-1) ->M1
             *-----------------------------------------------------------------------------*/
            mess_direct_clear(&sol);
            mess_direct_init(&sol);
            ret = mess_direct_create_lapack_lu(M3,sol);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_direct_create_lapack_lu);
            ret=mess_direct_inverse(sol,M1);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_direct_inverse);

            if(Err>0.1){
                switch (scale){
                    case MESS_SIGN_SCALE_NONE:
                        d=1;
                        break;
                    case MESS_SIGN_SCALE_FRO:
                        ret=mess_matrix_normf(M3,&nrm1); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);
                        ret=mess_matrix_normf(M1,&nrm2); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf);
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
            ret=mess_matrix_add(d/2.0,M1,1.0/(2.0*d),M3); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);

            /*-----------------------------------------------------------------------------
             * compute Z =[Z; d*Z*A^-(1)]/sqrt(2*d)= [Z; d*Z*M1]/sqrt(2*d) -> [Z/sqrt(2*d); M2/sqrt(2*d)]->M1
             *-----------------------------------------------------------------------------*/

            ret = mess_matrix_multiply(MESS_OP_NONE, Z, MESS_OP_NONE, M1, M2);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_scale(d,M2);                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
            ret = mess_matrix_cat(Z,NULL,M2,NULL,MESS_DENSE,M1);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);
            ret = mess_matrix_scale(1.0/sqrt(2.0*d), M1);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
            //
            // mess_matrix_printshort(M1);

            /*-----------------------------------------------------------------------------
             * compute  [~,Z,p]=qr(full(Z),0);
             *-----------------------------------------------------------------------------*/
            //Spaltenanzahl von Z andert sich nicht daher die groesse von perm nicht
            //in mess_direct_dgeqp3 wird jeder Wert von perm neu gesetzt daher brauch ich nicht zu initialisieren
            ret=mess_direct_dgeqp3 (M1, NULL, M2, perm); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_direct_dgeqp3);

            /*-----------------------------------------------------------------------------
             * r = sum(diag(Z)>r_tol*abs(Z(1,1))) = sum(diag(M2)>r_tol*(abs(M2(1,1))))
             *-----------------------------------------------------------------------------*/
            m11 = fabs(M2->values[0])*r_tol;
            r=0;
            for(j=0; j < MESS_MIN(M2->rows,M2->cols);++j){
                diagz = fabs(M2->values[j+j*M2->rows]);
                //  MSG_INFO("%20.15e \t %20.15e\n",diagz,m11);
                if(diagz>m11 ) ++r;
            }
            // MSG_INFO("rank = " MESS_PRINTF_INT "\n", r);

            /*-----------------------------------------------------------------------------
             *  for j=1:rc,  q(j) = find(p==j);  end / inverse permutation of p
             *-----------------------------------------------------------------------------*/
            for(j=0; j<n; ++j){
                invperm[perm[j]]=j;
            }

            /*-----------------------------------------------------------------------------
             *  compute Z = Z(1:r,q)
             *-----------------------------------------------------------------------------*/
            ret= mess_matrix_rowsub(M2,0,r-1,Z); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rowsub);
            ret= mess_matrix_colperm(Z,invperm); FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_colperm);
            /*------------------------------------------------------------------------
             * compute  Err  = norm(A + I,'fro')
             *-----------------------------------------------------------------------*/
            Err=.0;
            diagind=0;
            for (j = 0; j < M3->rows*M3->cols; ++j){
                if(j==diagind){
                    Err += (M3->values[j]+1) *(M3->values[j]+1);
                    diagind+=(n+1);
                }else{
                    Err += (M3->values[j] *M3->values[j]);
                }
            }
            Err=sqrt(Err);

            // MSG_INFO("Step " MESS_PRINTF_INT "  conv. crit. = %20.15e\n", iter, Err);

            // mess_matrix_printinfo(Z);
            if(Err<=_tol){ convergence=1; ++onemore;}else{convergence=0;}
        } else {
            break;
        }
    }

    /*------------------------------------------------------------------------
     * finalize computations
     *-----------------------------------------------------------------------*/
    // MSG_INFO("Step " MESS_PRINTF_INT "  conv. crit. = %20.15e\n", iter, Err);


    //check
    if(iter==_maxit){
        if (nrm1 > _tol){
            MSG_INFO("LYAP_SGN_FAC: No convergence in " MESS_PRINTF_INT  " iterations.\n", _maxit);
        }
    }
    ret=mess_matrix_scale(1.0/sqrt(2),Z);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scal);
    *maxit = iter;
    *tol = Err;
    // mess_matrix_printinfo(Z);

    /*------------------------------------------------------------------------
     * clear storage
     *-----------------------------------------------------------------------*/
    mess_free(perm);
    mess_free(invperm);
    MESS_CLEAR_MATRICES(&M1,&M2,&M3);
    ret=mess_direct_clear(&sol);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_clear);
    return (0);
}



