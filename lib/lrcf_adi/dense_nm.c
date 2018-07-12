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
 * @file lib/lrcf_adi/dense_nm.c
 * @brief Various implementations of the dense Newton method for the algebraic Riccati Equation.
 * @author @koehlerm
 * @author @dykstra
 */


#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "mess/mess.h"
#include "mess/error_macro.h"


// if stepsize of linesearch to small skip linesearch
#define MIN_LINESEARCH_PARA 1e-8

// if coefficients for linesearch are too small or large, skip linesearch
#define LINESEARCH_BREAKDOWN 1e-8


/**
 * @brief Implementation of the dense Newton method to solve the Riccati Equation \f$ A^T X E + E X A +/- X G X + Q = 0 \f$.
 * @param[in] X0           input starting guess symmetric matrix \f$ X0 \f$, @c NULL if zero initual guess
 * @param[in] A            input matrix \f$ A \f$ of the Riccati Equation
 * @param[in] E            input matrix \f$ E \f$ of the Riccati Equation, @c NULL if \f$ E \f$ is identity
 * @param[in] Q            input symmetric matrix \f$ Q \f$ of the Riccati Equation
 * @param[in] G            input symmetric matrix \f$ G \f$ of the Riccati Equation
 * @param[in] plus         input if zero then Riccati Equation with \f$ - \f$ is solved, otherwise \f$ + \f$
 * @param[in] linesearch   input it nonzero linesearch is used
 * @param[in] trans        input parameter to use the dual equation
 * @param[in] maxit        input maximum number of iterations
 * @param[in] nrm_t        input @ref mess_norm_t the desired type of norm
 * @param[in] absres_tol   input tolerance for absolute residual
 * @param[in] relres_tol   input tolerance for relative residual
 * @param[out] absres      output achieved absolute residual, @c NULL if not wanted
 * @param[out] relres      output achieved relative residual, @c NULL if not wanted
 * @param[out] RES         output residual matrix, @c NULL if not wanted
 * @param[out] X           solution of the Riccati Equation
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_dense_nm_gmpare function solves the algebraic Riccati Equation (ARE) with dense real symmetric positive definite square matrices @p G and @p Q and
 * dense real stable square matrix @p A.
 * If @p plus is nonzero the positive ARE with \f$ + \f$ in front of the quadratic term is solved, otherwise the ARE with \f$ - \f$ is solved.
 * If @p E points to @c NULL, then @p E is assumed to be the identity.
 *
 * Depending on @p trans the following AREs are solved:
 *
 * If @p trans is @ref MESS_OP_NONE
 *
 * \f[ A X E^T + E X A^T +/- E X G X E^T + Q = 0 \f]
 *
 * otherwise
 *
 * \f[ A^T X E + E^T X A +/- E^T X G X E + Q = 0 \f]
 *
 * The function uses a Newton method and the resulting Lyapunov Equations are solved using @ref mess_direct_create_generalized_lyapunov.
*
 * If @p X0 does not points to @c NULL, this matrix is used as an initual guess, otherwise a zero inital guess is used.
 *
 * The newton iteration solves for the next Newton step, instead of the directly for the next Newton iteration.
 * This is due to the numerical benfits, see @cite Ben97b and @cite BenB98.
 *
 * The polynomial root finding problem is solved via an generalized eigenvalue problem see @cite JonV06.
 *
 *
 * \attention This function completely relies on real and dense computations. It returns an error if it is called with
 * complex or sparse matrices.
 * \attention The definitenes, symmetry and stability conditions of the matrices are not checked.
 *
 * \sa mess_direct_create_generalized_lyapunov
 * \sa mess_lrcfadi_nm
 *
 * See @cite Ben97b, @cite BenB98 and @cite JonV06 for references.
 */
int mess_dense_nm_gmpare(mess_matrix X0, mess_matrix A, mess_matrix E, mess_matrix Q, mess_matrix G, mess_int_t plus, mess_int_t linesearch, mess_operation_t trans, mess_int_t maxit, mess_norm_t nrm_t,
                            double absres_tol, double relres_tol, double *absres, double *relres, mess_matrix RES, mess_matrix X) {

    MSG_FNAME(__func__);
    mess_int_t it, res_count = 0, ret = 0, i;
    double _res0, _absres, oldres, t = 1.0, alpha, beta, gamma, s, plus_scale=0;
    mess_operation_t dual_trans;
    mess_vector relres_hist, crit_points;
    mess_matrix Atmp, Btmp, tmp, tmp2, NGN, C1, C2;
    mess_direct lyap;


    /*-----------------------------------------------------------------------------
     *  solution of the Riccati Equation
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(G);
    mess_check_nullpointer(Q);
    mess_check_nullpointer(X);
    mess_check_real(A);
    mess_check_real(G);
    mess_check_real(Q);
    mess_check_dense(A);
    mess_check_dense(G);
    mess_check_dense(Q);
    mess_check_square(A);
    mess_check_square(Q);
    mess_check_square(G);
    mess_check_same_size(A,Q);
    mess_check_same_size(A,G);

    if(E){
        mess_check_real(E);
        mess_check_square(E);
        mess_check_dense(E);
        mess_check_same_size(A,E);
    }
    if(X0){
        mess_check_real(X0);
        mess_check_square(X0);
        mess_check_dense(X0);
        mess_check_same_size(A,X0);
    }

    if(!(trans == MESS_OP_NONE || trans == MESS_OP_TRANSPOSE)){
        trans = MESS_OP_TRANSPOSE;
    }

    mess_check_positive(absres_tol);
    mess_check_positive(relres_tol);
    mess_check_positive(maxit);

    MESS_MATRIX_RESET(X);


    /*-----------------------------------------------------------------------------
     *  alloc and prepare matrices for Newton Itertion
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&Atmp, &Btmp, &tmp, &tmp2, &NGN, &C1, &C2);
    MESS_INIT_VECTORS(&relres_hist,&crit_points);

    ret = mess_matrix_alloc(X, A->rows, A->cols, A->nnz, MESS_DENSE, MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_zeros(X);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);

    if(linesearch){
        ret = mess_matrix_alloc(C1, 3, 3, 9, MESS_DENSE, MESS_REAL);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_zeros(C1);                                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);
        ret = mess_matrix_setelement(C1, 0, 1, 1);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setelement);
        ret = mess_matrix_setelement(C1, 1, 2, 1);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setelement);
        ret = mess_matrix_eye(C2, 3, 3, MESS_DENSE);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
        ret = mess_vector_alloc(crit_points, 3, MESS_COMPLEX);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    }

    //set scalar for plus
    plus_scale = plus ? 1.0: -1.0;


    /*-----------------------------------------------------------------------------
     *  vector for residual logging
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_alloc(relres_hist,maxit,MESS_REAL);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_alloc);
    ret = mess_vector_zeros(relres_hist);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_zeros);


    /*-----------------------------------------------------------------------------
     *  dual transpose type
     *-----------------------------------------------------------------------------*/
    dual_trans = trans==MESS_OP_NONE? MESS_OP_TRANSPOSE: MESS_OP_NONE;


    /*-----------------------------------------------------------------------------
     *  compute norm of right hand side Q
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_norm(Q, nrm_t, &_res0);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_normf);


    /*-----------------------------------------------------------------------------
     *  handle initial guess
     *-----------------------------------------------------------------------------*/
    if(X0){
        // compute E^T*X0*G*X0*E -> tmp2
        ret = mess_matrix_multiply(MESS_OP_NONE, G, MESS_OP_NONE, X0, tmp);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, X0, MESS_OP_NONE, tmp, tmp2);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

        if(E){
            ret = mess_matrix_multiply(trans, E, MESS_OP_NONE, tmp2, tmp);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_NONE, tmp, dual_trans, E, tmp2);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        }

        // compute +/- E^T*X0*G*X0*E + Q -> Btmp
        ret = mess_matrix_copy(Q,Btmp);                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        ret = mess_matrix_add(plus_scale, tmp2, 1, Btmp);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);

        // compute E^T*X0*A -> tmp
        ret = mess_matrix_multiply(MESS_OP_NONE, X0, dual_trans, A, tmp);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

        // compute A^T*X0*E + E^T*X0*A +/- E^T X0 G X0 E + Q
        if(E){
            ret = mess_matrix_multiply(trans, E, MESS_OP_NONE, tmp, tmp2);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_add(1,tmp2,1,Btmp);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
            ret = mess_matrix_transpose(tmp2,tmp);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_transpose);
            ret = mess_matrix_add(1,tmp,1,Btmp);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        } else {
            ret = mess_matrix_add(1,tmp,1,Btmp);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
            ret = mess_matrix_transpose(tmp,tmp2);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_transpose);
            ret = mess_matrix_add(1,tmp2,1,Btmp);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        }

        /*-----------------------------------------------------------------------------
         *  make Btmp symmetric
         *-----------------------------------------------------------------------------*/
        //ret = mess_matrix_transpose(Btmp,tmp2);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_transpose);
        //ret = mess_matrix_add(0.5,tmp2,0.5,Btmp);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_proj_sym(Btmp);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_proj_sym);
    }


    /*-----------------------------------------------------------------------------
     *  main iteration
     *-----------------------------------------------------------------------------*/
    for ( it =0; it < maxit; it++){

        //
        // Here we solve for the Newton step and add it to X, due to numerical advantages.
        // See the @cite Ben97b and @cite BenB98 papers in doc for references.
        //
        if ( it == 0 && X0 == NULL) {
            ret = mess_matrix_copy(Q,Btmp);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            ret = mess_direct_init(&lyap);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
            ret = mess_direct_create_generalized_lyapunov(A, E, lyap);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_generalized_lyapunov);
            ret = mess_direct_solvem(trans, lyap, Btmp, tmp);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
            mess_direct_clear(&lyap);

            //
            // X0 is zero because of MESS_MATRIX_RESET therefore, absres = res0
            //
            _absres = _res0;
        } else {
            // compute A-G*X*E
            if(trans == MESS_OP_NONE){
                ret = mess_matrix_transpose(A,Atmp);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_transpose);
            } else {
                ret = mess_matrix_copy(A,Atmp);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            }
            ret = mess_matrix_multiply(MESS_OP_NONE, G, MESS_OP_NONE, X, tmp);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            if(E){
                ret = mess_matrix_multiply(MESS_OP_NONE, tmp, dual_trans, E, tmp2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
                ret = mess_matrix_add(plus_scale, tmp2, 1, Atmp);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
            } else {
                ret = mess_matrix_add(plus_scale, tmp, 1, Atmp);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
            }

            ret = mess_direct_init(&lyap);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
            ret = mess_direct_create_generalized_lyapunov(Atmp, E,lyap);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_generalized_lyapunov);
            ret = mess_direct_solvem(MESS_OP_TRANSPOSE, lyap, Btmp, tmp);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
            mess_direct_clear(&lyap);
        }


        /*-----------------------------------------------------------------------------
         *  handle linesearch
         *-----------------------------------------------------------------------------*/
        /* Set Stepsize to 1.0 before Linesearch process*/
        t = 1.0;
        if(linesearch){

            // NGN = E^T*NewtonStep*G*Newtonstep*E, NewtonStep<-tmp
            ret = mess_matrix_multiply(MESS_OP_NONE, G, MESS_OP_NONE, tmp, tmp2);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_NONE, tmp, MESS_OP_NONE, tmp2, NGN);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            if(E){
                ret = mess_matrix_multiply(MESS_OP_NONE, NGN, trans, E, tmp2);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
                ret = mess_matrix_multiply(dual_trans, E, MESS_OP_NONE, tmp2, NGN);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            }

            // alpha = trace(Btmp*Btmp)=FrobeniusNorm(Btmp)^2
            ret = mess_matrix_normf2(Btmp, &alpha);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_normf2);

            // beta = trace(Btmp, NGN)
            ret = mess_matrix_fro_inn(MESS_OP_NONE, Btmp, NGN, (void*)&beta);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_fro_inn);

            // modifiy beta in case of plus
            beta = plus? -beta:beta;

            // gamma = trace(NGN^2)=FrobeniusNorm(NGN)
            ret = mess_matrix_normf2(NGN, &gamma);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_normf2);


            /*-----------------------------------------------------------------------------
             *  If coefficients of polynomial are extremely small
             *  skip linesearch.
             *
             *  Matrix pencil approach to compute roots of polynomial, see @cite JonV06
             *  in doc.

             *-----------------------------------------------------------------------------*/
            if(fabs(alpha)>LINESEARCH_BREAKDOWN && fabs(beta)>LINESEARCH_BREAKDOWN && fabs(gamma)>LINESEARCH_BREAKDOWN){

                s = 4*gamma + 6*beta + 2*alpha - 4*beta - 2*alpha;

                ret = mess_matrix_setelement(C1, 2, 0, 2*alpha/s);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setelement);
                ret = mess_matrix_setelement(C1, 2, 1, (-2*alpha+4*beta)/s);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setelement);
                ret = mess_matrix_setelement(C1, 2, 2, plus_scale*6*(beta/s));              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setelement);

                ret = mess_matrix_setelement(C2, 2, 2, 4*(gamma/s));                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setelement);
                ret = mess_eigen_eigg(C1, C2, crit_points, NULL, NULL, NULL);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eigg);


                /*-----------------------------------------------------------------------------
                 *  check critical points
                 *-----------------------------------------------------------------------------*/
                for(i=0; i<3; i++){

                    // complex  roots can be skipped
                    s = cimag(crit_points->values_cpx[i]);
                    if(fabs(s) > 1e-12){
                        continue;
                    }

                    // real roots must be in interval [0, 2]
                    s = creal(crit_points->values_cpx[i]);
                    if(s < 0 || s > 2){
                        continue;
                    }

                    // find the minimizer
                    if(alpha*pow((1-t),2) - 2*beta*(1-t)*pow(t,2) + gamma*pow(t,4) > alpha*pow((1-s),2) - 2*beta*(1-s)*pow(s,2) + gamma*pow(s,4)) t = s;
                }
            }

            // if linesearch parameter to small set it to 1.0
            if(t < MIN_LINESEARCH_PARA) t = 1.0;

            MSG_INFO("it="MESS_PRINTF_INT", Line Search Parameter = %e\n", it, t);
        }

        /*-----------------------------------------------------------------------------
         *  Update Iterate with Newton Step
         *-----------------------------------------------------------------------------*/
        ret = mess_matrix_add(t, tmp, 1, X);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);


        /*-----------------------------------------------------------------------------
         *  compute residual
         *-----------------------------------------------------------------------------*/
        // compute E^T*X0*G*X0*E -> tmp2
        ret = mess_matrix_multiply(MESS_OP_NONE, G, MESS_OP_NONE, X, tmp);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, X, MESS_OP_NONE, tmp, tmp2);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        if(E){
            ret = mess_matrix_multiply(trans, E, MESS_OP_NONE, tmp2, tmp);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_NONE, tmp, dual_trans, E, tmp2);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        }

        // compute +/- E^T*X0*G*X0*E + Q -> Btmp
        ret = mess_matrix_copy(Q, Btmp);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        ret = mess_matrix_add(plus?1:-1, tmp2, 1, Btmp);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);

        // compute E^T*X0*A -> tmp
        ret = mess_matrix_multiply(MESS_OP_NONE, X, dual_trans, A, tmp);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

        // compute A^T*X0*E + E^T*X0*A +/- E^T X0 G X0 E + Q
        if(E){
            ret = mess_matrix_multiply(trans, E, MESS_OP_NONE, tmp, tmp2);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_add(1,tmp2,1,Btmp);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
            ret = mess_matrix_transpose(tmp2,tmp);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_transpose);
            ret = mess_matrix_add(1,tmp,1,Btmp);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        } else {
            ret = mess_matrix_add(1,tmp,1,Btmp);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
            ret = mess_matrix_transpose(tmp,tmp2);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_transpose);
            ret = mess_matrix_add(1,tmp2,1,Btmp);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        }

        /*-----------------------------------------------------------------------------
         *  make sure Btmp is symmetric
         *-----------------------------------------------------------------------------*/
        //ret = mess_matrix_transpose(Btmp,tmp2);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_transpose);
        //ret = mess_matrix_add(0.5,tmp2,0.5,Btmp);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_proj_sym(Btmp);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_proj_sym);

        oldres = _res0;
        ret = mess_matrix_norm(Btmp, nrm_t, &_absres);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm);


        /*-----------------------------------------------------------------------------
         *  print INFO about residual
         *-----------------------------------------------------------------------------*/
        MSG_INFO("it="MESS_PRINTF_INT", abs./rel. %s residual = %e / %e\n", it, mess_norm_t_str(nrm_t), _absres, _absres/_res0);


        /*-----------------------------------------------------------------------------
         * log residual
         *-----------------------------------------------------------------------------*/
        relres_hist->values[it]=_absres/_res0;
        if (_absres/_res0 < relres_tol && _absres < absres_tol){
            res_count++;
        }

        /*-----------------------------------------------------------------------------
         *  if res_count is large enough or change in residual is to small exit
         *  newton iteration
         *-----------------------------------------------------------------------------*/
        if ( res_count >= 3 ) break;
        if (it>0 && fabs(oldres-_absres) < mess_eps()){
            MSG_INFO("No change in absolute residual, therfore stop iteration\n");
            break;
        }
    }//* End of Main Loop *//


    /*-----------------------------------------------------------------------------
     *  set output arguments if wanted
     *-----------------------------------------------------------------------------*/
    if(absres) *absres = _absres;
    if(relres) *relres = _absres / _res0;
    if(RES){
        ret = mess_matrix_copy(Btmp,RES);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    }

    /*-----------------------------------------------------------------------------
     *  check if Newton iteration has convergend
     *-----------------------------------------------------------------------------*/
    if ( it == maxit) {
        MSG_ERROR("The Newton iteration did not converge within "MESS_PRINTF_INT" iterations.\n", maxit);
        for(i=0; i<maxit; ++i){
            MSG_INFO("it= " MESS_PRINTF_INT " resF/resF0 = %e\n",i,relres_hist->values[i]);
        }
        ret = MESS_ERROR_CONVERGE;
    }

    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&Atmp, &Btmp, &tmp, &tmp2, &NGN, &C1, &C2);
    MESS_INIT_VECTORS(&relres_hist,&crit_points);

    return ret;
}




#define MAXIT 50
#define RELRES_TOL sqrt(mess_eps())
#define ABSRES_TOL 1e-10


/**
 * @brief Implementation of the dense Newton-method to solve the generalized Riccati Equation \f$ A^T X E+E^T X A - E^T XG X E + Q = 0 \f$.
 * @param[in] A      input matrix \f$ A \f$ of the Riccati Equation
 * @param[in] E      input matrix \f$ E \f$ of the Riccati Equation, @c NULL if @p E is identity
 * @param[in] Q      input matrix \f$ Q \f$ of the Riccati Equation
 * @param[in] G      input matrix \f$ G \f$ of the Riccati Equation
 * @param[out] X    solution of the Riccati Equation
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_dense_nm_gmare function solves the generalized algebraic Riccati Equation (GARE)
 * \f[ A^T X E + E^T X A - E^T X G X E + Q = 0 \f]
 * using the Newton method. \n
 * The resulting internal generalized Lyapunov Equation
 * \f[ (A-GX_iE)^T X_{i+1}E + E^TX_{i+1} (A-GX_iE) - \left(Q + E^TX_iGX_iE \right) =0  \f]
 * is solved using a dense generalized  Lyapunov solver.
 *
 * \attention This function completely relies on dense computations. It returns an error if it is called with sparse
 * matrices.
 *
 * \sa mess_direct_create_generalized_lyapunov
 * \sa mess_lrcfadi_nm
 * \sa mess_dense_nm_gmpare
 *
 */
int mess_dense_nm_gmare(mess_matrix A, mess_matrix E,  mess_matrix Q, mess_matrix G, mess_matrix X) {

    MSG_FNAME(__func__);
    int ret = 0;
    ret = mess_dense_nm_gmpare(NULL, A, E, Q, G, 0, 0, MESS_OP_TRANSPOSE, MAXIT, MESS_2_NORM, ABSRES_TOL, RELRES_TOL, NULL, NULL, NULL, X);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dense_nm_gmpare);
    return 0;

}


/**
 * @brief Implementation of the dense Newton-method to solve the generalized Riccati Equation \f$ A^T X E+E^T X A + E^T XG X E + Q = 0 \f$.
 * @param[in] A      input matrix \f$ A \f$ of the Riccati Equation
 * @param[in] E      input matrix \f$ E \f$ of the Riccati Equation, @c NULL if @p E is identity
 * @param[in] Q      input matrix \f$ Q \f$ of the Riccati Equation
 * @param[in] G      input matrix \f$ G \f$ of the Riccati Equation
 * @param[out] X    solution of the Riccati Equation
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_dense_nm_gpare function solves the generalized algebraic Riccati Equation (GARE)
 * \f[ A^T X E + E^T X A + E^T X G X E + Q = 0 \f]
 * using the Newton method. \n
 * The resulting internal generalized Lyapunov Equation
 * \f[ (A-GX_iE)^T X_{i+1}E + E^TX_{i+1} (A-GX_iE) + \left(Q + E^TX_iGX_iE \right) =0  \f]
 * is solved using a dense generalized  Lyapunov solver.
 *
 * \attention This function completely relies on dense computations. It returns an error if it is called with sparse
 * matrices.
 *
 * \sa mess_direct_create_generalized_lyapunov
 * \sa mess_lrcfadi_nm
 * \sa mess_dense_nm_gmpare
 *
 */
int mess_dense_nm_gpare(mess_matrix A, mess_matrix E,  mess_matrix Q, mess_matrix G, mess_matrix X) {

    MSG_FNAME(__func__);
    int ret = 0;
    ret = mess_dense_nm_gmpare(NULL, A, E, Q, G, 1, 0, MESS_OP_TRANSPOSE, MAXIT, MESS_2_NORM, ABSRES_TOL, RELRES_TOL, NULL, NULL, NULL, X);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dense_nm_gmpare);
    return 0;

}


/**
 * @brief Compute the residual of the algebraic (generalized) Riccati (or Lyapunov) equation.
 * @param[in] A         input matrix
 * @param[in] E         input matrix \f$ E \f$ of the Riccati Equation, @c NULL if \f$ E \f$ is identity
 * @param[in] Q         input symmetric matrix
 * @param[in] G         input symmetric matrix (@p NULL if \f$ G \f$ is zero -> Lyapunov Equation)
 * @param[in] X         input symmetric solution of the Riccati Equation
 * @param[in] plus      input if zero then Riccati Equation with \f$ - \f$ is considered, otherwise \f$ + \f$
 * @param[in] trans     input parameter to use the dual equation
 * @param[in] nrm_t     input parameter to control desired norm type
 * @param[out] res      output asbolute residual
 * @param[out] rel      output relative residual ( absolute residual divided by norm of \f$Q\f$)
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_dense_res_gmpare  function computes residual of dense algebraic Riccati Equation (ARE).
 *
 * If @p plus is nonzero the positive ARE with \f$ + \f$ in front of the quadratic term is considered, otherwise the ARE with \f$ - \f$ is considered.
 *
 * If @p E points to @c NULL, then @p E is assumed to be the identity.
 *
 * Depending on @p trans the following reisduals are computed
 *
 * If @p trans is @ref MESS_OP_NONE
 *
 * \f[ \Vert A X E^T + E X A^T +/- E X G X E^T + Q \Vert \f]
 *
 * otherwise
 *
 * \f[ \Vert A^T X E + E^T X A +/- E^T X G X E + Q \Vert \f]

 *
 * If \f$ G \f$ is @c NULL the residual of a Lyapunov Equation is computed.
 *
 * @p nrm_t is used to define the desired norm.
 *
 * \f$ rel \f$ gives the relative residual \f[ rel=res/\Vert Q \Vert \f], if not wanted call the function with @p rel equals @c NULL.
 *
 *
 */

int mess_dense_res_gmpare(mess_matrix A, mess_matrix E,  mess_matrix Q, mess_matrix G, mess_matrix X, mess_int_t plus, mess_operation_t trans, mess_norm_t nrm_t, double *res, double *rel){
    MSG_FNAME(__func__);
    int ret =0;
    double normQ;
    mess_matrix TMP, TMP2, R;

    /*-----------------------------------------------------------------------------
     *  check nullpointer
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(Q);
    mess_check_nullpointer(X);
    mess_check_nullpointer(res);
    mess_check_real(A);
    mess_check_real(Q);
    mess_check_real(X);
    mess_check_same_size(A,Q);
    mess_check_same_size(A,X);

    if(E){
        mess_check_real(E);
        mess_check_same_size(A,E);
    }

    if(G){
        mess_check_real(G);
        mess_check_same_size(A,G);
    }

    if(!(trans == MESS_OP_NONE || trans == MESS_OP_TRANSPOSE)){
        MSG_ERROR("Operation Type not supported\n");
        return MESS_ERROR_ARGUMENTS;
    }

    mess_int_t dual_trans = trans==MESS_OP_NONE ? MESS_OP_TRANSPOSE : MESS_OP_NONE;

    /*-----------------------------------------------------------------------------
     *  init
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&TMP,&TMP2,&R);

    /*-----------------------------------------------------------------------------
     *  compute residual matrix
     *-----------------------------------------------------------------------------*/

    //compute A^TXE
    ret = mess_matrix_multiply(trans, A, MESS_OP_NONE, X, R);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    if(E){
        ret = mess_matrix_multiply(MESS_OP_NONE, R, dual_trans, E, TMP);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_ctranspose(TMP,R);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
    }else{
        ret = mess_matrix_ctranspose(R,TMP);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
    }

    //compute A^TXE + E^TXA
    ret = mess_matrix_add(1,TMP,1,R);                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);

    // compute A^TXE + E^TXA + Q
    ret = mess_matrix_add(1,Q,1,R);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);

    if(G){
        //compute E^TXGXE
        ret = mess_matrix_multiply(MESS_OP_NONE, G, MESS_OP_NONE, X, TMP);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_NONE, X, MESS_OP_NONE, TMP, TMP2);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        if(E){
            ret = mess_matrix_multiply(trans, E, MESS_OP_NONE, TMP2, TMP);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_NONE, TMP, dual_trans, E, TMP2);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        }
        //compute A^TX + E^TXA + Q  +/- E^TXGXE
        if(plus){
            ret = mess_matrix_add(1,TMP2,1,R);                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        } else {
            ret = mess_matrix_add(-1,TMP2,1,R);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        }
    }

    /*-----------------------------------------------------------------------------
     *  compute the 2 / Frobenius norm
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_norm(R, nrm_t, res);                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm);

    if(rel){
        ret = mess_matrix_norm(Q, nrm_t, &normQ);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_norm);
        *rel = *res/normQ;
    }

    /*-----------------------------------------------------------------------------
     * clear memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&TMP, &TMP2, &R);

    return 0;
}       /* -----  end of function mess_dense_res_gmpare  ----- */

