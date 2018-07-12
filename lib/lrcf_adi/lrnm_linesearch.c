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
 * @file lib/lrcf_adi/lrnm_linesearch.c
 * @brief Implementation of the linesearch function for the Newton method.
 * @author @mbehr
 *
 * @cite  BenHSetal16
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "mess/mess.h"
#include "mess/error_macro.h"

/**
 *
 * @brief Compute the coefficients of the polynomial for linesearch.
 * @param[in] W     input last residual factor from @ref mess_lradi
 * @param[in] Wnew  input actual residual factor from @ref mess_lradi
 * @param[in] deltaK    input last change in feedback.
 * @param[in] deltaKnew input actual change in feedback.
 * @param[out] alpha    output coefficient \f$ \alpha   \f$ of the polynomial.
 * @param[out] beta     output coefficient \f$ \beta    \f$ of the polynomial.
 * @param[out] delta    output coefficient \f$ \delta   \f$ of the polynomial.
 * @param[out] gamma    output coefficient \f$ \gamma   \f$ of the polynomial.
 * @param[out] epsilon  output coefficient \f$ \epsilon \f$ of the polynomial.
 * @param[out] zeta     output coefficient \f$ \zeta    \f$ of the polynomial.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_nm_linesearch_poly function computes the coefficients
 * of the polynomial needed by linesearch. More Information can be found in
 * @cite  BenHSetal16.
 *
 */
int mess_lrcfadi_nm_linesearch_poly(mess_matrix W, mess_matrix Wnew, mess_matrix deltaK, mess_matrix deltaKnew,
        double *alpha, double *beta, double *delta, double *gamma, double *epsilon, double *zeta){
    MSG_FNAME(__func__);
    mess_int_t ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(alpha);
    mess_check_nullpointer(beta);
    mess_check_nullpointer(delta);
    mess_check_nullpointer(gamma);
    mess_check_nullpointer(epsilon);
    mess_check_nullpointer(zeta);

    mess_check_nullpointer(W);          mess_check_real(W);         mess_check_dense(W);
    mess_check_nullpointer(Wnew);       mess_check_real(Wnew);      mess_check_dense(Wnew);
    mess_check_nullpointer(deltaKnew);  mess_check_real(deltaKnew); mess_check_dense(deltaKnew);
    mess_check_nullpointer(deltaK);     mess_check_real(deltaK);    mess_check_dense(deltaK);
    mess_check_same_rows(W,Wnew);

    if(W->rows < W->cols){
        MSG_ERROR("W has more cols than rows.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if(Wnew->rows < Wnew->cols){
        MSG_ERROR("Wnew has more cols than rows.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  transpose matrices if necessary
     *-----------------------------------------------------------------------------*/
    mess_matrix deltaK_, deltaKnew_;
    MESS_INIT_MATRICES(&deltaK_,&deltaKnew_);

    if(deltaK->rows > deltaK->cols){
        ret = mess_matrix_copy(deltaK,deltaK_);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    }else{
        ret = mess_matrix_ctranspose(deltaK,deltaK_);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    }

    if(deltaKnew->rows > deltaKnew->cols){
        ret = mess_matrix_copy(deltaKnew,deltaKnew_);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
    }else{
        ret = mess_matrix_ctranspose(deltaKnew,deltaKnew_);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
    }

    //now all matrices should have the right format
    mess_check_same_rows(W,deltaKnew_);
    mess_check_same_rows(deltaKnew_,deltaK_);


    /*-----------------------------------------------------------------------------
     *  compute coefficients of polynomial
     *-----------------------------------------------------------------------------*/
    //compute alpha = ||R(X_K)||_F^2 = || W*W'-deltaK*deltaK'||_F^2
    *alpha=0;
    if(deltaK){
        ret = mess_matrix_indefinite_dynormf2(W,deltaK_,alpha);                                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_indefinite_dynormf2);
    }else{
        ret = mess_matrix_dynormf2(W,alpha);                                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_dynormf2);
    }

    //compute beta = ||Wnew*Wnew'||_F^T
    *beta=0;
    ret = mess_matrix_dynormf2(Wnew,beta);                                                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_dynormf2);

    //compute delta = ||deltaKnew*deltaKnew'||_F^T
    *delta=0;
    ret = mess_matrix_dynormf2(deltaKnew_,delta);                                                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_dynormf2);

    //compute gamma = gamma1-gamma2=||W_old'*W_new||_F^2-||deltaK*deltaK'||_F^2
    *gamma=0;
    double gamma1=0,gamma2=0;
    //compute gamma1
    ret = mess_matrix_mulnormf2(MESS_OP_TRANSPOSE,W,MESS_OP_NONE,Wnew,&gamma1);                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mulnormf2);
    if(deltaK){
        ret = mess_matrix_mulnormf2(MESS_OP_TRANSPOSE,deltaK_,MESS_OP_NONE,Wnew,&gamma2);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mulnormf2);
    }
    *gamma = gamma1-gamma2;

    //compute epsilon
    *epsilon=0;
    double epsilon1=0, epsilon2=0;
    //compute epsilon1
    ret = mess_matrix_mulnormf2(MESS_OP_TRANSPOSE,W,MESS_OP_NONE, deltaKnew_, &epsilon1);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mulnormf2);
    //compute epsilon2
    if(deltaK){
        ret = mess_matrix_mulnormf2(MESS_OP_TRANSPOSE,deltaK_,MESS_OP_NONE,deltaKnew_, &epsilon2);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mulnormf2);
    }
    *epsilon = epsilon1 - epsilon2;

    //compute zeta
    *zeta=0;
    ret = mess_matrix_mulnormf2(MESS_OP_TRANSPOSE,Wnew,MESS_OP_NONE,deltaKnew_,zeta);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mulnormf2);


    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&deltaK_,&deltaKnew_);

    return 0;

}

static int __filter_linesearch_real(const double * val){return 0<(*val) && (*val)<1;};
static int __filter_linesearch_cpx(const mess_double_cpx_t * val){return ((cimag(*val)==0) && (0<creal(*val) && creal(*val)<1)) ;};
static double linesearch_poly_eval(double t, double *alpha, double *beta, double *delta, double *gamma, double *epsilon, double *zeta){
    return pow(1-t,2)*(*alpha)+pow(t,2)*(*beta)+pow(t,4)*(*delta)+2*t*(1-t)*(*gamma)-2*pow(t,2)*(1-t)*(*epsilon)-2*pow(t,3)*(*zeta);
}

/**
 *
 * @brief Compute the optimal step size for linesearch.
 * @param[in] alpha     input coefficient \f$ \alpha    \f$ of the polynomial.
 * @param[in] beta  input coefficient \f$ \beta     \f$ of the polynomial.
 * @param[in] delta     input coefficient \f$ \delta    \f$ of the polynomial.
 * @param[in] gamma input coefficient \f$ \gamma    \f$ of the polynomial.
 * @param[in] epsilon   input coefficient \f$ \epsilon \f$ of the polynomial.
 * @param[in] zeta  input coefficient \f$ \zeta     \f$ of the polynomial.
 * @param[out] lambda   output optimal step length \f$ \lambda  \f$ of the polynomial.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_nm_linesearch_lambda function computes the optimal
 * feasible step lenght  \f$ \lambda \f$ based on the coefficients of
 * the polynomial.  More Information can be found in
 * @cite  BenHSetal16.
 *
 */
int mess_lrcfadi_nm_linesearch_lambda(double *alpha, double* beta, double *delta, double *gamma, double *epsilon, double *zeta, double *lambda){
    MSG_FNAME(__func__);
    mess_int_t ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(alpha);
    mess_check_nullpointer(beta);
    mess_check_nullpointer(delta);
    mess_check_nullpointer(gamma);
    mess_check_nullpointer(epsilon);
    mess_check_nullpointer(zeta);
    mess_check_nullpointer(lambda);

    /*-----------------------------------------------------------------------------
     *  compute lambda via eigenproblem DOI. 10.1137/S0895479899365720
     *-----------------------------------------------------------------------------*/
    // setup entries
    double a0 = 2*((*gamma)-(*alpha));
    double a1=2*((*alpha)+(*beta)-2*((*gamma)+(*epsilon)));
    double a2=6*((*epsilon)-(*zeta));
    double a3=4*(*delta);
    double anrm = sqrt(pow(a0,2)+pow(a1,2)+pow(a2,2)+pow(a3,2));

    //set up matrices
    mess_matrix A,B;
    MESS_INIT_MATRICES(&A,&B);
    ret = mess_matrix_alloc(A,3,3,3*3,MESS_DENSE,MESS_REAL);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
    A->values[2+0*(A->ld)]=-a0/anrm;
    A->values[0+1*(A->ld)]= 1;
    A->values[2+1*(A->ld)]=-a1/anrm;
    A->values[1+2*(A->ld)]= 1;
    A->values[2+2*(A->ld)]=-a2/anrm;
    ret = mess_matrix_eye(B,3,3,MESS_DENSE);                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_eye);
    B->values[2+2*(B->ld)]=a3/anrm;

    //solve eigenvalue problem
    mess_vector ev;
    ret = mess_vector_init(&ev);                                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
    ret = mess_vector_alloc(ev,3,MESS_COMPLEX);                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_alloc);
    ret = mess_eigen_eigg(A,B,ev,NULL,NULL,NULL);                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_eigen_eigg);

    /*-----------------------------------------------------------------------------
     *  find optimal feasible value
     *-----------------------------------------------------------------------------*/
    //filter complex values and eigenvalues outside (0,1)
    ret = mess_vector_filter(ev,__filter_linesearch_real,__filter_linesearch_cpx);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_filter);
    ret = mess_vector_toreal_nowarn(ev);                                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal_nowarn);

    if(ev->dim==0){
        //no feasible lambda found
        MSG_INFO("No feasible lambda found. Set lambda to 1.0\n");
        *lambda=1;
    }else if(ev->dim==1){
        //exactly one feasible lambda found
        *lambda = ev->values[0];
    }else{
        //more than one feasible lambda found, find minimizer, by evaluating the polynomial
        mess_int_t i =1;
        mess_int_t minarg = 0;
        double minval = linesearch_poly_eval(ev->values[0],alpha,beta,delta,gamma,epsilon,zeta);
        double tempval;
        for(i=1;i<ev->dim;++i){
            tempval =  linesearch_poly_eval(ev->values[i],alpha,beta,delta,gamma,epsilon,zeta);
            if(tempval<minval){
                minval = tempval;
                minarg = i;
            }
        }
        *lambda = ev->values[minarg];
    }

    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&B);
    MESS_CLEAR_VECTORS(&ev);

    return 0;
}

/**
 *
 * @brief Compute the optimal step size for linesearch.
 * @param[in] W     input last residual factor from @ref mess_lradi
 * @param[in] Wnew  input actual residual factor from @ref mess_lradi
 * @param[in] deltaK    input last change in feedback.
 * @param[in] deltaKnew input actual change in feedback.
 * @param[out] lambda   output optimal step length \f$ \lambda  \f$ of the polynomial.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_nm_linesearch function computes the optimal
 * feasible step lenght  \f$ \lambda \f$.  More Information can be found in
 * @cite  BenHSetal16.
 *
 * @sa mess_lrcfadi_nm_linesearch_lambda
 * @sa mess_lrcfadi_nm_linesearch_poly
 *
 */
int mess_lrcfadi_nm_linesearch(mess_matrix W, mess_matrix Wnew, mess_matrix deltaK, mess_matrix deltaKnew,double *lambda){
    MSG_FNAME(__func__);
    mess_int_t ret = 0;
    double alpha, beta, delta, gamma, epsilon, zeta;

    /*-----------------------------------------------------------------------------
     *  compute coefficient of the polynomial
     *-----------------------------------------------------------------------------*/
    ret = mess_lrcfadi_nm_linesearch_poly(W, Wnew, deltaK, deltaKnew, &alpha, &beta, &delta, &gamma, &epsilon, &zeta);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_lrcfadi_nm_linesearch_poly);

    /*-----------------------------------------------------------------------------
     *  minimize the polynomial (compute lambda as minimizer in (0,1))
     *-----------------------------------------------------------------------------*/
    ret = mess_lrcfadi_nm_linesearch_lambda(&alpha, &beta, &delta, &gamma, &epsilon, &zeta, lambda);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_lrcfadi_nm_linesearch_lambda);

    MSG_INFO("Set Step Size lambda = %e\n",*lambda);

    return ret;
}

