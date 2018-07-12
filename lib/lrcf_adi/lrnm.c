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
 * @file lib/lrcf_adi/lrnm.c
 * @brief Implementation of the Newton method.
 * @author @koehlerm
 */


#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "mess/mess.h"
#include "mess/error_macro.h"


#define IS_REAL(X) (cimag((X))==0.0)

#ifdef MESS_MATLAB
#define SET_COLOR       MSG_PRINT("");
#define UNSET_COLOR     MSG_PRINT("");
#else
#define SET_COLOR       MSG_PRINT("\033[35m");
#define UNSET_COLOR     MSG_PRINT("\n\033[0m");
#endif

/**
 * @internal
 * @brief Perform  \f$ A = [\alpha A, \beta B] \f$ or \f$ A = [\alpha A; \beta B] \f$.
 * @param[in] op        input operation type
 * @param[in] alpha     input @c double scalar \f$ \alpha \f$
 * @param[in] A         input @ref mess_matrix A
 * @param[in] beta      input @c double scalar \f$ \beta \f$
 * @param[in,out] B     input @ref mess_matrix B
 * @return zero on succes or a non-zero error value otherwise
 *
 * The @ref _scalecat function performs  \f$ A = [\alpha A, \beta B] \f$
 * or \f$ A = [\alpha A; \beta B] \f$.\n
 * If @p op is equal to @ref MESS_OP_NONE the first operation will be performed, otherwise the second one.
 *
 * The output @p B will be @ref MESS_DENSE.
 *
 * @attention Internal use only.
 */

static int _scalecat(mess_operation_t op ,double alpha, mess_matrix A, double beta, mess_matrix B){
    MSG_FNAME(__func__);
    mess_int_t ret =0;
    mess_matrix tmp1, tmp2, tmp3;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_operation_type(op);

    /*-----------------------------------------------------------------------------
     *  perform operation
     *-----------------------------------------------------------------------------*/
    if(B){

        MESS_INIT_MATRICES(&tmp1,&tmp2,&tmp3);
        ret = mess_matrix_copy(A,tmp1);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_copy(B,tmp2);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);

        ret = mess_matrix_scale(alpha,tmp1);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
        ret = mess_matrix_scale(beta,tmp2);                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);

        if(op==MESS_OP_NONE){
            ret = mess_matrix_cat(tmp1,tmp2,NULL,NULL,MESS_DENSE, tmp3);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);
        }else{
            ret = mess_matrix_cat(tmp1,NULL,tmp2,NULL,MESS_DENSE, tmp3);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);
        }

        ret = mess_matrix_copy(tmp3,A);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);

        MESS_CLEAR_MATRICES(&tmp1,&tmp2,&tmp3);
    }else{
        //perform only A = [alpha*A]
        ret = mess_matrix_scale(alpha,A);                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
    }

    /*-----------------------------------------------------------------------------
     *  clear additional data
     *-----------------------------------------------------------------------------*/
    return ret;

}


/**
 * @brief Wrapper/Alternative name for \ref mess_lrcfadi_nm.
 * @param[in] eqn    input equation object defining Riccati Equation
 * @param[in] opt    input options object to configure the Newton iteration and the internal ADI
 * @param[out] status   status object containing all details about the iteration
 * @param[out] Zout     factor \f$ Z \f$ of the solution \f$ X \approx ZZ^T \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_lrnm function is a wrapper for \ref mess_lrcfadi_nm.
 *
 * @sa mess_lrcfadi_nm
 */
int mess_lrnm(mess_equation eqn, mess_options opt, mess_status status, mess_matrix Zout)
{
    return mess_lrcfadi_nm(eqn, opt, status, Zout);
}

/**
 * @brief Solve an algebraic Riccati Equation using the low-rank-Cholesky factor Newton method.
 * @param[in] eqn    input equation object defining Riccati Equation
 * @param[in] opt    input options object to configure the Newton iteration and the internal ADI
 * @param[out] status   status object containing all details about the iteration
 * @param[out] Zout     factor \f$ Z \f$ of the solution \f$ X \approx ZZ^T \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_lrcfadi_nm function solves algebraic Riccati Equations (ARE)
 * \f[ A^T X + X A - X B B^T X + C^T C = 0 \f]
 * or
 * \f[ A X + X A^T - X C^T C X + B B^T = 0 \f]
 * using the low-rank-Cholesky factor Newton method (LRCF-NM). \n
 * The solution of generalized algebraic Riccati Equations (GARE)
 * \f[ A^T X E + E^T X A - E^T X B B^T X E + C^T C = 0 \f]
 * or
 * \f[ A X E^T + E X A^T - E X C^T C X E^T + B B^T = 0 \f]
 * is also supported. \n
 * This function works using a previously created @ref mess_equation object which defines all necessary operations
 * for matrices. Internally a special @ref mess_equation object is created to handle the internal ADI iteration. \n
 * Details about the algorithm can be found in \cite Saa09 and \cite BenS10 .
 *
 * The necessary equation objects are generated using \ref mess_equation_riccati or \ref mess_equation_griccati
 * for example.
 *
 * The iteration is controlled by the corresponding entries in the given options object.
 *
 * \sa mess_equation_st
 * \sa mess_lrcfadi_adi
 * \sa mess_lrcfadi_lrcfnm
 */
int mess_lrcfadi_nm(mess_equation eqn, mess_options opt, mess_status status, mess_matrix Zout){
    MSG_FNAME(__func__);
    int ret = 0;
    int stop_res2 = 0, stop_rel = 0, stop_user = 0, have_mass_matrix = 0;
    mess_int_t it, i, projected=0;
    mess_matrix G, Kold, K, W, Zold, Z, deltaK, deltaKnew;
    mess_equation lyap;
    double res2, res2_old, te, ts;
    mess_options lyapopt;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(opt);
    mess_check_positive(opt->nm_maxit);
    mess_check_nullpointer(status);
    mess_check_nullpointer(Zout);
    mess_check_positive(opt->nm_maxit);

    if ( eqn->B == NULL || eqn->C == NULL) {
        MSG_ERROR("B or C matrix missing in the equation.\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( mess_equation_has_E(eqn)) have_mass_matrix = 1;

    /*-----------------------------------------------------------------------------
     *  Prepare
     *-----------------------------------------------------------------------------*/
    ts = mess_wtime();
    MESS_MATRIX_RESET(Zout);
    MESS_INIT_MATRICES(&G, &Kold, &K, &W, &Zold, &Z, &deltaK, &deltaKnew);

    ret = mess_options_init(&lyapopt);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_options_init);
    ret = mess_options_copy(opt,lyapopt);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_options_copy);

    //relative changes of low rank change are not tracked
    //ret = mess_vector_resize(status->rel_changes, 0);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    ret = mess_vector_resize(status->rel_changes, opt->nm_maxit+1);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    ret = mess_vector_resize(status->res2_norms, opt->nm_maxit+1);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    mess_try_alloc(status->internal_status, mess_status *, sizeof(mess_status) * (opt->nm_maxit+1));
    status->n_internal_status = opt->nm_maxit+1;
    for ( i = 0; i < opt->nm_maxit+1; i++) status->internal_status[i]=NULL;


    ret = mess_matrix_dynorm2((opt->type==MESS_OP_NONE)?eqn->B:eqn->C, &(status->res2_0));                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_dynorm2);
    MSG_INFO("using relative norm, res0 = %.e\n", status->res2_0);

    /*-----------------------------------------------------------------------------
     *  Setup the Lyapunov operator
     *-----------------------------------------------------------------------------*/
    ret = mess_equation_init(&lyap);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_init);
    if ( opt->type == MESS_OP_NONE) {
        ret =  mess_equation_stable(lyap, lyapopt , eqn, NULL, eqn->C );        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_equation_stable);
    } else {
        ret =  mess_equation_stable(lyap, lyapopt , eqn, eqn->B, NULL );        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_equation_stable);
    }

    if (opt->nm_singleshifts != 0 ){
        if (lyapopt->adi_shifts_p == NULL || lyapopt->adi_shifts_p->dim == 0 ) {
            ret = mess_parameter(lyap, lyapopt, status);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_parameter);
            if (opt->adi_output){
                MSG_INFO("ADI parameter.\n");
                mess_vector_printshort(lyapopt->adi_shifts_p);
            }
            lyap->ApEINV.clear(lyap);
        }
    }

    /*-----------------------------------------------------------------------------
     *  Begin the Iteration
     *-----------------------------------------------------------------------------*/
    res2_old = res2 = 0;

    /*-----------------------------------------------------------------------------
     *  Main Iteration loop
     *-----------------------------------------------------------------------------*/
    for(it = 0 ; it < opt->nm_maxit; it++){
        ret = mess_status_init(&(status->internal_status[it]));                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_status_init);

        /*-----------------------------------------------------------------------------
         *  Build G = [ C' K ] and set right hand side of lyapunov equation
         *-----------------------------------------------------------------------------*/
        //ret = mess_matrix_init(&G);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        MESS_MATRIX_RESET(G);
        {
            if ( opt->type == MESS_OP_NONE ) {
                ret = mess_matrix_convert (eqn->B , G, MESS_DENSE);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
            } else {
                mess_matrix Gtmp;
                ret = mess_matrix_init(&Gtmp);                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_ctranspose(eqn->C, Gtmp);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
                ret = mess_matrix_convert(Gtmp, G, MESS_DENSE);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
                mess_matrix_clear(&Gtmp);
            }
        }

        if ( it > 0 ) {
            //update feedback
            if ( opt->type == MESS_OP_NONE) {
                ret = mess_matrix_addcols(G, K);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_addcols);
                lyap->B = K;
            } else {
                mess_matrix tmp;
                mess_matrix_init(&tmp);
                mess_matrix_ctranspose(K,tmp);
                ret = mess_matrix_addcols(G, tmp);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_addcols);
                mess_matrix_clear(&tmp);
                lyap->K = K;
            }
        } else if(opt->K0){
            //first step with initial feedback
            MSG_INFO("Use given initial feedback.\n");
            ret = mess_matrix_copy(opt->K0,K);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            if ( opt->type == MESS_OP_NONE) {
                ret = mess_matrix_addcols(G, K);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_addcols);
                lyap->B = K;
            } else {
                lyap->K = K;
            }
        } else{
            //first step no initial feedback
        }

        ret = mess_equation_set_rhs2 (lyap, G);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_set_rhs2);

        if ( it >0 ) {
            ret = mess_matrix_clear(&Zold);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_clear);
            Zold = Z;
            Z = NULL;
            ret = mess_matrix_init(&Z);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        }

        /*-----------------------------------------------------------------------------
         *  store Kold
         *-----------------------------------------------------------------------------*/
        {
            if(it>0){
                ret = mess_matrix_copy(K,Kold);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            }else if(opt->K0 && it==0){
                //copy initial feedback, if there is one
                //ret = mess_matrix_copy(opt->K0,deltaKnew);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
                ret = mess_matrix_copy(opt->K0,Kold);                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            }else{
                //copy initial feedback, if there is one
                //initialize matrix of same format as K
                if (opt->type == MESS_OP_TRANSPOSE){
                    ret = mess_matrix_alloc(Kold, eqn->B->cols, eqn->B->rows,eqn->B->nnz, MESS_DENSE, MESS_REAL);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
                }else{
                    ret = mess_matrix_alloc(Kold, eqn->C->cols, eqn->C->rows,eqn->C->nnz, MESS_DENSE, MESS_REAL);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
                }
            }
        }

        /*-----------------------------------------------------------------------------
         *  Solve the internal lyapunov equation
         *-----------------------------------------------------------------------------*/
        if ( opt->nm_singleshifts == 0) {
            ret = mess_parameter(lyap, lyapopt, status->internal_status[it]);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_parameter);

            if (opt->adi_output) {
                if(lyapopt->adi_shifts_p){
                    MSG_INFO("ADI parameter.\n");
                    mess_vector_printshort(lyapopt->adi_shifts_p);
                }else{
                    MSG_INFO("No shift parameter available.\n");
                }
            }
            lyap->ApEINV.clear(lyap);
        }
        ret = mess_lrcfadi_adi(lyap, lyapopt, status->internal_status[it], Z);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_lrcfadi_adi);

        /*-----------------------------------------------------------------------------
         * Projection
         *-----------------------------------------------------------------------------*/
        projected = 0;
        if (opt->nm_gpStep > 0 && it % opt->nm_gpStep == 0 ){
            ret = mess_lrcfadi_galerkin(eqn, opt, have_mass_matrix ? MESS_EQN_GRICCATI:MESS_EQN_RICCATI, Z);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_lrcfadi_galerkin);
            projected = 1;
        }

        /*-----------------------------------------------------------------------------
         *  compute new feedback K for the next step
         *-----------------------------------------------------------------------------*/
        if ( have_mass_matrix  == 0) {
            mess_matrix tmp =NULL;
            ret = mess_matrix_init(&tmp);                                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
            if ( opt->type == MESS_OP_NONE ) {
                ret = mess_matrix_multiply(MESS_OP_HERMITIAN, Z, MESS_OP_HERMITIAN, eqn->C, tmp);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
                ret = mess_matrix_multiply(MESS_OP_NONE, Z, MESS_OP_NONE, tmp, K);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            } else {
                ret = mess_matrix_multiply(MESS_OP_HERMITIAN, eqn->B, MESS_OP_NONE, Z, tmp);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
                ret = mess_matrix_multiply(MESS_OP_NONE, tmp, MESS_OP_HERMITIAN, Z, K);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
            }
            mess_matrix_clear(&tmp);
        } else {
            mess_matrix tmp  = NULL;
            mess_matrix tmp2 = NULL;
            MESS_INIT_MATRICES(&tmp,&tmp2);
            if (opt->type == MESS_OP_NONE ) {
                ret = mess_matrix_multiply(MESS_OP_HERMITIAN, Z, MESS_OP_HERMITIAN, eqn->C, tmp);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
                ret = mess_matrix_multiply(MESS_OP_NONE, Z, MESS_OP_NONE, tmp, tmp2);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
                ret = mess_equation_E_apply(eqn,MESS_OP_NONE, tmp2, K);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply);
            } else {
                ret = mess_matrix_multiply(MESS_OP_HERMITIAN, Z, MESS_OP_NONE, eqn->B, tmp);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
                ret = mess_matrix_multiply(MESS_OP_NONE, Z, MESS_OP_NONE, tmp, tmp2);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
                ret = mess_equation_E_apply(eqn, MESS_OP_HERMITIAN, tmp2,tmp);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply);
                ret = mess_matrix_ctranspose(tmp,K);                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
            }
            MESS_CLEAR_MATRICES(&tmp,&tmp2);
        }

        /*-----------------------------------------------------------------------------
         *  use line search  if equation was not projected
         *-----------------------------------------------------------------------------*/
        if(it==0){
            //set initial value to W
            ret = mess_matrix_copy(lyap->RHS, W);                                                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);
            //set old feedback change to 0 or -opt->nm_K0 in the first iteration
            ret = mess_matrix_copy(Kold,deltaK);                                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            ret = mess_matrix_scale(-1.0,deltaK);                                                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_scale);
        }

        //compute new feedback change
        ret = mess_matrix_copy(Kold,deltaKnew);                                                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_add(1.0,K,-1.0,deltaKnew);                                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);

        //use line search only if no projection was done and residual is large
        if(opt->nm_linesearch && !projected && (it==2)){
            double lambda;

            //perform line search, first iteration use right hand side, not G, because init_rhs was maybe called by lradi
            ret = mess_lrcfadi_nm_linesearch(W, lyapopt->W, deltaK, deltaKnew,&lambda);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_lrcfadi_nm_linesearch);

            //lambda is 1.0 means perform no line search
            if( fabs(lambda-1.0)<sqrt(mess_eps())) goto nolinesearch;

            //compute new Z low rank solution factor
            ret = _scalecat(MESS_OP_NONE,sqrt(lambda),Z,sqrt(1-lambda),it?Zold:NULL);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),_scalecat);

            //compute new deltaK (Feedback Change)
            ret = _scalecat(opt->type,sqrt(1-lambda), deltaK, lambda, deltaKnew);                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),_scalecat);

            //compute new W
            ret = _scalecat(MESS_OP_NONE,sqrt(1-lambda), W, sqrt(lambda), lyapopt->W);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),_scalecat);

            //compute new Feedback Matrix K
            ret = mess_matrix_add(1-lambda,Kold,lambda,K);                                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
        }else{
nolinesearch:
            //store old feedbackchange if no linesearch is used
            ret = mess_matrix_copy(deltaKnew,deltaK);                                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            ret = mess_matrix_copy(lyapopt->W,W);                                                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        }

        /*-----------------------------------------------------------------------------
         *  Calculate 2 norm residual
         *  - use Riccati Operator if projection was done in this step
         *  - use accumulated feedback if no projection was done in this step
         *-----------------------------------------------------------------------------*/
        res2_old = res2;
        if( projected == 1 ){
            //use riccati operator to compute residual, because of projection
            if (opt->residual_method == MESS_RESIDUAL_INDEFINITE){
                MSG_WARN("Projection was done during Newton-Iteration, MESS_RESIDUAL_INDEFINITE method is not feasible.\n");
            }
            ret = mess_lrcfadi_residual(eqn,opt,Z,&res2);                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_lrcfadi_residual);
        }else{
            if (opt->residual_method == MESS_RESIDUAL_INDEFINITE){
                //use feedback change, applicable because no projection was done
                ret = mess_matrix_indefinite_dynorm2(W, deltaK, &res2);                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_indefinite_dynorm2);
            }else if (opt->residual_method == MESS_RESIDUAL_SPECTRAL){
                ret = mess_lrcfadi_residual(eqn,opt,Z,&res2);                                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_lrcfadi_residual);
            }else{
                MSG_ERROR("Unsupported residual method\n");
                return MESS_ERROR_ARGUMENTS;
            }
        }

#if 0
        /*-----------------------------------------------------------------------------
         *  THIS IS AN INTERNAL DOUBLE RESIDUAL CHECK. THIS PART IS ONLY COMPILED
         *  IF MESS_DEBUG IS ON. IT IS USEFULL BECAUSE
         *  LINESEARCH + INDEFINITE RICCATI RESIDUAL + GALERKIN PROJECTION
         *  + UPDATE OF THE LOW RANK FACTOR CAN EASILY LEAD TO WRONG
         *  RESIDUAL COMPUTATION DUE TO IMPLEMENTATION ERRORS.
         *  THEREFORE WE DO HERE A DOUBLE RESIDUAL COMPUTATION CHECK.
         *-----------------------------------------------------------------------------*/
        double myres2, rel_err=1, tol = 5e-2, diff_exp;
        ret = mess_lrcfadi_residual(eqn,opt,Z,&myres2);                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_lrcfadi_residual);
        rel_err = fabs(myres2-res2)/fabs(myres2);
        diff_exp = fabs(log10(myres2)-log10(res2));
        MSG_INFO("\n--------------------------------------\n"
                "DOUBLE RESIDUAL CHECK:\n"
                "                   it      ="MESS_PRINTF_INT"\n"
                "                   tol     =%e\n"
                "                   rel_err =%e\n"
                " log error exponents       =%e\n"
                "(indefinite)       res2    =%e\n"
                "(lrcfadi_residual) res2    =%e\n"
                "--------------------------------------\n",it,tol,rel_err,diff_exp,res2,myres2);

        if(rel_err>tol && diff_exp>=1 ){
            MSG_ERROR("\ntol = %e\t rel_err = %e \n"
                    "Iteration = "MESS_PRINTF_INT " Residual computed by mess_lrcfadi_residual = %e.\n"
                    "Iteration = "MESS_PRINTF_INT " Residual computed by indefinite   residual = %e.\n",tol,rel_err,it,myres2,it,res2);
            return MESS_ERROR_MISC;
        }
#endif

        //status->res2_norms->values[it] = res2/status->res2_0;
        status->res2_norms->values[it] = res2;
        status->rel_changes->values[it] = fabs(res2 - res2_old)/res2;
        if ( res2/status->res2_0 < opt->nm_res2_tol ) stop_res2 ++;

        /*-----------------------------------------------------------------------------
         *  Call the Step_debug function
         *-----------------------------------------------------------------------------*/
        if ( opt->nm_stepfunction ) {
            mess_lrcfadi_step step;
            step.it = it;
            step.res2 = res2;
            step.res2_change = status->rel_changes->values[it];
            step.rel_change = 0;
            step.stop_user = 0;
            step.Z = Z;
            opt->nm_stepfunction(opt->nm_stepfunction_aux, eqn, lyap, opt, &step);
            stop_user = step.stop_user;
        }

        if (opt->nm_output) {
            SET_COLOR;
            MSG_PRINT("\nLRCF-NM:\n");
            MSG_PRINT(" it   |   ||R(X_{k+1})||_2/||RHS||_2    | (rel.) Res-2 Change    | Inner ADI Steps\n");
            MSG_PRINT("------+---------------------------------+------------------------+----------------\n");

            for ( i = 0;  i <= it ; i++ ) {
                MSG_PRINT( "%5d |      %.8e             |     %.8e     | %5d\n",
                        (int) i ,
                        status->res2_norms->values[i]/status->res2_0,
                        status->rel_changes->values[i],
                        (int) status->internal_status[i]->it);
            }
            UNSET_COLOR;
        }

        if (stop_rel || stop_user || stop_res2 ){
            //MSG_PRINT("stop_rel  = %d\n",stop_rel);
            //MSG_PRINT("stop_user = %d\n",stop_user);
            //MSG_PRINT("stop_res  = %d\n",stop_res2);
            break;
        }

    }

    te = mess_wtime();

    ret = mess_matrix_clear(&G);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    mess_equation_clear(&lyap);

    /*-----------------------------------------------------------------------------
     *  post process
     *-----------------------------------------------------------------------------*/
    it = MESS_MIN(it,opt->nm_maxit-1);
    status->it        = it+1;
    //status->it        = it;
    status->stop_rel  = stop_rel;
    status->stop_res2 = stop_res2;
    status->stop_user = stop_user;
    status->time_all  = te - ts;
    ret = mess_vector_resize(status->rel_changes,it+1);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    ret = mess_vector_resize(status->res2_norms, it+1);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    mess_try_realloc(status->internal_status, mess_status *, sizeof(mess_status)*(it+1));
    status->n_internal_status=it+1;
    status-> res2_norm = status->res2_norms->values[it];
    if ( it > 0)
        status->res2_change = fabs(status->res2_norms->values[it]-status->res2_norms->values[it-1]);
    else
        status->res2_change = 0;

    *Zout = *Z;
    Z->values=NULL;

    mess_options_clear(&lyapopt);
    MESS_CLEAR_MATRICES(&G, &Kold, &K, &W, &Zold, &Z, &deltaK, &deltaKnew);


    /*-----------------------------------------------------------------------------
     *  print message if iteration has not convergend
     *-----------------------------------------------------------------------------*/
    if ( !(stop_rel || stop_user || stop_res2) ){
        MSG_WARN_NOBT("\n***********************************************************************\n"
                      "*    The NM-ADI iteration has not converged to the desired residual.  *\n"
                      "***********************************************************************\n");

    }


    return 0;

}

