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
 * @file lib/lrcf_adi/lradi.c
 * @brief The template formulation of the LRCF-ADI.
 * @author @koehlerm
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#define RES2_COUNT 3
#define RES_ABORT 1e+12

#define IS_REAL(X) (cimag((X))==0.0)

/**
 * @brief Wrapper/Alternative name for @ref mess_lrcfadi_adi.
 * @param[in]   eqn input   Equation object defining Lyapunov Equation
 * @param[in]   opt input   options for iteration
 * @param[out]  stat    status information about iteration
 * @param[out]   Z  output factor \f$ Z \f$ of the solution \f$ X \approx ZZ^T \f$
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_lradi function is a wrapper around \ref mess_lrcfadi_adi.
 *
 * @sa mess_lrcfadi_adi
 */
int mess_lradi(mess_equation eqn, mess_options opt, mess_status stat, mess_matrix Z)
{
    return mess_lrcfadi_adi(eqn, opt, stat, Z);
}



/**
 * @brief Solve a (generalized) Lyapunov Equation using the low-rank-Cholesky-factor ADI.
 * @param[in]   eqn      input Equation object defining Lyapunov Equation
 * @param[in]   opt  input options for iteration
 * @param[out]  stat    status information about iteration
 * @param[out]   Z  output factor \f$ Z \f$ of the solution \f$ X \approx ZZ^T \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_adi function solves a standard Lyapunov Equation
 * \f[ A X + X A^T + B B^T = 0 \f]
 * or a generalized Lyapunov Equation
 * \f[ A X E^T + E X A^T + B B^T = 0 \f]
 * using the low-rank Cholesky factor alternate direction implicit (LRCF-ADI) method. \n
 * The implementation produces a real solution factor even if complex shift parameters are used. \n
 * Details about the algorithm and its enhancements can be found in \cite BenKS12, \cite Saa09 and \cite KoeS09 . \n
 * The matrices are not passed in directly. They need to be part of a @ref mess_equation object before which can be
 * setup using \ref mess_equation_lyap or \ref mess_equation_glyap before. The selection of the
 * equation, normal or transposed, is done when this object is created.
 *
 * Before the function is called a set of iteration parameters need to be computed. Normally
 * this is done by calling \ref mess_lrcfadi_parameter. If other shifts have been computed before they need to be
 * put into the @ref mess_options.adi_shifts_p vector.
 *
 * The options for the ADI iteration like the maximum iteration number and the cancellation criteria
 * are configured by the @p opt argument. This argument is of type
 * \ref mess_options.
 *
 * If relative and absoulte residual is large than @c 1e+12 the algorithm return immediately an error.
 *
 * \sa mess_lrcfadi_lrcfadi
 * \sa mess_options
 * \sa mess_status
 *
 */
int mess_lrcfadi_adi(mess_equation eqn, mess_options opt, mess_status stat, mess_matrix Z)
{
    MSG_FNAME(__func__);

    int ret = 0 ;
    int adaptive_shifts = 0;
    int have_mass_matrix = 0;
    int stop_rel_change = 0;
    int stop_res2=0;
    int stop_res2c=0;
    int stop_user=0;
    int init_vreal_vimag = 0;
    mess_int_t it=0, ip=0, dim=0, stat_idx=0;

    mess_double_cpx_t pc=0;
    double res2 = 0, res2old = 0, res2c =0;
    double time_all_start= 0;
    double time_all_end= 0;
    double fnormz = 0, fnormv= 0;

    mess_matrix V1 = NULL;
    mess_matrix Vtmp = NULL;
    mess_matrix vreal = NULL;
    mess_matrix vimag = NULL;
    mess_operation_t op;


    /*-----------------------------------------------------------------------------
     *  Check Input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(stat);
    mess_check_nullpointer(Z);
    dim = mess_equation_dim(eqn);
    mess_check_positive(dim);

    /*-----------------------------------------------------------------------------
     *  check for feedback
     *-----------------------------------------------------------------------------*/

    if ( !mess_equation_has_ApEs(eqn)) {
        MSG_ERROR("No Shifted solve function defined.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if ( eqn->RHS == NULL) {
        MSG_ERROR("The equation does not provide a right hand side.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if ( mess_equation_has_E(eqn) ) {
        have_mass_matrix = 1;
    }

    /*-----------------------------------------------------------------------------
     *  Prepare shifts
     *-----------------------------------------------------------------------------*/

    adaptive_shifts = (opt->adi_shifts_paratype==MESS_LRCFADI_PARA_ADAPTIVE_V)|| (opt->adi_shifts_paratype==MESS_LRCFADI_PARA_ADAPTIVE_Z);
    if (!adaptive_shifts && opt->adi_shifts_p != NULL){
        mess_check_positive(opt->adi_shifts_p->dim);
        mess_check_positive(opt->adi_maxit);
        //try to convert if all entries are real
        ret = mess_vector_convert_if_real(opt->adi_shifts_p);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_convert_if_real);
        mess_vector_sort(opt->adi_shifts_p);
    }

    if (!adaptive_shifts && opt->adi_shifts_p == NULL) {
        MSG_INFO("Shift parameters for the ADI are missing. mess_parameter is called now.\n");
        ret = mess_parameter(eqn, opt, stat);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_parameter);
        //try to convert if all entries are real
        ret = mess_vector_convert_if_real(opt->adi_shifts_p);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_convert_if_real);
        //try to convert if all entries are real
        ret = mess_vector_sort(opt->adi_shifts_p);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_sort);
    }

    /*-----------------------------------------------------------------------------
     *  Get the operation mode
     *-----------------------------------------------------------------------------*/
    op = opt->type == MESS_OP_NONE ?  MESS_OP_NONE : MESS_OP_TRANSPOSE;

    /*-----------------------------------------------------------------------------
     *  Prepare the Iteration
     *-----------------------------------------------------------------------------*/
    time_all_start = mess_wtime();
    MESS_MATRIX_RESET(Z);
    if(opt->W){
        MESS_MATRIX_RESET(opt->W);
    }else{
        MESS_INIT_MATRICES(&(opt->W));
    }
    MESS_INIT_MATRICES(&V1,&Vtmp);
    ret = mess_matrix_alloc(Z,0,0,0,MESS_DENSE, MESS_REAL);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);


    //init eqn->RHS
    if(eqn->init_rhs){
        ret = eqn->init_rhs(eqn,opt);                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), eqn->init_rhs);
    }

    ret = mess_matrix_convert(eqn->RHS, opt->W,MESS_DENSE);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    ret = mess_matrix_dynorm2(opt->W, &(stat->res2_0));                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_dynorm2);

    // reset the residual history vectors
    ret = mess_vector_toreal_nowarn(stat->rel_changes);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_vector_toreal_nowarn(stat->res2_norms);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_vector_resize(stat->rel_changes,opt->adi_maxit+2);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    ret = mess_vector_resize(stat->res2_norms,opt->adi_maxit+2);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    ret = mess_vector_zeros(stat->res2_norms);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_zeros);

    /*-----------------------------------------------------------------------------
     *  Build the Multisolver and Setup other stuff
     *-----------------------------------------------------------------------------*/
    if( !adaptive_shifts){
        ret = mess_equation_ApEs_pre(eqn, opt->adi_shifts_p);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_ApEs_pre);
    }

    if ( have_mass_matrix ) {
        ret = mess_equation_E_pre(eqn);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_pre);
    }

    if (opt->adi_output){
        // MSG_PRINT(" it   |  ||V||_F/||Z||_F  |  ||L(ZZ')||_2/||RHS||_2  |   ||L(ZZ')||_2  | (rel.) Res-2 Change \n");
        MSG_PRINT(" it   |  ||V||_F/||Z||_F  |  ||L(ZZ')||_2/||RHS||_2  |   ||L(ZZ')||_2  | (rel.) Res-2 Change \n");
        MSG_PRINT("------+-------------------+--------------------------+-----------------+---------------------\n");
    }
    /*-----------------------------------------------------------------------------
     *  Main iteration
     *-----------------------------------------------------------------------------*/
    for ( it = 0; it < opt->adi_maxit; it++) {

        //compute only new shifts in the first iteration and if all shifts are used
        if ((adaptive_shifts && it==0) ||(adaptive_shifts && it>0 && ip+1>=opt->adi_shifts_p->dim )){
            // set the correct projection space
            if(it==0 ){
                ret = mess_options_set_proj_space(opt, opt->W);                                                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_options_set_proj_space);
            }else{
                if(IS_REAL(pc)){
                    if ( opt->adi_shifts_paratype == MESS_LRCFADI_PARA_ADAPTIVE_V ){
                        ret = mess_options_set_proj_space(opt, V1);                                                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_options_set_proj_space);
                    }else if (opt->adi_shifts_paratype == MESS_LRCFADI_PARA_ADAPTIVE_Z ){
                        ret = mess_options_set_proj_space(opt, Z);                                                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_options_set_proj_space);
                    }else{
                        return MESS_ERROR_NOSUPPORT;
                    }
                }else{
                    ret = mess_matrix_cat(vreal,vimag,NULL,NULL,MESS_DENSE,V1);                                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_cat);
                    if ( opt->adi_shifts_paratype == MESS_LRCFADI_PARA_ADAPTIVE_V ){
                        ret = mess_options_set_proj_space(opt, V1);                                                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_options_set_proj_space);
                    }else if (opt->adi_shifts_paratype == MESS_LRCFADI_PARA_ADAPTIVE_Z ){
                        ret = mess_options_set_proj_space(opt, Z);                                                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_options_set_proj_space);
                    }else{
                        return MESS_ERROR_NOSUPPORT;
                    }
                }
            }

            ret = mess_parameter(eqn,opt,stat);                                                                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_parameter);
            ret = mess_options_unset_proj_space(opt);                                                                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_options_unset_proj_space);

            //clear and set new shift parameters
            ret = mess_equation_ApEs_post(eqn);                                                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_ApEs_post);
            ret = mess_equation_ApEs_pre(eqn, opt->adi_shifts_p);                                                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_ApEs_pre);
            ip=-1;
        } /* Adaptive Shifts  */

        if(adaptive_shifts){
            ip++;
        }else{
            ip = it%(opt->adi_shifts_p->dim);
        }

        // select next shift
        ret = mess_vector_get(opt->adi_shifts_p,ip,&pc);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_get);

        /*-----------------------------------------------------------------------------
         *  Solve for V1
         *-----------------------------------------------------------------------------*/

        if ( op == MESS_OP_NONE){
            ret = mess_equation_ApEs_apply(eqn, MESS_OP_NONE, pc, ip, opt->W, V1);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_ApEs_apply);
        } else if(!IS_REAL(pc)){
            ip++;
            ret = mess_equation_ApEs_apply(eqn, MESS_OP_HERMITIAN, conj(pc), ip, opt->W, V1);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_ApEs_apply);
        }else{
            ret = mess_equation_ApEs_apply(eqn, MESS_OP_HERMITIAN, pc, ip, opt->W, V1);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_ApEs_apply);
        }

        if (IS_REAL(pc)) {
            /*-----------------------------------------------------------------------------
             * Real Shift parameter
             *-----------------------------------------------------------------------------*/
            ret = mess_matrix_toreal(V1);                                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_toreal);
            if (have_mass_matrix) {
                ret = mess_equation_E_apply(eqn, op, V1, Vtmp);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply);
                ret = mess_matrix_add(-2.0*creal(pc), Vtmp, 1, opt->W);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
            } else {
                ret = mess_matrix_add(-2.0*creal(pc), V1, 1, opt->W);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
            }

            /*-----------------------------------------------------------------------------
             *  Add cols to Z
             *-----------------------------------------------------------------------------*/
            ret = mess_matrix_addcols1(Z,sqrt(-2.0*creal(pc)),V1);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addcols1);

            /*-----------------------------------------------------------------------------
             *  calculate relative change by the new column block
             *-----------------------------------------------------------------------------*/
            {
                double fac= fabs(-2.0*creal(pc));
                ret = mess_matrix_normf2(V1,&fnormv);                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_normf2);
                fnormv*=fac;
                fnormz += fnormv;
                stat->rel_changes->values[stat_idx] = sqrt(fnormv/fnormz); ++stat_idx;
            }
        } else {

            /*-----------------------------------------------------------------------------
             *  Process a Complex Shift
             *-----------------------------------------------------------------------------*/
            double alpha = 2*sqrt(-creal(pc));
            double beta  = creal(pc)/cimag(pc);
            pc = conj(pc);

            if (init_vreal_vimag == 0) {
                MESS_INIT_MATRICES(&vreal,&vimag);
                init_vreal_vimag = 1;
            }
            ret = mess_matrix_realpart(V1, vreal);                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realpart);
            ret = mess_matrix_imagpart(V1, vimag);                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_imagpart);

            // Update W
            ret = mess_matrix_add(alpha*beta, vimag, alpha, vreal);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
            if ( have_mass_matrix ) {
                ret = mess_equation_E_apply(eqn, op, vreal, Vtmp);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply);
                ret = mess_matrix_add(alpha, Vtmp, 1.0, opt->W);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
            } else {
                ret = mess_matrix_add(alpha, vreal, 1.0,opt->W);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
            }

            /*-----------------------------------------------------------------------------
             *  calculate relative change by the new column block
             *-----------------------------------------------------------------------------*/
            {
                double f2norm_vreal, f2norm_vimag;
                mess_matrix_normf2(vreal,&f2norm_vreal);
                mess_matrix_normf2(vimag,&f2norm_vimag);
                fnormv = f2norm_vreal + alpha*alpha*(beta*beta+1)*f2norm_vimag;
                fnormz += fnormv;
                stat->rel_changes->values[stat_idx] = sqrt(fnormv/fnormz);  ++stat_idx;
                stat->rel_changes->values[stat_idx] = sqrt(fnormv/fnormz);  ++stat_idx;
            }


            /*-----------------------------------------------------------------------------
             *  Update Z
             *-----------------------------------------------------------------------------*/
            ret = mess_matrix_addcols(Z,vreal);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addcols);
            ret = mess_matrix_addcols1(Z,alpha*sqrt(beta*beta+1),vimag);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addcols1);

            // Skip next iteration step
            if ( opt ->adi_output ) {
#ifdef MESS64
                MSG_PRINT("%5"PRId64" |   -------------   |      --------------      |  -------------  |  ------------ \n", it);
#else
                MSG_PRINT("%5d |   -------------   |      -------------       |  ------------   |     ------------ \n", it);
#endif
            }

            it++;
        }


        /*-----------------------------------------------------------------------------
         *  Calculate Residuals
         *-----------------------------------------------------------------------------*/
        res2old = res2;
        if (opt->residual_method == MESS_RESIDUAL_INDEFINITE){
            ret = mess_matrix_dynorm2(opt->W,&res2);                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_dynorm2);
        }else if (opt->residual_method == MESS_RESIDUAL_SPECTRAL){

            /*-----------------------------------------------------------------------------
             *  Use Arpack or Arnoldi to calculate the 2-Norm residual.
             *-----------------------------------------------------------------------------*/
            mess_equation_t eqn_type    = eqn->eqn_type;
            eqn->eqn_type               = MESS_EQN_GLYAP;
            ret                         = mess_lrcfadi_residual(eqn,opt,Z,&res2);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_lrcfadi_residual);
            eqn->eqn_type               = eqn_type;
        }else{
            MSG_ERROR("Unsupported residual method\n");
            return MESS_ERROR_ARGUMENTS;
        }


        /*-----------------------------------------------------------------------------
         *  save residual to stat structure and compute stopping criteria
         *-----------------------------------------------------------------------------*/

        /*NOTE: the stat_idx counter was incremented to the next "unwritten" entry
         * during relative change compution after processing the last shift parameter*/
        stat->res2_norms->values[stat_idx-1] = res2;
        if(!IS_REAL(pc)){
            //for a complex shift log also the value before
            stat->res2_norms->values[stat_idx-2] = res2;
        }

        res2c = fabs(res2old-res2)/res2;
        if ( it > 0 && res2/stat->res2_0  < opt->adi_res2_tol){
            stop_res2++;
        }
        if ( it > 0 && res2c < opt->adi_res2c_tol){
            stop_res2c++;
        }

        /*-----------------------------------------------------------------------------
         *  Call the Step_debug function
         *-----------------------------------------------------------------------------*/
        if ( opt->stepfunction ) {
            mess_lrcfadi_step step;
            step.it = it;
            step.res2 = res2;
            step.res2_change = res2c;
            /*NOTE: the stat_idx counter was incremented to the next "unwritten" entry
             *during relative change compution after processing the last shift parameter*/
            step.rel_change = stat->rel_changes->values[stat_idx-1];
            step.stop_user = 0;
            step.Z = Z;
            opt->stepfunction(opt->stepfunction_aux, eqn, opt, &step);
            stop_user = step.stop_user;
        }


        /*-----------------------------------------------------------------------------
         *  Output
         *-----------------------------------------------------------------------------*/

        if (opt ->adi_output ){
#ifdef MESS64
            MSG_PRINT("%5" PRId64 " |  %.9e  |     %.9e      | %.9e | %.9e \n", it, stat->rel_changes->values[stat_idx-1], res2/stat->res2_0, res2, res2c);
#else
            MSG_PRINT("%5d |  %.9e  |     %.9e      | %.9e |   %.9e \n", (int) it, stat->rel_changes->values[stat_idx-1], res2/stat->res2_0, res2, res2c);
#endif
        }

        /*-----------------------------------------------------------------------------
         *  emergency stop, if residual is huge
         *-----------------------------------------------------------------------------*/
        if(res2/stat->res2_0 > RES_ABORT && res2 > RES_ABORT){
            MSG_ERROR("Abort, no convergence.\n")
                return MESS_ERROR_CONVERGE;
        }
        /*-----------------------------------------------------------------------------
         *  Check Criterias
         *-----------------------------------------------------------------------------*/
        if ( stat->rel_changes->values[stat_idx-1] < opt->adi_rel_change_tol) {
            stop_rel_change++;
        }

        if ( stop_rel_change>2 || stop_res2>0 || stop_res2c > RES2_COUNT || stop_user) {
            break;
        }
    }
    time_all_end = mess_wtime();

    //print new line
    if (opt ->adi_output ){MSG_PRINT("\n");}



    /*-----------------------------------------------------------------------------
     *  post process state information
     *-----------------------------------------------------------------------------*/
    //iteration until end or double step at the end

    stat->it = stat_idx;
    stat->stop_rel = stop_rel_change;
    stat->stop_res2 = stop_res2;
    stat->stop_res2c = stop_res2c;
    stat->stop_user = stop_user;
    stat->time_adi = time_all_end - time_all_start;
    ret = mess_vector_resize(stat->rel_changes,stat_idx);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    ret = mess_vector_resize(stat->res2_norms, stat_idx);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    stat-> res2_norm = stat->res2_norms->values[stat_idx-1];
    stat->res2_change = res2c ;
    stat->rel_change = stat->rel_changes->values[stat_idx-1];


    if ( vreal ) mess_matrix_clear(&vreal);
    if ( vimag ) mess_matrix_clear(&vimag);

    MESS_CLEAR_MATRICES(&Vtmp,&V1);

    /*-----------------------------------------------------------------------------
     *  print message if iteration has not convergend
     *-----------------------------------------------------------------------------*/
    if (!(stop_rel_change>2 || stop_res2>0 || stop_res2c > RES2_COUNT || stop_user)){
        MSG_WARN_NOBT("\n***********************************************************************\n"
                      "*    The ADI iteration has not converged to the desired residual.     *\n"
                      "***********************************************************************\n");
    }

    return 0;
}

