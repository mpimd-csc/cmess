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
 * @file lib/lrcf_adi/para.c
 * @brief Parameter compuation for ADI processes.
 * @author @koehlerm
 */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"

#ifdef _OPENMP_H
#include <omp.h>
#endif



/**
 * @brief Compute the maximum magnitude of the rational ADI function.
 * @param[in] p  input vector of ADI parameters
 * @param[in] lp  input length of the \f$ p \f$ vector
 * @param[in] set  input vector representing the discrete set
 * @param[in] lset  input length of the \f$ set \f$ vector
 * @param[out] max_r maximal magnitude of the rational ADI function over the set
 * @param[out] ind   index - maximum is attained for \f$ set(ind) \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_lp_s function computes the maximum magnitude of the rational ADI function over
 * a discrete subset of the left complex half plane. \n
 * The function is only an auxilary routine for \ref mess_lrcfadi_lp_mnmx.
 *
 * Details can be found in \cite Pen98 or \cite Saa09 .
 *
 * If a @c NULL-pointer is given @ref MESS_ERROR_NULLPOINTER is returned.
 *
 * \sa mess_lrcfadi_lp_mnmx
 */
int mess_lrcfadi_lp_s(mess_double_cpx_t *p, mess_int_t lp, mess_double_cpx_t *set, mess_int_t lset, double *max_r, mess_int_t *ind)
{
    MSG_FNAME(__func__);
    mess_int_t i, j;
    mess_double_cpx_t x;
    double rr;
    if ( p == NULL || set == NULL || max_r == NULL || ind == NULL) {
        MSG_ERROR("p or set or max_r or ind points to NULL");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( lp < 0 || lset < 0) {
        MSG_ERROR("lp < 0 \n");
        return MESS_ERROR_ARGUMENTS;
    }

    *max_r = -1;
    *ind = 0;
    for ( i = 0; i < lset; i++){
        x = set[i];
        rr = 1;
        for ( j = 0; j < lp; j++){
            rr = rr * cabs(p[j]-x)/cabs(p[j]+x);
        }
        if ( rr > *max_r){
            *max_r = rr;
            *ind = i;
        }
    }
    return 0;
}


/**
 * @brief Compute ADI parameters solving the ADI minimax problem suboptimal.
 * @param[in] rw  input vector containing numbers in the open left half plane, which
 *               approximate the spectrum of the corresponding matrix, e.g. a set of Ritz values \n
 *               (The set must be closed with respect to complex conjugation.)
 * @param[in] lrw  input length of the \f$ rw \f$ vector
 * @param[in,out] adi_shifts_l0 desired number of shift parameters (length \f$ (rw)  >= adi_shifts_l0 \f$) \n
 *               (The algorithm delivers either \f$ adi_shifts_l0 \f$ or \f$ adi_shifts_l0+1 \f$ parameters!)
 * @param[out] p   a \f$ adi_shifts_l0 \f$- or \f$ adi_shifts_l0+1 \f$-vector of suboptimal ADI parameters
 * @param[out] lp length of \f$ p \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_lp_mnmx function computes ADI parameters \f$ p \f$ by solving the rational
 * min-max problem
 * \f[ \min\limits_{p_j \in \mathbb{R}, j=1 , \cdots , lrw} \max\limits_{ \lambda \in rw}
 *     \Bigg \vert \prod\limits_{j=1}^{lrw} \frac{p_j - \lambda}{p_j + \lambda} \Bigg \vert\f]
 * suboptimal. \n
 * The delivered parameter set is closed under complex conjugation.
 *
 * Details can be found in \cite Pen98 or \cite Saa09 .
 *
 */
int mess_lrcfadi_lp_mnmx(mess_double_cpx_t *rw, mess_int_t lrw,  mess_int_t adi_shifts_l0, mess_double_cpx_t *p, mess_int_t *lp)
{
    MSG_FNAME(__func__);
    double max_rr = 0;
    double max_r;
    mess_double_cpx_t p0 = 0;
    mess_int_t i, j=0;
    mess_int_t ppos = 0;
    int ret = 0;

    if ( rw == NULL || p == NULL) {
        MSG_ERROR("rw or p or lp points to NULL.\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( lrw < 0) {
        MSG_ERROR("lrw < 0\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( lrw < adi_shifts_l0) {
        MSG_ERROR("length(rw) must be at least adi_shifts_l0\n");
        return MESS_ERROR_ARGUMENTS;
    }

    for (i=0; i < lrw; i++){
        ret = mess_lrcfadi_lp_s(&rw[i], 1, rw, lrw, &max_r, &j);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_lrcfadi_lp_s);
        if ( i == 0){
            max_rr = max_r;
            p0 = rw[i];
        }
        if ( max_r < max_rr ){
            p0 = rw [i];
            max_rr =max_r;
        }
    }
    if ( cimag(p0) != 0.0){
        p [ppos++] = p0;
        p [ppos++] = conj(p0);
    } else {
        p [ppos++] = p0;
    }

    mess_lrcfadi_lp_s(p, ppos, rw, lrw, &max_r, &i);

    while (ppos < adi_shifts_l0) {
        p0 = rw[i];
        if ( cimag(p0) != 0.0){
            p [ppos++] = p0;
            p [ppos++] = conj(p0);
        } else {
            p [ppos++] = p0;
        }
        mess_lrcfadi_lp_s(p, ppos, rw, lrw, &max_r, &i);
    }
    *lp = ppos;
    return 0;
}



/* Structure to pass the parameters to the Arnoldi procedures */
typedef struct {
    mess_equation eqn;
    mess_options opt;
    int have_mass;
} parameter_mvp_call ;

/**
 * @brief Operator times vector for the Arnoldi process.
 * @param[in] aux   input pointer to the data of the operator
 * @param[in] op    input (not used)
 * @param[in] in    input vector
 * @param[out] out  output vector
 *
 * The @ref mvp_parameter function computes out = Operator*in.  The in is not guaranteed
 * to be used read only.
 *
 */
static int mvp_parameter(void *aux, mess_operation_t op, mess_vector in, mess_vector out){
    parameter_mvp_call * argn = ( parameter_mvp_call *) aux;
    mess_equation eqn = argn->eqn;
    mess_options opt = argn->opt;

    if ( argn->have_mass  ) {
        mess_equation_A_apply_vector (eqn,    opt->type, in, out);
        mess_equation_Es_apply_vector(eqn,    opt->type, out, in);
        mess_vector_copy(in,out);
    }  else {
        mess_equation_A_apply_vector (eqn,opt->type, in,out);
    }
    return 0;
}

/**
 * @brief Inverse Operator times vector for the Arnoldi process.
 * @param[in] aux   input pointer to the data of the operator
 * @param[in] op    input (not used)
 * @param[in] in    input vector
 * @param[out] out  output vector
 *
 * The @ref mvp_parameterinv function computes Operator*out = in.  The in is not guaranteed
 * to be used read only.
 *
 */
static int mvp_parameterinv(void *aux, mess_operation_t op, mess_vector in, mess_vector out){
    parameter_mvp_call * argn = ( parameter_mvp_call *) aux;
    mess_equation eqn = argn->eqn;
    mess_options opt = argn->opt;

    if ( argn->have_mass ) {
        mess_equation_E_apply_vector (eqn, opt->type, in, out);
        mess_equation_As_apply_vector(eqn, opt->type, out, in);
        mess_vector_copy(in,out);
    }  else {
        mess_equation_As_apply_vector(eqn, opt->type, in, out);
    }
    return 0;
}


/**
 * @brief Compute a set of Ritz values from a given equation.
 * @param[in]   eqn      input equation object from which the Ritz values are computed
 * @param[in]   opt       input options for solver and shift parameter computation
 * @param[out]  ritz     vector containing Ritz values of the operator
 * @return zero on success or a non zero error code
 *
 * The @ref mess_lrcfadi_getritzvalues function computes a set of Ritz values of a given equation object. \n
 * Depending on the parameters @ref mess_options.adi_shifts_arp_p  and @ref mess_options.adi_shifts_arp_m  two Krylov subspaces are computed
 * using an Arnoldi method. These parameters set the dimension of the desired subspaces with respect to the operator
 * (arp_p) and its inverse (arp_m). \n
 * The Ritz values are combined in the output vector \f$ ritz \f$.
 *
 * In general this function is similar to \ref mess_eigen_eigs but works with a more abstract equation object and its
 * defined operators instead of matrices.
 *
 * \sa mess_lrcfadi_parameter
 * \sa mess_eigen_arnoldi_template
 *
 */
int mess_lrcfadi_getritzvalues ( mess_equation eqn, mess_options opt,  mess_vector ritz )
{
    MSG_FNAME(__func__);
    mess_int_t have_adi_shifts_b0 = 0, have_mass = 0;
    double * wr1, *wr2, *wi1, *wi2;
    mess_matrix Vp,Hp,Vm,Hm;
    mess_vector sv1,sv2;
    mess_mvpcall call1, call2;
    parameter_mvp_call mvp_para1,mvp_para2;
    mess_int_t n, lwork, info, one = 1, ritz_pos = 0, i = 0;
    double dummy = 0;
    double *work;
    int ret = 0;
    mess_int_t dim;

    /*-----------------------------------------------------------------------------
     *  Check Input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(ritz);
    dim = mess_equation_dim(eqn);
    mess_check_positive(dim);

    if ( opt->adi_shifts_b0 !=NULL ){
        have_adi_shifts_b0 = 1;
    }
    if ( opt->adi_shifts_arp_p >= dim) {
        MSG_ERROR("opt->adi_shifts_arp_p must be smaller than the dimension of the equation.\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( opt->adi_shifts_arp_m >= dim) {
        MSG_ERROR("opt->adi_shifts_arp_m must be smaller than the dimension of the equation!\n");
        return MESS_ERROR_ARGUMENTS;
    }

    have_mass = (mess_equation_has_E(eqn) && mess_equation_has_Es(eqn));

    /*-----------------------------------------------------------------------------
     *  generate Solver for the A operator
     *-----------------------------------------------------------------------------*/
    ret = mess_equation_A_pre(eqn);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_pre);
    ret = mess_equation_As_pre(eqn);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_As_pre);
    if ( have_mass ){
        MSG_INFO("compute solver for the mass matrix\n");
        ret = mess_equation_E_pre(eqn);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_pre);
        ret = mess_equation_Es_pre(eqn);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_Es_pre);
    }

    /*-----------------------------------------------------------------------------
     *  prepare eigenvalue computation
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(wr1,double *,sizeof(double)*(MESS_MAX(opt->adi_shifts_arp_p,opt->adi_shifts_arp_m)+1));
    mess_try_alloc(wi1, double *,sizeof(double)*(MESS_MAX(opt->adi_shifts_arp_p,opt->adi_shifts_arp_m)+1));
    mess_try_alloc(wr2, double *,sizeof(double)*(MESS_MAX(opt->adi_shifts_arp_p,opt->adi_shifts_arp_m)+1));
    mess_try_alloc(wi2, double *,sizeof(double)*(MESS_MAX(opt->adi_shifts_arp_p,opt->adi_shifts_arp_m)+1));
    ret = mess_matrix_init(&Vp);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&Hp);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&Vm);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
    ret = mess_matrix_init(&Hm);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);

    /*-----------------------------------------------------------------------------
     *  Arnoldi on the operator
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_init(&sv1);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sv1, dim, MESS_REAL);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    if (have_adi_shifts_b0){
        ret = mess_vector_copy(opt->adi_shifts_b0, sv1);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
    }else {
        ret = mess_vector_ones(sv1);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);
    }

    // call1.mvp = mvp_parameter;
    mvp_para1.eqn = eqn;
    mvp_para1.opt = opt;
    mvp_para1.have_mass = have_mass ;

    // call1.data = &mvp_para1;
    ret = mess_mvpcall_operator(&call1, dim, MESS_REAL, mvp_parameter, &mvp_para1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator );

    ret = mess_eigen_arnoldi_template(call1, opt->adi_shifts_arp_p, sv1, Hp, Vp );  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template);

    n =  Hp->rows-1;
    lwork = (MESS_MAX(opt->adi_shifts_arp_p,opt->adi_shifts_arp_m)+1) * 4;
    mess_try_alloc(work, double *, sizeof(double) * lwork);
    if ( opt->adi_shifts_arp_p!=0) {
        F77_GLOBAL(dgeev,DGEEV)("N","N", &n, Hp->values, &(Hp->rows), wr1, wi1, &dummy,&one, &dummy, &one, work, &lwork, &info);
        if ( info != 0) {
            MSG_ERROR("dgeev returned with info = " MESS_PRINTF_INT "\n", info);
            ret = info;
        }
    }
    mess_free(work);
    mess_vector_clear(&sv1);

    /*-----------------------------------------------------------------------------
     *  Arnoldi on the inverse operator
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_init(&sv2);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sv2, dim, MESS_REAL);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    if (have_adi_shifts_b0){
        ret = mess_vector_copy(opt->adi_shifts_b0, sv2);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
    }else {
        ret = mess_vector_ones(sv2);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);
    }

    // Call the Arnoldi Process
    mvp_para2.eqn = eqn;
    mvp_para2.opt = opt;
    mvp_para2.have_mass = have_mass ;


    // call2.data = &mvp_para2;
    ret = mess_mvpcall_operator(&call2, dim, MESS_REAL, mvp_parameterinv, &mvp_para2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator );

    ret = mess_eigen_arnoldi_template(call2, opt->adi_shifts_arp_m, sv2, Hm, Vm );      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template);

    n =  Hm->rows-1;
    lwork = (MESS_MAX(opt->adi_shifts_arp_p,opt->adi_shifts_arp_m)+1) * 4;
    mess_try_alloc(work, double *, sizeof(double) * lwork);
    if ( opt->adi_shifts_arp_m!=0){
        F77_GLOBAL(dgeev,DGEEV)("N","N", &n, Hm->values, &(Hm->rows), wr2, wi2, &dummy,&one, &dummy, &one, work, &lwork, &info);
        if ( info != 0) {
            MSG_ERROR("dgeev returned with info = " MESS_PRINTF_INT "\n", info);
            ret = info;
        }
    }
    mess_free(work);
    mess_vector_clear(&sv2);

    mess_mvpcall_clear(&call1);
    mess_mvpcall_clear(&call2);

    /*-----------------------------------------------------------------------------
     *  Collect the Eigenvalues
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_tocomplex(ritz);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
    ret = mess_vector_resize(ritz, (Hp->rows)+(Hm->rows)-2);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_vector_resize );
    ritz_pos = 0 ;
    for (i = 0; i < Hp->rows-1; i++) {
        ritz->values_cpx[ritz_pos++] = wr1[i] + wi1[i]*I;
    }
    for (i = 0; i < Hm->rows-1; i++) {
        ritz->values_cpx[ritz_pos++] = 1.0/(wr2[i] + wi2[i]*I);
    }

    /*-----------------------------------------------------------------------------
     *  cleanup the arnoldi process
     *-----------------------------------------------------------------------------*/
    mess_free(wr1);
    mess_free(wi1);
    mess_free(wr2);
    mess_free(wi2);
    MESS_CLEAR_MATRICES(&Hm,&Vm,&Hp,&Vp);


    return 0;
}       /* -----  end of function mess_lrcfadi_getritzvalues  ----- */

/**
 * @brief Compute parameters for the low-rank Cholesky factor ADI iteration.
 * @param[in] eqn    input equation object defining the desired equation
 * @param[in,out] opt   options controlling the process
 * @param[out] stat  status object for returning some states
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_parameter function computes parameters for the LRCF-ADI process. \n
 * Depending on the given mess_options.adi_shifts_paratype field in @p opt minmax, real minmax or Wachpress parameters are computed.
 * The type is defined using the @ref mess_parameter_t enumeration. \n
 * More details about shift parameters and their computation can be found in \cite Saa09  and \cite Pen98 .
 *
 * The parameter computation is controlled by 4 values from the Options object:
 * * mess_options.adi_shifts_l0 field in @p opt the desired number of shift parameters,
 * * mess_options.adi_shifts_arp_p field in @p opt sets the number of steps in the Arnoldi process with respect to the operator,
 * * mess_options.adi_shifts_arp_m field in @p opt sets the number of steps in the Arnoldi process with respect to the inverse of the operator,
 * * mess_options.adi_shifts_paratype field in @p opt sets the parameter type.
 * The necessary Ritz-values are computed via \ref mess_lrcfadi_getritzvalues.
 *
 * The @ref mess_options.adi_shifts_paratype is one of following enumeration values:
 * <ul>
 * <li> \ref MESS_LRCFADI_PARA_MINMAX - compute (complex) parameters out of the rational min-max problem.
 * <li> \ref MESS_LRCFADI_PARA_MINMAX_REAL - compute (real) parameters out of the rational min-max problem.\n
 * (They are basically the same  as \ref MESS_LRCFADI_PARA_MINMAX but the imaginary part is neglected.)
 * <li> \ref MESS_LRCFADI_PARA_WACHSPRESS - compute Wachspress parameters (\cite Wac95 ).
 * <li> \ref MESS_LRCFADI_PARA_ADAPTIVE_V - compute adaptive shift parameters (\cite BenKS14).
 * <li> \ref MESS_LRCFADI_PARA_ADAPTIVE_Z - compute adaptive shift parameters (\cite BenKS14).
 * </ul>
 * \sa mess_lrcfadi_para
 * \sa mess_lrcfadi_getritzvalues
 *
 */
int mess_lrcfadi_parameter(mess_equation eqn, mess_options opt, mess_status stat) {
    MSG_FNAME(__func__);
    mess_double_cpx_t *rw2;
    mess_int_t i;
    mess_int_t lp = 0;
    mess_double_cpx_t *p;
    mess_int_t rw2_pos = 0;
    int ret = 0;
    mess_vector ritz;
    mess_int_t dim;


    /*-----------------------------------------------------------------------------
     *  Check Input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(stat);

    mess_check_true(mess_equation_has_A(eqn));
    mess_check_positive((dim=mess_equation_dim(eqn)));


    /*-----------------------------------------------------------------------------
     *  Adaptive Shifts
     *-----------------------------------------------------------------------------*/
    if (opt->adi_shifts_paratype == MESS_LRCFADI_PARA_ADAPTIVE_V || opt->adi_shifts_paratype == MESS_LRCFADI_PARA_ADAPTIVE_Z ){


        /*-----------------------------------------------------------------------------
         *  some specifics checks for adaptive shifts
         *-----------------------------------------------------------------------------*/

        if(!mess_equation_has_A(eqn)){
            MSG_ERROR("AX.apply is necessary for adaptive shift strategy.\n");
            return MESS_ERROR_ARGUMENTS;
        }

        if(opt->proj_space==NULL){
            MSG_INFO("Please provide an project space by mess_options opt for adaptive shift parameter computation.\n");
            return 0;
        }



        /*-----------------------------------------------------------------------------
         *  get projection space and computed projected shfits
         *-----------------------------------------------------------------------------*/
        mess_vector proj_shifts;
        mess_matrix proj_orth, proj_temp, proj_A, proj_E, rand_rhs;
        MESS_INIT_VECTORS(&proj_shifts);
        mess_vector_alloc(proj_shifts,0,MESS_COMPLEX);
        MESS_INIT_MATRICES(&proj_orth,&proj_temp,&proj_A,&rand_rhs);

        mess_int_t have_mass_matrix=0, retry_random=0;

        ret = mess_equation_A_pre(eqn);                                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_pre);
        if ( mess_equation_has_E(eqn) ) {
            have_mass_matrix = 1;
            MESS_INIT_MATRICES(&proj_E);
            ret = mess_equation_E_pre(eqn);                                                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_equation_E_pre);
        }

retry_randomrhs_projection:
        //compute projected spectrum (independent of op type)
        ret = mess_matrix_orth(opt->proj_space,proj_orth);                                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_orth);
        ret = mess_equation_A_apply(eqn,MESS_OP_NONE,proj_orth,proj_temp);                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_equation_A_apply);
        ret = mess_matrix_multiply(MESS_OP_TRANSPOSE,proj_orth,MESS_OP_NONE,proj_temp,proj_A);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        if(have_mass_matrix){
            ret = mess_equation_E_apply(eqn,MESS_OP_NONE,proj_orth,proj_temp);                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_equation_E_apply);
            ret = mess_matrix_multiply(MESS_OP_TRANSPOSE,proj_orth,MESS_OP_NONE,proj_temp,proj_E);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
            ret = mess_eigen_eigg(proj_A, proj_E,proj_shifts,NULL,NULL,NULL);                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_eigen_eigg);
        }else{
            ret = mess_eigen_eig(proj_A,proj_shifts,NULL);                                             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_eigen_eig);
        }

        /*-----------------------------------------------------------------------------
         *  postprocess shift parameters
         *-----------------------------------------------------------------------------*/
        //filter instable shift values
        ret = mess_vector_filter_stable(proj_shifts);                                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_filter_stable);
        ret = mess_vector_convert_if_real(proj_shifts);                                                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_convert_if_real);
        ret = mess_vector_sort(proj_shifts);                                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_sort);
        if(proj_shifts->dim==0){
            MSG_WARN("Adaptive Shift parameter strategy found no stable shifts. opt->adi_shifts_p will be not changed.\n");
            if(opt->adi_shifts_p==NULL || opt->adi_shifts_p->dim==0){
                //retry one time with random rhs, increase number of cols every time
                if(retry_random<10){
                    mess_int_t rows = opt->proj_space->rows, cols = opt->proj_space->cols;
                    MESS_MATRIX_RESET(rand_rhs);
                    ret = mess_matrix_rand_dense(rand_rhs,rows,cols*(retry_random+1),MESS_REAL);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_rand_dense);
                    ret = mess_options_unset_proj_space(opt);                                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_options_unset_proj_space);
                    ret = mess_options_set_proj_space(opt,rand_rhs);                                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_options_set_proj_space);
                    MSG_WARN("Retry to get adaptive shifts with random projection space. rows= "MESS_PRINTF_INT" cols= "MESS_PRINTF_INT".\n",rand_rhs->rows,rand_rhs->cols);
                    retry_random++;
                    goto retry_randomrhs_projection;
                }
                MESS_CLEAR_VECTORS(&proj_shifts);
                MESS_CLEAR_MATRICES(&proj_A,&proj_temp,&proj_orth,&rand_rhs);
                if(have_mass_matrix){MESS_CLEAR_MATRICES(&proj_E);}
                MSG_ERROR("No shifts available for adaptive shift parameter strategy.\n");
                return MESS_ERROR_CONVERGE;
            }else{
                MSG_WARN("No stable shifts found by adaptive shift parameter strategy. Take the last ones.\n");
                MESS_CLEAR_VECTORS(&proj_shifts);
                MESS_CLEAR_MATRICES(&proj_A,&proj_temp,&proj_orth,&rand_rhs); if(have_mass_matrix){MESS_CLEAR_MATRICES(&proj_E);}
                return 0;
            }
        }else{
            if(opt->adi_shifts_p==NULL){
                ret = mess_vector_init(&(opt->adi_shifts_p));                                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
                ret = mess_vector_alloc((opt->adi_shifts_p),0,MESS_COMPLEX);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
            }
            ret = mess_vector_copy(proj_shifts,opt->adi_shifts_p);                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy);
        }
        MESS_CLEAR_VECTORS(&proj_shifts);
        MESS_CLEAR_MATRICES(&proj_A,&proj_temp,&proj_orth,&rand_rhs); if(have_mass_matrix){MESS_CLEAR_MATRICES(&proj_E);}
        return 0;
    }


    /*-----------------------------------------------------------------------------
     *  checks for other shift parameter strategy
     *-----------------------------------------------------------------------------*/
    mess_check_true(mess_equation_has_As(eqn));
    if ((opt->adi_shifts_l0 < 0) || (opt->adi_shifts_arp_p < 0) || (opt->adi_shifts_arp_m < 0)) {
        MSG_ERROR("opt->adi_shifts_l0 < 0 or opt->adi_shifts_arp_p < 0 or opt->adi_shifts_arp_m < 0\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( opt->adi_shifts_arp_p + opt->adi_shifts_arp_m  < 2 * opt->adi_shifts_l0) {
        MSG_ERROR("opt->adi_shifts_arp_p + opt->adi_shifts_arp_m >= 2*opt->adi_shifts_l0 required.\n");
        return MESS_ERROR_ARGUMENTS;
    }
    /*-----------------------------------------------------------------------------
     *  Adjust options
     *-----------------------------------------------------------------------------*/
    if (opt->type == MESS_OP_HERMITIAN ) opt->type = MESS_OP_TRANSPOSE;

    /*-----------------------------------------------------------------------------
     *  Compute the Ritz Values
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_init(&ritz);                                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(ritz,opt->adi_shifts_arp_p+opt->adi_shifts_arp_m,MESS_COMPLEX);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_lrcfadi_getritzvalues(eqn, opt, ritz);                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_lrcfadi_getritzvalues);
    mess_try_alloc(p, mess_double_cpx_t* ,sizeof(mess_double_cpx_t) * (opt->adi_shifts_l0+1));
    mess_try_alloc(rw2, mess_double_cpx_t* ,sizeof(mess_double_cpx_t) * (ritz->dim));
    mess_vector_tocomplex(ritz);

    /*-----------------------------------------------------------------------------
     *  Clean up is done by the @ref mess_equation_clear function
     *-----------------------------------------------------------------------------*/
    stat->unstable = 0;
    rw2_pos=0;
    for ( i = 0; i < ritz->dim; i++){
        mess_double_cpx_t rw = ritz->values_cpx[i];
        if ( creal (rw) < 0) {
            rw2[rw2_pos++] = rw;
        } else {
            MSG_WARN_NOBT("Non-Stable Ritzvalue: %lg + %lg i\n", creal(rw), cimag(rw));
            //MSG_PRINT("Non-Stable Ritzvalue: %lg + %lg i\n", creal(rw), cimag(rw));
            stat->unstable = 1;
        }
    }
    mess_vector_clear(&ritz);

    /*-----------------------------------------------------------------------------
     *  compute the parameters
     *-----------------------------------------------------------------------------*/
    if ( opt->adi_shifts_paratype == MESS_LRCFADI_PARA_MINMAX || opt->adi_shifts_paratype == MESS_LRCFADI_PARA_MINMAX_REAL ){

        /*-----------------------------------------------------------------------------
         *  Heuristic Parameters
         *-----------------------------------------------------------------------------*/
        ret = mess_lrcfadi_lp_mnmx(rw2, rw2_pos, opt->adi_shifts_l0, p, &lp);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_lrcfadi_lp_mnmx);
        // copy final values.
        if ( opt->adi_shifts_paratype == MESS_LRCFADI_PARA_MINMAX_REAL ) {
            if ( opt ->adi_shifts_p == NULL) {
                ret = mess_vector_init(&(opt->adi_shifts_p));                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_alloc((opt->adi_shifts_p), lp, MESS_REAL);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
            } else {
                ret = mess_vector_toreal(opt->adi_shifts_p);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
                ret = mess_vector_resize(opt->adi_shifts_p, lp);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
            }
            for (i = 0; i < lp; i++) {
                opt->adi_shifts_p->values[i] = creal(p[i]);
            }
        }else {
            if (opt ->adi_shifts_p == NULL) {
                ret = mess_vector_init(&(opt->adi_shifts_p));                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_alloc((opt->adi_shifts_p), lp, MESS_COMPLEX);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
            } else {
                ret = mess_vector_tocomplex(opt->adi_shifts_p);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
                ret = mess_vector_resize(opt->adi_shifts_p, lp);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
            }
            for (i = 0; i < lp; i++) {
                opt->adi_shifts_p->values_cpx[i] = p[i];
            }
        }
    } else if ( opt->adi_shifts_paratype  == MESS_LRCFADI_PARA_WACHSPRESS){
        /*-----------------------------------------------------------------------------
         *  Wachspress Parameters
         *-----------------------------------------------------------------------------*/
        double a, b, alpha, min, max, amax;
        double paras[50];
        mess_int_t l;

        min = -creal(rw2[0]);
        max = -creal(rw2[0]);
        amax = cimag(rw2[0])/creal(rw2[0]);
        for ( i = 1 ; i < rw2_pos; i++){
            if ( min > -creal(rw2[i]))
                min = -creal(rw2[i]);
            if ( max < -creal(rw2[i]))
                max = -creal(rw2[i]);
            if ( amax < cimag(rw2[i])/creal(rw2[i]))
                amax = cimag(rw2[i])/creal(rw2[i]);
        }
        a = min;
        b = max;
        alpha = atan (amax);
        double lx = opt->adi_shifts_l0;
        ret = mess_lrcfadi_para_wachspress(a, b, alpha, lx, paras, &l);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_lrcfadi_para_wachspress);
        if ( opt->adi_shifts_p == NULL ) {
            ret = mess_vector_init(&(opt->adi_shifts_p));                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
            ret = mess_vector_alloc((opt->adi_shifts_p),l, MESS_REAL);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        } else {
            ret = mess_vector_toreal(opt->adi_shifts_p);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
            ret = mess_vector_resize(opt->adi_shifts_p,l);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        }
        for (i = 0; i < l; i++){
            opt->adi_shifts_p->values[i] = paras[i];
        }
    } else {
        MSG_ERROR("unknown parameter type\n");
        ret = MESS_ERROR_ARGUMENTS;
    }

    mess_free(rw2);
    mess_free(p);

    return ret;
}


/**
 * @brief Wrapper around @ref mess_lrcfadi_parameter except that mess_equation.parameter is preferred (if available).
 * @param[in] eqn       input equation object defining the desired equation
 * @param[in] opt       options controlling the process
 * @param[out] stat     status object for returning some states
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_parameter function calls mess_equation.parameter, if it not points to @c NULL.
 * Otherwise the @ref mess_lrcfadi_parameter is called.
 * @attention Do not call @ref mess_parameter inside mess_equation.parameter.
 *
 * @sa mess_lrcfadi_para
 * @sa mess_lrcfadi_getritzvalues
 * @sa mess_equation.parameter
 *
 */
int mess_parameter(mess_equation eqn, mess_options opt, mess_status stat) {
    MSG_FNAME(__func__);
    if(eqn->parameter!=NULL){
        return eqn->parameter(eqn,opt,stat);
    }else{
        MSG_ERROR("equation does not provide a parameter function, this is mandatory!\n");
        return MESS_ERROR_NULLPOINTER;
    }
}

