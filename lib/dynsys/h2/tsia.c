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
 * @file lib/dynsys/h2/tsia.c
 * @brief TSIA algorithm.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "irka_common.c"

#define MAT_INIT(A) ret = mess_matrix_init(&(A)); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);

/**
 * @brief Compute a \f$ \mathcal{H}_2 \f$ reduced model via the TSIA algorithm.
 * @param[in] orig   input original system
 * @param[in] sigmaEX    input vector of initial interpolation points (NULL if not needed)
 * @param[in] opt    input options
 * @param[out] redu     reduced system
 * @param[in] V  input state space projection matrix
 * @param[in] W      input left projection matrix
 * @param[out] status   output status
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_h2_tsia function computes a \f$ \mathcal{H}_2 \f$ reduced model via the TSIA algorithm. \n
 * If sigmaEX!=NULL and the input system is a SISO system, the first projection will be build from this
 * interpolation points. \n
 * If no sigmaEX is given \f$ V \f$ and \f$ W \f$  are used as projection. If \f$ V \f$ and \f$ W \f$ are empty,
 * random matrices are used.
 *
 */
int mess_h2_tsia (  mess_dynsys orig, mess_vector sigmaEX, mess_h2_options opt,
            mess_dynsys redu, mess_matrix V, mess_matrix W, mess_h2_status status)

{
    MSG_FNAME(__func__);
    int usesigma = 0;
    mess_matrix Ar,ArT,Br,Cr;
    mess_matrix A,B,C;
    mess_matrix BBr;
    mess_matrix CCr;
    mess_matrix Wk,Vk;
    mess_direct sylv;
    int ret = 0;
    double tol;
    double h2err=0;
    double h2orig = 0;
    mess_int_t maxit;
    mess_int_t rdim;
    mess_vector sigma, sigmaold;
    double nrmsigma = 0;
    mess_int_t i, j , us, unstable;
    double ds ;
    double ts, te;
    mess_direct Abase;
    mess_int_t count_unstable = 0 ;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(orig);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(redu);
    mess_check_nullpointer(V);
    mess_check_nullpointer(W);
    mess_check_nullpointer(status);
    if ( sigmaEX != NULL && MESS_IS_DYNSYS_SISO(orig)) {
        usesigma =1 ;
        rdim = sigmaEX->dim;
    } else {
        usesigma = 0;
        rdim = opt->rdim;
    }
    if ( !MESS_IS_DYNSYS_LTI(orig)){
        MSG_ERROR("only standard state space LTI systems allowed.\n");
        return MESS_ERROR_DYNSYS;
    }
    tol = opt->tol;
    maxit = opt->maxit;
    if ( tol < 0 ){
        MSG_ERROR("tol must be at least zero.\n");
        return MESS_ERROR_ARGUMENTS;
    }
    if ( maxit < 0){
        MSG_ERROR("maxit (= " MESS_PRINTF_INT " ) have to be non negative.\n", maxit);
        return MESS_ERROR_ARGUMENTS;
    }
    if ( rdim <=0){
        MSG_ERROR("reduced dimension doesn't fit: rdim=" MESS_PRINTF_INT "\n", rdim);
        return MESS_ERROR_ARGUMENTS;
    }
    A=orig->A;
    B=orig->B;
    C=orig->C;


    ts = mess_wtime();
    /*-----------------------------------------------------------------------------
     *  prepare, compute intial reduced order model
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_resize(status->sigmadiff, maxit+1);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
    ret = mess_vector_zeros(status->sigmadiff);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_zeros);
    ret = mess_vector_resize(status->h2err, maxit+1);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
    ret = mess_vector_zeros(status->h2err);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_zeros);

    ret = mess_direct_init(&Abase);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_direct_init);
    ret = mess_direct_lu(A, Abase);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);

    MESS_INIT_VECTORS(&sigma,&sigmaold);
    ret = mess_vector_alloc(sigma, rdim, MESS_COMPLEX);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(sigmaold, rdim, MESS_COMPLEX);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_alloc);
    MAT_INIT(Ar);
    MAT_INIT(Br);
    MAT_INIT(Cr);
    MAT_INIT(ArT);
    MAT_INIT(Wk);
    MAT_INIT(Vk);
    MAT_INIT(CCr);
    MAT_INIT(BBr);
    if ( usesigma == 1) {
        ret = mess_vector_norm2(sigmaEX,&nrmsigma);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
        MSG_INFO("Use given interpolation points\n");
        ret = __constructVWmat(A, Abase, NULL ,B,C,sigmaEX, Vk,Wk);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0), __constructVWmat);
        ret = mess_matrix_biorth(Vk,Wk,V,W);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_biorth);
        ret = mess_vector_copy(sigmaEX, sigmaold);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy);
    } else {
        if ( V->cols == rdim && V->rows == A->rows ) {
            mess_check_same_size(V,W);
            MSG_INFO("Use given projection matrices\n");
            ret = mess_matrix_copy(V,Vk);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            ret = mess_matrix_copy(W,Wk);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        } else {
            MSG_INFO("Use random projection matrices\n");
            ret = mess_matrix_rand_dense(Wk,A->rows, rdim,MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_rand_dense);
            ret = mess_matrix_rand_dense(Vk,A->rows, rdim,MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_rand_dense);
        }
        ret = mess_matrix_biorth(Vk,Wk,V,W);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_biorth);
    }
    // project the system
    ret = __project_A(W,A,V,Ar);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), __project_A);
    ret = __project_Bmat(W,B,Br);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), __project_Bmat);
    ret = __project_Cmat(V,C,Cr);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), __project_Cmat);
    ret = mess_eigen_eig(Ar,sigma,NULL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eig);
    mess_vector_scalee(-1.0,sigma);         // FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_scalee);
    ret = mess_vector_sort(sigma);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_sort);
    if ( usesigma == 1)  {
        ret=  mess_vector_diffnorm(sigma,sigmaold, &ds);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_diffnorm);
        ds =ds /nrmsigma;
        status->sigmadiff->values[0] = ds;
    } else {
        ret = mess_vector_norm2(sigma,&nrmsigma);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
        status->sigmadiff->values[0] = 0;
        ds = 0;
    }

    if (opt->output ){
        printf(" step |  sigma-sigmaold  |   h2err          | stable               \n");
        printf("------+------------------+------------------+----------------------\n");
    }


    /*-----------------------------------------------------------------------------
     *  inital H2 Error
     *-----------------------------------------------------------------------------*/
    if (opt->calc_h2err == MESS_H2NORM_FULL ){
        ret = mess_h2_error_internal(A,B,C,NULL,Ar,Br,Cr,NULL,&h2err  );        FUNCTION_FAILURE_HANDLE(ret,(ret!=0 && ret!=MESS_ERROR_CONVERGE), mess_h2_error);
        status->h2err->values[0] = h2err;
    } else if ( opt->calc_h2err == MESS_H2NORM_UPDATE ) {
        ret = mess_h2_norm(orig,&h2orig);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_h2_norm);
        h2orig = h2orig*h2orig;
        h2err = 0;
        status->h2err->values[0] = h2err;
    }
    /*-----------------------------------------------------------------------------
     *  check stability
     *-----------------------------------------------------------------------------*/
    us = 0 ; unstable = 0;
    for ( j = 0 ; j < sigma->dim; j++) {    if ( creal ( sigma->values_cpx[j]) <=0) us++;   }
    if ( us > 0 ){
        MSG_INFO("reduced system is unstable.\n");
        unstable = 1;
    }

    if ( opt->output) {
#ifdef MESS64
        printf(" %4" PRId64 " | %16.10e | %16.10e | %d \n", (mess_int_t)0, ds,h2err, unstable?0:1);
#else
        printf(" %4d | %16.10e | %16.10e | %d \n", (mess_int_t)0, ds,h2err, unstable?0:1);
#endif
    }


    /*-----------------------------------------------------------------------------
     *  main loop
     *-----------------------------------------------------------------------------*/
    i = 1;
    while ( i <=maxit ){
    // for ( i = 1; i<=maxit ; i++ ){
        h2err = 0 ;
        ret = mess_matrix_ctranspose(Ar,ArT);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);

        // init sylvester solver
        ret = mess_direct_init(&sylv);                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
        ret = mess_direct_create_sylvester_sparsedense(A,NULL,NULL,ArT,sylv);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sylvester_sparsedense);
        // solve
        ret = mess_matrix_multiply(MESS_OP_NONE, B, MESS_OP_HERMITIAN, Br, BBr);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_HERMITIAN, C, MESS_OP_NONE, Cr, CCr);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_scale(-1,CCr);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);
        ret = mess_direct_solvem(MESS_OP_NONE,sylv,BBr,Vk);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
        ret = mess_direct_solvem(MESS_OP_HERMITIAN,sylv,CCr,Wk);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
        ret = mess_matrix_biorth(Vk,Wk,V,W);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_biorth);

        mess_direct_clear(&sylv);

        if ( opt->calc_h2err == MESS_H2NORM_UPDATE ){
            double h2min = 0;
            double h2kopp  =0;
            mess_matrix Q12BR;
            mess_matrix BQ12BR;
            MAT_INIT(Q12BR);
            MAT_INIT(BQ12BR);
            ret = mess_h2_norm_internal(Ar,Br,Cr,NULL, &h2min );        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_h2_norm_internal);
            h2min = h2min*h2min;

            ret = mess_matrix_multiply(MESS_OP_NONE, Wk, MESS_OP_NONE, Br, Q12BR);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_multiply(MESS_OP_HERMITIAN, B, MESS_OP_NONE, Q12BR, BQ12BR);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_trace(BQ12BR, &h2kopp);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_trace);

            h2err = h2orig + h2min + 2*h2kopp;
            h2err = sqrt(h2err);
            if ( isnan(h2err)) h2err=0;

            // printf("h2orig: %lg,\t h2min: %lg,\t h2kopp: %lg\n", h2orig, h2min, h2kopp);

            if ( opt->output ) printf("      |                  | %.10e\n", h2err);
            status->h2err->values[i-1] = h2err;

            h2err = 0;
            mess_matrix_clear(&Q12BR);
            mess_matrix_clear(&BQ12BR);
        }

        // project
        ret = __project_A(W,A,V,Ar);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_A);
        ret = __project_Bmat(W,B,Br);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_Bmat);
        ret = __project_Cmat(V,C,Cr);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_Cmat);


        // calc sigma
        ret = mess_vector_copy(sigma, sigmaold);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
        ret = mess_eigen_eig(Ar,sigma,NULL);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eig);
        mess_vector_scalee(-1.0,sigma);
        ret = mess_vector_sort(sigma);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_sort);

        // calc full h2 error if wanted

        if (opt->calc_h2err == MESS_H2NORM_FULL ){
            ret = mess_h2_error_internal(A,B,C,NULL,Ar,Br,Cr,NULL,&h2err);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0 && ret!=MESS_ERROR_CONVERGE), mess_h2_error);
        }

        /*-----------------------------------------------------------------------------
         *  check stability
         *-----------------------------------------------------------------------------*/
        us = 0 ;
        unstable = 0;
        for ( j = 0 ; j < sigma->dim; j++) {
            if ( creal ( sigma->values_cpx[j]) <=0) us++;
        }
        if ( us > 0 ){
            MSG_INFO("reduced system is unstable.\n");
            unstable = 1;
        }


        // pole change critera
        ret =mess_vector_diffnorm(sigma, sigmaold,&ds);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_diffnorm);
        ds = ds / nrmsigma;
        status->sigmadiff->values[i] = ds;
        status->h2err->values[i] = h2err;

        if (opt->output) {
#ifdef MESS64
            printf(" %4" PRId64 " | %16.10e | %16.10e | %d \n", i, ds,h2err, unstable?0:1 );
#else
            printf(" %4d | %16.10e | %16.10e | %d \n", i, ds,h2err, unstable?0:1 );
#endif
        }


        status->it=i;

        // printf("it = " MESS_PRINTF_INT " \t maxit = % ld \n", i, maxit );
        if ( unstable && i >= maxit) {
            count_unstable++;
            if ( count_unstable < 10 ) {
                maxit ++ ;
                // resize status vectors
                //
                ret = mess_vector_resize(status->sigmadiff, maxit+1);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);
                ret = mess_vector_resize(status->h2err, maxit+1);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);
            }
            i++;
            continue;
        }
        i++;

        if ( ds < tol ){
            status->cancel_sigma = 1;
            break;
        }
    }
    if ( status->cancel_sigma !=1){
        status->cancel_it = 1;
    }

    /*-----------------------------------------------------------------------------
     *  final h2 error calculation
     *-----------------------------------------------------------------------------*/
    if ( opt->calc_h2err == MESS_H2NORM_UPDATE ){
        ret = mess_matrix_ctranspose(Ar,ArT);                                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_ctranspose);
        ret = mess_direct_init(&sylv);                                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
        ret = mess_direct_create_sylvester_sparsedense(A,NULL,NULL, ArT,sylv);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_direct_create_sylvester_sparsedense);
        // solve
        ret = mess_matrix_multiply(MESS_OP_NONE, B, MESS_OP_HERMITIAN, Br, BBr);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_HERMITIAN, C, MESS_OP_NONE, Cr, CCr);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_scale(-1,CCr);                                                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_scale);
        ret = mess_direct_solvem(MESS_OP_NONE,sylv,BBr,Vk);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
        ret = mess_direct_solvem(MESS_OP_HERMITIAN,sylv,CCr,Wk);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
        mess_direct_clear(&sylv);

        double h2min = 0;
        double h2kopp  =0;
        mess_matrix Q12BR;
        mess_matrix BQ12BR;
        MAT_INIT(Q12BR);
        MAT_INIT(BQ12BR);
        ret = mess_h2_norm_internal(Ar,Br,Cr,NULL, &h2min);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_h2_norm_internal);
        h2min = h2min*h2min;

        ret = mess_matrix_multiply(MESS_OP_NONE, Wk, MESS_OP_NONE, Br, Q12BR);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_multiply);
        ret = mess_matrix_multiply(MESS_OP_HERMITIAN, B, MESS_OP_NONE, Q12BR, BQ12BR);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_multiply);
        // mess_matrix_print(BQ12BR);
        ret = mess_matrix_trace(BQ12BR, &h2kopp);                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_trace);

        h2err = h2orig + h2min + 2*h2kopp;
        h2err = sqrt(h2err);
        if ( isnan(h2err)) h2err=0;
        if ( opt->output ) printf("      |                  | %.10e\n", h2err);
        status->h2err->values[i-1]=h2err;
        status->finalh2 = h2err;
        mess_matrix_clear(&Q12BR);
        mess_matrix_clear(&BQ12BR);

    } else if ( !opt->calc_h2err && opt->calc_finalh2){
        MSG_INFO("calculate finale h2err\n");
        ret = mess_h2_error_internal(A,B,C,NULL,Ar,Br,Cr,NULL,&h2err );
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0 && ret!=MESS_ERROR_CONVERGE), mess_h2_error);
        if ( opt->output) printf("final h2-err: %lg\n", h2err);
        status->finalh2 = h2err;
    } else {
        status->finalh2 = status->h2err->values[maxit] ;
    }
    status->finalsigma = status->sigmadiff->values[maxit];

    ret = mess_dynsys_lti(redu, Ar,Br,Cr);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_dynsys_lti);

    te = mess_wtime();
    status->time = te-ts;

    /*-----------------------------------------------------------------------------
     *  finalize
     *-----------------------------------------------------------------------------*/
    mess_vector_clear(&sigma);
    mess_vector_clear(&sigmaold);
    mess_matrix_clear(&Wk);
    mess_matrix_clear(&Vk);
    mess_matrix_clear(&ArT);
    mess_matrix_clear(&CCr);
    mess_matrix_clear(&BBr);
    mess_direct_clear(&Abase);
    return 0;
}       /* -----  end of function mess_h2_tsia  ----- */

