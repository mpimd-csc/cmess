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
 * @file lib/dynsys/h2/irka_biorth.c
 * @brief IRKA algorithm (biorthonormal version).
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

#define STARTTIME(X) {MSG_WARN(X); ts=mess_wtime();}
#define ENDTIME {te=mess_wtime(); MSG_WARN("end: %lg s\n", te-ts);}

/**
 * @brief  Compute a \f$ \mathcal{H}_2 \f$ reduced model via the IRKA algorithm (biorthonormal version).
 * @param[in] orig   input original system
 * @param[in] sigma       input vector of initial interpolation points
 * @param[in] opt    input options
 * @param[out] reduced   reduced system
 * @param[out] V    state space projection matrix
 * @param[out] W    left projection matrix
 * @param[out] status   output status
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_h2_irka_biorth function computes a \f$ \mathcal{H}_2 \f$ reduced model with the IRKA algorithm. \n
 * The correction equation is solved via the biorthonormalization.
 *
 */
int  mess_h2_irka_biorth (  mess_dynsys orig, mess_vector sigma, mess_h2_options opt,
        mess_dynsys reduced, mess_matrix V, mess_matrix W, mess_h2_status status)
{
    MSG_FNAME(__func__);
    mess_vector b = NULL, c=NULL;
    mess_vector br = NULL, cr=NULL;
    mess_vector sigmaold;
    mess_int_t n = 0;
    int ret = 0;
    mess_int_t nsigma =0 ;
    mess_matrix Vk,Wk;
    mess_matrix A,B,C, E;
    mess_matrix Ar, Br, Cr;
    mess_int_t it = 0;
    mess_int_t i = 0 ;
    int conv = 0;
    double ds;
    double h2err = 0;
    mess_direct Asolver;
    double te,ts;
    mess_int_t maxit;
    double tol;
    int stepdebug = 0;
    double nrmsigma = 0;
    int girka=0;

    // double *vds;
    // double *vh2;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(orig);
    mess_check_nullpointer(sigma);
    mess_check_nullpointer(reduced);
    mess_check_nullpointer(V);
    mess_check_nullpointer(W);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(status);
    tol = opt->tol;
    if (tol<0) {
        MSG_ERROR("tol must at least be zero.\n");
        return MESS_ERROR_ARGUMENTS;
    }
    maxit = opt->maxit;

    if (mess_error_level > 2) {
        mess_matrix_printinfo(orig->A);
        mess_matrix_printinfo(orig->B);
        mess_matrix_printinfo(orig->C);
    }

    if ( maxit < 0){
        MSG_ERROR("maxit (= " MESS_PRINTF_INT " ) have to be non negative.\n", maxit);
        return MESS_ERROR_ARGUMENTS;
    }
    if ( MESS_IS_DYNSYS_LTI(orig)) {
        MSG_INFO("using standard state space system.\n");
        girka =0;
    } else if (MESS_IS_DYNSYS_GLTI(orig)){
        MSG_INFO("using generalized state space system.\n");
        girka =1;
    } else {
        MSG_ERROR("input system type not supported.\n");
        return MESS_ERROR_DYNSYS;
    }
    if (!MESS_IS_DYNSYS_SISO(orig)){
        MSG_ERROR("SISO system required.\n");
        return MESS_ERROR_DIMENSION;
    }
    A=orig->A;
    B=orig->B;
    C=orig->C;
    if (girka) E=orig->E; else E =  NULL;
    n = A->rows;

    if ( opt->stepdebug != NULL){
        MSG_INFO("using step debug function.\n");
        stepdebug = 1;
    }
    ret = mess_vector_norm2(sigma, &nrmsigma);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ts = mess_wtime();
    MESS_INIT_VECTORS(&b,&c);
    ret = mess_vector_alloc(b,n,MESS_REAL);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(c,n,MESS_REAL);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_matrix_getcol(B,0,b);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
    ret = mess_matrix_getrow(C,0,c);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getrow);
    ret = mess_vector_sort(sigma);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_sort);
    ret = mess_vector_tocomplex(sigma);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);

    nsigma=sigma->dim;

    MESS_INIT_VECTORS(&sigmaold,&br,&cr);
    MESS_INIT_MATRICES(&Ar,&Br,&Cr,&Vk,&Wk);
    ret = mess_vector_alloc(sigmaold, sigma->dim, MESS_COMPLEX);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(br,nsigma,MESS_REAL);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(cr,nsigma,MESS_REAL);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_matrix_alloc(Br,nsigma, 1, nsigma, MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(Cr,1, nsigma, nsigma, MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    // basic solver for A, only in the case of standard state space systems
    Asolver = NULL;
    // output status
    ret = mess_vector_resize(status->h2err, maxit+1);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    ret = mess_vector_resize(status->sigmadiff, maxit+1);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);


    /*-----------------------------------------------------------------------------
     *  build initial V and W
     *-----------------------------------------------------------------------------*/
    ret = __constructVW(A,Asolver, E ,b,c,sigma, Vk,Wk);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), __constructVW);

    // MSG_WARN("start biorth\n");ts=mess_wtime();
    if ( girka) {
        ret = mess_matrix_gbiorth(E, Vk,Wk,V,W);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_gbiorth);
    } else {
        ret = mess_matrix_biorth(Vk,Wk,V,W);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_biorth);
    }
    // te=mess_wtime();MSG_WARN("end biorth: %lg s\n", te-ts);

    it = 0;
    conv = 0;

    if (opt->output ){
        MSG_PRINT(" step |  sigma-sigmaold  |   h2err\n");
        MSG_PRINT("------+------------------+------------------\n");
    }

    do{
        // project A
        // STARTTIME("start project A\n");
        ret  = __project_A(W,A,V,Ar);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), __project_A);
        // ENDTIME;

        // compute new poles
        ret = mess_vector_copy(sigma, sigmaold);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_copy);
        ret = mess_eigen_eig(Ar,sigma,NULL);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eig);
        mess_vector_scalee(-1.0,sigma);

        /*-----------------------------------------------------------------------------
         *  testcode
         *-----------------------------------------------------------------------------*/
        mess_vector_tocomplex( sigma) ;
        for ( i = 0; i < sigma->dim; i++) {
            if ( creal(sigma->values_cpx[i]) < 0){
                // MSG_PRINT("Change the signs. \n");
                sigma->values_cpx[i] = - sigma->values_cpx[i];
            }
        }

        ret = mess_vector_sort(sigma);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_sort);
        // pole change critera
        ret =mess_vector_diffnorm(sigma, sigmaold,&ds);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_diffnorm);
        ds = ds / nrmsigma;



        /*-----------------------------------------------------------------------------
         *  H2-Error
         *-----------------------------------------------------------------------------*/
        if (opt->calc_h2err) {
            ret = __project_B(W,b,br);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0),__project_B);
            ret = __project_C(V,c,cr);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0),__project_C);
            ret = mess_matrix_setcol(Br, 0, br);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_setcol);
            ret = mess_matrix_setrow(Cr,0, cr);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setrow);

            ret = mess_h2_error_internal(A,B,C,E,Ar,Br,Cr,NULL,&h2err );        FUNCTION_FAILURE_HANDLE(ret,(ret!=0 && ret!=MESS_ERROR_CONVERGE), mess_h2_error);
        } else {
            h2err = 0;
        }
        status->h2err->values[it] = h2err;
        status->sigmadiff->values[it] = ds;

        if (opt->output) {
#ifdef MESS64
            MSG_PRINT(" %4" PRId64 "  | %.10e | %.10e\n", it, ds,h2err);
#else
            MSG_PRINT(" %4d | %.10e | %.10e\n", it, ds,h2err);
#endif

        }
        if ( stepdebug ){
            opt->stepdebug(opt->stepdebug_data, it, ds, h2err);
        }

        if ( it >=maxit ) {
            status->cancel_it  =1;
            break;
        }

        if ( ds < tol ) {
            status->cancel_sigma = 1;
            conv = 1;
        } else {
            ret = __constructVW(A,Asolver,E,b,c,sigma, Vk,Wk);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0), __constructVW);

            if ( girka ){
                ret = mess_matrix_gbiorth(E,Vk,Wk,V,W);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_gbiorth);
            } else {
                ret = mess_matrix_biorth(Vk,Wk,V,W);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_biorth);
            }
        }
        it++;
    }while (!conv);

    /*-----------------------------------------------------------------------------
     *  finalize
     *-----------------------------------------------------------------------------*/
    ret = __project_B(W,b,br);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_B);
    ret = __project_C(V,c,cr);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __project_C);
    ret = mess_matrix_setcol(Br, 0, br);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_setcol);
    ret = mess_matrix_setrow(Cr,0, cr);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setrow);

    if ( !opt->calc_h2err && opt->calc_finalh2){
        MSG_INFO("calculate finale h2err\n");
        ret = mess_h2_error_internal(A,B,C,E,Ar,Br,Cr,NULL,&h2err); FUNCTION_FAILURE_HANDLE(ret,(ret!=0 && ret!=MESS_ERROR_CONVERGE), mess_h2_error);
        if ( opt->output) MSG_PRINT("final h2-err: %lg\n", h2err);
        status->finalh2 = h2err;
    } else {
        status->finalh2 = status->h2err->values[it-conv];
    }
    status->finalsigma = status->sigmadiff->values[it-conv];
    te = mess_wtime();
    status->time = te -ts;

    ret = mess_dynsys_lti(reduced, Ar, Br, Cr);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dynsys_lti);

    if ( Asolver != NULL) mess_direct_clear(&Asolver);
    mess_vector_clear(&b);
    mess_vector_clear(&c);
    mess_vector_clear(&br);
    mess_vector_clear(&cr);
    // mess_matrix_clear(&V);
    mess_matrix_clear(&Vk);
    // mess_matrix_clear(&W);
    mess_matrix_clear(&Wk);
    mess_vector_clear(&sigmaold);
    //    mess_free(vds);
    //    mess_free(vh2);
    return 0;
}       /* -----  end of function mess_h2_irka_biorth  ----- */

