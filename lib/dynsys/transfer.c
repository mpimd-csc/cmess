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
 * @file lib/dynsys/transfer.c
 * @brief Evalute the transfer function.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "cscutils/hardware.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * @brief Invert a permutation.
 * @param[in] p  input permutation
 * @param[in] n  input length of the permutation
 * @return inverted permutation
 *
 * The @ref pinv function returns the inverted permution of p.
 *
 */
static mess_int_t *pinv (mess_int_t const *p,mess_int_t n)
{
    MSG_FNAME(__func__);
    mess_int_t k, *pinv ;
    if (!p) return (NULL) ;
    mess_try_alloc2(pinv, mess_int_t *, n*sizeof (mess_int_t)) ;
    if (!pinv) return (NULL) ;
    for (k = 0 ; k < n ; k++) pinv [p [k]] = k ;
    return (pinv) ;
}



/**
 * @brief Evaluate a transfer function.
 * @param[in] sys    input LTI System
 * @param[in] aa     input lower bound for the interval
 * @param[in] bb     input upper bound for the interval
 * @param[in] nsample  input number of samples
 * @param[in,out] omega output interval
 * @param[out] G    output \f$ \max(svd(G(s))) \f$
 * @param[out] Gout     pointer to a mess_matrix array of size nsample
 * @return zero on succes or a non-zero error value
 *
 * The @ref mess_dynsys_evaltransfer function evaluates the transfer function of a
 * LTI or a GLTI system on a set of logarithmic equidistant points in \f$ \left[ 10^{aa} , 10^{bb}\right]\f$. \n
 * The memory usage can be limited using the runtime configuration is OpenMP is enabled.
 * If Gout is not @c NULL, the G(s) at the interpolation points is returned, too.
 *
 *
 */
int  mess_dynsys_evaltransfer ( mess_dynsys sys, double aa, double bb, mess_int_t nsample, mess_vector omega, mess_vector G, mess_matrix *Gout )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i,j;
    mess_matrix shift;
    mess_matrix A, B, C,E;
    int use_gernalized = 0;
    mess_matrix work;
    mess_direct base;
    mess_int_t *dpos;
    mess_int_t *p,*q, *qinv;
    int perm_before = 0;
    int secondorder = 0;
    int gout = 0;
    mess_dynsys tmpsys = NULL;
    int gp = 0;
    // mess_vector b,c;

    // MSG_PRINT("nsample = " MESS_PRINTF_INT "\n",nsample);
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sys);
    mess_check_nullpointer(omega);
    mess_check_nullpointer(G);
    // mess_check_positive(nsample);
    if ( nsample == 0){
        MSG_INFO("use given points");
        gp = 1;
        nsample = omega->dim;
    }


    /*-----------------------------------------------------------------------------
     *  workaround for second order systems
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_DYNSYS_2ND(sys)) {
        MSG_INFO("use a second order system. convert it to first order.\n");
        secondorder  =1;
        ret = mess_dynsys_init(&tmpsys); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dynsys_init);
        ret = mess_dynsys_2nd_to_1st(sys, tmpsys); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_dynsys_2nd_to_1st);
        A=tmpsys->A;
        B=tmpsys->B;
        C=tmpsys->C;
        E=tmpsys->E;
    } else {
        if (!( MESS_IS_DYNSYS_LTI(sys)|| MESS_IS_DYNSYS_GLTI(sys)) ){
            MSG_ERROR("only for LTI and generalized LTI systems or 2nd order systems\n");
            return MESS_ERROR_DYNSYS;
        }
        secondorder = 0;
        A=sys->A;
        B=sys->B;
        C=sys->C;
        E=sys->E;
    }

    if ( Gout == NULL) {
        gout = 0;
    } else {
        gout = 1;
    }

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_init(&shift); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&work); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    if ( E != NULL ){
        use_gernalized = 1;
    }

    if ( gp ==0) {
        MSG_INFO("using logarithmic x scaling\n");
        ret = mess_vector_logspace10(omega,aa,bb,nsample); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_logspace10);
    }
    ret = mess_matrix_convert(A,work,MESS_CSR); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
    ret  = mess_matrix_scale(-1,work); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_scale);
    ret = mess_vector_resize(G,omega->dim); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_init);
    ret = mess_vector_toreal(G); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_toreal);

    /*-----------------------------------------------------------------------------
     *  check for diagonal elements
     *-----------------------------------------------------------------------------*/
    if ( use_gernalized == 0) {
        int problems = 0;
        mess_try_alloc(dpos, mess_int_t *, sizeof ( mess_int_t) * work->rows);
        ret =   mess_matrix_diagpos(work , dpos);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_diagpos);
        for ( i = 0; i < work->rows; i++) {
            if ( dpos[i] < 0) {
                MSG_WARN("diagonal entry in line " MESS_PRINTF_INT " doesn't exist.\n", i);
                problems ++;
            }
        }
        if ( problems > 0 ) {
            MSG_WARN("some diagonal entries don't exists. creating identity mass matrix.\n");
            perm_before = 1;
            use_gernalized = 1;
        }
        mess_free(dpos);
    }
    if ( use_gernalized != 0 ){
        if (E==NULL){
            ret = mess_matrix_eye(shift, A->rows, A->cols, MESS_CSR);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eyec);
        } else {
            ret = mess_matrix_copy(E, shift);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_copy);
        }
    }

    /*-----------------------------------------------------------------------------
     *  generating base solver
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_init(&base);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_direct_init);
    if (use_gernalized == 0){
        mess_try_alloc(dpos, mess_int_t*, sizeof(mess_int_t)*A->rows);
        ret =   mess_matrix_diagpos(work , dpos); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_diagpos);
        ret = mess_direct_lu(work, base); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);
    } else {
        mess_matrix work2;
        ret = mess_matrix_init(&work2);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_copy(work, work2);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        ret = mess_matrix_add(0.0,shift, 1, work);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_add(1.0,shift, 1, work2);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_add(0,work, 1, shift);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
        ret = mess_matrix_sort(work2);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
        ret = mess_matrix_sort(work);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
        ret = mess_matrix_sort(shift);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sort);
        ret = mess_direct_lu(work2, base); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);

        mess_matrix_clear(&work2);
        perm_before = 1;
    }

    /*-----------------------------------------------------------------------------
     *  get pattern
     *-----------------------------------------------------------------------------*/
    mess_matrix L=NULL, U=NULL;
    ret = mess_matrix_init(&L);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&U);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

    // get L
    ret  = mess_direct_getL(base, L);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_getL);
    if ( !MESS_IS_CSR(L) ) {
        MSG_ERROR ("at the moment only CSR in base\n");
        return MESS_ERROR_STORAGETYPE;
    }

    // get U
    ret = mess_direct_getU(base,U);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_getU);
    if ( !MESS_IS_CSR(U) ) {
        MSG_ERROR ("at the moment only CSR in base\n");
        return MESS_ERROR_STORAGETYPE;
    }

    // get permutations
    mess_try_alloc ( p, mess_int_t *, sizeof(mess_int_t) * work->rows);
    mess_try_alloc ( q, mess_int_t *, sizeof(mess_int_t) * work->cols);
    ret = mess_direct_getpermp(base, p);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0) , mess_direct_getpermp);
    ret = mess_direct_getpermq(base, q);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_direct_getpermq);
    qinv=pinv(q,work->cols);
    mess_direct_clear(&base);


    // perm_before = 1;
    if (perm_before) {
        MSG_INFO("permute system\n");
        ret = mess_matrix_perm(work, p, q);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_perm);
        if ( use_gernalized == 0) {
            ret =   mess_matrix_diagpos(work , dpos);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_diagpos);
        } else {
            // MSG_PRINT("perm shift...\n");
            ret= mess_matrix_perm(shift, p,q);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_perm);
        }
        // if ( MESS_IS_DENSE(B) ) MSG_PRINT("B is dense.\n");
        // else MSG_PRINT("B is sparse\n");
        // ret  = mess_matrix_perm(B, p, NULL);
        // FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_perm);
    }


    /*-----------------------------------------------------------------------------
     *  adjust OpenMP things
     *-----------------------------------------------------------------------------*/
#ifdef _OPENMP
    {
        size_t total_ram, free_ram;
        mess_int_t size;
        int max_threads;
        int os_num_threads;
        if( getenv("OMP_NUM_THREADS") != NULL){
            os_num_threads = atoi(getenv("OMP_NUM_THREADS"));
        } else {
            os_num_threads = omp_get_num_procs();
        }
        // MSG_INFO("os_num_threads: %d\n", os_num_threads);
        int percent = 80;
        csc_hardware_memory(&total_ram, &free_ram, NULL, NULL);

        if ( total_ram !=0 ) {
            max_threads  = os_num_threads +1;
            /* MSG_PRINT("MAX_THREADS: %d\n", max_threads);
               MSG_PRINT("MEM_TOTAL: %lu\n", (unsigned mess_int_t ) minfo.memtotal);
               MSG_PRINT("MEM_FREE: %lu\n", (unsigned mess_int_t)  minfo.memfree); */
            do {
                max_threads = max_threads - 1;
                size = sizeof(mess_int_t) * (L->nnz+U->nnz+2*(L->rows+1)) + max_threads*sizeof(mess_double_cpx_t)*(L->nnz+U->nnz)+max_threads*sizeof(mess_double_cpx_t)*work->nnz;
                // MSG_PRINT("Threads: %d memory usage: ", max_threads); mess_print_bytes(size); MSG_PRINT("\n");
            } while ( size > (total_ram*percent)/100 && max_threads> 1);
            if ( max_threads < 1 ) max_threads = 1;
            MSG_INFO("adjust number of threads to %d ( use: " MESS_PRINTF_INT " bytes ) \n", max_threads, size);
            omp_set_num_threads(max_threads);
        } else {
            MSG_WARN("can not get useful memory information. number of threads untouched.\n");
        }
    }
#endif
    /*-----------------------------------------------------------------------------
     *  main loop
     *-----------------------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for private(i,j) default(shared)
#endif
    for (i=0; i < nsample; i++){
        if ( mess_error_level > 2) {
            MSG_PRINT("%s: step " MESS_PRINTF_INT " / " MESS_PRINTF_INT "\r",__func__, i,nsample);
            fflush(stdout);
        }
        mess_double_cpx_t *lvalues;
        mess_double_cpx_t *uvalues;
        mess_double_cpx_t *vwork;
        mess_matrix T;
        mess_matrix CT;
        mess_vector S;

        mess_try_alloc2(vwork, mess_double_cpx_t *,work->nnz* sizeof(mess_double_cpx_t));
        mess_try_alloc2(lvalues, mess_double_cpx_t *, L->nnz*sizeof(mess_double_cpx_t));
        mess_try_alloc2(uvalues, mess_double_cpx_t *, U->nnz*sizeof(mess_double_cpx_t));
        if (use_gernalized == 0) {
            // memcpy(vwork, work->values_cpx, work->nnz * sizeof(mess_double_cpx_t));
            for ( j = 0 ; j < work->nnz; j++) vwork[j] = work->values[j];
            for ( j = 0 ; j < work->rows; j++)  vwork[dpos[j]] += omega->values[i]*I;
        } else {
            for ( j = 0; j < work->nnz; j++) vwork[j] = work->values[j] + I*(omega->values[i]*shift->values[j]);
        }

        /*-----------------------------------------------------------------------------
         *  factorize
         *-----------------------------------------------------------------------------*/
        if ( perm_before == 0)
            mess_decomp_lureuse_kernelcsr_complex(work->rows, work->cols,   vwork,  work->colptr, work->rowptr,
                    L->colptr, L->rowptr,                       // L pattern
                    U->colptr, U->rowptr,                       // U pattern
                    p, qinv,                                    // permuatations
                    lvalues, uvalues);
        else
            mess_decomp_lureuse_kernelcsr_complex(work->rows, work->cols,   vwork,  work->colptr, work->rowptr,
                    L->colptr, L->rowptr,                       // L pattern
                    U->colptr, U->rowptr,                       // U pattern
                    NULL, NULL,                                 // permuatations
                    lvalues, uvalues);


        /*-----------------------------------------------------------------------------
         *  solve
         *-----------------------------------------------------------------------------*/
        ret = mess_matrix_init(&T);                             FUNCTION_FAILURE_HANDLE_OMP(ret,(ret!=0),mess_matrix_init);
        ret = mess_matrix_init(&CT);                            FUNCTION_FAILURE_HANDLE_OMP(ret,(ret!=0),mess_matrix_init);
        ret = mess_vector_init(&S);                             FUNCTION_FAILURE_HANDLE_OMP(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_alloc(S, C->rows, MESS_REAL);         FUNCTION_FAILURE_HANDLE_OMP(ret,(ret!=0),mess_vector_alloc);
        ret = mess_solver_lusolvem_kernelcsr_complex(work->rows, lvalues, L->colptr, L->rowptr, uvalues, U->colptr, U->rowptr, /*  NULL */ p ,q,B,T);
        FUNCTION_FAILURE_HANDLE_OMP(ret,(ret!=0),mess_solver_lusolvem_kernelcsr_complex);
        ret = mess_matrix_multiply(MESS_OP_NONE, C, MESS_OP_NONE, T, CT);           FUNCTION_FAILURE_HANDLE_OMP(ret,(ret!=0),mess_matrix_multiply);
        if ( CT->rows == 1 && CT->cols == 1) {
            G->values[i] = cabs(CT->values_cpx[0]);
        } else {
            ret = mess_eigen_svd(CT,S,NULL,NULL);           FUNCTION_FAILURE_HANDLE_OMP(ret,(ret!=0),mess_eigen_svd);
            G->values[i] = fabs(S->values[0]);
        }
        if (gout){
            ret = mess_matrix_copy(CT,Gout[i]);     FUNCTION_FAILURE_HANDLE_OMP(ret,(ret!=0),mess_matrix_copy);
        }

        // mess_vector_clear(&x);
        mess_matrix_clear(&T);
        mess_matrix_clear(&CT);
        mess_vector_clear(&S);
        mess_free(vwork);
        mess_free(lvalues);
        mess_free(uvalues);
    }

    if ( secondorder == 1) {
        mess_dynsys_clear(&tmpsys);
    }
    if ( use_gernalized == 0)    mess_free(dpos);
    // ret = mess_vector_toreal(omega);
    // FUNCTION_FAILURE_HANDLE(ret,(ret=0), mess_vector_toreal);
    mess_free(qinv);
    mess_free(p);
    mess_free(q);
    mess_matrix_clear(&shift);
    mess_matrix_clear(&work);
    mess_matrix_clear(&L);
    mess_matrix_clear(&U);
    // mess_vector_clear(&b);
    // mess_vector_clear(&c);
    return 0;
}       /* -----  end of function mess_dynsys_evaltransfer  ----- */



/**
 * @brief Compute difference and relative difference between two transfer functions.
 * @param[in] nsample    input number of samples contained in the Go and Gr array
 * @param[in] Go     input array of G(i omega) of size nsample
 * @param[in] Gr     input array of Gr(i omega) of size nsample
 * @param[out] err2 two norm error
 * @param[out] rel2 relative two norm error (w.r.t. Go)
 * @param[out] errinf   infinity error
 * @param[out] relinf   relative infinity error (w.r.t. Go)
 * @param[out] diffvec  vector of the differences
 * @return zero on succes or a non-zero error value
 *
 * The @ref mess_dynsys_transfer_diff function computes various differences between two
 * discretized transfer functions. \n
 * The transfer functions are given as two array of size \c nsample. \n
 * Each array entry must be of the type \ref mess_matrix. If \f$ err2, rel2, errinf,
 * relinf \f$ or \f$ diffvec \f$ is set to @c NULL, the value will not be computed.
 *
 */
int mess_dynsys_transfer_diff ( mess_int_t nsample, mess_matrix *Go, mess_matrix *Gr, double *err2, double *rel2, double *errinf, double *relinf, mess_vector diffvec )
{
    MSG_FNAME(__func__);
    int calc_err2 = 0;
    int calc_rel2 = 0;
    int calc_errinf = 0;
    int calc_relinf = 0;
    int copy_diff = 0;
    mess_vector diff;
    mess_vector evGo;
    mess_vector S;
    mess_matrix tmp;
    mess_int_t i;
    int ret = 0;
    int siso = 0;
    double int_err2,int_errinf,ng,m;

    /*-----------------------------------------------------------------------------
     *  check inputs
     *----------------------------------------------------------------------------*/
    mess_check_positive(nsample);
    mess_check_nullpointer(Go);
    mess_check_nullpointer(Gr);
    if ( err2 != NULL) calc_err2  =1;
    if ( rel2 != NULL) calc_rel2  =1;
    if ( errinf != NULL) calc_errinf  =1;
    if ( relinf != NULL) calc_relinf  =1;
    if ( diffvec != NULL) copy_diff = 1;

    if ( Go[0]->rows != Gr[0]->rows) {
        MSG_ERROR("The matrices in the arrays Go and Gr must have the same number of rows.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( Go[0]->cols != Gr[0]->cols) {
        MSG_ERROR("The matrices in the arrays Go and Gr must have the same number of cols.\n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&diff,&evGo,&S);
    ret = mess_vector_alloc(diff, nsample, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(evGo, nsample, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(S, 1, MESS_REAL);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_alloc);
    ret = mess_matrix_init(&tmp);                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);

    if ( Go[0] ->rows == 1 && Go[0]->cols == 1) siso=1;


    /*-----------------------------------------------------------------------------
     *   evaluate
     *-----------------------------------------------------------------------------*/
    for (i=0; i < nsample; i++){
        // eval | Go |
        if (siso) {
            evGo->values[i] = cabs(Go[i]->values_cpx[0]);
            diff->values[i] = cabs(Go[i]->values_cpx[0]-Gr[i]->values_cpx[0]);
        } else {
            ret = mess_eigen_svd(Go[i],S,NULL,NULL);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_svd);
            evGo->values[i] = fabs(S->values[0]);
            ret = mess_matrix_copy(Gr[i], tmp);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
            ret = mess_matrix_addc(-1,Go[i], 1, tmp);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
            ret = mess_eigen_svd(tmp,S,NULL, NULL);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_svd);
            diff->values[i] = fabs(S->values[0]);
        }
    }

    if ( calc_err2 || calc_rel2 ){
        ret = mess_vector_norm2(diff,&int_err2);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    }
    if ( calc_err2 ) *err2 = int_err2;
    if ( calc_rel2 ){
        ret = mess_vector_norm2(evGo,&ng);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
        *rel2 = int_err2/ng;
    }
    if ( calc_errinf || calc_relinf){
        ret = mess_vector_norminf(diff,&int_errinf);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norminf);
    }
    if ( calc_errinf) *errinf = int_errinf;
    if ( calc_relinf ){
        ret = mess_vector_norminf(evGo,&m);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norminf);
        *relinf = int_errinf/m;
    }

    if ( copy_diff ) {
        ret = mess_vector_copy(diff, diffvec);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_copy);
    }
    mess_vector_clear(&S);
    mess_matrix_clear(&tmp);
    mess_vector_clear(&diff);
    mess_vector_clear(&evGo);
    return 0;
}       /* -----  end of function mess_dynsys_transfer_diff  ----- */
