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
 * @file lib/dynsys/h2/irka_para_init.c
 * @brief Initialize parameters for IRKA.
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
#include "blas_defs.h"

#ifdef _OPENMP_H
#include <omp.h>
#endif


/**
 * @brief Compare two complex numbers for qsort.
 * @param[in] pc1  input pointer to the first number
 * @param[in] pc2  input pointer to the second number
 * @return -1, 0 , 1
 * The @ref compare_complex function compares two complex numbers by the absolute value (@c cabs).
 * If the absolute value is equal, the comparison is done by the complex argument (@c carg)
 * The function is used for stdlibc qsort function.
 */
static int compare_complex(const void *pc1, const void *pc2){
    mess_double_cpx_t *c1 = (mess_double_cpx_t *) pc1;
    mess_double_cpx_t *c2 = (mess_double_cpx_t *) pc2;

    if ( cabs(*c1) < cabs(*c2)){
        return -1;
    } else if ( cabs(*c1) > cabs(*c2)){
        return 1;
    } else {
        if (carg(*c1) < carg(*c2)){
            return -1;
        }else if ( carg(*c1) > carg(*c2)){
            return 1;
        } else {
            return 0;
        }
    }
    return 0;
}

/**
 * @brief Initialize parameters for IRKA.
 * @param[in] A        input system matrix
 * @param[in] Asolver  input solver for \f$ A \f$
 * @param[in] E         input mass matrix (optional)
 * @param[in,out] r0    number of parameters
 * @param[in] kp     input number of arnoldi steps w.r.t. \f$ A \f$
 * @param[in] km     input number of arnoldi steps w.r.t. \f$ inv(A) \f$
 * @param[out] sigma    vector with the parameters
 * @return zero on succes or a non-zero error value
 *
 * The @ref mess_h2_irka_init function initializes the parameters for the IRKA algorithm.
 *
 */
int mess_h2_irka_init ( mess_matrix A, mess_direct Asolver, mess_matrix E, mess_int_t* r0, mess_int_t kp, mess_int_t km, mess_vector sigma )
{
    MSG_FNAME(__func__);
    int have_mass=0;
    mess_direct mass_solver;
    mess_vector b0 = NULL;
    mess_vector sv1=NULL;
    mess_vector sv2=NULL;
    mess_matrix Hp,Hm;
    mess_matrix Vp,Vm;
    int ret=0;
    mess_int_t lwork = 0;
    double *work = NULL;
    mess_double_cpx_t *rw;
    mess_double_cpx_t *rw2;
    double *wr1, *wi1;
    double *wr2, *wi2;
    mess_int_t r = 0,n = 0;
    int have_b0 = 1;
    mess_int_t one = 1;
    double dummy = 0;
    mess_int_t info =  0;
    mess_int_t i = 0;
    mess_int_t rw_pos= 0;
    mess_int_t rw2_pos = 0;
    double eps = mess_eps();
    mess_int_t end ;
    mess_int_t rl, rs;
    mess_double_cpx_t *sigmaL = NULL;
    mess_double_cpx_t *sigmaS = NULL;

    /*-----------------------------------------------------------------------------
     *  chech input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(Asolver);
    mess_check_nullpointer(sigma);
    mess_check_nullpointer(r0);
    r = *r0;
    if ( 2*(kp+km) < r) {
        MSG_ERROR("2*(kp+km) < r, please choose better values.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if ( E != NULL){
        MSG_INFO("generalized parameters\n");
        have_mass= 1;
        MSG_INFO("compute solver for the mass matrix\n");
        ret = mess_direct_init(&mass_solver);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
        ret = mess_direct_lu(E, mass_solver);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);
    }


    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_init(&b0);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(b0, A->rows, MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(b0);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_ones);

    MESS_INIT_MATRICES(&Vp,&Hp,&Vm,&Hm);

    mess_try_alloc(rw, mess_double_cpx_t* ,sizeof(mess_double_cpx_t) * (kp+km));
    mess_try_alloc(rw2,mess_double_cpx_t*, sizeof(mess_double_cpx_t) * (kp+km));
    mess_try_alloc(wr1,double *,sizeof(double)*(MESS_MAX(kp,km)+1));
    mess_try_alloc(wi1, double *,sizeof(double)*(MESS_MAX(kp,km)+1));
    mess_try_alloc(wr2, double *,sizeof(double)*(MESS_MAX(kp,km)+1));
    mess_try_alloc(wi2, double *,sizeof(double)*(MESS_MAX(kp,km)+1));
    // mess_try_alloc(p, mess_double_cpx_t* ,sizeof(mess_double_cpx_t) * (l0+1));

#ifdef _OPENMP_H
#pragma omp parallel sections private(work,sv1, n, lwork, sv2, info) shared (kp,km,wi1, wr1, wr2, wi2)
#endif
    {
#ifdef _OPENMP_H
#pragma omp section
#endif
        {
            // Arnoldi on A
            mess_vector_init(&sv1);
            mess_vector_alloc(sv1, A->rows, MESS_REAL);
            if (have_b0){
                mess_vector_copy(b0, sv1);
            }else {
                mess_vector_ones(sv1);
            }
            if ( have_mass ) {
                mess_eigen_arnoldig(A, NULL, NULL, mass_solver, kp, sv1, Hp, Vp);
            } else {
                mess_eigen_arnoldi(A, NULL, NULL, kp, sv1, Hp, Vp);
            }
            n =  Hp->rows-1;
            lwork = (MESS_MAX(kp,km)+1) * 4;
            mess_try_alloc( work, double *, sizeof(double) *lwork);
            F77_GLOBAL(dgeev,DGEEV)("N","N", &n, Hp->values, &(Hp->rows), wr1, wi1, &dummy,&one, &dummy, &one, work, &lwork, &info);
            if ( info != 0) {
                MSG_ERROR("dgeev returned with info = " MESS_PRINTF_INT "\n", info);
                ret = info;
            }
            mess_free(work);
            mess_vector_clear(&sv1);

        }
#ifdef _OPENMP_H
#pragma omp section
#endif
        {
            // Arnoldi on inv(A)
            mess_vector_init(&sv2);
            mess_vector_alloc(sv2, A->rows, MESS_REAL);

            if (have_b0){
                mess_vector_copy(b0, sv2);
            }else {
                mess_vector_ones(sv2);
            }
            if ( have_mass ) {
                MSG_WARN("arnoldi wrt inv(A)*E\n");
                mess_eigen_arnoldig_inv(Asolver, NULL, NULL,E , km, sv2, Hm, Vm);
            } else {
                mess_eigen_arnoldi_inv(Asolver, NULL, NULL,km, sv2, Hm, Vm);
            }
            n =  Hm->rows-1;
            lwork = (MESS_MAX(kp,km)+1) * 4;
            mess_try_alloc( work, double *, sizeof(double) *lwork);
            F77_GLOBAL(dgeev,DGEEV)("N","N", &n, Hm->values, &(Hm->rows), wr2, wi2, &dummy,&one, &dummy, &one, work, &lwork, &info);
            if ( info != 0) {
                MSG_ERROR("dgeev returned with info = " MESS_PRINTF_INT "\n", info);
                ret = info;
            }
            mess_free(work);
            mess_vector_clear(&sv2);

        }
    } // parallel section

    if (ret !=0 ) {
        return ret;
    }
    ret = 0;

    MSG_INFO("Arnoldi on A\n");
    for (i = 0; i < Hp->rows-1; i++) {
        MSG_INFO(" - eigenvalue [ " MESS_PRINTF_INT " ] = %lg + %lg i\n", i, wr1[i], wi1[i]);
        if ( wr1[i] < 0.0) {
            rw[rw_pos++] = wr1[i] + wi1[i]*I;
        } else {
            MSG_WARN("found unstable ritz value: %lg + %lg i, ignoring.\n",wr1[i],wi1[i]);
        }
    }

    MSG_INFO("Arnoldi on inv(A)\n");
    for (i = 0; i < Hm->rows-1; i++) {
        if ( wr2[i]<0.0) {
            rw[rw_pos++] = 1.0/(wr2[i] + wi2[i]*I);
        } else {
            MSG_WARN("found unstable ritz value: %lg + %lg i, ignoring.\n",wr2[i],wi2[i]);
        }
        MSG_INFO(" - eigenvalue [ " MESS_PRINTF_INT "  ] = %lg + %lg i\n", i, creal(rw[rw_pos-1]), cimag(rw[rw_pos-1]));
    }


    /*-----------------------------------------------------------------------------
     *  cleanup arnoldi
     *-----------------------------------------------------------------------------*/
    if ( have_b0 == 1){
        mess_vector_clear(&b0);
    }
    mess_matrix_clear(&Hp);
    mess_matrix_clear(&Vm);
    mess_matrix_clear(&Vp);
    mess_matrix_clear(&Hm);
    mess_free(wr1);
    mess_free(wi1);
    mess_free(wr2);
    mess_free(wi2);


    /*-----------------------------------------------------------------------------
     *  select points
     *-----------------------------------------------------------------------------*/
    qsort(rw,rw_pos,sizeof(mess_double_cpx_t), compare_complex);

    if ( ( cimag(rw[0])>eps) && (fabs(creal(rw[0])-creal(rw[1]))> 1e-13)){
        for (i=1;i<rw_pos;i++){
            rw[i-1]=rw[i];
        }
        rw_pos--;
    }
    end = rw_pos-1;
    if ( ( cimag(rw[end])>eps) && (fabs(creal(rw[end])-creal(rw[end-1]))> 1e-13)){
        rw_pos--;
    }

    rl = ceil(r/2);
    rs = r - rl;

    mess_try_alloc(sigmaL,mess_double_cpx_t *, sizeof(mess_double_cpx_t)*(rl+1));
    mess_try_alloc(sigmaS,mess_double_cpx_t *, sizeof(mess_double_cpx_t)*(rs+1));


    /*-----------------------------------------------------------------------------
     *  large ev
     *-----------------------------------------------------------------------------*/
    rw2_pos = 0;
    for (i = rw_pos-rl; i < rw_pos; i++){
        sigmaL[rw2_pos++] = rw[i];
    }
    if ( ( cimag(sigmaL[0])>eps) && (fabs(creal(sigmaL[0])-creal(sigmaL[1]))> 1e-13)){
        rl++;
        rs--;
        rw2_pos = 0;
        for (i = rw_pos-rl; i < rw_pos; i++){
            sigmaL[rw2_pos++] = rw[i];
        }
    }

    /*-----------------------------------------------------------------------------
     *  small eigen values
     *-----------------------------------------------------------------------------*/
    rw2_pos = 0;
    for (i = 0; i < rs; i++){
        sigmaS[rw2_pos++] = rw[i];
    }
    if ( ( cimag(sigmaS[rs-1])>eps) && (fabs(creal(sigmaS[rs-1])-creal(sigmaS[rs-2]))> 1e-13)){
        rs++;
        rw2_pos = 0;
        for (i = 0; i < rs; i++){
            sigmaS[rw2_pos++] = rw[i];
        }
    }
    *r0 = rs+rl;
    ret = mess_vector_resize(sigma,*r0);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_resize);
    ret = mess_vector_tocomplex(sigma);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);

    rw2_pos = 0;
    for (i=0;i<rl;i++){
        sigma->values_cpx[rw2_pos++]=-sigmaL[i];
    }
    for (i=0;i<rs;i++){
        sigma->values_cpx[rw2_pos++]=-sigmaS[i];
    }

    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
    if ( have_mass ==1 ){
        mess_direct_clear(&mass_solver);
    }
    mess_free(rw);
    mess_free(rw2);
    mess_free(sigmaS);
    mess_free(sigmaL);

    return 0;
}       /* -----  end of function mess_h2_irka_init  ----- */

