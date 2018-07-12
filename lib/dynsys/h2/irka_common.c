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

/** @addtogroup dynsys_mor_h2
 * @{
 */

/**
 * @file lib/dynsys/h2/irka_common.c
 * @brief Common functions for IRKA.
 * @author @koehlerm
 */

#ifdef _OPENMP
#include <omp.h>
#endif
#include <string.h>

// define to control the behavior
#define AR_W_VA 1   // A=W'(AV) else A=(W'A)V
// #define ITSOLVER 1
/**
 * @brief Project \f$ A \f$ with \f$ W \f$ and \f$ V \f$.
 * @param[in] W  input left projector matrix
 * @param[in] A  input system matrix
 * @param[in] V  input right projector matrix
 * @param[out] Ar output \f$ Ar=W^T A V \f$
 * @return always zero
 *
 * The @ref __project_A function projects \f$ A \f$, that means
 * \f[ Ar= W^T A V.\f]
 *
 */
__attribute__((unused))  static int __project_A ( mess_matrix W, mess_matrix A, mess_matrix V, mess_matrix Ar )
{
    MSG_FNAME(__func__);
    mess_matrix tmp;
    int ret =0 ;

    /*-----------------------------------------------------------------------------
     *   check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(W);
    mess_check_nullpointer(A);
    mess_check_nullpointer(V);
    mess_check_nullpointer(Ar);
    ret  = mess_matrix_init(&tmp);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

    /*-----------------------------------------------------------------------------
     *  project A
     *  tmp = A*V;
     *  Ar= W'*tmp;
     *-----------------------------------------------------------------------------*/
#ifdef AR_W_VA
    ret = mess_matrix_multiply(MESS_OP_NONE, A, MESS_OP_NONE, V, tmp);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mul);
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, W, MESS_OP_NONE, tmp, Ar);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mul);
#else
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, W, MESS_OP_NONE, A, tmp);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mul);
    ret = mess_matrix_multiply(MESS_OP_NONE, tmp, MESS_OP_NONE, V, Ar);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_mul);
#endif

    mess_matrix_clear(&tmp);
    return 0;
}       /* -----  end of function __project_A  ----- */


/**
 * @brief Project \f$ b \f$.
 * @param[in] W  input left projector matrix
 * @param[in] b  input vector
 * @param[out] br output \f$ br = W^T B \f$
 * @return always zero
 * The @ref __project_B function projects \f$ b \f$, that means
 * \f[ br=W^T b. \f]
 *
 */
__attribute__((unused)) static int  __project_B ( mess_matrix W, mess_vector b, mess_vector br )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(W);
    mess_check_nullpointer(b);
    mess_check_nullpointer(br);

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,W,b,br);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_mvpt);


    return 0;
}       /* -----  end of function __project_B  ----- */

/**
 * @brief Project \f$ B \f$ (matrix version).
 * @param[in] W  input left projector matrix
 * @param[in] B  input matrix
 * @param[out] Br output \f$ Br = W^T B \f$
 * @return always zero
 * The @ref __project_B function projects \f$ B \f$, that means
 * \f[ Br=W^T B. \f]
 *
 */
__attribute__((unused)) static int  __project_Bmat ( mess_matrix W, mess_matrix B, mess_matrix Br )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(W);
    mess_check_nullpointer(B);
    mess_check_nullpointer(Br);

    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, W, MESS_OP_NONE, B, Br);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_mul);


    return 0;
}       /* -----  end of function __project_Bmat  ----- */


/**
 * @brief Project \f$ c \f$.
 * @param[in] V  input left projector matrix
 * @param[in] c  input vector
 * @param[out] cr output \f$ cr = V^T*c \f$
 * @return always zero
 * The @ref __project_C function projects  \f$ c \f$, that means
 * \f[ cr=V^T C. \f]
 *
 */
__attribute__((unused))  static int  __project_C ( mess_matrix V, mess_vector c, mess_vector cr )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(V);
    mess_check_nullpointer(c);
    mess_check_nullpointer(cr);

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,V,c,cr);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_mvpt);


    return 0;
}       /* -----  end of function __project_C  ----- */

/**
 * @brief Project \f$ C \f$ (matrix version).
 * @param[in] V  input right projector matrix
 * @param[in] C  input matrix
 * @param[in] Cr  input output \f$ Cr = CV \f$
 * @return always zero
 * The @ref __project_Cmat function projects  \f$ C \f$, that means
 * \f[ Cr=CV .\f]
 *
 */
__attribute__((unused))  static int  __project_Cmat ( mess_matrix V, mess_matrix C, mess_matrix Cr )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(V);
    mess_check_nullpointer(C);
    mess_check_nullpointer(Cr);

    ret = mess_matrix_multiply(MESS_OP_NONE, C, MESS_OP_NONE, V, Cr);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_mul);


    return 0;
}       /* -----  end of function __project_Cmat  ----- */



/**
 * @brief Construct \f$ V \f$ and \f$ W \f$ for IRKA (sequential).
 * @param[in] A  input system matrix
 * @param[in] Abase  input direct solver object
 * @param[in] E  input mass matrix
 * @param[in] b input vector
 * @param[out] c output vector
 * @param[in] sigma  input set of interpolation points
 * @param[out] V    output span
 * @param[out] W    output span
 * @return always zero
 * The @ref __constructVW function constructs the span \f$ V \f$ and \f$ W \f$ for IRKA.
 *
 */
static int __attribute__((unused)) __constructVW_seq ( mess_matrix A, mess_direct Abase, mess_matrix E, mess_vector b, mess_vector c, mess_vector sigma, mess_matrix V, mess_matrix W)
{
    MSG_FNAME(__func__);
    mess_int_t i, r,n;
    int ret = 0;
    double eps = mess_eps();
    mess_multidirect Asolver = NULL;
    mess_vector bc, cc;
    // mess_vector Vt,Wt,Vtc,Wtc;
    mess_matrix Atmp;
    // int * done;
    // int isdone;
    // int ignore = 0;
    mess_int_t *msolverid,j;
    unsigned short *datatypes;
    mess_int_t new_r;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(b);
    mess_check_nullpointer(c);
    mess_check_nullpointer(sigma);
    mess_check_nullpointer(V);
    mess_check_nullpointer(W);
    n = A->rows;
    r = sigma->dim;
    i = 0;

    // MSG_WARN("VW start\n");
    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_alloc(V,n,r,n*r,MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(W,n,r,n*r,MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    MESS_INIT_VECTORS(&bc,&cc);
    ret = mess_vector_alloc(bc, n, MESS_COMPLEX);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(cc, n, MESS_COMPLEX);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);



    /*-----------------------------------------------------------------------------
     *  shifts anpassen
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(msolverid, mess_int_t*, r*sizeof(mess_int_t));
    mess_try_alloc(datatypes, unsigned short*, sizeof(unsigned short)*r);
    for (j = 0; j<r;j++) msolverid[j] = j;
    j=0;
    for (i = 0; i < r; i++){
        if ( fabs(cimag(sigma->values_cpx[i])) < eps){
            msolverid[j]=i;
            datatypes[j] = MESS_REAL;
            j++;
        }else{
            msolverid[j] = i;
            datatypes[j]= MESS_COMPLEX;
            i++;
            j++;
        }
    }
    new_r = j ;

    ret = mess_multidirect_init(&Asolver);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_multidirect_init);
    ret = mess_matrix_init(&Atmp);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_copy(A,Atmp);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_scale(-1.0,Atmp);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_scale);

    ret = mess_multidirect_create(Atmp,NULL, sigma,Asolver,Abase,E);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_multidirect_create);

    mess_matrix_clear(&Atmp);

    ret = mess_vector_copy_tocomplex(b,bc);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy_tocomplex);
    ret = mess_vector_copy_tocomplex(c,cc);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy_tocomplex);

    for (i=0;i < new_r; i++){
        mess_vector Vt=NULL,Wt=NULL,Vtc=NULL,Wtc=NULL;

        if (datatypes[i] == MESS_COMPLEX){
            // MSG_INFO("shift " MESS_PRINTF_INT " is complex\n", i);
            MESS_INIT_VECTORS(&Vt,&Wt,&Vtc,&Wtc);
            mess_vector_alloc(Vt,n,MESS_REAL);
            mess_vector_alloc(Wt,n,MESS_REAL);
            mess_vector_alloc(Vtc,n,MESS_COMPLEX);
            mess_vector_alloc(Wtc,n,MESS_COMPLEX);

            mess_multidirect_solve(MESS_OP_NONE,Asolver,msolverid[i],bc,Vtc);
            mess_vector_realpart(Vtc,Vt);
            mess_matrix_setcol(V,msolverid[i],Vt);
            mess_vector_imagpart(Vtc,Vt);
            mess_matrix_setcol(V,msolverid[i]+1,Vt);

            mess_multidirect_solve(MESS_OP_HERMITIAN,Asolver,msolverid[i],cc,Wtc);
            mess_vector_realpart(Wtc,Wt);
            mess_matrix_setcol(W,msolverid[i],Wt);
            mess_vector_imagpart(Wtc,Wt);
            mess_matrix_setcol(W,msolverid[i]+1,Wt);
        } else {
            MESS_INIT_VECTORS(&Vt,&Wt);
            mess_vector_alloc(Vt,n,MESS_REAL);
            mess_vector_alloc(Wt,n,MESS_REAL);

            mess_multidirect_solve(MESS_OP_NONE,Asolver, msolverid[i], b,Vt);
            mess_matrix_setcol(V,msolverid[i],Vt);

            mess_multidirect_solve(MESS_OP_HERMITIAN,Asolver,msolverid[i], c, Wt);
            mess_matrix_setcol(W,msolverid[i],Wt);
        }
        mess_vector_clear(&Vt);
        mess_vector_clear(&Vtc);
        mess_vector_clear(&Wt);
        mess_vector_clear(&Wtc);
    }


    // MSG_WARN("VW complete\n");
    mess_multidirect_clear(&Asolver);
    mess_free(msolverid);
    mess_free(datatypes);
    mess_vector_clear(&bc);
    mess_vector_clear(&cc);
    return 0;
}       /* -----  end of function __constructVW  ----- */




/**
 * @brief Construct \f$ V \f$ and \f$ W \f$ for IRKA (OpenMP version).
 * @param[in] A         input system matrix
 * @param[in] Abase     input direct solver object
 * @param[in] E         input mass matrix
 * @param[in] b         input vector
 * @param[out] c        output vector
 * @param[in] sigma     input set of interpolation points
 * @param[out] V        output span
 * @param[out] W        output span
 * @return always zero
 * The @ref __constructVW_sec function constructs the span \f$ V \f$ and \f$ W \f$ for IRKA (OpenMP version).
 *
 */
static int __attribute__((unused))  __constructVW_sec ( mess_matrix A, mess_direct Abase, mess_matrix E, mess_vector b, mess_vector c, mess_vector sigma, mess_matrix V, mess_matrix W)
{
    MSG_FNAME(__func__);
    mess_int_t i, r,n;
    int ret = 0;
    mess_multidirect Asolver = NULL;
    mess_vector Vt,Wt,Vtc,Wtc;
    mess_matrix Atmp;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(b);
    mess_check_nullpointer(c);
    mess_check_nullpointer(sigma);
    mess_check_nullpointer(V);
    mess_check_nullpointer(W);
    n = A->rows;
    r = sigma->dim;
    i = 0;

    // MSG_WARN("VW start\n");
    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_alloc(V,n,r,n*r,MESS_DENSE, MESS_REAL);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(W,n,r,n*r,MESS_DENSE, MESS_REAL);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    ret = mess_multidirect_init(&Asolver);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_multidirect_init);
    ret=mess_matrix_init(&Atmp);                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_copy(A,Atmp);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_scale(-1.0,Atmp);                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_scale);

    ret = mess_multidirect_create(Atmp,NULL, sigma,Asolver,Abase,E);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_multidirect_create);
    mess_matrix_clear(&Atmp);

    ret = mess_vector_tocomplex(b);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
    ret = mess_vector_tocomplex(c);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
    {
        MESS_INIT_VECTORS(&Vt,&Wt,&Vtc,&Wtc);
        ret = mess_vector_alloc(Vt,n,MESS_REAL);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_alloc(Wt,n,MESS_REAL);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_alloc(Vtc,n,MESS_COMPLEX);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_alloc(Wtc,n,MESS_COMPLEX);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);

        for (i=0;i < r; i++){
            if (cimag(sigma->values_cpx[i]) != 0.0){
                MSG_INFO("shift " MESS_PRINTF_INT " is complex\n", i);
                // complex
#ifdef _OPENMP
#pragma omp parallel sections
                {
#pragma omp section
                    {
#endif
                        mess_multidirect_solve(MESS_OP_NONE,Asolver,i,b,Vtc);
                        mess_vector_realpart(Vtc,Vt);
                        mess_matrix_setcol(V,i,Vtc);
                        mess_vector_imagpart(Vtc,Vt);
                        mess_matrix_setcol(V,i+1,Vt);
#ifdef _OPENMP
                    }
#pragma omp section
                    {
#endif

                        mess_multidirect_solve(MESS_OP_HERMITIAN,Asolver,i,c,Wtc);
                        mess_vector_realpart(Wtc,Wt);
                        mess_matrix_setcol(W,i,Wt);
                        mess_vector_imagpart(Wtc,Wt);
                        mess_matrix_setcol(W,i+1,Wt);
#ifdef _OPENMP
                    }
                }
#endif
                i++;


            } else {
                MSG_INFO("shift " MESS_PRINTF_INT " is real\n", i);
                // real
#ifdef _OPENMP
#pragma omp parallel sections
                {
#pragma omp section
                    {
#endif

                        mess_multidirect_solve(MESS_OP_NONE,Asolver, i, b,Vt);
                        mess_matrix_setcol(V,i,Vt);
#ifdef _OPENMP
                    }
#pragma omp section
                    {
#endif

                        mess_multidirect_solve(MESS_OP_HERMITIAN,Asolver,i, c, Wt);
                        mess_matrix_setcol(W,i,Wt);
#ifdef _OPENMP
                    }
                }
#endif

                // i++;
            }

        }
        mess_vector_clear(&Vt);
        mess_vector_clear(&Wt);
        mess_vector_clear(&Vtc);
        mess_vector_clear(&Wtc);

    }
    mess_multidirect_clear(&Asolver);
    return 0;
}       /* -----  end of function __constructVW  ----- */

/* Structure for Pthread arguments   */
struct VWjob {
    mess_double_cpx_t sigma;
    mess_int_t i;
    mess_multidirect Asolver;
    mess_vector b;
    mess_vector bc;
    mess_vector c;
    mess_vector cc;
    mess_matrix V;
    mess_matrix W;
    mess_int_t msolverid;
    mess_datatype_t datatype;
    mess_int_t n;
};

/* PThread Worker  */
static void VW_worker(void * data){
    MSG_FNAME(__func__);
    struct VWjob * job = (struct VWjob *) data;
    mess_vector Vt=NULL,Wt=NULL,Vtc=NULL,Wtc=NULL;

    if (job->datatype == MESS_COMPLEX){
        MSG_INFO("shift " MESS_PRINTF_INT " is complex\n", job->i);
        MESS_INIT_VECTORS(&Vt,&Wt,&Vtc,&Wtc);
        mess_vector_alloc(Vt,job->n,MESS_REAL);
        mess_vector_alloc(Wt,job->n,MESS_REAL);
        mess_vector_alloc(Vtc,job->n,MESS_COMPLEX);
        mess_vector_alloc(Wtc,job->n,MESS_COMPLEX);

        mess_multidirect_solve(MESS_OP_NONE, job->Asolver,job->msolverid,job->bc,Vtc);
        mess_vector_realpart(Vtc,Vt);
        mess_matrix_setcol(job->V,job->i,Vt);
        mess_vector_imagpart(Vtc,Vt);
        mess_matrix_setcol(job->V,job->i+1,Vt);

        mess_multidirect_solve(MESS_OP_HERMITIAN,job->Asolver,job->msolverid,job->cc,Wtc);
        mess_vector_realpart(Wtc,Wt);
        mess_matrix_setcol(job->W,job->i,Wt);
        mess_vector_imagpart(Wtc,Wt);
        mess_matrix_setcol(job->W,job->i+1,Wt);
    } else {
        MSG_INFO("shift " MESS_PRINTF_INT " is real\n", job->i);
        MESS_INIT_VECTORS(&Vt,&Wt);
        mess_vector_alloc(Vt,job->n,MESS_REAL);
        mess_vector_alloc(Wt,job->n,MESS_REAL);

        mess_multidirect_solve(MESS_OP_NONE,job->Asolver, job->msolverid, job->b,Vt);
        mess_matrix_setcol(job->V,job->i,Vt);

        mess_multidirect_solve(MESS_OP_HERMITIAN,job->Asolver,job->msolverid, job->c, Wt);
        mess_matrix_setcol(job->W,job->i,Wt);
    }

    mess_vector_clear(&Vt);
    mess_vector_clear(&Wt);
    mess_vector_clear(&Vtc);
    mess_vector_clear(&Wtc);
    mess_free(job);
}

/**
 * @brief Construct \f$ V \f$ and \f$ W \f$ for IRKA (OpenMP Preprocessing).
 * @param[in] A  input system matrix
 * @param[in] Abase  input direct solver object
 * @param[in] E  input mass matrix
 * @param[in] b input vector
 * @param[out] c output vector
 * @param[in] sigma  input set of interpolation points
 * @param[out] V    output span
 * @param[out] W    output span
 *
 * The @ref __constructVW function constructs the span \f$ V \f$ and \f$ W \f$ for IRKA (OpenMP Preprocessing).
 *
 */
static int __attribute__((unused)) __constructVW_pre( mess_matrix A, mess_direct Abase, mess_matrix E, mess_vector b, mess_vector c, mess_vector sigma, mess_matrix V, mess_matrix W)
{
        MSG_FNAME(__func__);
    mess_int_t i, r,n;
    int ret = 0;
    double eps = mess_eps();
    mess_multidirect Asolver = NULL;
    mess_vector bc, cc;
    mess_matrix Atmp;
    mess_int_t *msolverid,j;
    unsigned short *datatypes;
    mess_int_t new_r;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(b);
    mess_check_nullpointer(c);
    mess_check_nullpointer(sigma);
    mess_check_nullpointer(V);
    mess_check_nullpointer(W);
    n = A->rows;
    r = sigma->dim;
    i = 0;

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_alloc(V,n,r,n*r,MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(W,n,r,n*r,MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    MESS_INIT_VECTORS(&bc,&cc);
    ret = mess_vector_alloc(bc, n, MESS_COMPLEX);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(cc, n, MESS_COMPLEX);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);



    /*-----------------------------------------------------------------------------
     *  shifts anpassen
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(msolverid, mess_int_t*, r*sizeof(mess_int_t));
    mess_try_alloc(datatypes, unsigned short*, sizeof(unsigned short)*r);
    for (j = 0; j<r;j++) msolverid[j] = j;
    j=0;
    for (i = 0; i < r; i++){
        if ( fabs(cimag(sigma->values_cpx[i])) < eps){
            msolverid[j]=i;
            datatypes[j] = MESS_REAL;
            j++;
        }else{
            msolverid[j] = i;
            datatypes[j]= MESS_COMPLEX;
            i++;
            j++;
        }
    }
    new_r = j ;

    ret = mess_multidirect_init(&Asolver);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_multidirect_init);
    ret = mess_matrix_init(&Atmp);                                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_copy(A,Atmp);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_scale(-1.0,Atmp);                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_scale);

    ret = mess_multidirect_create(Atmp,NULL, sigma,Asolver,Abase,E);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_multidirect_create);
    mess_matrix_clear(&Atmp);

    ret = mess_vector_copy_tocomplex(b,bc);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy_tocomplex);
    ret = mess_vector_copy_tocomplex(c,cc);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy_tocomplex);

#ifdef _OPENMP
#pragma omp parallel for private(i) default(shared) schedule(static,1)
#endif
        for (i=0;i < new_r; i++){
            mess_vector Vt=NULL,Wt=NULL,Vtc=NULL,Wtc=NULL;

            if (datatypes[i] == MESS_COMPLEX){
                MSG_INFO("shift " MESS_PRINTF_INT " is complex\n", i);
                MESS_INIT_VECTORS(&Vt,&Wt,&Vtc,&Wtc);
                mess_vector_alloc(Vt,n,MESS_REAL);
                mess_vector_alloc(Wt,n,MESS_REAL);
                mess_vector_alloc(Vtc,n,MESS_COMPLEX);
                mess_vector_alloc(Wtc,n,MESS_COMPLEX);

                mess_multidirect_solve(MESS_OP_NONE,Asolver,msolverid[i],bc,Vtc);
                mess_vector_realpart(Vtc,Vt);
                mess_matrix_setcol(V,msolverid[i],Vt);
                mess_vector_imagpart(Vtc,Vt);
                mess_matrix_setcol(V,msolverid[i]+1,Vt);

                mess_multidirect_solve(MESS_OP_HERMITIAN,Asolver,msolverid[i],cc,Wtc);
                mess_vector_realpart(Wtc,Wt);
                mess_matrix_setcol(W,msolverid[i],Wt);
                mess_vector_imagpart(Wtc,Wt);
                mess_matrix_setcol(W,msolverid[i]+1,Wt);
            } else {
                MSG_INFO("shift " MESS_PRINTF_INT " is real\n", i);
                MESS_INIT_VECTORS(&Vt,&Wt);
                mess_vector_alloc(Vt,n,MESS_REAL);
                mess_vector_alloc(Wt,n,MESS_REAL);

                mess_multidirect_solve(MESS_OP_NONE,Asolver, msolverid[i], b,Vt);
                mess_matrix_setcol(V,msolverid[i],Vt);

                mess_multidirect_solve(MESS_OP_HERMITIAN,Asolver,msolverid[i], c, Wt);
                mess_matrix_setcol(W,msolverid[i],Wt);
            }
            mess_vector_clear(&Vt);
            mess_vector_clear(&Vtc);
            mess_vector_clear(&Wt);
            mess_vector_clear(&Wtc);
        }


    mess_multidirect_clear(&Asolver);
    mess_free(msolverid);
    mess_free(datatypes);
    mess_vector_clear(&bc);
    mess_vector_clear(&cc);
    return 0;
}       /* -----  end of function __constructVW  ----- */


/**
 * @brief Construct \f$ V \f$ and \f$ W \f$ for IRKA (OpenMP Preprocessing).
 * @param[in] A  input system matrix
 * @param[in] Abase  input direct solver object
 * @param[in] E  input mass matrix
 * @param[in] b input vector
 * @param[out] c output vector
 * @param[in] sigma  input set of interpolation points
 * @param[out] V    output span
 * @param[out] W    output span
 *
 * The @ref __constructVW function constructs the span \f$ V \f$ and \f$ W \f$ for IRKA (PThread enabled version).
 *
 */
static int __attribute__((unused)) __constructVW_pthread( mess_matrix A, mess_direct Abase, mess_matrix E, mess_vector b, mess_vector c, mess_vector sigma, mess_matrix V, mess_matrix W)
{
        MSG_FNAME(__func__);
    mess_int_t i, r,n;
    int ret = 0;
    double eps = mess_eps();
    mess_multidirect Asolver = NULL;
    mess_vector bc, cc;
    mess_matrix Atmp;
    mess_vector sigma_tmp;
    mess_threadpool pool;
    mess_int_t *ids;
    mess_int_t jobcount=0;
    mess_int_t *msolverid,j;
    unsigned short *datatypes;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(b);
    mess_check_nullpointer(c);
    mess_check_nullpointer(sigma);
    mess_check_nullpointer(V);
    mess_check_nullpointer(W);
    n = A->rows;
    r = sigma->dim;
    i = 0;

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_alloc(V,n,r,n*r,MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(W,n,r,n*r,MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    MESS_INIT_VECTORS(&sigma_tmp,&bc,&cc);
    ret = mess_vector_alloc(sigma_tmp, sigma->dim, MESS_COMPLEX);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_vector_alloc);
    ret = mess_vector_alloc(bc, n, MESS_COMPLEX);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(cc, n, MESS_COMPLEX);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);



    /*-----------------------------------------------------------------------------
     *  shifts anpassen
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(msolverid, mess_int_t*, r*sizeof(mess_int_t));
    mess_try_alloc(datatypes, unsigned short*, sizeof(unsigned short)*r);
    for (j = 0; j<r;j++) msolverid[j] = j;
    j=0;
    for (i = 0; i < r; i++){
        if ( fabs(cimag(sigma->values_cpx[i])) < eps){
            sigma_tmp->values_cpx[j] = sigma->values_cpx[i];
            // msolverid[i] = j;
            msolverid[i]=i;
            datatypes[i] = MESS_REAL;
            j++;
        }else{
            sigma_tmp->values_cpx[j] = sigma->values_cpx[i];
            msolverid[i] = i;
            msolverid[i+1] = i+1;
            datatypes[i]= MESS_COMPLEX;
            datatypes[i+1]=MESS_COMPLEX;
            i++;
            j++;
        }
    }
    ret = mess_vector_resize(sigma_tmp,j);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vetor_resize);
    // mess_vector_print(sigma_tmp);


    ret = mess_multidirect_init(&Asolver);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_multidirect_init);
    ret = mess_matrix_init(&Atmp);                                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
    ret = mess_matrix_copy(A,Atmp);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_scale(-1.0,Atmp);                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_scale);
    ret = mess_multidirect_create(Atmp,NULL, sigma,Asolver,Abase,E);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_multidirect_create);

    mess_matrix_clear(&Atmp);

    ret = mess_vector_copy_tocomplex(b,bc);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy_tocomplex);
    ret = mess_vector_copy_tocomplex(c,cc);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy_tocomplex);

#ifdef _OPENMP
        int th = omp_get_max_threads();
#else
        int th = 2;
#endif

    ret = mess_threadpool_init(&pool, th , 256);

    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_threadpool_init);
    jobcount = 0;
    mess_try_alloc(ids,  mess_int_t *, sizeof(mess_int_t)*r);

    for (i=0;i < r; i++){
        struct VWjob *vwjob;
        mess_threadpooljob job;

        mess_threadpooljob_init(&job);
        job->worker = VW_worker;
        mess_try_alloc(vwjob, struct VWjob*, sizeof(struct VWjob));

        vwjob->V=V;
        vwjob->W=W;
        vwjob->Asolver = Asolver;
        vwjob->b=b;
        vwjob->c=c;
        vwjob->bc=bc;
        vwjob->cc=cc;
        vwjob->n=n;
        vwjob->i=i;
        vwjob->msolverid = msolverid[i];
        vwjob->datatype=datatypes[i];
        vwjob->sigma=sigma->values_cpx[i];
        job->data = vwjob;

        mess_threadpool_insert(pool, job,&ids[jobcount++]);

        if (cimag(sigma->values_cpx[i]) != 0.0){
            i++;
        }
    }


    for (i=0; i<jobcount;i++){
        MSG_INFO("wait until " MESS_PRINTF_INT " is done\n", ids[i]);
        mess_threadpool_waitdone(pool, ids[i]);

    }

    mess_threadpool_clear(&pool);
    mess_multidirect_clear(&Asolver);
    mess_vector_clear(&sigma_tmp);
    mess_free(msolverid);
    mess_free(datatypes);
    mess_free(ids);
    mess_vector_clear(&bc);
    mess_vector_clear(&cc);
    return 0;
}       /* -----  end of function __constructVW  ----- */



static char * get_char(char *c1, char* c2)
{
    return "pthread";
}

/**
 * @brief Construct \f$ V \f$ and \f$ W \f$ for IRKA (wrapper).
 * @param[in] A         input system matrix
 * @param[in] Abase     input direct solver object
 * @param[in] E         input mass matrix
 * @param[in] b         input vector
 * @param[out] c        output vector
 * @param[in] sigma     input set of interpolation points
 * @param[out] V        output span
 * @param[out] W        output span
 *
 * The @ref __constructVW function contructs \f$ V \f$ and \f$ W \f$ for IRKA (wrapper).
 *
 */
static int __attribute__((unused))  __constructVW ( mess_matrix A, mess_direct Abase, mess_matrix E, mess_vector b, mess_vector c, mess_vector sigma, mess_matrix V, mess_matrix W)
{
    MSG_FNAME(__func__);

    if ( strcmp(get_char("mess.irka.constructVW","pthread"),"sequential") == 0 ) {
        MSG_INFO("use sequential V,W construction\n");
        return __constructVW_seq(A,Abase,E,b,c,sigma,V,W);
    } else if ( strcmp ( get_char("mess.irka.constructVW","pthread"), "omp_section") ==  0 ) {
        MSG_INFO("use OpenMP Section based V,W construction\n");
        return __constructVW_sec(A,Abase,E,b,c,sigma,V,W);
    } else if ( strcmp( get_char("mess.irka.constructVW", "pthread"), "omp_pre") == 0) {
        MSG_INFO("use OpenMP Preprocessing based V,W construction\n");
        return __constructVW_sec(A,Abase,E,b,c,sigma,V,W);
    } else {
        MSG_INFO("use PThreads based V,W construction\n");
        return __constructVW_pthread(A,Abase,E,b,c,sigma,V,W);
    }
    return 0;
}

/**
 * @brief Construct \f$ V \f$ and \f$ W \f$ for IRKA (wrapper, matrix inputs).
 * @param[in] A         input system matrix
 * @param[in] Abase     input direct solver object
 * @param[in] E         input mass matrix
 * @param[in] B         input matrix
 * @param[out] C        output matrix
 * @param[in] sigma     input set of interpolation points
 * @param[out] V        output span
 * @param[out] W        output span

 * The @ref __constructVWmat function contructs \f$ V \f$ and \f$ W \f$ for IRKA.
 *
 */
static int __attribute__((unused))  __constructVWmat ( mess_matrix A, mess_direct Abase, mess_matrix E, mess_matrix B, mess_matrix C, mess_vector sigma, mess_matrix V, mess_matrix W)
{
    MSG_FNAME(__func__);
    mess_vector b, c;
    mess_int_t n;
    int ret;

    n=A->rows;
    MESS_INIT_VECTORS(&b,&c);
    ret = mess_vector_alloc(b, n, MESS_REAL);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(c, n, MESS_REAL);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_alloc);
    ret = mess_matrix_getcol(B,0,b);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_getcol);
    ret = mess_matrix_getrow(C,0,c);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_getrow);

    if ( strcmp(get_char("mess.irka.constructVW","pthread"),"sequential") == 0 ) {
        MSG_INFO("use sequential V,W construction\n");
        ret =  __constructVW_seq(A,Abase,E,b,c,sigma,V,W);
    } else if ( strcmp ( get_char("mess.irka.constructVW","pthread"), "omp_section") ==  0 ) {
        MSG_INFO("use OpenMP Section based V,W construction\n");
        ret = __constructVW_sec(A,Abase,E,b,c,sigma,V,W);
    } else if ( strcmp( get_char("mess.irka.constructVW", "pthread"), "omp_pre") == 0) {
        MSG_INFO("use OpenMP Preprocessing based V,W construction\n");
        ret=  __constructVW_sec(A,Abase,E,b,c,sigma,V,W);
    } else {
        MSG_INFO("use PThreads based V,W construction\n");
        ret= __constructVW_pthread(A,Abase,E,b,c,sigma,V,W);
    }

    mess_vector_clear(&b);
    mess_vector_clear(&c);
    return ret;
}


/**
 * \}@
 */

