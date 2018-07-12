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
 * @file lib/eigen/arnoldi.c
 * @brief Arnoldi iteration.
 * @author @koehlerm
 *
 *   This file contains different variants of the Arnoldi process.
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


/*-----------------------------------------------------------------------------
 *  internal stuff for @ref mess_eigen_arnoldi
 *-----------------------------------------------------------------------------*/
struct mvpdata {
    mess_matrix A;
    mess_matrix U;
    mess_matrix W;
    mess_vector h1;
    mess_vector h2;
    int haveUW;
};

static int mvp_arnoldi(void *data, mess_operation_t op, mess_vector in, mess_vector out){
    struct mvpdata * d = (struct mvpdata*) data;
    if ( d->haveUW == 0){
        return mess_matrix_mvp(MESS_OP_NONE,d->A, in, out);
    } else {
        mess_matrix_mvp(MESS_OP_NONE,d->A, in, out);
        mess_matrix_mvp(MESS_OP_NONE,d->W, in, d->h1);
        mess_matrix_mvp(MESS_OP_NONE,d->U, d->h1, d->h2);
        mess_vector_axpy(-1, d->h2, out);
        return 0;
    }
}

/**
 * @brief Compute matrices \f$ H \f$ and \f$ V \f$ with Arnoldis method with respect to  \f$ A  \f$ or \f$ (A-UW)  \f$.
 * @param[in] A  input matrix for which the Arnoldi algorithm is supposed to run
 * @param[in] U  input matrix to work with  \f$ (A-UW)  \f$
 * @param[in] W  input matrix to work with  \f$ (A-UW)  \f$
 * @param[in] k  input number of Arnoldi steps (usually  \f$ k \ll n  \f$ )
 * @param[in] sv   input start vector
 * @param[out] H  \f$(k+1 \times k) \f$ upper Hessenberg matrix
 * @param[out] V  \f$ (n \times k+1)  \f$ orthogonal matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_arnoldi function computes matrices \f$ H  \f$ and  \f$ V  \f$ (if  \f$ V  \ne \f$ @c NULL) fullfilling
 * \f[
 * \begin{array}{ccc}
 *      V(:,1) &=& \frac{sv}{ \Vert sv \Vert }, \\
 *      V^T V  &=& I, \\
 *   AV(:,1:k) &=& VH .
 * \end{array}
 * \f]
 *
 */
int mess_eigen_arnoldi(mess_matrix A, mess_matrix U, mess_matrix W,mess_int_t k, mess_vector sv, mess_matrix H, mess_matrix V)
{
    MSG_FNAME(__func__);
    int ret = 0;
    struct mvpdata mvpdat;
    mess_mvpcall call;

    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_square(A);
    mess_check_real(A);
    mess_check_nullpointer(sv);
    mess_check_nullpointer(H);
    if ( k == 0) return 0;
    mess_check_positive(k);

    mvpdat.A = A;
    mvpdat.haveUW = 0;
    if ( U != NULL && W != NULL ) {
        if (!( U->rows == 0 || U->cols == 0 || W->rows == 0 || W->cols == 0)){
            MSG_INFO("Arnoldi process w.r.t A-UW\n");
            mvpdat.haveUW = 1;
            mvpdat.U = U;
            mvpdat.W = W;
            MESS_INIT_VECTORS(&(mvpdat.h1),&(mvpdat.h2));
            ret = mess_vector_alloc(mvpdat.h1, W->rows,MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
            ret = mess_vector_alloc(mvpdat.h2, U->rows,MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_allocc);
        }
    }
    ret = mess_mvpcall_operator(&call, A->rows, A->data_type, mvp_arnoldi, &mvpdat); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator );

    ret = mess_eigen_arnoldi_template(call, k, sv, H, V); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template);
    if ( mvpdat.haveUW == 1){
        mess_vector_clear(&(mvpdat.h1));
        mess_vector_clear(&(mvpdat.h2));
    }
    mess_mvpcall_clear(&call);
    return 0;
}


/*-----------------------------------------------------------------------------
 *  internal stuff for @ref mess_eigen_arnoldig
 *-----------------------------------------------------------------------------*/
struct mvpdatag {
    mess_matrix A;
    mess_matrix E;
    mess_matrix U;
    mess_matrix W;
    mess_vector h1;
    mess_vector h2;
    mess_vector y;
    mess_direct Esolver;
    mess_direct Asolver;
    mess_direct JaySolver;
    mess_matrix Utilde;
    mess_vector z,u,v, vj;
    int haveUW;
    int Atilde;
};

static int mvp_arnoldig(void *data, mess_operation_t op, mess_vector in, mess_vector out){
    struct mvpdatag * d = (struct mvpdatag*) data;
    if ( d->haveUW == 0){
        mess_matrix_mvp(MESS_OP_NONE,d->A, in, d->y);
    } else {
        mess_matrix_mvp(MESS_OP_NONE,d->A, in, d->y);
        mess_matrix_mvp(MESS_OP_NONE,d->W, in, d->h1);
        mess_matrix_mvp(MESS_OP_NONE,d->U, d->h1, d->h2);
        mess_vector_axpy(-1, d->h2, d->y);
    }
    mess_direct_solve(MESS_OP_NONE,d->Esolver, d->y, out);
    return 0;
}



/**
 * @brief Compute matrices \f$ H \f$ and \f$ V  \f$ with Arnoldis method with respect to  \f$E^{-1} A  \f$ or  \f$ E^{-1} (A-UW) \f$.
 * @param[in] A   input matrix for which the Arnoldi algorithm is supposed to run
 * @param[in] E  input solver to solve \f$Ex=b\f$
 * @param[in] U  input matrix to work with  \f$ (A-UW) \f$
 * @param[in] W  input matrix to work with  \f$(A-UW) \f$
 * @param[in] k  input number of Arnoldi steps (usually  \f$ k \ll n \f$ )
 * @param[in] sv   input start vector
 * @param[out] H  \f$(k+1 \times k) \f$ upper Hessenberg matrix
 * @param[out] V  \f$ (n \times k+1)  \f$ orthogonal matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_arnoldig function computes matrices \f$ H  \f$ and  \f$ V  \f$ (if  \f$ V  \ne \f$ @c NULL) fullfilling
 * \f[
 * \begin{array}{ccc}
 *      V(:,1) &=& \frac{sv}{ \Vert sv \Vert }, \\
 *      V^T V  &=& I, \\
 *   AV(:,1:k) &=& VH .
 * \end{array}
 * \f]
 *
 */
int mess_eigen_arnoldig(mess_matrix A, mess_matrix U, mess_matrix W, mess_direct E, mess_int_t k, mess_vector sv, mess_matrix H, mess_matrix V)
{
    MSG_FNAME(__func__);
    int ret = 0;
    struct mvpdatag mvpdat;
    mess_mvpcall call;

    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(E);
    mess_check_nullpointer(sv);
    mess_check_nullpointer(H);
    mess_check_real(sv);
    mess_check_real(A);
    // mess_check_real(E);
    if ( k == 0) return 0;
    mess_check_positive(k);

    mvpdat.A = A;
    mvpdat.Esolver = E;
    mvpdat.haveUW = 0;
    ret = mess_vector_init(&(mvpdat.y));                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(mvpdat.y, A->rows, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);

    if ( U != NULL && W != NULL ) {
        if (!( U->rows == 0 || U->cols == 0 || W->rows == 0 || W->cols == 0)){
            MSG_INFO("Arnoldi process w.r.t A-UW\n");
            mvpdat.haveUW = 1;
            mvpdat.U = U;
            mvpdat.W = W;
            MESS_INIT_VECTORS(&(mvpdat.h1),&(mvpdat.h2));
            ret = mess_vector_alloc((mvpdat.h1), W->rows,MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
            ret = mess_vector_alloc((mvpdat.h2), U->rows,MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        }
    }
    ret = mess_mvpcall_operator(&call, A->rows, A->data_type, mvp_arnoldig, &mvpdat); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);
    ret = mess_eigen_arnoldi_template(call, k, sv, H, V); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template);

    if ( mvpdat.haveUW == 1) {
        mess_vector_clear(&(mvpdat.h1));
        mess_vector_clear(&(mvpdat.h2));
    }
    mess_vector_clear(&(mvpdat.y));
    mess_mvpcall_clear(&call);
    return 0;
}


static int mvp_arnoldi_inv(void *data,mess_operation_t op,  mess_vector in, mess_vector out){
    struct mvpdatag * d = (struct mvpdatag*) data;
    if ( d->Atilde == 0){
        mess_direct_solve(MESS_OP_NONE,d->Asolver, in, out);
    } else {
        mess_direct_solve(MESS_OP_NONE,d->Asolver, in, d->z);
        mess_vector_copy(d->z,out);
        mess_matrix_mvp(MESS_OP_NONE,d->W, d->z, d->u);
        mess_direct_solve(MESS_OP_NONE,d->JaySolver, d->u,d->v);
        mess_matrix_mvp(MESS_OP_NONE,d->Utilde, d->v, d->z);
        mess_vector_axpy(1.0, d->z, out);
    }
    return 0;

}


/**
 * @brief Compute matrices \f$ H  \f$ and \f$ V  \f$ with Arnoldis method on \f$ A^{-1} \f$  or \f$  (A-UW)^{-1} \f$.
 * @param[in] A  input decomposed matrices for which the Arnoldi algorithm is supposed to run
 * @param[in] U  input matrix for \f$ (A-UW)^{-1} \f$
 * @param[in] W  input matrix for \f$ (A-UW)^{-1} \f$
 * @param[in] k  input number of Arnoldi steps (usually \f$  k \ll n \f$ )
 * @param[in] sv   input start vector
 * @param[out] H  \f$(k+1 \times k) \f$ upper Hessenberg matrix
 * @param[out] V  \f$ (n \times k+1)  \f$ orthogonal matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_arnoldi_inv function computes matrices \f$ H  \f$ and  \f$ V  \f$ (if  \f$ V  \ne \f$ @c NULL) fullfilling
 * \f[
 * \begin{array}{ccc}
 *      V(:,1) &=& \frac{sv}{ \Vert sv \Vert }, \\
 *      V^T V  &=& I, \\
 *   A^{-1}V(:,1:k) &=& VH .
 * \end{array}
 * \f]
 *
 */
int mess_eigen_arnoldi_inv(mess_direct A, mess_matrix U, mess_matrix W, mess_int_t k, mess_vector sv, mess_matrix H, mess_matrix V)
{
    MSG_FNAME(__func__);
    mess_int_t n,i, p = 0;
    int Atilde = 0;
    mess_matrix Utilde;
    mess_matrix Jay;
    mess_direct JaySolver;
    mess_vector z, u,v;
    int ret = 0;
    struct mvpdatag mvpdat;
    mess_mvpcall call;



    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(sv);
    mess_check_real(sv);
    mess_check_nullpointer(H);
    if ( k == 0) return 0;
    mess_check_positive(k);

    mvpdat.Atilde = 0;
    if ( U != NULL && W != NULL ){
        if (!( U->rows == 0 || U->cols == 0 || W->rows == 0 || W->cols == 0)){
            MSG_INFO("Arnoldi process w.r.t inv(A-UW)\n");
            Atilde = 1;
            mvpdat.U=U;
            mvpdat.W=W;
            mvpdat.Atilde = 1;
        }
    }
    n = A->rows;
    if ( k>n) {
        MSG_ERROR("k must be smaller than the order of A!\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  prepare the iteration
     *-----------------------------------------------------------------------------*/
    if ( Atilde == 1) {
        mess_matrix tmp;
        p = W->rows;
        MESS_INIT_MATRICES(&Utilde,&Jay,&tmp);
        MESS_INIT_VECTORS(&z,&u,&v);
        ret = mess_vector_alloc(z, n, MESS_REAL);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_alloc(u, p, MESS_REAL);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_alloc(v, p, MESS_REAL);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);

        // Btilde = A\B
        ret =mess_direct_solvem(MESS_OP_NONE,A, U, Utilde);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
        ret = mess_matrix_multiply(MESS_OP_NONE, W, MESS_OP_NONE, Utilde, tmp);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mul);

        ret = mess_matrix_alloc(Jay, p, p, p*p, MESS_DENSE, MESS_REAL);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);
        for ( i = 0; i  < p; i++){
            Jay->values[i+i*p] = 1.0;
        }
        mess_matrix_add(-1.0, tmp, 1, Jay);
        ret =mess_direct_init(&JaySolver);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
        ret = mess_direct_create_lapack_lu(Jay, JaySolver);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_lapack_lu);
        mess_matrix_clear(&tmp);
        mvpdat.JaySolver = JaySolver;
        mvpdat.z = z;
        mvpdat.v = v;
        mvpdat.u = u;
        mvpdat.Utilde = Utilde;
    }


    mvpdat.Asolver = A;
    ret = mess_mvpcall_operator(&call, A->rows, A->data_type, mvp_arnoldi_inv, &mvpdat); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);
    ret = mess_eigen_arnoldi_template(call, k, sv, H, V); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template);
    if ( Atilde == 1) {
        mess_vector_clear(&z);
        mess_vector_clear(&u);
        mess_vector_clear(&v);
        mess_matrix_clear(&Utilde);
        mess_matrix_clear(&Jay);
        mess_direct_clear(&JaySolver);
    }
    mess_mvpcall_clear(&call);
    return 0;
}


static int mvp_arnoldig_inv(void *data, mess_operation_t op, mess_vector in, mess_vector out){
    struct mvpdatag * d = (struct mvpdatag*) data;
    mess_matrix_mvp(MESS_OP_NONE,d->E, in, d->vj);
    if ( d->Atilde == 0){
        mess_direct_solve(MESS_OP_NONE,d->Asolver, d->vj, out);
    } else {
        mess_direct_solve(MESS_OP_NONE,d->Asolver, d->vj, d->z);
        mess_vector_copy(d->z,out);
        mess_matrix_mvp(MESS_OP_NONE,d->W, d->z, d->u);
        mess_direct_solve(MESS_OP_NONE,d->JaySolver, d->u,d->v);
        mess_matrix_mvp(MESS_OP_NONE,d->Utilde, d->v, d->z);
        mess_vector_axpy(1.0, d->z, out);
    }
    return 0;

}



/**
 * @brief Compute matrices \f$ H  \f$ and \f$ V  \f$ with Arnoldis method on \f$ A^{-1} E \f$ or \f$(A-UW)^{-1} E \f$.
 * @param[in] A  input decomposed matrices for which the Arnoldi algorithm is supposed to run
 * @param[in] E  input mass matrix
 * @param[in] U   input matrix for \f$(A-UW)^{-1} \f$
 * @param[in] W   input matrix for \f$(A-UW)^{-1} \f$
 * @param[in] k  input number of Arnoldi steps (usually \f$ k \ll n \f$);
 * @param[in] sv   input start vector
 * @param[out] H  \f$(k+1 \times k) \f$ upper Hessenberg matrix
 * @param[out] V  \f$ (n \times k+1)  \f$ orthogonal matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_arnoldig_inv function computes matrices \f$ H  \f$ and  \f$ V  \f$ (if  \f$ V  \ne \f$ @c NULL) fullfilling
 * \f[
 * \begin{array}{ccc}
 *      V(:,1) &=& \frac{sv}{ \Vert sv \Vert }, \\
 *      V^T V  &=& I, \\
 *   A^{-1}V(:,1:k) &=& VH .
 * \end{array}
 * \f]
 *
 */
int mess_eigen_arnoldig_inv(mess_direct A, mess_matrix U, mess_matrix W, mess_matrix E, mess_int_t k, mess_vector sv, mess_matrix H, mess_matrix V)
{
    MSG_FNAME(__func__);
    mess_int_t n,i, p = 0;
    // double normsv;
    // double beta, gamma;
    // mess_vector w;
    // mess_vector *Vi;
    // int retV = 0;
    int Atilde = 0;
    mess_matrix Utilde;
    mess_matrix Jay;
    mess_direct JaySolver;
    mess_vector z, u,v;
    mess_vector vj;
    int ret = 0;

    struct mvpdatag mvpdat;
    mess_mvpcall call;


    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(E);
    mess_check_real(E);
    mess_check_nullpointer(sv);
    mess_check_nullpointer(H);
    mess_check_real(sv);
    if ( k == 0) return 0;

    mess_check_positive(k);

    mvpdat.Atilde = 0;
    if ( U != NULL && W != NULL ){
        if (!( U->rows == 0 || U->cols == 0 || W->rows == 0 || W->cols == 0)){
            MSG_INFO("Arnoldi process w.r.t inv(A-UW) * E \n");
            Atilde = 1;
            mvpdat.U=U;
            mvpdat.W=W;
            mvpdat.Atilde = 1;
        }
    }
    n = A->rows;
    if ( k>n) {
        MSG_ERROR("k must be smaller than the order of A!\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  prepare the iteration
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_init(&vj);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(vj,n, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);

    if ( Atilde == 1) {
        mess_matrix tmp;
        p = W->rows;
        MESS_INIT_VECTORS(&z,&u,&v);
        ret = mess_vector_alloc(z, n, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_alloc(u, p, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_alloc(v, p, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        MESS_INIT_MATRICES(&Utilde,&Jay,&tmp);

        // Btilde = A\B
        ret = mess_direct_solvem(MESS_OP_NONE,A, U, Utilde);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solvem);
        ret = mess_matrix_multiply(MESS_OP_NONE, W, MESS_OP_NONE, Utilde, tmp);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
        ret = mess_matrix_alloc(Jay, p, p, p*p, MESS_DENSE, MESS_REAL);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        for ( i = 0; i  < p; i++){
            Jay->values[i+i*p] = 1.0;
        }
        mess_matrix_add(-1.0, tmp, 1, Jay);
        ret = mess_direct_init(&JaySolver);                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_direct_init);
        ret = mess_direct_create_lapack_lu(Jay, JaySolver);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_lapack_lu);
        mess_matrix_clear(&tmp);
        mvpdat.JaySolver = JaySolver;
        mvpdat.z = z;
        mvpdat.v = v;
        mvpdat.u = u;
        mvpdat.Utilde = Utilde;

    }

    mvpdat.Asolver = A;
    mvpdat.E = E ;
    mvpdat.vj = vj ;

    ret = mess_mvpcall_operator(&call, A->rows, A->data_type, mvp_arnoldig_inv, &mvpdat);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_mvpcall_operator );
    ret = mess_eigen_arnoldi_template(call, k, sv, H, V); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template);
    if ( Atilde == 1) {
        mess_vector_clear(&z);
        mess_vector_clear(&u);
        mess_vector_clear(&v);
        mess_matrix_clear(&Utilde);
        mess_matrix_clear(&Jay);
        mess_direct_clear(&JaySolver);
    }

    mess_vector_clear(&vj);
    mess_mvpcall_clear(&call);
    return 0;
}




/**
 * @brief Compute some eigenvalues with the Arnoldi process.
 * @param[in] A   input matrix
 * @param[in] kp     input number of Arnoldi steps with respect to \f$ A \f$
 * @param[in] km     input number of Arnoldi steps with respect to \f$ A^{-1} \f$
 * @param[in] r      input start vector (ones if @c NULL)
 * @param[out] ev   vector containing computed eigenvalues
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref mess_eigen_eigs function computes some eigenvalue approximations stored in \c ev for the eigenvalue problem
 * \f[ Ax = \lambda x \f]
 * with the Arnoldi process.
 *
 */
int mess_eigen_eigs ( mess_matrix A, mess_int_t kp, mess_int_t km, mess_vector r, mess_vector ev )
{
    MSG_FNAME(__func__);
    mess_direct solA = NULL;
    int ret = 0;
    mess_vector intR;
    mess_matrix Hp, Hm;
    mess_vector evp, evm;
    mess_int_t i, j;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(ev);
    mess_check_square(A);
    mess_check_real(A);
    mess_check_nonnegative(kp);
    mess_check_nonnegative(km);


    if ( kp >= A->rows) { kp = A->rows;}
    if ( km >= A->rows) { km = A->rows;}

    if ( kp + km >= A->rows) {
        km = A->rows - kp;
    }

    ret = mess_vector_init(&intR);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(intR, A->rows, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    if ( r ){
        ret = mess_vector_copy(r, intR);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_copy);
    } else {
        ret = mess_vector_ones(intR);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);
    }
    MESS_INIT_MATRICES(&Hp,&Hm);
    MESS_INIT_VECTORS(&evp,&evm);
    ret = mess_vector_alloc(evp, kp, MESS_REAL);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_alloc);
    ret = mess_vector_alloc(evm, km, MESS_REAL);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_alloc);



    /*-----------------------------------------------------------------------------
     *  perform kp steps w.r.t. A
     *-----------------------------------------------------------------------------*/
    if ( kp > 0) {
        ret = mess_eigen_arnoldi(A, NULL, NULL, kp, intR, Hp, NULL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi);
        ret = mess_matrix_resize(Hp, kp, kp);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_resize);
        ret = mess_eigen_eig(Hp, evp, NULL);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eig);
    }


    /*-----------------------------------------------------------------------------
     *  perform km steps w.r.t. inv(A)
     *-----------------------------------------------------------------------------*/
    if ( km > 0) {
        ret = mess_direct_init(&solA);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
        ret = mess_direct_lu(A, solA);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);
        ret = mess_eigen_arnoldi_inv(solA, NULL, NULL, km, intR, Hm, NULL);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_inv);
        ret = mess_matrix_resize(Hm, km, km);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_resize);
        ret = mess_eigen_eig(Hm, evm, NULL);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eig);
        mess_direct_clear(&solA);
    }

    /*-----------------------------------------------------------------------------
     *  finalize
     *-----------------------------------------------------------------------------*/
    for ( i = 0; i < evm->dim; i++) {
        if (MESS_IS_REAL(evm)){
            evm->values[i] = 1.0/evm->values[i];
        } else {
            evm->values_cpx[i] = 1.0/evm->values_cpx[i];
        }
    }

    if ( MESS_IS_COMPLEX(evm) || MESS_IS_COMPLEX(evp)){
        ret = mess_vector_tocomplex(evp);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
        ret = mess_vector_tocomplex(evm);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
        ret = mess_vector_tocomplex(ev);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
    }

    j = 0;
    ret = mess_vector_resize(ev, evp->dim + evm->dim);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);

    if ( MESS_IS_COMPLEX(ev) ){
        for (i = 0; i < evp->dim; i++){
            ev->values_cpx[j++] = evp->values_cpx[i];
        }
        for (i = 0; i < evm->dim; i++){
            ev->values_cpx[j++] = evm->values_cpx[i];
        }
    } else {
        for (i = 0; i < evp->dim; i++){
            ev->values[j++] = evp->values[i];
        }
        for (i = 0; i < evm->dim; i++){
            ev->values[j++] = evm->values[i];
        }
    }

    mess_vector_clear(&intR);
    mess_matrix_clear(&Hp);
    mess_matrix_clear(&Hm);
    mess_vector_clear(&evp);
    mess_vector_clear(&evm);
    return 0;
}       /* -----  end of function mess_eigen_eigs  ----- */

/**
 * @brief Compute some generalized eigenvalues with the Arnoldi process.
 * @param[in] A       input matrix
 * @param[in] E       input matrix
 * @param[in] kp     input number of Arnoldi steps with respect to \f$ A \f$
 * @param[in] km     input number of Arnoldi steps with respect to \f$ A^{-1} \f$
 * @param[in] r      input start vector (ones if @c NULL)
 * @param[out] ev   vector containing computed eigenvalues
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref mess_eigen_eiggs function computes some eigenvalue approximations stored in \c ev for the generalized eigenvalue
 * problem
 * \f[ A x = \lambda E x\f]
 * with the Arnoldi process.
 *
 */
int mess_eigen_eiggs ( mess_matrix A, mess_matrix E, mess_int_t kp, mess_int_t km, mess_vector r, mess_vector ev )
{
    MSG_FNAME(__func__);
    mess_direct solA = NULL;
    int ret = 0;
    mess_vector intR;
    mess_matrix Hp, Hm;
    mess_vector evp, evm;
    mess_int_t i, j;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(E);
    mess_check_nullpointer(ev);
    mess_check_square(A);
    mess_check_real(A);
    mess_check_real(E);
    mess_check_same_size(A,E);
    mess_check_positive(kp);
    mess_check_positive(km);


    if ( kp >= A->rows) { kp = A->rows;}
    if ( km >= A->rows) { km = A->rows;}

    if ( kp + km >= A->rows) {
        km = A->rows - kp;
    }

    ret = mess_vector_init(&intR);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(intR, A->rows, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    if ( r != NULL ){
        ret = mess_vector_copy(r, intR);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_copy);
    } else {
        ret = mess_vector_ones(intR);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);
    }
    MESS_INIT_MATRICES(&Hp,&Hm);
    MESS_INIT_VECTORS(&evp,&evm);
    mess_vector_alloc(evp, kp, MESS_REAL);
    mess_vector_alloc(evm, km, MESS_REAL);


    /*-----------------------------------------------------------------------------
     *  perform kp steps w.r.t. A
     *-----------------------------------------------------------------------------*/
    if ( kp > 0) {
        mess_direct Esolver;
        ret = mess_direct_init(&Esolver);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
        ret = mess_direct_lu(E, Esolver);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);
        ret = mess_eigen_arnoldig(A, NULL, NULL, Esolver, kp, intR, Hp, NULL);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi);
        ret = mess_matrix_resize(Hp, kp, kp);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_resize);
        ret = mess_eigen_eig(Hp, evp, NULL);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eig);
        mess_direct_clear(&Esolver);

    }


    /*-----------------------------------------------------------------------------
     *  perform km steps w.r.t. inv(A)
     *-----------------------------------------------------------------------------*/
    if ( km > 0) {
        ret = mess_direct_init(&solA);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
        ret = mess_direct_lu(A, solA);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_lu);

        ret = mess_eigen_arnoldig_inv(solA, NULL, NULL, E, km, intR, Hm, NULL);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_inv);
        ret = mess_matrix_resize(Hm, km, km);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_resize);
        ret = mess_eigen_eig(Hm, evm, NULL);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eig);
        mess_direct_clear(&solA);
    }

    /*-----------------------------------------------------------------------------
     *  finalize
     *-----------------------------------------------------------------------------*/
    for ( i = 0; i < evm->dim; i++) {
        if (MESS_IS_REAL(evm)){
            evm->values[i] = 1.0/evm->values[i];
        } else {
            evm->values_cpx[i] = 1.0/evm->values_cpx[i];
        }
    }

    if ( MESS_IS_COMPLEX(evm) || MESS_IS_COMPLEX(evp)){
        ret = mess_vector_tocomplex(evp);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
        ret = mess_vector_tocomplex(evm);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
        ret = mess_vector_tocomplex(ev);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_tocomplex);
    }

    j = 0;
    ret = mess_vector_resize(ev, evp->dim + evm->dim);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);

    if ( MESS_IS_COMPLEX(ev) ){
        for (i = 0; i < evp->dim; i++){
            ev->values_cpx[j++] = evp->values_cpx[i];
        }
        for (i = 0; i < evm->dim; i++){
            ev->values_cpx[j++] = evm->values_cpx[i];
        }
    } else {
        for (i = 0; i < evp->dim; i++){
            ev->values[j++] = evp->values[i];
        }
        for (i = 0; i < evm->dim; i++){
            ev->values[j++] = evm->values[i];
        }
    }

    mess_vector_clear(&intR);
    mess_matrix_clear(&Hp);
    mess_matrix_clear(&Hm);
    mess_vector_clear(&evp);
    mess_vector_clear(&evm);
    return 0;
}       /* -----  end of function mess_eigen_eigs  ----- */


