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
 * @file lib/eigen/lanczos_template_nrm.c
 * @brief The Lanczos iteration for \f$ 2 \f$-norm computation (real symmetric case).
 * @author @mbehr
 *
 *  This file contains the template for the lanczos iteration (hermitian case).
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>
#include <math.h>



/**
 * @brief Lanczos method for \f$ 2 \f$-norm computations (real symmetric case).
 * @param[in] A  input operator applying the \f$ A \f$ matrix for which the Lanczos algorithm is supposed to run
 * @param[in] k  input number of Lanczos steps (usually \f$ k \ll n \f$)
 * @param[in] sv   input start vector
 * @param[out] diag \f$ k \f$-vector containing diagonal elements of \f$ T \f$
 * @param[out] subdiag \f$ (k-1) \f$-vector containing subdiagonal elements of \f$ T \f$
 * @param[out] V \f$(n \times k)\f$ orthogonal matrix
 * @param[out] abseig absolute value of the largest eigenvalue
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_lanczos_template_nrm function produces matrices \f$ T \f$ (symmetric, tridiagonal) and \f$ V \f$
 * (if \f$ V \ne \f$ @c NULL) fullfilling
 * \f[
 * \begin{array}{ccc}
 *      V(:,1) &=& \frac{sv}{ \Vert sv \Vert }, \\
 *      V^T V  &=& I, \\
 *   AV(:,1:k) &=& VT(:,1:k).
 * \end{array}
 * \f]
 * \f$ T \f$ is stored as  \f$ k \f$-vector \c diag and \f$ (k-1) \f$-vector \c subdiag. \n
 * The algorithms stops if the maximum iteration number \f$ k \f$ is reached or the relative
 * difference between the largest eigenvalues of two step is smaller than \f$ \sqrt{eps} \f$.
 * \n
 *
 */

int mess_eigen_lanczos_template_nrm(mess_mvpcall A, mess_int_t k, mess_vector sv, mess_vector diag, mess_vector subdiag, mess_matrix V, double *abseig)
{
    MSG_FNAME(__func__);
    mess_int_t n,i;
    double normsv;
    double alpha, beta;

    mess_vector t2,t1;

    int retV = 0;
    int ret = 0;
    double me, meold;
    double *diagwork, *subdiagwork;
    mess_int_t dummy =1;
    mess_int_t info;
    double tol = sqrt(mess_eps());
    int err;
    // tol  = IT_TOL;

    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sv);
    mess_check_nullpointer(diag);
    mess_check_nullpointer(subdiag);
    if ( k == 0) return 0;
    mess_check_positive(k);

    if ( V != NULL ){
        retV = 1;
    }

    ret = mess_vector_norm2(sv,&normsv);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    n = A->dim;

    if(normsv == 0.0){
        ret = mess_vector_resize(sv,n);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        ret = mess_vector_ones(sv);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);
        ret = mess_vector_norm2(sv,&normsv);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_norm2);
    }

    if ( k>n) {
        MSG_ERROR("k must be smaller than the order of A!\n");
        return MESS_ERROR_ARGUMENTS;
    }
    // k = 5;
    /*-----------------------------------------------------------------------------
     *  prepare the iteration
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(diagwork,double * ,sizeof(double) * (k+1));
    mess_try_alloc(subdiagwork,double * ,sizeof(double) * (k));

    if( retV==1 ){
        ret = mess_matrix_init(&V);                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
        ret = mess_matrix_alloc(V, n, k, n*k, MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    }

    ret = mess_vector_toreal(diag);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
    ret = mess_vector_toreal(subdiag);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
    ret = mess_vector_resize(diag, k);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    ret = mess_vector_resize(subdiag, k);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);

    MESS_INIT_VECTORS(&t2,&t1);
    ret = mess_vector_alloc(t2, n, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(t1, n, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_scale(1.0/normsv, sv);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_scale);

    if( retV==1 ){
        ret = mess_matrix_setcol(V,0,sv); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
    }

    beta = 1.0;
    meold = 0;
    me = 1;
    err = 0;

    /*-----------------------------------------------------------------------------
     *  begin the iteration algorithm 9.2.1
     *-----------------------------------------------------------------------------*/
    mess_int_t j;
    for(j=0;j < k;++j){
        mess_int_t j2;
        if( j>0){
            /*-----------------------------------------------------------------------------
             *  sv -> t
             *  w/beta -> sv
             *  -beta*t -> w
             *----------------------------------------------------------------------------*/
            ret = mess_vector_copy(sv,t1);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
            ret = mess_vector_scale( 1.0/beta, t2);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_scale);
            ret = mess_vector_copy(t2,sv);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
            ret = mess_vector_scale( -1*beta, t1);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_scale);
            ret = mess_vector_copy(t1,t2);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);

            if( retV==1 ){
                ret = mess_matrix_setcol(V, j, sv); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
            }
        }

        /*-----------------------------------------------------------------------------
         *  perform Ax -> w
         *-----------------------------------------------------------------------------*/
        ret = mess_mvpcall_apply(A,MESS_OP_NONE,  sv, t1);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_apply);
        ret = mess_vector_axpy(1.0, t1, t2);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

        /*-----------------------------------------------------------------------------
         *  compute the values of T
         *-----------------------------------------------------------------------------*/
        //alpha = mess_vector_dot(sv,t2);
        mess_vector_dot(sv,t2,&alpha);
        diag->values[j]=alpha;
        ret = mess_vector_axpy(-1*alpha,sv, t2);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
        ret = mess_vector_norm2(t2,&beta);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_norm2);
        subdiag->values[j]= beta;

        if (beta == 0.0) {
            MSG_INFO("breakdown beta = 0.0 at j = " MESS_PRINTF_INT "\n", j);
            err = MESS_ERROR_CONVERGE;
            break;
        }

        /*-----------------------------------------------------------------------------
         *  the eigenvalue related part
         *-----------------------------------------------------------------------------*/
        meold=me;
        //diag and upperdiag kopieren da inhalt zerstört wird;
        for ( i =0; i <= j; ++i){
            diagwork[i]=diag->values[i];
            subdiagwork[i]=subdiag->values[i];
        }
        j2 = j +1;
        F77_GLOBAL(dstev,DSTEV)("N", &j2, diagwork, subdiagwork, NULL,&dummy,NULL, &info);

        //diaghelp contains ev in ascending order (dstev_ LAPACK documentation)
        me = diagwork[j];

        // MSG_PRINT("it = " MESS_PRINTF_INT " \t me = %.15e \t meold = %.15e \t abs(me-meold) = %.15e\n", j, me,meold, fabs(me-meold)/me);
        if ( fabs(me-meold) < me*tol && j>0) {
            break;
        }

    }

    mess_free(diagwork);
    mess_free(subdiagwork);
    mess_vector_clear(&t1);
    mess_vector_clear(&t2);

    if ( j == 0 && beta == 0.0) {
        *abseig = 0;
    } else {
        *abseig=me;
    }
    // MSG_INFO("Number of Iterations Lanczos: " MESS_PRINTF_INT " \n",j);
    return err;
}

