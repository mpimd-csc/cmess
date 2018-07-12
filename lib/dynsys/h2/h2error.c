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
 * @file lib/dynsys/h2/h2error.c
 * @brief Compute the \f$ \mathcal{H}_2 \f$-error between two dynamical systems.
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

#ifndef INFINITY
#define INFINITY 1.0/0.0
#endif


/**
 * @brief Compute the \f$ \mathcal{H}_2 \f$ error of two LTI systems.
 * @param[in] A      input system matrix
 * @param[in] B     input matrix
 * @param[in] C      input output matrix
 * @param[in] E      input mass Matrix ( @c NULL if not given )
 * @param[in] Ar     input system matrix of the second system
 * @param[in] Br    input matrix of the second system
 * @param[in] Cr     input output matrix of the second system
 * @param[in] Er              input reduced mass matrix. (NULL if not given)
 * @param[out] norm pointer to the double value where the norm is put in
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_h2_error function computes the \f$ \mathcal{H}_2 \f$-error of two SISO LTI systems using
 * Lyapunov Equations. \n
 * The \f$ \mathcal{H}_2 \f$ error is computed as in \ref mess_h2_norm_internal with input matrices
 * \f[ \begin{array}{ccc}
 * A_{big} &=&
 * \left[
 * \begin{array}{cc}
 * A & 0 \\ 0 & Ar
 * \end{array} \right], \\
 * B_{big} &=&
 * \left[
 * \begin{array}{c}
 * B  \\  Br
 * \end{array} \right] , \\
 * C_{big} &=&
 * \left[
 * \begin{array}{cc}
 * C & -Cr
 * \end{array} \right] , \\
 * E_{big} &=&
 * \left[
 * \begin{array}{cc}
 * E & 0 \\ 0 & Er
 * \end{array} \right] .
 * \end{array}
 * \f]
 *
 */
int mess_h2_error_internal ( mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix E, mess_matrix Ar, mess_matrix Br, mess_matrix Cr , mess_matrix Er, double *norm)
{
    MSG_FNAME(__func__);
    mess_matrix Abig=NULL;
    mess_matrix Bbig=NULL;
    mess_matrix Cbig=NULL;
    mess_matrix Ebig=NULL;
    mess_matrix CrM = NULL;
    int ret =0;
    double nrm = 0;
    int haveE=0, haveEr=0;
    int clearEbig = 0, clearE = 0, clearEr=0;
    int isinf = 0;



    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);
    mess_check_nullpointer(Ar);
    mess_check_nullpointer(Br);
    mess_check_nullpointer(Cr);
    mess_check_real(A);
    mess_check_real(B);
    mess_check_real(C);
    mess_check_real(Ar);
    mess_check_real(Br);
    mess_check_real(Cr);
    mess_check_square(A);
    mess_check_square(Ar);
    mess_check_nullpointer(norm);
    if ( E != NULL){
        mess_check_square(E);
        mess_check_real(E);
        haveE=1;
    }
    if ( Er != NULL){
        mess_check_square(Er);
        mess_check_real(Er);
        haveEr=1;
    }


    /*-----------------------------------------------------------------------------
     *  construct matrices
     *-----------------------------------------------------------------------------*/
    ret =mess_matrix_init(&Abig);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret =mess_matrix_init(&Bbig);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret =mess_matrix_init(&Cbig);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

    ret = mess_matrix_cat(A,NULL, NULL,Ar, MESS_CSR,Abig);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);
    ret = mess_matrix_cat(B,NULL,Br,NULL, MESS_DENSE, Bbig);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);

    ret = mess_matrix_init(&CrM);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_copy(Cr,CrM);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_scale(-1.0, CrM);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);

    ret = mess_matrix_cat(C, CrM, NULL, NULL, MESS_DENSE, Cbig);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);

    if ( haveE==0 && haveEr==0){
        Ebig=NULL;
    } else {
        ret = mess_matrix_init(&Ebig);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        if ( haveE==0){
            ret = mess_matrix_init (&E);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
            ret = mess_matrix_eye(E,A->rows, A->cols, A->store_type);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
            clearE=1;
        }
        if (haveEr==0){
            ret = mess_matrix_init (&Er);
            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_init);
            ret = mess_matrix_eye(Er,Ar->rows, Ar->cols, Ar->store_type);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
            clearEr=1;
        }
        ret = mess_matrix_cat(E,NULL, NULL,Er, MESS_CSR,Ebig);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_cat);
        clearEbig = 1;
    }

    if ( MESS_IS_DENSE (Ar)){

        if ( Er != NULL  ){
            mess_vector ew;
            mess_int_t i;
            ret = mess_vector_init(&ew);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
            ret = mess_vector_alloc(ew, Ar->rows, MESS_COMPLEX);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
            ret = mess_eigen_eigg(Ar,Er,ew,NULL,NULL,NULL);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eigg);
            ret = mess_vector_tocomplex(ew);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            // mess_vector_print(ew);
            for ( i=0 ;i< ew->dim ; i++ ){
                if ( creal(ew->values_cpx[i])>0 ) isinf = 1;
            }
            mess_vector_clear(&ew);
        } else {
            mess_vector ew;
            mess_int_t i;
            ret = mess_vector_init(&ew);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
            ret = mess_vector_alloc(ew, Ar->rows, MESS_COMPLEX);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
            ret = mess_eigen_eig(Ar,ew,NULL);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eigg);
            ret = mess_vector_tocomplex(ew);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            // mess_vector_print(ew);
            for ( i=0 ;i< ew->dim ; i++ ){
                if ( creal(ew->values_cpx[i])>0 ) isinf = 1;
            }
            mess_vector_clear(&ew);
        }
    }

    if ( isinf == 0 ) {
        ret = mess_h2_norm_internal(Abig,Bbig,Cbig,Ebig,&nrm );
    } else {
        nrm = INFINITY;
        ret = 0;
    }
    *norm=nrm;
    mess_matrix_clear(&Abig);
    mess_matrix_clear(&Bbig);
    mess_matrix_clear(&Cbig);
    mess_matrix_clear(&CrM);
    if ( clearEbig == 1){
        mess_matrix_clear(&Ebig);
    }
    if ( clearE==1) {
        mess_matrix_clear(&E);
    }
    if (clearEr==1){
        mess_matrix_clear(&Er);
    }
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_h2_norm_internal);
    *norm = nrm;
    return 0;
}


/**
 * @brief Compute the \f$ \mathcal{H}_2 \f$-error between two LTI systems.
 * @param[in] ltiA   input LTI system one
 * @param[in] ltiB   input LTI system two
 * @param[out] err  pointer to the error
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_h2_error function computes the \f$ \mathcal{H}_2 \f$-error between
 * \c ltiA and \c ltiB. \n
 * This function is only a wrapper around \ref mess_h2_error_internal.
 *
 */
int  mess_h2_error ( mess_dynsys ltiA, mess_dynsys ltiB, double *err )
{
    MSG_FNAME(__func__);
    int ret  =0;
    mess_check_nullpointer(ltiA);
    mess_check_nullpointer(ltiB);
    mess_check_nullpointer(err);

    if (!(MESS_IS_DYNSYS_LTI(ltiA) || MESS_IS_DYNSYS_GLTI(ltiA))){
        MSG_ERROR("the first system must be an LTI or a GLTI system.\n");
        return MESS_ERROR_DYNSYS;
    }
    if (!(MESS_IS_DYNSYS_LTI(ltiB) || MESS_IS_DYNSYS_GLTI(ltiB))){
        MSG_ERROR("the second system must be an LTI or a GLTI system.\n");
        return MESS_ERROR_DYNSYS;
    }


    ret = mess_h2_error_internal(ltiA->A, ltiA->B, ltiA->C,ltiA->E, ltiB->A, ltiB->B, ltiB->C, ltiB->E, err);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_h2_error_internal);
    return 0;
}       /* -----  end of function mess_h2_error  ----- */

