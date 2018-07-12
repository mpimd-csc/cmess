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
 * @file lib/lrcf_adi/norms/rcFnorm.c
 * @brief Frobenius-norm difference between factored solutions.
 * @author @koehlerm
 *
 * This file contains various functions to perform matrix-matrix additions and
 * similar problems.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/eigenvalue.h"
#include "mess/error_macro.h"
#include <complex.h>


/**
 * @brief Compute the relative Frobenius-norm difference between two factored solutions of a matrix equation.
 * @param[in] Zold  input old \f$ Z \f$ factor
 * @param[in] Z  input new \f$ Z \f$ factor
 * @param[out] chg relative Frobenius-norm difference
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_facFdiff function computes the relative Frobenius-norm difference
 * \f[ chg = \Vert \underbrace{Z_{old}Z_{old}^T}_{X_{old}} - \underbrace{ZZ^T}_{X} \Vert_F \cdot \Vert X \Vert_{F^{-1}} \f]
 * without forming \f$ X_{old}\f$ and \f$ X \f$.
 *
 *
 */
int  mess_lrcfadi_facFdiff ( mess_matrix Zold, mess_matrix Z, double *chg )
{
    MSG_FNAME(__func__);
    mess_matrix M1,M1t, M2t, M2, M3,M3t;
    double nrm, diffnrm, tr2, tr3;
    mess_int_t i, j;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(Zold);
    mess_check_nullpointer(Z);
    mess_check_nullpointer(chg);

    ret = mess_matrix_init(&M1);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&M2);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&M1t);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&M2t);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&M3);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_init(&M3t);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, Z, MESS_OP_NONE, Z, M1);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, Zold, MESS_OP_NONE, Zold, M2);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    ret = mess_matrix_multiply(MESS_OP_HERMITIAN, Zold, MESS_OP_NONE, Z, M3);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);


    if ( MESS_IS_REAL(M1)) {
        nrm=0;
#ifdef _OPENMP
#pragma omp parallel for private(i,j) reduction(+:nrm)
#endif
        for(i=0; i < M1->rows; i++){
            for (j=0; j < M1->rows; j++) {  nrm+= M1->values[i*M1->rows+j]*M1->values[j*M1->rows+i];    }
        }

    } else {
        nrm = 0;
#ifdef _OPENMP
#pragma omp parallel for private(i,j) reduction(+:nrm)
#endif
        for(i=0; i < M1->rows; i++){
            for (j=0; j < M1->rows; j++) {  nrm+= creal(M1->values_cpx[i*M1->rows+j]*M1->values_cpx[j*M1->rows+i]);}
        }

    }

    if ( MESS_IS_REAL(M2)) {
        tr2=0;
#ifdef _OPENMP
#pragma omp parallel for private(i,j) reduction(+:tr2)
#endif
        for(i=0; i < M2->rows; i++){
            for (j=0; j < M2->rows; j++) {  tr2+= M2->values[i*M2->rows+j]*M2->values[j*M2->rows+i];    }
        }

    } else {
        tr2 = 0;
#ifdef _OPENMP
#pragma omp parallel for private(i,j) reduction(+:tr2)
#endif
        for(i=0; i < M2->rows; i++){
            for (j=0; j < M2->rows; j++) {  tr2+= creal(M2->values_cpx[i*M2->rows+j]*M2->values_cpx[j*M2->rows+i]);}
        }

    }

    // mess_matrix_printinfo(M3);
    if ( MESS_IS_REAL(M3)) {
        tr3=0;
#ifdef _OPENMP
#pragma omp parallel for private(i,j) reduction(+:tr3)
#endif
        for(i=0; i < M3->rows; i++){
            for (j=0; j < M3->cols; j++) {  tr3+= M3->values[j*M3->rows+i]*M3->values[j*M3->rows+i];    }
        }

    } else {
        tr3 = 0;
#ifdef _OPENMP
#pragma omp parallel for private(i,j) reduction(+:tr3)
#endif
        for(i=0; i < M3->rows; i++){
            for (j=0; j < M3->cols; j++) {  tr3+= creal(M3->values_cpx[j*M3->rows+i]*conj(M3->values_cpx[j*M3->rows+i]));}
        }

    }

    mess_matrix_clear(&M1);
    mess_matrix_clear(&M2);
    mess_matrix_clear(&M3);
    mess_matrix_clear(&M1t);
    mess_matrix_clear(&M2t);
    mess_matrix_clear(&M3t);


    diffnrm = nrm+tr2-2*tr3;
    if ( diffnrm < 0) MSG_WARN("diffnrm = %.17e !!!!\n", diffnrm);
    if(fabs(nrm)< mess_eps()){
        MSG_ERROR("nrm =%e \t is smaller than eps",nrm);
        return MESS_ERROR_SINGULAR;
    }
    *chg = sqrt(fabs(diffnrm/nrm));


    return 0;
}       /* -----  end of function mess_lrcfadi_rcFnorm  ----- */

