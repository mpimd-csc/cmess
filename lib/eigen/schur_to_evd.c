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
 * @file lib/eigen/schur_to_evd.c
 * @brief Computes Spectral Decomposition from Schur decomposition.
 * @author @koehlerm
 * @author @mbehr
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
#include "blas_defs.h"
#include <complex.h>


/**
 * @brief Compute eigenvectors from a given Schur decomposition-.
 * @param[in] T     input quasi upper triangular matrix of Schur decomposition
 * @param[in] Q     input unitary transformation matrix of Schur decomposition
 * @param[in] V     input matrix with eigenvectors as columns
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_schur_to_evd function computes the eigenvectors
 * from a given Schur decomposition \f$ A = Q * T * Q^H \f$.
 * The eigenvectors are the columns of the matrix @p V at return.
 * Pairs of complex conjugated eigenvectors are stored in consecutive
 * columns.
 * The function calls @c DTREVC / @c ZTRVEC from @lapack.
 *
 * @sa mess_eigen_hessenberg_schur
 *
 */
int mess_eigen_schur_to_evd(mess_matrix T, mess_matrix Q, mess_matrix V){
    MSG_FNAME(__func__);

    int ret = 0;
    char SIDE[] = "R", HOWMNY[] = "B";
    mess_int_t M, info=0, i=0, j=0, cpx=0;
    double *WORK=NULL;
    mess_double_cpx_t *WORK_CPX=NULL;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(T);
    mess_check_nullpointer(Q);
    mess_check_nullpointer(V);

    mess_check_dense(T); mess_check_real_or_complex(T); mess_check_square(T);
    mess_check_dense(Q); mess_check_real_or_complex(Q); mess_check_square(Q);

    mess_check_same_size(T,Q);
    mess_check_same_datatype(T,Q);

    /*-----------------------------------------------------------------------------
     *  prepare dtrevc/ztrvec
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_copy(Q,V);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);

    if(MESS_IS_REAL(T)){
        //alloc workspace for dtrvec
        mess_try_alloc(WORK, double*, 3*sizeof(double)*T->rows);
        //call
        F77_GLOBAL(dtrevc,DTREVC)(SIDE, HOWMNY, NULL, &(T->rows), T->values, &(T->ld), NULL, &(V->cols) /*dummy*/, V->values, &(V->ld), &(V->cols), &M, WORK, &info);
        if(info){goto clear;}
    }else{
        //alloc workspace for ztrvec
        mess_try_alloc(WORK_CPX, mess_double_cpx_t*, 2*sizeof(mess_double_cpx_t)*T->rows);
        mess_try_alloc(WORK, double*, sizeof(double)*T->rows);
        //call
        F77_GLOBAL(ztrevc,ZTREVC)(SIDE, HOWMNY, NULL, &(T->rows), T->values_cpx, &(T->ld), NULL, &(V->cols) /*dummy*/, V->values_cpx, &(V->ld), &(V->cols), &M, WORK_CPX, WORK, &info);
        if(info){goto clear;}
    }


    /*-----------------------------------------------------------------------------
     *  postprocess eigenvectors in real case
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(T)){
        // check subdiagonal for complex eigenvalues
        for(i=0; i<T->cols-1; ++i){
            if(T->values[i+1 + i*T->ld]){
                cpx=1;
                break;
            }
        }

        //convert if necessary and postprocess eigenvectors
        if(cpx){
            ret = mess_matrix_tocomplex(V);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_tocomplex);
            for(i=0; i<T->cols-1; ++i){
                if(T->values[i+1 + i*T->ld]){
                    // 2-x-2 block postprocess the eigenvectors
                    for(j=0;j<V->rows;++j){
                        V->values_cpx[j+i*T->ld]        += V->values_cpx[j+(i+1)*T->ld]*I;
                        V->values_cpx[j+(i+1)*T->ld]    = conj(V->values_cpx[j+i*T->ld]);
                    }
                }
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  clear stuff
     *-----------------------------------------------------------------------------*/
clear:
    if(WORK)        mess_free(WORK);
    if(WORK_CPX)    mess_free(WORK_CPX);

    /*-----------------------------------------------------------------------------
     *  error if info !=0
     *-----------------------------------------------------------------------------*/
    if(info){
        MSG_ERROR("dtrevc/ztrevc returned with error = " MESS_PRINTF_INT "\n", info);
        return MESS_ERROR_LAPACK;
    }

    return 0;
}       /* -----  end of function mess_eigen_schur_to_evd  ----- */

