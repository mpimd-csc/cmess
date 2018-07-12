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
 * @file lib/matrix/mgsadd.c
 * @brief Add columns to an orthogonal matrix.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>
#include <math.h>



/**
 * @brief Use the modified Gram-Schmidt orthogonalization to add new columns to an orthogonal matrix.
 * @param[in,out]   Z   input/output orthogonal matrix \f$Z\f$
 * @param[in]       V   input new columns to add \f$V\f$
 * @param[in]       E   input mass matrix \f$E\f$, @c NULL if not given
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_mgs_add function computes an orthonormal basis for
 * \f$ \left[ \begin{array}{cc} Z & V \end{array}\right]\f$ with
 * the modified Gram-Schmidt process under the assumption that \f$ Z \f$ is already
 * orthogonal. \n
 * If \f$ E \f$ is given, we use the \f$ E \f$ inner product. In this case \f$ E \f$ need to be symmetric positive definite.
 * The result is written to \f$ Z \f$.
 *
 * @sa mess_matrix_mgs
 * @sa mess_matrix_mgs_inplace
 * @sa mess_matrix_orth
 */
int  mess_matrix_mgs_add ( mess_matrix Z, mess_matrix V, mess_matrix E )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t k,j;
    int haveE = 0;
    mess_int_t nz,nv;
    double kappa = 0.1;
    mess_int_t refine = 1, refinecnt;
    double nq;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(Z);
    mess_check_nullpointer(V);
    mess_check_dense(Z);
    mess_check_dense(V);
    mess_check_real_or_complex(Z);
    mess_check_real_or_complex(V);
    if ( E!=NULL){
        mess_check_real(E);
        mess_check_square(E);
        haveE=1;
    }
    nz = Z->cols;
    nv = V->cols;

    if ( nz == 0 ) {
        ret = mess_matrix_copy(V, Z);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    } else {
        ret = mess_matrix_addcols(Z,V); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addcols);
    }
    // MSG_INFO("mz= %d \t nz=%d \t mv = %d \t nv = %d\n", mz, nz, mv, nv);
    // if ( haveE ) refine = 2;

    if(MESS_IS_REAL(Z) &&  MESS_IS_REAL(V)){
        //real case
        double rw;
        if (haveE){
            for (k = 0 ; k < nv ; k++){
                ret = mess_matrix_colnorm(Z, nz + k, &nq );                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnormE);

                for (j=0; j < nz ; j++) {
                    ret = mess_matrix_coldotE(Z,E,j,nz+k, &rw);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldotE);
                    ret = mess_matrix_colaxpy2(Z,-rw,j,nz+k);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
                }

                ret = mess_matrix_colnorm(Z,nz+k, &rw);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);
                if ( rw  < kappa *nq) {
                    // MSG_WARN("refine\n");
                    for ( refinecnt =0 ; refinecnt < refine; refinecnt++){
                        for (j=0; j < nz ; j++) {
                            ret = mess_matrix_coldotE(Z,E,j,nz+k, &rw);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldotE);
                            ret = mess_matrix_colaxpy2(Z,-rw,j,nz+k);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_colaxpy2);
                        }
                    }
                }

                // modify the new columns
                ret = mess_matrix_colnormE(Z,E,nz+k, &rw);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnormE);
                // MSG_INFO("normEi k = % ld  \t w = %g\n",k, w);
                ret =   mess_matrix_colscale(Z, nz+k, 1.0/rw);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colscale);

                for (j=k+1; j <nv ; j++){
                    ret = mess_matrix_coldotE(Z,E,nz+k, nz+j, &rw);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldotE);
                    ret = mess_matrix_colaxpy2(Z,-rw,nz+k, nz+j);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
                }

            }
        }else{
            for (k = 0 ; k < nv ; k++){
                ret = mess_matrix_colnorm(Z, nz + k, &nq );                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);

                for (j=0; j < nz ; j++) {
                    ret = mess_matrix_coldot(Z,j,nz+k, &rw);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldot);
                    ret = mess_matrix_colaxpy2(Z,-rw,j,nz+k);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
                }

                ret = mess_matrix_colnorm(Z,nz+k, &rw);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);
                if ( rw  < kappa *nq) {
                    // MSG_WARN("refine\n");
                    for ( refinecnt =0 ; refinecnt < refine; refinecnt++){
                        for (j=0; j < nz ; j++) {
                            ret = mess_matrix_coldot(Z,j,nz+k, &rw);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldot);
                            ret = mess_matrix_colaxpy2(Z,-rw,j,nz+k);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_colaxpy2);
                        }
                    }
                }

                // modify the new columns

                ret = mess_matrix_colnorm(Z, nz+k, &rw);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);
                // MSG_INFO("normEi k = % ld  \t w = %g\n",k, w);
                ret =   mess_matrix_colscale(Z, nz+k, 1.0/rw);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colscale);

                for (j=k+1; j <nv ; j++){
                    ret = mess_matrix_coldot(Z,nz+k, nz+j, &rw);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldot);
                    ret = mess_matrix_colaxpy2(Z,-rw,nz+k, nz+j);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
                }

            }
        }
    }else{
        //complex case
        mess_double_cpx_t cw;
        double nrm;
        if (haveE){
            for (k = 0 ; k < nv ; k++){
                ret = mess_matrix_colnormE(Z,E, nz + k, &nq );                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnormE);

                for (j=0; j < nz ; j++) {
                    ret = mess_matrix_coldotcE(Z,E,j,nz+k, &cw);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldotE);
                    ret = mess_matrix_colaxpy2(Z,-cw,j,nz+k);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
                }

                ret = mess_matrix_colnormE(Z,E,nz+k, &nrm);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);
                if ( nrm  < kappa *nq) {
                    // MSG_WARN("refine\n");
                    for ( refinecnt =0 ; refinecnt < refine; refinecnt++){
                        for (j=0; j < nz ; j++) {
                            ret = mess_matrix_coldotcE(Z,E,j,nz+k, &cw);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldotE);
                            ret = mess_matrix_colaxpy2(Z,-cw,j,nz+k);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_colaxpy2);
                        }
                    }
                }

                // modify the new columns
                ret = mess_matrix_colnormE(Z,E,nz+k, &nrm);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnormE);
                // MSG_INFO("normEi k = % ld  \t w = %g\n",k, w);
                ret =   mess_matrix_colscale(Z, nz+k, 1.0/nrm);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colscale);

                for (j=k+1; j <nv ; j++){
                    ret = mess_matrix_coldotcE(Z,E,nz+k, nz+j, &cw);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldotE);
                    ret = mess_matrix_colaxpy2(Z,-cw,nz+k, nz+j);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
                }

            }
        }else{

            for (k = 0 ; k < nv ; k++){
                ret = mess_matrix_colnorm(Z, nz + k, &nq );                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);

                for (j=0; j < nz ; j++) {
                    ret = mess_matrix_coldotc(Z,j,nz+k, &cw);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldot);
                    ret = mess_matrix_colaxpy2(Z,-cw,j,nz+k);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
                }

                ret = mess_matrix_colnorm(Z,nz+k, &nrm);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);
                if ( nrm  < kappa *nq) {
                    // MSG_WARN("refine\n");
                    for ( refinecnt =0 ; refinecnt < refine; refinecnt++){
                        for (j=0; j < nz ; j++) {
                            ret = mess_matrix_coldotc(Z,j,nz+k, &cw);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldot);
                            ret = mess_matrix_colaxpy2(Z,-cw,j,nz+k);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_colaxpy2);
                        }
                    }
                }

                // modify the new columns

                ret = mess_matrix_colnorm(Z, nz+k, &nrm);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colnorm);
                // MSG_INFO("normEi k = % ld  \t w = %g\n",k, w);
                ret =   mess_matrix_colscale(Z, nz+k, 1.0/nrm);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colscale);
                for (j=k+1; j <nv ; j++){
                    ret = mess_matrix_coldotc(Z,nz+k, nz+j, &cw);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_coldot);
                    ret = mess_matrix_colaxpy2(Z,-cw,nz+k, nz+j);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy2);
                }

            }
        }


    }
    return (0);
}       /* -----  end of function mess_matrix_mgsadd  ----- */

