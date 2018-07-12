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
 * @file lib/matrix/condest.c
 * @brief Condition number estimaton.
 * @author @koehlerm
 * @author @mbehr
 *
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


#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif


/**
 * @brief Estimate \f$ 1 \f$-norm  of the inverse of a real square matrix.
 * @param[in] decomp    input contains decomposition of real matrix \f$ A\f$
 * @param[out] nrm      output estimated \f$1\f$-norm condition of \f$A^{-1}\f$
 * @return zero on success or a non zero error code
 *
 * The @ref mess_matrix_norminvest function estimates the \f$ 1 \f$-norm of the inverse
 * of a real square matrix \cite Dav06.
 *
 * @sa mess_matrix_condest
 */
int mess_matrix_norminvest(mess_direct decomp,double *nrm){
    MSG_FNAME(__func__);
    int ret =0;
    mess_int_t n=0,k, j=0,jold=0, m;
    double infnrm=.0;
    double nrm_old=.0;
    mess_vector x,y;
    double est = 0 ;
    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(decomp);
    mess_check_nullpointer(nrm);
    mess_check_nullpointer(decomp->data);

    if(decomp->cols!=decomp->rows){
        MSG_ERROR("Number of Rows and Number of Cols of A must be the same!\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  prepare computation
     *-----------------------------------------------------------------------------*/
    n = decomp->rows;
    MESS_INIT_VECTORS(&x,&y);

    ret = mess_vector_alloc(x,n, MESS_REAL);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_alloc);
    ret = mess_vector_alloc(y,n,MESS_REAL);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_alloc);

    /*-----------------------------------------------------------------------------
     * start computation
     *-----------------------------------------------------------------------------*/
    for(k =1;k<=5;++k){
        if(k==1){
            est=0;
            ret = mess_vector_ones(x);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_ones);
            ret = mess_vector_scale((double)1.0/n,x);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_scale);
            jold=-1;
        }else{
            //compute   j = min(find(abs(x)==norm(x,inf)));
            ret = mess_vector_norminf(x,&infnrm);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norminf)
                for(m=0 ; m<x->dim ; ++m){
                    if(fabs(x->values[m])==infnrm){
                        j=m;
                        break;
                    }
                }

            if(j==jold){
                break;
            }

            ret=mess_vector_zeros(x);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_zeros);
            x->values[j]=1;
            jold=j;
        }

        /*-----------------------------------------------------------------------------
         *  x=Q*(U\(L\(P*x)));
         *-----------------------------------------------------------------------------*/
        ret = mess_direct_solve(MESS_OP_NONE, decomp,x, y); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_direct_solve);
        nrm_old = est;
        ret = mess_vector_norm1(y,&est);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm);
        if(k>1 && est<=nrm_old){
            break ;
        }

        for(m=0; m<n;++m){
            if(y->values[m]<0)
                y->values[m]=-1;
            else
                y->values[m]=1;
        }

        /*-----------------------------------------------------------------------------
         *   x= P' *(L' \(U' \(Q'*x)));
         *-----------------------------------------------------------------------------*/
        ret = mess_direct_solve(MESS_OP_TRANSPOSE,decomp,y,x);  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_solve);

    }

    *nrm = est;
    mess_vector_clear(&x);
    mess_vector_clear(&y);
    return (0);
}


/**
 * @brief Compute a lower bound for the \f$ 1 \f$-norm condition number of a real square matrix.
 * @param[in] A     input real square  matrix \f$A\f$
 * @param[out] nrm  output estimated \f$1\f$-norm condition, if \f$A\f$ is not invertible @c nrm is @c INFINITY
 * @return  zero on success or a non zero error code
 *
 * The @ref mess_matrix_condest function estimates the \f$ 1 \f$-norm condition number of a real square
 * matrix \cite Dav06.
 *
 * @sa mess_matrix_norminvest
 */

int mess_matrix_condest(mess_matrix A, double *nrm){
    MSG_FNAME(__func__);
    int ret=0;
    double nrm1A=.0,nrm1Ainv=.0;
    mess_int_t n;
    mess_direct decomp;
    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_square(A);
    mess_check_real(A);
    n=A->rows;
    if(n==0){
        nrm=0;
        MSG_ERROR("A has no rows!\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     * start computation
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_init(&decomp);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_init);
    ret = mess_direct_lu(A, decomp);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_direct_lu);
    if ( ret != 0 ) {
        // matrix ist singulaer
#ifdef INFINITY
        *nrm =INFINITY;
#else
        *nrm = HUGE_VAL * HUGE_VAL;
#endif
        mess_direct_clear(&decomp);
        return (0);
    }

    ret = mess_matrix_norm1(A,&nrm1A);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_norm1);
    ret = mess_matrix_norminvest(decomp, &nrm1Ainv);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_norminvest);


    //printf("1-Norm A: %lg \t 1-Norm Ainv: %lg \n", nrm1A, nrm1Ainv);
    *nrm = nrm1A*nrm1Ainv;

    mess_direct_clear(&decomp);
    return(0);
}
