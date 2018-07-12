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
 * @file lib/matrix/twoside_mgs.c
 * @brief Skew projection.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/eigenvalue.h"
#include "mess/error_macro.h"
#include <complex.h>

#define REFINE 1
#define KAPPA  0.1


/**
 * @brief Skew projection using  \f$ Q \f$ and \f$ QT \f$.
 * @param[in] in    input vector
 * @param[in] Q      input first subspace to project on
 * @param[in] QT              input second subspace to project on
 * @param[out] out  output vector
 * @param[in] colcount   input column counter for \f$ Q \f$ and \f$ QT \f$\n
 *                       ( number of columns to use from \f$ Q \f$ and \f$ QT \f$)
 * @param[out] nrmout   norm of the output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_proj_mgs function computes the skew projection
 * \f[ out = \Pi_i^{colcount}(I-Q_iQT_i^T)u \f]
 * where \f$ Q_i \f$ is the \f$ i \f$-th column of \f$ Q \f$ and \f$ QT_i \f$ is the \f$ i \f$-th column of \f$ QT \f$.
 *
 */
int  mess_matrix_proj_mgs ( mess_vector in, mess_matrix Q, mess_matrix QT, mess_int_t colcount, mess_vector out, double *nrmout )
{
    MSG_FNAME(__func__);
    mess_vector tmp=NULL;
    int ret = 0;
    mess_int_t i;
    double hu=0;
    double nrm= 0;
    double nrmin = 0;



    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_nullpointer(Q);
    mess_check_nullpointer(QT);
    mess_check_nullpointer(nrmout);
    mess_check_real(in);
    mess_check_real(out);
    mess_check_real(QT);
    mess_check_real(Q);
    if ( Q->rows != in->dim ){
        MSG_ERROR("Q and in have the wrong dimension.\n");
        return ( MESS_ERROR_DIMENSION );
    }
    if ( QT->rows != in->dim ){
        MSG_ERROR("Q and in have the wrong dimension.\n");
        return (MESS_ERROR_DIMENSION);
    }
    if ( colcount > Q->cols || colcount > QT->cols || colcount <0){
        MSG_ERROR("colcount has a wrong value: " MESS_PRINTF_INT "\n", colcount);
        return (MESS_ERROR_DIMENSION);
    }

    ret = mess_vector_norm2(in,&nrmin);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    ret = mess_vector_copy(in,out);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_copy)

        ret = mess_vector_init(&tmp);
        ret = mess_vector_alloc(tmp, Q->rows, MESS_REAL);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_init);

    for ( i = 0; i< colcount; i++){
        ret = mess_matrix_getcol(QT,i,tmp);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_getcol);
        ret = mess_vector_dot(tmp,out,&hu);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_dot);
        ret = mess_matrix_getcol(Q,i,tmp);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
        ret =   mess_vector_axpy(-1.0*hu,tmp, out);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_axpy);
    }


    /*-----------------------------------------------------------------------------
     *  refinement
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_norm2(out,&nrm);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    if ( (REFINE) && (nrm < KAPPA*nrmin) ) {
        // MSG_INFO("refine...\n");
        for ( i = 0; i< colcount; i++){
            ret = mess_matrix_getcol(QT,i,tmp);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_getcol);
            ret = mess_vector_dot(tmp,out,&hu);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_dot);
            ret = mess_matrix_getcol(Q,i,tmp);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
            ret =   mess_vector_axpy(-1.0*hu,tmp, out);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_axpy);
        }
        ret = mess_vector_norm2(out,&nrm);           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    }
    *nrmout = nrm;
    mess_vector_clear(&tmp);
    return(0);
}       /* -----  end of function mess_matrix_proj_mgs  ----- */

/**
 * @brief Generalized skew projection using \f$ Q \f$ and \f$ QT \f$.
 * @param[in] in    input vector
 * @param[in] Q      input first subspace to project on
 * @param[in] QT         input second subspace to project on
 * @param[in] E          input mass matrix
 * @param[in] Eflag      input operation flag for \f$ E \f$
 * @param[out] out  output vector
 * @param[in] colcount   input column counter for \f$ Q \f$ and \f$ QT \f$ \n
 *                     ( number of columns to use from \f$ Q \f$ and \f$ QT \f$)
 * @param[out] nrmout   norm of the output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_proj_mgs function computes the skew projection
 * \f[ out = \Pi_i^{colcount}(I-Q_i op(E) QT_i^T)u \f]
 * where \f$ Q_i \f$ is the \f$ i \f$-th column of \f$ Q \f$ and \f$ QT_i \f$ is the \f$ i \f$-th column of \f$ QT \f$. \n
 * The @c Eflag argument specifies whether \f$ op(E)=E \f$ or \f$ op(E)=E^T \f$ is used.
 *
 */
int  mess_matrix_proj_gmgs ( mess_vector in, mess_matrix Q, mess_matrix QT, mess_matrix E, char *Eflag, mess_int_t colcount, mess_vector out, double *nrmout )
{
    MSG_FNAME(__func__);
    mess_vector tmp=NULL;
    mess_vector tmp2=NULL;
    int ret = 0;
    mess_int_t i;
    double hu=0;
    double nrm= 0;
    double nrmin = 0;
    int Etranspose = 0;



    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_nullpointer(Q);
    mess_check_nullpointer(QT);
    mess_check_nullpointer(nrmout);
    mess_check_real(in);
    mess_check_real(out);
    mess_check_real(QT);
    mess_check_real(Q);
    mess_check_nullpointer(E);
    mess_check_nullpointer(Eflag);
    mess_check_real(E);

    if ( Q->rows != in->dim ){
        MSG_ERROR("Q and in have the wrong dimension.\n");
        return ( MESS_ERROR_DIMENSION);
    }
    if ( QT->rows != in->dim ){
        MSG_ERROR("Q and in have the wrong dimension.\n");
        return (MESS_ERROR_DIMENSION);
    }
    if ( colcount > Q->cols || colcount > QT->cols || colcount <0){
        MSG_ERROR("colcount has a wrong value: " MESS_PRINTF_INT "\n", colcount);
        return (MESS_ERROR_DIMENSION);
    }

    if ( toupper(Eflag[0]) != 'T' && toupper(Eflag[0])!='N'){
        MSG_ERROR("unknown Eflag.\n");
        return (MESS_ERROR_ARGUMENTS);
    } else {
        if (toupper(Eflag[0])=='T'){
            Etranspose = 1;
        } else {
            Etranspose = 0;
        }
    }

    ret = mess_vector_norm2(in,&nrmin);                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    ret = mess_vector_copy(in,out);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_copy)
    MESS_INIT_VECTORS(&tmp,&tmp2);
    ret = mess_vector_alloc(tmp, Q->rows, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_init);
    ret = mess_vector_alloc(tmp2, Q->rows, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_init);

    for ( i = 0; i< colcount; i++){
        ret = mess_matrix_getcol(QT,i,tmp2);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_getcol);
        if (Etranspose == 1){
            ret = mess_matrix_mvp(MESS_OP_HERMITIAN,E,tmp2,tmp);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
        } else {
            ret = mess_matrix_mvp(MESS_OP_NONE,E,tmp2,tmp);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
        }

        ret = mess_vector_dot(tmp,out,&hu);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_dot);
        ret = mess_matrix_getcol(Q,i,tmp);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
        ret = mess_vector_axpy(-1.0*hu,tmp, out);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
    }


    /*-----------------------------------------------------------------------------
     *  refinement
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_norm2(out,&nrm);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    if ( (REFINE) && (nrm < KAPPA*nrmin) ) {
        // MSG_INFO("refine...\n");
        for ( i = 0; i< colcount; i++){
            ret = mess_matrix_getcol(QT,i,tmp2);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_getcol);
            if (Etranspose == 1){
                ret = mess_matrix_mvp(MESS_OP_HERMITIAN,E,tmp2,tmp);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
            } else {
                ret = mess_matrix_mvp(MESS_OP_NONE,E,tmp2,tmp);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
            }

            ret = mess_vector_dot(tmp, out,&hu);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_dot);
            ret = mess_matrix_getcol(Q,i,tmp);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
            ret =   mess_vector_axpy(-1.0*hu,tmp, out);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
        }
        ret = mess_vector_norm2(out,&nrm);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    }
    *nrmout = nrm;
    mess_vector_clear(&tmp);
    mess_vector_clear(&tmp2);
    return (0);
}       /* -----  end of function mess_matrix_proj_mgs  ----- */

