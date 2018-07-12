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
 * @file lib/matrix/biorth.c
 * @brief Biorthonormalize two matrices.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


/**
 * @brief Construct two biorthonormal matrices from two given matrices.
 * @param[in] Vin       input matrix \f$ V_{in}\f$
 * @param[in] Win       input matrix \f$ W_{in}\f$
 * @param[out] Vout     output matrix \f$V_{out}\f$
 * @param[out] Wout     output matrix \f$W_{out}\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_biorth function constructs a biorthonormal basis
 * \f[ W_{out}^TV_{out}=I \f]
 * from input matrices \f$V_{in}\f$ and \f$W_{in}\f$.
 *
 * @see mess_matrix_gbiorth
 * @attention  The function only works for real and dense matrices.
 *
 */
int  mess_matrix_biorth ( mess_matrix Vin, mess_matrix Win, mess_matrix Vout, mess_matrix Wout )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i=0, m=0;
    mess_vector v,w, tmp, tmp2;
    double nrmv, nrmw, nvw;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(Win);
    mess_check_nullpointer(Wout);
    mess_check_nullpointer(Vin);
    mess_check_nullpointer(Vout);
    mess_check_real(Vin);
    mess_check_real(Win);
    mess_check_dense(Vin);
    mess_check_dense(Win);


    if ( Vin->cols != Win->cols ){
        MSG_ERROR("matrices Vin and Win  must have the same number of columns.\n");
        return(MESS_ERROR_DIMENSION);
    }
    if ( Vin->rows != Win->rows ){
        MSG_ERROR("matrices Vin and Win  must have the same number of columns.\n");
        return(MESS_ERROR_DIMENSION);
    }
    if ( Win->cols > Win->rows ) {
        MSG_ERROR("Biorth only works if rows(" MESS_PRINTF_INT ") >= cols(" MESS_PRINTF_INT ")\n", Win->rows, Win->cols);
        return( MESS_ERROR_DIMENSION);
    }

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    MESS_MATRIX_RESET(Vout);
    MESS_MATRIX_RESET(Wout);
    ret = mess_matrix_alloc(Vout, Vin->rows, Vin->cols, 0, MESS_DENSE, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(Wout, Win->rows, Win->cols, 0, MESS_DENSE, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    MESS_INIT_VECTORS(&v,&w,&tmp,&tmp2);
    ret = mess_vector_alloc(v, Vin->rows, MESS_REAL);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(w, Vin->rows, MESS_REAL);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(tmp, Vin->rows, MESS_REAL);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(tmp2, Vin->rows, MESS_REAL);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);

    m = Vin->cols;
    for ( i = 0 ; i < m ; i++){
        ret = mess_matrix_getcol(Vin, i, tmp);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
        ret = mess_matrix_proj_mgs(tmp, Vout, Wout, i, v,&nrmv);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_proj_mgs);
        mess_vector_scalee(1.0/nrmv,v);

        ret = mess_matrix_getcol(Win, i, tmp2);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
        ret = mess_matrix_proj_mgs(tmp2, Wout, Vout, i, w, &nrmw);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_proj_mgs);
        mess_vector_scalee(1.0/nrmw,w);

        ret = mess_vector_dot(w,v,&nvw);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_dot);
        mess_vector_scalee(1.0/nvw, v);
        ret = mess_matrix_setcol(Vout, i, v);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
        ret = mess_matrix_setcol(Wout, i, w);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
    }
    mess_vector_clear(&w);
    mess_vector_clear(&v);
    mess_vector_clear(&tmp);
    mess_vector_clear(&tmp2);

    return(0);
}       /* -----  end of function mess_matrix_biorth  ----- */


/**
 * @brief Construct two biorthonormal matrices from two given matrices with respect to a weighted scalar product.
 * @param[in] E         input matrix \f$E\f$ symmetric positive definite
 * @param[in] Vin       input matrix \f$V_{in}\f$
 * @param[in] Win       input matrix \f$W_{in}\f$
 * @param[out] Vout     output matrix \f$V_{out}\f$
 * @param[out] Wout     output matrix \f$W_{out}\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_gbiorth function constructs a biorthonormal basis
 * \f[ W_{out}^TEV_{out}=I \f]
 * from the input matrices \f$V_{in}\f$ and \f$W_{in}\f$. \n
 * \f$ E \f$ need to be a symmetric positive definite matrix defining a proper scalar product.
 *
 * @see mess_matrix_biorth
 * @attention The function only works for real and dense matrices.
 *
 */
int  mess_matrix_gbiorth ( mess_matrix E, mess_matrix Vin, mess_matrix Win, mess_matrix Vout, mess_matrix Wout )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i=0, m=0;
    mess_vector v,w, tmp, tmp2, Ev;
    double nrmv, nrmw, nvw;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(E);
    mess_check_nullpointer(Win);
    mess_check_nullpointer(Wout);
    mess_check_nullpointer(Vin);
    mess_check_nullpointer(Vout);
    mess_check_real(Vin);
    mess_check_real(Win);
    mess_check_dense(Vin);
    mess_check_dense(Win);
    mess_check_real(E);
    mess_check_square(E);

    if ( Vin->cols != Win->cols ){
        MSG_ERROR("matrices Vin and Win  must have the same number of columns.\n");
        return (MESS_ERROR_DIMENSION);
    }
    if ( Vin->rows != Win->rows ){
        MSG_ERROR("matrices Vin and Win  must have the same number of columns.\n");
        return (MESS_ERROR_DIMENSION);
    }
    if ( E->rows != Win->rows){
        MSG_ERROR("dimension of the E matrix doesn't match.\n");
        return (MESS_ERROR_DIMENSION);
    }

    if ( Win->cols > Win->rows ) {
        MSG_ERROR("Biorth only works if rows(" MESS_PRINTF_INT ") >= cols(" MESS_PRINTF_INT ")\n", Win->rows, Win->cols);
        return ( MESS_ERROR_DIMENSION);
    }


    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    MESS_MATRIX_RESET(Vout);
    MESS_MATRIX_RESET(Wout);
    ret = mess_matrix_alloc(Vout, Vin->rows, Vin->cols, 0, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    ret = mess_matrix_alloc(Wout, Win->rows, Win->cols, 0, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    MESS_INIT_VECTORS(&v,&w,&tmp,&tmp2,&Ev);
    ret = mess_vector_alloc(v, Vin->rows, MESS_REAL);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(w, Vin->rows, MESS_REAL);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(tmp, Vin->rows, MESS_REAL);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(tmp2, Vin->rows, MESS_REAL);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(Ev, Vin->rows, MESS_REAL);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);

    m = Vin->cols;
    for ( i = 0 ; i < m ; i++){
        ret = mess_matrix_getcol(Vin, i, tmp);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
        ret = mess_matrix_proj_gmgs(tmp, Vout, Wout, E,"T",i, v,&nrmv);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_proj_mgs);
        mess_vector_scalee(1.0/nrmv,v);

        ret = mess_matrix_getcol(Win, i, tmp2);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
        ret = mess_matrix_proj_gmgs(tmp2, Wout, Vout,E,"N", i, w, &nrmw);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_proj_mgs);
        mess_vector_scalee(1.0/nrmw,w);

        ret = mess_matrix_mvp(MESS_OP_NONE,E,v,Ev);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_mvp);
        ret = mess_vector_dot(w,Ev,&nvw);                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_dot);
        mess_vector_scalee(1.0/nvw, v);
        ret = mess_matrix_setcol(Vout, i, v);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
        ret = mess_matrix_setcol(Wout, i, w);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
    }
    mess_vector_clear(&w);
    mess_vector_clear(&v);
    mess_vector_clear(&tmp);
    mess_vector_clear(&tmp2);
    mess_vector_clear(&Ev);

    return 0;
}       /* -----  end of function mess_matrix_gbiorth  ----- */

