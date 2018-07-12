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
 * @file lib/matrix/fdm_matrix.c
 * @brief Generate simple FDM matrices on the unit square.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>

static double zero_func(double x, double y) {
    return 0.0;
}

/**
 * @brief Generate a stiffness matrix \f$ A \f$ for the finite difference discretization (equidistant grid) of the PDE.
 * @param[out] A  \f$ (n \times n) \f$ sparse stiffness matrix, where \f$ n = n0^2 \f$
 * @param[in] n0  input number of inner grid points in each dimension
 * @param[in] fnfx  input pointer to the function describing \f$ f_x \f$ \n
 *                  (if @c NULL, \f$ f_x = 0 \f$ )
 * @param[in] fnfy  input pointer to the function describing \f$ f_y \f$ \n
 *                  (if @c NULL, \f$ f_y = 0 \f$ )
 * @param[in] fng   input pointer to the function describing \f$ g \f$  \n
 *                  (if @c NULL, \f$ g = 0 \f$)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matgen_fdmmatrix function generates the stiffness matrix \f$ A \f$ for the finite difference
 * discretization (equidistant grid) of the PDE
 * \f[ \Delta u - f_x \frac{du}{dx} - f_y \frac{du}{dy} - gu = RHS  \f]
 * on \f$ \Omega \f$ with boundary condition
 * \f[ u  =  0 \f]
 * on \f$ d \Omega \f$, where @f$ \Omega = (0,1) \times (0,1) @f$.
 *
 * This function is just used as an easy way to generate test problems rather than to solve PDEs.
 *
 */
int mess_matgen_fdmmatrix (mess_matrix A, mess_int_t n0, mess_matgen_fdm_function fnfx, mess_matgen_fdm_function fnfy, mess_matgen_fdm_function fng){
    MSG_FNAME(__func__);
    mess_int_t n2;
    double h, h2;
    mess_matrix intA;
    int ret = 0;
    mess_int_t len, ptr = 0, i = 0;
    double t1,t2,t3;
    mess_int_t iy,ix;
    double x,y;
    double fxv, fyv, gv;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_positive(n0);
    if ( fnfx == NULL) fnfx = zero_func;
    if ( fnfy == NULL) fnfy = zero_func;
    if ( fng  == NULL) fng  = zero_func;


    /*-----------------------------------------------------------------------------
     *  setup the matrix
     *-----------------------------------------------------------------------------*/
    n2 = n0*n0;
    h  = 1.0/((double)n0+1.0);
    h2 = h*h;
    t1 = 4.0/h2;
    t2 = -1.0/h2;
    t3 = 1.0/(2.0*h);
    len = 5*n2-4*n0;

    ret = mess_matrix_init(&intA);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_alloc(intA, n2,n2,len,MESS_COORD, MESS_REAL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    for ( iy = 1; iy <= n0; iy++) {
        y = iy * h;
        for (ix =1; ix<=n0; ix++) {
            x = ix * h;
            i++;

            fxv = fnfx(x,y);
            fyv = fnfy(x,y);
            gv = fng(x,y);
            if (iy > 1 ) {
                intA->rowptr[ptr] = i;
                intA->colptr[ptr] = i-n0;
                intA->values[ptr] = t2-fyv*t3;
                ptr ++;
            }
            if (ix > 1) {
                intA->rowptr[ptr] = i;
                intA->colptr[ptr] = i-1;
                intA->values[ptr] = t2-fxv*t3;
                ptr ++;
            }

            intA->rowptr[ptr] = i;
            intA->colptr[ptr] = i;
            intA->values[ptr] = t1+gv;
            ptr ++;

            if ( ix < n0) {
                intA->rowptr[ptr] = i;
                intA->colptr[ptr] = i+1;
                intA->values[ptr] = t2+fxv*t3;
                ptr ++;
            }

            if ( iy < n0) {
                intA->rowptr[ptr] = i;
                intA->colptr[ptr] = i+n0;
                intA->values[ptr] = t2+fyv*t3;
                ptr ++;
            }
        }
    }

    for ( iy = 0 ; iy < len ; iy++  ) {
        intA->rowptr[iy] = intA->rowptr[iy]-1;
        intA->colptr[iy] = intA->colptr[iy]-1;
    }

    ret = mess_matrix_convert(intA, A, MESS_CSR);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    ret = mess_matrix_dupl(A);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_dupl);
    ret = mess_matrix_scale(-1.0, A);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);

    mess_matrix_clear(&intA);
    return 0;
}       /* -----  end of function mess_matgen_fdmmatrix  ----- */

/**
 * @brief Generate a vector \f$ v \f$ containing values of a function \f$ f(x,y) \f$ on an equidistant grid in the interior of the unit square.
 * @param[out] v vector of order \f$ n = n_0^2 \f$ containing the values of \f$ f(x,y) \f$
 * @param[in]  n0  input number of inner grid points in each dimension
 * @param[in]  func  input the function \f$ f \f$ in the space variables \f$ x \f$ and \f$ y \f$ \n
 *             (if @c NULL, \f$ func = 0 \f$)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matgen_fdmvector function generates a vector \f$ v \f$ which contains the values of a function \f$ f(x,y) \f$
 * on an equidistant grid in the interior of the unit square. \n
 * The grid points are numbered consistently with those used in the function
 * \ref mess_matgen_fdmmatrix .
 *
 * This function is just used as an easy way to generate test problems rather than to solve PDEs.
 *
 */
int  mess_matgen_fdmvector ( mess_vector v, mess_int_t n0, mess_matgen_fdm_function func )
{
    MSG_FNAME(__func__);
    double h, y, x;
    mess_int_t n2, i, iy, ix;
    int ret = 0;

    mess_check_positive(n0);
    mess_check_nullpointer(v);
    if ( func == NULL ) func = zero_func;
    h = 1.0 /((double) n0+1);
    n2 = n0*n0;

    ret = mess_vector_toreal(v);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
    ret = mess_vector_resize(v, n2);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    ret = mess_vector_zeros(v);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_zeros);
    i = 0;
    for ( iy = 1 ; iy <= n0; iy++) {
        y = iy *h;
        for ( ix = 1; ix <= n0; ix++){
            x = ix * h ;
            v->values[i] = func(x,y);
            i++;
        }
    }


    return 0;
}       /* -----  end of function mess_matgen_fdmvector  ----- */


/**
 * @brief Generate a one-column matrix containing values of a function \f$ f(x,y) \f$ on an equidistant grid in the interior of the unit square.
 * @param[out] B \f$ (n \times 1) \f$ matrix of order \f$ n = n_0^2 \f$ containing the values of \f$ f(x,y) \f$
 * @param[in]  n0 input number of inner grid points in each dimension
 * @param[in]  func input the function \f$ f \f$ in the space variables \f$ x \f$ and \f$ y \f$ \n
 *                  (if @c NULL, \f$ func = 0 \f$)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matgen_fdmcolumn function generates a \f$ (n \times 1) \f$ matrix \f$ B \f$ which contains the values of a
 * function \f$ f(x,y) \f$ on an equidistant grid in the interior of the unit square. \n
 * The grid points are numbered consistently with those used in the function
 * \ref mess_matgen_fdmmatrix .
 *
 * This function is just used as an easy way to generate test problems rather than to solve PDEs.
 *
 * It is the same as \ref mess_matgen_fdmvector but with an other output format.
 *
 */
int mess_matgen_fdmcolumn ( mess_matrix B, mess_int_t n0, mess_matgen_fdm_function func  )
{
    MSG_FNAME(__func__);
    double h, y, x;
    mess_int_t n2, i, iy, ix;
    int ret = 0;

    mess_check_positive(n0);
    mess_check_nullpointer(B);
    if ( func == NULL ) func = zero_func;
    h = 1.0 /((double) n0+1);
    n2 = n0*n0;

    ret = mess_matrix_alloc(B,n2,1,n2,MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    i = 0;
    for ( iy = 1 ; iy <= n0; iy++) {
        y = iy *h;
        for ( ix = 1; ix <= n0; ix++){
            x = ix * h ;
            B->values[i] = func(x,y);
            i++;
        }
    }

    return 0;
}       /* -----  end of function mess_matgen_fdmcolumn  ----- */

/**
 * @brief Generates a one-row matrix containing values of a function \f$ f(x,y) \f$ on an equidistant grid in the interior of the unit square.
 * @param[out] C \f$ (1 \times n) \f$ matrix of order \f$ n = n_0^2 \f$ containing the values of \f$ f(x,y) \f$
 * @param[in]  n0 input number of inner grid points in each dimension
 * @param[in]  func input the function \f$ f \f$ in the space variables \f$ x \f$ and \f$ y \f$ \n
 *                  (if @c NULL, \f$ func = 0 \f$ )
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matgen_fdmrow function generates a \f$ (1 \times n) \f$ matrix \f$ C \f$ which contains the values of a
 * function \f$ f(x,y) \f$ on an equidistant grid in the interior of the unit square. \n
 * The grid points are numbered consistently with those used in the function
 * \ref mess_matgen_fdmmatrix .
 *
 * This function is just used as an easy way to generate test problems rather than to solve PDEs.
 *
 * It is the same as \ref mess_matgen_fdmvector but with an other output format.
 *
 */
int mess_matgen_fdmrow ( mess_matrix C, mess_int_t n0, mess_matgen_fdm_function func  )
{
    MSG_FNAME(__func__);
    mess_matrix B;
    int ret = 0;

    mess_check_positive(n0);
    mess_check_nullpointer(C);
    if ( func == NULL ) func = zero_func;

    ret = mess_matrix_init(&B);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matgen_fdmcolumn(B, n0, func);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matgen_fdmcolumn);
    ret = mess_matrix_ctranspose(B, C);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
    mess_matrix_clear(&B);

    return 0;
}       /* -----  end of function mess_matgen_fdmcolumn  ----- */
