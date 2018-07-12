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
 * @file lib/direct/splsolve.c
 * @brief Efficient solution with sparse lower trinagular matrices.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include "mess/mess.h"
#include "mess/error_macro.h"

/**
 * @brief Solve \f$ x=G^{-1} B(:,k) \f$ with all sparse.
 * @param[in] G   input lower triangular matrix
 * @param[in] pinv  input inverse row permutation of \f$ G \f$
 * @param[in] B   input right hand side matrix
 * @param[in] k   input number of the column in \f$ B \f$
 * @param[out] x  solution
 * @param[out] top  top of the chi stack
 * @param[out] xi   the reach(L,B(:,col)) stack
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_splsolve fucntion solves  \f$ x=G^{-1} B(:,k) \f$.\n
 * Ported from @csparse. For more details look at @cite Dav06.
 *
 */
int mess_direct_splsolve (mess_matrix G, mess_matrix B, mess_int_t k, mess_int_t *top, mess_int_t *xi, double *x, const mess_int_t *pinv)
{
    MSG_FNAME(__func__);
    mess_int_t j,J, p,q, px, n;
    int ret  =0 ;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(G);
    mess_check_nullpointer(B);
    mess_check_real(G);
    mess_check_real(B);
    mess_check_square(G);
    mess_check_csc(G);
    mess_check_csc(B);
    mess_check_nullpointer(top);
    mess_check_nullpointer(xi);
    mess_check_nullpointer(x);
    if(k < 0 || k >= B->cols) {
        MSG_ERROR("k = " MESS_PRINTF_INT " is out of range\n",k);
        return MESS_ERROR_DIMENSION;
    }
    n = G->rows;
    ret = mess_graph_reach(G,B,k,top,xi, pinv);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_graph_reach);
    for (p = *top ; p < n ; p++) x [xi [p]] = 0 ;    /* clear x */
    for (p = B->colptr[k] ; p < B->colptr[k+1] ; p++) x [B->rowptr[p]] = B->values[p] ; /* scatter B */

    /*-----------------------------------------------------------------------------
     *  solve
     *-----------------------------------------------------------------------------*/
    for (px = *top ; px < n ; px++) {
        j = xi [px] ;                               /* x(j) is nonzero */
        J = pinv ? (pinv [j]) : j ;                 /* j maps to col J of G */
        if (J < 0) continue ;                       /* column J is empty */
        x [j] /= G->values[(G->colptr [J])] ;/* x(j) /= G(j,j) */
        p = G->colptr[J]+1;            /* lo: L(j,j) 1st entry */
        q = G->colptr[J+1] ;        /* up: U(j,j) last entry */
        for ( ; p < q ; p++) {
            x [G->rowptr[p]] -= G->values [p] * x [j] ;          /* x(i) -= G(i,j) * x(j) */
        }
    }
    return 0;                                  /* return top of stack */
}


/**
 * @brief Solve \f$ x=G^{-1} B(:,k) \f$ with all sparse with complex data.
 * @param[in] G  input lower triangular matrix
 * @param[in] pinv  input inverse row permutation of \f$ G \f$
 * @param[in] B  input right hand side matrix
 * @param[in] k  input number of the column in \f$ B \f$
 * @param[out] x    solution
 * @param[out] top  top of the chi stack
 * @param[out] xi   the reach(L,B(:,col)) stack
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_splsolvec function solves \f$ x=G^{-1} B(:,k) \f$ with all sparse with complex data.\n
 * Ported from @csparse. For more details look at @cite Dav06.
 *
 */
int mess_direct_splsolvec (mess_matrix G, mess_matrix B, mess_int_t k, mess_int_t *top, mess_int_t *xi, mess_double_cpx_t *x, const mess_int_t *pinv)
{
    MSG_FNAME(__func__);
    mess_int_t j,J, p,q, px, n;
    int ret  =0 ;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(G);
    mess_check_nullpointer(B);
    mess_check_complex(G);
    mess_check_complex(B);
    mess_check_square(G);
    mess_check_csc(G);
    mess_check_csc(B);
    mess_check_nullpointer(top);
    mess_check_nullpointer(xi);
    mess_check_nullpointer(x);
    if(k < 0 || k >= B->cols) {
        MSG_ERROR("k = " MESS_PRINTF_INT " is out of range\n",k);
        return MESS_ERROR_DIMENSION;
    }
    n = G->rows;
    ret = mess_graph_reach(G,B,k,top,xi, pinv);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_graph_reach);
    for (p = *top ; p < n ; p++) x [xi [p]] = 0 ;    /* clear x */
    for (p = B->colptr[k] ; p < B->colptr[k+1] ; p++) x [B->rowptr[p]] = B->values_cpx[p] ; /* scatter B */

    /*-----------------------------------------------------------------------------
     *  solve
     *-----------------------------------------------------------------------------*/
    for (px = *top ; px < n ; px++) {
        j = xi [px] ;                               /* x(j) is nonzero */
        J = pinv ? (pinv [j]) : j ;                 /* j maps to col J of G */
        if (J < 0) continue ;                       /* column J is empty */
        x [j] /= G->values_cpx[(G->colptr [J])] ;/* x(j) /= G(j,j) */
        p = G->colptr[J]+1;            /* lo: L(j,j) 1st entry */
        q = G->colptr[J+1] ;        /* up: U(j,j) last entry */
        for ( ; p < q ; p++) {
            x [G->rowptr[p]] -= G->values_cpx [p] * x [j] ;          /* x(i) -= G(i,j) * x(j) */
        }
    }
    return 0;                                  /* return top of stack */
}

