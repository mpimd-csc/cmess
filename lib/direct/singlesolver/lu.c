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
 * @file lib/direct/singlesolver/lu.c
 * @brief Reuse kernel for sparse LU decomposition.
 * @author @koehlerm
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"

#ifdef MESS_HAVE_UMFPACK
#include <umfpack.h>
#endif

#ifdef _OPENMP_H
#include <omp.h>
#endif





#define USE_FORTRAN_KERNELS

#ifdef USE_FORTRAN_KERNELS
void  F77_GLOBAL(lurcsrr,LURCSRR)(mess_int_t *rows, mess_int_t * cols, double *values, mess_int_t *colptr, mess_int_t *rowptr, mess_int_t *lcolptr, mess_int_t *lrowptr, mess_int_t *ucolptr, mess_int_t *urowptr, mess_int_t *havep, mess_int_t *p,
        mess_int_t *haveq, mess_int_t *q, double *lvalues, double *uvalues, double *w);
void  F77_GLOBAL(lurcsrc,LURCSRC)(mess_int_t *rows, mess_int_t * cols, mess_double_cpx_t *values, mess_int_t *colptr, mess_int_t *rowptr, mess_int_t *lcolptr, mess_int_t *lrowptr, mess_int_t *ucolptr, mess_int_t *urowptr, mess_int_t *havep, mess_int_t *p,
        mess_int_t *haveq, mess_int_t *q, mess_double_cpx_t *lvalues, mess_double_cpx_t *uvalues, mess_double_cpx_t *w);

#endif


int mess_decomp_lureuse_kernelcsr(mess_int_t rows, mess_int_t cols,                                 // basic
        double *values, mess_int_t *colptr, mess_int_t *rowptr,         // input matrix
        mess_int_t *lcolptr, mess_int_t *lrowptr,                       // L pattern
        mess_int_t *ucolptr, mess_int_t *urowptr,                       // U pattern
        mess_int_t *p, mess_int_t *q,                                   // permuatations
        double *lvalues, double *uvalues)
{
#ifndef USE_FORTRAN_KERNELS
    MSG_FNAME(__func__);
    double *w;
    mess_int_t i, j, k;
    mess_int_t col;
    mess_int_t row;
    mess_int_t *qinv;
    register double alpha;


    mess_try_alloc(w, double *, sizeof (double) *cols);
    qinv = q;

    // UN(1,:) = A (1,:);
    row = ( p != NULL ) ? p[0] : 0;
    for ( i = rowptr[row];  i < rowptr[row+1]; i++) {
        col = (qinv != NULL) ? qinv[colptr[i]]: colptr[i];
        w [col] = values[i];
    }
    for ( i = urowptr[0];  i < urowptr[1]; i ++ ){
        uvalues[i] = w[ucolptr[i]];
    }
    for ( i = 0; i < cols; i++) {
        w[i] = 0.0;
    }
    lvalues[0] = 1.0;

    // main loop
    for ( i = 1; i < rows; i++) {
        // MSG_INFO("i = " MESS_PRINTF_INT " \n", i);
        /* if ( i & 0x7 ) {
           MSG_INFO(" progress: " MESS_PRINTF_INT " / " MESS_PRINTF_INT "     \r", i, rows);
           } */
        /* % copy a row to the working array
         *  for j = 1:n
         *     if A(i,j) ~= 0
         *         w(j) = A(i,j);
         *     end
         *  end
         */
        row = ( p != NULL ) ? p[i] : i;
        for (j=rowptr[row]; j < rowptr[row+1]; j++){
            col = (qinv != NULL) ? qinv[colptr[j]]: colptr[j];
            // MSG_PRINT("w [ %4ld ] = %lg \n", col, values[j]);
            w[ col ] = values[j];
        }
        /*
         *  for j=1:i-1
         *  if LP(i,j)~=0
         *     alpha = w(j)/UN(j,j);
         *     w = w - alpha*UN(j,:);
         *     w(j) = 0;
         *     LN(i,j) = alpha;
         *  end
         * end
         */
        for (j=lrowptr[i]; j < lrowptr[i+1]-1; j++){
            col = lcolptr[j];
            // MSG_PRINT("w(col) = %lg  \t uvalues = %lg \n", w[col], uvalues[urowptr[col]]);
            alpha = w[col] / uvalues[urowptr[col]];
            if (  uvalues[urowptr[col]] == 0 ) {
                // MSG_PRINT("i=" MESS_PRINTF_INT ", j=" MESS_PRINTF_INT ", alpha = %lg \n",i,j,alpha);
                abort();
            }
            for ( k = urowptr[col]; k < urowptr[col+1]; k++){
                w[ucolptr[k]] -= alpha * uvalues[k];
                // if (w[ucolptr[k]] == 0 ) {
                // MSG_PRINT(" i = " MESS_PRINTF_INT "  , j = " MESS_PRINTF_INT ", k = " MESS_PRINTF_INT " col = " MESS_PRINTF_INT ", ucolptr = " MESS_PRINTF_INT "\n", i , j, k, col, ucolptr[k]);
                // abort();
                // }
            }
            w[col] = 0;
            lvalues[j] = alpha;
            if ( isnan(alpha) || isinf(alpha) ) abort();
        }
        lvalues[j] = 1.0;

        /*
         *  for j = i:n
         *  if UP(i,j) ~= 0
         *     UN(i,j) = w(j);
         *     w(j) = 0;
         *   end
         * end
         */
        for (j=urowptr[i]; j < urowptr[i+1]; j++){
            uvalues[j] = w[ucolptr[j]];
            w[ucolptr[j]] = 0;
        }
    }
    mess_free(w);
#else
    MSG_FNAME(__func__);
    mess_int_t TRUE = 1;
    mess_int_t FALSE = 0;
    mess_int_t havep, haveq;
    double * w;
    mess_int_t j;
    // MSG_INFO("use fortran kernel\n");
    mess_try_alloc(w, double *, sizeof(double) * cols);

    for ( j = 0 ; j < cols ; j++ ) {
        w[j]=0;
    }
    if ( p != NULL) havep = TRUE; else havep = FALSE;
    if ( q != NULL) haveq = TRUE; else haveq = FALSE;
    F77_GLOBAL(lurcsrr,LURCSRR)(&rows, &cols, values, colptr, rowptr, lcolptr, lrowptr, ucolptr, urowptr, &havep, p, &haveq, q, lvalues, uvalues, w);
    mess_free(w);
#endif
    return 0;
}

int mess_decomp_lureuse_kernelcsr_complex(mess_int_t rows, mess_int_t cols,                                 // basic
        mess_double_cpx_t *values, mess_int_t *colptr, mess_int_t *rowptr,         // input matrix
        mess_int_t *lcolptr, mess_int_t *lrowptr,                       // L pattern
        mess_int_t *ucolptr, mess_int_t *urowptr,                       // U pattern
        mess_int_t *p, mess_int_t *q,                                   // permuatations
        mess_double_cpx_t *lvalues, mess_double_cpx_t *uvalues)
{
#ifndef USE_FORTRAN_KERNELS
    MSG_FNAME(__func__);
    mess_double_cpx_t *w;
    mess_int_t i, j, k;
    mess_int_t col;
    mess_int_t row;
    mess_int_t *qinv;
    mess_double_cpx_t alpha;


    mess_try_alloc(w, mess_double_cpx_t *, sizeof (mess_double_cpx_t) *cols);
    qinv = q;

    // UN(1,:) = A (1,:);
    row = ( p != NULL ) ? p[0] : 0;
    for ( i = rowptr[row];  i < rowptr[row+1]; i++) {
        col = (qinv != NULL) ? qinv[colptr[i]]: colptr[i];
        w [col] = values[i];
    }
    for ( i = urowptr[0];  i < urowptr[1]; i ++ ){
        uvalues[i] = w[ucolptr[i]];
    }
    for ( i = 0; i < cols; i++) {
        w[i] = 0.0;
    }
    lvalues[0] = 1.0;

    // main loop
    for ( i = 1; i < rows; i++) {
        /* if ( i & 0x7 ) {
           MSG_INFO(" progress: " MESS_PRINTF_INT " / " MESS_PRINTF_INT "     \r", i, rows);
           } */
        /* % copy a row to the working array
         *  for j = 1:n
         *     if A(i,j) ~= 0
         *         w(j) = A(i,j);
         *     end
         *  end
         */
        row = ( p != NULL ) ? p[i] : i;
        for (j=rowptr[row]; j < rowptr[row+1]; j++){
            col = (qinv != NULL) ? qinv[colptr[j]]: colptr[j];
            w[ col ] = values[j];
        }
        /*
         *  for j=1:i-1
         *  if LP(i,j)~=0
         *     alpha = w(j)/UN(j,j);
         *     w = w - alpha*UN(j,:);
         *     w(j) = 0;
         *     LN(i,j) = alpha;
         *  end
         * end
         */
        for (j=lrowptr[i]; j < lrowptr[i+1]-1; j++){
            col = lcolptr[j];
            alpha = w[col] / uvalues[urowptr[col]];
            for ( k = urowptr[col]; k < urowptr[col+1]; k++){
                w[ucolptr[k]] -= alpha * uvalues[k];
            }
            w[col] = 0;
            lvalues[j] = alpha;
            if ( isnan(alpha) || isinf(alpha) ) abort();

        }
        lvalues[j] = 1.0;

        /*
         *  for j = i:n
         *  if UP(i,j) ~= 0
         *     UN(i,j) = w(j);
         *     w(j) = 0;
         *   end
         * end
         */
        for (j=urowptr[i]; j < urowptr[i+1]; j++){
            uvalues[j] = w[ucolptr[j]];
            w[ucolptr[j]] = 0;
        }
    }
    mess_free(w);
#else
    MSG_FNAME(__func__);
    mess_int_t TRUE = 1;
    mess_int_t FALSE = 0;
    mess_int_t havep, haveq;
    mess_double_cpx_t * w;
    mess_int_t j;
    // MSG_INFO("use fortran kernel\n");
    mess_try_alloc(w, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * cols);
    for ( j = 0 ; j < cols; j++) w[j] = 0;
    if ( p != NULL) havep = TRUE; else havep = FALSE;
    if ( q != NULL) haveq = TRUE; else haveq = FALSE;
    F77_GLOBAL(lurcsrc,LURCSRC)(&rows, &cols, values, colptr, rowptr, lcolptr, lrowptr, ucolptr, urowptr, &havep, p, &haveq, q, lvalues, uvalues, w);
    mess_free(w);
#endif

    return 0;
}


