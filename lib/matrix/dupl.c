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
 * @file lib/matrix/dupl.c
 * @brief Remove duplicate entries and zero entries from sparse storage schemes.
 * @author @koehlerm
 * @author @mbehr
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>

/**
 * @brief Detect and sum up duplicate entries in a sparse matrix.
 * @param[in,out] mat matrix from which the duplicate entries should be removed
 * @return 0 on success or a non zero error code otherwise.
 *
 * The @ref mess_matrix_dupl function detects and sums up duplicate entries in CSR/CSC
 * matrices, that means if there exists two or more entries in one row or column which
 * have the same column or row indicies.\n
 * The function assumes that the data structures are sorted. If this
 * is not fulfilled the matrix can be sorted using \ref mess_matrix_sort before.
 *
 * @sa mess_matrix_sort
 */
int  mess_matrix_dupl ( mess_matrix mat )
{
    MSG_FNAME(__func__);
    mess_int_t *Ap, *Ai, *w;
    double *Ax;
    mess_double_cpx_t *Az;
    mess_int_t i,j,p,nz = 0, m, n,q;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(mat);
    mess_check_real_or_complex(mat);

    if ( MESS_IS_DENSE(mat)) {
        return 0;
    }
    if ( !(MESS_IS_CSR(mat) || MESS_IS_CSC(mat))) {
        MSG_ERROR("The matrix has the wrong storage type. Only CSR or CSC allowed.\n");
        return MESS_ERROR_STORAGETYPE;
    }
    m = mat->rows;
    n = mat->cols;

    /*-----------------------------------------------------------------------------
     *  CSC variant
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_CSC(mat)) {
        if(MESS_IS_REAL(mat)){
            Ap = mat->colptr;
            Ai = mat->rowptr;
            Ax = mat->values;
            mess_try_alloc(w, mess_int_t *, sizeof(mess_int_t) * m);
            for ( i = 0; i < m; i++ ) w[i] = -1;
            for ( j = 0; j < n; j++){
                q = nz;
                for ( p = Ap[j]; p < Ap[j+1]; p++){
                    i = Ai[p];
                    if ( w[i] >= q) {
                        Ax[w[i]] += Ax[p];
                    } else {
                        w[i] = nz;
                        Ai[nz] = i;
                        Ax[nz++] = Ax[p];
                    }
                }
                Ap[j] = q;
            }
            Ap[n] = nz;
            mess_free(w);
            if ( nz > 0 ) {
                mess_try_realloc(mat->values, double *, nz * sizeof(double));
                mess_try_realloc(mat->rowptr, mess_int_t*, nz*sizeof(mess_int_t));
            } else {
                mess_free(mat->values);
                mess_free(mat->rowptr);
                mat->values = NULL;
                mat->rowptr = NULL;
            }
            mat->nnz = nz;
        }else{
            Ap = mat->colptr;
            Ai = mat->rowptr;
            Az = mat->values_cpx;
            mess_try_alloc(w, mess_int_t *, sizeof(mess_int_t) * m);
            for ( i = 0 ; i < m ; i++ ) w[i] = -1;
            for ( j = 0; j < n ; j++){
                q = nz;
                for ( p = Ap[j]; p < Ap[j+1]; p++){
                    i = Ai[p];
                    if ( w[i] >= q) {
                        Az[w[i]] += Az[p];
                    } else {
                        w[i] = nz;
                        Ai[nz] = i;
                        Az[nz++] = Az[p];
                    }
                }
                Ap[j] = q;
            }
            Ap[n] = nz;
            mess_free(w);
            if ( nz > 0 ) {
                mess_try_realloc(mat->values_cpx, mess_double_cpx_t *, nz * sizeof(mess_double_cpx_t));
                mess_try_realloc(mat->rowptr, mess_int_t*, nz*sizeof(mess_int_t));
            } else {
                mess_free(mat->values_cpx);
                mess_free(mat->rowptr);
                mat->values = NULL;
                mat->rowptr = NULL;
            }
            mat->nnz = nz;

        }
    } else if ( MESS_IS_CSR(mat)){
        if(MESS_IS_REAL(mat)){
            Ap = mat->rowptr;
            Ai = mat->colptr;
            Ax = mat->values;
            mess_try_alloc(w, mess_int_t *, sizeof(mess_int_t) * n);
            for ( i = 0 ; i < n ; i++ ) w[i] = -1;
            for ( j = 0; j < m ; j++){
                q = nz;
                for ( p = Ap[j]; p < Ap[j+1]; p++){
                    i = Ai[p];
                    if ( w[i] >= q) {
                        Ax[w[i]] += Ax[p];
                    } else {
                        w[i] = nz;
                        Ai[nz] = i;
                        Ax[nz++] = Ax[p];
                    }
                }
                Ap[j] = q;
            }
            Ap[m] = nz;
            mess_free(w);
            if ( nz > 0 ) {
                mess_try_realloc(mat->values, double *, nz * sizeof(double));
                mess_try_realloc(mat->colptr, mess_int_t*, nz*sizeof(mess_int_t));
            } else {
                mess_free(mat->values);
                mess_free(mat->colptr);
                mat->values = NULL;
                mat->colptr = NULL;
            }
            mat->nnz = nz;

        }else{
            Ap = mat->rowptr;
            Ai = mat->colptr;
            Az = mat->values_cpx;
            mess_try_alloc(w, mess_int_t *, sizeof(mess_int_t) * n);
            for ( i = 0 ; i < n ; i++ ) w[i] = -1;
            for ( j = 0; j < m ; j++){
                q = nz;
                for ( p = Ap[j]; p < Ap[j+1]; p++){
                    i = Ai[p];
                    if ( w[i] >= q) {
                        Az[w[i]] += Az[p];
                    } else {
                        w[i] = nz;
                        Ai[nz] = i;
                        Az[nz++] = Az[p];
                    }
                }
                Ap[j] = q;
            }
            Ap[m] = nz;
            mess_free(w);
            if ( nz > 0 ) {
                mess_try_realloc(mat->values_cpx, mess_double_cpx_t *, nz * sizeof(mess_double_cpx_t));
                mess_try_realloc(mat->colptr, mess_int_t*, nz*sizeof(mess_int_t));
            } else {
                mess_free(mat->values_cpx);
                mess_free(mat->colptr);
                mat->values = NULL;
                mat->colptr = NULL;
            }
            mat->nnz = nz;

        }
    }
    return 0;
}       /* -----  end of function mess_matrix_dupl  ----- */




/**
 * @brief Remove zero entries in a sparse matrix.
 * @param[in,out] mat matrix from which the zero entries should be removed
 * @return 0 on success or a non zero error code otherwise.
 *
 * The @ref mess_matrix_eliminate_zeros function detects and removes zero entries in CSR/CSC/COORD
 * matrices.\n
 * Please note that a dense matrix is immediately returned.
 *
 */
int  mess_matrix_eliminate_zeros ( mess_matrix mat )
{
    MSG_FNAME(__func__);

    mess_int_t i=0 , j=0, end=0, idx=0, nnz=0;
    double val = 0;
    mess_double_cpx_t val_cpx = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(mat);
    mess_check_real_or_complex(mat);

    if ( MESS_IS_DENSE(mat)) {
        return 0;
    }

    if ( MESS_IS_CSC(mat)) {
        /*-----------------------------------------------------------------------------
         *  CSC variant
         *-----------------------------------------------------------------------------*/
        if(MESS_IS_REAL(mat)){
            /*-----------------------------------------------------------------------------
             *  MESS_REAL
             *-----------------------------------------------------------------------------*/
            for(i=0; i<mat->cols; ++i){
                j = end;
                end = mat->colptr[i+1];
                while(j<end){
                    idx = mat->rowptr[j];
                    val = mat->values[j];
                    if(val!=0.0){
                        mat->rowptr[nnz]=idx;
                        mat->values[nnz]=val;
                        ++nnz;
                    }
                    ++j;
                }
                mat->colptr[i+1]=nnz;
            }
            mat->nnz = nnz;
            mess_try_realloc(mat->values, double*, sizeof(double)*nnz);
            mess_try_realloc(mat->rowptr, mess_int_t*, sizeof(mess_int_t)*nnz);
        }else{
            /*-----------------------------------------------------------------------------
             *  MESS_COMPLEX
             *-----------------------------------------------------------------------------*/
            for(i=0; i<mat->cols; ++i){
                j = end;
                end = mat->colptr[i+1];
                while(j<end){
                    idx = mat->rowptr[j];
                    val_cpx = mat->values_cpx[j];
                    if(val_cpx!=0.0){
                        mat->rowptr[nnz]=idx;
                        mat->values_cpx[nnz]=val_cpx;
                        ++nnz;
                    }
                    ++j;
                }
                mat->colptr[i+1]=nnz;
            }
            mat->nnz = nnz;
            mess_try_realloc(mat->values_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*nnz);
            mess_try_realloc(mat->rowptr, mess_int_t*, sizeof(mess_int_t)*nnz);
        }
    } else if ( MESS_IS_CSR(mat)){
        /*-----------------------------------------------------------------------------
         *  CSR  variant
         *-----------------------------------------------------------------------------*/
        if(MESS_IS_REAL(mat)){
            /*-----------------------------------------------------------------------------
             *  MESS_REAL
             *-----------------------------------------------------------------------------*/
            for(i=0; i<mat->rows; ++i){
                j = end;
                end = mat->rowptr[i+1];
                while(j<end){
                    idx = mat->colptr[j];
                    val = mat->values[j];
                    if(val!=0.0){
                        mat->colptr[nnz]=idx;
                        mat->values[nnz]=val;
                        ++nnz;
                    }
                    ++j;
                }
                mat->rowptr[i+1]=nnz;
            }
            mat->nnz = nnz;
            mess_try_realloc(mat->values, double*, sizeof(double)*nnz);
            mess_try_realloc(mat->colptr, mess_int_t*, sizeof(mess_int_t)*nnz);
        }else{
            /*-----------------------------------------------------------------------------
             *  MESS_COMPLEX
             *-----------------------------------------------------------------------------*/
            for(i=0; i<mat->rows; ++i){
                j = end;
                end = mat->rowptr[i+1];
                while(j<end){
                    idx = mat->colptr[j];
                    val_cpx = mat->values_cpx[j];
                    if(val_cpx!=0.0){
                        mat->colptr[nnz]=idx;
                        mat->values_cpx[nnz]=val_cpx;
                        ++nnz;
                    }
                    ++j;
                }
                mat->rowptr[i+1]=nnz;
            }
            mat->nnz = nnz;
            mess_try_realloc(mat->values_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*nnz);
            mess_try_realloc(mat->colptr, mess_int_t*, sizeof(mess_int_t)*nnz);
        }
    } else if (MESS_IS_COORD(mat)){
        if(MESS_IS_REAL(mat)){
            /*-----------------------------------------------------------------------------
             *  MESS_REAL
             *-----------------------------------------------------------------------------*/
            for(j=0;j<mat->nnz;++j){
                if(mat->values[j]!=0.0){
                    mat->rowptr[nnz]=mat->rowptr[j];
                    mat->colptr[nnz]=mat->colptr[j];
                    mat->values[nnz]=mat->values[j];
                    ++nnz;
                }
            }
            mat->nnz = nnz;
            mess_try_realloc(mat->rowptr, mess_int_t*, sizeof(mess_int_t)*nnz);
            mess_try_realloc(mat->colptr, mess_int_t*, sizeof(mess_int_t)*nnz);
            mess_try_realloc(mat->values, double*, sizeof(double)*nnz);
        }else{
            /*-----------------------------------------------------------------------------
             *  MESS_COMPLEX
             *-----------------------------------------------------------------------------*/
            for(j=0;j<mat->nnz;++j){
                if(mat->values_cpx[j]!=0.0){
                    mat->rowptr[nnz]=mat->rowptr[j];
                    mat->colptr[nnz]=mat->colptr[j];
                    mat->values_cpx[nnz]=mat->values_cpx[j];
                    ++nnz;
                }
            }
            mat->nnz = nnz;
            mess_try_realloc(mat->rowptr, mess_int_t*, sizeof(mess_int_t)*nnz);
            mess_try_realloc(mat->colptr, mess_int_t*, sizeof(mess_int_t)*nnz);
            mess_try_realloc(mat->values_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*nnz);
        }
    }else{
        MSG_ERROR("Unkown storage type.\n");
        return MESS_ERROR_STORAGETYPE;
    }
    return 0;
}       /* -----  end of function mess_matrix_eliminiate_zeros  ----- */




