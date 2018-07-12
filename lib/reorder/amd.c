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
 * @file lib/reorder/amd.c
 * @brief Compute AMD reordering of a matrix.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>

#ifdef MESS_HAVE_AMD
#include <amd.h>
#endif

#ifdef MESS_HAVE_COLAMD
#include <colamd.h>
#endif

#ifdef MESS64
#define AMD_VALID amd_l_valid
#define AMD_ORDER amd_l_order
#define AMD_PINFO  amd_l_info
#else
#define AMD_VALID amd_valid
#define AMD_ORDER amd_order
#define AMD_PINFO  amd_info
#endif


/**
 * @brief Compute @amd reordering of a matrix.
 * @param[in]  A        input matrix \f$A\f$
 * @param[out] p        output permutation array of length  rows of \f$ A\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_reorder_amd function reorders a matrix with approximate minimum degree (@amd) ordering.\n
 * It uses the @amd reordering subroutines from @suitesparse.
 *
 * @sa mess_matrix_reorder_rcm
 * @sa mess_matrix_reorder_colamd
 * @sa mess_matrix_reorder
 */
int mess_matrix_reorder_amd(mess_matrix A, mess_int_t *p){
    MSG_FNAME(__func__);
#ifdef MESS_HAVE_AMD
    mess_matrix csc;
    int conv = -1;
    int ret = 0;
    double Info [AMD_INFO];

    mess_check_nullpointer(A);
    mess_check_nullpointer(p);

    if ( MESS_IS_CSR(A)) {
        ret = AMD_VALID(A->rows, A->cols, A->rowptr, A->colptr) ;
        if ( ret != AMD_OK && ret != AMD_OK_BUT_JUMBLED) {
            MSG_ERROR("matrix is invalid as input of amd_order\n");
            return MESS_ERROR_GENERAL;
        }
        ret = AMD_ORDER (A->rows, A->rowptr, A->colptr, p, (double*)NULL, Info );
    } else if ( MESS_IS_CSC(A)){
        ret = AMD_VALID(A->rows, A->cols, A->colptr, A->rowptr) ;
        if ( ret != AMD_OK && ret != AMD_OK_BUT_JUMBLED) {
            MSG_ERROR("matrix is invalid as input of amd_order\n");
            return MESS_ERROR_GENERAL;
        }
        ret = AMD_ORDER(A->cols, A->colptr, A->rowptr, p, (double*)NULL, Info );
    } else {
        MESS_MATRIX_CHECKFORMAT(A, csc, conv, MESS_CSC);
        ret = AMD_ORDER(csc->cols, csc->colptr, csc->rowptr, p, (double*)NULL, Info );      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),AMD_ORDER);
        if (conv == 0){
            mess_matrix_clear(&csc);
        }
    }
    if ( mess_error_level > 2)
        AMD_PINFO(Info);

    return 0;
#else
    MSG_ERROR("no AMD support available.\n");
    return MESS_ERROR_NOSUPPORT;

#endif
}


#ifdef MESS64
#define COLAMD_RECOMMENDED colamd_l_recommended
#define COLAMD_SET_DEFAULTS colamd_l_set_defaults
#define COLAMD colamd_l
#define COLAMD_REPORT colamd_l_report
#else
#define COLAMD_RECOMMENDED colamd_recommended
#define COLAMD_SET_DEFAULTS colamd_set_defaults
#define COLAMD colamd
#define COLAMD_REPORT colamd_report

#endif

/**
 * @brief Compute the @colamd reordering of a matrix.
 * @param[in] A     input matrix \f$A\f$
 * @param[out] p    input/output permutation of length columns of  \f$ A \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_reorder_colamd function computes the column approximate minimum degree (@colamd) reordering of a given
 * matrix. \n
 * It uses the @colamd package from @suitesparse.
 *
 * @sa mess_matrix_reorder_rcm
 * @sa mess_matrix_reorder_amd
 * @sa mess_matrix_reorder
 *
 */
int mess_matrix_reorder_colamd ( mess_matrix A, mess_int_t *p  )
{
    MSG_FNAME(__func__);
#ifdef MESS_HAVE_COLAMD
    size_t alen = 0;
    double knobs [COLAMD_KNOBS] ;
    mess_int_t stats [COLAMD_STATS];
    mess_int_t *Aval;
    mess_int_t *pi;
    int ret = 0 ;

    mess_check_nullpointer(A);
    mess_check_nullpointer(p);

    COLAMD_SET_DEFAULTS(knobs) ;
    alen = COLAMD_RECOMMENDED (A->nnz, A->rows, A->cols) ;

    if (MESS_IS_CSR(A)) {
        mess_try_alloc(Aval, mess_int_t*, sizeof(mess_int_t )*alen);
        memcpy(Aval, A->colptr, sizeof(mess_int_t) *A->nnz);
        mess_try_alloc(pi, mess_int_t *, sizeof(mess_int_t)*(A->rows+1));
        memcpy(pi, A->rowptr, sizeof(mess_int_t)*(A->rows+1));
        ret = COLAMD (A->rows, A->cols, alen, Aval, pi , knobs ,stats) ;
        if ( ret < 0 ) {
            mess_free(Aval);
            mess_free(pi);
            MSG_ERROR("colamd returned with an error ( %d )\n", ret );
            COLAMD_REPORT (stats) ;
            return MESS_ERROR_GENERAL;
        }
        memcpy(p, pi, sizeof(mess_int_t)*A->rows);
        mess_free(pi);
        mess_free(Aval);
    } else if ( MESS_IS_CSC (A)) {
        mess_try_alloc(Aval, mess_int_t*, sizeof(mess_int_t )*alen);
        memcpy(Aval, A->rowptr, sizeof(mess_int_t) *A->nnz);
        mess_try_alloc(pi, mess_int_t *, sizeof(mess_int_t)*(A->rows+1));
        memcpy(pi, A->colptr, sizeof(mess_int_t)*(A->rows+1));
        ret = COLAMD (A->rows, A->cols, alen, Aval, pi , knobs ,stats) ;
        if ( ret < 0 ) {
            mess_free(Aval);
            mess_free(pi);
            MSG_ERROR("colamd returned with an error ( %d )\n", ret );
            COLAMD_REPORT (stats) ;
            return MESS_ERROR_GENERAL;
        }
        memcpy(p, pi, sizeof(mess_int_t)*A->rows);
        mess_free(pi);
        mess_free(Aval);
    } else {
        MSG_ERROR("datatype not supported.\n");
        return MESS_ERROR_DATATYPE;
    }
    if ( mess_error_level > 2)
        COLAMD_REPORT(stats);
    return 0;
#else
    MSG_ERROR("no COLAMD support available.\n");
    return MESS_ERROR_NOSUPPORT;
#endif
}       /* -----  end of function mess_matrix_reorder_colamd  ----- */

