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
 * @file lib/matrix/fill.c
 * @brief Fill up symmetry in coordinate matrices.
 *
 * @author @koehlerm
 *
 *
 */
#include    <stdlib.h>
#include    <stdio.h>
#include    <math.h>
#include    "mess/mess.h"
#include    "mess/error_macro.h"

/**
 * @brief Fill up the symmetry in a coordinate matrix.
 * @param[in,out] mat matrix to fill up
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_symfillup function fills up the symmetry in a given symmetric matrix.
 * It does not check if the not stored triangle exists. \n
 * It only works if the matrix is stored as coordinate one and the SYMMETRY flag is set in the data structure.\n
 * It is mostly called from the \ref mess_matrix_read functions to output a general matrices without
 * any stored symmetry.
 *
 * @sa mess_matrix_read
 *
 */
int mess_matrix_symfillup(mess_matrix mat) {
    MSG_FNAME(__func__);
    mess_int_t nnz = 0;
    mess_int_t oldnnz = 0;
    mess_int_t i = 0;

    mess_check_nullpointer(mat);
    mess_check_coord(mat);

    nnz = mat->nnz;
    oldnnz = nnz;

    if ( MESS_IS_COMPLEX(mat) ){
        mess_try_realloc( mat->values_cpx , mess_double_cpx_t *, nnz * 2 * sizeof ( mess_double_cpx_t ) );
        mess_try_realloc( mat->rowptr, mess_int_t *, nnz * 2 * sizeof (mess_int_t));
        mess_try_realloc( mat->colptr, mess_int_t *, nnz * 2 * sizeof (mess_int_t));
        for ( i = 0; i < oldnnz ; i++ ) {
            if ( mat->colptr[i] == mat->rowptr[i] ) continue;
            mat->rowptr[nnz] = mat->colptr[i];
            mat->colptr[nnz] = mat->rowptr[i];
            if(MESS_IS_SKEWSYMMETRIC(mat)){
                mat->values_cpx[nnz] = -conj(mat->values_cpx[i]);
            } else {
                mat->values_cpx[nnz] = mat->values_cpx[i];
            }
            nnz++;
        }
        mess_try_realloc( mat->values_cpx , mess_double_cpx_t *, nnz * sizeof ( mess_double_cpx_t ) );
        mess_try_realloc( mat->rowptr, mess_int_t *, nnz  * sizeof (mess_int_t));
        mess_try_realloc( mat->colptr, mess_int_t *, nnz  * sizeof (mess_int_t));
        mat->nnz = nnz;
        mat->symmetry = MESS_GENERAL;
    } else if ( MESS_IS_REAL ( mat ) ) {
        mess_try_realloc( mat->values , double *, nnz * 2 * sizeof ( double  ) );
        mess_try_realloc( mat->rowptr, mess_int_t *, nnz * 2 * sizeof (mess_int_t));
        mess_try_realloc( mat->colptr, mess_int_t *, nnz * 2 * sizeof (mess_int_t));
        for ( i = 0; i < oldnnz ; i++ ) {
            if ( mat->colptr[i] == mat->rowptr[i] ) continue;
            mat->rowptr[nnz] = mat->colptr[i];
            mat->colptr[nnz] = mat->rowptr[i];
            if ( MESS_IS_SKEWSYMMETRIC(mat)){
                mat->values[nnz] = -mat->values[i];
            } else {
                mat->values[nnz] = mat->values[i];
            }
            nnz++;
        }
        mess_try_realloc( mat->values , double *, nnz * sizeof ( double  ) );
        mess_try_realloc( mat->rowptr, mess_int_t *, nnz  * sizeof (mess_int_t));
        mess_try_realloc( mat->colptr, mess_int_t *, nnz  * sizeof (mess_int_t));
        mat->nnz = nnz;
        mat->symmetry = MESS_GENERAL;
    } else {
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}

