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
 * @file lib/matrix/lift.c
 * @brief Add a lower trailing block zeros to a matrix.
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>



/**
 * @brief Add lower trailing block of n rows to matrix.
 * @param[in]       in input matrix
 * @param[in]        n input positive number of rows of zeros to add
 * @param[out]     out matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_lift function adds a trailing block of @c n rows zeros at
 * the lower of @ref mess_matrix @c in. The result is in @ref mess_matrix @c out.
 *
 * @attention The function only works for dense matrices.
 *
 * @see mess_matrix_cat
 * @see mess_matrix_zeros
 *
 */
int mess_matrix_lift(mess_matrix in, mess_int_t n,  mess_matrix out ){
    MSG_FNAME(__func__);
    mess_int_t ret=0;

    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_positive(n);
    mess_check_nullpointer(in);
    mess_check_real_or_complex(in);
    mess_check_dense(in);
    mess_check_nullpointer(out);

    /*-----------------------------------------------------------------------------
     *  perform [in;zeros(n,in->cols)]->out
     *-----------------------------------------------------------------------------*/
    ret = mess_matrix_alloc(out,n+in->rows,in->cols,(n+in->rows)*(in->cols),MESS_DENSE,in->data_type);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
    if(MESS_IS_REAL(in)){
        F77_GLOBAL(dlacpy,DLACPY)("A", &in->rows, &in->cols, in->values, &in->ld, out->values, &out->ld);
    }else{
        F77_GLOBAL(zlacpy,ZLACPY)("A", &in->rows, &in->cols, in->values_cpx, &in->ld, out->values_cpx, &out->ld);
    }

    return ret;
}


