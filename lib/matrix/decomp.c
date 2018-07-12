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
 * @file lib/matrix/decomp.c
 * @brief Decompose a Matrix in symmetric (hermitian) and skewsymmetric (skewhermitian) part.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>

/**
 * @brief Decompose a matrix in a symmetric (hermitian) and skewsymmetric (skewhermitian) part.
 * @param[in]   A           input matrix
 * @param[out]  Asym        output symmetric (hermitian) part of \f$ A\f$ (@c NULL if not wanted)
 * @param[out]  Askewsym    output skewsymmetric (skewhermitian) oart of \f$ A\f$ (@c NULL if not wanted)
 * @return zero on success or a non-zero error value otherwise

 * The @ref mess_matrix_decomp function decomposes a matrix \f$A\f$ into
 * \f[
 * \begin{aligned}
 *  A_{sym}\leftarrow \frac{1}{2}(A+A^H)  \\
 *  A_{skewsym} \leftarrow \frac{1}{2}(A-A^H)
 *  \end{aligned}.
 * \f]
 *  Note that @p Asym and @p Aksewsym have the same storage type as @p A.
 */

int mess_matrix_decomp ( mess_matrix A, mess_matrix Asym, mess_matrix Askewsym){
    MSG_FNAME(__func__);
    mess_int_t ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_real_or_complex(A);
    mess_check_square(A)


        /*-----------------------------------------------------------------------------
         * compute hermitian and skewhermitian part
         *-----------------------------------------------------------------------------*/
        if(Asym==NULL && Askewsym==NULL){
            return 0;
        }

    if(Asym){
        MESS_MATRIX_RESET(Asym);
        ret = mess_matrix_ctranspose(A,Asym);                        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
        ret = mess_matrix_add(0.5,A,0.5,Asym);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
    }

    if(Askewsym){
        MESS_MATRIX_RESET(Askewsym);
        ret = mess_matrix_ctranspose(A,Askewsym);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
        ret = mess_matrix_add(0.5,A,-0.5,Askewsym);                 FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_add);
    }

    return ret ;
}




















