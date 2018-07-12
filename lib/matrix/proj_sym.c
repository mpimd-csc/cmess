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
 * @file lib/matrix/proj_sym.c
 * @brief Projection of a matrix on symmetric matrix space.
 * @author @dykstra
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"

#include <complex.h>


/**
 * @brief Get projection of a matrix on symmetric matrix space.
 * @param[in,out] matrix    input matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_proj_sym function projects the given matrix on the space of symmetric matrices
 * by computing \f$A \mapsto \mathcal{P}(A):=frac{1}{2}(A+A^T)\f$, where \f$\mathcal(P):\mathbb(R)^{n \times n}\to \mathbb(S)^{n \times n}\f$.
 *
 */
int mess_matrix_proj_sym(mess_matrix matrix){
    MSG_FNAME(__func__);
    mess_int_t i,j;
    int ret = 0 ;

    mess_check_nullpointer(matrix);
    mess_check_square(matrix);


    if(MESS_IS_DENSE(matrix)){
        if(MESS_IS_REAL(matrix)){
            for(j = 0; j < matrix->cols; j++){
                for(i = j+1; i < matrix->cols; i++){
                    matrix->values[i+j*matrix->ld] = 0.5*(matrix->values[i+j*matrix->ld] + matrix->values[j+i*matrix->ld]);
                }
            }
            for(j = 0; j < matrix->cols; j++){
                for(i = 0; i < j; i++){
                    matrix->values[i+j*matrix->ld] = matrix->values[j+i*matrix->ld];
                }
            }
        } else if(MESS_IS_COMPLEX(matrix)){
            for(j = 0; j < matrix->cols; j++){
                for(i = j+1; i < matrix->cols; i++){
                    matrix->values_cpx[i+j*matrix->ld] = 0.5*(matrix->values_cpx[i+j*matrix->ld] + matrix->values_cpx[j+i*matrix->ld]);
                }
            }
            for(j = 0; j < matrix->cols; j++){
                for(i = 0; i < j; i++){
                    matrix->values_cpx[i+j*matrix->ld] = matrix->values_cpx[j+i*matrix->ld];
                }
            }
        } else {
            MSG_ERROR("Unsupported Data type: %s \n", mess_datatype_t_str(matrix->data_type));
            return(MESS_ERROR_DATATYPE);
        }
    } else {
        MSG_WARN("Computing the symmetric projection of a sparse matrix may lead to big fill-in.\n");
        mess_matrix trans;
        ret = mess_matrix_init(&trans);
        ret = mess_matrix_transpose(matrix,trans);      FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_transpose);
        ret = mess_matrix_add(0.5,trans,0.5,matrix);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_add);
        mess_matrix_clear(&trans);
    }


    return(0);
}

