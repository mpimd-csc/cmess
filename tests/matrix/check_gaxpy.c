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
 * @addtogroup test_matrix
 * @{
 * @file tests/matrix/check_gaxpy.c
 * @brief Check the matrix-vector product with an update of a vector.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_gaxpy function defined in gaxpy.c that means it checks if the matrix-vector
 * product with an update of a vector \f$ y \f$
 * \f[ y \leftarrow y + op(A)x\f]
 * is computed correctly. \n
 * Operations \f$ op(.) \f$ on matrix \f$ A \f$ can be
 *  * \ref MESS_OP_NONE (\f$ op(A)= A \f$),
 *  * \ref MESS_OP_TRANSPOSE (\f$ op(A)= A^T \f$),
 *  * \ref MESS_OP_HERMITIAN (\f$ op(A)= A^H \f$).
 *
 *
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mess/mess.h"

#include "../call_macro.h"

/*
   ADENSE                   - dense input matrix for comparision
   TEMP                     - arbitrary matrix
   STORETYPEMAT         - storetype of the matrix for comparision
   OP                       - operation for gaxpy
   VIN                      - input vector for gaxpy (x)
   VTEMP1, VTEMP2, VTEMP3   - arbitrary vectors
   DATATYPEOUT              - datatype of output vector (y)
   DIFF                 - double diff between dense solution and computed
   ERR                      - mess_int_t inkrements if error occurs
   EPS                      - machine epsilon precision bound for DIFF
   */
#define CHECKGAXPY(ADENSE, TEMP, STORETYPEMAT, OP, VIN, VTEMP1, VTEMP2, VTEMP3, DATATYPEOUT, DIFF, ERR, EPS){                                   \
    CALL(mess_vector_totype(VTEMP2,DATATYPEOUT));                                                                                           \
    CALL(mess_vector_totype(VTEMP3,DATATYPEOUT));                                                                                           \
    CALL(mess_vector_copy(VIN,VTEMP1));                                                                                                     \
    CALL(mess_matrix_convert(ADENSE, TEMP, STORETYPEMAT));                                                                                  \
    CALL(mess_matrix_gaxpy(OP, ADENSE, VTEMP1, VTEMP2)); /*solution*/                                                                       \
    CALL(mess_vector_copy(VIN,VTEMP1));                                                                                                     \
    CALL(mess_matrix_gaxpy(OP,TEMP,VTEMP1,VTEMP3));                                                                                         \
    /* @ref mess_vector_print(VTEMP2); mess_vector_print(VTEMP3); */                                                                            \
    CALL(mess_vector_diffnorm(VTEMP2,VTEMP3,&DIFF));                                                                                        \
    if(DIFF>10*EPS){++ERR;                                                                                                                  \
        printf("Diff: %e\n",DIFF);                                                                                                              \
    };/*reset data*/                                                                                                                        \
    CALL(mess_vector_zeros(VTEMP1));                                                                                                        \
    CALL(mess_vector_zeros(VTEMP2));                                                                                                        \
    CALL(mess_vector_zeros(VTEMP3));                                                                                                        \
}


int main ( int argc, char ** argv) {
    mess_init();
    mess_int_t dim=10;
    int ret;
    int err = 0;
    double eps = mess_eps(), diff;
    mess_matrix Adense, Atemp;
    mess_vector  vin, vtemp1, vtemp2, vtemp3;

    //init real matrices
    CALL(mess_matrix_init(&Adense));
    CALL(mess_matrix_init(&Atemp));

    //init vector
    MESS_INIT_VECTORS(&vtemp1,&vtemp2,&vtemp3,&vin);
    CALL(mess_vector_alloc(vtemp1, dim, MESS_REAL));
    CALL(mess_vector_alloc(vtemp2, dim, MESS_REAL));
    CALL(mess_vector_alloc(vtemp3, dim, MESS_REAL));
    CALL(mess_vector_alloc(vin, dim, MESS_REAL));


    CALL(mess_matrix_rand(Adense,dim,dim,MESS_DENSE,MESS_REAL,1));
    CALL(mess_vector_rand(vin));

    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_toreal(Adense));
    CALL(mess_vector_toreal(vin));
    //real matrix, real input vector real output vector
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);

    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);

    //real matrix, real input vector complex output vector
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);

    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);

    CALL(mess_vector_tocomplex(vin));

    //real matrix, complex input vector real output vector
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);

    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);

    //real matrix, complex input vector complex output vector
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);

    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);


    CALL(mess_matrix_tocomplex(Adense));
    CALL(mess_vector_toreal(vin));
    //real matrix, real input vector real output vector
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);

    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);

    //real matrix, real input vector complex output vector
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);

    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);

    CALL(mess_vector_tocomplex(vin));

    //real matrix, complex input vector real output vector
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);

    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_REAL, diff, err, eps);

    //real matrix, complex input vector complex output vector
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSC, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);

    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_NONE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_TRANSPOSE, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);
    CHECKGAXPY(Adense, Atemp, MESS_CSR, MESS_OP_HERMITIAN, vin, vtemp1, vtemp2, vtemp3, MESS_COMPLEX, diff, err, eps);

    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/

    mess_matrix_clear(&Adense);
    mess_matrix_clear(&Atemp);
    mess_vector_clear(&vin);
    mess_vector_clear(&vtemp1);
    mess_vector_clear(&vtemp2);
    mess_vector_clear(&vtemp3);

    return (err>0)?(1):(0);
}
