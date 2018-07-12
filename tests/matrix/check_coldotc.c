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
 * @file tests/matrix/check_coldotc.c
 * @brief Check the column dot product of a complex matrix.
 * @author @koehlerm
 * @test
 * This function check the @ref mess_matrix_coldotc function defined in colops.c, that means it checks if the dot product
 * \f[ Q(:,col_1)^HQ(:,col_2)\f]
 * is computed correctly for a given complex matrix \f$ Q \f$ and columns \f$ col_1\f$ and \f$ col_2 \f$.
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mess/mess.h"

#include "../call_macro.h"


mess_double_cpx_t matrix_data[16] = {11,0,31,0,12-3*I,22+I,0,42,0,0,33,43,14,24,0,0};

int main ( int argc, char ** argv) {
    mess_init();
    int ret;
    int err = 0;
    double eps = mess_eps();
    mess_matrix A,Acsr,Acsc;
    mess_double_cpx_t dot1,dot2,dot3;

    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_init(&Acsc));
    CALL(mess_matrix_init(&Acsr));

    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(A,4,4,0,NULL,matrix_data));
    CALL(mess_matrix_convert(A,Acsr,MESS_CSR));
    CALL(mess_matrix_convert(A,Acsc,MESS_CSC));

    CALL(mess_matrix_coldotc(A,1,3,&dot1));
    CALL(mess_matrix_coldotc(Acsr,1,3,&dot2));
    CALL(mess_matrix_coldotc(Acsc,1,3,&dot3));

    if ( cabs(dot1-dot2)/cabs(dot1) > 10 * eps) err++;
    if ( cabs(dot1-dot3)/cabs(dot1) > 10 * eps) err++;
    printf("dot1 = %lg %lg\n", creal(dot1),cimag(dot1));
    printf("dot2 = %lg %lg\n", creal(dot2),cimag(dot2));
    printf("dot3 = %lg %lg\n", creal(dot3),cimag(dot3));


    printf("err = %d \n",err);
    mess_matrix_clear(&A);
    mess_matrix_clear(&Acsr);
    mess_matrix_clear(&Acsc);
    return (err>0)?(1):(0);
}

