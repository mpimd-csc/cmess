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
 * @file tests/matrix/check_coldot.c
 * @brief Check the column dot product.
 * @author @koehlerm
 * @test
 * This function check the @ref mess_matrix_coldot function defined in colops.c, that means
 * it checks if the dot product
 * \f[ Q(:,col_1)^TQ(:,col_2)\f]
 * is computed correctly for a given matrix \f$ Q \f$ and columns \f$ col_1\f$ and \f$ col_2 \f$.
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mess/mess.h"

#include "../call_macro.h"
#define  CALL2(X) ret = (X); tests++; if ( ret == 1) { err ++; fprintf(stderr,"================  test %35s failed =================\n", #X);  } else {fprintf(stderr, "test %35s passed.\n",#X);  }


double matrix_data[16] = {11,0,31,0,12,22,0,42,0,0,33,43,14,24,0,0};

int main ( int argc, char ** argv) {
    mess_init();
    int ret;
    int err = 0;
    double eps = mess_eps();
    mess_matrix A,Acsr,Acsc;
    double dot1,dot2,dot3;

    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_init(&Acsc));
    CALL(mess_matrix_init(&Acsr));

    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(A,4,4,0,matrix_data,NULL));
    CALL(mess_matrix_convert(A,Acsr,MESS_CSR));
    CALL(mess_matrix_convert(A,Acsc,MESS_CSC));

    CALL(mess_matrix_coldot(A,1,3,&dot1));
    CALL(mess_matrix_coldot(Acsr,1,3,&dot2));
    CALL(mess_matrix_coldot(Acsc,1,3,&dot3));

    if ( fabs(dot1-dot2)/fabs(dot1) > 10 * eps) err++;
    if ( fabs(dot1-dot3)/fabs(dot1) > 10 * eps) err++;

    printf("err = %d \nnorm1=%lg\nnorm2=%lg\nnorm3=%lg\n",err,dot1,dot2,dot3);
    mess_matrix_clear(&A);
    mess_matrix_clear(&Acsr);
    mess_matrix_clear(&Acsc);
    return (err>0)?(1):(0);
}

