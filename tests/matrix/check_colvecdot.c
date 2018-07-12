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
 * @file tests/matrix/check_colvecdot.c
 * @brief Check the column vector product.
 * @author @koehlerm
 * @test
 * This function check the @ref mess_matrix_colvecdot function defined in colops.c that means it checks if the
 * vector product
 * \f[Q(:,col)^Tv \f]
 * is computed correctly for a given matrix \f$ Q \f$, a vector \f$ v \f$ and a column \f$ col \f$.
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


double matrix_data[16] = {11,0,31,0,12,22,0,42,0,0,33,43,14,24,0,0};

int main ( int argc, char ** argv) {
    mess_init();
    int ret;
    int err = 0;
    double eps = mess_eps();
    mess_matrix A,Acsr,Acsc;
    mess_vector v;
    double dot1=0,dot2=0,dot3=0;

    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_init(&Acsc));
    CALL(mess_matrix_init(&Acsr));


    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(A,4,4,0,matrix_data,NULL));
    CALL(mess_matrix_convert(A,Acsr,MESS_CSR));
    CALL(mess_matrix_convert(A,Acsc,MESS_CSC));
    MESS_INIT_VECTORS(&v);
    CALL(mess_vector_alloc(v,A->rows,MESS_REAL));
    CALL(mess_vector_ones(v));

    CALL(mess_matrix_colvecdot(A,3,v,&dot1));
    CALL(mess_matrix_colvecdot(Acsr,3,v,&dot2));
    CALL(mess_matrix_colvecdot(Acsc,3,v,&dot3));

    if ( fabs(dot1-dot2)/fabs(dot1) > 10 * eps) err++;
    if ( fabs(dot1-dot3)/fabs(dot1) > 10 * eps) err++;
    printf("dot1 = %lg\n", dot1);
    printf("dot2 = %lg\n", dot2);
    printf("dot3 = %lg\n", dot3);

    printf("err = %d \n",err);
    mess_matrix_clear(&A);
    mess_matrix_clear(&Acsr);
    mess_matrix_clear(&Acsc);
    mess_vector_clear(&v);
    return (err>0)?(1):(0);
}

