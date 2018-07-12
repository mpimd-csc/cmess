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
 * @file tests/matrix/check_diag.c
 * @brief Check the output of diagonal elements of a matrix.
 * @author @koehlerm
 * @test
 * This function checks the @ref mess_matrix_diag function defined in basic/diag.c that means it checks if the
 * output is equal to the main diagonal elements of a given matrix.
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
#define  CALL2(X) ret = (X); tests++; if ( ret == 1) { err ++; fprintf(stderr,"================  test %35s failed =================\n", #X);  } else {fprintf(stderr, "test %35s passed.\n",#X);  }


double matrix_data[16] = {11,0,31,0,12,22,0,42,0,0,33,43,14,24,0,0};
mess_double_cpx_t matrix_datac[16] = {11,0,31,0,12-3*I,22+I,0,42,0,0,33,43,14,24,0,0};


int main ( int argc, char ** argv) {
    mess_init();
    int ret;
    int err = 0;
    double eps = mess_eps();
    mess_matrix A,Acsr,Acsc;
    mess_matrix Ai,Acsri,Acsci;
    mess_vector d1,d2,d3;
    mess_vector di1,di2,di3;
    double nrm1,nrm2,nrm3,nrm4;


    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_init(&Acsc));
    CALL(mess_matrix_init(&Acsr));
    CALL(mess_matrix_init(&Ai));
    CALL(mess_matrix_init(&Acsci));
    CALL(mess_matrix_init(&Acsri));


    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(A,4,4,0,matrix_data,NULL));
    CALL(mess_matrix_convert(A,Acsr,MESS_CSR));
    CALL(mess_matrix_convert(A,Acsc,MESS_CSC));

    CALL(mess_matrix_dense_from_farray(Ai,4,4,0,NULL,matrix_datac));
    CALL(mess_matrix_convert(Ai,Acsri,MESS_CSR));
    CALL(mess_matrix_convert(Ai,Acsci,MESS_CSC));

    MESS_INIT_VECTORS(&d1,&d2,&d3,&di1,&di2,&di3);
    CALL(mess_vector_alloc(d1,A->rows,MESS_REAL));
    CALL(mess_vector_alloc(d2,A->rows,MESS_REAL));
    CALL(mess_vector_alloc(d3,A->rows,MESS_REAL));
    CALL(mess_vector_alloc(di1,A->rows,MESS_COMPLEX));
    CALL(mess_vector_alloc(di2,A->rows,MESS_COMPLEX));
    CALL(mess_vector_alloc(di3,A->rows,MESS_COMPLEX));

    CALL(mess_matrix_diag(A,d1));
    CALL(mess_matrix_diag(Acsr,d2));
    CALL(mess_matrix_diag(Acsc,d3));
    CALL(mess_matrix_diag(Ai,di1));
    CALL(mess_matrix_diag(Acsri,di2));
    CALL(mess_matrix_diag(Acsci,di3));

    CALL(mess_vector_diffnorm(d1,d2,&nrm1));
    CALL(mess_vector_diffnorm(d1,d3,&nrm2));
    CALL(mess_vector_diffnorm(di1,di2,&nrm3));
    CALL(mess_vector_diffnorm(di1,di3,&nrm4));

    if ( nrm1 > 10*eps) err++;
    if ( nrm2 > 10*eps) err++;
    if ( nrm3 > 10*eps) err++;
    if ( nrm4 > 10*eps) err++;


    printf("err = %d \n",err);
    mess_matrix_clear(&A);
    mess_matrix_clear(&Acsr);
    mess_matrix_clear(&Acsc);
    mess_matrix_clear(&Ai);
    mess_matrix_clear(&Acsri);
    mess_matrix_clear(&Acsci);
    mess_vector_clear(&d1);
    mess_vector_clear(&d2);
    mess_vector_clear(&d3);
    mess_vector_clear(&di1);
    mess_vector_clear(&di2);
    mess_vector_clear(&di3);

    return (err>0)?(1):(0);
}

