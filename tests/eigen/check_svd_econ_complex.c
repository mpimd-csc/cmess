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
 * @addtogroup test_eigen
 * @{
 * @file tests/eigen/check_svd_econ_complex.c
 * @brief Check the economy size singular value decomposition (complex version).
 * @author @koehlerm
 * @test
 * This function checks the @ref mess_eigen_svd_econ function defined in svd.c that means it checks if the economy size
 * singular value decomposition
 * \f[ A =  U S V^T\f]
 * is computed correctly.
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

int main ( int argc, char ** argv) {
    mess_init();
    mess_matrix A,U1,V1,U2,V2,tmp,tmp2, eye;
    mess_vector sigma_ist, sigma_soll;
    mess_int_t i,n,m ;
    int err = 0;
    int ret = 0;
    double r, nrm;

    if (argc!=3) {
        fprintf(stderr, "usage: %s rows cols\n",argv[0]);
        return -1;
    }
    n = atoi(argv[1]);
    m = atoi(argv[2]);

    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_init(&U1));
    CALL(mess_matrix_init(&U2));
    CALL(mess_matrix_init(&V1));
    CALL(mess_matrix_init(&V2));
    CALL(mess_matrix_init(&tmp));
    CALL(mess_matrix_init(&tmp2));
    CALL(mess_matrix_init(&eye));
    MESS_INIT_VECTORS(&sigma_ist,&sigma_soll);
    CALL(mess_vector_alloc(sigma_ist,m,MESS_REAL));
    CALL(mess_vector_alloc(sigma_soll,m,MESS_REAL));
    CALL(mess_matrix_eyec(eye,n,n,MESS_DENSE));

    CALL(mess_matrix_rand_dense(tmp, n, m,MESS_REAL));
    CALL(mess_matrix_rand_dense(tmp2, n, m,MESS_REAL));
    CALL(mess_matrix_addc(3,tmp2,I,tmp));
    CALL(mess_matrix_orth(tmp,U1));
    CALL(mess_matrix_rand_dense(tmp, m, m,MESS_REAL));
    CALL(mess_matrix_rand_dense(tmp2, m,m,MESS_REAL ));
    CALL(mess_matrix_addc(3,tmp2,-I,tmp));
    CALL(mess_matrix_orth(tmp,V1));

    for ( i = 0; i < m;  i++) {
        sigma_soll->values[i] = m-i;
        CALL(mess_matrix_colscale(U1, i, (double) m-i));
    }
    CALL( mess_matrix_multiply(MESS_OP_NONE, U1, MESS_OP_HERMITIAN, V1, A));

    CALL(mess_eigen_svd_econ(A,sigma_ist, U2, V2));

    CALL(mess_vector_sort(sigma_ist));
    CALL(mess_vector_sort(sigma_soll));
    CALL(mess_vector_norm2(sigma_soll, &nrm));
    CALL(mess_vector_diffnorm(sigma_soll, sigma_ist, &r));
    printf("diff = %lg\n", r /nrm);
    if ( r > nrm * 1000 * mess_eps()) err++;

    CALL( mess_matrix_multiply(MESS_OP_HERMITIAN, U2, MESS_OP_NONE, U2, tmp));
    CALL(mess_matrix_eyec(eye, tmp->rows, tmp->cols,MESS_DENSE));
    CALL(mess_matrix_diffnorm(tmp, eye, &r));
    if ( r > n * 100 * mess_eps()) err++;
    if ( tmp->rows != m ) err ++;
    if ( U2->cols != m ) err ++;

    printf("diffU2 = %lg\n", r );

    CALL( mess_matrix_multiply(MESS_OP_HERMITIAN, V2, MESS_OP_NONE, V2, tmp));
    CALL(mess_matrix_eyec(eye, tmp->rows, tmp->cols,MESS_DENSE));
    CALL(mess_matrix_diffnorm(tmp, eye, &r));
    printf("diffV2 = %lg\n", r );
    if ( r > n * 100 * mess_eps()) err++;
    if ( tmp->rows != m ) err ++;
    if ( V2->rows != m ) err ++;

    printf("err = %d\n", err);


    mess_matrix_clear(&A);
    mess_matrix_clear(&V1);
    mess_matrix_clear(&V2);
    mess_matrix_clear(&U1);
    mess_matrix_clear(&U2);
    mess_matrix_clear(&eye);
    mess_matrix_clear(&tmp);
    mess_matrix_clear(&tmp2);
    mess_vector_clear(&sigma_ist);
    mess_vector_clear(&sigma_soll);
    return err;
}

