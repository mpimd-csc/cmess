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
 * @file tests/eigen/check_arnoldi_template_nrm.c
 * @brief Check the @ref mess_eigen_arnoldi_template_nrm function.
 * @author @mbehr
 * @test
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
    mess_error_level=3;
    mess_int_t ret = 0;
    mess_int_t k = 10;
    double diff = 0, tol = sqrt(mess_eps()), eigA, testeigA;

    mess_vector sv, svc, ewA;
    mess_matrix A, Ac, Tmp1;
    mess_mvpcall mvpA, mvpAc;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 3) {
        printf("usage: %s A.mtx k \n", argv[0]);
        return 1;
    }

    k = atoi(argv[2]);

    /*-----------------------------------------------------------------------------
     *  init matrices and vectors
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&Ac,&Tmp1);
    MESS_INIT_VECTORS(&sv,&svc, &ewA);

    CALL(mess_matrix_read_formated(argv[1], A, MESS_CSC));
    CALL(mess_mvpcall_matrix(&mvpA, MESS_OP_NONE, A));

    CALL(mess_matrix_copy(A,Ac));
    CALL(mess_matrix_scalec(1+2*I,Ac));
    CALL(mess_mvpcall_matrix(&mvpAc, MESS_OP_NONE, Ac));

    CALL(mess_vector_alloc(sv, A->rows, MESS_REAL));
    CALL(mess_vector_rand(sv));

    CALL(mess_vector_alloc(svc, A->rows, MESS_COMPLEX));
    CALL(mess_vector_rand(svc));

    /*-----------------------------------------------------------------------------
     * call mess_arnoldi_template_nrm and check results
     *-----------------------------------------------------------------------------*/
    // real variant
    CALL(mess_eigen_arnoldi_template_nrm(mvpA, k, sv, &eigA));

    CALL(mess_matrix_convert(A, Tmp1, MESS_DENSE));
    CALL(mess_eigen_eig(Tmp1, ewA, NULL));
    CALL(mess_vector_norminf(ewA, &testeigA));

    diff = fabs(eigA-testeigA);
    printf("diff = %e\n", diff);
    if(diff >= tol){
        printf("Failed:\n");
        printf("Got         = %e\n", eigA);
        printf("Expected    = %e\n", testeigA);
        ret = 1;
        goto clear;
    }

    //complex variant
    CALL(mess_eigen_arnoldi_template_nrm(mvpAc, k, svc, &eigA));

    CALL(mess_matrix_convert(Ac, Tmp1, MESS_DENSE));
    CALL(mess_eigen_eig(Tmp1, ewA, NULL));
    CALL(mess_vector_norminf(ewA, &testeigA));

    diff = fabs(eigA-testeigA);
    printf("diff = %e\n", diff);
    if(diff >= tol){
        printf("Failed:\n");
        printf("Got         = %e\n", eigA);
        printf("Expected    = %e\n", testeigA);
        ret = 1;
        goto clear;
    }


    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
clear:
    MESS_CLEAR_MATRICES(&A,&Ac,&Tmp1);
    MESS_CLEAR_VECTORS(&sv,&svc,&ewA);
    CALL(mess_mvpcall_clear(&mvpA));
    CALL(mess_mvpcall_clear(&mvpAc));

    return ret;
}

