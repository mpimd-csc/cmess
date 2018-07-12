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
 * @file tests/eigen/check_eig_hessenberg_schur_largest_ev.c
 * @brief Check the @ref mess_eigen_hessenberg_abs_largest.
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



int main(int argc, char ** argv){
    mess_init();
    mess_error_level = 3;
    mess_int_t ret = 0, seed = 1, dim = 50, k_triu=-1, maxind;
    double diff = 0, tol = sqrt(mess_eps()), dummy;
    mess_double_cpx_t largest_ev, largest_ev_test;
    mess_datatype_t dt = MESS_REAL;
    mess_vector ev, evector, tmp;
    mess_matrix hessA;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    if ( argc != 4 ) {
        fprintf (stderr, "Usage: %s dim cpx k_triu\n", argv[0]) ;
        exit(1);
    }

    dim     = atoi(argv[1]);
    dt      = atoi(argv[2])?MESS_COMPLEX:MESS_REAL;
    k_triu  = atoi(argv[3]);

    if(k_triu<-1){
        printf("ktriu must be larger or equal to -1, otherwise no Hessenberg matrix is generated\n");
        ret = 1;
        goto clear;
    }

    /*-----------------------------------------------------------------------------
     *  init matrices and vectors
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&hessA);
    MESS_INIT_VECTORS(&ev, &evector, &tmp);

    /*-----------------------------------------------------------------------------
     *  create upper hessenberg matrix
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_rand_init(&seed));
    CALL(mess_matrix_rand(hessA, dim, dim, MESS_DENSE, dt, 1.0));
    CALL(mess_matrix_triu(hessA, k_triu));

    /*-----------------------------------------------------------------------------
     *  call mess_eigen_hessenberg_abs_largest and check and check function
     *-----------------------------------------------------------------------------*/
    CALL(mess_eigen_hessenberg_abs_largest(hessA, evector, &largest_ev));

    //check residual of eigenvalue and eigenvector
    CALL(mess_matrix_mvp(MESS_OP_NONE, hessA, evector, tmp));
    CALL(mess_vector_scalec(largest_ev, evector));
    CALL(mess_vector_diffnorm(evector, tmp, &diff));

    printf("|| Av - lambda v ||_2 = %e\n", diff);
    if(diff>=tol){
        printf("Failed:\n");
        ret = 1;
        goto clear;
    }

    //check eigenvalue/eigenvector property
    CALL(mess_eigen_eig(hessA, ev, NULL));
    CALL(mess_vector_max(ev, &dummy, &maxind));
    CALL(mess_vector_get(ev, maxind, &largest_ev_test));

    //check if largest eigenvalue was computed
    diff = MESS_MIN(cabs(largest_ev-largest_ev_test),cabs(largest_ev-conj(largest_ev_test)));
    printf("| largest_ev - largest_ev_test|  = %e\n", diff);
    if(diff>=tol){
        printf("Failed:\n");
        printf("Got                     = %e + %eI\n",creal(largest_ev),cimag(largest_ev));
        printf("Got abs. value          = %e\n",cabs(largest_ev));
        printf("Expected                = %e + %eI\n",creal(largest_ev_test),cimag(largest_ev_test));
        printf("Expected abs. value     = %e\n",cabs(largest_ev_test));
        mess_vector_print(ev);
        ret = 1;
        goto clear;
    }

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
clear:
    MESS_CLEAR_MATRICES(&hessA);
    MESS_CLEAR_VECTORS(&ev, &evector, &tmp);

    return ret;
}

