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
 * @file tests/eigen/check_eig.c
 * @brief Check the solution of the eigenvalue problem using xGEEV from @lapack.
 * @author @koehlerm
 * @test
 * This function checks the @ref mess_eigen_eig function defined in eig.c that means it checks if the eigenvalue problem
 * \f[ Ax = \lambda x \f]
 * is solved correctly .
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

/**
 * @brief Check if the eigenvalue problem is solved correctly.
 * @param[in] A     input matrix
 * @param[in] Z     input matrix containing eigenvectors
 * @param[in] ev    input eigenvalue
 * @param[in] col   input column of eigenvector matrix to check
 * @param[in] tol   tolerance for error
 * @return @p diff > @p tol
 *
 * The @ref check_eigen_vector function checks if the eigenvalue problem
 * \f[ AV(:,col)= ev V(:,col) \f]
 * is solved correctly for a given eigenvalue \f$ ev \f$ and an eigenvector \f$ Z(:,col) \f$. \n
 * If
 * \f[\Vert A V(:,col)- ev V(:,col) \Vert_2 > \sqrt(\epsilon) \f]
 * the error value is increased by \f$ 1 \f$. \n
 * \f$ n \f$ is the number of rows of \f$ Z \f$.
 *
 */



int check_eigen_vector(mess_matrix A, mess_matrix Z, double ev , mess_int_t col, double tol) {
    double diff;
    int ret = 0;
    mess_vector EV, x;
    MESS_INIT_VECTORS(&EV,&x);
    CALL(mess_vector_alloc(EV, Z->rows, MESS_REAL));
    CALL(mess_vector_alloc(x,  Z->rows, MESS_REAL));
    CALL(mess_matrix_getcol(Z, col, EV));
    CALL(mess_matrix_mvp(MESS_OP_NONE,A, EV, x));
    CALL(mess_vector_axpy(-ev,EV,x));
    CALL(mess_vector_norm2(x,&diff));
    MESS_CLEAR_VECTORS(&EV,&x);
    printf("tol = %e \t diff = %lg\n",tol, diff);
    return diff>tol;
}



int main ( int argc, char ** argv) {
    mess_init();
    int ret = 0;
    mess_int_t i, n;
    double r=0, r2=0, tol = sqrt(mess_eps());
    mess_vector ev_soll, ev_ist;
    mess_matrix A, Z, Zorth, Zorth2;

    /*-----------------------------------------------------------------------------
     *  get input arguments
     *-----------------------------------------------------------------------------*/
    if (argc != 2 ) {
        fprintf(stderr, "usage: %s dim \n", argv[0]);
        return -1;
    }
    n = atoi(argv[1]);

    /*-----------------------------------------------------------------------------
     *  init variables
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&ev_ist,&ev_soll);
    CALL(mess_vector_alloc(ev_ist, n, MESS_REAL));
    CALL(mess_vector_alloc(ev_soll, n, MESS_REAL));
    MESS_INIT_MATRICES(&A,&Z,&Zorth,&Zorth2);

    /*-----------------------------------------------------------------------------
     *  generate data
     *-----------------------------------------------------------------------------*/
    for (i = 0; i < n; i++) {
        ev_soll->values[i]= i +1 ;
    }
    CALL(mess_matrix_rand_init(0));
    CALL(mess_matrix_rand_dense(Z, n, n,MESS_REAL ));
    CALL(mess_matrix_orth(Z, Zorth));
    CALL(mess_matrix_copy(Zorth, Zorth2));
    for (i = 0; i < n; i++) {
        CALL(mess_matrix_colscale(Zorth, i, ev_soll->values[i]));
    }
    CALL( mess_matrix_multiply(MESS_OP_NONE, Zorth, MESS_OP_HERMITIAN, Zorth2, A));

    /*-----------------------------------------------------------------------------
     *  compute eigenvalues/eigenvectors and test
     *-----------------------------------------------------------------------------*/
    CALL(mess_eigen_eig(A, ev_ist, Z));
    for (i = 0; i < n; i++) {
        CALL(check_eigen_vector(A,Z,ev_ist->values[i], i, tol));
    }

    CALL(mess_vector_sort(ev_ist));
    CALL(mess_vector_norm2(ev_soll, &r2));
    CALL(mess_vector_diffnorm(ev_ist, ev_soll, &r));

    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_VECTORS(&ev_soll,&ev_ist);
    MESS_CLEAR_MATRICES(&A,&Z,&Zorth,&Zorth2);

    /*-----------------------------------------------------------------------------
     * return
     *-----------------------------------------------------------------------------*/
    if(r/r2>tol){
        printf("tol=%e \t diff=%lg\n",tol, r/r2);
        return 1;
    }
    return 0;
}

