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
 * @file tests/eigen/check_eig2.c
 * @brief Check the solution of the eigenvalue problem using xGEEV from @lapack (complex version).
 * @author @koehlerm
 * @test
 *
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
 * @param[in] A         input matrix
 * @param[in] Z         input matrix containing eigenvectors
 * @param[in] ev        input eigenvalue
 * @param[in] col       input column of eigenvector matrix to check
 * @param[in] tol       tolerance for error
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
int check_eigen_vector(mess_matrix A, mess_matrix Z, mess_double_cpx_t ev, mess_int_t col, double tol) {
    int ret = 0;
    double diff;
    mess_vector EV, x;
    MESS_INIT_VECTORS(&EV,&x);
    CALL(mess_vector_alloc(EV, Z->rows, MESS_REAL));
    CALL(mess_vector_alloc(x,  Z->rows, MESS_REAL));
    CALL(mess_matrix_getcol(Z, col, EV));
    CALL(mess_matrix_mvp(MESS_OP_NONE,A, EV, x));
    CALL(mess_vector_axpyc(-ev,EV,x));
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
    mess_vector ev_soll, ev_ist, h;
    mess_matrix A, Z, Zorth, Zorth2;

    /*-----------------------------------------------------------------------------
     *  get input arguments
     *-----------------------------------------------------------------------------*/
    if (argc != 2 ) {
        fprintf(stderr, "usage: %s dim \n", argv[0]);
        return -1;
    }
    n = atoi(argv[1]);
    n += n%2;

    /*-----------------------------------------------------------------------------
     *  init variables
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&ev_ist);
    MESS_INIT_VECTORS(&ev_soll);
    CALL(mess_vector_alloc(ev_ist, n-1, MESS_COMPLEX));
    CALL(mess_vector_alloc(ev_soll, n, MESS_COMPLEX));
    MESS_INIT_MATRICES(&A,&Z,&Zorth,&Zorth2);

    /*-----------------------------------------------------------------------------
     *  generate data
     *-----------------------------------------------------------------------------*/
    for (i = 0; i < n; i=i+2) {
        ev_soll->values_cpx[i]= i+1 + (n-i) * I ;
        ev_soll->values_cpx[i+1]= i+1 - (n-i) * I ;
    }

    /*-----------------------------------------------------------------------------
     *  Assemble eigenvector matrix
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_rand_init(NULL));
    CALL(mess_matrix_rand_dense(Zorth, n, n, MESS_REAL));
    CALL(mess_matrix_orth(Zorth,Z));
    CALL(mess_matrix_rand_dense(Zorth, n, n,MESS_REAL));
    CALL(mess_matrix_orth(Zorth,Zorth2));
    CALL(mess_matrix_addc(1,Zorth2,I, Z));
    MESS_INIT_VECTORS(&h);
    CALL(mess_vector_alloc(h, n, MESS_COMPLEX));
    for (i = 0; i < n; i=i+2) {
        Z->values_cpx[i+i*Z->ld] = 1;
        CALL(mess_matrix_getcol(Z, i, h));
        CALL(mess_vector_conj(h));
        CALL(mess_matrix_setcol(Z, i+1, h));
    }
    mess_vector_clear(&h);

    /*-----------------------------------------------------------------------------
     *  Setup final matrix
     *-----------------------------------------------------------------------------*/
    mess_direct sol;
    CALL(mess_direct_init(&sol));
    CALL(mess_direct_lu(Z,sol));
    CALL(mess_direct_inverse(sol, Zorth2));
    mess_direct_clear(&sol);
    for (i = 0; i < n; i++) {
        CALL(mess_matrix_colscale(Z, i, ev_soll->values_cpx[i]));
    }
    CALL( mess_matrix_multiply(MESS_OP_NONE, Z, MESS_OP_NONE, Zorth2, A));
    CALL(mess_matrix_toreal(A));

    CALL(mess_eigen_eig(A, ev_ist, Z));
    for (i = 0; i < n; i++) {
        CALL(check_eigen_vector(A,Z,ev_ist->values_cpx[i], i, tol));
    }

    CALL(mess_vector_sort(ev_soll));
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

