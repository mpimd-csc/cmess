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
 * @addtogroup test_direct
 * @{
 * @file tests/direct/inverse.c
 * @brief Check the computation of the inverse of a matrix.
 * @author @koehlerm
 * @test
 * This function checks the @ref mess_direct_inverse function defined in direct.c that means it checks if the computed
 * inverse of a matrix \f$ A \f$ fullfills
 * \f[ A A^{-1}=I, \f]
 * where \f$ I \f$ denotes the identity matrix. \n
 * The solver can be selected by the third argument for the second input:
 * <ul>
 * <li> 0 - \ref mess_direct_create_lapack_lu,
 * <li> 1 - \ref mess_direct_create_csparse_lu,
 * <li> 2 - \ref mess_direct_create_umfpack,
 * <li> 3 - \ref mess_direct_create_sparse_lu,
 * <li> 4 - \ref mess_direct_create_cholesky,
 * <li> 5 - \ref mess_direct_create_csparse_cholesky.
 * </ul>
 *
 * @}
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"

#include "../call_macro.h"

/**
 * @brief Select a solver.
 * @param[in] i      input solver to choose
 * @param[in] A          input matrix
 * @param[out] sol  generated solver
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref select_solver function generates a solver for a matrix \f$ A \f$. \n
 * Solvers are choosen with \f$ i \f$. Possible values are:
 * <ul>
 * <li> 0 - \ref mess_direct_create_lapack_lu,
 * <li> 1 - \ref mess_direct_create_csparse_lu,
 * <li> 2 - \ref mess_direct_create_umfpack,
 * <li> 3 - \ref mess_direct_create_sparse_lu,
 * <li> 4 - \ref mess_direct_create_cholesky,
 * <li> 5 - \ref mess_direct_create_csparse_cholesky.
 * </ul>
 *
 */

int select_solver(int i, mess_matrix A, mess_direct sol) {
    int ret =0;
    switch(i){
        case 0:
            CALL(mess_direct_create_lapack_lu(A,sol));
            break;
        case 1:
#ifdef MESS_HAVE_CSPARSE
            CALL(mess_direct_create_csparse_lu(A,sol));
            break;
#endif
        case 2:
#ifdef MESS_HAVE_UMFPACK
            CALL(mess_direct_create_umfpack(A,sol));
            break;
#endif
        case 3:
            CALL(mess_direct_create_sparse_lu(A,sol));
            break;
        case 4:
            CALL(mess_direct_create_cholesky(A,sol));
            break;
        case 5:
#ifdef MESS_HAVE_CSPARSE
            CALL(mess_direct_create_csparse_cholesky(A,sol));
            break;
#endif
        case 6:
#ifdef MESS_HAVE_CHOLMOD
            CALL(mess_direct_create_cholmod_cholesky(A,sol));
            break;
#endif
        default:

            return 1;

    }
    return 0;
}

/**
 * @brief Test inverting a matrix.
 * @param[in] i      input solver to choose
 * @param[in] A          input matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref testmat function generates an inverse for a matrix \f$ A \f$ and checks if
 * \f[ A A^{-1}= I, \f]
 * where \f$ I \f$ denotes the identity matrix.
 * Parameter \f$ i \f$ is only used to choose the solver. Possible values are:
 * <ul>
 * <li> 0 - \ref mess_direct_create_lapack_lu,
 * <li> 1 - \ref mess_direct_create_csparse_lu,
 * <li> 2 - \ref mess_direct_create_umfpack,
 * <li> 3 - \ref mess_direct_create_sparse_lu,
 * <li> 4 - \ref mess_direct_create_cholesky,
 * <li> 5 - \ref mess_direct_create_csparse_cholesky.
 * </ul>
 *
 */
int testmat(int i, mess_matrix A) {
    mess_direct sol;
    mess_matrix inv, eye, t;
    double res1;
    double nrm;
    int ret = 0 ;
    int err = 0;

    CALL(mess_direct_init(&sol));
    CALL(mess_matrix_init(&eye));
    CALL(mess_matrix_init(&inv));
    CALL(mess_matrix_init(&t));

    CALL(mess_matrix_norm2(A,&nrm));
    CALL(select_solver(i,A,sol));
    CALL(mess_direct_inverse(sol,inv));
    CALL(mess_matrix_eye(eye, A->rows, A->cols, MESS_CSR));
    CALL( mess_matrix_multiply(MESS_OP_NONE, A, MESS_OP_NONE, inv, t));
    mess_matrix_printinfo(inv);
    mess_matrix_printshort(t);
    CALL(mess_matrix_diffnorm(t, eye, &res1));

    printf("res1 = %lg nrm=%lg\n", res1, nrm);

    if ( res1 > 1e-10*MESS_MAX(nrm,1)) err++;
    if (err ) {
        mess_matrix_printshort(t);
        mess_matrix_printshort(eye);
    }
    mess_direct_clear(&sol);
    mess_matrix_clear(&eye);
    mess_matrix_clear(&inv);
    mess_matrix_clear(&t);
    return err ;
}


int main ( int argc, char ** argv) {
    mess_matrix mat_coord;
    mess_matrix mat_csr;
    mess_matrix mat_csc;
    mess_matrix mat_dense;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if ( argc != 2 && argc != 3) {
        printf("check the direct linear solvers\n");
        printf("usage: %s op matrix.mtx\n", argv[0]);
        printf(" op: 0 -> lapacklu\n");
        return 1;
    }
    /*-----------------------------------------------------------------------------
     *  init matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_init(&mat_coord));
    CALL(mess_matrix_init(&mat_csr));
    CALL(mess_matrix_init(&mat_dense));
    CALL(mess_matrix_init(&mat_csc));


    /*-----------------------------------------------------------------------------
     *  read matrices
     *-----------------------------------------------------------------------------*/
    if (argc == 3 ) {
        CALL(mess_matrix_read (argv[2], mat_coord));
    } else {
        CALL(mess_matrix_rand(mat_dense, 50,50, MESS_DENSE, MESS_REAL,0.5));
        CALL(mess_matrix_convert(mat_dense, mat_coord,MESS_COORD));
    }
    CALL(mess_matrix_convert(mat_coord, mat_csr, MESS_CSR));
    CALL(mess_matrix_convert(mat_coord, mat_csc, MESS_CSC));
    CALL(mess_matrix_convert(mat_coord, mat_dense, MESS_DENSE));


    /*-----------------------------------------------------------------------------
     *  test with CSR
     *-----------------------------------------------------------------------------*/
    CALL(testmat(atoi(argv[1]),mat_csr));
    CALL(testmat(atoi(argv[1]),mat_csc));
    CALL(testmat(atoi(argv[1]),mat_coord));
    CALL(testmat(atoi(argv[1]),mat_dense));


    mess_matrix_clear(&mat_coord);
    mess_matrix_clear(&mat_csr);
    mess_matrix_clear(&mat_dense);
    mess_matrix_clear(&mat_csc);

    return 0;
}

