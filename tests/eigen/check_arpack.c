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
 *
 * @file tests/eigen/check_arpack.c
 * @brief Check the solution of the eigenvalue problem using the Arnoldi process from @arpack.
 * @test
 * This function checks the @ref mess_eigen_arpack and @ref mess_eigen_arpack_template functions defined in arpack.c that means it
 * checks if the eigenvalue problem \f[ Ax = \lambda x \f]
 * is solved correctly for a given matrix \f$ A \f$ using the Arnoldi process from @arpack.
 *
 *
 * @}
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/mess.h"

#include "../call_macro.h"

/**
 * @brief Check if the eigenvalue problem is solved correctly.
 * @param[in] A  input matrix
 * @param[in] ev     input vector containing eigenvalues
 * @param[in] V      input matrix containing eigenvectors
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref check_ev function checks if the eigenvalue problem
 * \f[ AV(:,i)= ev(i) V(:,i) \f]
 * is solved correctly for \f$ i= 1, \cdots , dim(ev) \f$. \n
 * If
 * \f[ \frac{\Vert A V(:,i)- ev(i) V(:,i) \Vert_2}{ \Vert V(:,i) \Vert_2} > \sqrt{ \epsilon n }\f]
 * the error value is increased by \f$ 1 \f$. \n
 * \f$ n \f$ is the number of rows of \f$ V \f$.
 *
 */
int check_ev(mess_matrix A, mess_vector ev, mess_matrix V ) {
    mess_int_t i;
    mess_vector EV,R;
    double norm, res;
    int err = 0;
    int ret = 0;

    MESS_INIT_VECTORS(&EV,&R);
    CALL(mess_vector_alloc(EV,A->rows,V->data_type));
    CALL(mess_vector_alloc(R,A->rows,V->data_type));

    for (i = 0; i < ev->dim; i++) {
        CALL(mess_matrix_getcol(V,i,EV));
        CALL(mess_vector_norm2(EV,&norm));
        CALL(mess_matrix_mvp(MESS_OP_NONE, A, EV, R));
        if ( MESS_IS_REAL(ev)) {
            CALL(mess_vector_axpyc(-ev->values[i],EV,R));
        } else {
            CALL(mess_vector_axpyc(-ev->values_cpx[i],EV,R));
        }
        CALL(mess_vector_norm2(R,&res));
        printf("%ld \t %lg\n",(long) i,res/norm);
        if ( res / norm > sqrt(mess_eps()*A->rows)) err ++;
    }

    mess_vector_clear(&EV);
    mess_vector_clear(&R);
    return err;
}

/**
 * @brief Check if the eigenvalue problem is solved correctly.
 * @param[in] A      input matrix to compute its eigenvalues from
 * @param[in] nev    input number of desired eigenvalues
 * @param[in] which      input eigenvalues to select in the Arnoldi process
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref call_arnoldi function computes the eigenvalues and eigenvector for the eigenvalue problem
 * \f[ A x= \lambda x\f]
 * and checks if the eigenvalue problem is fulfilled.
 *
 */
int call_arnoldi(mess_matrix A, mess_int_t nev, mess_eigen_arpack_which_t which) {
    mess_vector ev;
    mess_matrix V;
    mess_eigen_arpack_options_t opt = MESS_EIGEN_ARPACK_DEFAULT;
    mess_mvpcall mvpcall;
    int err = 0;
    int ret = 0;

    CALL(mess_matrix_init(&V));
    MESS_INIT_VECTORS(&ev,&opt.b0);
    CALL(mess_vector_alloc(ev, A->rows, MESS_REAL));
    opt.which = which;
    CALL(mess_vector_alloc(opt.b0,A->rows,MESS_REAL));
    CALL(mess_vector_ones(opt.b0));

    CALL(mess_eigen_arpack( A, nev, opt, ev,V));
    mess_vector_print(ev);
    err += check_ev(A,ev,V);

    CALL(mess_matrix_zeros(V));
    CALL(mess_mvpcall_matrix(&mvpcall, MESS_OP_NONE, A));
    CALL(mess_eigen_arpack_template( mvpcall , nev, opt, ev,V));
    CALL(mess_mvpcall_clear(&mvpcall));
    mess_vector_print(ev);
    err += check_ev(A,ev,V);

    mess_matrix_clear(&V);
    mess_vector_clear(&ev);
    mess_vector_clear(&opt.b0);
    return err;

}

int main (int argc, char *argv[])
{
    mess_init();
    mess_matrix matrix;
    mess_matrix input;
    int ret = 0 ;

    mess_version();
    if ( argc < 2 ) {
        fprintf (stderr, "eigenvalues with LAPACK\n");
        fprintf (stderr, "Usage: %s file.mtx\n", argv[0]) ;
        exit(1);
    }
    mess_matrix_init(&input);
    mess_matrix_init(&matrix);

    mess_matrix_read(argv[1], input);
    mess_matrix_convert(input, matrix, MESS_CSR);
    mess_matrix_clear(&input);


    if (matrix->cols != matrix->rows){
        return -1;
    }
    int nev = 3;
    printf("Check MESS_EIGEN_ARPACK_LARGE_MAGNITUDE:\n:");
    CALL(call_arnoldi(matrix,nev,MESS_EIGEN_ARPACK_LARGE_MAGNITUDE));

    printf("Check MESS_EIGEN_ARPACK_SMALL_MAGNITUDE:\n:");
    CALL(call_arnoldi(matrix,nev,MESS_EIGEN_ARPACK_SMALL_MAGNITUDE));

    printf("Check MESS_EIGEN_ARPACK_LARGE_REALPART:\n:");
    CALL(call_arnoldi(matrix,nev,MESS_EIGEN_ARPACK_LARGE_REALPART));

    printf("Check MESS_EIGEN_ARPACK_SMALL_REALPART:\n:");
    CALL(call_arnoldi(matrix,nev,MESS_EIGEN_ARPACK_SMALL_REALPART));

    if ( MESS_IS_COMPLEX(matrix)) {
        printf("Check MESS_EIGEN_ARPACK_LARGE_IMAGPART:\n:");
        CALL(call_arnoldi(matrix,nev,MESS_EIGEN_ARPACK_LARGE_IMAGPART));

        printf("Check MESS_EIGEN_ARPACK_SMALL_IMAGPART:\n:");
        CALL(call_arnoldi(matrix,nev,MESS_EIGEN_ARPACK_SMALL_IMAGPART));
    }
    mess_matrix_clear(&matrix);
    return ret;
}
