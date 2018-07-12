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
 *
 * @file tutorials/eig/tutorial_eigs.c
 * @brief Demonstrate how to approximate some eigenvalues of a matrix.
 *
 * # Tutorial: Eigenvalue Computation for Sparse Matrices
 *
 * This function demonstrates how to approximate some eigenvalues of the eigenvalue problem
 * \f[ A x = \lambda x, \f]
 * where \f$ A \f$ is a matrix, \f$ \lambda \f$ an eigenvalue and \f$ x \f$ its corresponding eigenvector
 * using the Arnoldi process.
 *
 * @sa mess_eigen_arpack_template
 * @sa mess_eigen_eigs
 *
 * @snippet "tutorials/eig/tutorial_eigs.c" CODE
 *
 */

///@cond
///[CODE]
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/mess.h"

#ifdef MESS_HAVE_ARPACK

/**
 * @brief Compute a matrix-vector product.
 * @param[in] data  input matrix data
 * @param[in] x     input vector
 * @param[out] y        output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mvp function computes the matrix-vector product between \f$ data \f$ and an input vector \f$ x \f$.
 *
 * @sa mess_matrix_mvp
 *
 */
int mvp(void * data, mess_vector x, mess_vector y) {
    mess_matrix A = (mess_matrix) data;
    mess_matrix_mvp(MESS_OP_NONE, A, x, y);
    return 0;
}
#endif

int main (int argc, char *argv[])
{
    mess_matrix matrix;
    mess_vector ev;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 2 ) {
        printf ("Usage: %s file.mtx\n", argv[0]) ;
        exit(1);
    }

    /*-----------------------------------------------------------------------------
     *  init and read matrics
     *-----------------------------------------------------------------------------*/
    mess_matrix_init(&matrix);
    mess_matrix_read_formated(argv[1], matrix, MESS_CSR);
    if(matrix->rows!=matrix->cols){
        MESS_CLEAR_MATRICES(&matrix);
        printf("Failed:Matrix has to be square.\n");
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  compute and print eigenvalues
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&ev);
    mess_vector_alloc(ev, matrix->rows, MESS_REAL);

    mess_eigen_eigs(matrix,30, 0, NULL, ev);
    printf("Eigenvalues:\n");
    mess_vector_print(ev);


    /*-----------------------------------------------------------------------------
     *  compute eigenvalues with arpack if available
     *-----------------------------------------------------------------------------*/
#ifdef MESS_HAVE_ARPACK
    mess_matrix V;
    mess_vector EV,R;
    mess_int_t i;
    mess_int_t err = 0;
    double norm, res;
    mess_eigen_arpack_options_t opt = MESS_EIGEN_ARPACK_DEFAULT;
    mess_mvpcall mvpcall;
    mess_mvpcall_matrix(&mvpcall, MESS_OP_NONE, matrix);
    opt.which = MESS_EIGEN_ARPACK_LARGE_MAGNITUDE;

    mess_matrix_init(&V);
    mess_vector_init(&EV);
    mess_vector_init(&R);
    mess_vector_alloc(EV,matrix->rows,MESS_REAL);
    mess_vector_alloc(R,matrix->rows,MESS_REAL);

    printf("ARPACK:\n");
    mess_eigen_arpack_template( mvpcall , 5, opt, ev,V);
    mess_mvpcall_clear(&mvpcall);
    mess_vector_print(ev);
    for (i = 0; i < ev->dim; i++) {
        mess_matrix_getcol(V,i,EV);
        mess_vector_norm2(EV,&norm);
        mess_matrix_mvp(MESS_OP_NONE, matrix, EV, R);
        if ( MESS_IS_REAL(ev)) {
            mess_vector_axpyc(-ev->values[i],EV,R);
        } else {
            mess_vector_axpyc(-ev->values_cpx[i],EV,R);
        }
        mess_vector_norm2(R,&res);
        printf("%ld \t %lg\n",(long) i,res/norm);
        if ( res / norm > sqrt(mess_eps()*matrix->rows)) err ++;
    }
    mess_matrix_clear(&V);
    mess_vector_clear(&EV);
    mess_vector_clear(&R);
    if ( err > 0 ) {
        printf("Eigenvalue Computations may be wrong\n");
        return 1;
    } else {
        printf("Eigenvalue Computations are correct.\n");
    }
#endif

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    mess_matrix_clear(&matrix);
    mess_vector_clear(&ev);

    return 0;
}
///[CODE]
///@endcond
