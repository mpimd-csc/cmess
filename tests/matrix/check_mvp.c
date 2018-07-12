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
 * @file tests/matrix/check_mvp.c
 * @brief Check matrix-vector product computation.
 * @test
 * This function checks the @ref mess_matrix_mvp function defined in mvp.c that means it checks if the matrix-vector
 * product
 * \f[ op(A) x \f]
 * is computed correctly for a given matrix \f$ A \f$ and a given vector \f$ x \f$.
 * Operations \f$ op (.)\f$ on \f$ A \f$ can be
 * <ul>
 * <li> \f$ 0 \f$: \ref MESS_OP_NONE (\f$ op(A)=A \f$),
 * <li> \f$ 1 \f$: \ref MESS_OP_TRANSPOSE (\f$ op(A)= A^T \f$),
 * <li> \f$ 2 \f$: \ref MESS_OP_HERMITIAN (\f$ op(A)=A^H \f$).
 * </ul>
 * \f$ 1 \f$ is default value.
 * @}
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/mess.h"

#include "../call_macro.h"
int main (int argc, char *argv[])
{
    mess_init();
    mess_matrix mat_coord;
    mess_matrix mat_csr;
    mess_matrix mat_csc;
    mess_matrix mat_dense;
    mess_vector x, ycoord, ycsr, ycsc, ydense;
    mess_operation_t op;
    double nrmA;
    double eps = mess_eps();
    double tol;

    mess_int_t i;
    mess_int_t n,m;
    int ret;
    int err = 0;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if ( argc != 3) {
        printf("compute the norm of a matrix\n");
        printf("usage: %s matrix.mtx op\n", argv[0]);
        printf(" op: 0 -> y=Ax , 1-> y=A^Hx , 2-> y = A^Tx");
        return 1;
    }
    switch(atoi(argv[2])){
        case 0:
            op = MESS_OP_NONE; break;
        case 1:
            op = MESS_OP_HERMITIAN; break;
        case 2:
            op = MESS_OP_TRANSPOSE; break;
        default: return 1;
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
    CALL(mess_matrix_read (argv[1], mat_coord));
    CALL(mess_matrix_convert(mat_coord, mat_csr, MESS_CSR));
    CALL(mess_matrix_convert(mat_coord, mat_csc, MESS_CSC));
    CALL(mess_matrix_convert(mat_coord, mat_dense, MESS_DENSE));
    CALL(mess_matrix_norm2(mat_dense, &nrmA));
    printf("norm2 = %lg \t norm2*eps*dim=%lg\n",nrmA, nrmA*eps*MESS_MAX(mat_csr->rows, mat_csr->cols));
    tol = 10.0*nrmA*eps*MESS_MAX(mat_csr->rows,mat_csr->cols);


    /*-----------------------------------------------------------------------------
     *  init vectors
     *-----------------------------------------------------------------------------*/
    if ( op == MESS_OP_NONE ) n = mat_coord->cols; else n = mat_coord->rows;
    if ( op == MESS_OP_NONE ) m = mat_coord->rows; else m = mat_coord->cols;

    MESS_INIT_VECTORS(&x,&ycoord,&ycsr,&ycsc,&ydense);
    CALL(mess_vector_alloc(x, n, mat_coord->data_type));
    if ( MESS_IS_COMPLEX(mat_coord)) {
        for ( i = 0; i < n; i++) x->values_cpx[i] = 1;
    } else {
        for ( i = 0; i < n; i++) x->values[i] = 1;
    }
    CALL(mess_vector_alloc(ycoord, m, mat_coord->data_type));
    CALL(mess_vector_alloc(ycsr, m, mat_coord->data_type));
    CALL(mess_vector_alloc(ycsc, m, mat_coord->data_type));
    CALL(mess_vector_alloc(ydense, m, mat_coord->data_type));



    /*-----------------------------------------------------------------------------
     * Compute
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_mvp( op ,mat_coord, x, ycoord));
    CALL(mess_matrix_mvp( op ,mat_csr, x, ycsr));
    CALL(mess_matrix_mvp( op ,mat_csc, x, ycsc));
    CALL(mess_matrix_mvp( op ,mat_dense, x, ydense));


    /*-----------------------------------------------------------------------------
     *  check output
     *-----------------------------------------------------------------------------*/
    double diff_csr, diff_csc, diff_coord;
    CALL(mess_vector_diffnorm(ydense, ycsr, &diff_csr));
    CALL(mess_vector_diffnorm(ydense, ycsc, &diff_csc));
    CALL(mess_vector_diffnorm(ydense, ycoord, &diff_coord));

    printf("diff csr:    %lg\n", diff_csr);
    printf("diff csc:    %lg\n", diff_csc);
    printf("diff coord:  %lg\n", diff_coord);

    if ( diff_csr > tol) {
        err ++;
        printf("ERROR DENSE - CSR\n");
        mess_matrix_printinfo(mat_dense);
        mess_matrix_printinfo(mat_csr);
    }
    if ( diff_csc > tol){
        err ++;
        printf("ERROR DENSE - CSC\n");
        mess_matrix_printinfo(mat_dense);
        mess_matrix_printinfo(mat_csc);
    }
    if ( diff_coord > tol){
        err ++;
        printf("ERROR DENSE - COORD\n");
        mess_matrix_printinfo(mat_dense);
        mess_matrix_printinfo(mat_coord);
    }



    if ( err ) {
        printf("FAILED\n");
    } else {
        printf("PASSED\n");
    }
    mess_vector_clear(&x);
    mess_vector_clear(&ycsr);
    mess_vector_clear(&ycsc);
    mess_vector_clear(&ycoord);
    mess_vector_clear(&ydense);

    mess_matrix_clear(&mat_coord);
    mess_matrix_clear(&mat_csr);
    mess_matrix_clear(&mat_dense);
    mess_matrix_clear(&mat_csc);

    return err;
}
/** \}@ */

