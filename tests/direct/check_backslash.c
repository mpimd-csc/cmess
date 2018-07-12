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
 * @file tests/direct/check_backslash.c
 * @brief Check the solution of a system of linear equations.
 * @author  @mbehr
 * @test
 *
 * This function checks the @ref mess_matrix_backslash and
 * @ref mess_matrix_backslashm functions defined in backslash.c that means
 * it checks if the solution of
 * \f[ op(A) x= b \f]
 * and
 * \f[ op(A) X= B \f]
 * is computed correctly for given matrices \f$ A \f$ and \f$ B \f$ and given vector \f$ b \f$. \n
 * Operations \f$ op (.)\f$ on \f$ A \f$ can be
 * * \ref MESS_OP_NONE (\f$ op(A) = A \f$),
 * * \ref MESS_OP_TRANSPOSE (\f$ op(A) =A^T \f$),
 * * \ref MESS_OP_HERMITIAN (\f$ op(A)= A^H \f$).
 *
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



#define CHECK_BACKSLASH(OP,A,B,X,TEMP,ERR,EPS,DIFF){                            \
    CALL(mess_matrix_backslash(OP,A,B,X));                                  \
    CALL(mess_matrix_mvp(OP,A,X,TEMP));                                     \
    CALL(mess_vector_diffnorm(B,TEMP,&DIFF));                               \
    if(DIFF>EPS){                                                           \
        printf("Failed:\n");                                                \
        printf("OP: %s\n",#OP);                                             \
        CALL(mess_matrix_printinfo(A));                                     \
        CALL(mess_vector_printinfo(B));                                     \
        printf("Diff: %e\n",DIFF);                                          \
        ++ERR;                                                              \
    }                                                                       \
}


#define CHECK_BACKSLASHM(OP,A,B,X,TEMP,ERR,EPS,DIFF){                           \
    CALL(mess_matrix_backslashm(OP,A,B,X));                                 \
    CALL( mess_matrix_multiply(OP, A, MESS_OP_NONE, X, TEMP));                  \
    CALL(mess_matrix_diffnorm(B,TEMP,&DIFF));                               \
    if(DIFF>EPS){                                                           \
        printf("BACKSLASHM Failed:\n");                                     \
        printf("OP: %s\n",#OP);                                             \
        CALL(mess_matrix_printinfo(A));                                     \
        CALL(mess_matrix_printinfo(B));                                     \
        printf("Diff: %e\n",DIFF);                                          \
        ++ERR;                                                              \
    }                                                                       \
}



int main ( int argc, char ** argv) {
    mess_init();
    int  ret=0, err=0, rectangular;
    double eps = sqrt(mess_eps()),  diff=.0;

    mess_vector  x, b, temp;
    mess_matrix  A, X, B, B2, Temp;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if ( argc != 3) {
        printf("usage: %s matrix.mtx [rectangular]\n", argv[0]);
        return 1;
    }

    rectangular = atoi(argv[2]);

    /*-----------------------------------------------------------------------------
     * Init/load matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&X,&B,&B2,&Temp);
    CALL(mess_matrix_read_formated(argv[1], A, MESS_DENSE ));

    if(rectangular){
        CALL(mess_matrix_sub(A,0,A->rows-1, 0, MESS_MAX((A->cols)/2,1),Temp));
        CALL(mess_matrix_copy(Temp,A));
    }

    CALL(mess_matrix_rand(B,A->rows,A->cols,MESS_DENSE,MESS_REAL,0.5));
    CALL(mess_matrix_alloc(Temp,A->cols,A->rows,A->cols*A->rows,MESS_DENSE,MESS_REAL));

    /*-----------------------------------------------------------------------------
     * Init/load vectors
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&b,&x,&temp);
    CALL(mess_vector_alloc(b,A->rows,MESS_REAL));
    CALL(mess_vector_alloc(x,A->cols,MESS_REAL));
    CALL(mess_vector_alloc(temp,A->rows,MESS_REAL));
    CALL(mess_vector_rand(b))

        /*-----------------------------------------------------------------------------
         *  Test different cases
         *-----------------------------------------------------------------------------*/
        int storetype=0;

    for(storetype=0;storetype<4;++storetype){

        switch(storetype){
            case 0:
                CALL(mess_matrix_convert(A,Temp,MESS_DENSE));
                CALL(mess_matrix_copy(Temp,A));
                break;
            case 1:
                CALL(mess_matrix_convert(A,Temp,MESS_CSC));
                CALL(mess_matrix_copy(Temp,A));
                break;
            case 2:
                CALL(mess_matrix_convert(A,Temp,MESS_CSR));
                CALL(mess_matrix_copy(Temp,A));
                break;
            case 3:
                CALL(mess_matrix_convert(A,Temp,MESS_COORD));
                CALL(mess_matrix_copy(Temp,A));
                break;
            default:
                break;
        }

        CHECK_BACKSLASH(MESS_OP_NONE,A,b,x,temp,err,eps,diff);
        CHECK_BACKSLASH(MESS_OP_TRANSPOSE,A,b,x,temp,err,eps,diff);
        CHECK_BACKSLASH(MESS_OP_HERMITIAN,A,b,x,temp,err,eps,diff);

        CHECK_BACKSLASHM(MESS_OP_NONE    ,A,B,X,Temp,err,eps,diff);
        CHECK_BACKSLASHM(MESS_OP_TRANSPOSE,A,B,X,Temp,err,eps,diff);
        CHECK_BACKSLASHM(MESS_OP_HERMITIAN,A,B,X,Temp,err,eps,diff);
    }


    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&B,&B2,&X,&Temp);
    MESS_CLEAR_VECTORS(&x,&b,&temp);


    return (err>0)?(1):(0);
}
