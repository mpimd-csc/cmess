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
 * @file tests/matrix/check_mulnormf2.c
 * @brief Check the Frobenius-norm computation of a matrix.
 * @test
 * This function checks the @ref mess_matrix_mulnormf2 function defined in norm.c that means it checks if the
 * Frobenius-norm of a matrix \f$ G \f$
 * \f[ \Vert op(A)*op(B) \Vert_F^2 \f]
 * is computed correctly.
 *
 * @}
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "mess/config.h"

#include "mess/mess.h"
#include "mess/error_macro.h"

#include "../call_macro.h"


int main ( int argc, char **argv){
    mess_init();
    mess_error_level = 0;

    mess_matrix A,B,TMP1;
    mess_int_t rowsA=6, colsA=6, rowsB=6, colsB=6;
    mess_int_t irowsA, icolsA, irowsB, icolsB;

    mess_operation_t  ops[] = {MESS_OP_NONE, MESS_OP_HERMITIAN, MESS_OP_TRANSPOSE};
    mess_datatype_t   dts[] = {MESS_REAL, MESS_COMPLEX};
    mess_storage_t    sts[] = {MESS_DENSE, MESS_CSR, MESS_CSC, MESS_COORD};
    mess_int_t iopsA, iopsB, idtsA, idtsB, istsA, istsB;

    double p = 0.5, nrm, ref, tol = sqrt(mess_eps());
    int ret;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 1 ) {
        fprintf(stderr, "usage: %s\n", argv[0]);
        return 1;
    }


    /*-----------------------------------------------------------------------------
     *  test function
     *-----------------------------------------------------------------------------*/

    //rows and cols of A
    for(irowsA=1;irowsA<rowsA;++irowsA){
        for(icolsA=1;icolsA<colsA;++icolsA){

            //rows and cols of B
            for(irowsB=1;irowsB<rowsB;++irowsB){
                for(icolsB=1;icolsB<colsB;++icolsB){

                    //datatype of A and B
                    for(idtsA=0;idtsA<2;++idtsA){
                        for(idtsB=0;idtsB<2;++idtsB){

                            //storage type of A and B
                            for(istsA=0;istsA<4;++istsA){
                                for(istsB=0;istsB<4;++istsB){

                                    //get storage and datatype
                                    mess_storage_t stA =  sts[istsA];
                                    mess_storage_t stB =  sts[istsB];
                                    mess_datatype_t dtA = dts[idtsA];
                                    mess_datatype_t dtB = dts[idtsB];

                                    //operation type for A and B
                                    for(iopsA=0;iopsA<3;iopsA++){
                                        for(iopsB=0;iopsB<3;iopsB++){

                                            mess_operation_t opA = ops[iopsA];
                                            mess_operation_t opB = ops[iopsB];

                                            //check if product is possible
                                            mess_int_t l1 = opA==MESS_OP_NONE?colsA:rowsA;
                                            mess_int_t l2 = opB==MESS_OP_NONE?rowsB:colsB;

                                            if(l1==l2){

                                                // create random matrices
                                                MESS_INIT_MATRICES(&A,&B,&TMP1);
                                                CALL(mess_matrix_rand(A, rowsA, colsA, stA, dtA, p));
                                                CALL(mess_matrix_rand(B, rowsB, colsB, stB, dtB, p));
                                                CALL(mess_matrix_mulnormf2(opA,A,opB,B,&nrm));

                                                CALL(mess_matrix_multiply(opA,A,opB,B,TMP1));
                                                CALL(mess_matrix_normf2(TMP1,&ref));

                                                if (fabs(nrm-ref)>tol){
                                                    printf("FAILED\n");
                                                    printf("opA=%s\n",mess_operation_t_str(opA));
                                                    printf("opB=%s\n",mess_operation_t_str(opB));

                                                    printf("Matrix A\n");
                                                    mess_matrix_printinfo(A);

                                                    printf("Matrix B\n");
                                                    mess_matrix_printinfo(B);

                                                    return 1;
                                                }
                                                MESS_CLEAR_MATRICES(&A,&B,&TMP1);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&B,&TMP1);


    return 0;
}

