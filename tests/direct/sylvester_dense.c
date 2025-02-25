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
 * @file tests/direct/sylvester_dense.c
 * @brief Check the solution of dense Sylvester equations.
 * @test
 * This function checks if the solver generated by @ref mess_direct_create_sylvester_dense solves the dense Sylvester equations:
 * <ul>
 * <li>  \f[ A X   +   X H   + M = 0 \f], \f[ A^T X     + X H^T     + M = 0 \f], \f[ A^H X     + X H^H     + M = 0 \f]
 * <li>  \f[ A X F + E X H + M = 0 \f],   \f[ A^T X F^T + E^T X H^T + M = 0 \f], \f[ A^T X F^H + E^H X H^H + M = 0 \f]
 * </ul>
 * correctly.
 * \n
 * All matrix \f$ A, E, F, H, M \f$ must be dense.
 *
 * @author @mbehr
 *
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "../call_macro.h"
#include <complex.h>

/**
 * @internal
 * @brief Print information about test instance.$
 *
 * @attention Internal use only.
 */
void test_sylvester_info(mess_direct sylv, mess_operation_t op, mess_matrix A,mess_matrix F, mess_matrix E, mess_matrix H, mess_matrix M, double res2, double res2_0){

    /*-----------------------------------------------------------------------------
     *  print informations about data
     *-----------------------------------------------------------------------------*/
    printf("Solver datatype = %s\n",mess_datatype_t_str(sylv->data_type));
    printf("Solver Name     = %s\n",sylv->name);

    printf("Matrix A        = %s, %s\n",mess_storage_t_str(A->store_type), mess_datatype_t_str(A->data_type));
    if(F){
        printf("Matrix F        = %s, %s\n",mess_storage_t_str(F->store_type), mess_datatype_t_str(F->data_type));
    }else{
        printf("Matrix F        = NULL \n");
    }
    if(E){
        printf("Matrix E        = %s, %s\n",mess_storage_t_str(E->store_type), mess_datatype_t_str(E->data_type));
    }else{
        printf("Matrix E        = NULL \n");
    }
    printf("Matrix H        = %s, %s\n",mess_storage_t_str(H->store_type), mess_datatype_t_str(H->data_type));
    printf("Matrix M        = %s, %s\n",mess_storage_t_str(M->store_type), mess_datatype_t_str(M->data_type));
    printf("op              = %s\n",mess_operation_t_str(op));
    printf("size(X)         = %d-by-%d\n",A->rows,H->rows);

    /*-----------------------------------------------------------------------------
     *  print residual informations
     *-----------------------------------------------------------------------------*/
    if(op==MESS_OP_NONE){
        if(E){
            if(F){
                printf("|| AXF + EXH + M ||_2                       = %e\n", res2);
                printf("|| AXF + EXH + M ||_2 / || M ||_2           = %e\n", res2/res2_0);
            }else{
                printf("|| AX + EXH + M ||_2                        = %e\n", res2);
                printf("|| AX + EXH + M ||_2 / || M ||_2            = %e\n", res2/res2_0);
            }
        }else{
            printf("|| AX +  XH + M ||_2                        = %e\n", res2);
            printf("|| AX +  XH + M ||_2 / || M ||_2            = %e\n", res2/res2_0);
        }
    }else if (op==MESS_OP_TRANSPOSE){
        if(E){
            if(F){
                printf("|| A.'XF.' + E.'XH.' + M ||_2               = %e\n", res2);
                printf("|| A.'XF.' + E.'XH.' + M ||_2 / || M ||_2   = %e\n", res2/res2_0);
            }else{
                printf("|| A.'X + E.'XH.' + M ||_2                  = %e\n", res2);
                printf("|| A.'X + E.'XH.' + M ||_2 / || M ||_2      = %e\n", res2/res2_0);
            }
        }else{
            printf("|| A.'X +  XH.' + M ||_2                    = %e\n", res2);
            printf("|| A.'X +  XH.' + M ||_2 / || M ||_2        = %e\n", res2/res2_0);
        }
    }else{
        if(E){
            if(F){
                printf("|| A'XF' + E'XH' + M ||_2                   = %e\n", res2);
                printf("|| A'XF' + E'XH' + M ||_2 / || M ||_2       = %e\n", res2/res2_0);
            }else{
                printf("|| A'X + E'XH' + M ||_2                     = %e\n", res2);
                printf("|| A'X + E'XH' + M ||_2 / || M ||_2         = %e\n", res2/res2_0);
            }
        }else{
            printf("|| A'X +  XH' + M ||_2                      = %e\n", res2);
            printf("|| A'X +  XH' + M ||_2 / || M ||_2          = %e\n", res2/res2_0);
        }
    }
    printf("|| M ||_2                                   = %e\n", res2_0);
    printf("\n");
}


int main ( int argc, char **argv){

    int ret = 0, Xdatatype_test_failed = 0;
    const int verbose = 1;
    int generalized = 0;
    double res2 = 0, res2_0 = 0, tol = sqrt(mess_eps());
    mess_int_t seed = 1234;
    mess_int_t dimA, dimH;
    mess_matrix A, H, M, X, F=NULL, E=NULL;
    mess_direct sylv;
    mess_error_level = 3;

    mess_operation_t ops [] = {MESS_OP_NONE, MESS_OP_TRANSPOSE, MESS_OP_HERMITIAN};
    const int length_ops = sizeof(ops)/sizeof(mess_operation_t);
    int i_op;

    mess_datatype_t dts [] = {MESS_REAL, MESS_COMPLEX};
    //mess_datatype_t dts [] = {MESS_REAL};
    const int length_dts = sizeof(dts)/sizeof(mess_datatype_t);
    int i_dtA, i_dtE, i_dtH, i_dtM, i_dtF;

    mess_int_t dimsA [] = {1,2,3,5,10,33,40};
    const int length_dimsA = sizeof(dimsA)/sizeof(mess_int_t);

    mess_int_t dimsH [] = {1,2,3,5,10,33,40};
    const int length_dimsH = sizeof(dimsH)/sizeof(mess_int_t);


    /*-----------------------------------------------------------------------------
     *  check input arguments and read the first argument
     *-----------------------------------------------------------------------------*/
    if (argc != 2){
        printf("usage: %s [0|1](standard/generalized) \n", argv[0]);
        return 1;
    }

    generalized = atoi(argv[1]);

    /*-----------------------------------------------------------------------------
     *  init and and create random matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&H,&M,&X);
    if(generalized){
        MESS_INIT_MATRICES(&E,&F);
    }

    CALL(mess_matrix_rand_init(&seed));

    /*-----------------------------------------------------------------------------
     *  test mess_direct_sylvester_sparsedense with different cases
     *-----------------------------------------------------------------------------*/
    /* iterate over all possible datatypes for E and F, if generalized == 0 this loop is encountered exactly ones*/
    for(i_dtF=generalized?0:length_dts-1, i_dtE=generalized?0:length_dts-1; i_dtF<length_dts; i_dtF += (i_dtE+1)/length_dts, i_dtE = (i_dtE+1)%length_dts){

        /* iterate over all possible datatypes for A */
        for(i_dtA=0;i_dtA<length_dts;++i_dtA){

            /* iterate over all possible dimension for A */
            for(dimA=0;dimA<length_dimsA;++dimA){

                /* iterate over all possible datatypes for H*/
                for(i_dtH=0;i_dtH<length_dts;++i_dtH){

                    /* iterate over all possible dimension for H */
                    for(dimH=0;dimH<length_dimsH;++dimH){

                        /* iterate over all possible datatypes for M*/
                        for(i_dtM=0;i_dtM<length_dts;++i_dtM){


                            /*-----------------------------------------------------------------------------
                             *  create random matrices
                             *----------------------------------------------------------------------------*/
                            CALL(mess_matrix_rand_dense(A,dimsA[dimA],dimsA[dimA],dts[i_dtA]));
                            if(generalized){
                                CALL(mess_matrix_rand_dense(F,dimsH[dimH],dimsH[dimH],dts[i_dtF]));
                                CALL(mess_matrix_rand_dense(E,dimsA[dimA],dimsA[dimA],dts[i_dtE]));
                            }
                            CALL(mess_matrix_rand_dense(H,dimsH[dimH],dimsH[dimH],dts[i_dtH]));
                            CALL(mess_matrix_rand_dense(M,dimsA[dimA],dimsH[dimH],dts[i_dtM]));
                            CALL(mess_matrix_norm2(M,&res2_0));

                            /*-----------------------------------------------------------------------------
                             *  create sparse dense solver and solve
                             *-----------------------------------------------------------------------------*/
                            CALL(mess_direct_init(&sylv));
                            CALL(mess_direct_create_sylvester_dense(A, F, E, H, sylv));

                            /*-----------------------------------------------------------------------------
                             * solve dense sylvester equation and  print residual
                             *-----------------------------------------------------------------------------*/
                            /* iterate over all possible operation types*/
                            for(i_op=0;i_op<length_ops;++i_op){

                                /*-----------------------------------------------------------------------------
                                 *  solve sylvester equation and compute residual
                                 *-----------------------------------------------------------------------------*/
                                CALL(mess_direct_solvem(ops[i_op],sylv,M,X));
                                if(generalized){
                                    CALL(mess_direct_generalized_sylvester_res2(ops[i_op],A,F,E,H,M,X,&res2));
                                }else{
                                    CALL(mess_direct_sylvester_res2(ops[i_op],A,H,M,X,&res2));
                                }

                                /*-----------------------------------------------------------------------------
                                 *  print infos and check residual and datatype
                                 *-----------------------------------------------------------------------------*/
                                if(verbose){
                                    test_sylvester_info(sylv, ops[i_op], A, F, E, H, M, res2, res2_0);
                                }

                                /*-----------------------------------------------------------------------------
                                 *  check data_type and residual
                                 *-----------------------------------------------------------------------------*/
                                if(MESS_IS_REAL(A) && MESS_IS_REAL(H) && MESS_IS_REAL(M)){
                                    if(!F){
                                        Xdatatype_test_failed =(!generalized || MESS_IS_REAL(E))?(MESS_IS_REAL(X)? 0:1):(MESS_IS_COMPLEX(X)? 0:1);
                                    }else{
                                        if(MESS_IS_COMPLEX(F)){
                                            Xdatatype_test_failed = MESS_IS_COMPLEX(X)? 0:1;
                                        }else{
                                            Xdatatype_test_failed =(!E || MESS_IS_REAL(E))?(MESS_IS_REAL(X)? 0:1):(MESS_IS_COMPLEX(X)? 0:1);
                                        }
                                    }
                                }else{
                                    Xdatatype_test_failed = MESS_IS_COMPLEX(X)? 0:1;
                                }

                                if(res2 > tol || res2/res2_0 > tol || Xdatatype_test_failed ){
                                    printf("Failed:\n");
                                    printf("tol             = %e\n", tol);
                                    printf("Residual in 2-Norm\n");
                                    printf("abs. residual   = %e\n", res2);
                                    printf("rel. residual   = %e\n", res2/res2_0);
                                    printf("ops             = %s\n", mess_operation_t_str(ops[i_op]));
                                    printf("Datatype X      = %s\n", mess_datatype_t_str(X->data_type));
                                    printf("Matrix A:\n"); mess_matrix_printshort(A);
                                    if(F){
                                        printf("Matrix F:\n"); mess_matrix_printshort(F);
                                    }
                                    if(E){
                                        printf("Matrix E:\n"); mess_matrix_printshort(E);
                                    }
                                    printf("Matrix H:\n"); mess_matrix_printshort(H);
                                    printf("Matrix M:\n"); mess_matrix_printshort(M);
                                    printf("Matrix X:\n"); mess_matrix_printshort(X);
                                    printf("\n");
                                    return 1;
                                }

                            }

                            /*-----------------------------------------------------------------------------
                             *  clear solver for the next round
                             *-----------------------------------------------------------------------------*/
                            mess_direct_clear(&sylv);
                        }
                    }
                }
            }
        }
    }
    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&H,&M,&X);
    if(F){MESS_CLEAR_MATRICES(&F);}
    if(E){MESS_CLEAR_MATRICES(&E);}

    return 0;
}

