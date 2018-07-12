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
 * @file tests/direct/sylvester_sparsedense.c
 * @brief Check the solution of a sparse-dense Sylvester equations with a small second matrix.
 * @test
 * This function checks the @ref mess_direct_solvem, @ref mess_direct_solvemt, @ref mess_direct_solvemh functions for the sparse-dense Sylvester equation solver that
 * means it checks if Sylvester equations
 * <ul>
 * <li>  \f[ A X   +   X H   + M = 0 \f], \f[ A^T X     + X H^T     + M = 0 \f], \f[ A^H X     + X H^H     + M = 0 \f]
 * <li>  \f[ A X   + E X H + M = 0 \f],   \f[ A^T X     + E^T X H^T + M = 0 \f], \f[ A^T X     + E^H X H^H + M = 0 \f]
 * <li>  \f[ A X F + E X H + M = 0 \f],   \f[ A^T X F^T + E^T X H^T + M = 0 \f], \f[ A^T X F^H + E^H X H^H + M = 0 \f]
 * </ul>
 * are solved correctly with a solver created by the @ref mess_direct_create_sylvester_sparsedense.
 * \n
 * \f$ A, E \f$ is a sparse or dense matrix and \f$ F, H \f$ is a small and dense matrix.
 *
 * @author @dykstra
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

    int ret=0;
    const int verbose=1;
    double res2 = 0, res2_0 = 0;
    double tol = 1e-6;
    int Xdatatype_test_failed=0;
    mess_int_t seed = 1234;
    mess_matrix A=NULL, F=NULL, E=NULL, H=NULL, M=NULL, X=NULL, Temp=NULL;
    mess_direct sylv;
    mess_error_level = 1;

    mess_operation_t ops [] = {MESS_OP_NONE, MESS_OP_TRANSPOSE, MESS_OP_HERMITIAN};
    const int length_ops = sizeof(ops)/sizeof(mess_operation_t);
    int i_op;

    mess_datatype_t dts [] = {MESS_REAL, MESS_COMPLEX};
    const int length_dts = sizeof(dts)/sizeof(mess_datatype_t);
    int i_dtA, i_dtE, i_dtH, i_dtM, i_dtF;

    //storage type for sparse matrices
    mess_storage_t sts [] = {MESS_CSR, MESS_CSC, MESS_COORD};
    const int length_sts = sizeof(sts)/sizeof(mess_storage_t);
    int i_stA, i_stE;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if (!(3<= argc  && argc <= 5)){
        printf("usage: %s A.mtx H.mtx\n", argv[0]);
        printf("usage: %s A.mtx E.mtx H.mtx\n", argv[0]);
        printf("usage: %s A.mtx F.mtx E.mtx H.mtx\n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  init and read matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&H,&M,&X,&Temp);
    if(argc == 3){
        CALL(mess_matrix_read_formated(argv[1], A, MESS_CSR));
        CALL(mess_matrix_read_formated(argv[2], H, MESS_DENSE));
    }else if(argc == 4){
        MESS_INIT_MATRICES(&E);
        CALL(mess_matrix_read_formated(argv[1], A, MESS_CSR));
        CALL(mess_matrix_read_formated(argv[2], E, MESS_CSR));
        CALL(mess_matrix_read_formated(argv[3], H, MESS_DENSE));
    }else if (argc == 5){
        MESS_INIT_MATRICES(&E,&F);
        CALL(mess_matrix_read_formated(argv[1], A, MESS_CSR));
        CALL(mess_matrix_read_formated(argv[2], F, MESS_DENSE));
        CALL(mess_matrix_read_formated(argv[3], E, MESS_CSR));
        CALL(mess_matrix_read_formated(argv[4], H, MESS_DENSE));
    }

    /*-----------------------------------------------------------------------------
     *  create random rhs
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_rand_init(&seed));
    CALL(mess_matrix_rand_dense(M,A->rows, H->cols,MESS_REAL));

    /*-----------------------------------------------------------------------------
     *  test mess_direct_sylvester_sparsedense with different cases
     *-----------------------------------------------------------------------------*/

    /* iterate over all possible datatypes for F, if F is NULL this loop is encountered exactly ones*/
    for(i_dtF=F?0:length_dts-1; i_dtF<length_dts; ++i_dtF){

        /* iterate over all possible datatypes / storagetypes for E, if E is NULL this loop is encountered exactly ones*/
        for(i_dtE=E?0:length_dts-1, i_stE=E?0:length_sts-1; i_stE<length_sts; i_stE += (i_dtE+1)/length_dts, i_dtE = (i_dtE+1)%length_dts){

            /* iterate over all possible datatypes for A */
            for(i_dtA=0;i_dtA<length_dts;++i_dtA){

                /* iterate over all possible storagetypes for A */
                for(i_stA=0;i_stA<length_sts;++i_stA){

                    /* iterate over all possible datatypes for H, H must be dense */
                    for(i_dtH=0;i_dtH<length_dts;++i_dtH){

                        /* iterate over all possible datatypes for M, M must be dense */
                        for(i_dtM=0;i_dtM<length_dts;++i_dtM){

                            /*-----------------------------------------------------------------------------
                             *  convert Matrix A
                             *-----------------------------------------------------------------------------*/
                            CALL(mess_matrix_totype(A,dts[i_dtA]));
                            if(MESS_IS_COMPLEX(A)){
                                CALL(mess_matrix_scalec(1+I*1, A));
                            }
                            CALL(mess_matrix_convert(A,Temp,sts[i_stA]));
                            CALL(mess_matrix_copy(Temp,A));

                            /*-----------------------------------------------------------------------------
                             *  convert Matrix H
                             *-----------------------------------------------------------------------------*/
                            CALL(mess_matrix_totype(H,dts[i_dtH]));
                            if(MESS_IS_COMPLEX(H)){
                                CALL(mess_matrix_scalec(1+I*1, H));
                            }

                            /*-----------------------------------------------------------------------------
                             *  convert Matrix M
                             *-----------------------------------------------------------------------------*/
                            CALL(mess_matrix_totype(M,dts[i_dtM]));
                            if(MESS_IS_COMPLEX(M)){
                                CALL(mess_matrix_scalec(1+I*1, M));
                            }
                            CALL(mess_matrix_norm2(M, &res2_0));

                            /*-----------------------------------------------------------------------------
                             *  convert Matrix E if E is not NULL
                             *-----------------------------------------------------------------------------*/
                            if(E){
                                CALL(mess_matrix_totype(E,dts[i_dtE]));
                                if(MESS_IS_COMPLEX(E)){
                                    CALL(mess_matrix_scalec(1+I*1, E));
                                }
                                CALL(mess_matrix_convert(E,Temp,sts[i_stE]));
                                CALL(mess_matrix_copy(Temp,E));
                            }

                            /*-----------------------------------------------------------------------------
                             *  convert Matrix F if F is not NULL
                             *-----------------------------------------------------------------------------*/
                            if(F){
                                CALL(mess_matrix_totype(F,dts[i_dtF]));
                                if(MESS_IS_COMPLEX(F)){
                                    CALL(mess_matrix_scalec(1+I*1, F));
                                }
                            }

                            /*-----------------------------------------------------------------------------
                             *  create sparse dense solver and solve
                             *-----------------------------------------------------------------------------*/
                            CALL(mess_direct_init(&sylv));
                            if(E){
                                if(F){
                                    CALL(mess_direct_create_sylvester_sparsedense(A, F, E, H, sylv));
                                }else{
                                    CALL(mess_direct_create_sylvester_sparsedense(A, NULL, E, H, sylv));
                                }
                            }else{
                                CALL(mess_direct_create_sylvester_sparsedense(A, NULL, NULL, H, sylv));
                            }

                            /*-----------------------------------------------------------------------------
                             * solve sparse-dense Sylvester equation and  print residual
                             *-----------------------------------------------------------------------------*/
                            /* iterate over all possible operation types*/
                            for(i_op=0;i_op<length_ops;++i_op){

                                /*-----------------------------------------------------------------------------
                                 *  solve sylvester equation and compute residual
                                 *-----------------------------------------------------------------------------*/
                                CALL(mess_direct_solvem(ops[i_op],sylv,M,X));
                                if(E){
                                    if(F){
                                        CALL(mess_direct_generalized_sylvester_res2(ops[i_op],A,F,E,H,M,X,&res2));
                                    }else{
                                        CALL(mess_direct_sylvestersg_res2(ops[i_op],A,E,H,M,X,&res2));
                                    }
                                }else{
                                    CALL(mess_direct_sylvester_res2(ops[i_op],A,H,M,X,&res2));
                                }

                                /*-----------------------------------------------------------------------------
                                 *  print infos and check residual and datatype
                                 *-----------------------------------------------------------------------------*/
                                if(verbose){
                                    test_sylvester_info(sylv, ops[i_op], A, F,  E, H, M, res2, res2_0);
                                }

                                /*-----------------------------------------------------------------------------
                                 *  check data_type and residual
                                 *-----------------------------------------------------------------------------*/
                                if(MESS_IS_REAL(A) && MESS_IS_REAL(H) && MESS_IS_REAL(M)){
                                    if(!F){
                                        Xdatatype_test_failed =(!E || MESS_IS_REAL(E))?(MESS_IS_REAL(X)? 0:1):(MESS_IS_COMPLEX(X)? 0:1);
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
                                    printf("Matrix A:\n"); mess_matrix_printinfo(A);
                                    if(F){
                                        printf("Matrix F:\n"); mess_matrix_printinfo(F);
                                    }
                                    if(E){
                                        printf("Matrix E:\n"); mess_matrix_printinfo(E);
                                    }
                                    printf("Matrix H:\n"); mess_matrix_printinfo(H);
                                    printf("Matrix M:\n"); mess_matrix_printinfo(M);
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

    if(!mess_direct_create_sylvester_sparsedense(A,H,NULL,H,sylv)){
        printf("didn't catch wrong case AXF + XH + M = 0.\n");
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&H,&M,&X,&Temp);
    if(F){MESS_CLEAR_MATRICES(&F);}
    if(E){MESS_CLEAR_MATRICES(&E);}

    return 0;
}

