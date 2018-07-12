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
 * @addtogroup test_format
 * @{
 * @file tests/formats/check_mvpcall_matrix.c
 * @brief Check the @ref mess_mvpcall_matrix function.
 * @author  @mbehr
 * @test
 *
 * This function checks the @ref mess_mvpcall_matrix function defined in mvpcall.c that means it uses
 * @ref mess_matrix_mvp to compare the results.
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

int check_mvp_call_matrix(mess_matrix mat, mess_operation_t op1, mess_operation_t op2,mess_vector v1, mess_vector v2){
    int ret =0;
    switch(op1){
        case MESS_OP_NONE:
            CALL(mess_matrix_mvp(op2, mat, v1, v2));
            break;
        case MESS_OP_TRANSPOSE:
            switch(op2){
                case MESS_OP_NONE:
                    CALL(mess_matrix_mvp(MESS_OP_TRANSPOSE, mat, v1,v2));
                    break;
                case MESS_OP_TRANSPOSE:
                    CALL(mess_matrix_mvp(MESS_OP_NONE, mat, v1,v2));
                    break;
                case MESS_OP_HERMITIAN:
                    CALL(mess_matrix_mvp(MESS_OP_NONE,mat,v1,v2));
                    CALL(mess_vector_conj(v2));
                    break;
            }
            break;
        case MESS_OP_HERMITIAN:
            switch(op2){
                case MESS_OP_NONE:
                    CALL(mess_matrix_mvp(MESS_OP_HERMITIAN,mat,v1,v2));
                    break;
                case MESS_OP_TRANSPOSE:
                    CALL(mess_matrix_mvp(MESS_OP_NONE, mat, v1,v2));
                    CALL(mess_vector_conj(v2));
                    break;
                case MESS_OP_HERMITIAN:
                    CALL(mess_matrix_mvp(MESS_OP_NONE, mat, v1,v2));
                    break;
            }
            break;
    }
    return ret;
}


int main ( int argc, char ** argv) {

    mess_init();
    double eps=sqrt(mess_eps()), diff, p = 0.4;
    int ret=0, err=0;
    mess_matrix mat, matc, temp;
    mess_vector v1,v2,v3;
    mess_storage_t  storage_types [] = {MESS_DENSE, MESS_CSR, MESS_CSC, MESS_COORD};
    mess_datatype_t data_types [] = {MESS_REAL, MESS_COMPLEX};
    mess_operation_t op_types [] = {MESS_OP_NONE, MESS_OP_TRANSPOSE, MESS_OP_HERMITIAN};

    mess_int_t st=0, dt=0, op1=0, op2=0;
    mess_int_t rows, cols;
    mess_mvpcall mvp, mvpc;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if (argc !=3) {
        printf("usage: %s rows cols\n", argv[0]);
        return 1;
    }

    rows = atoi(argv[1]);
    cols = atoi(argv[2]);

    /*-----------------------------------------------------------------------------
     *  generate matrix
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&mat,&matc,&temp);
    CALL(mess_matrix_rand(mat, rows, cols,  MESS_CSR, MESS_REAL,p));
    CALL(mess_matrix_rand(matc, rows, cols, MESS_CSR, MESS_COMPLEX, p));

    /*-----------------------------------------------------------------------------
     *  init vectors
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&v1,&v2,&v3);
    CALL(mess_vector_alloc(v1, 1, MESS_REAL));
    CALL(mess_vector_alloc(v2, 1, MESS_REAL));
    CALL(mess_vector_alloc(v3, 1, MESS_REAL));

    /*-----------------------------------------------------------------------------
     *  check different storage and datatypes
     *-----------------------------------------------------------------------------*/

    for(op1=0;op1<3;++op1){
        for(op2=0;op2<3;++op2){
            for(st=0;st<3;++st){
                for(dt=0;dt<2;++dt){

                    /*-----------------------------------------------------------------------------
                     *  convert to datatype and storage type
                     *-----------------------------------------------------------------------------*/
                    CALL(mess_matrix_convert(mat, temp, storage_types[st]));
                    CALL(mess_matrix_copy(temp,mat));

                    CALL(mess_matrix_convert(matc, temp, storage_types[st]));
                    CALL(mess_matrix_copy(temp,matc));

                    /*-----------------------------------------------------------------------------
                     *  generate mvpcall operator
                     *-----------------------------------------------------------------------------*/
                    CALL(mess_mvpcall_matrix(&mvp, op_types[op1], mat));
                    CALL(mess_mvpcall_matrix(&mvpc, op_types[op1], matc));

                    /*-----------------------------------------------------------------------------
                     *  generate random vector of data_type and needed size
                     *-----------------------------------------------------------------------------*/
                    CALL(mess_vector_totype(v1,data_types[dt]));
                    switch(op_types[op1]){
                        case MESS_OP_NONE:
                            CALL(mess_vector_resize(v1,(op_types[op2]==MESS_OP_NONE)?cols:rows));
                            break;
                        case MESS_OP_TRANSPOSE:
                            CALL(mess_vector_resize(v1,(op_types[op2]==MESS_OP_NONE)?rows:cols));
                            break;
                        case MESS_OP_HERMITIAN:
                            CALL(mess_vector_resize(v1,(op_types[op2]==MESS_OP_NONE)?rows:cols));
                            break;
                    }
                    CALL(mess_vector_rand(v1));

                    /*-----------------------------------------------------------------------------
                     *  perform multiplication with mvp call and mess_matrix_mvp (MESS_REAL operator)
                     *-----------------------------------------------------------------------------*/
                    CALL(mess_mvpcall_apply(mvp,op_types[op2],v1,v2));
                    CALL(check_mvp_call_matrix(mat, op_types[op1], op_types[op2], v1, v3));

                    /*-----------------------------------------------------------------------------
                     *  compare results
                     *-----------------------------------------------------------------------------*/
                    CALL(mess_vector_diffnorm(v2,v3,&diff));
                    if(diff>eps){
                        printf("FAILED REAL operator mvp:\n");
                        mess_matrix_printinfo(matc);
                        mess_vector_printinfo(v1);
                        printf("op1=%s\n", mess_operation_t_str(op_types[op1]));
                        printf("op2=%s\n", mess_operation_t_str(op_types[op2]));
                        printf("DIFF=%e\n",diff);
                        return 1;
                    }

                    /*-----------------------------------------------------------------------------
                     *  perform multiplication with mvp call and mess_matrix_mvp (MESS_COMPLEX operator)
                     *-----------------------------------------------------------------------------*/
                    CALL(mess_mvpcall_apply(mvpc,op_types[op2],v1,v2));
                    CALL(check_mvp_call_matrix(matc, op_types[op1], op_types[op2], v1, v3));

                    /*-----------------------------------------------------------------------------
                     *  compare results
                     *-----------------------------------------------------------------------------*/
                    CALL(mess_vector_diffnorm(v2,v3,&diff));
                    if(diff>eps){
                        printf("FAILED COMPLEX operator mvpc:\n");
                        mess_matrix_printinfo(matc);
                        mess_vector_printinfo(v1);
                        printf("op1=%s\n", mess_operation_t_str(op_types[op1]));
                        printf("op2=%s\n", mess_operation_t_str(op_types[op2]));
                        printf("DIFF=%e\n",diff);
                        return 1;
                    }

                    /*-----------------------------------------------------------------------------
                     *  clear operators
                     *-----------------------------------------------------------------------------*/
                    CALL(mess_mvpcall_clear(&mvp));
                    CALL(mess_mvpcall_clear(&mvpc));

                }
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&mat, &matc, &temp);
    MESS_CLEAR_VECTORS(&v1,&v2,&v3);

    return (err>0)?(1):(0);
}

