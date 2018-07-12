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
 * @file tests/matrix/check_matrix_null.c
 * @brief Check the orthogonal basis construction of the nullspace.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_null function defined in orth.c that means it checks if an orthogonal basis
 * of the nullspace is constructed correctly.
 *
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>
#include "mess/mess.h"

#include "../call_macro.h"


#define CHECK_MESS_MATRIX_NULL(Q,ZNULL,ZORTH,TEMP1,TEMP2,EYE,RANK,NRM,ERR,EPS){                 \
    CALL(mess_matrix_null( Q, ZNULL ));                                                         \
    CALL(mess_matrix_orth(Q,ZORTH));                                                            \
    /*check if dim(ker(Q)) + dim(range(Q)) == Q->cols*/                                         \
    if(ZNULL->cols+ZORTH->cols==Q->cols){                                                       \
        /*check Q*ZNULL==0 if dim(ker(Q))>0 */                                                  \
        if(ZNULL->cols>0){                                                                      \
            CALL( mess_matrix_multiply(MESS_OP_NONE, Q, MESS_OP_NONE, ZNULL , TEMP1));                                  \
            CALL(mess_matrix_norm1(TEMP1,&NRM));                                                \
            if(NRM>EPS){                                                                        \
                printf("Failed! Basis of nullspace seems to be wrong\n");                       \
                printf("Diff: %e\n",NRM);                                                       \
                ++ERR;                                                                          \
            }                                                                                   \
            \
            /*compute RANK of ZNULL, check if colums are linear independent */                  \
            CALL(mess_matrix_rank(ZNULL,&RANK));                                                \
            if(RANK!=ZNULL->cols){                                                              \
                printf("Failed! Basis of nullspace seems to be linear dependent\n");            \
                ++ERR;                                                                          \
            }                                                                                   \
            \
            /*check if basis of nullspace is orthogonal*/                                       \
            CALL(mess_matrix_copy(ZNULL,TEMP1));                                                \
            CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE, ZNULL, MESS_OP_NONE, TEMP1, TEMP2));      \
            CALL(mess_matrix_eye(EYE,TEMP2->rows,TEMP2->cols,MESS_DENSE));                      \
            CALL(mess_matrix_diffnorm(EYE,TEMP2,&NRM));                                             \
            if(NRM>EPS){                                                                        \
                printf("Failed! Basis of nullspace seems to be not orthogonal\n");              \
                printf("Diff: %e\n",NRM);                                                       \
                ++ERR;                                                                          \
            }                                                                                   \
        }                                                                                       \
    }else{                                                                                      \
        printf("Failed\n");                                                                     \
        printf("Either number of basis vector of nullspace or range is wrong\n");               \
        ++ERR;                                                                                  \
        \
    }                                                                                           \
}


int main ( int argc, char ** argv) {
    mess_init();
    int rows=10, cols=5,i=0, ret=0, err=0;
    double p =1 ;                   //density for random matrix generation negelected in dense case
    double eps = 1e-11;
    double nrm=0;

    mess_int_t rank =0;

    mess_matrix Q, Znull, Zorth,Temp1, Temp2, Eye, Qcol;


    //init real matrices
    CALL(mess_matrix_init(&Q));
    CALL(mess_matrix_init(&Znull));
    CALL(mess_matrix_init(&Zorth))

        //init temp matrices;
        CALL(mess_matrix_init(&Temp1));
    CALL(mess_matrix_init(&Temp2));
    CALL(mess_matrix_init(&Eye));

    CALL(mess_matrix_init(&Qcol));



    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    //load  matrices
    CALL(mess_matrix_rand(Q,rows,cols,MESS_DENSE,MESS_REAL,p));


    /*-----------------------------------------------------------------------------
     *  Test
     *-----------------------------------------------------------------------------*/

    CHECK_MESS_MATRIX_NULL(Q,Znull,Zorth,Temp1,Temp2,Eye,rank,nrm,err,eps);

    for(i=0;i<cols;++i){
        mess_matrix_colsub(Q,i,i,Qcol);
        mess_matrix_addcols(Q,Qcol);
        CHECK_MESS_MATRIX_NULL(Q,Znull,Zorth,Temp1,Temp2,Eye,rank,nrm,err,eps);
    }


    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/
    mess_matrix_clear(&Q);
    mess_matrix_clear(&Zorth);
    mess_matrix_clear(&Znull);
    mess_matrix_clear(&Temp1);
    mess_matrix_clear(&Temp2);
    mess_matrix_clear(&Eye);
    mess_matrix_clear(&Qcol);

    return (err>0)?(1):(0);
}
