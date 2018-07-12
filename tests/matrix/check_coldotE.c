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
 * @file tests/matrix/check_coldotE.c
 * @brief Check the \ref mess_matrix_coldotE function.
 * @author  @mbehr
 * @test
 * This function check the @ref mess_matrix_coldotE function defined in
 * colops.c, that means it checks if the dot product
 * \f[Q(:,col_1)^T*E*Q(:,col_2)\f]
 * is computed correctly for a given matrix \f$ Q \f$ and a
 * symmetric positive definite matrix \f$ E \f$ and columns \f$ col_1\f$ and \f$ col_2 \f$.
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

#define CHECKCOLVDOTE(Q,E,COL1,COL2,SOL,DOT,ERR,EPS,DIFF){                  \
    CALL(mess_matrix_coldotE(Q, E, COL1, COL2, &DOT ))                  \
    DIFF=fabs(SOL->values[COL1*SOL->ld+COL2]-DOT);                      \
    printf("Dot: %e\n",DOT);                                                \
    printf("Diff%e\n\n",DIFF);                                          \
    if(DIFF>EPS){++ERR;}                                                \
}

int main ( int argc, char ** argv) {
    mess_init();
    int rows=6, cols=12,i=0, ret=0, err=0;
    double p = 0.6;                     //density for random matrix generation
    double eps = 1e-12;
    double dot=0,diff=0,diagonaladd=0;  //needed for creating a symmetric positiv definit Matrix E
    mess_int_t col1,col2;

    mess_matrix Q, Qcsc, Qcsr, E, Ecsc, Ecsr, temp, sol;

    //init Q matrices
    CALL(mess_matrix_init(&Q));
    CALL(mess_matrix_init(&Qcsc));
    CALL(mess_matrix_init(&Qcsr));

    //init E matrices
    CALL(mess_matrix_init(&E));
    CALL(mess_matrix_init(&Ecsc));
    CALL(mess_matrix_init(&Ecsr));

    //init additional Matrices
    CALL(mess_matrix_init(&temp));
    CALL(mess_matrix_init(&sol));

    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    //load real matrices
    CALL(mess_matrix_rand(Q,rows,cols,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_convert(Q,Qcsr,MESS_CSR));
    CALL(mess_matrix_convert(Q,Qcsc,MESS_CSC));

    //first create a positiv definite matrix <-> symmetric diagonaldominant matrix
    CALL(mess_matrix_rand(E,rows,rows,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_ctranspose(E,temp));
    CALL(mess_matrix_add(0.5,temp,0.5,E)); //E is symmetric and all entries e holds 0<=e<=1
    diagonaladd=Q->rows*Q->cols+1;
    for (i = 0; i < E->cols; ++i) {E->values[i*E->ld+i]+=diagonaladd;} //E is symmetric and diagonaldominant-> E is positiv definite

    CALL(mess_matrix_convert(E,Ecsr,MESS_CSR));
    CALL(mess_matrix_convert(E,Ecsc,MESS_CSC));

    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/

    //compute solution for comparison
    mess_matrix_multiply(MESS_OP_TRANSPOSE, Q, MESS_OP_NONE, E, temp);
    mess_matrix_multiply(MESS_OP_NONE, temp, MESS_OP_NONE, Q, sol); //sol(i,j)=Q(:,i)'E*Q(:,j)

    //test function
    for(col1=0;col1<cols;++col1){
        for(col2=0;col2<cols;++col2){
            CHECKCOLVDOTE(Q,E,col1,col2,sol,dot,err,eps,diff);
            CHECKCOLVDOTE(Qcsc,E,col1,col2,sol,dot,err,eps,diff);
            CHECKCOLVDOTE(Qcsr,E,col1,col2,sol,dot,err,eps,diff);
            CHECKCOLVDOTE(Q,Ecsc,col1,col2,sol,dot,err,eps,diff);
            CHECKCOLVDOTE(Qcsc,Ecsc,col1,col2,sol,dot,err,eps,diff);
            CHECKCOLVDOTE(Qcsr,Ecsc,col1,col2,sol,dot,err,eps,diff);
            CHECKCOLVDOTE(Q,Ecsr,col1,col2,sol,dot,err,eps,diff);
            CHECKCOLVDOTE(Qcsc,Ecsr,col1,col2,sol,dot,err,eps,diff);
            CHECKCOLVDOTE(Qcsr,Ecsr,col1,col2,sol,dot,err,eps,diff);
        }
    }


    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/


    mess_matrix_clear(&Q);
    mess_matrix_clear(&Qcsr);
    mess_matrix_clear(&Qcsc);
    mess_matrix_clear(&E);
    mess_matrix_clear(&Ecsr);
    mess_matrix_clear(&Ecsc);
    mess_matrix_clear(&temp);
    mess_matrix_clear(&sol);


    return (err>0)?(1):(0);
}
