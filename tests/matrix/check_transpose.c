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
 * @file tests/matrix/check_transpose.c
 * @brief Check the transpose function.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_ctranspose function defined in transpose.c that means it checks if the transpose
 * of a matrix is computed correctly.
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


#define CHECKTRANSPOSE(A,SOL,DIFF,EPS, ERR){                        \
    mess_matrix_ctranspose(A,TEMP);                                 \
    mess_matrix_diffnorm(TEMP,SOL,&(DIFF));                         \
    if((DIFF)>(EPS))    {                                           \
        ++(ERR);                                                    \
        mess_matrix_printinfo(A);                                   \
    }                                                               \
}


double matrix_A[9]  = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
double matrix_AT[9] = { 1, 4, 7, 2, 5, 8, 3, 6, 9};



int main ( ) {
    mess_init();

    double diff=.0, eps=mess_eps();
    int ret=0, err=0;
    // int rows=30,cols=40;

    //prepare dense matrices
    mess_matrix Adenser, Adensec, Acscr, Acscc, Acsrr, Acsrc, Acoordr, Acoordc, SOL,TEMP;

    //init matrices
    CALL(mess_matrix_init(&SOL));
    CALL(mess_matrix_init(&TEMP));

    CALL(mess_matrix_init(&Adenser));
    CALL(mess_matrix_init(&Adensec));

    CALL(mess_matrix_init(&Acsrr));
    CALL(mess_matrix_init(&Acsrc));

    CALL(mess_matrix_init(&Acscr));
    CALL(mess_matrix_init(&Acscc));

    CALL(mess_matrix_init(&Acoordr));
    CALL(mess_matrix_init(&Acoordc));



    //generate matrices
    CALL(mess_matrix_dense_from_farray(SOL,3,3,3,matrix_AT,NULL));

    CALL(mess_matrix_dense_from_farray(Adenser,3,3,3,matrix_A,NULL));
    CALL(mess_matrix_convert(Adenser,Acsrr,MESS_CSR));
    CALL(mess_matrix_convert(Adenser,Acscr,MESS_CSC));
    CALL(mess_matrix_convert(Adenser,Acoordr,MESS_COORD));


    CALL(mess_matrix_convert(Adenser,Adensec,MESS_DENSE));  CALL(mess_matrix_tocomplex(Adensec));
    CALL(mess_matrix_convert(Adenser,Acsrc,MESS_CSR));      CALL(mess_matrix_tocomplex(Acsrc));
    CALL(mess_matrix_convert(Adenser,Acscc,MESS_CSC));      CALL(mess_matrix_tocomplex(Acscc));
    CALL(mess_matrix_convert(Adenser,Acoordc,MESS_CSC));    CALL(mess_matrix_tocomplex(Acoordc));


    //test cases
    CHECKTRANSPOSE(Adenser,SOL,diff,eps, err);
    CHECKTRANSPOSE(Adensec,SOL,diff,eps, err);
    CHECKTRANSPOSE(Acsrr,SOL,diff,eps, err);
    CHECKTRANSPOSE(Acsrc,SOL,diff,eps, err);
    CHECKTRANSPOSE(Acscr,SOL,diff,eps, err);
    CHECKTRANSPOSE(Acscc,SOL,diff,eps, err);
    CHECKTRANSPOSE(Acoordr,SOL,diff,eps, err);
    CHECKTRANSPOSE(Acoordc,SOL,diff,eps, err);



    mess_matrix_clear(&Adenser);
    mess_matrix_clear(&Adensec);
    mess_matrix_clear(&Acscr);
    mess_matrix_clear(&Acscc);
    mess_matrix_clear(&Acsrr);
    mess_matrix_clear(&Acsrc);
    mess_matrix_clear(&Acoordr);
    mess_matrix_clear(&Acoordc);
    mess_matrix_clear(&SOL);
    mess_matrix_clear(&TEMP);

    return (err>0)?(1):(0);
}

