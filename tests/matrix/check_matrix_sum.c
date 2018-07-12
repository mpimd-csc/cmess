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
 * @file tests/matrix/check_matrix_sum.c
 * @brief Check the column and row sum computation of a matrix.
 * @author  @mbehr
 * @test
 * This function checks the mess_matrix_colsums and @ref mess_matrix_rowsums function defined in colrowsums.c that means
 * it checks if a column and row sum of a given matrix is computed correctly.
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


//checks the mess_matrix_rowsum against the mess_matrix_colsum function
#define CHECKMATRIXSUM(A,TEMP,VSUM1,VSUM2,ERR,EPS,DIFF){                    \
    CALL(mess_matrix_rowsums(A,VSUM1));                                 \
    CALL(mess_matrix_ctranspose(A,TEMP));                               \
    CALL(mess_matrix_colsums(TEMP,VSUM2));                              \
    CALL(mess_vector_diffnorm(VSUM1,VSUM2,&DIFF));                  \
    if(DIFF>10*EPS){                                                \
        printf("Failed:\n");                                        \
        printf("A: %s\n", #A);                                      \
        printf("TEMP: %s\n", #TEMP);                                \
        printf("VSUM1: %s\n", #VSUM1);                              \
        printf("VSUM2: %s\n", #VSUM2);                              \
        printf("Diff: %e\n",DIFF);                                  \
        ++ERR;                                                      \
    }                                                               \
}



int main ( int argc, char ** argv) {
    mess_init();
    int rows=10, cols=5, ret=0, err=0;
    double p = 1;                   //density for random matrix generation negelected in dense case
    double eps = 1e-12;
    double diff=0;

    mess_matrix Adenser, Acscr, Acsrr, Adensec, Acscc, Acsrc, Temp1;
    mess_vector rowsum1, rowsum2;

    //init matrices
    CALL(mess_matrix_init(&Adenser));
    CALL(mess_matrix_init(&Acscr));
    CALL(mess_matrix_init(&Acsrr));
    CALL(mess_matrix_init(&Adensec));
    CALL(mess_matrix_init(&Acsrc));
    CALL(mess_matrix_init(&Acscc));
    CALL(mess_matrix_init(&Temp1));

    MESS_INIT_VECTORS(&rowsum1,&rowsum2);
    CALL(mess_vector_alloc(rowsum1,rows,MESS_REAL));
    CALL(mess_vector_alloc(rowsum2,rows,MESS_REAL));


    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    //load real matrices
    CALL(mess_matrix_rand(Adenser,rows,cols,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_convert(Adenser,Acscr,MESS_CSC));
    CALL(mess_matrix_convert(Adenser,Acsrr,MESS_CSR));

    //load complex matrices
    CALL(mess_matrix_rand(Adensec,rows,cols,MESS_DENSE,MESS_COMPLEX,p));
    CALL(mess_matrix_convert(Adensec,Acscc,MESS_CSC));
    CALL(mess_matrix_convert(Adensec,Acsrc,MESS_CSR));

    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/

    CHECKMATRIXSUM(Adenser,Temp1,rowsum1,rowsum2,err,eps,diff);
    CHECKMATRIXSUM(Acscr,Temp1,rowsum1,rowsum2,err,eps,diff);
    CHECKMATRIXSUM(Acsrr,Temp1,rowsum1,rowsum2,err,eps,diff);
    CHECKMATRIXSUM(Adensec,Temp1,rowsum1,rowsum2,err,eps,diff);
    CHECKMATRIXSUM(Acscc,Temp1,rowsum1,rowsum2,err,eps,diff);
    CHECKMATRIXSUM(Acsrc,Temp1,rowsum1,rowsum2,err,eps,diff);

    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_clear(&Adenser));
    CALL(mess_matrix_clear(&Adensec));
    CALL(mess_matrix_clear(&Acscr));
    CALL(mess_matrix_clear(&Acscc));
    CALL(mess_matrix_clear(&Acsrr));
    CALL(mess_matrix_clear(&Acsrc));
    CALL(mess_matrix_clear(&Temp1));

    CALL(mess_vector_clear(&rowsum1));
    CALL(mess_vector_clear(&rowsum2));

    return (err>0)?(1):(0);
}


