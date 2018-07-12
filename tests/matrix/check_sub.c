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
 * @file tests/matrix/check_sub.c
 * @brief Check the sub matrix computation.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_sub function defined in sub.c that means it checks if the returned sub matrix
 * of a given matrix is correct.
 *
 * @}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mess/mess.h"

#include "../call_macro.h"

#define CHECKSUB(A,TEMP,ROW1,ROW2,COL1,COL2,SOL,ERR,EPS,DIFF){              \
    mess_matrix_sub(A,ROW1,ROW2, COL1,COL2, TEMP);                          \
    mess_matrix_diffnorm(SOL,TEMP,&DIFF);                                   \
    /*mess_matrix_print(TEMP);*/                                            \
    printf("%e\n",DIFF);                                                    \
    if(DIFF>EPS){++ERR;}                                                    \
}

int main ( ) {
    mess_init();
    int rows=6, cols=12,ret=0, err=0;
    double p = 0.6;                     //density for random matrix generation
    double eps = 1e-12;
    double diff=0;
    mess_int_t row1, row2, col1, col2;

    mess_matrix Ar, Acscr, Acsrr, Ac, Acscc, Acsrc, temp, solr,solc;

    //init real matrices
    CALL(mess_matrix_init(&Ar));
    CALL(mess_matrix_init(&Acscr));
    CALL(mess_matrix_init(&Acsrr));

    //init complex matrices
    CALL(mess_matrix_init(&Ac));
    CALL(mess_matrix_init(&Acscc));
    CALL(mess_matrix_init(&Acsrc));

    //init additional Matrices
    CALL(mess_matrix_init(&temp));
    CALL(mess_matrix_init(&solr));
    CALL(mess_matrix_init(&solc));

    //load real matrices
    CALL(mess_matrix_rand(Ar,rows,cols,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_convert(Ar,Acsrr,MESS_CSR));
    CALL(mess_matrix_convert(Ar,Acscr,MESS_CSC));

    //load complex matrices
    CALL(mess_matrix_rand(Ac,rows,cols,MESS_DENSE,MESS_COMPLEX,p));
    CALL(mess_matrix_convert(Ac,Acsrc,MESS_CSR));
    CALL(mess_matrix_convert(Ac,Acscc,MESS_CSC));
    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/

    for(row1=0;row1<rows;++row1){
        for(row2=row1;row2<rows;++row2){
            for(col1=0;col1<cols;++col1){
                for(col2=col1;col2<cols;++col2){

                    //  row1=1; row2=3; col1=0;col2=3;
                    //compute real solution for comparison
                    CALL(mess_matrix_rowsub(Ar,row1,row2,temp));
                    CALL(mess_matrix_colsub(temp,col1,col2,solr));

                    //compute complex solution for comparison
                    CALL(mess_matrix_rowsub(Ac,row1,row2,temp));
                    CALL(mess_matrix_colsub(temp,col1,col2,solc));

                    //start testing
                    CHECKSUB(Ar, temp, row1, row2, col1, col2, solr, err, eps, diff);
                    CHECKSUB(Ac, temp, row1, row2, col1, col2, solc, err, eps, diff);
                    CHECKSUB(Acsrr, temp, row1, row2, col1, col2, solr, err, eps, diff);
                    CHECKSUB(Acscr, temp, row1, row2, col1, col2, solr, err, eps, diff);
                    CHECKSUB(Acscc, temp, row1, row2, col1, col2, solc, err, eps, diff);
                    CHECKSUB(Acsrc, temp, row1, row2, col1, col2, solc, err, eps, diff);
                }
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/

    mess_matrix_clear(&Ar);
    mess_matrix_clear(&Acsrr);
    mess_matrix_clear(&Acscr);
    mess_matrix_clear(&Ac);
    mess_matrix_clear(&Acsrc);
    mess_matrix_clear(&Acscc);
    mess_matrix_clear(&solr);
    mess_matrix_clear(&solc);
    mess_matrix_clear(&temp);

    return (err>0)?(1):(0);
}
