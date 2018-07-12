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
 * @file tests/matrix/check_rowsub.c
 * @brief Check the row submatrix computation.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_rowsub function defined in sub.c that means it checks if the returned row
 * submatrix is computed correctly from a given matrix.
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

#define CHECKROWSUB(A,TEMP,COL1,COL2,SOL,ERR,EPS,DIFF){                     \
    mess_matrix_rowsub(A, COL1, COL2, TEMP);                                \
    mess_matrix_diffnorm(SOL,TEMP,&DIFF);                                   \
    /*mess_matrix_print(TEMP);*/                                            \
    printf("%e\n",DIFF);                                                    \
    if(DIFF>EPS){++ERR;}                                                    \
}

int main ( ) {
    mess_init();
    int rows=6, cols=12, ret=0, err=0;
    double p = 0.6;                     //density for random matrix generation
    double eps = 1e-12;
    double diff=0;
    mess_int_t ind1,ind2;

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

    for(ind1=0;ind1<rows;++ind1){
        for(ind2=ind1;ind2<rows;++ind2){
            //compute real solution for comparison
            CALL(mess_matrix_rowsub(Ar,ind1,ind2,solr));
            /*mess_matrix_print(solr);*/

            //compute complex solution for comparison
            CALL(mess_matrix_rowsub(Ac,ind1,ind2,solc));

            //start testing
            CHECKROWSUB(Acsrr, temp, ind1, ind2, solr, err, eps, diff);
            CHECKROWSUB(Acscr, temp, ind1, ind2, solr, err, eps, diff);
            CHECKROWSUB(Acscc, temp, ind1, ind2, solc, err, eps, diff);
            CHECKROWSUB(Acsrc, temp, ind1, ind2, solc, err, eps, diff);
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
