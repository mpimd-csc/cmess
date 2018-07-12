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
 * @file tests/matrix/check_orth.c
 * @brief Check the orthogonal basis construction using SVD.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_orth function defined in orth.c that means it checks if an orthogonal
 * basis of a matrix is computed using singular value decomposition.
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

#define CHECKORTH(MAT,Q,EYE,TEMP1,TEMP2,ERR,EPS,DIFF){                                      \
    mess_matrix_orth(MAT,Q);                                                            \
    mess_matrix_copy(Q,TEMP1);                                                          \
    if(MESS_IS_COMPLEX(MAT)){                                                           \
        CALL( mess_matrix_multiply(MESS_OP_HERMITIAN, TEMP1, MESS_OP_NONE, Q, TEMP2));          \
    }else{                                                                              \
        CALL( mess_matrix_multiply(MESS_OP_HERMITIAN, TEMP1, MESS_OP_NONE, Q, TEMP2));          \
    }                                                                                   \
    CALL(mess_matrix_diffnorm(TEMP2,EYE,&DIFF));                                            \
    if(DIFF>EPS){                                                                       \
        ++ERR;                                                                          \
        printf("mess_matrix_mgsadd failed \n");                                         \
        printf("Error:%e\n",DIFF);                                                      \
    }                                                                                   \
}






int main ( int argc, char ** argv) {
    mess_init();
    int err=0, ret=0;
    double eps = 1e-12;
    double diff=0;

    mess_int_t cols = 10;
    mess_int_t rows = 10;

    mess_matrix Matr, Matc, Qr, Qc, Eyer, Temp1, Temp2;

    //init  matrices
    CALL(mess_matrix_init(&Matr));
    CALL(mess_matrix_init(&Matc));
    CALL(mess_matrix_init(&Qr));
    CALL(mess_matrix_init(&Qc));
    CALL(mess_matrix_init(&Eyer));
    CALL(mess_matrix_init(&Temp1));
    CALL(mess_matrix_init(&Temp2));



    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    //load  matrices
    CALL(mess_matrix_rand(Matr,rows,cols,MESS_DENSE,MESS_REAL,1));
    CALL(mess_matrix_rand(Matc,rows,cols,MESS_DENSE,MESS_COMPLEX,1));
    CALL(mess_matrix_eye(Eyer,rows,cols,MESS_DENSE));


    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/
    CHECKORTH(Matr,Qr,Eyer,Temp1,Temp2,err,eps,diff);
    CHECKORTH(Matc,Qc,Eyer,Temp1,Temp2,err,eps,diff);

    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/

    mess_matrix_clear(&Matr);
    mess_matrix_clear(&Matc);
    mess_matrix_clear(&Qr);
    mess_matrix_clear(&Qc);
    mess_matrix_clear(&Eyer);
    mess_matrix_clear(&Temp1);
    mess_matrix_clear(&Temp2);

    return (err>0)?(1):(0);
}
