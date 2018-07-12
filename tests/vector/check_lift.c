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
 * @file tests/vector/check_lift.c
 * @brief Check the \ref mess_matrix_lift and \ref mess_vector_lift function.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_lift function defined in lift.c and the mess_vector_lift function defined in
 * vector.c. \n
 * It checks if the addition of a zero block or a zero vector at the end of a marix or a vector is computed correctly.
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


#define CHECKVECTORLIFT(X,LIFT,XLIFT,XZEROS,XTMP,DIFF,EPS,ERR)          \
    CALL(mess_vector_lift((X),(LIFT),(XLIFT)));                     \
CALL(mess_vector_resize((XZEROS),(LIFT)));                      \
CALL(mess_vector_zeros((XZEROS)));                              \
CALL(mess_vector_cat((X),(XZEROS),(XTMP)));                     \
CALL(mess_vector_diffnorm((XLIFT),(XTMP),&(DIFF)));             \
if((DIFF)>(EPS)){                                               \
    printf("mess_vector_lift failed with n=%d\n",(LIFT));       \
    CALL(mess_vector_printinfo((X)));                           \
    ++(ERR);                                                    \
}

#define CHECKMATRIXLIFT(A,LIFT,ALIFT,AZEROS,ATMP,DIFF,EPS,ERR)          \
    CALL(mess_matrix_resize((AZEROS),(LIFT),cols));                     \
CALL(mess_matrix_cat((A),NULL,(AZEROS),NULL,MESS_DENSE,(ALIFT)));   \
CALL(mess_matrix_lift((A),(LIFT), (ATMP)));                         \
mess_matrix_printinfo(A);                                           \
mess_matrix_print(A);                                               \
mess_matrix_printinfo(AZEROS);                                      \
mess_matrix_print(AZEROS);                                          \
mess_matrix_printinfo(ALIFT);                                       \
mess_matrix_print(ALIFT);                                           \
mess_matrix_print(ATMP);                                            \
CALL(mess_matrix_add(-1,(ALIFT),1,(ATMP)));                         \
CALL(mess_matrix_normf((ATMP),&(DIFF)));                            \
if((DIFF)>(EPS)){                                                   \
    printf("mess_matrix_lift failed with n=%d\n",(LIFT));           \
    CALL(mess_matrix_printinfo(A));                                 \
    printf("diff=%e\n",(DIFF));                                     \
    ++(ERR);                                                        \
}




int main ( int argc, char ** argv) {
    mess_init();

    double diff=.0, eps=1e-8;
    int ret=0, err=0;
    int rows=2,cols=3;
    int lift=0;
    mess_matrix A, Alift, Atmp, Azeros;
    mess_vector x, xlift,xzeros,xtmp;

    /*-----------------------------------------------------------------------------
     *  init matrices and vectors
     *-----------------------------------------------------------------------------*/

    CALL(mess_vector_init(&x,rows,MESS_REAL));
    CALL(mess_vector_init(&xlift,rows,MESS_REAL));
    CALL(mess_vector_init(&xzeros,rows,MESS_REAL));
    CALL(mess_vector_init(&xtmp,rows,MESS_REAL));

    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_init(&Alift));
    CALL(mess_matrix_init(&Atmp));
    CALL(mess_matrix_init(&Azeros));

    //CALL(mess_matrix_alloc(A,rows,cols,rows*cols,MESS_DENSE,MESS_REAL));
    CALL(mess_matrix_alloc(Alift,rows,cols,rows*cols,MESS_DENSE,MESS_REAL));
    CALL(mess_matrix_alloc(Atmp,rows,cols,rows*cols,MESS_DENSE,MESS_REAL));
    CALL(mess_matrix_alloc(Azeros,rows,cols,rows*cols,MESS_DENSE,MESS_REAL));

    /*-----------------------------------------------------------------------------
     *  check @ref mess_vector_lift
     *-----------------------------------------------------------------------------*/
    CALL(mess_vector_rand(x));
    for(lift=1;lift<10;++lift){
        CHECKVECTORLIFT(x,lift,xlift,xzeros,xtmp,diff,eps,err);
    }

    CALL(mess_vector_tocomplex(x));
    CALL(mess_vector_scalec(1+2*I,x));
    for(lift=1;lift<10;++lift){
        CHECKVECTORLIFT(x,lift,xlift,xzeros,xtmp,diff,eps,err);
    }

    /*-----------------------------------------------------------------------------
     *  check @ref mess_matrix_lift
     *-----------------------------------------------------------------------------*/

    //test DENSE-REAL
    CALL(mess_matrix_rand(A,rows,cols,MESS_DENSE,MESS_REAL,1.0));
    for(lift=1;lift<3;++lift){
        CHECKMATRIXLIFT(A,lift,Alift,Azeros,Atmp,diff,eps,err);
    }


    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/

    CALL(mess_matrix_clear(&A));
    CALL(mess_matrix_clear(&Alift));
    CALL(mess_matrix_clear(&Atmp));
    CALL(mess_matrix_clear(&Azeros));

    CALL(mess_vector_clear(&x));
    CALL(mess_vector_clear(&xlift));
    CALL(mess_vector_clear(&xzeros));
    CALL(mess_vector_clear(&xtmp));

    return (err>0)?(1):(0);
}

