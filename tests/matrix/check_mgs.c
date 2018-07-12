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
 * @file tests/matrix/check_mgs.c
 * @brief Check the ortohogonal basis construction using the modified Gram-Schmidt process.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_mgs function defined in mgs.c that means it checks if an orthogonal basis for
 * a given matrix is computed using the modified Gram-Schmidt process.
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

#define CHECKMGS(A,Q,R,EYE,TEMP1,TEMP2,ERR,EPS,DIFF){                       \
    mess_matrix_mgs (A, Q, R);                                              \
    if(R!=NULL){                                                            \
        mess_matrix_multiply(MESS_OP_NONE, Q, MESS_OP_NONE, R, TEMP1);              \
        mess_matrix_diffnorm(TEMP1,A,&DIFF);                                        \
        if(DIFF>EPS){++ERR;}                                                    \
    }                                                                       \
    mess_matrix_copy(Q,TEMP1);                                              \
    mess_matrix_multiply(MESS_OP_HERMITIAN, Q, MESS_OP_NONE, TEMP1, TEMP2); \
    mess_matrix_diffnorm(TEMP2,EYE,&DIFF);                                      \
    if(DIFF>eps){++ERR;}                                                    \
}

int main ( int argc, char ** argv) {
    mess_init();
    int n=4, ret=0, err=0;
    double p = 0.6;                     //density for random matrix generation
    double eps = 1e-12;
    double diff=0;

    mess_matrix Ar,Ac,Qr,Qc,Rr,Rc, temp1,temp2, eyer,eyec;

    //init A matrices
    CALL(mess_matrix_init(&Ar));
    CALL(mess_matrix_init(&Ac));

    //init Q matrices
    CALL(mess_matrix_init(&Qr));
    CALL(mess_matrix_init(&Qc));

    //init R matrices
    CALL(mess_matrix_init(&Rr));
    CALL(mess_matrix_init(&Rc));

    CALL(mess_matrix_init(&temp1));
    CALL(mess_matrix_init(&temp2));
    CALL(mess_matrix_init(&eyer));
    CALL(mess_matrix_init(&eyec));


    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    //load A matrices
    CALL(mess_matrix_rand(Ar,n,n,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_rand(Ac,n,n,MESS_DENSE,MESS_COMPLEX,p));

    //load eye matrices
    CALL(mess_matrix_eye(eyer,n,n,MESS_DENSE));
    CALL(mess_matrix_eyec(eyec,n,n,MESS_DENSE));

    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/


    CHECKMGS(Ar,Qr,Rr,eyer,temp1,temp2,err,eps,diff)
        CHECKMGS(Ar,Qr,NULL,eyer,temp1,temp2,err,eps,diff)
        CHECKMGS(Ac,Qc,Rc,eyec,temp1,temp2,err,eps,diff)
        CHECKMGS(Ac,Qc,NULL,eyec,temp1,temp2,err,eps,diff)

        /*-----------------------------------------------------------------------------
         *  Clear Memory
         *-----------------------------------------------------------------------------*/


        mess_matrix_clear(&Ar);
    mess_matrix_clear(&Ac);
    mess_matrix_clear(&Qr);
    mess_matrix_clear(&Qc);
    mess_matrix_clear(&Rr);
    mess_matrix_clear(&Rc);
    mess_matrix_clear(&temp1);
    mess_matrix_clear(&temp2);
    mess_matrix_clear(&eyer);
    mess_matrix_clear(&eyec);

    return (err>0)?(1):(0);
}
