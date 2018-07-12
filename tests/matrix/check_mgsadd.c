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
 * @file tests/matrix/check_mgsadd.c
 * @brief Check the orthogonal basis construction using the modified Gram-Schmidt process (column adding version).
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_mgs_add function defined in mgsadd.c that means it checks if an orthogonal basis
 * for
 * \f[ \left[ \begin{array}{cc} Z & V \end{array}\right] \f]
 * is computed using the modified Gram-Schmidt process under the assumption that the matrix \f$ Z \f$ is already
 * orthogonal and \f$ V \f$ are new columns to add.
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

//checkout without Mass Matrix
#define CHECKMGSADD(Z,V,EYE,TEMP1,TEMP2,TEMP3,ERR,EPS,DIFF){                                \
        CALL(mess_matrix_mgs(Z,TEMP1,NULL));                                                \
        CALL(mess_matrix_mgs_add (TEMP1, V, NULL ));                                        \
        CALL(mess_matrix_copy(TEMP1,TEMP2));                                                \
        if(MESS_IS_COMPLEX(Z)||MESS_IS_COMPLEX(V)){                                         \
        CALL( mess_matrix_multiply(MESS_OP_NONE, TEMP2, MESS_OP_HERMITIAN, TEMP1, TEMP3));      \
        }else{                                                                              \
            CALL( mess_matrix_multiply(MESS_OP_NONE, TEMP2, MESS_OP_TRANSPOSE, TEMP1, TEMP3));  \
        }                                                                                   \
        CALL(mess_matrix_diffnorm(TEMP3,EYE,&DIFF));                                            \
        if(DIFF>EPS){                                                                       \
            ++ERR;                                                                          \
            printf("mess_matrix_mgsadd failed \n");                                         \
            printf("Error:%e\n",DIFF);                                                      \
        }                                                                                   \
}

#define CHECKMGSADDE(Z,V,EYE,E,CHOLSOLVER,TEMP1,TEMP2,TEMP3,ERR,EPS,DIFF){                  \
        CALL(mess_direct_init(&CHOLSOLVER));                                                \
        CALL(mess_direct_create_cholesky(E,CHOLSOLVER));                                            \
        CALL(mess_direct_getL(CHOLSOLVER,TEMP1));                                           \
        CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE, TEMP1, MESS_OP_NONE, Z, TEMP2));          \
        CALL(mess_matrix_mgs(TEMP2,TEMP3,NULL));                                            \
        CALL(mess_matrix_backslashm(MESS_OP_TRANSPOSE,TEMP1,TEMP3,TEMP2));                  \
        CALL(mess_matrix_mgs_add(TEMP2,V,E));                                               \
        CALL(mess_matrix_copy(TEMP2,TEMP1));                                                \
        if(MESS_IS_COMPLEX(Z)||MESS_IS_COMPLEX(V)){                                         \
        CALL( mess_matrix_multiply(MESS_OP_HERMITIAN, TEMP1, MESS_OP_NONE, E, TEMP3));          \
        }else{                                                                              \
        CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE, TEMP1, MESS_OP_NONE, E, TEMP3));          \
        }                                                                                   \
        CALL( mess_matrix_multiply(MESS_OP_NONE, TEMP3, MESS_OP_NONE, TEMP2, TEMP1));           \
        CALL(mess_matrix_diffnorm(TEMP1,EYE,&DIFF));                                            \
        if(DIFF>EPS){                                                                       \
            ++ERR;                                                                          \
            printf("mess_matrix_mgsadd failed \n");                                         \
            printf("Error:%e\n",DIFF);                                                      \
        }                                                                                   \
        CALL(mess_direct_clear(&CHOLSOLVER));                                               \
}






//generate diagonal matrix with positiv elements (spd)
void generatespd(mess_matrix E, int dim){
    mess_matrix_eye(E,dim,dim, MESS_DENSE);
    int i;
    srand(time(NULL));
    for ( i=0 ; i < dim; ++i) {
        E->values[i+i*E->ld]=1+ rand()/((double) RAND_MAX);
    }
}




int main ( int argc, char ** argv) {
    int zrows=10,zcols=4,vrows=10, vcols=6,ret=0, err=0;
    double p = 1.0;                     //density for random matrix generation
    double eps = 1e-12;
    double diff=0;

    mess_matrix Zr, Vr, Eyer, Qr, Zc,Vc, E, Temp1, Temp2, Temp3;
    mess_direct chol;


    //init real matrices
    CALL(mess_matrix_init(&Zr));
    CALL(mess_matrix_init(&Vr));
    CALL(mess_matrix_init(&Qr));
    CALL(mess_matrix_init(&Zc));
    CALL(mess_matrix_init(&Vc));
    CALL(mess_matrix_init(&Eyer));
    CALL(mess_matrix_init(&E));

    //init temp matrices;
    CALL(mess_matrix_init(&Temp1));
    CALL(mess_matrix_init(&Temp2));
    CALL(mess_matrix_init(&Temp3));

    //init chol solver
    //CALL(mess_direct_init(&chol));




    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    //load  matrices
    CALL(mess_matrix_rand(Zr,zrows,zcols,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_rand(Vr,vrows,vcols,MESS_DENSE,MESS_REAL,p));


    CALL(mess_matrix_rand(Zc,vrows,zcols,MESS_DENSE,MESS_COMPLEX,p));
    CALL(mess_matrix_rand(Vc,vrows,vcols,MESS_DENSE,MESS_COMPLEX,p));


    //load eye matrices
    CALL(mess_matrix_eye(Eyer,zrows,zcols+vcols,MESS_DENSE));

    //load spd mass matrix
    generatespd(E,zrows);

    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/
    //check without mass matrix
    CHECKMGSADD(Zr,Vr,Eyer,Temp1,Temp2,Temp3,err,eps,diff);
    CHECKMGSADD(Zc,Vc,Eyer,Temp1,Temp2,Temp3,err,eps,diff);



    //check with mass matrix
    CHECKMGSADDE(Zr,Vr,Eyer,E,chol,Temp1,Temp2,Temp3,err,eps,diff);
    CHECKMGSADDE(Zc,Vc,Eyer,E,chol,Temp1,Temp2,Temp3,err,eps,diff);


    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/

    mess_matrix_clear(&Zr);
    mess_matrix_clear(&Vr);
    mess_matrix_clear(&Qr);
    mess_matrix_clear(&Zc);
    mess_matrix_clear(&Vc);
    mess_matrix_clear(&Eyer);
    mess_matrix_clear(&E);
    mess_matrix_clear(&Temp1);
    mess_matrix_clear(&Temp2);
    mess_matrix_clear(&Temp3);
    mess_direct_clear(&chol);

    return (err>0)?(1):(0);
}
