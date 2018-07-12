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
 * @file tests/matrix/check_matrix_gbiorth.c
 * @brief Check the biorthonormal basis construction with respect to a weighted scalar product.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_gbiorth function defined in biorth.c that means it checks if a biorthonormal
 * basis
 * \f[ W_{out}^T E V_{out}=I \f]
 * is constructed correctly from input matrices \f$ V_{in} \f$, \f$ W_{in} \f$ and a symmetric positiv definite matrix
 * \f$ E \f$ defining a weighted scalar product.
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


#define CHECKMARTIXGBIORTH(E,VIN,WIN,VOUT,WOUT,EYE,TEMP1,TEMP2,ERR,EPS,DIFF){                   \
    CALL(mess_matrix_gbiorth(E,VIN,WIN,VOUT,WOUT))                                      \
    CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE, WOUT, MESS_OP_NONE, E, TEMP1));           \
    CALL( mess_matrix_multiply(MESS_OP_NONE, TEMP1, MESS_OP_NONE, VOUT, TEMP2));                \
    CALL(mess_matrix_diffnorm(TEMP2,EYE,&DIFF));                                            \
    if(DIFF>EPS){                                                                       \
        ++ERR;                                                                          \
        printf("mess_matrix_mgsadd failed \n");                                         \
        printf("Error:%e\n",DIFF);                                                      \
    }                                                                                   \
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
    mess_init();
    int Wrows=10,Wcols=4, ret=0, err=0;
    double p = 1.0;                     //density for random matrix generation
    double eps = 1e-12;
    double diff=0;

    mess_matrix Win,Vin,Wout,Vout,E,Temp1, Temp2, Eye;


    //init real matrices
    CALL(mess_matrix_init(&Win));
    CALL(mess_matrix_init(&Wout));
    CALL(mess_matrix_init(&Vin));
    CALL(mess_matrix_init(&Vout));
    CALL(mess_matrix_init(&E));
    CALL(mess_matrix_init(&Eye));

    //init temp matrices;
    CALL(mess_matrix_init(&Temp1));
    CALL(mess_matrix_init(&Temp2));


    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    //load  matrices
    CALL(mess_matrix_rand(Win,Wrows,Wcols,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_rand(Vin,Wrows,Wcols,MESS_DENSE,MESS_REAL,p));


    //load eye matrices
    CALL(mess_matrix_eye(Eye,Wcols,Wcols,MESS_DENSE));

    //load spd mass matrix
    generatespd(E,Wrows);

    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/

    CHECKMARTIXGBIORTH(E,Vin,Win,Vout,Wout,Eye,Temp1,Temp2,err,eps,diff)

        /*-----------------------------------------------------------------------------
         *  Clear Memory
         *-----------------------------------------------------------------------------*/

        mess_matrix_clear(&Win);
    mess_matrix_clear(&Wout);
    mess_matrix_clear(&Vin);
    mess_matrix_clear(&Vout);
    mess_matrix_clear(&E);
    mess_matrix_clear(&Eye);
    mess_matrix_clear(&Temp1);
    mess_matrix_clear(&Temp2);

    return (err>0)?(1):(0);
}
