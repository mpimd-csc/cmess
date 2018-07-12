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
 * @file tests/matrix/check_matrix_biorth.c
 * @brief Check the biorthonormal basis construction.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_biorth function defined in biorth.c that means it checks if a biorthonormal
 * basis
 * \f[ W_{out}^T V_{out}=I \f]
 * is constructed correctly from input matrices \f$ V_{in} \f$ and \f$ W_{in} \f$.
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


#define CHECKMARTIXBIORTH(VIN,WIN,VOUT,WOUT,EYE,TEMP1,ERR,EPS,DIFF){                        \
    CALL(mess_matrix_biorth(VIN,WIN,VOUT,WOUT))                                         \
    CALL( mess_matrix_multiply(MESS_OP_TRANSPOSE, WOUT, MESS_OP_NONE, VOUT, TEMP1));            \
    CALL(mess_matrix_diffnorm(TEMP1,EYE,&DIFF));                                            \
    if(DIFF>EPS){                                                                       \
        ++ERR;                                                                          \
        printf("mess_matrix_mgsadd failed \n");                                         \
        printf("Error:%e\n",DIFF);                                                      \
    }                                                                                   \
}




int main ( int argc, char ** argv) {
    mess_init();
    int Wrows=10,Wcols=4, ret=0, err=0;
    double p = 1.0;                     //density for random matrix generation
    double eps = 1e-12;
    double diff=0;

    mess_matrix Win,Vin,Wout,Vout,Temp1, Eye;


    //init real matrices
    CALL(mess_matrix_init(&Win));
    CALL(mess_matrix_init(&Wout));
    CALL(mess_matrix_init(&Vin));
    CALL(mess_matrix_init(&Vout));
    CALL(mess_matrix_init(&Eye));

    //init temp matrices;
    CALL(mess_matrix_init(&Temp1));


    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    //load  matrices
    CALL(mess_matrix_rand(Win,Wrows,Wcols,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_rand(Vin,Wrows,Wcols,MESS_DENSE,MESS_REAL,p));


    //load eye matrices
    CALL(mess_matrix_eye(Eye,Wcols,Wcols,MESS_DENSE));


    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/

    CHECKMARTIXBIORTH(Vin,Win,Vout,Wout,Eye,Temp1,err,eps,diff)

        /*-----------------------------------------------------------------------------
         *  Clear Memory
         *-----------------------------------------------------------------------------*/

        mess_matrix_clear(&Win);
    mess_matrix_clear(&Wout);
    mess_matrix_clear(&Vin);
    mess_matrix_clear(&Vout);
    mess_matrix_clear(&Eye);
    mess_matrix_clear(&Temp1);

    return (err>0)?(1):(0);
}
