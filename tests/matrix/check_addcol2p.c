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
 * @file tests/matrix/check_addcol2p.c
 * @brief Check the \ref mess_matrix_addcols2p function.
 * @author  @mbehr
 * @test
 * This function checks if the @ref mess_matrix_addcols2p function works, that means it checks if
 * \f[ Z= [Z,s_1 \cdot V_1= s_2 \cdot V_2] \f]
 * is concatenated correcly, where \f$ s_1\f$ and \f$ s_2 \f$ are given scaling factors for given matrices \f$ V_1\f$ and \f$ V_2 \f$.
 * @}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mess/mess.h"

#include "../call_macro.h"

#define CHECKADDCOL2P(Z,S1,V1,S2,V2,TEMP1,TEMP2,TEMP3,ERR,EPS,DIFF){        \
    CALL(mess_matrix_copy(Z,TEMP1));                                    \
    CALL(mess_matrix_copy(V1,TEMP2));                                   \
    CALL(mess_matrix_add(S2,V2,S1,TEMP2));                              \
    CALL(mess_matrix_cat(TEMP1,TEMP2,NULL,NULL,MESS_DENSE,TEMP3));      \
    CALL(mess_matrix_addcols2p(Z,S1,V1,S2,V2));                         \
    CALL(mess_matrix_diffnorm(Z,TEMP3,&DIFF));                              \
    if(DIFF>EPS){                                                       \
        printf("Failed\n");                                             \
        printf("Diff:%e\n",diff);                                       \
        ++ERR;                                                          \
    }                                                                   \
}

int main ( int argc, char ** argv) {
    mess_init();
    int Zrows=10, Zcols=5, Vcols=4, err=0;
    double eps = 1e-12;
    double diff=0;
    double p  = 1;  //neglected for dense matrices
    double s1=rand(), s2=rand();
    int ret=0;

    mess_matrix Z,V1,V2,Temp1,Temp2, Temp3;

    /*-----------------------------------------------------------------------------
     *  Init matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&Z,&V1,&V2,&Temp1,&Temp2,&Temp3)

        /*-----------------------------------------------------------------------------
         *  Load matrices
         *-----------------------------------------------------------------------------*/
        CALL(mess_matrix_rand(Z,Zrows,Zcols,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_rand(V1,Zrows,Vcols,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_rand(V2,Zrows,Vcols,MESS_DENSE,MESS_REAL,p));
    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/
    CHECKADDCOL2P(Z,s1,V1,s2,V2,Temp1,Temp2,Temp3,err,eps,diff)

        /*-----------------------------------------------------------------------------
         *  Clear Memory
         *-----------------------------------------------------------------------------*/

        MESS_CLEAR_MATRICES(&Z,&V1,&V2,&Temp1,&Temp2,&Temp3)

        return (err>0)?(1):(0);
}
