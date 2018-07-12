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
 * @file tests/vector/check_kron.c
 * @brief Check the Kronecker product computation.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_vector_kron function defined in vector.c that means it checks if the Kronecker product
 * of two random vectors \f$ v_1 \f$ and \f$ v_2 \f$
 * \f[ v_1 \otimes v_2 \f]
 * is computed correctly.
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


//very naive test function for kronecker produkt
int testkron(mess_vector v1, mess_vector v2, mess_vector out){
    int ret;
    mess_vector temp1,temp2,temp3,temp4,temp5,temp6;

    CALL(mess_vector_init(&temp1,v1->dim,MESS_COMPLEX));
    CALL(mess_vector_init(&temp2,v2->dim,MESS_COMPLEX));
    CALL(mess_vector_init(&temp3,v2->dim,MESS_COMPLEX));
    CALL(mess_vector_init(&temp4,v2->dim,MESS_COMPLEX));
    CALL(mess_vector_init(&temp5,v2->dim,MESS_COMPLEX));
    CALL(mess_vector_init(&temp6,v2->dim,MESS_COMPLEX));

    CALL(mess_vector_copy_tocomplex(v1,temp1));
    CALL(mess_vector_copy_tocomplex(v2,temp2));

    CALL(mess_vector_copy(temp1,temp3));
    CALL(mess_vector_copy(temp2,temp4));

    mess_int_t i=0;

    CALL(mess_vector_scalec(temp1->values_cpx[i],temp4));
    for(i=1;i<v1->dim;++i){
        CALL(mess_vector_copy(temp2,temp5));
        CALL(mess_vector_scalec(temp1->values_cpx[i],temp5));
        CALL(mess_vector_cat(temp4,temp5,temp6));
        CALL(mess_vector_copy(temp6,temp4));
    }

    CALL(mess_vector_copy(temp4,out));

    CALL(mess_vector_clear(&temp1));
    CALL(mess_vector_clear(&temp2));
    CALL(mess_vector_clear(&temp3));
    CALL(mess_vector_clear(&temp4));
    CALL(mess_vector_clear(&temp5));
    CALL(mess_vector_clear(&temp6));

    return 0;

}

#define CHECK_KRON(V1,V2,V3,VTEST,DIFF,ERR,EPS)                         \
    CALL(mess_vector_kron(V1,V2,V3));                               \
CALL(testkron(V1,V2,VTEST));                                    \
CALL(mess_vector_diffnorm(V3,VTEST,&DIFF));                     \
if(DIFF>EPS){                                                   \
    printf("failed\tdiff=%e\n",diff);                           \
    mess_vector_printinfo(V1);                                  \
    mess_vector_printinfo(V2);                                  \
    mess_vector_printinfo(V3);                                  \
    ERR++;                                                      \
}else{                                                          \
    printf(" passed with diff=%e.\n",diff);                     \
}                                                               \



int main ( int argc, char ** argv) {
    mess_init();
    mess_int_t dim1=2,dim2=3;
    int ret;
    int err = 0;
    double eps = 1e-10, diff;
    mess_vector  v1,v2,v3,vtest;

    /*-----------------------------------------------------------------------------
     *  init vectors
     *-----------------------------------------------------------------------------*/

    //init vector
    CALL(mess_vector_init(&v1, dim1, MESS_REAL));
    CALL(mess_vector_init(&v2, dim2, MESS_REAL));
    CALL(mess_vector_init(&v3, dim1*dim2, MESS_REAL));
    CALL(mess_vector_init(&vtest, dim1*dim2, MESS_REAL));
    CALL(mess_vector_rand(v1));
    CALL(mess_vector_rand(v2));

    /*-----------------------------------------------------------------------------
     *  test different cases
     *-----------------------------------------------------------------------------*/
    for(dim2=2;dim2<10;++dim2){
        for(dim1=2;dim1<5;++dim1){
            CALL(mess_vector_resize(v1,dim1));
            CALL(mess_vector_rand(v1));
            CALL(mess_vector_resize(v2,dim2));
            CALL(mess_vector_rand(v2));

            CHECK_KRON(v1,v2,v3,vtest,diff,err,eps);
            CHECK_KRON(v2,v1,v3,vtest,diff,err,eps);
            mess_vector_tocomplex(v1);
            mess_vector_toreal_nowarn(v2);

            CHECK_KRON(v1,v2,v3,vtest,diff,err,eps);
            CHECK_KRON(v2,v1,v3,vtest,diff,err,eps);
            mess_vector_toreal_nowarn(v1);
            mess_vector_tocomplex(v2);

            CHECK_KRON(v1,v2,v3,vtest,diff,err,eps);
            CHECK_KRON(v2,v1,v3,vtest,diff,err,eps);
            mess_vector_tocomplex(v1);
            mess_vector_tocomplex(v2);
        }

    }





    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/
    CALL(mess_vector_clear(&v1));
    CALL(mess_vector_clear(&v2));
    CALL(mess_vector_clear(&v3));
    CALL(mess_vector_clear(&vtest));










    return (err>0)?(1):(0);
}
