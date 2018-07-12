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
 * @file tests/matrix/check_colscale.c
 * @brief Check the column scaling function of a matrix.
 * @author  @mbehr
 * @test
 * This function check the @ref mess_matrix_colscale function defined in colops.c, that means it checks if the scaling
 * \f[Q(:,col) \leftarrow scale \cdot Q(:,col)  \f]
 * is computed correctly for a given matrix \f$ Q \f$, a column \f$ col \f$ and a scaling factor \f$ scale \f$.
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

double real_matrix_data[9] = {1,2,3,4,5,6,7,8,9};
mess_double_cpx_t complex_matrix_data[9] = {1,2*I,3,4*I,5,6*I,7,8*I,9};

#define CHECKCOLSCALE(A,COL,SCALAR, DIFF, EPS,ERR){                             \
    mess_vector check1, check2;                                             \
    MESS_INIT_VECTORS(&check1,&check2);                                     \
    mess_vector_alloc(check1, (A)->cols,(A)->data_type);                    \
    mess_vector_alloc(check2, (A)->cols,(A)->data_type);                    \
    mess_matrix_getcol((A),(COL),check1);                                   \
    if(MESS_IS_REAL(check1)){                                               \
        mess_vector_scale((SCALAR),check1);                                     \
    }else{                                                                  \
        mess_vector_scalec((SCALAR),check1);                                        \
    }                                                                       \
    mess_matrix_colscale((A),(COL),(SCALAR));                               \
    mess_matrix_getcol((A),(COL),check2);                                   \
    mess_vector_diffnorm(check1,check2,&(DIFF));                            \
    if((DIFF)>10*(EPS)){++(ERR); printf("diff = %lg\n", DIFF);};                                            \
    mess_vector_clear(&check1);                                             \
    mess_vector_clear(&check2);                                             \
}



int main ( int argc, char ** argv) {
    mess_init();

    double diff=.0, eps=mess_eps();
    int ret=0, err=0;
    mess_double_cpx_t scalac=3-2*I;
    double scalar = 3;

    //prepare dense matrices
    mess_matrix denser, densec;
    CALL(mess_matrix_init(&denser));
    CALL(mess_matrix_init(&densec));
    CALL(mess_matrix_dense_from_farray(denser,3,3,0,real_matrix_data,NULL));
    CALL(mess_matrix_dense_from_farray(densec,3,3,0,NULL, complex_matrix_data));


    //prepare csr matrices
    mess_matrix csrr,csrc;
    CALL(mess_matrix_init(&csrr));
    CALL(mess_matrix_init(&csrc));
    CALL(mess_matrix_convert(denser,csrr,MESS_CSR));
    CALL(mess_matrix_convert(densec,csrc,MESS_CSR));

    //prepare csr matrices
    mess_matrix cscr,cscc;
    CALL(mess_matrix_init(&cscr));
    CALL(mess_matrix_init(&cscc));
    CALL(mess_matrix_convert(denser,cscr,MESS_CSC));
    CALL(mess_matrix_convert(densec,cscc,MESS_CSC));

    mess_int_t i=0;
    for(i=0;i<3;++i){
        //check dense
        CHECKCOLSCALE(densec,i,scalac,diff,eps,err);
        CHECKCOLSCALE(denser,i,scalar,diff,eps,err);

        //check csc
        CHECKCOLSCALE(cscc,i,scalac,diff,eps,err);
        CHECKCOLSCALE(cscr,i,scalar,diff,eps,err);

        //check csr
        CHECKCOLSCALE(csrc,i,scalac,diff,eps,err);
        CHECKCOLSCALE(csrr,i,scalar,diff,eps,err);
    }
    mess_matrix_clear(&denser);
    mess_matrix_clear(&densec);
    mess_matrix_clear(&cscr);
    mess_matrix_clear(&cscc);
    mess_matrix_clear(&csrr);
    mess_matrix_clear(&csrc);

    return (err>0)?(1):(0);
}

