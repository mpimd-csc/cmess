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
 * @file tests/matrix/check_colaxpy2.c
 * @brief Check the \ref mess_matrix_colaxpy2 function.
 * @author  @mbehr
 * @test
 *  This function checks the @ref mess_matrix_colaxpy2 function defined in colops.c, that means it checks if
 *  \f[ Q(:,col_1) \leftarrow  Q(:,col_1)+coeff*Q(:,col_c)\f]
 *  is computed correctly for a given matrix \f$ Q \f$ and given columns \f$ col_1, col_c \f$
 *  and a coefficient \f$ coeff \f$.
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



int main ( int argc, char ** argv) {

    mess_init();
    double diff=.0, eps=mess_eps();
    int ret=0, err=0;
    mess_int_t n=104;


    //prepare dense matrices
    mess_matrix denser, densec;
    CALL(mess_matrix_init(&denser));
    CALL(mess_matrix_init(&densec));
    CALL(mess_matrix_alloc(denser,n,n,n*n,MESS_DENSE,MESS_REAL));
    CALL(mess_matrix_alloc(densec,n,n,n*n,MESS_DENSE,MESS_COMPLEX));

    //prepare csr matrices
    mess_matrix csrr,csrc;
    CALL(mess_matrix_init(&csrr));
    CALL(mess_matrix_init(&csrc));
    CALL(mess_matrix_alloc(csrr,n,n,n*n,MESS_CSR,MESS_REAL));
    CALL(mess_matrix_alloc(csrr,n,n,n*n,MESS_CSR,MESS_COMPLEX));

    //prepare csc matrices
    mess_matrix cscr,cscc;
    CALL(mess_matrix_init(&cscr));
    CALL(mess_matrix_init(&cscc));
    CALL(mess_matrix_alloc(cscr,n,n,n*n,MESS_CSC,MESS_REAL));
    CALL(mess_matrix_alloc(cscr,n,n,n*n,MESS_CSC,MESS_COMPLEX));

    //test values
    double p [5]={0.1,0.2,0.5,0.7,1.0}; int k=0;
    mess_int_t col1[5]={0,5,0,4,2};     int l=0;
    mess_int_t colc[5]={1,6,2,1,8};     int m=0;

    for(k=0;k<5;++k){
        for(l=0;l<5;++l){
            for(m=0;m<5;++m){

                //test real csc case
                CALL(mess_matrix_rand_csc(cscr,n,n,p[k],MESS_REAL));
                CALL(mess_matrix_convert(cscr,denser,MESS_DENSE));
                CALL(mess_matrix_colaxpy2(denser,p[k]+I,colc[l],col1[m]));
                CALL(mess_matrix_colaxpy2(cscr,p[k]+I,colc[l],col1[m]));
                CALL(mess_matrix_add(-1.0,cscr,1.0,denser));
                CALL(mess_matrix_normf(denser,&diff));
                if(diff>10*eps){++err;};
                printf("%e\n",diff);
                if ( err ) { printf("A\n"); abort(); }

                //test complex csc case
                CALL(mess_matrix_rand_csc(cscc,n,n,p[k],MESS_REAL));
                CALL(mess_matrix_tocomplex(cscc));
                CALL(mess_matrix_scalec(1+I,cscc));
                CALL(mess_matrix_convert(cscc,densec,MESS_DENSE));
                CALL(mess_matrix_colaxpy2(densec,p[k]+I,colc[l],col1[m]));
                CALL(mess_matrix_colaxpy2(cscc,p[k]+I,colc[l],col1[m]));
                CALL(mess_matrix_add(-1.0,cscc,1.0,densec));
                CALL(mess_matrix_normf(densec,&diff));
                if(diff>10*eps){++err;};
                printf("%e\n",diff);
                if ( err ) { printf("B\n"); abort(); }


                //test real csr case
                CALL(mess_matrix_rand_csr(csrr,n,n,p[k],MESS_REAL));
                CALL(mess_matrix_convert(csrr,denser,MESS_DENSE));
                CALL(mess_matrix_colaxpy2(denser,p[k],colc[l],col1[m]));
                CALL(mess_matrix_colaxpy2(csrr,p[k],colc[l],col1[m]));
                CALL(mess_matrix_add(-1.0,csrr,1.0,denser));
                CALL(mess_matrix_normf(denser,&diff));
                if(diff>10*eps){++err;};
                printf("%e\n",diff);
                if ( err ) { printf("C\n"); abort(); }

                //test complex csr case
                CALL(mess_matrix_rand_csr(csrc,n,n,p[k],MESS_REAL));
                CALL(mess_matrix_tocomplex(csrc));
                CALL(mess_matrix_scalec(1+I,csrc));
                CALL(mess_matrix_convert(csrc,densec,MESS_DENSE));
                CALL(mess_matrix_colaxpy2(densec,p[k]+I,colc[l],col1[m]));
                CALL(mess_matrix_colaxpy2(csrc,p[k]+I,colc[l],col1[m]));
                CALL(mess_matrix_add(-1.0,csrc,1.0,densec));
                CALL(mess_matrix_normf(densec,&diff));
                if(diff>10*eps){++err;};
                printf("%e\n",diff);
                if ( err ) { printf("D\n"); abort(); }

            }

        }
    }

    mess_matrix_clear(&cscc);
    mess_matrix_clear(&cscr);
    mess_matrix_clear(&csrc);
    mess_matrix_clear(&csrr);
    mess_matrix_clear(&denser);
    mess_matrix_clear(&densec);

    return (err>0)?(1):(0);
}

