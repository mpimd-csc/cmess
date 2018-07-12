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
 * @file tests/formats/check_perm.c
 * @brief Check the permutation function for matrices.
 *
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_perm
 * for @ref MESS_REAL, @ref MESS_COMPLEX, @ref MESS_DENSE, @ref MESS_CSR, @ref MESS_CSC format.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mess/mess.h"
//#include "mess/interface_cholmod.h"
#include "../call_macro.h"



mess_int_t * permute(mess_int_t dim) {
    mess_int_t *p = malloc(sizeof(mess_int_t)*dim);
    mess_int_t k;
    for (k = 0; k < dim; k++)
        p[k] = k;
    for (k = dim-1; k > 0; k--) {
        mess_int_t j = rand() % (k+1);
        mess_int_t temp = p[j];
        p[j] = p[k];
        p[k] = temp;
    }
    return p;
}


#define CHECKPERM(A1,A2,DIFF,ERR,EPS)                               \
    CALL(mess_matrix_diffnorm(A1,A2,&DIFF));                    \
printf("diff=%e\n",DIFF);                                   \
if(DIFF>EPS){                                               \
    printf("Failed with diff=%e\n",DIFF);                   \
    mess_matrix_printinfo(A1);                              \
    mess_matrix_printinfo(A2);                              \
}



int main(){
    mess_init();
    int ret=0, err=0;
    int dim= 50;
    double diff, eps=1e-10;


    mess_matrix A, Adense, Acsc, Acsr,A_cpx, Adense_cpx, Acsc_cpx, Acsr_cpx;

    /*-----------------------------------------------------------------------------
     *  Init/Load matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_init(&Adense));
    CALL(mess_matrix_init(&Acsc));
    CALL(mess_matrix_init(&Acsr));
    CALL(mess_matrix_init(&A_cpx));
    CALL(mess_matrix_init(&Adense_cpx));
    CALL(mess_matrix_init(&Acsc_cpx));
    CALL(mess_matrix_init(&Acsr_cpx));

    /*
       CALL(mess_matrix_alloc(A,dim,dim,dim*dim,MESS_DENSE,MESS_REAL));
       CALL(mess_matrix_alloc(Adense,dim,dim,dim*dim,MESS_DENSE,MESS_REAL));
       CALL(mess_matrix_alloc(Acsc,dim,dim,dim*dim,MESS_CSC,MESS_REAL));
       CALL(mess_matrix_alloc(Acsr,dim,dim,dim*dim,MESS_CSR,MESS_REAL));
       CALL(mess_matrix_alloc(A_cpx,dim,dim,dim*dim,MESS_DENSE,MESS_COMPLEX));
       CALL(mess_matrix_alloc(Adense_cpx,dim,dim,dim*dim,MESS_DENSE,MESS_COMPLEX));
       CALL(mess_matrix_alloc(Acsc_cpx,dim,dim,dim*dim,MESS_CSC,MESS_COMPLEX));
       CALL(mess_matrix_alloc(Acsr_cpx,dim,dim,dim*dim,MESS_CSR,MESS_COMPLEX));
       */

    CALL(mess_matrix_rand(A,dim,dim,MESS_DENSE,MESS_REAL,1));

    CALL(mess_matrix_copy(A,A_cpx));
    CALL(mess_matrix_scalec(1+I,A_cpx));

    /*-----------------------------------------------------------------------------
     *  build a permutation
     *-----------------------------------------------------------------------------*/
    mess_int_t *p,*q;
    p = permute(dim);
    q = permute(dim);

    /*-----------------------------------------------------------------------------*
     * permute and check
     *-----------------------------------------------------------------------------*/
    //test with p and q, real matrices
    CALL(mess_matrix_convert(A,Adense,MESS_DENSE));
    CALL(mess_matrix_convert(A,Acsc,MESS_CSC));
    CALL(mess_matrix_convert(A,Acsr,MESS_CSR));
    CALL(mess_matrix_perm(Adense,p,q));
    CALL(mess_matrix_perm(Acsc,p,q));
    CALL(mess_matrix_perm(Acsr,p,q));
    CHECKPERM(Adense,Acsc,diff,err,eps);
    CHECKPERM(Adense,Acsr,diff,err,eps);
    CHECKPERM(Acsc,Acsr,diff,err,eps);

    //test with p and q=identity, real matrices
    CALL(mess_matrix_convert(A,Adense,MESS_DENSE));
    CALL(mess_matrix_convert(A,Acsc,MESS_CSC));
    CALL(mess_matrix_convert(A,Acsr,MESS_CSR));
    CALL(mess_matrix_perm(Adense,p,NULL));
    CALL(mess_matrix_perm(Acsc,p,NULL));
    CALL(mess_matrix_perm(Acsr,p,NULL));
    CHECKPERM(Adense,Acsc,diff,err,eps);
    CHECKPERM(Adense,Acsr,diff,err,eps);
    CHECKPERM(Acsc,Acsr,diff,err,eps);

    //test with p=identity and q, real matrices
    CALL(mess_matrix_convert(A,Adense,MESS_DENSE));
    CALL(mess_matrix_convert(A,Acsc,MESS_CSC));
    CALL(mess_matrix_convert(A,Acsr,MESS_CSR));
    CALL(mess_matrix_perm(Adense,NULL,q));
    CALL(mess_matrix_perm(Acsc,NULL,q));
    CALL(mess_matrix_perm(Acsr,NULL,q));
    CHECKPERM(Adense,Acsc,diff,err,eps);
    CHECKPERM(Adense,Acsr,diff,err,eps);
    CHECKPERM(Acsc,Acsr,diff,err,eps);

    //test with p and q, complex matrices
    CALL(mess_matrix_convert(A_cpx,Adense_cpx,MESS_DENSE));
    CALL(mess_matrix_convert(A_cpx,Acsc_cpx,MESS_CSC));
    CALL(mess_matrix_convert(A_cpx,Acsr_cpx,MESS_CSR));
    CALL(mess_matrix_perm(Adense_cpx,p,q));
    CALL(mess_matrix_perm(Acsc_cpx,p,q));
    CALL(mess_matrix_perm(Acsr_cpx,p,q));
    CHECKPERM(Adense_cpx,Acsc_cpx,diff,err,eps);
    CHECKPERM(Adense_cpx,Acsr_cpx,diff,err,eps);
    CHECKPERM(Acsc_cpx,Acsr_cpx,diff,err,eps);

    //test with p and q=identity, complex matrices
    CALL(mess_matrix_convert(A_cpx,Adense_cpx,MESS_DENSE));
    CALL(mess_matrix_convert(A_cpx,Acsc_cpx,MESS_CSC));
    CALL(mess_matrix_convert(A_cpx,Acsr_cpx,MESS_CSR));
    CALL(mess_matrix_perm(Adense_cpx,p,NULL));
    CALL(mess_matrix_perm(Acsc_cpx,p,NULL));
    CALL(mess_matrix_perm(Acsr_cpx,p,NULL));
    CHECKPERM(Adense_cpx,Acsc_cpx,diff,err,eps);
    CHECKPERM(Adense_cpx,Acsr_cpx,diff,err,eps);
    CHECKPERM(Acsc_cpx,Acsr_cpx,diff,err,eps);

    //test with p=identity and q, complex matrices
    CALL(mess_matrix_convert(A_cpx,Adense_cpx,MESS_DENSE));
    CALL(mess_matrix_convert(A_cpx,Acsc_cpx,MESS_CSC));
    CALL(mess_matrix_convert(A_cpx,Acsr_cpx,MESS_CSR));
    CALL(mess_matrix_perm(Adense_cpx,NULL,q));
    CALL(mess_matrix_perm(Acsc_cpx,NULL,q));
    CALL(mess_matrix_perm(Acsr_cpx,NULL,q));
    CHECKPERM(Adense_cpx,Acsc_cpx,diff,err,eps);
    CHECKPERM(Adense_cpx,Acsr_cpx,diff,err,eps);
    CHECKPERM(Acsc_cpx,Acsr_cpx,diff,err,eps);


    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/

    CALL(mess_matrix_clear(&A));
    CALL(mess_matrix_clear(&Acsc));
    CALL(mess_matrix_clear(&Acsr));
    CALL(mess_matrix_clear(&Adense));
    CALL(mess_matrix_clear(&A_cpx));
    CALL(mess_matrix_clear(&Acsc_cpx));
    CALL(mess_matrix_clear(&Acsr_cpx));
    CALL(mess_matrix_clear(&Adense_cpx));
    free(p);
    free(q);

    return err ;
}


