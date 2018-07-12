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
 * @addtogroup test_lrcfadi
 * @{
 *
 * @file tests/lrcf_adi/check_dae2_handles_para.c
 * @brief Check function handles for Hessenberg index 2 DAE systems, which are necessary for mess_lrcfadi_parameter.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_equation_glyap_dae2 function defined in equation_glyap_dae2.c.
 *
 * \sa mess_equation_glyap_dae2
 *
 * @}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mess/mess.h"
#include "../lib/lrcf_adi/equation_glyap_dae2.c"


#include "../call_macro.h"

#define CHECK_FUNCTION(INFO,FNAME,DIFFFUNC,Y1,Y2,DIFF,ERR,EPS)          \
    CALL(DIFFFUNC(Y1,Y2,&DIFF));                                    \
if(DIFF>EPS){                                                   \
    printf(INFO " " #FNAME " failed\tdiff=%e\n",diff);          \
    ERR++;                                                      \
}else{                                                          \
    printf(INFO " " #FNAME  " passed with diff=%e.\n",diff);    \
}                                                               \

int main ( int argc, char ** argv) {

    int ret=0, err=0;
    double diff=.0, eps=1e-6, delta=-0.02;
    mess_int_t nv, np, i, j, cols=2;

    mess_matrix M, A, G, B, Gt, fullA, fullM, Tmp1, Tmp2, X1, X2, Y1, Y2;
    mess_vector x1, x2, y1, y2;
    mess_equation eqn;
    mess_options opt;

    mess_version();
    mess_init();


    if ( argc != 6 ) {
        fprintf(stderr, "usage: %s M.mtx A.mtx G.mtx B.mtx memory_usage \n", argv[0]);
        return -1;
    }


    /*-----------------------------------------------------------------------------
     *  Init and read data
     *-----------------------------------------------------------------------------*/
    //init matrices
    CALL(mess_matrix_init(&M));
    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_init(&G));
    CALL(mess_matrix_init(&B));
    CALL(mess_matrix_init(&Gt));
    CALL(mess_matrix_init(&fullA));
    CALL(mess_matrix_init(&fullM));
    CALL(mess_matrix_init(&Tmp1));
    CALL(mess_matrix_init(&Tmp2));
    CALL(mess_matrix_init(&X1));
    CALL(mess_matrix_init(&X2));
    CALL(mess_matrix_init(&Y1));
    CALL(mess_matrix_init(&Y2));

    //read matrices
    CALL(mess_matrix_read_formated(argv[1], M, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[2], A, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[3], G, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[4], B, MESS_DENSE));

    //size of linearized navier stokes system for velocity and pressure
    nv = M->rows;
    np = G->cols;

    //random data for matrix
    CALL(mess_matrix_rand(X1,nv+np,cols,MESS_DENSE,MESS_REAL,1.0));
    CALL(mess_matrix_copy(X1,X2));

    //build fullM =  [M,delta*G;delta*G^T,0]- = [M,Tmp1;Tmp2,0]
    CALL(mess_matrix_copy(G,Tmp1));
    CALL(mess_matrix_scale(delta,Tmp1));
    CALL(mess_matrix_ctranspose(Tmp1,Tmp2));
    CALL(mess_matrix_cat(M,Tmp1,Tmp2,NULL,MESS_CSR,fullM));

    //build fullA = [A,G;G^T,0];
    CALL(mess_matrix_ctranspose(G,Tmp1));
    CALL(mess_matrix_cat(A,G,Tmp1,NULL,MESS_CSR,fullA));

    //init and alloc vector
    MESS_INIT_VECTORS(&x1,&x2,&y1,&y2);
    CALL(mess_vector_alloc(x1,nv+np,MESS_REAL));
    CALL(mess_vector_alloc(x2,nv+np,MESS_REAL));
    CALL(mess_vector_alloc(y1,nv+np,MESS_REAL));
    CALL(mess_vector_alloc(y2,nv+np,MESS_REAL));

    //random data for vector
    CALL(mess_vector_rand(x1));
    CALL(mess_vector_copy(x1,x2));

    //init equation structure
    CALL(mess_equation_init(&eqn));

    //init options structure
    CALL(mess_options_init(&opt));
    opt->memory_usage = atoi(argv[5]);

    //build equation structure for Index 2 DAE
    CALL(mess_equation_glyap_dae2(eqn,opt,M,A,G,B,delta));

    //set function handles for shifts parameter computation
    eqn->AX.apply           = AX_apply_shifts;
    eqn->EX.apply           = EX_apply_shifts;
    eqn->AINV.apply         = AINV_apply_shifts;
    eqn->EINV.apply         = EINV_apply_shifts;

    /*-----------------------------------------------------------------------------
     *  Test function handles for shift parameter computation, i.e. size == nv+np
     *-----------------------------------------------------------------------------*/
    int ops [] = {MESS_OP_NONE,MESS_OP_TRANSPOSE,MESS_OP_HERMITIAN};

    for(j=0;j<2;++j){
        for(i=0;i<3;++i){
            //test AX.apply
            eqn->AX.apply(eqn,ops[i],X1,Y1);
            CALL( mess_matrix_multiply(ops[i], fullA, MESS_OP_NONE, X2, Y2));
            CHECK_FUNCTION("Shift",AX.apply,mess_matrix_diffnorm,Y1,Y2,diff,err,eps);

            //test EX.apply
            eqn->EX.apply(eqn,ops[i],X1,Y1);
            CALL( mess_matrix_multiply(ops[i], fullM, MESS_OP_NONE, X2, Y2));
            CHECK_FUNCTION("Shift",EX.apply,mess_matrix_diffnorm,Y1,Y2,diff,err,eps);

            //test AINV.apply
            //eqn->AINV.generate(eqn);
            eqn->AINV.apply(eqn,ops[i],X1,Y1);
            CALL( mess_matrix_multiply(ops[i], fullA, MESS_OP_NONE, Y1, X2));
            CHECK_FUNCTION("Shift",AINV.apply,mess_matrix_diffnorm,X1,X2,diff,err,eps);
            //eqn->AINV.clear(eqn);

            //test EINV.apply
            //CALL(eqn->EINV.generate(eqn));
            eqn->EINV.apply(eqn,ops[i],X1,Y1);
            mess_matrix_multiply(ops[i], fullM, MESS_OP_NONE, Y1, X2);
            CHECK_FUNCTION("Shift",EINV.apply,mess_matrix_diffnorm,X1,X2,diff,err,eps);
            //eqn->EINV.clear(eqn);
        }

        //perform same tests with complex data
        CALL(mess_vector_tocomplex(x1));
        CALL(mess_vector_scalec(1+2*I,x1));
        CALL(mess_vector_copy(x1,x2));

        CALL(mess_matrix_tocomplex(X1));
        CALL(mess_matrix_scalec(1+2*I,X1));
        CALL(mess_matrix_copy(X1,X2));
    }


    /*-----------------------------------------------------------------------------
     *  Clear data
     *-----------------------------------------------------------------------------*/
    //clear matrices
    CALL(mess_matrix_clear(&M));
    CALL(mess_matrix_clear(&A));
    CALL(mess_matrix_clear(&G));
    CALL(mess_matrix_clear(&B));
    CALL(mess_matrix_clear(&Gt));
    CALL(mess_matrix_clear(&fullA));
    CALL(mess_matrix_clear(&fullM));
    CALL(mess_matrix_clear(&Tmp1));
    CALL(mess_matrix_clear(&Tmp2));
    CALL(mess_matrix_clear(&X1));
    CALL(mess_matrix_clear(&X2));
    CALL(mess_matrix_clear(&Y1));
    CALL(mess_matrix_clear(&Y2));

    //clear vectors
    CALL(mess_vector_clear(&x1));
    CALL(mess_vector_clear(&x2));
    CALL(mess_vector_clear(&y1));
    CALL(mess_vector_clear(&y2));

    //clear eqn structure
    CALL(mess_equation_clear(&eqn));

    //clear options structure
    CALL(mess_options_clear(&opt));

    return (err>0)?(1):(0);
}















