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
 * @file tests/lrcf_adi/check_dae2_handles_adi.c
 * @brief Check function handles for Hessenberg index 2 DAE systems, which are necessary for mess_lrcfadi_adi.
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


#include "../call_macro.h"

#define CHECK_FUNCTION(INFO,FNAME,DIFFFUNC,Y1,Y2,DIFF,ERR,EPS)          \
        CALL(DIFFFUNC(Y1,Y2,&DIFF));                                    \
        if(DIFF>EPS){                                                   \
            printf(INFO " " #FNAME " failed\tdiff=%e\n",diff);          \
            ERR++;                                                      \
        }else{                                                          \
            printf(INFO " " #FNAME " passed with diff=%e.\n",diff);     \
        }                                                               \

int main ( int argc, char ** argv) {

    int ret=0, err=0;
    int ops [] = {MESS_OP_NONE,MESS_OP_TRANSPOSE,MESS_OP_HERMITIAN};
    mess_int_t nv, np, i, j, k, cols=2, noshifts=3;
    double diff=.0, eps=1e-6, delta=-0.02;
    mess_vector x1, x2, y1, y2,shifts;
    mess_matrix M, A, G, B, Gt, fullA, MM, X1lift, ApMM, X1, X2, Y1, Y2,Tmp, Tmp1, Tmp2, Tmp3, Tmp4;
    mess_direct solver;
    mess_equation eqn;
    mess_options opt;

    mess_version();
    mess_init();

    if ( argc != 7 ) {
        fprintf(stderr, "usage: %s M.mtx A.mtx G.mtx B.mtx MM.mtx memory_usage\n", argv[0]);
        return 0;
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
    CALL(mess_matrix_init(&MM));
    CALL(mess_matrix_init(&X1lift));
    CALL(mess_matrix_init(&ApMM));
    CALL(mess_matrix_init(&X1));
    CALL(mess_matrix_init(&X2));
    CALL(mess_matrix_init(&Y1));
    CALL(mess_matrix_init(&Y2));
    CALL(mess_matrix_init(&Tmp));
    CALL(mess_matrix_init(&Tmp1));
    CALL(mess_matrix_init(&Tmp2));
    CALL(mess_matrix_init(&Tmp3));
    CALL(mess_matrix_init(&Tmp4));

    //read matrices
    CALL(mess_matrix_read_formated(argv[1], M, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[2], A, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[3], G, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[4], B, MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[5], MM, MESS_CSR));


    //size of linearized navier stokes system for velocity and pressure
    nv = M->rows;
    np = G->cols;

    //build FullA = [A,G;G',0]
    CALL(mess_matrix_ctranspose(G,Gt));
    CALL(mess_matrix_cat(A,G,Gt,NULL,MESS_CSR,fullA));


    //init and alloc vector
    MESS_INIT_VECTORS(&x1,&x2,&y1,&y2,&shifts);
    CALL(mess_vector_alloc(x1,nv,MESS_REAL));
    CALL(mess_vector_alloc(x2,nv,MESS_REAL));
    CALL(mess_vector_alloc(y1,nv,MESS_REAL));
    CALL(mess_vector_alloc(y2,nv,MESS_REAL));
    CALL(mess_vector_alloc(shifts,noshifts,MESS_REAL));

    //init equation structure
    CALL(mess_equation_init(&eqn));

    //init options structure
    CALL(mess_options_init(&opt));
    opt->memory_usage = atoi(argv[6]);


    //build equation structure for Index 2 DAE
    CALL(mess_equation_glyap_dae2(eqn,opt,M,A,G,B,delta));

    /*-----------------------------------------------------------------------------
     *  Test function handles for lrcfadi computation, i.e. size == nv
     *-----------------------------------------------------------------------------*/

    //----------------------test ApEINV----------------------//
    //generate test data
    CALL(mess_matrix_alloc(X1,nv,cols,nv*cols,MESS_DENSE,MESS_REAL));
    for(i=0;i<X1->rows;++i){
        for(j=0;j<X1->cols;++j){
            X1->values[i+X1->ld*j]*=(i%3)*(i%6<3?1:-1)+ (j%3)*(j%6<5?1:-1);
        }
    }

    //RHS for test ApEINV
    CALL(mess_matrix_lift(X1,np,X1lift));

    //generate some random shift values for testing
    CALL(mess_vector_ones(shifts));
    for(k=0;k<shifts->dim;++k){shifts->values[k]*=(k%3 + 1)*((k%2)?1:-1);};


    for(k=0;k<2;++k){
        for(j=0;j<shifts->dim;++j){
            //build A+pMM
            CALL(mess_matrix_copy(MM,ApMM));
            if(MESS_IS_REAL(shifts)){
                CALL(mess_matrix_add(1.0,fullA,shifts->values[j],ApMM));
            }else{
                CALL(mess_matrix_addc(1.0,fullA,shifts->values_cpx[j],ApMM));
            }

            //init solver for comparison
            CALL(mess_direct_init(&solver));
        #ifdef MESS_HAVE_UMFPACK
            CALL(mess_direct_create_umfpack(ApMM,solver));
        #else
            CALL(mess_direct_create_sparse_lu(ApMM,solver));
        #endif

            //init multisolver
            CALL(eqn->ApEINV.generate(eqn,shifts));

            //test for MESS_OP_NONE, MESS_OP_TRANSPOSE, MESS_OP_HERMITIAN
            for(i=0;i<3;++i){
                if(MESS_IS_REAL(shifts)){
                    CALL(eqn->ApEINV.apply(eqn, ops[i], shifts->values[j], j, X1, Y1));;
                }else{
                    CALL(eqn->ApEINV.apply(eqn, ops[i], shifts->values_cpx[j], j, X1, Y1));
                }
                CALL(mess_direct_solvem(ops[i],solver,X1lift,Y2));
                CALL(mess_matrix_sub(Y2,0,nv-1,0,Y2->cols-1,Tmp));
                CHECK_FUNCTION("LRCFADI",ApEINV.apply,mess_matrix_diffnorm,Tmp,Y1,diff,err,eps);
            }

            CALL(eqn->ApEINV.clear(eqn));
            CALL(mess_direct_clear(&solver));
        }

        //perform same tests with complex data
        CALL(mess_matrix_tocomplex(X1));
        CALL(mess_matrix_scalec(1+2*I,X1));

        //generate some complex shifts
        CALL(mess_vector_tocomplex(shifts));
        CALL(mess_vector_scalec(1+I,shifts));

        //RHS for test ApEINV
        CALL(mess_matrix_lift(X1,np,X1lift));
    }

    //----------------------Ex.apply EX.apply----------------------//
    //generate real input data
    CALL(mess_matrix_toreal(X1));
    CALL(mess_matrix_copy(X1,X2));

    CALL(mess_vector_ones(x1));
    for(j=0;j<shifts->dim;++j){x1->values[j]*=(j%3 + 1)*((j%2)?1:-1);};
    CALL(mess_vector_copy(x1,x2));

    //test Ex.apply and EX.apply
    for (j=0;j<2;++j){
        for(i=0;i<3;++i){
            //test EX.apply
            CALL(eqn->EX.apply(eqn,ops[i],X1,Y1));
            CALL( mess_matrix_multiply(ops[i], M, MESS_OP_NONE, X1, Y2));
            CHECK_FUNCTION("LRCFADI",EX.apply,mess_matrix_diffnorm,Y1,Y2,diff,err,eps);
        }

        //perform same tests with complex data
        CALL(mess_matrix_tocomplex(X1));
        CALL(mess_matrix_scalec(1+2*I,X1));
        CALL(mess_matrix_copy(X1,X2));

        CALL(mess_vector_tocomplex(x1));
        CALL(mess_vector_scalec(1+2*I,x1));
        CALL(mess_vector_copy(x1,x2));

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
    CALL(mess_matrix_clear(&MM));
    CALL(mess_matrix_clear(&ApMM));
    CALL(mess_matrix_clear(&X1));
    CALL(mess_matrix_clear(&X1lift));
    CALL(mess_matrix_clear(&X2));
    CALL(mess_matrix_clear(&Y1));
    CALL(mess_matrix_clear(&Y2));
    CALL(mess_matrix_clear(&Tmp));
    CALL(mess_matrix_clear(&Tmp1));
    CALL(mess_matrix_clear(&Tmp2));
    CALL(mess_matrix_clear(&Tmp3));
    CALL(mess_matrix_clear(&Tmp4));


    //clear vectors
    CALL(mess_vector_clear(&x1));
    CALL(mess_vector_clear(&x2));
    CALL(mess_vector_clear(&y1));
    CALL(mess_vector_clear(&y2));
    CALL(mess_vector_clear(&shifts));


    //clear eqn structure
    CALL(mess_equation_clear(&eqn));

    //clear options structure
    CALL(mess_options_clear(&opt));

    return (err>0)?(1):(0);
}















