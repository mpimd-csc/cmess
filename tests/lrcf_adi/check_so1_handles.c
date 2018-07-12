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
 * @file tests/lrcf_adi/check_so1_handles.c
 * @brief Check function handles for second order systems (so1).
 * @author  @mbehr
 * @test
 * This function uses the @ref mess_equation_glyap_so1 function defined in equation_glyap_so1.c that means it
 * generates a @ref mess_equation object from a second order system
 * \f[
 * \begin{array}{ccc}
 *     M\ddot{x} + D\dot{x} + K x &=& B u
 * \end{array},
 * \f]
 * which is first transformed to a first order system
 * \f[
 *    \underbrace{\left[
 *    \begin{array}{cc}
 *       -K &  0 \\  0 &  M
 *    \end{array}
 *    \right]}_{\mathcal{E}}
 *    \left[
 *    \begin{array}{c}
 *       \dot{x} \\ \ddot{x}
 *    \end{array}
 *    \right] =
 *    \underbrace{\left[
 *    \begin{array}{cc}
 *       0 & -K \\ -K & -D
 *    \end{array}
 *    \right]}_{\mathcal{A}}
 *    \left[
 *    \begin{array}{c}
 *       x  \\ \dot{x}
 *    \end{array}
 *    \right] +
 *    \underbrace{ \left[
 *    \begin{array}{c}
 *       0 \\ B
 *    \end{array}
 *    \right]}_{\mathcal{B}}u
 *   \f]
 * where \f$ M \f$, \f$ D \f$ and \f$ K \f$ are symmetric matrices of dimension \f$ n \f$
 * and checks if the resulting equation object is solved correctly by function handles defined in equation_glyap_so1.
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
    mess_int_t n, i, j, cols=2;
    double diff=.0, eps=1e-7;
    mess_vector x1, x2, y1, y2, shifts;
    mess_matrix M, D, K, B, fullA, fullM, X1, X2, Y1, Y2, Tmp;
    mess_equation eqn;
    mess_options opt;

    mess_version();
    mess_init();

    if ( argc != 4 ) {
        fprintf(stderr, "usage: %s M.mtx D.mtx K.mtx \n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  Init and read data
     *-----------------------------------------------------------------------------*/
    //init matrices
    CALL(mess_matrix_init(&M));
    CALL(mess_matrix_init(&D));
    CALL(mess_matrix_init(&K));
    CALL(mess_matrix_init(&B));
    CALL(mess_matrix_init(&fullA));
    CALL(mess_matrix_init(&fullM));
    CALL(mess_matrix_init(&X1));
    CALL(mess_matrix_init(&X2));
    CALL(mess_matrix_init(&Y1));
    CALL(mess_matrix_init(&Y2));
    CALL(mess_matrix_init(&Tmp));

    //read matrices
    CALL(mess_matrix_read_formated(argv[1], M, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[2], D, MESS_CSR));
    CALL(mess_matrix_read_formated(argv[3], K, MESS_CSR));

    n= M->rows;

    //init and alloc vector
    MESS_INIT_VECTORS(&x1,&x2,&y1,&y2,&shifts)
    CALL(mess_vector_alloc(x1,2*n,MESS_REAL));
    CALL(mess_vector_alloc(x2,2*n,MESS_REAL));
    CALL(mess_vector_alloc(y1,2*n,MESS_REAL));
    CALL(mess_vector_alloc(y2,2*n,MESS_REAL));
    CALL(mess_vector_alloc(shifts,4,MESS_REAL));

    //init mess_matrix
    CALL(mess_matrix_alloc(X1,2*n,cols, 2*n*cols,MESS_DENSE,MESS_REAL));


    //init equation structure
    CALL(mess_equation_init(&eqn));

    //init options structure
    CALL(mess_options_init(&opt));

    //build equation structure for Second Order Systems
    CALL(mess_matrix_rand(Tmp,n,cols,MESS_DENSE,MESS_REAL,0.5));
    CALL(mess_equation_glyap_so1(eqn,opt,M,D,K,Tmp,1e-6,1e+6));


    //build fullA=[0,-K;-K,-D]
    CALL(mess_matrix_cat(NULL,K,K,D,MESS_CSR,fullA));
    CALL(mess_matrix_scale(-1,fullA));

    //build fullE=[-K,0;0,M]
    CALL(mess_matrix_copy(K,Tmp));
    CALL(mess_matrix_scale(-1,Tmp));
    CALL(mess_matrix_cat(Tmp,NULL,NULL,M,MESS_CSR,fullM));


    //print condition numbers
    mess_matrix_condest(M,&diff);       printf("cond(M)=%e\n",diff);
    mess_matrix_condest(D,&diff);       printf("cond(D)=%e\n",diff);
    mess_matrix_condest(K,&diff);       printf("cond(K)=%e\n",diff);
    mess_matrix_condest(fullA,&diff);   printf("cond(fullA)=%e\n",diff);
    mess_matrix_condest(fullM,&diff);   printf("cond(fullM)=%e\n",diff);


    /*-----------------------------------------------------------------------------
     *  Test function handles
     *-----------------------------------------------------------------------------*/
    //generate input data
    CALL(mess_vector_ones(x1)); CALL(mess_vector_toreal_nowarn(x1));
    for(i=0;i<x1->dim;++i){x1->values[i]*=(i%3)*(i%6<3?1:-1);}
    CALL(mess_vector_copy(x1,x2));

    CALL(mess_matrix_ones(X1)); CALL(mess_matrix_toreal(X1));
    for(i=0;i<X1->rows;++i){
        for(j=0;j<X1->cols;++j){
            X1->values[i+X1->ld*j]*=(i%3)*(i%6<3?1:-1)+ (j%3)*(j%6<5?1:-1);
        }
    }

    CALL(mess_matrix_copy(X1,X2));

    int ops [] = {MESS_OP_NONE,MESS_OP_TRANSPOSE,MESS_OP_HERMITIAN};

    for(j=0;j<2;++j){
        for(i=0;i<3;++i){
            //test AX.apply
            CALL(eqn->AX.apply(eqn,ops[i],X1,Y1));
            CALL( mess_matrix_multiply(ops[i], fullA, MESS_OP_NONE, X2, Y2));
            CHECK_FUNCTION("AX.apply",AX.apply,mess_matrix_diffnorm,Y1,Y2,diff,err,eps);

            //test EX.apply
            CALL(eqn->EX.apply(eqn,ops[i],X1,Y1));
            CALL( mess_matrix_multiply(ops[i], fullM, MESS_OP_NONE, X2, Y2));
            CHECK_FUNCTION("EX.apply",EX.apply,mess_matrix_diffnorm,Y1,Y2,diff,err,eps);

            //test AINV.apply
            //CALL(eqn->AINV.generate(eqn));
            CALL(eqn->AINV.apply(eqn,ops[i],X1,Y1));
            CALL( mess_matrix_multiply(ops[i], fullA, MESS_OP_NONE, Y1, X2));
            CHECK_FUNCTION("AINV",AINV.apply,mess_matrix_diffnorm,X1,X2,diff,err,eps);
            //CALL(eqn->AINV.clear(eqn));

            //test EINV.apply
            //CALL(eqn->EINV.generate(eqn));
            CALL(eqn->EINV.apply(eqn,ops[i],X1,Y1));
            CALL( mess_matrix_multiply(ops[i], fullM, MESS_OP_NONE, Y1, X2));
            CHECK_FUNCTION("EINV",EINV.apply,mess_matrix_diffnorm,X1,X2,diff,err,eps);
            //CALL(eqn->EINV.clear(eqn));

            //test ApEINV
            mess_int_t jj;
            if(j==0){
                double p;
                CALL(mess_vector_ones(shifts));
                for(jj=0;jj<shifts->dim;++jj){shifts->values[jj]*=(jj%3 + 1)*((jj%2)?1:-1);};
                if(eqn->ApEINV.generate){CALL(eqn->ApEINV.generate(eqn,shifts));}//dummy shiftvector
                for(jj=0;jj<shifts->dim;++jj){
                    p = shifts->values[jj];
                    CALL(eqn->ApEINV.apply(eqn,ops[i],p,jj,X1,Y1));
                    CALL( mess_matrix_multiply(ops[i], fullA, MESS_OP_NONE, Y1, X2));
                    CALL( mess_matrix_multiply(ops[i], fullM, MESS_OP_NONE, Y1, Y2));
                    CALL(mess_matrix_addc(1.0,X2,p,Y2));
                    CHECK_FUNCTION("ApEINV",ApEINV.apply,mess_matrix_diffnorm,X1,Y2,diff,err,eps);
                }
                if(eqn->ApEINV.clear){CALL(eqn->ApEINV.clear(eqn))};
            }else{
                mess_double_cpx_t p;
                CALL(mess_vector_ones(shifts)); CALL(mess_vector_tocomplex(shifts));
                for(jj=0;jj<shifts->dim;++jj){shifts->values_cpx[jj]*=(jj%3 +1)*((jj%2)?1:-1) + 1.0*I;};
                if(eqn->ApEINV.generate){CALL(eqn->ApEINV.generate(eqn,shifts));}//dummy shiftvector
                for(jj=0;jj<shifts->dim;++jj){
                    p = shifts->values_cpx[jj];
                    CALL(eqn->ApEINV.apply(eqn,ops[i],p,jj,X1,Y1));
                    CALL( mess_matrix_multiply(ops[i], fullA, MESS_OP_NONE, Y1, X2));
                    CALL( mess_matrix_multiply(ops[i], fullM, MESS_OP_NONE, Y1, Y2));
                    CALL(mess_matrix_addc(1.0,X2,(ops[i]==MESS_OP_NONE ||ops[i]==MESS_OP_TRANSPOSE)? p:conj(p),Y2));
                    CHECK_FUNCTION("ApEINV",ApEINV.apply,mess_matrix_diffnorm,X1,Y2,diff,err,eps);
                }
                if(eqn->ApEINV.clear){CALL(eqn->ApEINV.clear(eqn));}
            }

            CALL(mess_matrix_copy(X1,X2));
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
    CALL(mess_matrix_clear(&D));
    CALL(mess_matrix_clear(&K));
    CALL(mess_matrix_clear(&B));
    CALL(mess_matrix_clear(&fullA));
    CALL(mess_matrix_clear(&fullM));
    CALL(mess_matrix_clear(&X1));
    CALL(mess_matrix_clear(&X2));
    CALL(mess_matrix_clear(&Y1));
    CALL(mess_matrix_clear(&Y2));
    CALL(mess_matrix_clear(&Tmp));

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















