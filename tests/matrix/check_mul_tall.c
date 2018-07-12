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
 * @file tests/matrix/check_mul_tall.c
 * @brief Check the multiplication of two matrices (version with tall matrices).
 * @author @koehlerm
 * @test
 * This function checks the @ref mess_matrix_multiply function defined in mult.c that means it checks if the matrix matrix
 * multiplication
 * \f[ op(A) op(B) \f]
 * is computed correctly.
 * Operations \f$ op (.)\f$ on matrices \f$ A \f$ and \f$ B \f$ can be
 * <ul>
 * <li> \ref MESS_OP_NONE (\f$ op(X)= X \f$),
 * <li> \ref MESS_OP_TRANSPOSE (\f$ op(X)= X^T \f$),
 * <li> \ref MESS_OP_HERMITIAN (\f$ op(X)= X^H \f$).
 * </ul>
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
#define  CALL2(X) ret = (X); tests++; if ( ret == 1) { err ++; fprintf(stderr,"================  test %35s failed =================\n", #X); abort();  } else {fprintf(stderr, "test %35s passed.\n",#X);  }


/*-----------------------------------------------------------------------------
 *  real test data
 *-----------------------------------------------------------------------------*/
double matrix_A[4*2] = { 11, 21, 0, 41, 0, 22, 32, 42};
double matrix_B[2*4] = { 11, 21,0,22,0,23,14,24};
double matrix_ANBN[16] = { 121, 693, 672, 1333, 0, 484, 704, 924, 0, 506, 736, 966, 154, 822, 768, 1582 };
double matrix_ANAT[16] = { 121, 231, 0, 451, 231, 925, 704, 1785, 0, 704, 1024, 1344, 451, 1785, 1344, 3445 };
double matrix_BTBN[16] = { 562, 462, 483, 658, 462, 484, 506, 528, 483, 506, 529, 552, 658, 528, 552, 772 };
double matrix_BTAT[16] = { 121, 0, 0, 154, 693, 484, 506, 822, 672, 704, 736, 768, 1333, 924, 966, 1582 };


/*-----------------------------------------------------------------------------
 * Mixed Test Data complex-real
 *-----------------------------------------------------------------------------*/
mess_double_cpx_t matrix_Ai[8] = { 11 + 11*I, 21 + 0*I, 0 + 0*I, 41 + 14*I, 0 + 21*I, 22 + 22*I, 32 + 23*I, 42 + 24*I };
mess_double_cpx_t matrix_Bi[8] = { 11 + -22*I, 21 + -0*I, 0 + -42*I, 22 + -44*I, 0 + -0*I, 23 + -64*I, 14 + -82*I, 24 + -84*I };
mess_double_cpx_t matrix_AiNBN[16] = { 121 + 562*I, 693 + 462*I, 672 + 483*I, 1333 + 658*I, 0 + 462*I, 484 + 484*I, 704 + 506*I, 924 + 528*I, 0 + 483*I, 506 + 506*I, 736 + 529*I, 966 + 552*I, 154 + 658*I, 822 + 528*I, 768 + 552*I, 1582 + 772*I };
mess_double_cpx_t matrix_AiNAT[16] = { 121 + 121*I, 231 + 0*I, 0 + 0*I, 451 + 154*I, 231 + 693*I, 925 + 484*I, 704 + 506*I, 1785 + 822*I, 0 + 672*I, 704 + 704*I, 1024 + 736*I, 1344 + 768*I, 451 + 1333*I, 1785 + 924*I, 1344 + 966*I, 3445 + 1582*I };
mess_double_cpx_t matrix_AiNAH[16] = { 121 + 121*I, 231 + 0*I, 0 + 0*I, 451 + 154*I, 231 + 693*I, 925 + 484*I, 704 + 506*I, 1785 + 822*I, 0 + 672*I, 704 + 704*I, 1024 + 736*I, 1344 + 768*I, 451 + 1333*I, 1785 + 924*I, 1344 + 966*I, 3445 + 1582*I };
mess_double_cpx_t matrix_BiTBN[16] = { 562 + -242*I, 462 + -1386*I, 483 + -1344*I, 658 + -2666*I, 462 + 0*I, 484 + -968*I, 506 + -1408*I, 528 + -1848*I, 483 + 0*I, 506 + -1012*I, 529 + -1472*I, 552 + -1932*I, 658 + -308*I, 528 + -1644*I, 552 + -1536*I, 772 + -3164*I };
mess_double_cpx_t matrix_BiTAT[16] = { 121 + -242*I, 0 + -462*I, 0 + 0*I, 154 + -902*I, 693 + -462*I, 484 + -1850*I, 506 + -1408*I, 822 + -3570*I, 672 + 0*I, 704 + -1408*I, 736 + -2048*I, 768 + -2688*I, 1333 + -902*I, 924 + -3570*I, 966 + -2688*I, 1582 + -6890*I };
mess_double_cpx_t matrix_BiTAH[16] = { 121 + -242*I, 0 + -462*I, 0 + 0*I, 154 + -902*I, 693 + -462*I, 484 + -1850*I, 506 + -1408*I, 822 + -3570*I, 672 + 0*I, 704 + -1408*I, 736 + -2048*I, 768 + -2688*I, 1333 + -902*I, 924 + -3570*I, 966 + -2688*I, 1582 + -6890*I };
mess_double_cpx_t matrix_BiHBN[16] = { 562 + 242*I, 462 + 1386*I, 483 + 1344*I, 658 + 2666*I, 462 + 0*I, 484 + 968*I, 506 + 1408*I, 528 + 1848*I, 483 + 0*I, 506 + 1012*I, 529 + 1472*I, 552 + 1932*I, 658 + 308*I, 528 + 1644*I, 552 + 1536*I, 772 + 3164*I };
mess_double_cpx_t matrix_BiHAT[16] = { 121 + 242*I, 0 + 462*I, 0 + 0*I, 154 + 902*I, 693 + 462*I, 484 + 1850*I, 506 + 1408*I, 822 + 3570*I, 672 + 0*I, 704 + 1408*I, 736 + 2048*I, 768 + 2688*I, 1333 + 902*I, 924 + 3570*I, 966 + 2688*I, 1582 + 6890*I };
mess_double_cpx_t matrix_BiHAH[16] = { 121 + 242*I, 0 + 462*I, 0 + 0*I, 154 + 902*I, 693 + 462*I, 484 + 1850*I, 506 + 1408*I, 822 + 3570*I, 672 + 0*I, 704 + 1408*I, 736 + 2048*I, 768 + 2688*I, 1333 + 902*I, 924 + 3570*I, 966 + 2688*I, 1582 + 6890*I };

/*-----------------------------------------------------------------------------
 *  Real Complex
 *-----------------------------------------------------------------------------*/
mess_double_cpx_t matrix_ANBiN[16] = { 121 + -242*I, 693 + -462*I, 672 + 0*I, 1333 + -902*I, 0 + -462*I, 484 + -1850*I, 704 + -1408*I, 924 + -3570*I, 0 + 0*I, 506 + -1408*I, 736 + -2048*I, 966 + -2688*I, 154 + -902*I, 822 + -3570*I, 768 + -2688*I, 1582 + -6890*I };
mess_double_cpx_t matrix_ANAiT[16] = { 121 + 121*I, 231 + 693*I, 0 + 672*I, 451 + 1333*I, 231 + 0*I, 925 + 484*I, 704 + 704*I, 1785 + 924*I, 0 + 0*I, 704 + 506*I, 1024 + 736*I, 1344 + 966*I, 451 + 154*I, 1785 + 822*I, 1344 + 768*I, 3445 + 1582*I };
mess_double_cpx_t matrix_ANAiH[16] = { 121 + -121*I, 231 + -693*I, 0 + -672*I, 451 + -1333*I, 231 + 0*I, 925 + -484*I, 704 + -704*I, 1785 + -924*I, 0 + 0*I, 704 + -506*I, 1024 + -736*I, 1344 + -966*I, 451 + -154*I, 1785 + -822*I, 1344 + -768*I, 3445 + -1582*I };
mess_double_cpx_t matrix_BTBiN[16] = { 562 + -242*I, 462 + 0*I, 483 + 0*I, 658 + -308*I, 462 + -1386*I, 484 + -968*I, 506 + -1012*I, 528 + -1644*I, 483 + -1344*I, 506 + -1408*I, 529 + -1472*I, 552 + -1536*I, 658 + -2666*I, 528 + -1848*I, 552 + -1932*I, 772 + -3164*I };
mess_double_cpx_t matrix_BTAiT[16] = { 121 + 562*I, 0 + 462*I, 0 + 483*I, 154 + 658*I, 693 + 462*I, 484 + 484*I, 506 + 506*I, 822 + 528*I, 672 + 483*I, 704 + 506*I, 736 + 529*I, 768 + 552*I, 1333 + 658*I, 924 + 528*I, 966 + 552*I, 1582 + 772*I };
mess_double_cpx_t matrix_BTAiH[16] = { 121 + -562*I, 0 + -462*I, 0 + -483*I, 154 + -658*I, 693 + -462*I, 484 + -484*I, 506 + -506*I, 822 + -528*I, 672 + -483*I, 704 + -506*I, 736 + -529*I, 768 + -552*I, 1333 + -658*I, 924 + -528*I, 966 + -552*I, 1582 + -772*I };
mess_double_cpx_t matrix_BHBiN[16] = { 562 + -242*I, 462 + 0*I, 483 + 0*I, 658 + -308*I, 462 + -1386*I, 484 + -968*I, 506 + -1012*I, 528 + -1644*I, 483 + -1344*I, 506 + -1408*I, 529 + -1472*I, 552 + -1536*I, 658 + -2666*I, 528 + -1848*I, 552 + -1932*I, 772 + -3164*I };
mess_double_cpx_t matrix_BHAiT[16] = { 121 + 562*I, 0 + 462*I, 0 + 483*I, 154 + 658*I, 693 + 462*I, 484 + 484*I, 506 + 506*I, 822 + 528*I, 672 + 483*I, 704 + 506*I, 736 + 529*I, 768 + 552*I, 1333 + 658*I, 924 + 528*I, 966 + 552*I, 1582 + 772*I };
mess_double_cpx_t matrix_BHAiH[16] = { 121 + -562*I, 0 + -462*I, 0 + -483*I, 154 + -658*I, 693 + -462*I, 484 + -484*I, 506 + -506*I, 822 + -528*I, 672 + -483*I, 704 + -506*I, 736 + -529*I, 768 + -552*I, 1333 + -658*I, 924 + -528*I, 966 + -552*I, 1582 + -772*I };

/*-----------------------------------------------------------------------------
 *  complex
 *-----------------------------------------------------------------------------*/
mess_double_cpx_t matrix_AiNBiN[16] = { 363 + 320*I, 693 + 0*I, 672 + 483*I, 1641 + -244*I, 1386 + 0*I, 1452 + -1366*I, 1716 + -902*I, 2568 + -3042*I, 1344 + 483*I, 1914 + -902*I, 2208 + -1519*I, 2502 + -2136*I, 2820 + -244*I, 2670 + -3042*I, 2700 + -2136*I, 4746 + -6118*I };
mess_double_cpx_t matrix_AiNAiT[16] = { -441 + 242*I, -231 + 693*I, -483 + 672*I, -207 + 1487*I, -231 + 693*I, 441 + 968*I, 198 + 1210*I, 1257 + 1746*I, -483 + 672*I, 198 + 1210*I, 495 + 1472*I, 792 + 1734*I, -207 + 1487*I, 1257 + 1746*I, 792 + 1734*I, 2673 + 3164*I };
mess_double_cpx_t matrix_AiNAiH[16] = { 683 + 0*I, 693 + -693*I, 483 + -672*I, 1109 + -1179*I, 693 + 693*I, 1409 + 0*I, 1210 + -198*I, 2313 + -102*I, 483 + 672*I, 1210 + 198*I, 1553 + 0*I, 1896 + -198*I, 1109 + 1179*I, 2313 + 102*I, 1896 + 198*I, 4217 + 0*I };
mess_double_cpx_t matrix_BiTBiN[16] = { 78 + -484*I, -462 + -1386*I, 483 + -1344*I, -1146 + -2974*I, -462 + -1386*I, -3216 + -1936*I, -2310 + -2420*I, -6612 + -3492*I, 483 + -1344*I, -2310 + -2420*I, -3567 + -2944*I, -4824 + -3468*I, -1146 + -2974*I, -6612 + -3492*I, -4824 + -3468*I, -13008 + -6328*I };
mess_double_cpx_t matrix_BiTAiT[16] = { 363 + 320*I, 1386 + 0*I, 1344 + 483*I, 2820 + -244*I, 693 + 0*I, 1452 + -1366*I, 1914 + -902*I, 2670 + -3042*I, 672 + 483*I, 1716 + -902*I, 2208 + -1519*I, 2700 + -2136*I, 1641 + -244*I, 2568 + -3042*I, 2502 + -2136*I, 4746 + -6118*I };
mess_double_cpx_t matrix_BiTAiH[16] = { -121 + -804*I, -1386 + -924*I, -1344 + -483*I, -2512 + -1560*I, 693 + -924*I, -484 + -2334*I, -902 + -1914*I, -1026 + -4098*I, 672 + -483*I, -308 + -1914*I, -736 + -2577*I, -1164 + -3240*I, 1025 + -1560*I, -720 + -4098*I, -570 + -3240*I, -1582 + -7662*I };
mess_double_cpx_t matrix_BiHBiN[16] = { 1046 + 0*I, 1386 + 1386*I, 483 + 1344*I, 2462 + 2358*I, 1386 + -1386*I, 4184 + 0*I, 3322 + 396*I, 7668 + 204*I, 483 + -1344*I, 3322 + -396*I, 4625 + 0*I, 5928 + 396*I, 2462 + -2358*I, 7668 + -204*I, 5928 + -396*I, 14552 + 0*I };
mess_double_cpx_t matrix_BiHAiT[16] = { -121 + 804*I, -1386 + 924*I, -1344 + 483*I, -2512 + 1560*I, 693 + 924*I, -484 + 2334*I, -902 + 1914*I, -1026 + 4098*I, 672 + 483*I, -308 + 1914*I, -736 + 2577*I, -1164 + 3240*I, 1025 + 1560*I, -720 + 4098*I, -570 + 3240*I, -1582 + 7662*I };
mess_double_cpx_t matrix_BiHAiH[16] = { 363 + -320*I, 1386 + 0*I, 1344 + -483*I, 2820 + 244*I, 693 + 0*I, 1452 + 1366*I, 1914 + 902*I, 2670 + 3042*I, 672 + -483*I, 1716 + 902*I, 2208 + 1519*I, 2700 + 2136*I, 1641 + 244*I, 2568 + 3042*I, 2502 + 2136*I, 4746 + 6118*I };




int check_mul(mess_matrix A, mess_operation_t opA, mess_matrix B, mess_operation_t opB,  mess_matrix C, mess_matrix Csoll){
    int ret;
    double res;
    double eps = mess_eps();
    CALL( mess_matrix_multiply(opA, A, opB, B, C));
    CALL(mess_matrix_diffnorm(C,Csoll, &res));
    if ( res > 10 * eps) {
        fprintf(stderr,"->res = %lg\n", res);
        mess_matrix_print(C);
        mess_matrix_print(Csoll);
        return 1;
    }
    else return 0;

}

int main ( int argc, char ** argv) {
    mess_init();
    int ret;
    int err = 0;
    int tests=0;
    mess_matrix A,Acsr,Acsc;
    mess_matrix B,Bcsr,Bcsc;
    mess_matrix Ai,Acsri,Acsci;
    mess_matrix Bi,Bcsri,Bcsci;

    mess_matrix ANBN, ANAT, BTBN, BTAT;
    mess_matrix AiNBN, AiNAT, AiNAH;
    mess_matrix BiTBN, BiTAT, BiTAH;
    mess_matrix BiHBN, BiHAT, BiHAH;
    mess_matrix ANBiN, ANAiT, ANAiH;
    mess_matrix BTBiN, BTAiT, BTAiH;
    mess_matrix BHBiN, BHAiT, BHAiH;
    mess_matrix AiNBiN, AiNAiT, AiNAiH;
    mess_matrix BiTBiN, BiTAiT, BiTAiH;
    mess_matrix BiHBiN, BiHAiT, BiHAiH;

    mess_matrix C;


    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_init(&Acsr));
    CALL(mess_matrix_init(&Acsc));
    CALL(mess_matrix_init(&B));
    CALL(mess_matrix_init(&Bcsr));
    CALL(mess_matrix_init(&Bcsc));
    CALL(mess_matrix_init(&Ai));
    CALL(mess_matrix_init(&Acsri));
    CALL(mess_matrix_init(&Acsci));
    CALL(mess_matrix_init(&Bi));
    CALL(mess_matrix_init(&Bcsri));
    CALL(mess_matrix_init(&Bcsci));

    CALL(mess_matrix_init(&C));
    CALL(mess_matrix_init(&ANBN));
    CALL(mess_matrix_init(&ANAT));
    CALL(mess_matrix_init(&BTBN));
    CALL(mess_matrix_init(&BTAT));

    CALL(mess_matrix_init(&AiNBN));
    CALL(mess_matrix_init(&AiNAT));
    CALL(mess_matrix_init(&AiNAH));
    CALL(mess_matrix_init(&BiTBN));
    CALL(mess_matrix_init(&BiTAT));
    CALL(mess_matrix_init(&BiTAH));
    CALL(mess_matrix_init(&BiHBN));
    CALL(mess_matrix_init(&BiHAT));
    CALL(mess_matrix_init(&BiHAH));

    CALL(mess_matrix_init(&ANBiN));
    CALL(mess_matrix_init(&ANAiT));
    CALL(mess_matrix_init(&ANAiH));
    CALL(mess_matrix_init(&BTBiN));
    CALL(mess_matrix_init(&BTAiT));
    CALL(mess_matrix_init(&BTAiH));
    CALL(mess_matrix_init(&BHBiN));
    CALL(mess_matrix_init(&BHAiT));
    CALL(mess_matrix_init(&BHAiH));

    CALL(mess_matrix_init(&AiNBiN));
    CALL(mess_matrix_init(&AiNAiT));
    CALL(mess_matrix_init(&AiNAiH));
    CALL(mess_matrix_init(&BiTBiN));
    CALL(mess_matrix_init(&BiTAiT));
    CALL(mess_matrix_init(&BiTAiH));
    CALL(mess_matrix_init(&BiHBiN));
    CALL(mess_matrix_init(&BiHAiT));
    CALL(mess_matrix_init(&BiHAiH));

    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(A,4,2,0,matrix_A,NULL));
    CALL(mess_matrix_dense_from_farray(B,2,4,0,matrix_B,NULL));
    CALL(mess_matrix_dense_from_farray(Ai,4,2,0,NULL, matrix_Ai));
    CALL(mess_matrix_dense_from_farray(Bi,2,4,0,NULL, matrix_Bi));

    /*-----------------------------------------------------------------------------
     *  Load real results
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(ANBN,4,4,0,matrix_ANBN,NULL));
    CALL(mess_matrix_dense_from_farray(ANAT,4,4,0,matrix_ANAT,NULL));
    CALL(mess_matrix_dense_from_farray(BTBN,4,4,0,matrix_BTBN,NULL));
    CALL(mess_matrix_dense_from_farray(BTAT,4,4,0,matrix_BTAT,NULL));

    /*-----------------------------------------------------------------------------
     *  load complex-real results
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(AiNBN,4,4,0,NULL,matrix_AiNBN));
    CALL(mess_matrix_dense_from_farray(AiNAT,4,4,0,NULL,matrix_AiNAT));
    CALL(mess_matrix_dense_from_farray(AiNAH,4,4,0,NULL,matrix_AiNAH));
    CALL(mess_matrix_dense_from_farray(BiTBN,4,4,0,NULL,matrix_BiTBN));
    CALL(mess_matrix_dense_from_farray(BiTAT,4,4,0,NULL,matrix_BiTAT));
    CALL(mess_matrix_dense_from_farray(BiTAH,4,4,0,NULL,matrix_BiTAH));
    CALL(mess_matrix_dense_from_farray(BiHBN,4,4,0,NULL,matrix_BiHBN));
    CALL(mess_matrix_dense_from_farray(BiHAT,4,4,0,NULL,matrix_BiHAT));
    CALL(mess_matrix_dense_from_farray(BiHAH,4,4,0,NULL,matrix_BiHAH));

    /*-----------------------------------------------------------------------------
     *  load real-complex results
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(ANBiN,4,4,0,NULL,matrix_ANBiN));
    CALL(mess_matrix_dense_from_farray(ANAiT,4,4,0,NULL,matrix_ANAiT));
    CALL(mess_matrix_dense_from_farray(ANAiH,4,4,0,NULL,matrix_ANAiH));
    CALL(mess_matrix_dense_from_farray(BTBiN,4,4,0,NULL,matrix_BTBiN));
    CALL(mess_matrix_dense_from_farray(BTAiT,4,4,0,NULL,matrix_BTAiT));
    CALL(mess_matrix_dense_from_farray(BTAiH,4,4,0,NULL,matrix_BTAiH));
    CALL(mess_matrix_dense_from_farray(BHBiN,4,4,0,NULL,matrix_BHBiN));
    CALL(mess_matrix_dense_from_farray(BHAiT,4,4,0,NULL,matrix_BHAiT));
    CALL(mess_matrix_dense_from_farray(BHAiH,4,4,0,NULL,matrix_BHAiH));


    /*-----------------------------------------------------------------------------
     *  load complex results
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(AiNBiN,4,4,0,NULL,matrix_AiNBiN));
    CALL(mess_matrix_dense_from_farray(AiNAiT,4,4,0,NULL,matrix_AiNAiT));
    CALL(mess_matrix_dense_from_farray(AiNAiH,4,4,0,NULL,matrix_AiNAiH));
    CALL(mess_matrix_dense_from_farray(BiTBiN,4,4,0,NULL,matrix_BiTBiN));
    CALL(mess_matrix_dense_from_farray(BiTAiT,4,4,0,NULL,matrix_BiTAiT));
    CALL(mess_matrix_dense_from_farray(BiTAiH,4,4,0,NULL,matrix_BiTAiH));
    CALL(mess_matrix_dense_from_farray(BiHBiN,4,4,0,NULL,matrix_BiHBiN));
    CALL(mess_matrix_dense_from_farray(BiHAiT,4,4,0,NULL,matrix_BiHAiT));
    CALL(mess_matrix_dense_from_farray(BiHAiH,4,4,0,NULL,matrix_BiHAiH));


    /*-----------------------------------------------------------------------------
     *  Setup the matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_convert(A,Acsr,MESS_CSR));
    CALL(mess_matrix_convert(A,Acsc,MESS_CSC));
    CALL(mess_matrix_convert(B,Bcsr,MESS_CSR));
    CALL(mess_matrix_convert(B,Bcsc,MESS_CSC));
    CALL(mess_matrix_convert(Ai,Acsri,MESS_CSR));
    CALL(mess_matrix_convert(Ai,Acsci,MESS_CSC));
    CALL(mess_matrix_convert(Bi,Bcsri,MESS_CSR));
    CALL(mess_matrix_convert(Bi,Bcsci,MESS_CSC));



    /*-----------------------------------------------------------------------------
     *  real checks
     *-----------------------------------------------------------------------------*/
    /*-----------------------------------------------------------------------------
     *  DENSE Checks REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(A,MESS_OP_NONE,B,MESS_OP_NONE,C,ANBN));
    CALL2(check_mul(A,MESS_OP_NONE,A,MESS_OP_TRANSPOSE,C,ANAT));
    CALL2(check_mul(A,MESS_OP_NONE,A,MESS_OP_HERMITIAN,C,ANAT));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,A,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,A,MESS_OP_HERMITIAN,C,BTAT));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,A,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,A,MESS_OP_HERMITIAN,C,BTAT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) *DENSE Checks REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsr,MESS_OP_NONE,B,MESS_OP_NONE,C,ANBN));
    CALL2(check_mul(Acsr,MESS_OP_NONE,A,MESS_OP_TRANSPOSE,C,ANAT));
    CALL2(check_mul(Acsr,MESS_OP_NONE,A,MESS_OP_HERMITIAN,C,ANAT));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,A,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,A,MESS_OP_HERMITIAN,C,BTAT));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,A,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,A,MESS_OP_HERMITIAN,C,BTAT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) *DENSE Checks REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsc,MESS_OP_NONE,B,MESS_OP_NONE,C,ANBN));
    CALL2(check_mul(Acsc,MESS_OP_NONE,A,MESS_OP_TRANSPOSE,C,ANAT));
    CALL2(check_mul(Acsc,MESS_OP_NONE,A,MESS_OP_HERMITIAN,C,ANAT));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,A,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,A,MESS_OP_HERMITIAN,C,BTAT));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,A,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,A,MESS_OP_HERMITIAN,C,BTAT));

    /*-----------------------------------------------------------------------------
     *  DENSE*Sparse(CSR) Checks REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(A,MESS_OP_NONE,Bcsr,MESS_OP_NONE,C,ANBN));
    CALL2(check_mul(A,MESS_OP_NONE,Acsr,MESS_OP_TRANSPOSE,C,ANAT));
    CALL2(check_mul(A,MESS_OP_NONE,Acsr,MESS_OP_HERMITIAN,C,ANAT));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Acsr,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Acsr,MESS_OP_HERMITIAN,C,BTAT));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Bcsr,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Acsr,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Acsr,MESS_OP_HERMITIAN,C,BTAT));

    /*-----------------------------------------------------------------------------
     *  DENSE*Sparse(CSC) Checks REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(A,MESS_OP_NONE,Bcsc,MESS_OP_NONE,C,ANBN));
    CALL2(check_mul(A,MESS_OP_NONE,Acsc,MESS_OP_TRANSPOSE,C,ANAT));
    CALL2(check_mul(A,MESS_OP_NONE,Acsc,MESS_OP_HERMITIAN,C,ANAT));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Acsc,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Acsc,MESS_OP_HERMITIAN,C,BTAT));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Bcsc,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Acsc,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Acsc,MESS_OP_HERMITIAN,C,BTAT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR)*Sparse(CSR) Checks REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsr,MESS_OP_NONE,C,ANBN));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Acsr,MESS_OP_TRANSPOSE,C,ANAT));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Acsr,MESS_OP_HERMITIAN,C,ANAT));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Acsr,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Acsr,MESS_OP_HERMITIAN,C,BTAT));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Bcsr,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Acsr,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Acsr,MESS_OP_HERMITIAN,C,BTAT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR)*Sparse(CSC) Checks REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsc,MESS_OP_NONE,C,ANBN));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Acsc,MESS_OP_TRANSPOSE,C,ANAT));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Acsc,MESS_OP_HERMITIAN,C,ANAT));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Acsc,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Acsc,MESS_OP_HERMITIAN,C,BTAT));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Bcsc,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Acsc,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Acsc,MESS_OP_HERMITIAN,C,BTAT));


    /*-----------------------------------------------------------------------------
     *  Sparse(CSC)*Sparse(CSR) Checks REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsr,MESS_OP_NONE,C,ANBN));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Acsr,MESS_OP_TRANSPOSE,C,ANAT));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Acsr,MESS_OP_HERMITIAN,C,ANAT));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Acsr,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Acsr,MESS_OP_HERMITIAN,C,BTAT));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Bcsr,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Acsr,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Acsr,MESS_OP_HERMITIAN,C,BTAT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC)*Sparse(CSC) Checks REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsc,MESS_OP_NONE,C,ANBN));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Acsc,MESS_OP_TRANSPOSE,C,ANAT));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Acsc,MESS_OP_HERMITIAN,C,ANAT));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Acsc,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Acsc,MESS_OP_HERMITIAN,C,BTAT));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Bcsc,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Acsc,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Acsc,MESS_OP_HERMITIAN,C,BTAT));


    /*-----------------------------------------------------------------------------
     *  complex-real checks
     *-----------------------------------------------------------------------------*/
    /*-----------------------------------------------------------------------------
     *  DENSE Checks complex real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Ai,MESS_OP_NONE,B,MESS_OP_NONE,C,AiNBN));
    CALL2(check_mul(Ai,MESS_OP_NONE,A,MESS_OP_TRANSPOSE,C,AiNAT));
    CALL2(check_mul(Ai,MESS_OP_NONE,A,MESS_OP_HERMITIAN,C,AiNAH));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,BiTBN));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,A,MESS_OP_TRANSPOSE,C,BiTAT));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,A,MESS_OP_HERMITIAN,C,BiTAH));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,BiHBN));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,A,MESS_OP_TRANSPOSE,C,BiHAT));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,A,MESS_OP_HERMITIAN,C,BiHAH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) *DENSE Checks complex real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsri,MESS_OP_NONE,B,MESS_OP_NONE,C,AiNBN));
    CALL2(check_mul(Acsri,MESS_OP_NONE,A,MESS_OP_TRANSPOSE,C,AiNAT));
    CALL2(check_mul(Acsri,MESS_OP_NONE,A,MESS_OP_HERMITIAN,C,AiNAH));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,BiTBN));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,A,MESS_OP_TRANSPOSE,C,BiTAT));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,A,MESS_OP_HERMITIAN,C,BiTAH));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,BiHBN));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,A,MESS_OP_TRANSPOSE,C,BiHAT));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,A,MESS_OP_HERMITIAN,C,BiHAH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) *DENSE Checks complex  REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsci,MESS_OP_NONE,B,MESS_OP_NONE,C,AiNBN));
    CALL2(check_mul(Acsci,MESS_OP_NONE,A,MESS_OP_TRANSPOSE,C,AiNAT));
    CALL2(check_mul(Acsci,MESS_OP_NONE,A,MESS_OP_HERMITIAN,C,AiNAH));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,BiTBN));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,A,MESS_OP_TRANSPOSE,C,BiTAT));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,A,MESS_OP_HERMITIAN,C,BiTAH));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,BiHBN));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,A,MESS_OP_TRANSPOSE,C,BiHAT));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,A,MESS_OP_HERMITIAN,C,BiHAH));

    /*-----------------------------------------------------------------------------
     *  DENSE*Sparse(CSR) Checks complex REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsr,MESS_OP_NONE,C,AiNBN));
    CALL2(check_mul(Ai,MESS_OP_NONE,Acsr,MESS_OP_TRANSPOSE,C,AiNAT));
    CALL2(check_mul(Ai,MESS_OP_NONE,Acsr,MESS_OP_HERMITIAN,C,AiNAH));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_NONE,C,BiTBN));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Acsr,MESS_OP_TRANSPOSE,C,BiTAT));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Acsr,MESS_OP_HERMITIAN,C,BiTAH));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Bcsr,MESS_OP_NONE,C,BiHBN));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Acsr,MESS_OP_TRANSPOSE,C,BiHAT));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Acsr,MESS_OP_HERMITIAN,C,BiHAH));

    /*-----------------------------------------------------------------------------
     *  DENSE*Sparse(CSC) Checks complex REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsc,MESS_OP_NONE,C,AiNBN));
    CALL2(check_mul(Ai,MESS_OP_NONE,Acsc,MESS_OP_TRANSPOSE,C,AiNAT));
    CALL2(check_mul(Ai,MESS_OP_NONE,Acsc,MESS_OP_HERMITIAN,C,AiNAH));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_NONE,C,BiTBN));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Acsc,MESS_OP_TRANSPOSE,C,BiTAT));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Acsc,MESS_OP_HERMITIAN,C,BiTAH));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Bcsc,MESS_OP_NONE,C,BiHBN));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Acsc,MESS_OP_TRANSPOSE,C,BiHAT));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Acsc,MESS_OP_HERMITIAN,C,BiHAH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR)*Sparse(CSR) Checks complex  REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsr,MESS_OP_NONE,C,AiNBN));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Acsr,MESS_OP_TRANSPOSE,C,AiNAT));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Acsr,MESS_OP_HERMITIAN,C,AiNAH));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_NONE,C,BiTBN));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Acsr,MESS_OP_TRANSPOSE,C,BiTAT));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Acsr,MESS_OP_HERMITIAN,C,BiTAH));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Bcsr,MESS_OP_NONE,C,BiHBN));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Acsr,MESS_OP_TRANSPOSE,C,BiHAT));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Acsr,MESS_OP_HERMITIAN,C,BiHAH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR)*Sparse(CSC) Checks complex REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsc,MESS_OP_NONE,C,AiNBN));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Acsc,MESS_OP_TRANSPOSE,C,AiNAT));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Acsc,MESS_OP_HERMITIAN,C,AiNAH));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_NONE,C,BiTBN));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Acsc,MESS_OP_TRANSPOSE,C,BiTAT));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Acsc,MESS_OP_HERMITIAN,C,BiTAH));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Bcsc,MESS_OP_NONE,C,BiHBN));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Acsc,MESS_OP_TRANSPOSE,C,BiHAT));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Acsc,MESS_OP_HERMITIAN,C,BiHAH));


    /*-----------------------------------------------------------------------------
     *  Sparse(CSC)*Sparse(CSR) Checks complex REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsr,MESS_OP_NONE,C,AiNBN));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Acsr,MESS_OP_TRANSPOSE,C,AiNAT));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Acsr,MESS_OP_HERMITIAN,C,AiNAH));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_NONE,C,BiTBN));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Acsr,MESS_OP_TRANSPOSE,C,BiTAT));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Acsr,MESS_OP_HERMITIAN,C,BiTAH));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Bcsr,MESS_OP_NONE,C,BiHBN));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Acsr,MESS_OP_TRANSPOSE,C,BiHAT));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Acsr,MESS_OP_HERMITIAN,C,BiHAH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC)*Sparse(CSC) Checks complex  REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsc,MESS_OP_NONE,C,AiNBN));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Acsc,MESS_OP_TRANSPOSE,C,AiNAT));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Acsc,MESS_OP_HERMITIAN,C,AiNAH));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_NONE,C,BiTBN));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Acsc,MESS_OP_TRANSPOSE,C,BiTAT));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Acsc,MESS_OP_HERMITIAN,C,BiTAH));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Bcsc,MESS_OP_NONE,C,BiHBN));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Acsc,MESS_OP_TRANSPOSE,C,BiHAT));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Acsc,MESS_OP_HERMITIAN,C,BiHAH));

    /*-----------------------------------------------------------------------------
     *  real complex  checks
     *-----------------------------------------------------------------------------*/
    /*-----------------------------------------------------------------------------
     *  DENSE Checks real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(A,MESS_OP_NONE,Bi,MESS_OP_NONE,C,ANBiN));
    CALL2(check_mul(A,MESS_OP_NONE,Ai,MESS_OP_TRANSPOSE,C,ANAiT));
    CALL2(check_mul(A,MESS_OP_NONE,Ai,MESS_OP_HERMITIAN,C,ANAiH));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Bi,MESS_OP_NONE,C,BTBiN));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Ai,MESS_OP_TRANSPOSE,C,BTAiT));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Ai,MESS_OP_HERMITIAN,C,BTAiH));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Bi,MESS_OP_NONE,C,BHBiN));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Ai,MESS_OP_TRANSPOSE,C,BHAiT));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Ai,MESS_OP_HERMITIAN,C,BHAiH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) *DENSE Checks real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bi,MESS_OP_NONE,C,ANBiN));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Ai,MESS_OP_TRANSPOSE,C,ANAiT));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Ai,MESS_OP_HERMITIAN,C,ANAiH));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Bi,MESS_OP_NONE,C,BTBiN));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Ai,MESS_OP_TRANSPOSE,C,BTAiT));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Ai,MESS_OP_HERMITIAN,C,BTAiH));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Bi,MESS_OP_NONE,C,BHBiN));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Ai,MESS_OP_TRANSPOSE,C,BHAiT));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Ai,MESS_OP_HERMITIAN,C,BHAiH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) *DENSE Checks real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bi,MESS_OP_NONE,C,ANBiN));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Ai,MESS_OP_TRANSPOSE,C,ANAiT));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Ai,MESS_OP_HERMITIAN,C,ANAiH));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Bi,MESS_OP_NONE,C,BTBiN));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Ai,MESS_OP_TRANSPOSE,C,BTAiT));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Ai,MESS_OP_HERMITIAN,C,BTAiH));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Bi,MESS_OP_NONE,C,BHBiN));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Ai,MESS_OP_TRANSPOSE,C,BHAiT));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Ai,MESS_OP_HERMITIAN,C,BHAiH));

    /*-----------------------------------------------------------------------------
     *  DENSE*Sparse(CSR) Checks real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(A,MESS_OP_NONE,Bcsri,MESS_OP_NONE,C,ANBiN));
    CALL2(check_mul(A,MESS_OP_NONE,Acsri,MESS_OP_TRANSPOSE,C,ANAiT));
    CALL2(check_mul(A,MESS_OP_NONE,Acsri,MESS_OP_HERMITIAN,C,ANAiH));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_NONE,C,BTBiN));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Acsri,MESS_OP_TRANSPOSE,C,BTAiT));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Acsri,MESS_OP_HERMITIAN,C,BTAiH));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Bcsri,MESS_OP_NONE,C,BHBiN));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Acsri,MESS_OP_TRANSPOSE,C,BHAiT));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Acsri,MESS_OP_HERMITIAN,C,BHAiH));

    /*-----------------------------------------------------------------------------
     *  DENSE*Sparse(CSC) Checks real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(A,MESS_OP_NONE,Bcsci,MESS_OP_NONE,C,ANBiN));
    CALL2(check_mul(A,MESS_OP_NONE,Acsci,MESS_OP_TRANSPOSE,C,ANAiT));
    CALL2(check_mul(A,MESS_OP_NONE,Acsci,MESS_OP_HERMITIAN,C,ANAiH));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_NONE,C,BTBiN));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Acsci,MESS_OP_TRANSPOSE,C,BTAiT));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,Acsci,MESS_OP_HERMITIAN,C,BTAiH));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Bcsci,MESS_OP_NONE,C,BHBiN));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Acsci,MESS_OP_TRANSPOSE,C,BHAiT));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,Acsci,MESS_OP_HERMITIAN,C,BHAiH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR)*Sparse(CSR) Checks real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsri,MESS_OP_NONE,C,ANBiN));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Acsri,MESS_OP_TRANSPOSE,C,ANAiT));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Acsri,MESS_OP_HERMITIAN,C,ANAiH));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_NONE,C,BTBiN));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Acsri,MESS_OP_TRANSPOSE,C,BTAiT));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Acsri,MESS_OP_HERMITIAN,C,BTAiH));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Bcsri,MESS_OP_NONE,C,BHBiN));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Acsri,MESS_OP_TRANSPOSE,C,BHAiT));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Acsri,MESS_OP_HERMITIAN,C,BHAiH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR)*Sparse(CSC) Checks real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsci,MESS_OP_NONE,C,ANBiN));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Acsci,MESS_OP_TRANSPOSE,C,ANAiT));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Acsci,MESS_OP_HERMITIAN,C,ANAiH));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_NONE,C,BTBiN));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Acsci,MESS_OP_TRANSPOSE,C,BTAiT));
    CALL2(check_mul(Bcsr,MESS_OP_TRANSPOSE,Acsci,MESS_OP_HERMITIAN,C,BTAiH));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Bcsci,MESS_OP_NONE,C,BHBiN));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Acsci,MESS_OP_TRANSPOSE,C,BHAiT));
    CALL2(check_mul(Bcsr,MESS_OP_HERMITIAN,Acsci,MESS_OP_HERMITIAN,C,BHAiH));


    /*-----------------------------------------------------------------------------
     *  Sparse(CSC)*Sparse(CSR) Checks real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsri,MESS_OP_NONE,C,ANBiN));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Acsri,MESS_OP_TRANSPOSE,C,ANAiT));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Acsri,MESS_OP_HERMITIAN,C,ANAiH));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_NONE,C,BTBiN));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Acsri,MESS_OP_TRANSPOSE,C,BTAiT));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Acsri,MESS_OP_HERMITIAN,C,BTAiH));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Bcsri,MESS_OP_NONE,C,BHBiN));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Acsri,MESS_OP_TRANSPOSE,C,BHAiT));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Acsri,MESS_OP_HERMITIAN,C,BHAiH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC)*Sparse(CSC) Checks real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsci,MESS_OP_NONE,C,ANBiN));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Acsci,MESS_OP_TRANSPOSE,C,ANAiT));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Acsci,MESS_OP_HERMITIAN,C,ANAiH));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_NONE,C,BTBiN));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Acsci,MESS_OP_TRANSPOSE,C,BTAiT));
    CALL2(check_mul(Bcsc,MESS_OP_TRANSPOSE,Acsci,MESS_OP_HERMITIAN,C,BTAiH));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Bcsci,MESS_OP_NONE,C,BHBiN));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Acsci,MESS_OP_TRANSPOSE,C,BHAiT));
    CALL2(check_mul(Bcsc,MESS_OP_HERMITIAN,Acsci,MESS_OP_HERMITIAN,C,BHAiH));

    /*-----------------------------------------------------------------------------
     *  real complex  checks
     *-----------------------------------------------------------------------------*/
    /*-----------------------------------------------------------------------------
     *  DENSE Checks complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Ai,MESS_OP_NONE,Bi,MESS_OP_NONE,C,AiNBiN));
    CALL2(check_mul(Ai,MESS_OP_NONE,Ai,MESS_OP_TRANSPOSE,C,AiNAiT));
    CALL2(check_mul(Ai,MESS_OP_NONE,Ai,MESS_OP_HERMITIAN,C,AiNAiH));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Bi,MESS_OP_NONE,C,BiTBiN));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Ai,MESS_OP_TRANSPOSE,C,BiTAiT));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Ai,MESS_OP_HERMITIAN,C,BiTAiH));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Bi,MESS_OP_NONE,C,BiHBiN));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Ai,MESS_OP_TRANSPOSE,C,BiHAiT));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Ai,MESS_OP_HERMITIAN,C,BiHAiH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) *DENSE Checks complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bi,MESS_OP_NONE,C,AiNBiN));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Ai,MESS_OP_TRANSPOSE,C,AiNAiT));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Ai,MESS_OP_HERMITIAN,C,AiNAiH));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Bi,MESS_OP_NONE,C,BiTBiN));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Ai,MESS_OP_TRANSPOSE,C,BiTAiT));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Ai,MESS_OP_HERMITIAN,C,BiTAiH));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Bi,MESS_OP_NONE,C,BiHBiN));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Ai,MESS_OP_TRANSPOSE,C,BiHAiT));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Ai,MESS_OP_HERMITIAN,C,BiHAiH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) *DENSE Checks complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bi,MESS_OP_NONE,C,AiNBiN));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Ai,MESS_OP_TRANSPOSE,C,AiNAiT));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Ai,MESS_OP_HERMITIAN,C,AiNAiH));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Bi,MESS_OP_NONE,C,BiTBiN));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Ai,MESS_OP_TRANSPOSE,C,BiTAiT));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Ai,MESS_OP_HERMITIAN,C,BiTAiH));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Bi,MESS_OP_NONE,C,BiHBiN));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Ai,MESS_OP_TRANSPOSE,C,BiHAiT));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Ai,MESS_OP_HERMITIAN,C,BiHAiH));

    /*-----------------------------------------------------------------------------
     *  DENSE*Sparse(CSR) Checks complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsri,MESS_OP_NONE,C,AiNBiN));
    CALL2(check_mul(Ai,MESS_OP_NONE,Acsri,MESS_OP_TRANSPOSE,C,AiNAiT));
    CALL2(check_mul(Ai,MESS_OP_NONE,Acsri,MESS_OP_HERMITIAN,C,AiNAiH));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_NONE,C,BiTBiN));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Acsri,MESS_OP_TRANSPOSE,C,BiTAiT));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Acsri,MESS_OP_HERMITIAN,C,BiTAiH));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Bcsri,MESS_OP_NONE,C,BiHBiN));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Acsri,MESS_OP_TRANSPOSE,C,BiHAiT));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Acsri,MESS_OP_HERMITIAN,C,BiHAiH));

    /*-----------------------------------------------------------------------------
     *  DENSE*Sparse(CSC) Checks complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsci,MESS_OP_NONE,C,AiNBiN));
    CALL2(check_mul(Ai,MESS_OP_NONE,Acsci,MESS_OP_TRANSPOSE,C,AiNAiT));
    CALL2(check_mul(Ai,MESS_OP_NONE,Acsci,MESS_OP_HERMITIAN,C,AiNAiH));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_NONE,C,BiTBiN));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Acsci,MESS_OP_TRANSPOSE,C,BiTAiT));
    CALL2(check_mul(Bi,MESS_OP_TRANSPOSE,Acsci,MESS_OP_HERMITIAN,C,BiTAiH));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Bcsci,MESS_OP_NONE,C,BiHBiN));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Acsci,MESS_OP_TRANSPOSE,C,BiHAiT));
    CALL2(check_mul(Bi,MESS_OP_HERMITIAN,Acsci,MESS_OP_HERMITIAN,C,BiHAiH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR)*Sparse(CSR) Checks  complex
     *----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsri,MESS_OP_NONE,C,AiNBiN));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Acsri,MESS_OP_TRANSPOSE,C,AiNAiT));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Acsri,MESS_OP_HERMITIAN,C,AiNAiH));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_NONE,C,BiTBiN));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Acsri,MESS_OP_TRANSPOSE,C,BiTAiT));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Acsri,MESS_OP_HERMITIAN,C,BiTAiH));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Bcsri,MESS_OP_NONE,C,BiHBiN));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Acsri,MESS_OP_TRANSPOSE,C,BiHAiT));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Acsri,MESS_OP_HERMITIAN,C,BiHAiH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR)*Sparse(CSC) Checks  complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsci,MESS_OP_NONE,C,AiNBiN));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Acsci,MESS_OP_TRANSPOSE,C,AiNAiT));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Acsci,MESS_OP_HERMITIAN,C,AiNAiH));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_NONE,C,BiTBiN));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Acsci,MESS_OP_TRANSPOSE,C,BiTAiT));
    CALL2(check_mul(Bcsri,MESS_OP_TRANSPOSE,Acsci,MESS_OP_HERMITIAN,C,BiTAiH));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Bcsci,MESS_OP_NONE,C,BiHBiN));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Acsci,MESS_OP_TRANSPOSE,C,BiHAiT));
    CALL2(check_mul(Bcsri,MESS_OP_HERMITIAN,Acsci,MESS_OP_HERMITIAN,C,BiHAiH));


    /*-----------------------------------------------------------------------------
     *  Sparse(CSC)*Sparse(CSR) Checks  complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsri,MESS_OP_NONE,C,AiNBiN));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Acsri,MESS_OP_TRANSPOSE,C,AiNAiT));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Acsri,MESS_OP_HERMITIAN,C,AiNAiH));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_NONE,C,BiTBiN));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Acsri,MESS_OP_TRANSPOSE,C,BiTAiT));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Acsri,MESS_OP_HERMITIAN,C,BiTAiH));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Bcsri,MESS_OP_NONE,C,BiHBiN));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Acsri,MESS_OP_TRANSPOSE,C,BiHAiT));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Acsri,MESS_OP_HERMITIAN,C,BiHAiH));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC)*Sparse(CSC) Checks  complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsci,MESS_OP_NONE,C,AiNBiN));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Acsci,MESS_OP_TRANSPOSE,C,AiNAiT));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Acsci,MESS_OP_HERMITIAN,C,AiNAiH));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_NONE,C,BiTBiN));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Acsci,MESS_OP_TRANSPOSE,C,BiTAiT));
    CALL2(check_mul(Bcsci,MESS_OP_TRANSPOSE,Acsci,MESS_OP_HERMITIAN,C,BiTAiH));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Bcsci,MESS_OP_NONE,C,BiHBiN));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Acsci,MESS_OP_TRANSPOSE,C,BiHAiT));
    CALL2(check_mul(Bcsci,MESS_OP_HERMITIAN,Acsci,MESS_OP_HERMITIAN,C,BiHAiH));


    fprintf(stderr," %d of %d tests passed.\nfailed: %d\n",tests-err,tests,err);


    mess_matrix_clear(&A);
    mess_matrix_clear(&Acsr);
    mess_matrix_clear(&Acsc);
    mess_matrix_clear(&B);
    mess_matrix_clear(&Bcsr);
    mess_matrix_clear(&Bcsc);
    mess_matrix_clear(&Ai);
    mess_matrix_clear(&Acsri);
    mess_matrix_clear(&Acsci);
    mess_matrix_clear(&Bi);
    mess_matrix_clear(&Bcsri);
    mess_matrix_clear(&Bcsci);


    mess_matrix_clear(&C);
    mess_matrix_clear(&ANBN);
    mess_matrix_clear(&ANAT);
    mess_matrix_clear(&BTBN);
    mess_matrix_clear(&BTAT);

    mess_matrix_clear(&AiNBN);
    mess_matrix_clear(&AiNAT);
    mess_matrix_clear(&AiNAH);
    mess_matrix_clear(&BiTBN);
    mess_matrix_clear(&BiTAT);
    mess_matrix_clear(&BiTAH);
    mess_matrix_clear(&BiHBN);
    mess_matrix_clear(&BiHAT);
    mess_matrix_clear(&BiHAH);

    mess_matrix_clear(&ANBiN);
    mess_matrix_clear(&ANAiT);
    mess_matrix_clear(&ANAiH);
    mess_matrix_clear(&BTBiN);
    mess_matrix_clear(&BTAiT);
    mess_matrix_clear(&BTAiH);
    mess_matrix_clear(&BHBiN);
    mess_matrix_clear(&BHAiT);
    mess_matrix_clear(&BHAiH);

    mess_matrix_clear(&AiNBiN);
    mess_matrix_clear(&AiNAiT);
    mess_matrix_clear(&AiNAiH);
    mess_matrix_clear(&BiTBiN);
    mess_matrix_clear(&BiTAiT);
    mess_matrix_clear(&BiTAiH);
    mess_matrix_clear(&BiHBiN);
    mess_matrix_clear(&BiHAiT);
    mess_matrix_clear(&BiHAiH);


    return (err>0)?(1):(0);
}

