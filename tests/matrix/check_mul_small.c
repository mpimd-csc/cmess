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
 * @file tests/matrix/check_mul_small.c
 * @brief Check the multiplication of two matrices (version with small matrices).
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
#define  CALL2(X) ret = (X); tests++; if ( ret == 1) { err ++; fprintf(stderr,"================  test %35s failed =================\n", #X); /* abort(); */  } else {fprintf(stderr, "test %35s passed.\n",#X);  }


/*-----------------------------------------------------------------------------
 *  real test data
 *-----------------------------------------------------------------------------*/
double matrix_B[4*2] = { 11, 21, 0, 41, 0, 22, 32, 42};
double matrix_A[2*4] = { 11, 21,0,22,0,23,14,24};

double matrix_ANBN[4] = { 695, 1677, 588, 2228 };
double matrix_ANAT[4] = { 317, 567, 567, 2030 };
double matrix_ANAH[4] = { 317, 567, 567, 2030 };

double matrix_BTBN[4] = { 2243, 2184, 2184, 3272 };
double matrix_BHBN[4] = { 2243, 2184, 2184, 3272 };
double matrix_BTAT[4] = { 695, 588, 1677, 2228 };
double matrix_BHAT[4] = { 695, 588, 1677, 2228 };
double matrix_BTAH[4] = { 695, 588, 1677, 2228 };
double matrix_BHAH[4] = { 695, 588, 1677, 2228 };


/*-----------------------------------------------------------------------------
 * Mixed Test Data complex-real
 *-----------------------------------------------------------------------------*/
mess_double_cpx_t matrix_Bi[8] = { 11 + 11*I, 21 + 0*I, 0 + 0*I, 41 + 14*I, 0 + 21*I, 22 + 22*I, 32 + 23*I, 42 + 24*I };
mess_double_cpx_t matrix_Ai[8] = { 11 + -22*I, 21 + -0*I, 0 + -42*I, 22 + -44*I, 0 + -0*I, 23 + -64*I, 14 + -82*I, 24 + -84*I };

mess_double_cpx_t matrix_AiNBN[4] = { 695 + -4486*I, 1677 + -4368*I, 588 + -4368*I, 2228 + -6544*I };
mess_double_cpx_t matrix_AiNAT[4] = { 317 + -1390*I, 567 + -1176*I, 567 + -3354*I, 2030 + -4456*I };
mess_double_cpx_t matrix_AiNAH[4] = { 317 + -1390*I, 567 + -1176*I, 567 + -3354*I, 2030 + -4456*I };
mess_double_cpx_t matrix_BiTBN[4] = { 2243 + 695*I, 2184 + 1677*I, 2184 + 588*I, 3272 + 2228*I };
mess_double_cpx_t matrix_BiTAT[4] = { 695 + 317*I, 588 + 567*I, 1677 + 567*I, 2228 + 2030*I };
mess_double_cpx_t matrix_BiTAH[4] = { 695 + 317*I, 588 + 567*I, 1677 + 567*I, 2228 + 2030*I };
mess_double_cpx_t matrix_BiHBN[4] = { 2243 + -695*I, 2184 + -1677*I, 2184 + -588*I, 3272 + -2228*I };
mess_double_cpx_t matrix_BiHAT[4] = { 695 + -317*I, 588 + -567*I, 1677 + -567*I, 2228 + -2030*I };
mess_double_cpx_t matrix_BiHAH[4] = { 695 + -317*I, 588 + -567*I, 1677 + -567*I, 2228 + -2030*I };

/*-----------------------------------------------------------------------------
 *  Real Complex
 *-----------------------------------------------------------------------------*/
mess_double_cpx_t matrix_ANBiN[4] = { 695 + 317*I, 1677 + 567*I, 588 + 567*I, 2228 + 2030*I };
mess_double_cpx_t matrix_ANAiT[4] = { 317 + -1390*I, 567 + -3354*I, 567 + -1176*I, 2030 + -4456*I };
mess_double_cpx_t matrix_ANAiH[4] = { 317 + 1390*I, 567 + 3354*I, 567 + 1176*I, 2030 + 4456*I };
mess_double_cpx_t matrix_BTBiN[4] = { 2243 + 695*I, 2184 + 588*I, 2184 + 1677*I, 3272 + 2228*I };
mess_double_cpx_t matrix_BTAiT[4] = { 695 + -4486*I, 588 + -4368*I, 1677 + -4368*I, 2228 + -6544*I };
mess_double_cpx_t matrix_BTAiH[4] = { 695 + 4486*I, 588 + 4368*I, 1677 + 4368*I, 2228 + 6544*I };
mess_double_cpx_t matrix_BHBiN[4] = { 2243 + 695*I, 2184 + 588*I, 2184 + 1677*I, 3272 + 2228*I };
mess_double_cpx_t matrix_BHAiT[4] = { 695 + -4486*I, 588 + -4368*I, 1677 + -4368*I, 2228 + -6544*I };
mess_double_cpx_t matrix_BHAiH[4] = { 695 + 4486*I, 588 + 4368*I, 1677 + 4368*I, 2228 + 6544*I };

/*-----------------------------------------------------------------------------
 *  complex
 *-----------------------------------------------------------------------------*/
mess_double_cpx_t matrix_AiNBiN[4] = { 2085 + -4169*I, 2853 + -3801*I, 3942 + -3801*I, 6684 + -4514*I };
mess_double_cpx_t matrix_AiNAiT[4] = { -8655 + -2780*I, -8169 + -4530*I, -8169 + -4530*I, -11058 + -8912*I };
mess_double_cpx_t matrix_AiNAiH[4] = { 9289 + 0*I, 9303 + 2178*I, 9303 + -2178*I, 15118 + 0*I };
mess_double_cpx_t matrix_BiTBiN[4] = { 1926 + 1390*I, 1617 + 2265*I, 1617 + 2265*I, 1242 + 4456*I };
mess_double_cpx_t matrix_BiTAiT[4] = { 2085 + -4169*I, 3942 + -3801*I, 2853 + -3801*I, 6684 + -4514*I };
mess_double_cpx_t matrix_BiTAiH[4] = { -695 + 4803*I, -2766 + 4935*I, 501 + 4935*I, -2228 + 8574*I };
mess_double_cpx_t matrix_BiHBiN[4] = { 2560 + 0*I, 2751 + -1089*I, 2751 + 1089*I, 5302 + 0*I };
mess_double_cpx_t matrix_BiHAiT[4] = { -695 + -4803*I, -2766 + -4935*I, 501 + -4935*I, -2228 + -8574*I };
mess_double_cpx_t matrix_BiHAiH[4] = { 2085 + 4169*I, 3942 + 3801*I, 2853 + 3801*I, 6684 + 4514*I };



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

    mess_matrix ANBN, ANAT, ANAH;
    mess_matrix BTBN, BTAT, BTAH;
    mess_matrix BHBN, BHAT, BHAH;

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
    CALL(mess_matrix_init(&ANAH));
    CALL(mess_matrix_init(&BTBN));
    CALL(mess_matrix_init(&BTAT));
    CALL(mess_matrix_init(&BTAH));
    CALL(mess_matrix_init(&BHBN));
    CALL(mess_matrix_init(&BHAT));
    CALL(mess_matrix_init(&BHAH));

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
    CALL(mess_matrix_dense_from_farray(A,2,4,0,matrix_A,NULL));
    CALL(mess_matrix_dense_from_farray(B,4,2,0,matrix_B,NULL));
    CALL(mess_matrix_dense_from_farray(Ai,2,4,0,NULL, matrix_Ai));
    CALL(mess_matrix_dense_from_farray(Bi,4,2,0,NULL, matrix_Bi));

    /*-----------------------------------------------------------------------------
     *  Load real results
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(ANBN,2,2,0,matrix_ANBN,NULL));
    CALL(mess_matrix_dense_from_farray(ANAT,2,2,0,matrix_ANAT,NULL));
    CALL(mess_matrix_dense_from_farray(ANAH,2,2,0,matrix_ANAT,NULL));
    CALL(mess_matrix_dense_from_farray(BTBN,2,2,0,matrix_BTBN,NULL));
    CALL(mess_matrix_dense_from_farray(BTAT,2,2,0,matrix_BTAT,NULL));
    CALL(mess_matrix_dense_from_farray(BTAH,2,2,0,matrix_BTAH,NULL));
    CALL(mess_matrix_dense_from_farray(BHBN,2,2,0,matrix_BHBN,NULL));
    CALL(mess_matrix_dense_from_farray(BHAT,2,2,0,matrix_BHAT,NULL));
    CALL(mess_matrix_dense_from_farray(BHAH,2,2,0,matrix_BHAH,NULL));

    /*-----------------------------------------------------------------------------
     *  load complex-real resultfs
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(AiNBN,2,2,0,NULL,matrix_AiNBN));
    CALL(mess_matrix_dense_from_farray(AiNAT,2,2,0,NULL,matrix_AiNAT));
    CALL(mess_matrix_dense_from_farray(AiNAH,2,2,0,NULL,matrix_AiNAH));
    CALL(mess_matrix_dense_from_farray(BiTBN,2,2,0,NULL,matrix_BiTBN));
    CALL(mess_matrix_dense_from_farray(BiTAT,2,2,0,NULL,matrix_BiTAT));
    CALL(mess_matrix_dense_from_farray(BiTAH,2,2,0,NULL,matrix_BiTAH));
    CALL(mess_matrix_dense_from_farray(BiHBN,2,2,0,NULL,matrix_BiHBN));
    CALL(mess_matrix_dense_from_farray(BiHAT,2,2,0,NULL,matrix_BiHAT));
    CALL(mess_matrix_dense_from_farray(BiHAH,2,2,0,NULL,matrix_BiHAH));

    /*-----------------------------------------------------------------------------
     *  load real-complex resultfs
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(ANBiN,2,2,0,NULL,matrix_ANBiN));
    CALL(mess_matrix_dense_from_farray(ANAiT,2,2,0,NULL,matrix_ANAiT));
    CALL(mess_matrix_dense_from_farray(ANAiH,2,2,0,NULL,matrix_ANAiH));
    CALL(mess_matrix_dense_from_farray(BTBiN,2,2,0,NULL,matrix_BTBiN));
    CALL(mess_matrix_dense_from_farray(BTAiT,2,2,0,NULL,matrix_BTAiT));
    CALL(mess_matrix_dense_from_farray(BTAiH,2,2,0,NULL,matrix_BTAiH));
    CALL(mess_matrix_dense_from_farray(BHBiN,2,2,0,NULL,matrix_BHBiN));
    CALL(mess_matrix_dense_from_farray(BHAiT,2,2,0,NULL,matrix_BHAiT));
    CALL(mess_matrix_dense_from_farray(BHAiH,2,2,0,NULL,matrix_BHAiH));


    /*-----------------------------------------------------------------------------
     *  load complex results
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(AiNBiN,2,2,0,NULL,matrix_AiNBiN));
    CALL(mess_matrix_dense_from_farray(AiNAiT,2,2,0,NULL,matrix_AiNAiT));
    CALL(mess_matrix_dense_from_farray(AiNAiH,2,2,0,NULL,matrix_AiNAiH));
    CALL(mess_matrix_dense_from_farray(BiTBiN,2,2,0,NULL,matrix_BiTBiN));
    CALL(mess_matrix_dense_from_farray(BiTAiT,2,2,0,NULL,matrix_BiTAiT));
    CALL(mess_matrix_dense_from_farray(BiTAiH,2,2,0,NULL,matrix_BiTAiH));
    CALL(mess_matrix_dense_from_farray(BiHBiN,2,2,0,NULL,matrix_BiHBiN));
    CALL(mess_matrix_dense_from_farray(BiHAiT,2,2,0,NULL,matrix_BiHAiT));
    CALL(mess_matrix_dense_from_farray(BiHAiH,2,2,0,NULL,matrix_BiHAiH));

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
    CALL2(check_mul(A,MESS_OP_NONE,A,MESS_OP_HERMITIAN,C,ANAH));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,BTBN));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,A,MESS_OP_TRANSPOSE,C,BTAT));
    CALL2(check_mul(B,MESS_OP_TRANSPOSE,A,MESS_OP_HERMITIAN,C,BTAH));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,BHBN));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,A,MESS_OP_TRANSPOSE,C,BHAT));
    CALL2(check_mul(B,MESS_OP_HERMITIAN,A,MESS_OP_HERMITIAN,C,BHAH));
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
    mess_matrix_printinfo(Bcsri);
    mess_matrix_printinfo(Acsc);
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
    mess_matrix_clear(&ANAH);
    mess_matrix_clear(&BTBN);
    mess_matrix_clear(&BTAT);
    mess_matrix_clear(&BTAH);
    mess_matrix_clear(&BHBN);
    mess_matrix_clear(&BHAT);
    mess_matrix_clear(&BHAH);


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

