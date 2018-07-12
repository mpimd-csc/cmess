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
 * @file tests/matrix/check_mul.c
 * @brief Check the multiplication of two matrices.
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
#define  CALL2(X) ret = (X); tests++; if ( ret == 1) { err ++; fprintf(stderr,"================  test %35s failed =================\n", #X);  } else {fprintf(stderr, "test %35s passed.\n",#X);  }


/*-----------------------------------------------------------------------------
 *  real test data
 *-----------------------------------------------------------------------------*/

double matrix_data[16] = {11,0,31,0,12,22,0,42,0,0,33,43,14,24,0,0};
double matrix_NN[16]  = {121,0,1364,1333,984,1492,372,924,602,1032,1089,1419,442,528,434,1008};
double matrix_NT[16]  = {461,600,341,504,600,1060,0,924,341,0,2050,1419,504,924,1419,3613};
double matrix_TN[16] = { 1082,132,1023,154,132,2392,1806,696,1023,1806,2938,0,154,696,0,772};
double matrix_TT[16] = {121,984,602,442,0,1492,1032,528,1364,372,1089,434,1333,924,1419,1008};

/*-----------------------------------------------------------------------------
 *
 *-----------------------------------------------------------------------------*/
mess_double_cpx_t matrixi_data[16]={11+11*I,0+12*I,31,14*I,12,22+22*I,0,42+24*I,31*I,0,33+33*I,43,14,24+42*I,43*I,I};

mess_double_cpx_t matrix_NNi[16]={121 +        461*I, 0 +        600*I, 1364 +        341*I, 1333 +        504*I,
    984 +        600*I, 1492 +       1060*I, 372 ,  924 +        924*I,
    602 +        341*I, 1032,   1089 +       2050*I, 1419 +       1419*I,
    442 +        518*I, 528 +        948*I,  434 +       1419*I, 1008 +       3613*I};
mess_double_cpx_t matrix_NTi[16]={ 461 +        121*I, 600              ,     341 +       1364*I, 504 +       1333*I,
    600 +        984*I, 1060 +       1492*I,  0 +        372*I, 924 +        924*I,
    341 +        602*I, 0 +       1032*I, 2050 +       1089*I, 1419 +       1419*I,
    504 +        456*I,    924 +        552*I,    1419 +        434*I,    3613 +       1008*I};
mess_double_cpx_t matrix_NHi[16]={  461 -        121*I,       600              ,           341 -       1364*I,        504 -       1333*I,
    600 -        984*I,         1060 -       1492*I,     0 -        372*I,     924 -        924*I,
    341 -        602*I,      0 -       1032*I,     2050 -       1089*I,      1419 -       1419*I,
    504 -        456*I,      924 -        552*I,  1419 -        434*I,  3613 -       1008*I };
mess_double_cpx_t matrix_TNi[16] = { 1082 +        121*I,  132 +        984*I, 1023 +        602*I, 154 +        442*I,
    132   ,  2392 +       1492*I,  1806 +       1032*I,  696 +        528*I,
    1023 +       1364*I,     1806 +        372*I,          2938 +       1089*I,     0 +        434*I,
    154 +       1333*I,  696 +        966*I,  0 +       1462*I,  772 +       1008*I };
mess_double_cpx_t matrix_TTi[16] = { 121 +       1082*I, 984 +        132*I, 602 +       1023*I, 442 +        154*I,
    0 +        132*I,  1492 +       2392*I, 1032 +       1806*I, 528 +        696*I,
    1364 +       1023*I, 372 +       1806*I, 1089 +       2938*I, 434              ,
    1333 +        154*I,924 +        738*I,  1419 +         43*I, 1008 +        772*I };
mess_double_cpx_t matrix_THi[16] = { 121 -       1082*I,      984 -        132*I,    602 -       1023*I,    442 -        154*I,
    0 -        132*I,  1492 -       2392*I,    1032 -       1806*I,    528 -        696*I,
    1364 -       1023*I, 372 -       1806*I,    1089 -       2938*I,    434              ,
    1333 -        154*I, 924 -        738*I,    1419 -         43*I,  1008 -        772*I};
mess_double_cpx_t matrix_NiN[16] ={   121 +       1082*I, 0 +        132*I, 1364 +       1023*I, 1333 +        154*I,
    984 +        132*I, 1492 +       2392*I, 372 +       1806*I, 924 +        738*I,
    602 +       1023*I, 1032 +       1806*I, 1089 +       2938*I, 1419 +         43*I,
    442 +        154*I, 528 +        696*I, 434              ,  1008 +        772*I};
mess_double_cpx_t matrix_NiT[16] ={  461 +        121*I, 600 +        984*I, 341 +        602*I, 504 +        456*I,
    600              ,     1060 +       1492*I,    0 +       1032*I, 924 +        552*I,
    341 +       1364*I,   0 +        372*I, 2050 +       1089*I,  1419 +        434*I,
    504 +       1333*I,  924 +        924*I, 1419 +       1419*I,3613 +       1008*I} ;
mess_double_cpx_t matrix_TiN[16] ={ 1082 +        121*I,132              ,  1023 +       1364*I, 154 +       1333*I,
    132 +        984*I, 2392 +       1492*I, 1806 +        372*I, 696 +        966*I,
    1023 +        602*I, 1806 +       1032*I,2938 +       1089*I, 0 +       1462*I,
    154 +        442*I,  696 +        528*I, 0 +        434*I,772 +       1008*I};
mess_double_cpx_t matrix_TiT[16] ={ 121 +        461*I,  984 +        600*I, 602 +        341*I,442 +        518*I,
    0 +        600*I,  1492 +       1060*I,  1032     , 528 +        948*I,
    1364 +        341*I, 372              , 1089 +       2050*I, 434 +       1419*I,
    1333 +        504*I, 924 +        924*I, 1419 +       1419*I, 1008 +       3613*I};
mess_double_cpx_t matrix_HiN[16] ={ 1082 -        121*I, 132              ,1023 -       1364*I, 154 -       1333*I,
    132 -        984*I, 2392 -       1492*I,   1806 -        372*I, 696 -        966*I,
    1023 -        602*I,1806 -       1032*I, 2938 -       1089*I, 0 -       1462*I,
    154 -        442*I, 696 -        528*I,  0 -        434*I,  772 -       1008*I};
mess_double_cpx_t matrix_HiT[16] ={ 121 -        461*I, 984 -        600*I,  602 -        341*I,442 -        518*I,
    0 -        600*I,  1492 -       1060*I, 1032 , 528 -        948*I,
    1364 -        341*I,  372              ,    1089 -       2050*I,    434 -       1419*I,
    1333 -        504*I,  924 -        924*I,   1419 -       1419*I,  1008 -       3613*I };

mess_double_cpx_t matrix_NiNi[16] = {   0 +       1543*I, -984 +        732*I,  762 +       1364*I,  877 +        658*I,
    984 +        732*I, 0 +       3452*I,   -660 +       1806*I,  372 +       1662*I,
    -762 +       1364*I, 660 +       1806*I, 0 +       4988*I,      985 +       1462*I,
    -891 +        672*I,-438 +       1644*I, -1028 +       1419*I, -1 +       4385*I, } ;
mess_double_cpx_t matrix_NiTi[16] = {   -621 +        242*I, 468 +        984*I, -682 +       1966*I,  350 +       1789*I,
    468 +        984*I, -1332 +       2984*I, -1806 +       1404*I,  186 +       1476*I,
    -682 +       1966*I, -1806 +       1404*I, -888 +       2178*I,  1376 +       1853*I,
    350 +       1789*I,  186 +       1476*I, 1376 +       1853*I,   2840 +       2016*I } ;
mess_double_cpx_t matrix_NiHi[16] = {  1543,    732 +        984*I,        1364 -        762*I,     658 -        877*I,
    732 -        984*I,   3452 ,  1806 +        660*I,      1662 -        372*I,
    1364 +        762*I,  1806 -        660*I,  4988,  1462 -        985*I,
    658 +        877*I,   1662 +        372*I,  1462 +        985*I, 4386  } ;
mess_double_cpx_t matrix_TiNi[16] = { 621 +        242*I,   -468 +        984*I,  682 +       1966*I,  -364 +       1775*I,
    -468 +        984*I,  1332 +       2984*I,  1806 +       1404*I, -252 +       1494*I,
    682 +       1966*I,  1806 +       1404*I,   888 +       2178*I,  -1419 +       1896*I,
    -364 +       1775*I, -252 +       1494*I,  -1419 +       1896*I, -2842 +       2016*I } ;
mess_double_cpx_t matrix_TiTi[16] = { 0 +       1543*I,   984 +        732*I, -762 +       1364*I, -891 +        672*I,
    -984 +        732*I,  0 +       3452*I,  660 +       1806*I, -438 +       1644*I,
    762 +       1364*I, -660 +       1806*I, 0 +       4988*I,  -1028 +       1419*I,
    877 +        658*I,  372 +       1662*I, 985 +       1462*I,  -1 +       4385*I};
mess_double_cpx_t matrix_TiHi[16] = {  242 -        621*I, 984 +        468*I, 1966 -        682*I,1775 +        364*I,
    984 +        468*I, 2984 -       1332*I, 1404 -       1806*I, 1494 +        252*I,
    1966 -        682*I, 1404 -       1806*I,  2178 -        888*I, 1896 +       1419*I,
    1789 +        350*I, 1476 +        186*I,  1853 +       1376*I,  2017 +       2841*I};
mess_double_cpx_t matrix_HiNi[16] = {  1543, 732 +        984*I, 1364 -        762*I, 672 -        891*I,
    732 -        984*I,   3452, 1806 +        660*I,1644 -        438*I,
    1364 +        762*I, 1806 -        660*I, 4988, 1419 -       1028*I,
    672 +        891*I, 1644 +        438*I, 1419 +       1028*I,  4386     } ;
mess_double_cpx_t matrix_HiTi[16] = {   242 +        621*I,  984 -        468*I,  1966 +        682*I, 1775 -        364*I,
    984 -        468*I,   2984 +       1332*I,  1404 +       1806*I, 1494 -        252*I,
    1966 +        682*I, 1404 +       1806*I,2178 +        888*I, 1896 -       1419*I,
    1789 -        350*I, 1476 -        186*I,1853 -       1376*I, 2017 -       2841*I};
mess_double_cpx_t matrix_HiHi[16] = {  0 -       1543*I, 984 -        732*I,-762 -       1364*I, -891 -        672*I,
    -984 -        732*I,  0 -       3452*I,  660 -       1806*I,  -438 -       1644*I,
    762 -       1364*I,  -660 -       1806*I,  0 -       4988*I, -1028 -       1419*I,
    877 -        658*I,  372 -       1662*I, 985 -       1462*I,   -1 -       4385*I};

int check_mul(mess_matrix A, mess_operation_t  opA, mess_matrix B, mess_operation_t  opB,  mess_matrix C, mess_matrix Csoll){
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
    mess_matrix CNN,CNT,CTN,CTT;
    mess_matrix CNNi,CNTi,CTNi,CTTi,C;
    mess_matrix CNHi,CTHi;
    mess_matrix CNiN,CNiT,CTiN,CTiT,CHiN,CHiT;
    mess_matrix CNiNi,CNiTi, CNiHi, CTiNi,CTiTi, CTiHi, CHiNi,CHiTi,CHiHi;


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
    CALL(mess_matrix_init(&CNN));
    CALL(mess_matrix_init(&CNT));
    CALL(mess_matrix_init(&CTN));
    CALL(mess_matrix_init(&CTT));
    CALL(mess_matrix_init(&CNNi));
    CALL(mess_matrix_init(&CNTi));
    CALL(mess_matrix_init(&CTNi));
    CALL(mess_matrix_init(&CTTi));
    CALL(mess_matrix_init(&CNHi));
    CALL(mess_matrix_init(&CTHi));
    CALL(mess_matrix_init(&CNiN));
    CALL(mess_matrix_init(&CNiT));
    CALL(mess_matrix_init(&CTiN));
    CALL(mess_matrix_init(&CTiT));
    CALL(mess_matrix_init(&CHiN));
    CALL(mess_matrix_init(&CHiT));
    CALL(mess_matrix_init(&CNiNi));
    CALL(mess_matrix_init(&CNiTi));
    CALL(mess_matrix_init(&CNiHi));
    CALL(mess_matrix_init(&CTiNi));
    CALL(mess_matrix_init(&CTiTi));
    CALL(mess_matrix_init(&CTiHi));
    CALL(mess_matrix_init(&CHiNi));
    CALL(mess_matrix_init(&CHiTi));
    CALL(mess_matrix_init(&CHiHi));



    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(A,4,4,0,matrix_data,NULL));
    CALL(mess_matrix_dense_from_farray(B,4,4,0,matrix_data,NULL));
    CALL(mess_matrix_dense_from_farray(Ai,4,4,0,NULL,matrixi_data));
    CALL(mess_matrix_dense_from_farray(Bi,4,4,0,NULL,matrixi_data));

    /*-----------------------------------------------------------------------------
     *  Load real results
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(CNN,4,4,0,matrix_NN,NULL));
    CALL(mess_matrix_dense_from_farray(CNT,4,4,0,matrix_NT,NULL));
    CALL(mess_matrix_dense_from_farray(CTN,4,4,0,matrix_TN,NULL));
    CALL(mess_matrix_dense_from_farray(CTT,4,4,0,matrix_TT,NULL));

    /*-----------------------------------------------------------------------------
     *  load real-complex reasults
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(CNNi,4,4,0,NULL,matrix_NNi));
    CALL(mess_matrix_dense_from_farray(CNTi,4,4,0,NULL,matrix_NTi));
    CALL(mess_matrix_dense_from_farray(CNHi,4,4,0,NULL,matrix_NHi));
    CALL(mess_matrix_dense_from_farray(CTNi,4,4,0,NULL,matrix_TNi));
    CALL(mess_matrix_dense_from_farray(CTTi,4,4,0,NULL,matrix_TTi));
    CALL(mess_matrix_dense_from_farray(CTHi,4,4,0,NULL,matrix_THi));

    /*-----------------------------------------------------------------------------
     *  Load complex -real results
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(CNiN,4,4,0,NULL,matrix_NiN));
    CALL(mess_matrix_dense_from_farray(CNiT,4,4,0,NULL,matrix_NiT));
    CALL(mess_matrix_dense_from_farray(CTiN,4,4,0,NULL,matrix_TiN));
    CALL(mess_matrix_dense_from_farray(CTiT,4,4,0,NULL,matrix_TiT));
    CALL(mess_matrix_dense_from_farray(CHiN,4,4,0,NULL,matrix_HiN));
    CALL(mess_matrix_dense_from_farray(CHiT,4,4,0,NULL,matrix_HiT));

    /*-----------------------------------------------------------------------------
     *  Load complex results
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(CNiNi,4,4,0,NULL,matrix_NiNi));
    CALL(mess_matrix_dense_from_farray(CNiTi,4,4,0,NULL,matrix_NiTi));
    CALL(mess_matrix_dense_from_farray(CNiHi,4,4,0,NULL,matrix_NiHi));
    CALL(mess_matrix_dense_from_farray(CTiNi,4,4,0,NULL,matrix_TiNi));
    CALL(mess_matrix_dense_from_farray(CTiTi,4,4,0,NULL,matrix_TiTi));
    CALL(mess_matrix_dense_from_farray(CTiHi,4,4,0,NULL,matrix_TiHi));
    CALL(mess_matrix_dense_from_farray(CHiNi,4,4,0,NULL,matrix_HiNi));
    CALL(mess_matrix_dense_from_farray(CHiTi,4,4,0,NULL,matrix_HiTi));
    CALL(mess_matrix_dense_from_farray(CHiHi,4,4,0,NULL,matrix_HiHi));

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
     *  DENSE Checks REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(A,MESS_OP_NONE,B,MESS_OP_NONE,C,CNN));
    CALL2(check_mul(A,MESS_OP_NONE,B,MESS_OP_TRANSPOSE,C,CNT));
    CALL2(check_mul(A,MESS_OP_NONE,B,MESS_OP_HERMITIAN,C,CNT));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,B,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,B,MESS_OP_HERMITIAN,C,CTT));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,B,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,B,MESS_OP_HERMITIAN,C,CTT));


    /*-----------------------------------------------------------------------------
     *  DENSE REAL-COMPLEX
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(A,MESS_OP_NONE,Bi,MESS_OP_NONE,C,CNNi));
    CALL2(check_mul(A,MESS_OP_NONE,Bi,MESS_OP_TRANSPOSE,C,CNTi));
    CALL2(check_mul(A,MESS_OP_NONE,Bi,MESS_OP_HERMITIAN,C,CNHi));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bi,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bi,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bi,MESS_OP_HERMITIAN,C,CTHi));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bi,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bi,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bi,MESS_OP_HERMITIAN,C,CTHi));

    /*-----------------------------------------------------------------------------
     *  DENSE COMPLEX-REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Ai,MESS_OP_NONE,B,MESS_OP_NONE,C,CNiN));
    CALL2(check_mul(Ai,MESS_OP_NONE,B,MESS_OP_TRANSPOSE,C,CNiT));
    CALL2(check_mul(Ai,MESS_OP_NONE,B,MESS_OP_HERMITIAN,C,CNiT));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,CTiN));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,B,MESS_OP_TRANSPOSE,C,CTiT));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,B,MESS_OP_HERMITIAN,C,CTiT));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,CHiN));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,B,MESS_OP_TRANSPOSE,C,CHiT));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,B,MESS_OP_HERMITIAN,C,CHiT));


    /*-----------------------------------------------------------------------------
     *  Dense Complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Ai,MESS_OP_NONE,Bi,MESS_OP_NONE,C,CNiNi));
    CALL2(check_mul(Ai,MESS_OP_NONE,Bi,MESS_OP_TRANSPOSE,C,CNiTi));
    CALL2(check_mul(Ai,MESS_OP_NONE,Bi,MESS_OP_HERMITIAN,C,CNiHi));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bi,MESS_OP_NONE,C,CTiNi));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bi,MESS_OP_TRANSPOSE,C,CTiTi));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bi,MESS_OP_HERMITIAN,C,CTiHi));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bi,MESS_OP_NONE,C,CHiNi));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bi,MESS_OP_TRANSPOSE,C,CHiTi));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bi,MESS_OP_HERMITIAN,C,CHiHi));


    /*-----------------------------------------------------------------------------
     *  Dense * Spase (CSR)  REAL
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(A,MESS_OP_NONE,Bcsr,MESS_OP_NONE,C,CNN));
    CALL2(check_mul(A,MESS_OP_NONE,Bcsr,MESS_OP_TRANSPOSE,C,CNT));
    CALL2(check_mul(A,MESS_OP_NONE,Bcsr,MESS_OP_HERMITIAN,C,CNT));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_HERMITIAN,C,CTT));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bcsr,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bcsr,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bcsr,MESS_OP_HERMITIAN,C,CTT));

    /*-----------------------------------------------------------------------------
     *  Dense * Spase (CSR)  real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(A,MESS_OP_NONE,Bcsri,MESS_OP_NONE,C,CNNi));
    CALL2(check_mul(A,MESS_OP_NONE,Bcsri,MESS_OP_TRANSPOSE,C,CNTi));
    CALL2(check_mul(A,MESS_OP_NONE,Bcsri,MESS_OP_HERMITIAN,C,CNHi));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_HERMITIAN,C,CTHi));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bcsri,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bcsri,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bcsri,MESS_OP_HERMITIAN,C,CTHi));

    /*-----------------------------------------------------------------------------
     *  Dense * Spase (CSR)  complex real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsr,MESS_OP_NONE,C,CNiN));
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsr,MESS_OP_TRANSPOSE,C,CNiT));
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsr,MESS_OP_HERMITIAN,C,CNiT));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_NONE,C,CTiN));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_TRANSPOSE,C,CTiT));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_HERMITIAN,C,CTiT));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bcsr,MESS_OP_NONE,C,CHiN));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bcsr,MESS_OP_TRANSPOSE,C,CHiT));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bcsr,MESS_OP_HERMITIAN,C,CHiT));

    /*-----------------------------------------------------------------------------
     *  Dense * Spase (CSR)  complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsri,MESS_OP_NONE,C,CNiNi));
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsri,MESS_OP_TRANSPOSE,C,CNiTi));
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsri,MESS_OP_HERMITIAN,C,CNiHi));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_NONE,C,CTiNi));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_TRANSPOSE,C,CTiTi));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_HERMITIAN,C,CTiHi));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bcsri,MESS_OP_NONE,C,CHiNi));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bcsri,MESS_OP_TRANSPOSE,C,CHiTi));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bcsri,MESS_OP_HERMITIAN,C,CHiHi));

    /*-----------------------------------------------------------------------------
     *  Dense * Spase (CSC)
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(A,MESS_OP_NONE,Bcsc,MESS_OP_NONE,C,CNN));
    CALL2(check_mul(A,MESS_OP_NONE,Bcsc,MESS_OP_TRANSPOSE,C,CNT));
    CALL2(check_mul(A,MESS_OP_NONE,Bcsc,MESS_OP_HERMITIAN,C,CNT));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_HERMITIAN,C,CTT));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bcsc,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bcsc,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bcsc,MESS_OP_HERMITIAN,C,CTT));

    /*-----------------------------------------------------------------------------
     *  Dense * Spase (CSC) real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(A,MESS_OP_NONE,Bcsci,MESS_OP_NONE,C,CNNi));
    CALL2(check_mul(A,MESS_OP_NONE,Bcsci,MESS_OP_TRANSPOSE,C,CNTi));
    CALL2(check_mul(A,MESS_OP_NONE,Bcsci,MESS_OP_HERMITIAN,C,CNHi));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(A,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_HERMITIAN,C,CTHi));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bcsci,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bcsci,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(A,MESS_OP_HERMITIAN,Bcsci,MESS_OP_HERMITIAN,C,CTHi));

    /*-----------------------------------------------------------------------------
     *  Dense * Spase (CSC) complex real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsc,MESS_OP_NONE,C,CNiN));
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsc,MESS_OP_TRANSPOSE,C,CNiT));
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsc,MESS_OP_HERMITIAN,C,CNiT));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_NONE,C,CTiN));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_TRANSPOSE,C,CTiT));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_HERMITIAN,C,CTiT));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bcsc,MESS_OP_NONE,C,CHiN));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bcsc,MESS_OP_TRANSPOSE,C,CHiT));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bcsc,MESS_OP_HERMITIAN,C,CHiT));

    /*-----------------------------------------------------------------------------
     *  Dense * Spase (CSC) omplex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsci,MESS_OP_NONE,C,CNiNi));
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsci,MESS_OP_TRANSPOSE,C,CNiTi));
    CALL2(check_mul(Ai,MESS_OP_NONE,Bcsci,MESS_OP_HERMITIAN,C,CNiHi));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_NONE,C,CTiNi));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_TRANSPOSE,C,CTiTi));
    CALL2(check_mul(Ai,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_HERMITIAN,C,CTiHi));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bcsci,MESS_OP_NONE,C,CHiNi));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bcsci,MESS_OP_TRANSPOSE,C,CHiTi));
    CALL2(check_mul(Ai,MESS_OP_HERMITIAN,Bcsci,MESS_OP_HERMITIAN,C,CHiHi));


    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) * Dense real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsr,MESS_OP_NONE,B,MESS_OP_NONE,C,CNN));
    CALL2(check_mul(Acsr,MESS_OP_NONE,B,MESS_OP_TRANSPOSE,C,CNT));
    CALL2(check_mul(Acsr,MESS_OP_NONE,B,MESS_OP_HERMITIAN,C,CNT));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,B,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,B,MESS_OP_HERMITIAN,C,CTT));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,B,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,B,MESS_OP_HERMITIAN,C,CTT));


    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) * Dense real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bi,MESS_OP_NONE,C,CNNi));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bi,MESS_OP_TRANSPOSE,C,CNTi));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bi,MESS_OP_HERMITIAN,C,CNHi));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bi,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bi,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bi,MESS_OP_HERMITIAN,C,CTHi));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bi,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bi,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bi,MESS_OP_HERMITIAN,C,CTHi));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) * Dense complex real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsri,MESS_OP_NONE,B,MESS_OP_NONE,C,CNiN));
    CALL2(check_mul(Acsri,MESS_OP_NONE,B,MESS_OP_TRANSPOSE,C,CNiT));
    CALL2(check_mul(Acsri,MESS_OP_NONE,B,MESS_OP_HERMITIAN,C,CNiT));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,CTiN));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,B,MESS_OP_TRANSPOSE,C,CTiT));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,B,MESS_OP_HERMITIAN,C,CTiT));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,CHiN));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,B,MESS_OP_TRANSPOSE,C,CHiT));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,B,MESS_OP_HERMITIAN,C,CHiT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) * Dense complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bi,MESS_OP_NONE,C,CNiNi));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bi,MESS_OP_TRANSPOSE,C,CNiTi));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bi,MESS_OP_HERMITIAN,C,CNiHi));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bi,MESS_OP_NONE,C,CTiNi));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bi,MESS_OP_TRANSPOSE,C,CTiTi));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bi,MESS_OP_HERMITIAN,C,CTiHi));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bi,MESS_OP_NONE,C,CHiNi));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bi,MESS_OP_TRANSPOSE,C,CHiTi));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bi,MESS_OP_HERMITIAN,C,CHiHi));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) * Dense real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsc,MESS_OP_NONE,B,MESS_OP_NONE,C,CNN));
    CALL2(check_mul(Acsc,MESS_OP_NONE,B,MESS_OP_TRANSPOSE,C,CNT));
    CALL2(check_mul(Acsc,MESS_OP_NONE,B,MESS_OP_HERMITIAN,C,CNT));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,B,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,B,MESS_OP_HERMITIAN,C,CTT));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,B,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,B,MESS_OP_HERMITIAN,C,CTT));


    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) * Dense real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bi,MESS_OP_NONE,C,CNNi));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bi,MESS_OP_TRANSPOSE,C,CNTi));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bi,MESS_OP_HERMITIAN,C,CNHi));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bi,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bi,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bi,MESS_OP_HERMITIAN,C,CTHi));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bi,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bi,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bi,MESS_OP_HERMITIAN,C,CTHi));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) * Dense complex real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsci,MESS_OP_NONE,B,MESS_OP_NONE,C,CNiN));
    CALL2(check_mul(Acsci,MESS_OP_NONE,B,MESS_OP_TRANSPOSE,C,CNiT));
    CALL2(check_mul(Acsci,MESS_OP_NONE,B,MESS_OP_HERMITIAN,C,CNiT));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,B,MESS_OP_NONE,C,CTiN));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,B,MESS_OP_TRANSPOSE,C,CTiT));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,B,MESS_OP_HERMITIAN,C,CTiT));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,B,MESS_OP_NONE,C,CHiN));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,B,MESS_OP_TRANSPOSE,C,CHiT));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,B,MESS_OP_HERMITIAN,C,CHiT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) * Dense complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bi,MESS_OP_NONE,C,CNiNi));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bi,MESS_OP_TRANSPOSE,C,CNiTi));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bi,MESS_OP_HERMITIAN,C,CNiHi));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bi,MESS_OP_NONE,C,CTiNi));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bi,MESS_OP_TRANSPOSE,C,CTiTi));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bi,MESS_OP_HERMITIAN,C,CTiHi));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bi,MESS_OP_NONE,C,CHiNi));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bi,MESS_OP_TRANSPOSE,C,CHiTi));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bi,MESS_OP_HERMITIAN,C,CHiHi));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) * Sparse(CSR) real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsr,MESS_OP_NONE,C,CNN));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsr,MESS_OP_TRANSPOSE,C,CNT));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsr,MESS_OP_HERMITIAN,C,CNT));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_HERMITIAN,C,CTT));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bcsr,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bcsr,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bcsr,MESS_OP_HERMITIAN,C,CTT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) * Sparse(CSR) real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsri,MESS_OP_NONE,C,CNNi));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsri,MESS_OP_TRANSPOSE,C,CNTi));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsri,MESS_OP_HERMITIAN,C,CNHi));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_HERMITIAN,C,CTHi));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bcsri,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bcsri,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bcsri,MESS_OP_HERMITIAN,C,CTHi));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) * Sparse(CSR) complex real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsr,MESS_OP_NONE,C,CNiN));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsr,MESS_OP_TRANSPOSE,C,CNiT));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsr,MESS_OP_HERMITIAN,C,CNiT));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_NONE,C,CTiN));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_TRANSPOSE,C,CTiT));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_HERMITIAN,C,CTiT));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bcsr,MESS_OP_NONE,C,CHiN));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bcsr,MESS_OP_TRANSPOSE,C,CHiT));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bcsr,MESS_OP_HERMITIAN,C,CHiT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) * Sparse(CSR) complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsri,MESS_OP_NONE,C,CNiNi));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsri,MESS_OP_TRANSPOSE,C,CNiTi));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsri,MESS_OP_HERMITIAN,C,CNiHi));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_NONE,C,CTiNi));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_TRANSPOSE,C,CTiTi));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_HERMITIAN,C,CTiHi));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bcsri,MESS_OP_NONE,C,CHiNi));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bcsri,MESS_OP_TRANSPOSE,C,CHiTi));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bcsri,MESS_OP_HERMITIAN,C,CHiHi));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) * Sparse(CSC) real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsc,MESS_OP_NONE,C,CNN));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsc,MESS_OP_TRANSPOSE,C,CNT));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsc,MESS_OP_HERMITIAN,C,CNT));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_HERMITIAN,C,CTT));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bcsc,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bcsc,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bcsc,MESS_OP_HERMITIAN,C,CTT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) * Sparse(CSC) real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsci,MESS_OP_NONE,C,CNNi));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsci,MESS_OP_TRANSPOSE,C,CNTi));
    CALL2(check_mul(Acsr,MESS_OP_NONE,Bcsci,MESS_OP_HERMITIAN,C,CNHi));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(Acsr,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_HERMITIAN,C,CTHi));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bcsci,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bcsci,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(Acsr,MESS_OP_HERMITIAN,Bcsci,MESS_OP_HERMITIAN,C,CTHi));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) * Sparse(CSC) complex real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsc,MESS_OP_NONE,C,CNiN));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsc,MESS_OP_TRANSPOSE,C,CNiT));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsc,MESS_OP_HERMITIAN,C,CNiT));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_NONE,C,CTiN));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_TRANSPOSE,C,CTiT));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_HERMITIAN,C,CTiT));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bcsc,MESS_OP_NONE,C,CHiN));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bcsc,MESS_OP_TRANSPOSE,C,CHiT));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bcsc,MESS_OP_HERMITIAN,C,CHiT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSR) * Sparse(CSC) complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsci,MESS_OP_NONE,C,CNiNi));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsci,MESS_OP_TRANSPOSE,C,CNiTi));
    CALL2(check_mul(Acsri,MESS_OP_NONE,Bcsci,MESS_OP_HERMITIAN,C,CNiHi));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_NONE,C,CTiNi));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_TRANSPOSE,C,CTiTi));
    CALL2(check_mul(Acsri,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_HERMITIAN,C,CTiHi));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bcsci,MESS_OP_NONE,C,CHiNi));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bcsci,MESS_OP_TRANSPOSE,C,CHiTi));
    CALL2(check_mul(Acsri,MESS_OP_HERMITIAN,Bcsci,MESS_OP_HERMITIAN,C,CHiHi));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) * Sparse(CSR) real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsr,MESS_OP_NONE,C,CNN));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsr,MESS_OP_TRANSPOSE,C,CNT));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsr,MESS_OP_HERMITIAN,C,CNT));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_HERMITIAN,C,CTT));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bcsr,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bcsr,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bcsr,MESS_OP_HERMITIAN,C,CTT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) * Sparse(CSR) real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsri,MESS_OP_NONE,C,CNNi));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsri,MESS_OP_TRANSPOSE,C,CNTi));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsri,MESS_OP_HERMITIAN,C,CNHi));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_HERMITIAN,C,CTHi));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bcsri,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bcsri,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bcsri,MESS_OP_HERMITIAN,C,CTHi));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) * Sparse(CSR) complex real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsr,MESS_OP_NONE,C,CNiN));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsr,MESS_OP_TRANSPOSE,C,CNiT));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsr,MESS_OP_HERMITIAN,C,CNiT));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_NONE,C,CTiN));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_TRANSPOSE,C,CTiT));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bcsr,MESS_OP_HERMITIAN,C,CTiT));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bcsr,MESS_OP_NONE,C,CHiN));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bcsr,MESS_OP_TRANSPOSE,C,CHiT));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bcsr,MESS_OP_HERMITIAN,C,CHiT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) * Sparse(CSR) complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsri,MESS_OP_NONE,C,CNiNi));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsri,MESS_OP_TRANSPOSE,C,CNiTi));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsri,MESS_OP_HERMITIAN,C,CNiHi));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_NONE,C,CTiNi));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_TRANSPOSE,C,CTiTi));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bcsri,MESS_OP_HERMITIAN,C,CTiHi));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bcsri,MESS_OP_NONE,C,CHiNi));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bcsri,MESS_OP_TRANSPOSE,C,CHiTi));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bcsri,MESS_OP_HERMITIAN,C,CHiHi));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) * Sparse(CSC) real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsc,MESS_OP_NONE,C,CNN));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsc,MESS_OP_TRANSPOSE,C,CNT));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsc,MESS_OP_HERMITIAN,C,CNT));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_HERMITIAN,C,CTT));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bcsc,MESS_OP_NONE,C,CTN));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bcsc,MESS_OP_TRANSPOSE,C,CTT));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bcsc,MESS_OP_HERMITIAN,C,CTT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) * Sparse(CSC) real complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsci,MESS_OP_NONE,C,CNNi));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsci,MESS_OP_TRANSPOSE,C,CNTi));
    CALL2(check_mul(Acsc,MESS_OP_NONE,Bcsci,MESS_OP_HERMITIAN,C,CNHi));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(Acsc,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_HERMITIAN,C,CTHi));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bcsci,MESS_OP_NONE,C,CTNi));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bcsci,MESS_OP_TRANSPOSE,C,CTTi));
    CALL2(check_mul(Acsc,MESS_OP_HERMITIAN,Bcsci,MESS_OP_HERMITIAN,C,CTHi));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) * Sparse(CSC) complex real
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsc,MESS_OP_NONE,C,CNiN));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsc,MESS_OP_TRANSPOSE,C,CNiT));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsc,MESS_OP_HERMITIAN,C,CNiT));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_NONE,C,CTiN));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_TRANSPOSE,C,CTiT));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bcsc,MESS_OP_HERMITIAN,C,CTiT));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bcsc,MESS_OP_NONE,C,CHiN));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bcsc,MESS_OP_TRANSPOSE,C,CHiT));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bcsc,MESS_OP_HERMITIAN,C,CHiT));

    /*-----------------------------------------------------------------------------
     *  Sparse(CSC) * Sparse(CSC) complex
     *-----------------------------------------------------------------------------*/
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsci,MESS_OP_NONE,C,CNiNi));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsci,MESS_OP_TRANSPOSE,C,CNiTi));
    CALL2(check_mul(Acsci,MESS_OP_NONE,Bcsci,MESS_OP_HERMITIAN,C,CNiHi));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_NONE,C,CTiNi));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_TRANSPOSE,C,CTiTi));
    CALL2(check_mul(Acsci,MESS_OP_TRANSPOSE,Bcsci,MESS_OP_HERMITIAN,C,CTiHi));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bcsci,MESS_OP_NONE,C,CHiNi));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bcsci,MESS_OP_TRANSPOSE,C,CHiTi));
    CALL2(check_mul(Acsci,MESS_OP_HERMITIAN,Bcsci,MESS_OP_HERMITIAN,C,CHiHi));







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
    mess_matrix_clear(&CNN);
    mess_matrix_clear(&CNT);
    mess_matrix_clear(&CTN);
    mess_matrix_clear(&CTT);
    mess_matrix_clear(&CNNi);
    mess_matrix_clear(&CNTi);
    mess_matrix_clear(&CTNi);
    mess_matrix_clear(&CTTi);
    mess_matrix_clear(&CNHi);
    mess_matrix_clear(&CTHi);
    mess_matrix_clear(&CNiN);
    mess_matrix_clear(&CNiT);
    mess_matrix_clear(&CTiN);
    mess_matrix_clear(&CTiT);
    mess_matrix_clear(&CHiN);
    mess_matrix_clear(&CHiT);
    mess_matrix_clear(&CNiNi);
    mess_matrix_clear(&CNiTi);
    mess_matrix_clear(&CNiHi);
    mess_matrix_clear(&CTiNi);
    mess_matrix_clear(&CTiTi);
    mess_matrix_clear(&CTiHi);
    mess_matrix_clear(&CHiNi);
    mess_matrix_clear(&CHiTi);
    mess_matrix_clear(&CHiHi);


    return (err>0)?(1):(0);
}

