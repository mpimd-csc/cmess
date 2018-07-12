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
 * @file tests/matrix/check_add.c
 * @brief Check matrix addition.
 * @author @koehlerm
 * @test
 * This function checks if @ref mess_matrix_add  and @ref mess_matrix_addc works correctly.
 * @}
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <math.h>

#include "../call_macro.h"
#include "mess/mess.h"

#define CHECK_ADD(TOL,P,ROWS,COLS,ALPHA,BETA,DTA,DTB,STA,STB){              \
    int ret;                                                                \
    double diff;                                                            \
    mess_matrix A, B, TMPA, TMPB, COPYA, COPYB;                             \
    MESS_INIT_MATRICES(&A,&B,&TMPA,&TMPB,&COPYA,&COPYB);                    \
                                                                            \
    /*generate random matrix and create reference solution*/                \
    CALL(mess_matrix_rand(A,ROWS,COLS,MESS_DENSE,DTA,P));                   \
    CALL(mess_matrix_rand(B,ROWS,COLS,MESS_DENSE,DTB,P));                   \
                                                                            \
    /* create a copy of the original matrices */                            \
    CALL(mess_matrix_copy(A,COPYA));                                        \
    CALL(mess_matrix_copy(B,COPYB));                                        \
    CALL(mess_matrix_convert(A,TMPA,STA));                                  \
    CALL(mess_matrix_convert(B,TMPB,STB));                                  \
    CALL(mess_matrix_addc(ALPHA,A,BETA,B));                                 \
    CALL(mess_matrix_addc(ALPHA,TMPA,BETA,TMPB));                           \
    CALL(mess_matrix_diffnorm(B,TMPB,&diff));                               \
                                                                            \
    if(diff>TOL){                                                           \
        printf("FAILED\n");                                                 \
        printf("Matrix A\n");                                               \
        mess_matrix_printinfo(TMPA);                                        \
        mess_matrix_print(TMPA);                                            \
        printf("Matrix B\n");                                               \
        mess_matrix_printinfo(TMPB);                                        \
        mess_matrix_print(TMPB);                                            \
        printf("Expected Result:\n");                                       \
        mess_matrix_print(B);                                               \
        printf("Got:\n");                                                   \
        mess_matrix_print(TMPB);                                            \
        printf("Original Matrix A\n");                                      \
        mess_matrix_printinfo(COPYA);                                       \
        mess_matrix_print(COPYA);                                           \
        printf("Original Matrix B\n");                                      \
        mess_matrix_printinfo(COPYB);                                       \
        mess_matrix_print(COPYB);                                           \
        printf("alpha= %e + %e I\n",creal(ALPHA),cimag(ALPHA));             \
        printf("beta = %e + %e I\n",creal(BETA),cimag(BETA));               \
        printf("diff = %e\n",diff);                                         \
        printf("tol  = %e\n",tol);                                          \
        return 1;                                                           \
    }                                                                       \
    MESS_CLEAR_MATRICES(&A,&B,&TMPA,&TMPB,&COPYA,&COPYB);                   \
}



int main (int argc, char *argv[]){
    mess_init();
    mess_error_level=0;
    mess_int_t  rows, cols;
    double tol = sqrt(mess_eps());

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 1) {
        printf("usage: %s \n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  check cases
     *-----------------------------------------------------------------------------*/
    mess_double_cpx_t alpha1=2.0,      beta1=3.0;
    mess_double_cpx_t alpha2=-1.0+2*I, beta2=2.5+3*I;

    mess_int_t idtA,idtB;
    mess_datatype_t dts [] = {MESS_REAL, MESS_COMPLEX};
    mess_int_t istA,istB;
    mess_storage_t sts [] = {MESS_DENSE, MESS_CSR, MESS_CSC, MESS_COORD};
    double ps [] = {0.1, 0.2, 0.5, 0.8, 1.0};
    mess_int_t ip;

    for(ip=0;ip<5;++ip){
        for(rows=1;rows<40;rows=rows+2){
            printf("rows="MESS_PRINTF_INT"\n",rows);
            for(cols=1;cols<40;cols*=2){
                for(istA=0;istA<4;istA++){
                    for(istB=0;istB<4;istB++){
                        for(idtA=0;idtA<2;idtA++){
                            for(idtB=0;idtB<2;idtB++){
                                CHECK_ADD(tol,ps[ip],rows,cols,alpha1,beta1,dts[idtA],dts[idtB],sts[istA],sts[istB]);
                                CHECK_ADD(tol,ps[ip],rows,cols,alpha2,beta2,dts[idtA],dts[idtB],sts[istA],sts[istB]);
                                CHECK_ADD(tol,ps[ip],rows,cols,alpha1,beta2,dts[idtA],dts[idtB],sts[istA],sts[istB]);
                                CHECK_ADD(tol,ps[ip],rows,cols,alpha2,beta1,dts[idtA],dts[idtB],sts[istA],sts[istB]);
                            }
                        }
                    }
                }
            }
        }
    }


    return 0;
}

