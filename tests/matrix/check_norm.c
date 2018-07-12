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
 * @file tests/matrix/check_norm.c
 * @brief Check norm functions.
 * @test
 * This function checks functions @ref mess_matrix_norm1, mess_matrix_norm2, mess_matrix_norm2inv, mess_matrix_normf and
 * @ref mess_matrix_norminf defined in norm.c that means it checks if the column-sum (\f$ 1 \f$), row-sum (\f$ \infty \f$),
 * Frobenius-norm of a matrix and the \f$ 2 \f$ norm of a matrix or its inverse.
 *
 * @}
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "mess/config.h"

#include "mess/mess.h"
#include "mess/error_macro.h"

#include "../call_macro.h"


#define CHECKNORM(NORMFUNC,A, NRM, SOL, EPS, ERR){  \
    NORMFUNC(A,&NRM);                           \
    if(fabs(NRM-SOL)>EPS){++ERR;                \
        printf("%s\n",#NORMFUNC);                   \
        printf("Matrix:%s\n",#A);                   \
        printf("Norm:%e\n",NRM);                    \
        printf("Sol:%e\n",SOL);                     \
        printf("Diff:%e\n",fabs(NRM-SOL));          \
        printf("-----------------------\n");        \
    };                                          \
}


const double rand50x50norm1     =   29.8018246341451;
const double rand50x50norminf   =   28.9667636414232;
const double rand50x50normf     =   28.8202850526756;
const double rand50x50norm2     =   25.2249050294185;
const double rand50x50norm2inv  =   29.5091653482019;

const double randc50x50norm1    =   41.6678409812798;
const double randc50x50norminf  =   43.6532647822569;
const double randc50x50normf    =   40.8594662754399;
const double randc50x50norm2    =   35.6331007126333;
const double randc50x50norm2inv =   12.4776860845471;

const double issAnorm2inv       =   2.572750289476989;


int main ( int argc, char **argv){

    mess_init();
    mess_matrix  Adenser,Adensec, Acscr,Acscc, Acsrr, Acsrc, Acoordr, Acoordc, issA;

    double nrm;
    double eps = 1e-12;
    int err=0;
    int ret = 0;

    mess_init();
    mess_error_level = 0;

    if(argc!=4){
        printf("Invalid Number of input arguments!\n");
        return 1;
    }


    CALL(mess_matrix_init(&Adenser));
    CALL(mess_matrix_init(&Adensec));

    CALL(mess_matrix_init(&Acscr));
    CALL(mess_matrix_init(&Acscc));

    CALL(mess_matrix_init(&Acsrr));
    CALL(mess_matrix_init(&Acsrc));

    CALL(mess_matrix_init(&Acoordr));
    CALL(mess_matrix_init(&Acoordc));

    CALL(mess_matrix_init(&issA));

    /*-----------------------------------------------------------------------------
     *  load matrices
     *-----------------------------------------------------------------------------*/

    CALL(mess_matrix_read(argv[1],Adenser));
    CALL(mess_matrix_read(argv[2],Adensec));
    CALL(mess_matrix_read(argv[3],issA));

    CALL(mess_matrix_convert(Adenser,Acsrr,MESS_CSR));
    CALL(mess_matrix_convert(Adenser,Acscr,MESS_CSC));
    CALL(mess_matrix_convert(Adenser,Acoordr,MESS_COORD));

    CALL(mess_matrix_convert(Adensec,Acsrc,MESS_CSR));
    CALL(mess_matrix_convert(Adensec,Acscc,MESS_CSC));
    CALL(mess_matrix_convert(Adensec,Acoordc,MESS_COORD));


    /*-----------------------------------------------------------------------------
     *  check values
     *-----------------------------------------------------------------------------*/

    //
    //  //check for real case dense
    CHECKNORM(mess_matrix_norm1,Adenser, nrm, rand50x50norm1, eps, err);
    CHECKNORM(mess_matrix_norminf,Adenser, nrm, rand50x50norminf, eps, err);
    CHECKNORM(mess_matrix_normf,Adenser, nrm, rand50x50normf, eps, err);
    CHECKNORM(mess_matrix_norm2,Adenser, nrm, rand50x50norm2, eps, err);
    CHECKNORM(mess_matrix_norm2inv,Adenser, nrm, rand50x50norm2inv, eps, err);


    //check for real case csc
    CHECKNORM(mess_matrix_norm1,Acscr, nrm, rand50x50norm1, eps, err);
    CHECKNORM(mess_matrix_norminf,Acscr, nrm, rand50x50norminf, eps, err);
    CHECKNORM(mess_matrix_normf,Acscr, nrm, rand50x50normf, eps, err);
    CHECKNORM(mess_matrix_norm2,Acscr, nrm, rand50x50norm2, eps, err);
    CHECKNORM(mess_matrix_norm2inv,Acscr, nrm, rand50x50norm2inv, eps, err);


    //check for real case csr
    CHECKNORM(mess_matrix_norm1,Acsrr, nrm, rand50x50norm1, eps, err);
    CHECKNORM(mess_matrix_norminf,Acsrr, nrm, rand50x50norminf, eps, err);
    CHECKNORM(mess_matrix_normf,Acsrr, nrm, rand50x50normf, eps, err);
    CHECKNORM(mess_matrix_norm2,Acsrr, nrm, rand50x50norm2, eps, err);
    CHECKNORM(mess_matrix_norm2inv,Acsrr, nrm, rand50x50norm2inv, eps, err);

    //check for real case coord
    CHECKNORM(mess_matrix_norm1,Acoordr, nrm, rand50x50norm1, eps, err);
    CHECKNORM(mess_matrix_norminf,Acoordr, nrm, rand50x50norminf, eps, err);
    CHECKNORM(mess_matrix_normf,Acoordr, nrm, rand50x50normf, eps, err);
    CHECKNORM(mess_matrix_norm2,Acoordr, nrm, rand50x50norm2, eps, err);
    CHECKNORM(mess_matrix_norm2inv,Acoordr, nrm, rand50x50norm2inv, eps, err);

    //check for complex case dense
    CHECKNORM(mess_matrix_norm1,Adensec, nrm, randc50x50norm1, eps, err);
    CHECKNORM(mess_matrix_norminf,Adensec, nrm, randc50x50norminf, eps, err);
    CHECKNORM(mess_matrix_normf,Adensec, nrm, randc50x50normf, eps, err);
    CHECKNORM(mess_matrix_norm2,Adensec, nrm, randc50x50norm2, eps, err);
    CHECKNORM(mess_matrix_norm2inv,Adensec, nrm, randc50x50norm2inv, eps, err);

    //check for complex case csr
    CHECKNORM(mess_matrix_norm1,Acsrc, nrm, randc50x50norm1, eps, err);
    CHECKNORM(mess_matrix_norminf,Acsrc, nrm, randc50x50norminf, eps, err);
    CHECKNORM(mess_matrix_normf,Acsrc, nrm, randc50x50normf, eps, err);
    CHECKNORM(mess_matrix_norm2,Acsrc, nrm, randc50x50norm2, eps, err);
    CHECKNORM(mess_matrix_norm2inv,Acsrc, nrm, randc50x50norm2inv, eps, err);

    //check for complex case csc
    CHECKNORM(mess_matrix_norm1,Acscc, nrm, randc50x50norm1, eps, err);
    CHECKNORM(mess_matrix_norminf,Acscc, nrm, randc50x50norminf, eps, err);
    CHECKNORM(mess_matrix_normf,Acscc, nrm, randc50x50normf, eps, err);
    CHECKNORM(mess_matrix_norm2,Acscc, nrm, randc50x50norm2, eps, err);
    CHECKNORM(mess_matrix_norm2inv,Acscc, nrm, randc50x50norm2inv, eps, err);

    //check for complex case coord
    CHECKNORM(mess_matrix_norm1,Acoordc, nrm, randc50x50norm1, eps, err);
    CHECKNORM(mess_matrix_norminf,Acoordc, nrm, randc50x50norminf, eps, err);
    CHECKNORM(mess_matrix_normf,Acoordc, nrm, randc50x50normf, eps, err);
    CHECKNORM(mess_matrix_norm2,Acoordc, nrm, randc50x50norm2, eps, err);
    CHECKNORM(mess_matrix_norm2inv,Acoordc, nrm, randc50x50norm2inv, eps, err);


    //check mess_matrix_norm2inv for larger matrix for better code coverage
    CHECKNORM(mess_matrix_norm2inv,issA, nrm, issAnorm2inv, eps, err);


    //clear matrices
    mess_matrix_clear(&Adenser);
    mess_matrix_clear(&Adensec);

    mess_matrix_clear(&Acscr);
    mess_matrix_clear(&Acscc);

    mess_matrix_clear(&Acsrr);
    mess_matrix_clear(&Acsrc);

    mess_matrix_clear(&Acoordr);
    mess_matrix_clear(&Acoordc);


    mess_matrix_clear(&issA);


    return err;
}

