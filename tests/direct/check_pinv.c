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
 * @addtogroup test_direct
 * @{
 * @file tests/direct/check_pinv.c
 * @brief Check the computation of the pseudoinverse of a matrix.
 * @author @mbehr
 * @test
 *
 * This function checks the @ref mess_matrix_pinv function defined in pinv.c.
 * It checks if
 * * \f$ AP_{inv}A=A\f$
 * * \f$ P_{inv}AP_{inv}=P_{inv}\f$
 * * \f$ (AP_{inv})^{H}=AP_{inv}\f$
 * * \f$ (P_{inv}A)^{H}=P_{inv}A \f$
 * where \f$ P_{inv} \f$ is the pseudoinverse of  \f$A\f$.
 *
 * @}
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"

#include "../call_macro.h"


/**
 * @internal
 * @brief Check if pseudoinverse \f$ P_{inv} \f$ of \f$ A \f$ fullfills \f$  AP_{inv}A=A \f$.
 * */
double check_pseudo1(mess_matrix A, mess_matrix Pinv){
    double err=0;
    mess_matrix M1,M2;
    mess_matrix_init(&M1);
    mess_matrix_init(&M2);
    mess_matrix_multiply(MESS_OP_NONE,A,MESS_OP_NONE,Pinv,M1);
    mess_matrix_multiply(MESS_OP_NONE,M1,MESS_OP_NONE,A,M2);
    mess_matrix_diffnorm(M2,A,&err);
    mess_matrix_clear(&M1);
    mess_matrix_clear(&M2);
    return err;
}

/**
 * @internal
 * @brief Check if pseudoinverse \f$ P_{inv} \f$ of \f$ A \f$ fullfills \f$  P_{inv}AP_{inv}=P_{inv} \f$.
 * */
double check_pseudo2(mess_matrix A, mess_matrix Pinv){
    double err=0;
    mess_matrix M1,M2;
    mess_matrix_init(&M1);
    mess_matrix_init(&M2);
    mess_matrix_multiply(MESS_OP_NONE,Pinv,MESS_OP_NONE,A,M1);
    mess_matrix_multiply(MESS_OP_NONE,M1,MESS_OP_NONE,Pinv,M2);
    mess_matrix_diffnorm(M2,Pinv,&err);
    mess_matrix_clear(&M1);
    mess_matrix_clear(&M2);
    return err;
}

/**
 * @internal
 * @brief Check if pseudoinverse \f$ P_{inv} \f$ of \f$ A \f$ fullfills \f$ (AP_{inv})^H=AP_{inv} \f$.
 * */
double check_pseudo3(mess_matrix A, mess_matrix Pinv){
    double err=0;
    mess_matrix M1,M2;
    mess_matrix_init(&M1);
    mess_matrix_init(&M2);
    mess_matrix_multiply(MESS_OP_NONE,A,MESS_OP_NONE,Pinv,M1);
    mess_matrix_ctranspose(M1,M2);
    mess_matrix_diffnorm(M1,M2,&err);
    mess_matrix_clear(&M1);
    mess_matrix_clear(&M2);
    return err;
}

/**
 * @internal
 * @brief Check if pseudoinverse \f$ P_{inv} \f$ of \f$ A \f$ fullfills \f$ (AP_{inv})^H=AP_{inv} \f$.
 * */
double check_pseudo4(mess_matrix A, mess_matrix Pinv){
    double err=0;
    mess_matrix M1,M2;
    mess_matrix_init(&M1);
    mess_matrix_init(&M2);
    mess_matrix_multiply(MESS_OP_NONE,Pinv,MESS_OP_NONE,A,M1);
    mess_matrix_ctranspose(M1,M2);
    mess_matrix_diffnorm(M1,M2,&err);
    mess_matrix_clear(&M1);
    mess_matrix_clear(&M2);
    return err;
}


int main ( int argc, char ** argv) {

    mess_init();
    int  err=0,ret;
    double eps,diff;
    mess_matrix A,Pinv;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if ( argc != 3) {
        printf("usage: %s  matrix.mtx op\n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     * Init/load matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_read (argv[1], A ));

    //make matrix rectangular
    mess_matrix_resize(A,A->rows,(mess_int_t)((A->cols+0.5)/2));
    if(atoi(argv[2])){
        //take hermitian of A
        mess_matrix AT;
        mess_matrix_init(&AT);
        mess_matrix_ctranspose(A,AT);
        mess_matrix_copy(AT,A);
        mess_matrix_clear(&AT);
    }

    eps = sqrt(mess_eps()*(A->rows)*(A->cols));
    CALL(mess_matrix_init(&Pinv));
    CALL(mess_matrix_pinv(A,Pinv));

    /*-----------------------------------------------------------------------------
     *  check properties of pseudoinverse
     *-----------------------------------------------------------------------------*/
    diff = check_pseudo1(A,Pinv);
    printf("check_pseudo1: eps=%e \t diff=%e\n",eps,diff);
    if(diff>eps){
        printf("check_pseudo1 failed with:\n");
        printf("diff=%e\n",diff);
        printf("eps=%e\n",eps);
        mess_matrix_printinfo(A);
        mess_matrix_printshort(A);
        err++;
    }

    diff = check_pseudo2(A,Pinv);
    printf("check_pseudo2: eps=%e \t diff=%e\n",eps,diff);
    if(diff>eps){
        printf("check_pseudo2 failed with:\n");
        printf("diff=%e\n",diff);
        printf("eps=%e\n",eps);
        mess_matrix_printinfo(A);
        mess_matrix_printshort(A);
        err++;
    }

    diff = check_pseudo3(A,Pinv);
    printf("check_pseudo3: eps=%e \t diff=%e\n",eps,diff);
    if(diff>eps){
        printf("check_pseudo3 failed with:\n");
        printf("diff=%e\n",diff);
        printf("eps=%e\n",eps);
        mess_matrix_printinfo(A);
        mess_matrix_printshort(A);
        err++;
    }

    diff = check_pseudo4(A,Pinv);
    printf("check_pseudo4: eps=%e \t diff=%e\n",eps,diff);
    if(diff>eps){
        printf("check_pseudo4 failed with:\n");
        printf("diff=%e\n",diff);
        printf("eps=%e\n",eps);
        mess_matrix_printinfo(A);
        mess_matrix_printshort(A);
        err++;
    }


    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_clear(&A));
    CALL(mess_matrix_clear(&Pinv));



    return err;
}

