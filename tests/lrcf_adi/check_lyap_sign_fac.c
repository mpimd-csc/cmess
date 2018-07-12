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
 *
 * @addtogroup test_lrcfadi
 * @{
 *
 * @file tests/lrcf_adi/check_lyap_sign_fac.c
 * @brief Check the computation of the solution of Lyapunov Equations using matrix sign function and its residual.
 * @test
 *
 * This function checks the solution of the generalized Lyapunov Equation
 * \f[ AX E^T + EXA^T + BB^T = 0 \f]
 * and the Lyapunov Equations
 * \f[ AX + XA^T + BB^T = 0 \f]
 * by computing the corresponding resiudal. \n
 * Time is measured for computing the solution of one Lyapunov Equation and its residual.
 *
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "../call_macro.h"

int  main(int argc , char **argv) {
    mess_init();
    mess_error_level = 3;
    mess_version();

    //define parameters for lypa_sgn_fac
    int ret;
    mess_int_t maxit = 50;
    double tol = 1e-8;
    double tstart=0, tend=0, trun=0, nrm, nrmB;
    mess_matrix A,AT,B,E,ET,Z, ZT,BT;
    int err = 0;

    if ( argc != 4  ) {
        printf("Usage: %s A.mtx B.mtx E.mtx \n", argv[0]);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  init and read data
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&BT,&A,&B,&E,&ET,&AT,&Z,&ZT);
    CALL(mess_matrix_read_formated(argv[1], A,MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[2], B,MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[3], E,MESS_DENSE));
    CALL(mess_matrix_ctranspose(B, BT));
    CALL(mess_matrix_ctranspose(A, AT));
    CALL(mess_matrix_ctranspose(E, ET));

    //compute rhs two norm
    ret = mess_matrix_dynorm2(B, &nrmB);
    printf("nrmB = %lg\n", nrmB);

    /*-----------------------------------------------------------------------------
     * computation
     *-----------------------------------------------------------------------------*/
    printf(" Computation of Lyapunov residual form Matrix Sign Function via Arnoldi\n");
    printf("========================================================================\n");
    printf("Instance: %s %s %s \n", argv[1], argv[2], argv[3]);

    mess_sign_scale_t scales [] ={MESS_SIGN_SCALE_NONE, MESS_SIGN_SCALE_FRO, MESS_SIGN_SCALE_DET};
    mess_int_t scale;
    for(scale=0;scale<3;scale++){
        maxit=50;
        tol=1e-8;
        printf("Scaling Method:%s\n",mess_sign_scale_t_str(scales[scale]));
        tstart=mess_wtime();
        CALL(mess_lrcf_gsignfac(AT, ET, BT,&maxit, &tol,ZT, scales[scale]));
        mess_matrix_ctranspose(ZT,Z);
        CALL(mess_lrcfadi_res2g(A, E, B, Z, &nrm));
        tend = mess_wtime();
        trun = tend-tstart;
        printf("Generalized Lyapunov residual from Matrix Sign Function arnoldi:\n\tIterations: "MESS_PRINTF_INT"\n\trelative: %.15e \n\tabsolute: %.15e\n",maxit, nrm/nrmB, nrm);
        printf("Time for computation: %.4e \n",trun);

        if (nrm/nrmB > 1e-8){
            err++;
            printf("sign test failed with scaling %s\n",mess_sign_scale_t_str(scales[scale]));
        }

        maxit = 50;
        tol = 1e-8;
        tstart=mess_wtime();
        CALL( mess_lrcf_signfac(AT, BT,&maxit,&tol ,ZT,scales[scale]));
        CALL(mess_matrix_ctranspose(ZT,Z));
        CALL(mess_lrcfadi_res2 ( A, B, Z, &nrm));
        tend = mess_wtime();
        trun = tend-tstart;
        printf("Lyapunov residual from Matrix Sign Function arnoldi:\n\tIterations: "MESS_PRINTF_INT"\n\trelative:%.15e \n\tabsolute: %.15e\n",maxit, nrm/nrmB, nrm);
        printf("Time for computation: %.4e \n",trun);
        if (nrm/nrmB > 1e-8){
            err++;
            printf("sign test failed with scaling %s\n",mess_sign_scale_t_str(scales[scale]));
        }
        printf("\n");
    }
    /*-----------------------------------------------------------------------------
     * clear Matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&B,&E,&Z,&AT,&BT,&ET,&ZT);

    if (err ) {
        printf("At least one sign function test failed.\n");
    }
    mess_exit();
    return (err);
}
