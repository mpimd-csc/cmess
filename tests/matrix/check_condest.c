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
 * @file tests/matrix/check_condest.c
 * @brief Check the condition number estimation.
 * @author  @mbehr
 * @test
 * This function check the @ref mess_matrix_condest function defined in condest.c that means it checks if the
 * \f$ 1 \f$-norm condition number of a real square matrix is computed correctly.
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

#define CHECKCONDEST(A,NRM,SOL,ERR,EPS){                                                \
    CALL(mess_matrix_condest(A,&NRM));                                              \
    if(fabs(NRM-SOL)>EPS){                                                          \
        printf("Failed Nrm: %e \t Sol: %e \t Err: %e\n",NRM,SOL,fabs(NRM-SOL));     \
        ++ERR;                                                                      \
    }                                                                               \
}


double rand50x50condest = 1.962107663241329e+03;


int main ( int argc, char ** argv) {
    mess_init();
    int ret=0, err=0;
    double eps = 1e-10;
    double nrm=0,sol=0;

    mess_matrix Ar, Acscr, Acsrr;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if ( argc != 3) {
        printf("Usage: %s matrix.mtx cond(matrix.mat,'fro')\n", argv[0]);
    }

    /*-----------------------------------------------------------------------------
     *  init matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_init(&Ar));
    CALL(mess_matrix_init(&Acscr));
    CALL(mess_matrix_init(&Acsrr));

    /*-----------------------------------------------------------------------------
     *  read matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_read (argv[1], Ar));
    sol=atof(argv[2]);
    CALL(mess_matrix_convert(Ar,Acsrr,MESS_CSR));
    CALL(mess_matrix_convert(Ar,Acscr,MESS_CSC));

    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/
    CHECKCONDEST(Ar,nrm,sol,err,eps);
    CHECKCONDEST(Acscr,nrm,sol,err,eps);
    CHECKCONDEST(Acsrr,nrm,sol,err,eps);

    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/
    mess_matrix_clear(&Ar);
    mess_matrix_clear(&Acsrr);
    mess_matrix_clear(&Acscr);

    return (err>0)?(1):(0);
}
