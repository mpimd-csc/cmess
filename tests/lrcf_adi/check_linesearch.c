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
 * @file tests/lrcf_adi/check_linesearch.c
 * @brief Check the linesearch function defined in lrnm_linesearch.c
 * @test
 *
 * This function checks the @ref mess_lrcfadi_nm_linesearch_poly and @ref mess_lrcfadi_nm_linesearch_lambda
 * functions defined in lrnm_linesearch.c that means it checks if the
 * scalars of the polynomial are computed correctly and the minimizer \f$ lambda \f$ of the
 * polynomial was found correct. The check is done with tests data generated from
 * the matlab implementation.
 *
 * @}
 */



#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "../call_macro.h"


int main ( int argc, char **argv){
    mess_version();
    mess_init();
    mess_error_level = 2;
    int ret =0, err=0;
    double tol = sqrt(mess_eps());

    //correct coefficent and minimizer computed by matlab
    double alpha, beta, delta, gamma, epsilon, zeta, lambda;

    //input matrices for linesearch
    mess_matrix W, Wnew, DeltaK, DeltaKnew;

    //computed coefficient  and minimizer by C
    double _alpha, _beta, _delta, _gamma, _epsilon, _zeta, _lambda;
    double abs_err, rel_err;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 12) {
        fprintf(stderr, "usage: %s W.mtx Wnew.mtx DeltaK.mtx DeltaKnew.mtx  alpha beta delta gamma epsilon zeta lambda\n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  read matrices, coefficient and minimizer
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&W,&Wnew,&DeltaK,&DeltaKnew);
    CALL(mess_matrix_read_formated(argv[1], W,          MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[2], Wnew,       MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[3], DeltaK,     MESS_DENSE));
    CALL(mess_matrix_read_formated(argv[4], DeltaKnew,  MESS_DENSE));

    alpha   = atof(argv[ 5]);
    beta    = atof(argv[ 6]);
    delta   = atof(argv[ 7]);
    gamma   = atof(argv[ 8]);
    epsilon = atof(argv[ 9]);
    zeta    = atof(argv[10]);
    lambda  = atof(argv[11]);

    /*-----------------------------------------------------------------------------
     *  compute coefficients and compare the results
     *-----------------------------------------------------------------------------*/
    CALL(mess_lrcfadi_nm_linesearch_poly(W, Wnew, DeltaK, DeltaKnew, &_alpha, &_beta, &_delta, &_gamma, &_epsilon, &_zeta));

#define SCALAR_ERROR_MACRO(NAME,SCALAR, MYSCALAR, ABS_ERR, REL_ERR, TOL, ERR)                                               \
    (ABS_ERR) = fabs((SCALAR)-(MYSCALAR));                                                                                  \
    (REL_ERR) = (ABS_ERR)/fabs(SCALAR);                                                                                     \
    printf(#NAME "\n");                                                                                                     \
    printf("correct=%+.10e\t computed=%+.10e\t abs_err=%+.10e\t rel_err=%+.10e\n",(SCALAR),(MYSCALAR),(ABS_ERR),(REL_ERR)); \
    if((REL_ERR)>(TOL)){printf("failed for: " #NAME"\n"); ++(ERR);}                                                         \
    printf("----------------------\n");

    //error for each scalar of the polynomial
    SCALAR_ERROR_MACRO(alpha,   alpha,      _alpha,     abs_err, rel_err, tol, err);
    SCALAR_ERROR_MACRO(beta,    beta,       _beta,      abs_err, rel_err, tol, err);
    SCALAR_ERROR_MACRO(delta,   delta,      _delta,     abs_err, rel_err, tol, err);
    SCALAR_ERROR_MACRO(gamma,   gamma,      _gamma,     abs_err, rel_err, tol, err);
    SCALAR_ERROR_MACRO(epsilon, epsilon,    _epsilon,   abs_err, rel_err, tol, err);
    SCALAR_ERROR_MACRO(zeta,    zeta,       _zeta,      abs_err, rel_err, tol, err);

    /*-----------------------------------------------------------------------------
     *  compute minimizer with correct (given) coefficients
     *-----------------------------------------------------------------------------*/
    CALL(mess_lrcfadi_nm_linesearch_lambda(&alpha, &beta, &delta, &gamma, &epsilon, &zeta, &_lambda));

    SCALAR_ERROR_MACRO(lambda, lambda, _lambda, abs_err, rel_err, tol, err);

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&W,&Wnew,&DeltaK,&DeltaKnew);

    return 0;

}


