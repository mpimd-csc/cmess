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
 * @file tests/matrix/check_dynormf.c
 * @brief Check the Frobenius-norm computation of a matrix.
 * @test
 * This function checks the @ref mess_matrix_dynormf function defined in dynorm.c that means it checks if the
 * Frobenius-norm of a matrix \f$ G \f$
 * \f[ \Vert G^H G \Vert_F \f]
 * is computed correctly.
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

/**
 * @brief Determine values used in a FDM matrix.
 * @param[in] x      input space variable \f$ x \f$
 * @param[in] y          input space variable \f$ y \f$
 * @return 1 if \f$ x \f$ is between \f$ 0.1 \f$ and \f$ 0.4 \f$, otherwise 0
 *
 * The @ref fdmcol function returns \f$ 1 \f$ if the space variable \f$ x \f$ fullfills
 * \f[ 0.1 < x < 0.4, \f]
 * otherwise \f$ 0 \f$.
 */
double fdmcol(double x, double y) {
    if ( x > 0.1 && x < 0.4 ) return 1;
    return 0;
}

int main ( int argc, char **argv){
    mess_init();
    mess_matrix  G;
    double nrm;
    double ref;
    int ret = 0;
    mess_int_t dim;
    if ( argc != 3 ) {
        fprintf(stderr, "usage: %s dim ref\n", argv[0]);
        return 0;
    }
    dim = atoi(argv[1]);
    ref = atof(argv[2]);
    mess_init();
    mess_error_level = 3;

    CALL( mess_matrix_init(&G));


    /*-----------------------------------------------------------------------------
     *  generate matrix
     *-----------------------------------------------------------------------------*/
    CALL(mess_matgen_fdmcolumn(G,dim, fdmcol));
    CALL(mess_matrix_dynormf(G,&nrm));
    printf("NNormest = %21.17e\n",nrm);
    printf("err = %21.17e\n", fabs(nrm-ref)/fabs(ref));
    mess_matrix_clear(&G);
    if (  fabs(nrm-ref)/fabs(ref) < dim*dim*mess_eps()) {
        return 0;
    }
    return 1;
}

