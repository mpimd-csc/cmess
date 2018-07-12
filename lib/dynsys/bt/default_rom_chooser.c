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
 * @file lib/dynsys/bt/default_rom_chooser.c
 * @brief Default adaptiv ROM dimension chooser.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"

#ifdef _OPENMP
#include <omp.h>
#endif



/**
 * @brief Choose dimension of reduced order model adaptively for balanced truncation (BT).
 * @param[in] SIGMA      input vector with the HANKEL singular values
 * @param[in] tol        input pointer to the cut off tolerance
 * @param[in,out] maxr      pointer to the reduced order
 * @param[in] FOM        input order of the original system (not used)
 * @param[in] aux        input auxillary data (not used)
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_bt_chooseorder_default function chooses the dimension of the reduced order
 * adaptively be using the criteria
 * \f[ \Vert H-H_r \Vert_{\infty} < 2 \sum\limits_{i= maxr+1}^{n} \sigma_i . \f]
 * It is the default for the first order balanced truncation algorithms.
 *
 */
int mess_bt_chooseorder_default(mess_vector SIGMA, double *tol, mess_int_t* maxr, mess_int_t FOM, void*aux) {
    MSG_FNAME(__func__);
    mess_int_t k0;
    double sum;

    mess_check_nullpointer(SIGMA);
    if ( mess_error_level >= 3 ) {
        MSG_INFO("Sigma: \n");
        mess_vector_print(SIGMA);
    }

    k0=SIGMA->dim;
    sum = 0 ;
    while ( sum < (*tol)/2 && k0 > 1){
        k0--;
        sum +=SIGMA->values[k0];
    }
    k0++; // adjust C index
    MSG_INFO("k0 = " MESS_PRINTF_INT " \t sum = %lg\n", k0, sum);
    *maxr = MESS_MIN(k0, *maxr);
    return 0;
}

/**
 * @brief Choose dimension of reduced order model adaptively for balanced truncation (BT) with cut of of nearly zero singular values.
 * @param[in] SIGMA      input vector with the HANKEL singular values
 * @param[in] tol        input pointer to the cut off tolerance
 * @param[in,out] maxr      pointer to the reduced order
 * @param[in] FOM        input order of the original system (not used)
 * @param[in] aux        input auxillary data (not used)
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_bt_chooseorder_minreal function chooses the order of the reduced system such that
 * \f[ \frac{\sigma_{maxr}}{\sigma_1} < tol\f]
 * is fulfilled. \n
 * It is the default for second order balanced truncation algorithms.
 *
 */
int mess_bt_chooseorder_minreal ( mess_vector SIGMA, double *tol, mess_int_t* maxr, mess_int_t FOM, void*aux)
{
    MSG_FNAME(__func__);
    mess_int_t k0 = 0;
    mess_check_nullpointer(SIGMA);
    mess_check_nullpointer(tol);
    mess_check_nullpointer(maxr);

    if ( mess_error_level >= 3 ) {
        MSG_INFO("Sigma: \n");
        mess_vector_print(SIGMA);
    }

    for (k0 = 0; k0 < SIGMA->dim; k0++){
        // MSG_INFO("rel sigma: %lg\n", SIGMA->values[k0]/SIGMA->values[0] );
        if ( SIGMA->values[k0]/SIGMA->values[0] <= *tol) break;
    }
    if (k0 == 0) k0 = 1;
    *maxr=MESS_MIN(k0, *maxr);
    MSG_INFO("k0 = " MESS_PRINTF_INT "\n", k0);


    return 0;
}       /* -----  end of function mess_bt_chooseorder_minreal  ----- */

