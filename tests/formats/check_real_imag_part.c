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
 * @addtogroup test_format
 * @{
 * @file tests/formats/check_real_imag_part.c
 * @brief Check the @ref mess_matrix_realpart @ref mess_matrix_imagpart and @ref mess_matrix_complex_from_parts functions.
 * @author  @mbehr
 * @test
 *
 * This function checks the @ref mess_matrix_realpart @ref mess_matrix_imagpart @ref mess_matrix_complex_from_parts function defined in
 * realimag.c
 * that means it uses these functions together to tests against each other.
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


int main ( int argc, char ** argv) {

    mess_init();
    double eps=mess_eps(), diff, p = 0.4;
    int ret=0, err=0;
    mess_matrix mat, realp, imagp, matp;
    mess_storage_t  storage_types [] = {MESS_DENSE, MESS_CSR, MESS_CSC, MESS_COORD};
    mess_datatype_t data_types [] = {MESS_REAL, MESS_COMPLEX};

    mess_int_t st=0, dt=0;
    mess_int_t rows, cols;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if (argc !=3) {
        printf("usage: %s rows cols\n", argv[0]);
        return 1;
    }

    rows = atoi(argv[1]);
    cols = atoi(argv[2]);

    /*-----------------------------------------------------------------------------
     *  generate matrix
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&mat,&realp,&imagp,&matp);
    mess_int_t seed = 1;
    CALL(mess_matrix_rand_init(&seed));
    CALL(mess_matrix_rand(mat, rows, cols,  MESS_CSR, MESS_REAL, p));

    /*-----------------------------------------------------------------------------
     *  check different storage and datatypes
     *-----------------------------------------------------------------------------*/
    for(st=0;st<3;++st){
        for(dt=0;dt<2;++dt){

            /*-----------------------------------------------------------------------------
             *  convert to datatype and storage type
             *-----------------------------------------------------------------------------*/
            if(data_types[dt]==MESS_REAL){
                CALL(mess_matrix_toreal(mat));
            }else{
                CALL(mess_matrix_tocomplex(mat));
            }

            CALL(mess_matrix_convert(mat, matp, storage_types[st]));
            CALL(mess_matrix_copy(matp,mat));

            /*-----------------------------------------------------------------------------
             *  get real part
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_realpart(mat,realp));
            /*-----------------------------------------------------------------------------
             *  get complex part
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_imagpart(mat,imagp));

            /*-----------------------------------------------------------------------------
             *  recover mat from parts
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_complex_from_parts(realp,imagp,matp));

            /*-----------------------------------------------------------------------------
             *  compare results
             *-----------------------------------------------------------------------------*/
            CALL(mess_matrix_diffnormf(mat,matp,&diff));
            if(diff>eps){
                printf("FAILED:\n");
                printf("Matrix:\n");    mess_matrix_printinfo(mat); mess_matrix_print(mat);
                printf("Real part:\n"); mess_matrix_printinfo(realp); mess_matrix_print(realp);
                printf("Imag part:\n"); mess_matrix_printinfo(imagp); mess_matrix_print(imagp);
                printf("Recover from parts:\n"); mess_matrix_printinfo(matp); mess_matrix_print(matp);
                printf("DIFF=%e\n",diff);
                return 1;
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&mat, &realp, &imagp, &matp);

    return (err>0)?(1):(0);
}

