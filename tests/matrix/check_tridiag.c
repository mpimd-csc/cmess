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
 * @file tests/matrix/check_tridiag.c
 * @brief Check @ref mess_matrix_tridiag function.
 * @author @mbehr
 * @test
 * This function checks if @ref mess_matrix_tridiag works correctly.
 * @}
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <math.h>

#include "../call_macro.h"
#include "mess/mess.h"



int main (int argc, char *argv[]){
    mess_init();
    mess_error_level=0;
    mess_double_cpx_t lower = -2.0 + 1.0*I, diag = -1.0 + 1.0*I, upper =  2.0 + 1.0*I;
    mess_int_t  rows = 0, cols = 0, i=0;
    double tol = sqrt(mess_eps());
    double diff = 0;
    mess_int_t err = 0, ret=0;

    /*-----------------------------------------------------------------------------
     *  check and read input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 3) {
        printf("usage: %s ROWS COLS \n", argv[0]);
        return 1;
    }

    rows = atoi(argv[1]);
    cols = atoi(argv[2]);

    /*-----------------------------------------------------------------------------
     *  create matrices
     *-----------------------------------------------------------------------------*/
    mess_matrix Adense, Acsc, Acsr, Acoord, Atest;
    MESS_INIT_MATRICES(&Adense, &Acsc, &Acsr, &Acoord, &Atest);

    /*-----------------------------------------------------------------------------
     *  test tridiag function for real
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_alloc(Atest, rows, cols, rows*cols, MESS_DENSE, MESS_COMPLEX));
    CALL(mess_matrix_zeros(Atest));

    for(i=0;i<MESS_MIN(rows,cols);++i){
        // create tridiagonl
        CALL(mess_matrix_setelement_complex(Atest, i, i, creal(diag)));
        if(i+1 < Atest->cols){
            CALL(mess_matrix_setelement_complex(Atest, i, i+1, creal(upper)));
        }
        if(i+1 < Atest->rows){
            CALL(mess_matrix_setelement_complex(Atest, i+1, i, creal(lower)));
        }
    }

    CALL(mess_matrix_tridiag(Adense,    rows,   cols,   MESS_DENSE,     MESS_REAL,  lower,  diag,   upper));
    CALL(mess_matrix_tridiag(Acsc,      rows,   cols,   MESS_CSC,       MESS_REAL,  lower,  diag,   upper));
    CALL(mess_matrix_tridiag(Acsr,      rows,   cols,   MESS_CSR,       MESS_REAL,  lower,  diag,   upper));
    CALL(mess_matrix_tridiag(Acoord,    rows,   cols,   MESS_COORD,     MESS_REAL,  lower,  diag,   upper));

    CALL(mess_matrix_diffnormf(Adense, Atest, &diff));      if(diff>tol) {printf("FAILED:\n"); mess_matrix_print(Adense);  mess_matrix_print(Atest);   return 1;}
    CALL(mess_matrix_diffnormf(Acsr, Atest, &diff));        if(diff>tol) {printf("FAILED:\n"); mess_matrix_print(Acsr);    mess_matrix_print(Acsr);    return 1;}
    CALL(mess_matrix_diffnormf(Acsc, Atest, &diff));        if(diff>tol) {printf("FAILED:\n"); mess_matrix_print(Acsc);    mess_matrix_print(Acsc);    return 1;}
    CALL(mess_matrix_diffnormf(Acoord, Atest, &diff));      if(diff>tol) {printf("FAILED:\n"); mess_matrix_print(Acoord);  mess_matrix_print(Acoord);  return 1;}

    /*-----------------------------------------------------------------------------
     *  test tridiag function for complex
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_alloc(Atest, rows, cols, rows*cols, MESS_DENSE, MESS_COMPLEX));
    CALL(mess_matrix_zeros(Atest));

    for(i=0;i<MESS_MIN(rows,cols);++i){
        // create tridiagonl
        CALL(mess_matrix_setelement_complex(Atest, i, i, diag));
        if(i+1 < Atest->cols){
            CALL(mess_matrix_setelement_complex(Atest, i, i+1, upper));
        }
        if(i+1 < Atest->rows){
            CALL(mess_matrix_setelement_complex(Atest, i+1, i, lower));
        }
    }

    CALL(mess_matrix_tridiag(Adense,    rows,   cols,   MESS_DENSE,     MESS_COMPLEX,  lower,  diag,   upper));
    CALL(mess_matrix_tridiag(Acsc,      rows,   cols,   MESS_CSC,       MESS_COMPLEX,  lower,  diag,   upper));
    CALL(mess_matrix_tridiag(Acsr,      rows,   cols,   MESS_CSR,       MESS_COMPLEX,  lower,  diag,   upper));
    CALL(mess_matrix_tridiag(Acoord,    rows,   cols,   MESS_COORD,     MESS_COMPLEX,  lower,  diag,   upper));

    CALL(mess_matrix_diffnormf(Adense, Atest, &diff));      if(diff>tol) {printf("FAILED:\n"); mess_matrix_print(Adense);  mess_matrix_print(Atest);   return 1;}
    CALL(mess_matrix_diffnormf(Acsr, Atest, &diff));        if(diff>tol) {printf("FAILED:\n"); mess_matrix_print(Acsr);    mess_matrix_print(Acsr);    return 1;}
    CALL(mess_matrix_diffnormf(Acsc, Atest, &diff));        if(diff>tol) {printf("FAILED:\n"); mess_matrix_print(Acsc);    mess_matrix_print(Acsc);    return 1;}
    CALL(mess_matrix_diffnormf(Acoord, Atest, &diff));      if(diff>tol) {printf("FAILED:\n"); mess_matrix_print(Acoord);  mess_matrix_print(Acoord);  return 1;}

    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&Adense, &Acsc, &Acsr, &Acoord, &Atest);


    return err;
}

