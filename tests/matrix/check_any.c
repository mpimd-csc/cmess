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
 * @file tests/matrix/check_any.c
 * @brief Check function @ref mess_matrix_any.
 * @author @mbehr
 * @test
 * This function checks if @ref mess_matrix_any  works correctly
 * @}
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <math.h>

#include "../call_macro.h"
#include "mess/mess.h"



static mess_int_t f_real(double val){ return fabs(val)<0.1; }
static mess_int_t f_cpx(mess_double_cpx_t val){ return cabs(val)<0.1; }


#define CHECK_ANY(TOL,ROWS,COLS,DTA,STA,DIM){                               \
    int ret;                                                                \
    double diff;                                                            \
    mess_matrix A, TMPA;                                                    \
    MESS_INIT_MATRICES(&A,&TMPA);                                           \
    \
    mess_vector v, tmpv;                                                    \
    MESS_INIT_VECTORS(&v,&tmpv);                                            \
    CALL(mess_vector_alloc(v,0,MESS_REAL));                                 \
    CALL(mess_vector_alloc(tmpv,0,MESS_REAL));                              \
    \
    /*generate random matrix and create reference solution*/                \
    CALL(mess_matrix_rand(A,ROWS,COLS,MESS_DENSE,DTA,0.5));                 \
    CALL(mess_matrix_convert(A,TMPA,STA));                                  \
    CALL(mess_matrix_any(A,f_real,f_cpx,DIM,v));                            \
    \
    /* compute solution and compare */                                      \
    CALL(mess_matrix_any(TMPA,f_real,f_cpx,DIM,tmpv));                      \
    CALL(mess_vector_diffnorminf(v,tmpv,&diff));                            \
    \
    if(diff>TOL){                                                           \
        printf("FAILED\n");                                                 \
        printf("Matrix\n");                                                 \
        mess_matrix_printinfo(TMPA);                                        \
        return 1;                                                           \
    }                                                                       \
    MESS_CLEAR_MATRICES(&A,&TMPA);                                          \
    MESS_CLEAR_VECTORS(&v,&tmpv);                                           \
}



int main (int argc, char *argv[]){
    mess_init();
    mess_error_level=0;
    mess_int_t rows, cols;
    double tol = mess_eps();

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
    mess_int_t idtA;
    mess_datatype_t dts [] = {MESS_REAL, MESS_COMPLEX};
    mess_int_t istA;
    mess_storage_t sts [] = {MESS_DENSE, MESS_CSR, MESS_CSC, MESS_COORD};
    mess_int_t dim;
    for(rows=1;rows<100;rows=rows+5){
        for(cols=1;cols<100;cols*=2){
            for(istA=0;istA<4;istA++){
                for(idtA=0;idtA<2;idtA++){
                    for(dim=0;dim<2;dim++){
                        CHECK_ANY(tol,rows,cols,dts[idtA],sts[istA],dim)
                    }
                }
            }
        }
    }
    return 0;
}

