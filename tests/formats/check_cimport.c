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
 * @file tests/formats/check_cimport.c
 * @brief Check the creation of different storage type matrices from given arrays.
 * @author @mbehr
 * @test
 * This function checks the @ref mess_matrix_csr, @ref mess_matrix_csc and @ref mess_matrix_coord functions defined in cimport.c
 * that means it checks if a @ref MESS_CSR, @ref MESS_CSC and @ref MESS_COORD matrix is created correctly.
 *
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mess/mess.h"
#include "../call_macro.h"



int main (int argc, char *argv[]){
    mess_init();
    mess_error_level=2;
    mess_int_t  rows, cols;
    int ret;
    double tol = sqrt(mess_eps()), diff;
    mess_matrix A, Atmp;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc != 1) {
        printf("usage: %s \n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  init matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&Atmp);


    /*-----------------------------------------------------------------------------
     *  check cases
     *-----------------------------------------------------------------------------*/
    mess_int_t idt;
    mess_datatype_t dts [] = {MESS_REAL, MESS_COMPLEX};
    mess_int_t ist;
    mess_storage_t sts [] = {MESS_DENSE, MESS_CSR, MESS_CSC, MESS_COORD};
    double ps [] = {0.1, 0.2, 0.5, 0.8, 1.0};
    mess_int_t ip;

    for(ip=0;ip<5;++ip){
        for(idt=0;idt<2;idt++){
            for(rows=1;rows<40;rows=rows+2){
                for(cols=1;cols<40;cols*=2){
                    for(ist=0;ist<4;ist++){

                        /*-----------------------------------------------------------------------------
                         *  - create random matrix
                         *  - call mess_matrix_csr, mess_matrix_csc, mess_matrix_coord
                         *  - check the results
                         *-----------------------------------------------------------------------------*/
                        mess_matrix_rand(A, rows, cols, sts[ist], dts[idt], ps[ip]);

                        if(MESS_IS_CSR(A)){
                            CALL(mess_matrix_csr(Atmp, A->rows, A->cols, A->rowptr, A->colptr, A->values, A->values_cpx));
                            CALL(mess_matrix_diffnormf(A,Atmp,&diff));
                            if(diff>tol || A->data_type != Atmp->data_type || A->store_type != Atmp->store_type){
                                printf("FAILED:\n");
                                printf("Matrix A:\t"); mess_matrix_printinfo(A);
                                MESS_CLEAR_MATRICES(&A,&Atmp);
                                return 1;
                            }

                        }else if (MESS_IS_CSC(A)){
                            CALL(mess_matrix_csc(Atmp, A->rows, A->cols, A->rowptr, A->colptr, A->values, A->values_cpx));
                            CALL(mess_matrix_diffnormf(A,Atmp,&diff));
                            if(diff>tol || A->data_type != Atmp->data_type || A->store_type != Atmp->store_type){
                                printf("FAILED:\n");
                                printf("Matrix A:\t"); mess_matrix_printinfo(A);
                                MESS_CLEAR_MATRICES(&A,&Atmp);
                                return 1;
                            }

                        }else if (MESS_IS_COORD(A)){
                            CALL(mess_matrix_coord(Atmp, A->rows, A->cols, A->nnz, 0,  A->rowptr, A->colptr, A->values, A->values_cpx));
                            CALL(mess_matrix_diffnormf(A,Atmp,&diff));
                            if(diff>tol || A->data_type != Atmp->data_type || A->store_type != Atmp->store_type){
                                printf("FAILED:\n");
                                printf("Matrix A:\t"); mess_matrix_printinfo(A);
                                MESS_CLEAR_MATRICES(&A,&Atmp);
                                return 1;
                            }
                        }
                        CALL(mess_matrix_reset(Atmp));
                    }
                }
            }
        }
    }


    /*-----------------------------------------------------------------------------
     *  clear and return
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&Atmp);

    return 0;
}
