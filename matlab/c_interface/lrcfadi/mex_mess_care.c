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
 * @file matlab/c_interface/lrcfadi/mex_mess_care.c
 * @brief
 * @author @mbehr
 */


#include "interface_matlab.h"


void mex_mess_care( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    init_mexmess();

    mess_matrix A = NULL;
    mess_matrix E = NULL;
    mess_matrix B = NULL;
    mess_matrix C = NULL;
    mess_matrix Z = NULL;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if (nrhs != 4 && nrhs!=3) {
        csc_error_message("The function needs three or four input arguments but %d are given.\n", (int) nrhs);
    }

    A = mess_matrix_from_mexmess(prhs[0]);
    if ( A == NULL ) {
        csc_error_message("The first argument is invalid.\n");
        return;
    }
    B = mess_matrix_from_mexmess(prhs[(nrhs==4)?2:1]);
    if ( B == NULL ) {
        MESS_CLEAR_MATRICES(&A);
        if(nrhs==4){
            csc_error_message("The third argument is invalid.\n");
        }else{
            csc_error_message("The second argument is invalid.\n");
        }
        return;
    }
    C = mess_matrix_from_mexmess(prhs[(nrhs==4)?3:2]);
    if ( C == NULL ) {
        MESS_CLEAR_MATRICES(&A, &B);
        if(nrhs==4){
            csc_error_message("The fourth argument is invalid.\n");
        }else{
            csc_error_message("The third argument is invalid.\n");
        }
        return;
    }


    /*-----------------------------------------------------------------------------
     *  Check Input
     *-----------------------------------------------------------------------------*/
    if ( A->rows != B->rows ) {
        csc_error_message("The number of rows of A ( %d, ...) and B( %d, ...) does not match.\n", A->rows, B->rows);
        MESS_CLEAR_MATRICES(&A, &B, &C);
        return;
    }
    if ( A->rows != C->cols ) {
        csc_error_message("The number of rows of A ( %d, ...) and the number of colums of C( ..., %d) does not match.\n", A->rows, C->cols);
        MESS_CLEAR_MATRICES(&A, &B, &C);
        return;
    }

    if ( !MESS_IS_REAL(A)) {
        csc_error_message("Matrix A is not real\n");
        MESS_CLEAR_MATRICES(&A, &B, &C);
        return;
    }
    if ( !MESS_IS_REAL(B)) {
        csc_error_message("Matrix B is not real\n");
        MESS_CLEAR_MATRICES(&A, &B, &C);
        return;
    }
    if ( !MESS_IS_REAL(C)) {
        csc_error_message("Matrix C is not real\n");
        MESS_CLEAR_MATRICES(&A, &B, &C);
        return;
    }


    /*-----------------------------------------------------------------------------
     *  Pick the optional E
     *-----------------------------------------------------------------------------*/
    if(nrhs==4){
        E = mess_matrix_from_mexmess(prhs[1]);
        if ( E != NULL ) {
            if ( E->rows != A->rows || E->cols != A->cols) {
                csc_error_message("Dimension of A ( %d, %d ) and E ( %d, %d ) does not match.\n", A->rows, A->cols, E->rows, E->rows);
                MESS_CLEAR_MATRICES(&A, &E, &B, &C);
                return;
            }
            if ( !MESS_IS_REAL(E)) {
                csc_error_message("Matrix E is not real\n");
                MESS_CLEAR_MATRICES(&A, &E, &B, &C);
                return;
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  Solve
     *-----------------------------------------------------------------------------*/
    mess_matrix_init(&Z);

    if(nrhs==4){
        ret = mess_care(A, E, B, C, Z);
    }else{
        ret = mess_care(A, NULL, B, C, Z);
    }

    mess_matrix_clear(&A);
    mess_matrix_clear(&B);
    mess_matrix_clear(&C);
    if (E) mess_matrix_clear(&E);

    if ( ret != 0 ) {
        mess_matrix_clear(&Z);
        csc_error_message("Failed to solve the Lyapunov Equation using the high-level interface.\n-> %d - %s \n", (int) ret, mess_get_error(ret));
        return;
    }

    plhs[0] = mess_matrix_to_mexmess(Z);
    mess_matrix_clear(&Z);
    return;
}


