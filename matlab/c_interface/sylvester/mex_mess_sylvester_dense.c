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
 * @file matlab/c_interface/sylvester/mex_mess_sylvester_sparsedense.c
 * @brief
 * @author @mbehr
 */


#include "interface_matlab.h"


void mex_mess_sylvester_dense( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    init_mexmess();

    mess_matrix A = NULL, F = NULL, E = NULL, H = NULL, M = NULL, X = NULL;
    mess_direct syl=NULL;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if (!( 3<=nrhs && nrhs<=5)) {
        csc_error_message("The function needs three, four or five input arguments but %d are given.\n", (int) nrhs);
    }

    if(nrhs==3){
        A = mess_matrix_from_mexmess(prhs[0]);
        H = mess_matrix_from_mexmess(prhs[1]);
        M = mess_matrix_from_mexmess(prhs[2]);
        if(!A || !H || !M){
            csc_error_message("Cannot transfer matrices from MEX-M.E.S.S. to C-M.E.S.S.\n");
            goto CLEAR;
        }
    }else{
        A = mess_matrix_from_mexmess(prhs[0]);
        F = mess_matrix_from_mexmess(prhs[1]);
        E = mess_matrix_from_mexmess(prhs[2]);
        H = mess_matrix_from_mexmess(prhs[3]);
        M = mess_matrix_from_mexmess(prhs[4]);
        if(!A || !F || !E || !H || !M){
            csc_error_message("Cannot transfer matrices from MEX-M.E.S.S. to C-M.E.S.S.\n");
            goto CLEAR;
        }
    }

    /*-----------------------------------------------------------------------------
     *  Solve
     *-----------------------------------------------------------------------------*/
    ret = mess_direct_init(&syl);
    ret = mess_direct_create_sylvester_dense(A,F,E,H,syl);
    if(ret){
        csc_error_message("Failed to create solver for dense Sylvester equation.\n--> %d - %s \n", (int) ret, mess_get_error(ret));
        goto CLEAR;
    }

    mess_matrix_init(&X);
    ret = mess_direct_solvem(MESS_OP_NONE, syl, M, X);
    if(ret){
        csc_error_message("Failed to solve for dense Sylvester equation.\n--> %d - %s \n", (int) ret, mess_get_error(ret));
        goto CLEAR;
    }

    /*-----------------------------------------------------------------------------
     *  transfer result to MEX-M.E.S.S.
     *-----------------------------------------------------------------------------*/
    plhs[0] = mess_matrix_to_mexmess(X);

    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
CLEAR:
    if(A)       mess_matrix_clear(&A);
    if(F)       mess_matrix_clear(&F);
    if(E)       mess_matrix_clear(&E);
    if(H)       mess_matrix_clear(&H);
    if(M)       mess_matrix_clear(&M);
    if(X)       mess_matrix_clear(&X);
    if(syl)     mess_direct_clear(&syl);

    return;
}


