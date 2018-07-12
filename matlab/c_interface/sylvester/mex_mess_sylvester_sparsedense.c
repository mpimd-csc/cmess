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


void mex_mess_sylvester_sparsedense( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    init_mexmess();

    mess_matrix A = NULL;
    mess_matrix F = NULL;
    mess_matrix E = NULL;
    mess_matrix H = NULL;
    mess_matrix M = NULL;
    mess_matrix X = NULL;
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
        if(!A){
            csc_error_message("The first argument is invalid.\n");
            return;
        }

        H = mess_matrix_from_mexmess(prhs[1]);
        if(!H){
            csc_error_message("The second argument is invalid.\n");
            MESS_CLEAR_MATRICES(&A);
            return;
        }

        M = mess_matrix_from_mexmess(prhs[2]);
        if(!M){
            csc_error_message("The third argument is invalid.\n");
            MESS_CLEAR_MATRICES(&A,&H);
            return;
        }
    }else if(nrhs==4){
        A = mess_matrix_from_mexmess(prhs[0]);
        if(!A){
            csc_error_message("The first argument is invalid.\n");
            return;
        }

        E = mess_matrix_from_mexmess(prhs[1]);
        if(!E){
            csc_error_message("The second argument is invalid.\n");
            MESS_CLEAR_MATRICES(&A);
            return;
        }

        H = mess_matrix_from_mexmess(prhs[2]);
        if(!H){
            csc_error_message("The third argument is invalid.\n");
            MESS_CLEAR_MATRICES(&A,&E);
            return;
        }

        M = mess_matrix_from_mexmess(prhs[3]);
        if(!M){
            csc_error_message("The fourth argument is invalid.\n");
            MESS_CLEAR_MATRICES(&A,&E,&H);
            return;
        }
    }else{
        A = mess_matrix_from_mexmess(prhs[0]);
        if(!A){
            csc_error_message("The first argument is invalid.\n");
            return;
        }

        F = mess_matrix_from_mexmess(prhs[1]);
        if(!F){
            csc_error_message("The second argument is invalid.\n");
            MESS_CLEAR_MATRICES(&A);
            return;
        }

        E = mess_matrix_from_mexmess(prhs[2]);
        if(!E){
            csc_error_message("The third argument is invalid.\n");
            MESS_CLEAR_MATRICES(&A,&F);
            return;
        }

        H = mess_matrix_from_mexmess(prhs[3]);
        if(!H){
            csc_error_message("The fourth argument is invalid.\n");
            MESS_CLEAR_MATRICES(&A,&F,&E);
            return;
        }

        M = mess_matrix_from_mexmess(prhs[4]);
        if(!M){
            csc_error_message("The fifth argument is invalid.\n");
            MESS_CLEAR_MATRICES(&A,&F,&E,&H);
            return;
        }
    }

    /*-----------------------------------------------------------------------------
     *  Solve
     *-----------------------------------------------------------------------------*/

    ret = mess_direct_init(&syl);
    ret = mess_direct_create_sylvester_sparsedense(A,F,E,H,syl);
    if(ret){
        mess_direct_clear(&syl);
        MESS_CLEAR_MATRICES(&A,&F,&E,&H,&M);
        csc_error_message("Failed to create solver for sparse dense Sylvester equation.\n--> %d - %s \n", (int) ret, mess_get_error(ret));
        return;
    }

    mess_matrix_init(&X);

    ret = mess_direct_solvem(MESS_OP_NONE, syl, M, X);
    if(ret){
        mess_direct_clear(&syl);
        MESS_CLEAR_MATRICES(&A,&F,&E,&H,&M,&X);
        csc_error_message("Failed to solve for sparse dense Sylvester equation.\n--> %d - %s \n", (int) ret, mess_get_error(ret));
        return;
    }

    MESS_CLEAR_MATRICES(&A,&H,&M);
    if (E) mess_matrix_clear(&E);
    if (F) mess_matrix_clear(&F);

    plhs[0] = mess_matrix_to_mexmess(X);
    mess_matrix_clear(&X);
    return;
}


