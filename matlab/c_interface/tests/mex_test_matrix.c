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
 * @file matlab/c_interface/tests/mex_test_matrix.c
 * @brief
 * @author @mbehr
 *
 */


#include "interface_matlab.h"


void mex_test_matrix( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    init_mexmess();


    /*-----------------------------------------------------------------------------
     *  check input/output arguments
     *-----------------------------------------------------------------------------*/
    if (nrhs != 1 ) {
        csc_error_message("The function needs no input argument but %d are given.\n", (int) nrhs);
    }

    if (nlhs != 1){
        csc_error_message("The function needs one out argument but %d are given.\n", (int) nlhs);
    }

    /*-----------------------------------------------------------------------------
     *  get matrix from mexmess
     *-----------------------------------------------------------------------------*/
    mess_matrix A  = mess_matrix_from_mexmess(prhs[0]);

    /*-----------------------------------------------------------------------------
     *  transfer matrix to mexmess
     *-----------------------------------------------------------------------------*/
    plhs[0] = mess_matrix_to_mexmess(A);

    /*-----------------------------------------------------------------------------
     *  clear matrix
     *-----------------------------------------------------------------------------*/
    mess_matrix_clear(&A);

    return;
}


