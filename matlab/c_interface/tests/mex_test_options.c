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
 * @file matlab/c_interface/tests/mex_test_options.c
 * @brief
 * @author @mbehr
 */

#include "interface_matlab.h"


void mex_test_options( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    init_mexmess();

    /*-----------------------------------------------------------------------------
     *  check input/output arguments
     *-----------------------------------------------------------------------------*/
    if (nrhs != 1 ) {
        csc_error_message("The function needs one input argument but %d are given.\n", (int) nrhs);
    }

    if (nlhs != 1){
        csc_error_message("The function needs one out argument but %d are given.\n", (int) nlhs);
    }

    /*-----------------------------------------------------------------------------
     *  check if class is correct
     *-----------------------------------------------------------------------------*/
    if (!mxIsClass(prhs[0], "mess_options")){
        csc_error_message("Expected mess_options class.\n");
        return ;
    }

    /*-----------------------------------------------------------------------------
     *  get options from matlab
     *-----------------------------------------------------------------------------*/
    mess_options opt =  mess_options_from_mexmess (prhs[0]);

    /*-----------------------------------------------------------------------------
     *  convert options to matlab
     *-----------------------------------------------------------------------------*/
    plhs[0] =  mess_options_to_mexmess(opt);

    /*-----------------------------------------------------------------------------
     *  clear options structure
     *-----------------------------------------------------------------------------*/
    mess_options_clear(&opt);


    return;
}


