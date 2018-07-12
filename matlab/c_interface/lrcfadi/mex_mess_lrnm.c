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
 * @file matlab/c_interface/lrcfadi/mex_mess_lrnm.c
 * @brief
 * @author @mbehr
 */

#include "interface_matlab.h"

void mex_mess_lrnm( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    init_mexmess();

    int ret = 0;
    mess_matrix Z      = NULL;
    mess_equation eqn  = NULL;
    mess_options opt   = NULL;
    mess_freelist memlist;
    mess_freelist_init(&memlist);

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if (nrhs != 2 ) {
        mess_freelist_clear(&memlist);
        csc_error_message("The function needs two input arguments but %d are given.\n", (int) nrhs);
    }

    if (nlhs != 2 && nlhs != 1 ) {
        mess_freelist_clear(&memlist);
        csc_error_message("The function needs one or two output arguments but %d are given.\n", (int) nlhs);
    }

    if ( ret ) {
        mess_freelist_clear(&memlist);
        csc_error_message("Failed to initialize the equation.");
        return;
    }

    /*-----------------------------------------------------------------------------
     *  get options structure from matlab
     *-----------------------------------------------------------------------------*/
    opt = mess_options_from_mexmess(prhs[1]);
    if ( opt == NULL ) {
        mess_freelist_clear(&memlist);
        csc_error_message("Failed to copy options from MATLAB.\n");
        return;
    }

    /*-----------------------------------------------------------------------------
     *  get equation from matlab
     *-----------------------------------------------------------------------------*/
    //cast as non const to remove compiler warnings
    mxArray* myprhs0 = (mxArray*)prhs[0];
    mxArray* myprhs1 = (mxArray*)prhs[1];
    eqn = mess_equation_from_mexmess(myprhs0, myprhs1, memlist, MESS_EQN_GRICCATI);

    if ( ret ) {
        mess_freelist_clear(&memlist);
        mess_options_clear(&opt);
        csc_error_message("Failed to get the equation from MATLAB.\n");
        return;
    }

    /*-----------------------------------------------------------------------------
     *  init Z and status and solve lyapunov equation
     *-----------------------------------------------------------------------------*/
    mess_status status = NULL;
    mess_matrix_init(&Z);
    mess_status_init(&status);

    ret = mess_lrnm(eqn, opt, status, Z);


    /*-----------------------------------------------------------------------------
     *  convert Z to matlab
     *-----------------------------------------------------------------------------*/
    if (nlhs == 2 ) {
        plhs[0] = mess_matrix_to_mexmess(Z);
        plhs[1] = mess_status_to_mexmess(status);
    }else{
        plhs[0] = mess_matrix_to_mexmess(Z);
    }


    /*-----------------------------------------------------------------------------
     *  clear everything
     *-----------------------------------------------------------------------------*/
    mess_matrix_clear(&Z);
    mess_options_clear(&opt);
    mess_status_clear(&status);
    mess_freelist_clear(&memlist);
    return;
}


