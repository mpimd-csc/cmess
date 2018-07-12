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
 * @file matlab/c_interface/tests/mex_test_status.c
 * @brief
 * @author @mbehr
 */


#include "interface_matlab.h"


void mex_test_status( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    init_mexmess();
    //   MSG_FNAME(__func__);
    int i;

    /*-----------------------------------------------------------------------------
     *  check input/output arguments
     *-----------------------------------------------------------------------------*/
    if (nrhs != 0 ) {
        csc_error_message("The function needs no input argument but %d are given.\n", (int) nrhs);
    }

    if (nlhs != 1){
        csc_error_message("The function needs one out argument but %d are given.\n", (int) nlhs);
    }

    /*-----------------------------------------------------------------------------
     *  create status
     *-----------------------------------------------------------------------------*/
    mess_status stat;
    mess_status_init(&stat);

    /*-----------------------------------------------------------------------------
     *  set values
     *-----------------------------------------------------------------------------*/
    mess_vector_resize(stat->res2_norms,3);
    mess_vector_resize(stat->rel_changes,3);

    stat->res2_norms->values[0]  = 1;
    stat->res2_norms->values[1]  = 2;
    stat->res2_norms->values[2]  = 3;

    stat->rel_changes->values[0] = -1;
    stat->rel_changes->values[1] = -2;
    stat->rel_changes->values[2] = -3;

    stat->it                = 2;
    stat->res2_change       = 2;
    stat->res2_0            = 2;
    stat->res2_norm         = 2;
    stat->rel_change        = 2;
    stat->stop_res2         = 1;
    stat->stop_res2c        = 1;
    stat->stop_rel          = 1;
    stat->stop_user         = 1;
    stat->time_all          = 2;
    stat->time_adi          = 2;
    stat->unstable          = 1;
    stat->n_internal_status = 2;


    /*-----------------------------------------------------------------------------
     *  create sub status
     *-----------------------------------------------------------------------------*/
    stat->internal_status = (mess_status*) __mess_malloc(sizeof(mess_status)*2);
    for ( i = 0; i < 2; i++){
        mess_status_init(&(stat->internal_status[i]));
        stat->internal_status[i]->it                = i;
        stat->internal_status[i]->res2_change       = i;
        stat->internal_status[i]->res2_norm         = i;
        stat->internal_status[i]->res2_0            = i;
        stat->internal_status[i]->rel_change        = i;
        stat->internal_status[i]->stop_res2         = i;
        stat->internal_status[i]->stop_res2c        = i;
        stat->internal_status[i]->stop_rel          = i;
        stat->internal_status[i]->stop_user         = i;
        stat->internal_status[i]->time_all          = i;
        stat->internal_status[i]->time_adi          = i;
        stat->internal_status[i]->unstable          = i;
        stat->internal_status[i]->n_internal_status = 0;
    }

    /*-----------------------------------------------------------------------------
     *  convert status to matlab
     *-----------------------------------------------------------------------------*/
    plhs[0] =  mess_status_to_mexmess(stat);

    /*-----------------------------------------------------------------------------
     *  clear options structure
     *-----------------------------------------------------------------------------*/
    mess_status_clear(&stat);

    return;
}


