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
 * @file matlab/c_interface/conv/mess_status.c
 * @brief
 * @author @mbehr
 */


#include "interface_matlab.h"

static void set_status_properties(mxArray* m_stat, mess_status stat){

    mxSetProperty(m_stat,0,"res2_norms",        mess_vector_to_mexmess(     stat->res2_norms        ));
    mxSetProperty(m_stat,0,"rel_changes",       mess_vector_to_mexmess(     stat->rel_changes       ));
    mxSetProperty(m_stat,0,"it",                mxCreateDoubleScalar(       stat->it                ));
    mxSetProperty(m_stat,0,"res2_change",       mxCreateDoubleScalar(       stat->res2_change       ));
    mxSetProperty(m_stat,0,"res2_0",            mxCreateDoubleScalar(       stat->res2_0            ));
    mxSetProperty(m_stat,0,"res2_norm",         mxCreateDoubleScalar(       stat->res2_norm         ));
    mxSetProperty(m_stat,0,"rel_change",        mxCreateDoubleScalar(       stat->rel_change        ));
    mxSetProperty(m_stat,0,"stop_res2",         mxCreateLogicalScalar(      stat->stop_res2         ));
    mxSetProperty(m_stat,0,"stop_res2c",        mxCreateLogicalScalar(      stat->stop_res2c        ));
    mxSetProperty(m_stat,0,"stop_rel",          mxCreateLogicalScalar(      stat->stop_rel          ));
    mxSetProperty(m_stat,0,"stop_user",         mxCreateLogicalScalar(      stat->stop_user         ));
    mxSetProperty(m_stat,0,"time_all",          mxCreateDoubleScalar(       stat->time_all          ));
    mxSetProperty(m_stat,0,"time_adi",          mxCreateDoubleScalar(       stat->time_adi          ));
    mxSetProperty(m_stat,0,"unstable",          mxCreateLogicalScalar(      stat->unstable          ));
    mxSetProperty(m_stat,0,"n_internal_status", mxCreateDoubleScalar(       stat->n_internal_status ));

}


mxArray* mess_status_to_mexmess (mess_status stat){

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    if(!stat){
        csc_error_message("input points to NULL\n");
        return NULL;
    }

    /*-----------------------------------------------------------------------------
     *  create mex status instance and set substructures
     *-----------------------------------------------------------------------------*/
    mxArray* m_stat = NULL;
    mexCallMATLAB(1,&m_stat,0,NULL,"mess_status");

    //set status of root
    set_status_properties(m_stat,stat);

    //get cell array of substructures
    mwSize dims[1]={stat->n_internal_status};
    mxArray* cell_array = mxCreateCellArray(1, dims);

    //set internal_status sub status
    int i;
    for(i=0;i<stat->n_internal_status;++i){
        mxArray* m_sub_stat = NULL;
        m_sub_stat = mess_status_to_mexmess(stat->internal_status[i]);
        mxSetCell(cell_array, i, m_sub_stat);
        mxSetProperty(m_stat,0,"internal_status", cell_array );
    }

    return m_stat;
}











