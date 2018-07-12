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
 * @file matlab/c_interface/conv/mess_options.c
 * @brief
 * @author @mbehr
 */


#include "interface_matlab.h"


mess_options mess_options_from_mexmess (const mxArray *instance){

    mess_options opt;

    /*-----------------------------------------------------------------------------
     *  check input pointer
     *-----------------------------------------------------------------------------*/
    if(!instance){
        csc_error_message("NULL Pointer Exception.\n");
        csc_error_message("An error occured return NULL Pointer.\n");
        return NULL;
    }

    /*-----------------------------------------------------------------------------
     *  check input class
     *-----------------------------------------------------------------------------*/
    if (!mxIsClass(instance, "mess_options")){
        csc_error_message("Expected mess_options class.\n");
        csc_error_message("An error occured return NULL Pointer.\n");
        return NULL;
    }

    /*-----------------------------------------------------------------------------
     *  get type property from options
     *-----------------------------------------------------------------------------*/
    mess_options_init(&opt);
    opt->type               = mess_operation_t_from_mexmess(mxGetProperty(instance,0,"type"));
    opt->residual_method    = mess_residual_t_from_mexmess(mxGetProperty(instance,0,"residual_method"));

    /*-----------------------------------------------------------------------------
     *  get options_newton properties from options
     *-----------------------------------------------------------------------------*/
    opt->nm_maxit           = (mess_int_t)  mxGetScalar(mxGetProperty(mxGetProperty(instance,0,"nm"),0,"maxit"));
    opt->nm_res2_tol        = (double)      mxGetScalar(mxGetProperty(mxGetProperty(instance,0,"nm"),0,"res2_tol"));
    opt->nm_gpStep          = (mess_int_t)  mxGetScalar(mxGetProperty(mxGetProperty(instance,0,"nm"),0,"gpStep"));
    opt->nm_output          = (mess_int_t)  mxGetScalar(mxGetProperty(mxGetProperty(instance,0,"nm"),0,"output"));
    opt->nm_singleshifts    = (mess_int_t)  mxGetScalar(mxGetProperty(mxGetProperty(instance,0,"nm"),0,"singleshifts"));
    opt->nm_linesearch      = (mess_int_t)  mxGetScalar(mxGetProperty(mxGetProperty(instance,0,"nm"),0,"linesearch"));
    opt->K0                 =               mess_matrix_from_mexmess(mxGetProperty(mxGetProperty(instance,0,"nm"),0,"K0"));

    /*-----------------------------------------------------------------------------
     * get options_adi properties from options
     *-----------------------------------------------------------------------------*/
    opt->adi_maxit          = (mess_int_t)   mxGetScalar(mxGetProperty(mxGetProperty(instance,0,"adi"),0,"maxit"));
    opt->adi_res2_tol       = (double)       mxGetScalar(mxGetProperty(mxGetProperty(instance,0,"adi"),0,"res2_tol"));
    opt->adi_res2c_tol      = (double)       mxGetScalar(mxGetProperty(mxGetProperty(instance,0,"adi"),0,"res2c_tol"));
    opt->adi_rel_change_tol = (double)       mxGetScalar(mxGetProperty(mxGetProperty(instance,0,"adi"),0,"rel_change_tol"));
    opt->adi_output         = (mess_int_t)   mxGetScalar(mxGetProperty(mxGetProperty(instance,0,"adi"),0,"output"));
    opt->memory_usage       = mess_memusage_t_from_mexmess(mxGetProperty(mxGetProperty(instance,0,"adi"),0,"memory_usage"));


    /*-----------------------------------------------------------------------------
     *  get options_adi_shifts properties from options_adi
     *-----------------------------------------------------------------------------*/
    opt->adi_shifts_p           = mess_vector_from_mexmess(                mxGetProperty(mxGetProperty(mxGetProperty(instance,0,"adi"),0,"shifts"),0,"p"));
    opt->adi_shifts_arp_p       = (mess_int_t)                 mxGetScalar(mxGetProperty(mxGetProperty(mxGetProperty(instance,0,"adi"),0,"shifts"),0,"arp_p"));
    opt->adi_shifts_arp_m       = (mess_int_t)                 mxGetScalar(mxGetProperty(mxGetProperty(mxGetProperty(instance,0,"adi"),0,"shifts"),0,"arp_m"));
    opt->adi_shifts_paratype    = mess_parameter_t_from_mexmess(           mxGetProperty(mxGetProperty(mxGetProperty(instance,0,"adi"),0,"shifts"),0,"paratype"));
    opt->adi_shifts_l0          = (mess_int_t)                 mxGetScalar(mxGetProperty(mxGetProperty(mxGetProperty(instance,0,"adi"),0,"shifts"),0,"l0"));
    opt->adi_shifts_b0          = mess_vector_from_mexmess(                mxGetProperty(mxGetProperty(mxGetProperty(instance,0,"adi"),0,"shifts"),0,"b0"));


    return opt;

}

mxArray* mess_options_to_mexmess (mess_options opt){

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    if(!opt){
        csc_error_message("input points to NULL\n");
        return NULL;
    }


    /*-----------------------------------------------------------------------------
     *  create mex options instance and get substructures
     *-----------------------------------------------------------------------------*/
    mxArray* m_opt = NULL;
    mexCallMATLAB(1,&m_opt,0,NULL,"mess_options");


    // set opt->type
    mxSetProperty(m_opt,0,"type",    mess_operation_t_to_mexmess(opt->type)   );

    // set residual_method
    mxSetProperty(m_opt,0,"residual_method",    mess_residual_t_to_mexmess(opt->residual_method)   );

    // get nm field
    mxArray* m_opt_nm  = mxGetProperty(m_opt,0,"nm");

    // get adi field
    mxArray* m_opt_adi = mxGetProperty(m_opt,0,"adi");

    // get adi shifts field
    mxArray* m_opt_adi_shifts = mxGetProperty(m_opt_adi,0,"shifts");


    /*-----------------------------------------------------------------------------
     *  set properties in mex options_adi_shifts instance
     *-----------------------------------------------------------------------------*/
    mxSetProperty(m_opt_adi_shifts,0,"p",           mess_vector_to_mexmess(      opt->adi_shifts_p           ));
    mxSetProperty(m_opt_adi_shifts,0,"arp_p",       mxCreateDoubleScalar(       opt->adi_shifts_arp_p       ));
    mxSetProperty(m_opt_adi_shifts,0,"arp_m",       mxCreateDoubleScalar(       opt->adi_shifts_arp_m       ));
    mxSetProperty(m_opt_adi_shifts,0,"paratype",    mess_parameter_t_to_mexmess(opt->adi_shifts_paratype    ));
    mxSetProperty(m_opt_adi_shifts,0,"l0",          mxCreateDoubleScalar(       opt->adi_shifts_l0          ));
    mxSetProperty(m_opt_adi_shifts,0,"b0",          mess_vector_to_mexmess(      opt->adi_shifts_b0          ));

    /*-----------------------------------------------------------------------------
     *  set properties in mex options_adi instance
     *-----------------------------------------------------------------------------*/
    mxSetProperty(m_opt_adi,0,"maxit",              mxCreateDoubleScalar(       opt->adi_maxit          ));
    mxSetProperty(m_opt_adi,0,"res2_tol",           mxCreateDoubleScalar(       opt->adi_res2_tol       ));
    mxSetProperty(m_opt_adi,0,"res2c_tol",          mxCreateDoubleScalar(       opt->adi_res2c_tol      ));
    mxSetProperty(m_opt_adi,0,"rel_change_tol",     mxCreateDoubleScalar(       opt->adi_rel_change_tol ));
    mxSetProperty(m_opt_adi,0,"output",             mxCreateLogicalScalar(      opt->adi_output         ));
    mxSetProperty(m_opt_adi,0,"memory_usage",       mess_memusage_t_to_mexmess( opt->memory_usage       ));
    mxSetProperty(m_opt_adi,0,"shifts",             m_opt_adi_shifts                                    );

    /*-----------------------------------------------------------------------------
     *  set properties in mex options_nm instance
     *-----------------------------------------------------------------------------*/
    mxSetProperty(m_opt_nm,0,"maxit",               mxCreateDoubleScalar(       opt->nm_maxit           ));
    mxSetProperty(m_opt_nm,0,"res2_tol",            mxCreateDoubleScalar(       opt->nm_res2_tol        ));
    mxSetProperty(m_opt_nm,0,"gpStep",              mxCreateDoubleScalar(       opt->nm_gpStep          ));
    mxSetProperty(m_opt_nm,0,"output",              mxCreateLogicalScalar(      opt->nm_output          ));
    mxSetProperty(m_opt_nm,0,"singleshifts",        mxCreateLogicalScalar(      opt->nm_singleshifts    ));
    mxSetProperty(m_opt_nm,0,"linesearch",          mxCreateLogicalScalar(      opt->nm_linesearch      ));
    mxSetProperty(m_opt_nm,0,"K0",                  mess_matrix_to_mexmess(      opt->K0                 ));

    return m_opt;
}


