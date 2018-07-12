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
 * @file matlab/c_interface/conv/eqn/dae2.c
 * @brief Equation conversion for so1 function handles.
 * @author @mbehr
 *
 */


#include "interface_matlab.h"



mess_equation eqn_conv_lyap_dae2(const mxArray* m_eqn, const mxArray* m_opt, mess_freelist mem) {
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input pointer m_eqn and class
     *-----------------------------------------------------------------------------*/
    if(!m_eqn){
        csc_error_message("m_eqn points to NULL.\n");
        return NULL;
    }

    if (!mxIsClass(m_eqn, "mess_equation_glyap_dae2")){
        csc_error_message("Expected mess_equation_glyap_dae2 class.\nAn error occured return NULL Pointer.");
        return NULL;
    }

    /*-----------------------------------------------------------------------------
     *  check input pointer m_opt and class
     *------------------------------------------------------------------------------*/
    if(!m_opt){
        csc_error_message("m_opt points to NULL.\nAn error occured return NULL Pointer.");
        return NULL;
    }

    if (!mxIsClass(m_opt, "mess_options")){
        csc_error_message("Expected mess_options class.\nAn error occured return NULL Pointer.");
        return NULL;
    }

    /*-----------------------------------------------------------------------------
     * check mem pointer
     *-----------------------------------------------------------------------------*/
    if(!mem){
        csc_error_message("mem points to NULL.\n");
        return NULL;
    }

    /*-----------------------------------------------------------------------------
     *  get matrices
     *-----------------------------------------------------------------------------*/
    mess_matrix M       = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"M"));
    mess_matrix A       = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"A"));
    mess_matrix G       = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"G"));
    mess_matrix B       = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"B"));
    mess_matrix K0      = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"K0"));
    double delta        = mxGetScalar(mxGetProperty(m_eqn,0,"delta"));

    /*-----------------------------------------------------------------------------
     *  get options
     *-----------------------------------------------------------------------------*/
    mess_options opt = mess_options_from_mexmess(m_opt);

    /*-----------------------------------------------------------------------------
     *  build equation
     *-----------------------------------------------------------------------------*/
    mess_equation eqn;
    mess_equation_init(&eqn);

    ret = mess_equation_glyap_dae2(eqn, opt, M, A, G, B, delta);

    /*-----------------------------------------------------------------------------
     *  stabilize system if K0 is given
     *-----------------------------------------------------------------------------*/
    mess_equation eqn_stable;
    if(K0){
        MSG_PRINT("Stabilize system with given initial feedback.");
        mess_equation_init(&eqn_stable);

        if(opt->type ==MESS_OP_NONE){
            ret = mess_equation_stable(eqn_stable,opt,eqn,eqn->B,K0);
        }else{
            ret = mess_equation_stable(eqn_stable,opt,eqn,K0,eqn->B);
        }
    }

    mess_options_clear(&opt);

    /*-----------------------------------------------------------------------------
     *  add everything we need to freelist
     *-----------------------------------------------------------------------------*/
    mess_freelist_add_mess_matrix(mem, M);
    mess_freelist_add_mess_matrix(mem, A);
    mess_freelist_add_mess_matrix(mem, G);
    mess_freelist_add_mess_matrix(mem, B);

    if(K0){
        //it is important that eqn_stable is inserted before eqn, because of the cleaning method afterwards
        mess_freelist_add_mess_equation(mem, eqn_stable);
        mess_freelist_add_mess_equation(mem, eqn);
        mess_freelist_add_mess_matrix(mem, K0);
        return  eqn_stable;
    }else{
        mess_freelist_add_mess_equation(mem, eqn);
        return eqn;
    }

}

mess_equation eqn_conv_ric_dae2(const mxArray* m_eqn, const mxArray* m_opt,  mess_freelist mem) {
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input pointer m_eqn and class
     *-----------------------------------------------------------------------------*/
    if(!m_eqn){
        csc_error_message("m_eqn points to NULL.\n");
        return NULL;
    }

    if (!mxIsClass(m_eqn, "mess_equation_griccati_dae2")){
        csc_error_message("Expected mess_equation_griccati_dae2 class.\nAn error occured return NULL Pointer.");
        return NULL;
    }

    /*-----------------------------------------------------------------------------
     *  check input pointer m_opt and class
     *------------------------------------------------------------------------------*/
    if(!m_opt){
        csc_error_message("m_opt points to NULL.\nAn error occured return NULL Pointer.");
        return NULL;
    }

    if (!mxIsClass(m_opt, "mess_options")){
        csc_error_message("Expected mess_options class.\nAn error occured return NULL Pointer.");
        return NULL;
    }

    /*-----------------------------------------------------------------------------
     * check mem pointer
     *-----------------------------------------------------------------------------*/
    if(!mem){
        csc_error_message("mem points to NULL.\n");
        return NULL;
    }

    /*-----------------------------------------------------------------------------
     *  get matrices
     *-----------------------------------------------------------------------------*/
    mess_matrix M       = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"M"));
    mess_matrix A       = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"A"));
    mess_matrix G       = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"G"));
    mess_matrix B       = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"B"));
    mess_matrix C       = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"C"));
    double delta        = mxGetScalar(mxGetProperty(m_eqn,0,"delta"));

    /*-----------------------------------------------------------------------------
     *  get options
     *-----------------------------------------------------------------------------*/
    mess_options opt = mess_options_from_mexmess(m_opt);

    /*-----------------------------------------------------------------------------
     *  build equation
     *-----------------------------------------------------------------------------*/
    mess_equation eqn;
    mess_equation_init(&eqn);
    ret = mess_equation_griccati_dae2(eqn, opt, M, A, G, B, C, delta);

    /*-----------------------------------------------------------------------------
     *  add everything we need to freelist
     *-----------------------------------------------------------------------------*/
    mess_freelist_add_mess_matrix(mem, M);
    mess_freelist_add_mess_matrix(mem, A);
    mess_freelist_add_mess_matrix(mem, G);
    mess_freelist_add_mess_matrix(mem, B);
    mess_freelist_add_mess_matrix(mem, C);
    mess_freelist_add_mess_equation(mem,eqn);

    /*-----------------------------------------------------------------------------
     *  clear options structure
     *-----------------------------------------------------------------------------*/
    mess_options_clear(&opt);


    return eqn;
}

