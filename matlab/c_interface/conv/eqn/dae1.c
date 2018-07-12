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
 * @file matlab/c_interface/conv/eqn/dae1.c
 * @brief Equation conversion for dae1 function handles.
 * @author @mbehr
 *
 */


#include "interface_matlab.h"


mess_equation eqn_conv_lyap_dae1(const mxArray* m_eqn, const mxArray* m_opt, mess_freelist mem) {
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input pointer m_eqn and class
     *-----------------------------------------------------------------------------*/
    if(!m_eqn){
        csc_error_message("m_eqn points to NULL.\n");
        return NULL;
    }

    if (!mxIsClass(m_eqn, "mess_equation_glyap_dae1")){
        csc_error_message("Expected mess_equation_glyap_dae1 class.\nAn error occured return NULL Pointer.");
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
    mess_matrix E11     = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"E11"));
    mess_matrix A11     = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"A11"));
    mess_matrix A12     = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"A12"));
    mess_matrix A21     = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"A21"));
    mess_matrix A22     = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"A22"));
    mess_matrix B       = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"B"));

    mess_freelist_add_mess_matrix(mem, E11);
    mess_freelist_add_mess_matrix(mem, A11);
    mess_freelist_add_mess_matrix(mem, A12);
    mess_freelist_add_mess_matrix(mem, A21);
    mess_freelist_add_mess_matrix(mem, A22);
    mess_freelist_add_mess_matrix(mem, B);

    /*-----------------------------------------------------------------------------
     *  get options
     *-----------------------------------------------------------------------------*/
    mess_options opt = mess_options_from_mexmess(m_opt);

    /*-----------------------------------------------------------------------------
     *  build equation
     *-----------------------------------------------------------------------------*/
    mess_equation eqn;
    mess_equation_init(&eqn);
    ret = mess_equation_glyap_dae1(eqn, opt, E11, A11, A12, A21, A22, B);

    /*-----------------------------------------------------------------------------
     *  add equation we need to freelist
     *-----------------------------------------------------------------------------*/
    mess_freelist_add_mess_equation(mem,eqn);

    /*-----------------------------------------------------------------------------
     *  clear options structure
     *-----------------------------------------------------------------------------*/
    mess_options_clear(&opt);


    return eqn;
}



mess_equation eqn_conv_ric_dae1(const mxArray* m_eqn, const mxArray* m_opt,  mess_freelist mem) {
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input pointer m_eqn and class
     *-----------------------------------------------------------------------------*/
    if(!m_eqn){
        csc_error_message("m_eqn points to NULL.\n");
        return NULL;
    }

    if (!mxIsClass(m_eqn, "mess_equation_griccati_dae1")){
        csc_error_message("Expected mess_equation_griccati_dae1 class.\nAn error occured return NULL Pointer.");
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
    mess_matrix E11     = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"E11"));
    mess_matrix A11     = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"A11"));
    mess_matrix A12     = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"A12"));
    mess_matrix A21     = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"A21"));
    mess_matrix A22     = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"A22"));
    mess_matrix B       = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"B"));
    mess_matrix C       = mess_matrix_from_mexmess(mxGetProperty(m_eqn,0,"C"));

    /*-----------------------------------------------------------------------------
     *  get options
     *-----------------------------------------------------------------------------*/
    mess_options opt = mess_options_from_mexmess(m_opt);

    /*-----------------------------------------------------------------------------
     *  build equation
     *-----------------------------------------------------------------------------*/
    mess_equation eqn;
    mess_equation_init(&eqn);
    ret = mess_equation_griccati_dae1(eqn, opt, E11, A11, A12, A21, A22, B, C);

    /*-----------------------------------------------------------------------------
     *  add everything we need to freelist
     *-----------------------------------------------------------------------------*/
    mess_freelist_add_mess_matrix(mem, E11);
    mess_freelist_add_mess_matrix(mem, A11);
    mess_freelist_add_mess_matrix(mem, A12);
    mess_freelist_add_mess_matrix(mem, A21);
    mess_freelist_add_mess_matrix(mem, A22);
    mess_freelist_add_mess_matrix(mem, B);
    mess_freelist_add_mess_matrix(mem, C);
    mess_freelist_add_mess_equation(mem,eqn);

    /*-----------------------------------------------------------------------------
     *  clear options structure
     *-----------------------------------------------------------------------------*/
    mess_options_clear(&opt);


    return eqn;
}


