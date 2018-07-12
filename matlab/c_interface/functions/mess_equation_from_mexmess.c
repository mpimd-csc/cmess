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
 * @file matlab/c_interface/functions/mess_equation_from_mexmess.c
 * @brief
 * @author @mbehr
 */


#include "interface_matlab.h"

/* Typedef for the convert functions  */
typedef mess_equation (*eqn_converter)(const mxArray* , const mxArray* m_opt,  mess_freelist );

/**
 * @internal
 * @brief Structure for all @pymess to @mess converters.
 * Structure for all @pymess to @mess converters.
 *
 * @attention Internal use only.
 *
 **/
typedef struct _equations_t {
    char * name;
    eqn_converter converter;
} equations_t;


// Register all converts
static equations_t eqns[] = {
    {"mess_equation_glyap",                  eqn_conv_lyap},
    {"mess_equation_griccati",               eqn_conv_ric},
    {"mess_equation_glyap_so1",              eqn_conv_lyap_so1},
    {"mess_equation_griccati_so1",           eqn_conv_ric_so1},
    {"mess_equation_glyap_so2",              eqn_conv_lyap_so2},
    {"mess_equation_griccati_so2",           eqn_conv_ric_so2},
    {"mess_equation_glyap_dae2",             eqn_conv_lyap_dae2},
    {"mess_equation_griccati_dae2",          eqn_conv_ric_dae2},
    {"mess_equation_glyap_dae1",             eqn_conv_lyap_dae1},
    {"mess_equation_griccati_dae1",          eqn_conv_ric_dae1},
    {NULL, NULL}
};

mess_equation mess_equation_from_mexmess(mxArray * m_eqn, mxArray* m_opt, mess_freelist  mem, mess_equation_t eqn_type) {

    init_mexmess();
    mess_equation eqn;
    int cnt = 0;

    if ( eqn_type == MESS_EQN_LYAP ) eqn_type = MESS_EQN_GLYAP;
    if ( eqn_type == MESS_EQN_RICCATI ) eqn_type = MESS_EQN_GRICCATI;

    /*-----------------------------------------------------------------------------
     * Check for special equations like Lyapunov, riccati, ... defined in C
     *-----------------------------------------------------------------------------*/
    cnt = 0;

    //MSG_PRINT("\n");
    MSG_PRINT("\n");
    while ( eqns[cnt].name != NULL ) {
        //MSG_PRINT("Check for.... %s\n", eqns[cnt].name);
        MSG_PRINT("Check for.... %s\n", eqns[cnt].name);


        if (mxIsClass(m_eqn, eqns[cnt].name)){
            //MSG_PRINT("OK\n");
            MSG_PRINT("OK\n");
            eqn = eqns[cnt].converter(m_eqn, m_opt, mem);
            return eqn;
        }
        cnt++;
    }


    /*-----------------------------------------------------------------------------
     *  callback equation
     *-----------------------------------------------------------------------------*/
    //check if equation is inherited from equation
    if(mxIsSubClass(m_eqn,"equation")){
        MSG_PRINT("CALLBACK equation assumed.\n");
        mess_equation_init(&eqn);
        mess_options opt = mess_options_from_mexmess(m_opt);
        mess_callback_equation_mexmess(eqn, opt, m_eqn);
        mess_options_clear(&opt);
        return eqn;
    }

    csc_error_message("Cannot deduce type of defined equation");
    return NULL;

}


