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
 * @file matlab/c_interface/tests/enum_conversions.c
 * @brief
 * @author @mbehr
 *
 */


#include "interface_matlab.h"


#undef MEXMESS_GENERATOR
#define MEXMESS_GENERATOR(CMESS_ENUM)                                                                       \
    void mex_test_##CMESS_ENUM( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){                \
                                                                                                            \
        init_mexmess();                                                                                     \
                                                                                                            \
         /* check input/output arguments */                                                                 \
        if (nrhs != 1 ) {                                                                                   \
            csc_error_message("The function needs one input argument but %d are given.\n", (int) nrhs);     \
        }                                                                                                   \
                                                                                                            \
         /* check if class is correct */                                                                    \
        if (!mxIsClass(prhs[0], #CMESS_ENUM )){                                                             \
            csc_error_message("Expected "#CMESS_ENUM" class.\n");                                           \
            return ;                                                                                        \
        }                                                                                                   \
                                                                                                            \
        /* call conversion function from MEX-M.E.S.S. to C-M.E.S.S. */                                      \
        CMESS_ENUM c_enum = CMESS_ENUM##_from_mexmess(prhs[0]);                                             \
                                                                                                            \
        /* call conversion function from C-M.E.S.S. to MEX-M.E.S.S. */                                      \
        plhs[0] = CMESS_ENUM##_to_mexmess(c_enum);                                                          \
                                                                                                            \
        return;                                                                                             \
    }

MEXMESS_GENERATOR(mess_direct_cholpackage_t);
MEXMESS_GENERATOR(mess_direct_lupackage_t);
MEXMESS_GENERATOR(mess_equation_t);
MEXMESS_GENERATOR(mess_memusage_t);
MEXMESS_GENERATOR(mess_multidirect_t);
MEXMESS_GENERATOR(mess_operation_t);
MEXMESS_GENERATOR(mess_parameter_t);
MEXMESS_GENERATOR(mess_residual_t);
MEXMESS_GENERATOR(mess_norm_t);






