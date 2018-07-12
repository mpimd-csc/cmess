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
 * @file matlab/c_interface/conv/mess_enum_converter.c
 * @brief
 * @author @mbehr
 */


#include "interface_matlab.h"

#define MESS_ENUM_CONVERTER(CMESS_ENUM_NAME, DEFAULT_ENUM)                                              \
    CMESS_ENUM_NAME CMESS_ENUM_NAME## _from_mexmess (const mxArray *instance){                          \
                                                                                                        \
        mxArray * temp = NULL;                                                                          \
                                                                                                        \
        /* check input pointer*/                                                                        \
        if(!instance){                                                                                  \
            csc_error_message("NULL Pointer Exception.\n");                                             \
            csc_error_message("An error occured return " #DEFAULT_ENUM "\n");                           \
            return DEFAULT_ENUM;                                                                        \
        }                                                                                               \
                                                                                                        \
        /*  check input class */                                                                        \
        if (!mxIsClass(instance, #CMESS_ENUM_NAME)){                                                    \
            csc_error_message("Expected " # CMESS_ENUM_NAME " class.\n");                               \
            csc_error_message("An error occured return "#DEFAULT_ENUM"\n");                             \
            return DEFAULT_ENUM;                                                                        \
        }                                                                                               \
                                                                                                        \
        /* convert to mess_parameter_t */                                                               \
        mxArray* myinstance = (mxArray*) instance;                                                      \
        mexCallMATLAB(1,&temp,1,&myinstance,"int8");                                                    \
        if(temp){                                                                                       \
            int8_t * val = (int8_t*) mxGetData(temp);                                                   \
            CMESS_ENUM_NAME para = *val;                                                                \
            return para;                                                                                \
        }                                                                                               \
                                                                                                        \
        csc_error_message("mexCallMATLAB int8 failed.\n");                                              \
        csc_error_message("An error occured return "#DEFAULT_ENUM".\n");                                \
                                                                                                        \
        return DEFAULT_ENUM;                                                                            \
    }                                                                                                   \
                                                                                                        \
    mxArray* CMESS_ENUM_NAME## _to_mexmess (CMESS_ENUM_NAME c_enum){                                    \
                                                                                                        \
        mxArray* temp = NULL;                                                                           \
        mxArray* m_enum = NULL;                                                                         \
                                                                                                        \
        /* create double scalar and parameter type */                                                   \
        temp = mxCreateDoubleScalar((int)c_enum);                                                       \
        if(temp==NULL){                                                                                 \
            csc_error_message("mxCreateDoubleScalar failed.\n");                                        \
            return NULL;                                                                                \
        }                                                                                               \
                                                                                                        \
        mexCallMATLAB(1,&m_enum,1,&temp,#CMESS_ENUM_NAME);                                              \
                                                                                                        \
        /*  check result and return value */                                                            \
        if(!m_enum){                                                                                    \
            csc_error_message("mexCallMATLAB "#CMESS_ENUM_NAME" failed\n");                             \
            return NULL;                                                                                \
        }                                                                                               \
        return m_enum;                                                                                  \
    }


MESS_ENUM_CONVERTER(mess_direct_cholpackage_t,  MESS_DIRECT_DEFAULT_CHOLESKY);
MESS_ENUM_CONVERTER(mess_direct_lupackage_t,    MESS_DIRECT_DEFAULT_LU);
MESS_ENUM_CONVERTER(mess_equation_t,            MESS_EQN_NONE);
MESS_ENUM_CONVERTER(mess_memusage_t,            MESS_MEMORY_LOW);
MESS_ENUM_CONVERTER(mess_multidirect_t,         MESS_MULTIDIRECT_SPARSE_LU);
MESS_ENUM_CONVERTER(mess_operation_t,           MESS_OP_NONE);
MESS_ENUM_CONVERTER(mess_parameter_t,           MESS_LRCFADI_PARA_MINMAX);
MESS_ENUM_CONVERTER(mess_residual_t,            MESS_RESIDUAL_INDEFINITE);
MESS_ENUM_CONVERTER(mess_norm_t,                MESS_2_NORM);





