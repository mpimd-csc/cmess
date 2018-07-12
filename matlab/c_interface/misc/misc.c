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
 * @file matlab/c_interface/misc/misc.c
 * @brief
 * @author @mbehr
 *
 */


#include "interface_matlab.h"


#undef MEXMESS_GENERATOR
#define MEXMESS_GENERATOR(CMESS_CALL, MXCREATECALL)                                             \
    void mex_##CMESS_CALL (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){         \
        init_mexmess();                                                                         \
        *plhs = MXCREATECALL( CMESS_CALL ());                                                   \
    }


MEXMESS_GENERATOR(mess_have_amd,        mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_have_arpack,     mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_have_bzip2,      mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_have_cholmod,    mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_have_colamd,     mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_have_csparse,    mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_have_matio,      mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_have_mess64,     mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_have_mklpardiso, mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_have_openmp,     mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_have_superlu,    mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_have_umfpack,    mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_have_zlib,       mxCreateLogicalScalar);
MEXMESS_GENERATOR(mess_is_debug,        mxCreateLogicalScalar);

MEXMESS_GENERATOR(mess_version_major,   mxCreateDoubleScalar);
MEXMESS_GENERATOR(mess_version_minor,   mxCreateDoubleScalar);
MEXMESS_GENERATOR(mess_version_patch,   mxCreateDoubleScalar);
MEXMESS_GENERATOR(mess_git_id,          mxCreateString);
MEXMESS_GENERATOR(mess_git_branch,      mxCreateString);


void mex_mess_version(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    init_mexmess();
    mess_version();
}

void mex_mess_version_verbose(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    init_mexmess();
    mess_version_verbose();
}


#undef MEXMESS_GENERATOR
#define MEXMESS_GENERATOR(CMESS_DIRECT_SELECT, CMESS_DIRECT_TYPE)                                                           \
    void mex_##CMESS_DIRECT_SELECT(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){                             \
        init_mexmess();                                                                                                     \
                                                                                                                            \
        int ret = 0;                                                                                                        \
        const int nargs = 1;                                                                                                \
                                                                                                                            \
        /* check input/output arguments */                                                                                  \
        if (nrhs != nargs ) {                                                                                               \
            csc_error_message("The function needs %d input argument but %d are given.\n", (int) nargs ,(int) nrhs);         \
            return;                                                                                                         \
        }                                                                                                                   \
                                                                                                                            \
        /* check if class is correct */                                                                                     \
        if (!mxIsClass(prhs[0], #CMESS_DIRECT_TYPE)){                                                                       \
            csc_error_message("Expected "#CMESS_DIRECT_TYPE" class.\n");                                                    \
            return ;                                                                                                        \
        }                                                                                                                   \
                                                                                                                            \
        /* call conversion function from MEX-M.E.S.S. to C-M.E.S.S. */                                                      \
        CMESS_DIRECT_TYPE c_direct = CMESS_DIRECT_TYPE##_from_mexmess(prhs[0]);                                             \
                                                                                                                            \
        /* call direct select function */                                                                                   \
        ret = CMESS_DIRECT_SELECT ( c_direct );                                                                             \
        if(ret){                                                                                                            \
            csc_error_message("C-M.E.S.S. returned an error %s.\n",mess_get_error(ret));                                    \
            return;                                                                                                         \
        }                                                                                                                   \
    return;                                                                                                                 \
    }


MEXMESS_GENERATOR(mess_direct_lu_select, mess_direct_lupackage_t);
MEXMESS_GENERATOR(mess_direct_chol_select, mess_direct_cholpackage_t);
MEXMESS_GENERATOR(mess_multidirect_select, mess_multidirect_t);
