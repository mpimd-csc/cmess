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
 * @file matlab/mexmess/mess_call.c
 * @brief
 * @author @mbehr
 *
 */


#include "interface_matlab.h"


typedef void (*mess_call_func)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]);


typedef struct _mess_call {
    char * name;
    mess_call_func func;
} mess_call_t;


/* Register all mess_calls  */
static mess_call_t mess_calls[] = {
    {"mess_care",                               mex_mess_care},
    {"mess_direct_chol_select",                 mex_mess_direct_chol_select},
    {"mess_direct_lu_select",                   mex_mess_direct_lu_select},
    {"mess_git_branch",                         mex_mess_git_branch},
    {"mess_git_id",                             mex_mess_git_id},
    {"mess_glyap",                              mex_mess_glyap},
    {"mess_gstein",                             mex_mess_gstein},
    {"mess_have_amd",                           mex_mess_have_amd},
    {"mess_have_arpack",                        mex_mess_have_arpack},
    {"mess_have_bzip2",                         mex_mess_have_bzip2},
    {"mess_have_cholmod",                       mex_mess_have_cholmod},
    {"mess_have_colamd",                        mex_mess_have_colamd},
    {"mess_have_csparse",                       mex_mess_have_csparse},
    {"mess_have_matio",                         mex_mess_have_matio},
    {"mess_have_mess64",                        mex_mess_have_mess64},
    {"mess_have_mklpardiso",                    mex_mess_have_mklpardiso},
    {"mess_have_openmp",                        mex_mess_have_openmp},
    {"mess_have_superlu",                       mex_mess_have_superlu},
    {"mess_have_umfpack",                       mex_mess_have_umfpack},
    {"mess_have_zlib",                          mex_mess_have_zlib},
    {"mess_is_debug",                           mex_mess_is_debug},
    {"mess_lradi",                              mex_mess_lradi},
    {"mess_dense_nm_gmpare",                    mex_mess_dense_nm_gmpare},
    {"mess_lrnm",                               mex_mess_lrnm},
    {"mess_lyap",                               mex_mess_lyap},
    {"mess_multidirect_select",                 mex_mess_multidirect_select},
    {"mess_sylvester_sparsedense",              mex_mess_sylvester_sparsedense},
    {"mess_sylvester_dense",                    mex_mess_sylvester_dense},
    {"mess_version",                            mex_mess_version},
    {"mess_version_major",                      mex_mess_version_major},
    {"mess_version_minor",                      mex_mess_version_minor},
    {"mess_version_patch",                      mex_mess_version_patch},
    {"mess_version_verbose",                    mex_mess_version_verbose},
    {"mex_test_matrix",                         mex_test_matrix},
    {"mex_test_mess_direct_cholpackage_t",      mex_test_mess_direct_cholpackage_t},
    {"mex_test_mess_direct_lupackage_t",        mex_test_mess_direct_lupackage_t},
    {"mex_test_mess_equation_t",                mex_test_mess_equation_t},
    {"mex_test_mess_memusage_t",                mex_test_mess_memusage_t},
    {"mex_test_mess_multidirect_t",             mex_test_mess_multidirect_t},
    {"mex_test_mess_operation_t",               mex_test_mess_operation_t},
    {"mex_test_mess_parameter_t",               mex_test_mess_parameter_t},
    {"mex_test_mess_residual_t",                mex_test_mess_residual_t},
    {"mex_test_mess_norm_t",                    mex_test_mess_norm_t},
    {"mex_test_options",                        mex_test_options},
    {"mex_test_status",                         mex_test_status},
    {"mex_test_vector",                         mex_test_vector},
    {NULL, NULL}
};


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    init_mexmess();

    char *input = NULL;
    int cnt = 0;

    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if (nrhs == 0 ){
        csc_print_message("Callable functions:\n");
        while(mess_calls[cnt].name != NULL){
            csc_print_message("%s\n",mess_calls[cnt].name);
            cnt++;
        }
        csc_error_message("The function needs at least one input argument to call a function.\n", (int) nrhs);
        return;
    }

    /*-----------------------------------------------------------------------------
     *  check if first argument is a string
     *-----------------------------------------------------------------------------*/
    if ( mxIsChar(prhs[0]) != 1 ){
        csc_error_message("First argument must be a string.\n");
        return;
    }

    if ( mxGetM(prhs[0]) != 1){
        csc_error_message("First argument must be a string.\n");
        return;
    }

    /*-----------------------------------------------------------------------------
     *  get string
     *-----------------------------------------------------------------------------*/
    input = mxArrayToString(prhs[0]);
    if(!input){
        csc_error_message("First argument must be a string.\n");
        return;
    }

    /*-----------------------------------------------------------------------------
     *  compare and call function
     *-----------------------------------------------------------------------------*/
    cnt = 0;
    while(mess_calls[cnt].name != NULL ){
        if(strcmp(mess_calls[cnt].name,input) == 0){
            //shift the function name away
            mess_calls[cnt].func(nlhs,plhs,nrhs-1,++prhs);
            mxFree(input);
            return;
        }
        cnt++;
    }

    /*-----------------------------------------------------------------------------
     *  function not implemented
     *-----------------------------------------------------------------------------*/
    csc_error_message("Requested function %s is not callable.\n",input);

    /*-----------------------------------------------------------------------------
     *  free input_len
     *-----------------------------------------------------------------------------*/
    mxFree(input);

    return;
}


