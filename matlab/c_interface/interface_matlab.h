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
 * @file matlab/c_interface/interface_matlab.h
 * @brief Interface to @matlab .
 *
 * This header file provides an interface to @matlab to convert matrices easily from @matlab  to @mess
 * and the way back.
 * Include this file after you included the mex.h header.
 */

#ifndef INTERFACE_MATLAB_H_
#define INTERFACE_MATLAB_H_

#ifndef MESS64
#define MESS64
#endif
#ifndef MESS_MATLAB
#define MESS_MATLAB
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mex.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "cscutils/error_message.h"
#include <complex.h>

/** @addtogroup interfaces_matlab
 * @{ */
/**
 * @brief Function failure handler for @matlab.
 *
 * Function failure handler for @matlab
 */
#define     FUNCTION_FAILURE_HANDLE_MEX( err, cond , fun ) if ( (cond)) { char text[1024]; snprintf(text,1023," %s returned with %d - %s\n", #fun, err, mess_get_error(err)); mexErrMsgTxt(text); }

/**
 * @brief Check if the given mxArray is a structure.
 *
 * The macro forces to be X an instance of a @matlab  structure otherwise an error is printed.
 **/
#define mess_check_struct(X,str) if (!mxIsStruct(X)) { csc_error_message(str); return MESS_ERROR_ARGUMENTS;  }

#ifdef __cplusplus
extern "C" {
#endif


    /*-----------------------------------------------------------------------------
     *  some helper functions
     *-----------------------------------------------------------------------------*/
    char* mess_classname_from_mexmess (const mxArray *instance);
    double mess_double_from_mexmess (const mxArray *Amatlab);
    int mxIsSubClass (const mxArray  *instance, char* class);
    mxArray * mxCreateComplexScalar(mess_double_cpx_t s );
    mess_double_cpx_t mess_complex_from_mexmess (const mxArray *Amatlab);
    void init_mexmess();

    /*-----------------------------------------------------------------------------
     *  conversion of C-M.E.S.S. enums
     *-----------------------------------------------------------------------------*/
    mess_direct_cholpackage_t mess_direct_cholpackage_t_from_mexmess(const mxArray *instance);
    mess_direct_lupackage_t mess_direct_lupackage_t_from_mexmess(const mxArray *instance);
    mess_equation_t mess_equation_t_from_mexmess(const mxArray *instance);
    mess_memusage_t mess_memusage_t_from_mexmess(const mxArray *instance);
    mess_multidirect_t mess_multidirect_t_from_mexmess(const mxArray *instance);
    mess_operation_t mess_operation_t_from_mexmess(const mxArray *instance);
    mess_parameter_t mess_parameter_t_from_mexmess(const mxArray *instance);
    mess_residual_t mess_residual_t_from_mexmess(const mxArray *instance);
    mess_norm_t mess_norm_t_from_mexmess(const mxArray *instance);

    mxArray* mess_direct_cholpackage_t_to_mexmess (mess_direct_cholpackage_t chol_t);
    mxArray* mess_direct_lupackage_t_to_mexmess (mess_direct_lupackage_t chol_t);
    mxArray* mess_equation_t_to_mexmess (mess_equation_t eqn_t);
    mxArray* mess_memusage_t_to_mexmess (mess_memusage_t mem);
    mxArray* mess_multidirect_t_to_mexmess (mess_multidirect_t mdirect_t);
    mxArray* mess_operation_t_to_mexmess (mess_operation_t op);
    mxArray* mess_parameter_t_to_mexmess (mess_parameter_t para);
    mxArray* mess_residual_t_to_mexmess (mess_residual_t res_t);
    mxArray* mess_norm_t_to_mexmess (mess_norm_t nrm_t);


    /*-----------------------------------------------------------------------------
     *  conversion of mess_options structure
     *-----------------------------------------------------------------------------*/
    mess_options mess_options_from_mexmess (const mxArray *instance);
    mxArray * mess_options_to_mexmess (mess_options opt);


    /*-----------------------------------------------------------------------------
     *  conversion of mess_status structure
     *-----------------------------------------------------------------------------*/
    mxArray* mess_status_to_mexmess (mess_status stat);

    /*-----------------------------------------------------------------------------
     *  conversion of mess_equation and predefined equations structure
     *-----------------------------------------------------------------------------*/
    int mess_callback_equation_mexmess(mess_equation eqn, mess_options opt , const  mxArray *meqn);
    mess_equation eqn_conv_lyap(const mxArray* m_eqn, const mxArray* m_opt, mess_freelist mem);
    mess_equation eqn_conv_lyap_dae1(const mxArray* m_eqn, const mxArray* m_opt, mess_freelist mem);
    mess_equation eqn_conv_lyap_dae2(const mxArray* m_eqn, const mxArray* m_opt, mess_freelist mem);
    mess_equation eqn_conv_lyap_so1(const mxArray* m_eqn, const mxArray* m_opt, mess_freelist mem);
    mess_equation eqn_conv_lyap_so2(const mxArray* m_eqn, const mxArray* m_opt, mess_freelist mem);

    mess_equation eqn_conv_ric(const mxArray* m_eqn, const mxArray* m_opt, mess_freelist mem);
    mess_equation eqn_conv_ric_dae1(const mxArray* m_eqn, const mxArray* m_opt, mess_freelist mem);
    mess_equation eqn_conv_ric_dae2(const mxArray* m_eqn, const mxArray* m_opt, mess_freelist mem);
    mess_equation eqn_conv_ric_so1(const mxArray* m_eqn, const mxArray* m_opt, mess_freelist mem);
    mess_equation eqn_conv_ric_so2(const mxArray* m_eqn, const mxArray* m_opt, mess_freelist mem);

    mess_equation mess_equation_from_mexmess(mxArray * m_eqn, mxArray * m_opt, mess_freelist  mem, mess_equation_t eqn_type);

    /*-----------------------------------------------------------------------------
     *  conversion of mess_matrix structure
     *-----------------------------------------------------------------------------*/
    mess_matrix mess_matrix_from_mexmess(const mxArray *Amatlab);
    mxArray * mess_matrix_to_mexmess(mess_matrix A);

    /*-----------------------------------------------------------------------------
     *  conversion of mess_vector structure
     *-----------------------------------------------------------------------------*/
    mess_vector mess_vector_from_mexmess (const mxArray *Amatlab);
    mxArray * mess_vector_to_mexmess (mess_vector A);

    /*-----------------------------------------------------------------------------
     *  direct select routines
     *-----------------------------------------------------------------------------*/
    void mex_mess_direct_lu_select(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_direct_chol_select(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_multidirect_select(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);

    /*-----------------------------------------------------------------------------
     *  informations about cmess
     *-----------------------------------------------------------------------------*/
    void mex_mess_have_amd(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_have_arpack(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_have_bzip2(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_have_cholmod(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_have_colamd(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_have_csparse(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_have_matio(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_have_mess64(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_have_mklpardiso(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_have_openmp(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_have_superlu(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_have_umfpack(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_have_zlib(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_is_debug(int nlhs, mxArray *plhs[], int nrsh, const mxArray *prhs[]);
    void mex_mess_version(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_version_major(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_version_minor(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_version_patch(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_version_verbose(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_git_id(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_git_branch(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

    /*-----------------------------------------------------------------------------
     *  lyap, care, lradi, lrnm, sylvester and glyap3 stuff
     *-----------------------------------------------------------------------------*/
    void mex_mess_lyap(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_care(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_dense_nm_gmpare(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_lrnm(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_lradi(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_sylvester_sparsedense(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_sylvester_dense(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_glyap(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_mess_gstein(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

    /*-----------------------------------------------------------------------------
     *  some test functions
     *-----------------------------------------------------------------------------*/
    void mex_test_mess_direct_cholpackage_t(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_test_mess_direct_lupackage_t(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_test_mess_equation_t(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_test_mess_memusage_t(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_test_mess_multidirect_t(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_test_mess_operation_t(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_test_mess_parameter_t(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_test_mess_residual_t(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_test_mess_norm_t(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_test_vector(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_test_matrix(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_test_status(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
    void mex_test_options(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);


#ifdef __cplusplus
};
#endif


/* @} */

#endif /* INTERFACE_MATLAB_H_ */

