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
 * @file matlab/c_interface/lrcfadi/mex_mess_dense_nm_gmpare.c
 * @brief
 * @author @mbehr
 * @author @dykstra
 */


#include "interface_matlab.h"

//    [X, absres, relres] = mess_call('mess_dense_nm_gmpare', p.Results.X0, A, E, Q, G, p.Results.plus, p.Results.linesearch, p.Results.op, p.Results.maxit, ...
//    p.Results.nrm, p.Results.absres_tol, p.Results.relres_tol);


void mex_mess_dense_nm_gmpare( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    init_mexmess();

    mess_matrix X0 = NULL;
    mess_matrix A = NULL;
    mess_matrix E = NULL;
    mess_matrix Q = NULL;
    mess_matrix G = NULL;
    mess_matrix X = NULL;
    mess_int_t plus, linesearch;
    mess_operation_t trans;
    mess_norm_t norm;
    mess_int_t maxit;
    double absres_tol, relres_tol;
    double absres, relres;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  read input arguments, input args are checked by matlab
     *-----------------------------------------------------------------------------*/
    if (nrhs != 12) {
        csc_error_message("The function needs exactly 12 input arguments but %d are given.\n", (int) nrhs);
    }


    X0 = mess_matrix_from_mexmess(prhs[0]);

    A = mess_matrix_from_mexmess(prhs[1]);
    if (!A) {
        if(X0) mess_matrix_clear(&X0);
        csc_error_message("The second argument is invalid.\n");
        return;
    }

    E = mess_matrix_from_mexmess(prhs[2]);

    Q = mess_matrix_from_mexmess(prhs[3]);
    if (!Q) {
        if(X0) mess_matrix_clear(&X0);
        if(E) mess_matrix_clear(&E);
        MESS_CLEAR_MATRICES(&A);
        csc_error_message("The fourth argument is invalid.\n");
        return;
    }

    G = mess_matrix_from_mexmess(prhs[4]);
    if (!G) {
        if(X0) mess_matrix_clear(&X0);
        if(E) mess_matrix_clear(&E);
        MESS_CLEAR_MATRICES(&A, &Q);
        csc_error_message("The fifth argument is invalid.\n");
        return;
    }

    plus        = (mess_int_t) mxGetScalar(prhs[5]);
    linesearch  = (mess_int_t) mxGetScalar(prhs[6]);
    trans       = mess_operation_t_from_mexmess(prhs[7]);
    maxit       = (mess_int_t) mxGetScalar(prhs[8]);
    norm        = mess_norm_t_from_mexmess(prhs[9]);
    absres_tol  = (double) mxGetScalar(prhs[10]);
    relres_tol  = (double) mxGetScalar(prhs[11]);


    /*-----------------------------------------------------------------------------
     *  call mess_dense_nm_gmpare
     *-----------------------------------------------------------------------------*/
    mess_matrix_init(&X);
    ret = mess_dense_nm_gmpare(X0, A, E, Q, G, plus, linesearch, trans, maxit, norm, absres_tol, relres_tol, &absres, &relres, NULL, X);

    MESS_CLEAR_MATRICES(&A,&Q,&G);
    if (E) mess_matrix_clear(&E);
    if (X0) mess_matrix_clear(&X0);

    if ( ret != 0 ) {
        mess_matrix_clear(&X);
        csc_error_message("Failed to solve the Riccati Equation.\n-> %d - %s \n", (int) ret, mess_get_error(ret));
        return;
    }

    plhs[0] = mess_matrix_to_mexmess(X);
    plhs[1] = mxCreateDoubleScalar(absres);
    plhs[2] = mxCreateDoubleScalar(relres);

    mess_matrix_clear(&X);
    return;
}


