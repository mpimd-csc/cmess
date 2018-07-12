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
 * @file matlab/c_interface/glyap3/mex_mess_glyap.c
 * @brief
 * @author @mbehr
 */

#include "interface_matlab.h"

void mex_mess_glyap( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    init_mexmess();

    mess_matrix A = NULL, E = NULL, Ahat = NULL, Ehat = NULL, QA = NULL, QE = NULL, Y = NULL, X = NULL;
    mess_operation_t op = MESS_OP_NONE;
    int ret = 0;

    MESS_INIT_MATRICES(&X, &Ahat, &Ehat, &QA, &QE);


    /*-----------------------------------------------------------------------------
     *  perform only solve call and solve and get schur decomposition call
     *-----------------------------------------------------------------------------*/
    if((nlhs == 1 || nlhs == 3 || nlhs == 5) && 3 == nrhs || (nrhs == 4 && mxIsClass(prhs[3], "mess_operation_t"))){

        /*-----------------------------------------------------------------------------
         *  read matrices and operation mode
         *-----------------------------------------------------------------------------*/
        A = mess_matrix_from_mexmess(prhs[0]);
        E = mess_matrix_from_mexmess(prhs[1]);
        Y = mess_matrix_from_mexmess(prhs[2]);

        if(nrhs == 4){
            op = mess_operation_t_from_mexmess(prhs[3]);
        }

        if(!A || !Y){
            csc_error_message("Cannot transfer Matrices from MEX-M.E.S.S. to C-M.E.S.S.\n");
            goto CLEAR;
        }

        if (E && nlhs == 3){
            csc_error_message("Invalid calling sequence.\n");
            goto CLEAR;
        }

        /*-----------------------------------------------------------------------------
         *  perform solve and get schur decomposition
         *-----------------------------------------------------------------------------*/
        ret = mess_glyap(op, A, E, Y, Ahat, QA, Ehat, QE, X);

        if (ret){
            csc_error_message("C-M.E.S.S.: mess_glyap returned = %d.\n",ret);
            goto CLEAR;
        }

    }else if(nlhs == 1 &&
            (nrhs == 4  || (nrhs == 5 && mxIsClass(prhs[4], "mess_operation_t"))  ||
             nrhs == 7  || (nrhs == 8 && mxIsClass(prhs[7], "mess_operation_t"))
            )
            ){

        if(nrhs == 4 || nrhs == 5){

            /*-----------------------------------------------------------------------------
             *  read matrices and operation mode
             *-----------------------------------------------------------------------------*/
            Ahat = mess_matrix_from_mexmess(prhs[1]);
            QA   = mess_matrix_from_mexmess(prhs[2]);
            Y    = mess_matrix_from_mexmess(prhs[3]);

            if(nrhs == 5){
                op = mess_operation_t_from_mexmess(prhs[4]);
            }

            if(!Ahat || !QA || !Y){
                csc_error_message("Cannot transfer Matrices from MEX-M.E.S.S. to C-M.E.S.S.\n");
                goto CLEAR;
            }

            /*-----------------------------------------------------------------------------
             *  perform solve and using schur decomposition
             *-----------------------------------------------------------------------------*/
            ret = mess_tglyap(op, Ahat, QA, NULL, NULL, Y, X);

            if (ret){
                csc_error_message("C-M.E.S.S.: mess_tglyap returned = %d.\n",ret);
                goto CLEAR;
            }

        }

        if(nrhs == 7 || nrhs == 8){

            /*-----------------------------------------------------------------------------
             *  read matrices and operation mode
             *-----------------------------------------------------------------------------*/
            Ahat = mess_matrix_from_mexmess(prhs[2]);
            Ehat = mess_matrix_from_mexmess(prhs[3]);
            QA   = mess_matrix_from_mexmess(prhs[4]);
            QE   = mess_matrix_from_mexmess(prhs[5]);
            Y    = mess_matrix_from_mexmess(prhs[6]);

            if(nrhs == 8){
                op = mess_operation_t_from_mexmess(prhs[7]);
            }

            if(!Ahat || !Ehat || !QA || !QE || !Y){
                csc_error_message("Cannot transfer Matrices from MEX-M.E.S.S. to C-M.E.S.S.\n");
                goto CLEAR;
            }

            /*-----------------------------------------------------------------------------
             *  perform solve and using schur decomposition
             *-----------------------------------------------------------------------------*/
            ret = mess_tglyap(op, Ahat, QA, Ehat, QE, Y, X);

        }

    }else{
        csc_error_message("Invalid calling sequence.\n");
        goto CLEAR;
    }

    /*-----------------------------------------------------------------------------
     *  write results to left hand sides
     *-----------------------------------------------------------------------------*/
    if( nlhs == 1 ){
        plhs[0] = mess_matrix_to_mexmess(X);
    }

    if( nlhs == 3 ){
        plhs[0] = mess_matrix_to_mexmess(X);
        plhs[1] = mess_matrix_to_mexmess(Ahat);
        plhs[2] = mess_matrix_to_mexmess(QA);
    }

    if( nlhs == 5 ){
        plhs[0] = mess_matrix_to_mexmess(X);
        plhs[1] = mess_matrix_to_mexmess(Ahat);
        plhs[2] = mess_matrix_to_mexmess(Ehat);
        plhs[3] = mess_matrix_to_mexmess(QA);
        plhs[4] = mess_matrix_to_mexmess(QE);
    }

CLEAR:
    if(A)       mess_matrix_clear(&A);
    if(E)       mess_matrix_clear(&E);
    if(Ahat)    mess_matrix_clear(&Ahat);
    if(Ehat)    mess_matrix_clear(&Ehat);
    if(QA)      mess_matrix_clear(&QA);
    if(QE)      mess_matrix_clear(&QE);
    if(Y)       mess_matrix_clear(&Y);
    if(X)       mess_matrix_clear(&X);
    return;
}


