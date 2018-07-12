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
 * @file lib/direct/singlesolver/sylvester_sparsedense.c
 * @brief Generate a solver for a sparse-dense Sylvester equation with a small second matrix.
 *
 * @author @koehlerm
 * @author @dykstra
 * @author @mbehr
 *
 * This code implements a sparse-dense Sylvester equation solver for the following Sylvester equations:
 * <ul>
 * <li> Standard:  \f$AX+XH+M=0\f$, \f$A^TX+XH^T+M=0\f$, \f$A^HX+XH^H+M=0\f$
 * <li> Half Generalized Sylvester: \f$AX+EXH+M=0\f$, \f$A^TX+E^TXH^T+M=0\f$, \f$A^HX+E^HXH^H+M=0\f$
 * <li> Generalized Sylvester: \f$AXF+EXH+M=0\f$, \f$A^TXF^T+E^TXH^T+M=0\f$, \f$A^HXF^H+E^HXH^H+M=0\f$
 * </ul>
 * with \f$ A \f$ and \f$ E \f$ sparse and large and \f$ F \f$ and \f$ H \f$ are small and dense.
 *
 * See @cite GarLAM92 and @cite morBenKS11 for references.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"



/**
 * @internal
 * @brief Enumeration to represent the type of a sparse-dense Sylvester equation.
 *
 * The @ref _sylv_equation_sd_t  enumeration is used to represent different sparse-dense Sylvester equations.
 *
 * @attention Internal use only.
 **/
typedef enum {
    /** Represents a Sylvester Equation \f$AX+XH+M=0\f$.*/
    SYLV_SPARSE_DENSE_STANDARD = 0,
    /** Represents a Sylvester Equation \f$AX+EXH+M=0\f$.*/
    SYLV_SPARSE_DENSE_HALF_GENERALIZED = 1,
    /** Represents a Sylvester Equation \f$AXF+EXH+M=0\f$.*/
    SYLV_SPARSE_DENSE_GENERALIZED = 2,
} _sylv_equation_sd_t;



/**
 * @internal
 * @brief Internal structure for sparse-dense Sylvester equation solver.
 *
 * @attention Internal use only.
 */
typedef struct _sylv_solver_sd_st {
    mess_matrix A;                          /** Matrix\f$ A \f$ of Sylvester Equation.*/
    mess_matrix Fhat;                       /** Schur Form of Matrix \f$ F \f$ of Sylvester Equation.*/
    mess_matrix Hhat;                       /** Schur Form of Matrix \f$ H \f$ of Sylvester Equation.*/
    mess_matrix Q;                          /** Schur Transformation Matrix of Matrix Pair \f$(F,H)\f$ of Sylvester Equation.*/
    mess_matrix Z;                          /** Schur Transformation Matrix of Matrix Pair \f$(F,H)\f$ of Sylvester Equation.*/
    mess_matrix E;                          /** Matrix \f$ E \f$ of Sylvester Equation. */
    mess_int_t rowsX;                       /** Rows of Solution \f$ X\f$ of Sylvester Equation.*/
    mess_int_t colsX;                       /** Columns of Solution \f$ X\f$ of Sylvester Equation.*/
    int isreal;                             /** @c 1 if all Matrices of Sylvester Equation are real, @c 0 otherwise.*/
    _sylv_equation_sd_t _sylv_eqn_sd_t;     /** Represents the Type of Sylvester Equation.*/
} _sylv_solver_sd;



/**
 * @internal
 * @brief Clear a @ref _sylv_solver_sd struct.
 * @param[in,out] solver solver data
 * @return always zero
 *
 * The @ref sylvester_sd_clear function is clears a @ref _sylv_solver_sd.
 *
 * @attention Internal use only.
 */
static int sylvester_sd_clear(void *solver){
    _sylv_solver_sd * sol = (_sylv_solver_sd*) solver;
    if ( sol != NULL) {
        if(sol->A) mess_matrix_clear(&(sol->A));
        if(sol->Fhat) mess_matrix_clear(&(sol->Fhat));
        if(sol->Hhat) mess_matrix_clear(&(sol->Hhat));
        if(sol->Q) mess_matrix_clear(&(sol->Q));
        if(sol->Z) mess_matrix_clear(&(sol->Z));
        if(sol->E) mess_matrix_clear(&(sol->E));
        mess_free( sol );
    }
    return 0;
}


/**
 * @internal
 * @brief Solve \f$ AX+XH+M=0 \f$, \f$AX+EXH+M=0\f$, \f$AXF+EXH+M=0\f$.
 * @param[in] datain  input solver data
 * @param[in] M  input right hand sides dense matrix.
 * @param[out] X solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref sylvester_sd_solvem function solves sparse-dense Sylvester equations:
 *  <ul>
 *  <li> \f$AX+XH+M=0\f$
 *  <li> \f$AX+EXH+M=0\f$
 *  <li> \f$AXF+EXH+M=0\f$
 *  </ul>
 *
 * @attention Internal use only.
 */
static int sylvester_sd_solvem(void* datain, mess_matrix M, mess_matrix X) {
    MSG_FNAME(__func__);
    int ret = 0, onebyoneblock=0;
    mess_int_t i,j,realM;
    mess_matrix Mtilde, Xtilde, E, CompEvFhat, CompEvHhat, F2, H2, Q2, Z2, HELP, MZ2;
    mess_vector rhs, xt;
    mess_direct LUsolver;
    _sylv_solver_sd  * data = (_sylv_solver_sd*) datain;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(M);
    mess_check_nullpointer(X);
    mess_check_dense(M);
    if ( M->rows != data->rowsX) {
        MSG_ERROR("number of rows doesn't match.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( M->cols != data->colsX){
        MSG_ERROR("number of columns doesn't match.\n");
        return MESS_ERROR_DIMENSION;
    }

    if (!(  data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_STANDARD          ||
                data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_HALF_GENERALIZED  ||
                data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_GENERALIZED)){
        MSG_ERROR("Unknown sparse-dense Sylvester equation.");
        return MESS_ERROR_NOSUPPORT;
    }


    /*-----------------------------------------------------------------------------
     *  prepare and init matrices
     *-----------------------------------------------------------------------------*/
    realM = MESS_IS_REAL(M);
    MESS_INIT_MATRICES(&Mtilde, &Xtilde, &E, &CompEvFhat, &CompEvHhat, &F2, &H2, &Q2, &Z2, &HELP, &MZ2);
    MESS_INIT_VECTORS(&rhs,&xt);

    /*-----------------------------------------------------------------------------
     *  handle the case of complex Eigenvalues for real matrices
     *-----------------------------------------------------------------------------*/
    if (data->isreal) {
        ret = mess_matrix_alloc(CompEvHhat, 2, 2, 2*2, MESS_DENSE, MESS_REAL);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(CompEvFhat, 2, 2, 2*2, MESS_DENSE, MESS_REAL);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(F2, 2, 2, 2*2, MESS_DENSE, MESS_COMPLEX);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(H2, 2, 2, 2*2, MESS_DENSE, MESS_COMPLEX);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(Q2, 2, 2, 2*2, MESS_DENSE, MESS_COMPLEX);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(Z2, 2, 2, 2*2, MESS_DENSE, MESS_COMPLEX);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(HELP, M->rows, 2, M->rows*2, MESS_DENSE, MESS_COMPLEX);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(MZ2, M->rows, 2, M->rows*2, MESS_DENSE, MESS_COMPLEX);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    }

    if(data->isreal && realM){
        ret = mess_matrix_alloc(Xtilde, M->rows, M->cols, M->rows*M->cols, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    } else {
        ret = mess_matrix_alloc(Xtilde, M->rows, M->cols, M->rows*M->cols, MESS_DENSE, MESS_COMPLEX);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    }

    if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_GENERALIZED){
        ret = mess_matrix_multiply(MESS_OP_NONE, M, MESS_OP_NONE, data->Z, Mtilde);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    }else{
        ret = mess_matrix_multiply(MESS_OP_NONE, M, MESS_OP_NONE, data->Q, Mtilde);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    }
    ret = mess_matrix_scale(-1.0,Mtilde);                                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);

    /*-----------------------------------------------------------------------------
     *  main loop
     *-----------------------------------------------------------------------------*/
    for (j = 0; j < data->colsX; j++){

        /*-----------------------------------------------------------------------------
         *  prepare next iteration
         *-----------------------------------------------------------------------------*/
        if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_STANDARD){
            if(data->isreal && realM){
                ret = mess_matrix_eye(E, data->rowsX, data->rowsX, MESS_CSC);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
            }else{
                ret = mess_matrix_eyec(E, data->rowsX, data->rowsX, MESS_CSC);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eyec);
            }
        }else{
            ret = mess_matrix_copy(data->E, E);                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        }

        ret = mess_matrix_getcol(Mtilde, j, rhs);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);

        /*-----------------------------------------------------------------------------
         *  check if 1-by-1 block is on diagonal
         *-----------------------------------------------------------------------------*/
        if(data->isreal){
            if(data->_sylv_eqn_sd_t != SYLV_SPARSE_DENSE_GENERALIZED){
                onebyoneblock = (j == data->colsX-1 || data->Hhat->values[j+1+j*data->Hhat->ld] == 0);
            }else{
                onebyoneblock = (j == data->colsX - 1 || (data->Fhat->values[j+1+j*data->Fhat->ld] == 0 && data->Hhat->values[j+1+j*data->Hhat->ld] == 0));
            }
        }else{
            //complex data 1-by-1 block on diagonal because of complex schur.
            onebyoneblock = 1;
        }

        if(onebyoneblock){
            /*-----------------------------------------------------------------------------
             * Real eigenvalue and data->isreal:
             * Solve (F(j,j)A + S(j,j)*E) xt = rhs, for xt.
             * or
             * Full complex Case. No 2x2 Block in Diagonal of S.
             *-----------------------------------------------------------------------------*/
            if(data->_sylv_eqn_sd_t != SYLV_SPARSE_DENSE_GENERALIZED){
                if(data->isreal){
                    ret = mess_matrix_addc(1,data->A,data->Hhat->values[j+j*data->Hhat->ld],E);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
                }else{
                    ret = mess_matrix_addc(1,data->A,data->Hhat->values_cpx[j+j*data->Hhat->ld],E);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
                }
            }else{
                if(data->isreal){
                    ret = mess_matrix_addc(data->Fhat->values[j+j*data->Fhat->ld], data->A, data->Hhat->values[j+j*data->Hhat->ld], E);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
                }else{
                    ret = mess_matrix_addc(data->Fhat->values_cpx[j+j*data->Fhat->ld], data->A, data->Hhat->values_cpx[j+j*data->Hhat->ld], E);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
                }
            }
            ret = mess_direct_init(&LUsolver);                                                                                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
            ret = mess_direct_create_sparse_lu(E,LUsolver);                                                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sparse_lu);
            ret = mess_direct_solve(MESS_OP_NONE,LUsolver,rhs,xt);                                                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solve);
        } else{
            /*-----------------------------------------------------------------------------
             * Complex eigenvalues and data->isreal:
             * The Idea is to solve
             * A*XTilde(:,j:j+1)*CompEvFhat + E*XTilde(:,j:j+1)*CompEvHhat = MTilde(:,j:j+1)
             * for XTilde.
             * We have a 2x2 Block on Diagonal of S.
             *-----------------------------------------------------------------------------*/
            if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_GENERALIZED){
                ret = mess_matrix_sub(data->Fhat,j,j+1,j,j+1,CompEvFhat);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sub);
                ret = mess_matrix_sub(data->Hhat,j,j+1,j,j+1,CompEvHhat);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sub);
                ret = mess_eigen_gschur_complex(CompEvFhat, CompEvHhat, F2, H2, Q2, Z2, NULL, NULL, NULL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_gschur_complex);
            }else{
                ret = mess_matrix_sub(data->Hhat,j,j+1,j,j+1,CompEvHhat);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sub);
                ret = mess_eigen_schur_complex(CompEvHhat, H2, Z2, NULL);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_schur_complex);
            }

            ret = mess_matrix_setcol(HELP, 0, rhs);                                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
            ret = mess_matrix_getcol(Mtilde, j+1, rhs);                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
            ret = mess_matrix_setcol(HELP, 1, rhs);                                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
            ret = mess_matrix_multiply(MESS_OP_NONE, HELP, MESS_OP_NONE, Z2, MZ2);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_getcol(MZ2, 0, rhs);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);

            /*-----------------------------------------------------------------------------
             * Complex eigenvalues and data->isreal:
             * Solve (F(0,0)*A + H2(0,0)*E) * xt = rhs2, for xt.
             *-----------------------------------------------------------------------------*/
            ret = mess_direct_init(&LUsolver);                                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
            if(data->_sylv_eqn_sd_t != SYLV_SPARSE_DENSE_GENERALIZED){
                ret = mess_matrix_addc(1,data->A,H2->values_cpx[0],E);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
            }else{
                ret = mess_matrix_addc(F2->values_cpx[0], data->A, H2->values_cpx[0], E);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
            }
            ret = mess_direct_create_sparse_lu(E,LUsolver);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sparse_lu);
            ret = mess_direct_solve(MESS_OP_NONE,LUsolver,rhs,xt);                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solve);
            ret = mess_matrix_setcol(HELP, 0, xt);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);

            /*-----------------------------------------------------------------------------
             *  update right hand side MZ2 with solution xt and get new right hand side
             *-----------------------------------------------------------------------------*/
            if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_STANDARD){
                ret = mess_matrix_colaxpy(-H2->values_cpx[0+1*H2->ld],xt,1,MZ2);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
            }else if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_HALF_GENERALIZED){
                ret = mess_matrix_mvp(MESS_OP_NONE, data->E, xt, rhs);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                ret = mess_matrix_colaxpy(-H2->values_cpx[0+1*H2->ld],rhs,1,MZ2);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
            }else{
                ret = mess_matrix_mvp(MESS_OP_NONE, data->A, xt, rhs);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                ret = mess_matrix_colaxpy(-(F2->values_cpx[0+1*F2->ld]), rhs, 1, MZ2);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                ret = mess_matrix_mvp(MESS_OP_NONE, data->E, xt, rhs);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                ret = mess_matrix_colaxpy(-(H2->values_cpx[0+1*H2->ld]), rhs, 1, MZ2);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
            }
            ret = mess_matrix_getcol(MZ2, 1, rhs);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);

            /*-----------------------------------------------------------------------------
             * Solve for second column
             * diagonal element is the conjugated EV for standard
             * and half generalized case.
             * F2(1,1)*A+H2(1,1)*E)Xhat_j+1= RHS - (F2(1,2)A+H2(1,2)E)Xhat_j
             *-----------------------------------------------------------------------------*/
            if(data->_sylv_eqn_sd_t != SYLV_SPARSE_DENSE_GENERALIZED){
                ret = mess_vector_conj(rhs);                                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_conj);
                ret = mess_direct_solve(MESS_OP_NONE,LUsolver,rhs,xt);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solve);
                ret = mess_vector_conj(xt);                                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_conj);
            }else{
                ret = mess_direct_clear(&LUsolver);                                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_clear);
                ret = mess_matrix_copy(data->E,E);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
                ret = mess_matrix_addc(F2->values_cpx[1+1*F2->ld], data->A, H2->values_cpx[1+1*H2->ld], E);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
                ret = mess_direct_init(&LUsolver);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
                ret = mess_direct_create_sparse_lu(E,LUsolver);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sparse_lu);
                ret = mess_direct_solve(MESS_OP_NONE,LUsolver,rhs,xt);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solve);
            }
            ret = mess_matrix_setcol(HELP, 1, xt);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);

            /*-----------------------------------------------------------------------------
             * get the resulting columns of Xtilde in MZ2
             *-----------------------------------------------------------------------------*/
            if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_GENERALIZED){
                ret = mess_matrix_multiply(MESS_OP_NONE, HELP, MESS_OP_HERMITIAN, Q2, MZ2);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            }else{
                ret = mess_matrix_multiply(MESS_OP_NONE, HELP, MESS_OP_HERMITIAN, Z2, MZ2);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            }

            ret = mess_matrix_getcol(MZ2, 0, xt);                                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);

            if(realM){
                ret = mess_vector_toreal_nowarn(xt);                                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
            }
            ret = mess_matrix_setcol(Xtilde, j, xt);                                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);

            /*-----------------------------------------------------------------------------
             * Here the columns of MTilde are updated with the new solution xt.
             *-----------------------------------------------------------------------------*/
            for (i=j+1; i < data->colsX ; i++){
                if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_STANDARD){
                    if(data->isreal){
                        ret = mess_matrix_colaxpy(-data->Hhat->values[j+i*data->Hhat->ld],xt,i,Mtilde);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    }else{
                        ret = mess_matrix_colaxpy(-data->Hhat->values_cpx[j+i*data->Hhat->ld],xt,i,Mtilde);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    }
                }else if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_HALF_GENERALIZED){
                    if(data->isreal){
                        ret = mess_matrix_mvp(MESS_OP_NONE, data->E, xt, rhs);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                        ret = mess_matrix_colaxpy(-data->Hhat->values[j+i*data->Hhat->ld],rhs,i,Mtilde);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    }else{
                        ret = mess_matrix_mvp(MESS_OP_NONE, data->E, xt, rhs);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                        ret = mess_matrix_colaxpy(-data->Hhat->values_cpx[j+i*data->Hhat->ld],rhs,i,Mtilde);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    }
                }else{
                    if(data->isreal){
                        ret = mess_matrix_mvp(MESS_OP_NONE, data->A, xt, rhs);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                        ret = mess_matrix_colaxpy(-data->Fhat->values[j+i*data->Fhat->ld],rhs,i,Mtilde);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                        ret = mess_matrix_mvp(MESS_OP_NONE, data->E, xt, rhs);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                        ret = mess_matrix_colaxpy(-data->Hhat->values[j+i*data->Hhat->ld],rhs,i,Mtilde);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    }else{
                        ret = mess_matrix_mvp(MESS_OP_NONE, data->A, xt, rhs);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                        ret = mess_matrix_colaxpy(-data->Fhat->values_cpx[j+i*data->Fhat->ld],rhs,i,Mtilde);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                        ret = mess_matrix_mvp(MESS_OP_NONE, data->E, xt, rhs);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                        ret = mess_matrix_colaxpy(-data->Hhat->values_cpx[j+i*data->Hhat->ld],rhs,i,Mtilde);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    }
                }
            }

            ret = mess_matrix_getcol(MZ2, 1, xt);                                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
            j++;
        }
        ret = mess_matrix_setcol(Xtilde, j, xt);                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);


        /*-----------------------------------------------------------------------------
         * Here the columns of MTilde are updated with the new solution xt.
         *-----------------------------------------------------------------------------*/
        for (i=j+1; i < data->colsX ; i++){
            if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_STANDARD){
                if(data->isreal){
                    ret = mess_matrix_colaxpy(-data->Hhat->values[j+i*data->Hhat->ld],xt,i,Mtilde);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                }else{
                    ret = mess_matrix_colaxpy(-data->Hhat->values_cpx[j+i*data->Hhat->ld],xt,i,Mtilde);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                }
            }else if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_HALF_GENERALIZED){
                if(data->isreal){
                    ret = mess_matrix_mvp(MESS_OP_NONE, data->E, xt, rhs);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_colaxpy(-data->Hhat->values[j+i*data->Hhat->ld],rhs,i,Mtilde);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                }else{
                    ret = mess_matrix_mvp(MESS_OP_NONE, data->E, xt, rhs);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_colaxpy(-data->Hhat->values_cpx[j+i*data->Hhat->ld],rhs,i,Mtilde);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                }
            }else{
                if(data->isreal){
                    ret = mess_matrix_mvp(MESS_OP_NONE, data->A, xt, rhs);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_colaxpy(-data->Fhat->values[j+i*data->Fhat->ld],rhs,i,Mtilde);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    ret = mess_matrix_mvp(MESS_OP_NONE, data->E, xt, rhs);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_colaxpy(-data->Hhat->values[j+i*data->Hhat->ld],rhs,i,Mtilde);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                }else{
                    ret = mess_matrix_mvp(MESS_OP_NONE, data->A, xt, rhs);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_colaxpy(-data->Fhat->values_cpx[j+i*data->Fhat->ld],rhs,i,Mtilde);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    ret = mess_matrix_mvp(MESS_OP_NONE, data->E, xt, rhs);                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_colaxpy(-data->Hhat->values_cpx[j+i*data->Hhat->ld],rhs,i,Mtilde);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                }
            }
        }
        mess_direct_clear(&LUsolver);
    }

    ret = mess_matrix_multiply(MESS_OP_NONE, Xtilde, MESS_OP_HERMITIAN, data->Q, X);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);

    /*-----------------------------------------------------------------------------
     *  clear memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&Mtilde, &Xtilde, &E, &CompEvFhat, &CompEvHhat, &F2, &H2, &Q2, &Z2, &HELP, &MZ2);
    MESS_CLEAR_VECTORS(&rhs,&xt);

    return 0;
}



/**
 * @internal
 * @brief Solve \f$ A^HX+XH^H+M=0 \f$, \f$A^HX+E^HXH^H+M=0\f$, \f$A^HXF^H+E^HXH^H+M=0\f$.
 * @param[in] datain  input solver data
 * @param[in] M  input right hand sides dense matrix.
 * @param[out] X solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref sylvester_sd_solvemh function solves sparse-dense Sylvester equations:
 *  <ul>
 *  <li> \f$A^HX+XH^H+M=0\f$
 *  <li> \f$A^HX+E^HXH^H+M=0\f$
 *  <li> \f$A^HXF^H+E^HXH^H+M=0\f$
 *  </ul>
 *
 * @attention Internal use only.
 */

static int sylvester_sd_solvemh(void* datain, mess_matrix M, mess_matrix X) {
    MSG_FNAME(__func__);
    int ret = 0, onebyoneblock=0;
    mess_int_t i,j,realM;
    mess_matrix Mtilde, Xtilde, E, CompEvFhat, CompEvHhat, F2, H2, Q2, Z2, HELP, MZ2;
    mess_vector rhs, xt;
    mess_direct LUsolver;
    _sylv_solver_sd  * data = (_sylv_solver_sd*) datain;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(data);
    mess_check_nullpointer(M);
    mess_check_nullpointer(X);
    mess_check_dense(M);
    if ( M->rows != data->rowsX) {
        MSG_ERROR("number of rows doesn't match.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( M->cols != data->colsX){
        MSG_ERROR("number of columns doesn't match.\n");
        return MESS_ERROR_DIMENSION;
    }

    if (!(  data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_STANDARD          ||
                data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_HALF_GENERALIZED  ||
                data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_GENERALIZED)){
        MSG_ERROR("Unknown sparse-dense Sylvester equation.");
        return MESS_ERROR_NOSUPPORT;
    }


    /*-----------------------------------------------------------------------------
     *  prepare and init matrices
     *-----------------------------------------------------------------------------*/
    realM = MESS_IS_REAL(M);
    MESS_INIT_MATRICES(&Mtilde, &Xtilde, &E, &CompEvFhat, &CompEvHhat, &F2, &H2, &Q2, &Z2, &HELP, &MZ2);
    MESS_INIT_VECTORS(&rhs,&xt);

    /*-----------------------------------------------------------------------------
     *  handle the case of complex Eigenvalues for real matrices
     *-----------------------------------------------------------------------------*/
    if (data->isreal) {
        ret = mess_matrix_alloc(CompEvHhat, 2, 2, 2*2, MESS_DENSE, MESS_REAL);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(CompEvFhat, 2, 2, 2*2, MESS_DENSE, MESS_REAL);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(F2, 2, 2, 2*2, MESS_DENSE, MESS_COMPLEX);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(H2, 2, 2, 2*2, MESS_DENSE, MESS_COMPLEX);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(Q2, 2, 2, 2*2, MESS_DENSE, MESS_COMPLEX);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(Z2, 2, 2, 2*2, MESS_DENSE, MESS_COMPLEX);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(HELP, M->rows, 2, M->rows*2, MESS_DENSE, MESS_COMPLEX);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(MZ2, M->rows, 2, M->rows*2, MESS_DENSE, MESS_COMPLEX);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    }

    if(data->isreal && realM){
        ret = mess_matrix_alloc(Xtilde, M->rows, M->cols, M->rows*M->cols, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    } else {
        ret = mess_matrix_alloc(Xtilde, M->rows, M->cols, M->rows*M->cols, MESS_DENSE, MESS_COMPLEX);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
    }

    if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_GENERALIZED){
        ret = mess_matrix_multiply(MESS_OP_NONE, M, MESS_OP_NONE, data->Q, Mtilde);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    }else{
        ret = mess_matrix_multiply(MESS_OP_NONE, M, MESS_OP_NONE, data->Q, Mtilde);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    }
    ret = mess_matrix_scale(-1.0,Mtilde);                                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_scale);

    /*-----------------------------------------------------------------------------
     *  main loop
     *-----------------------------------------------------------------------------*/
    for (j = data->colsX-1; j>=0; j--){

        /*-----------------------------------------------------------------------------
         *  prepare next iteration
         *-----------------------------------------------------------------------------*/
        if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_STANDARD){
            if(data->isreal && realM){
                ret = mess_matrix_eye(E, data->rowsX, data->rowsX, MESS_CSC);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
            }else{
                ret = mess_matrix_eyec(E, data->rowsX, data->rowsX, MESS_CSC);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eyec);
            }
        }else{
            ret = mess_matrix_copy(data->E, E);                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        }

        ret = mess_matrix_getcol(Mtilde, j, rhs);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);

        /*-----------------------------------------------------------------------------
         *  check if 1-by-1 block is on diagonal
         *-----------------------------------------------------------------------------*/
        if(data->isreal){
            if(data->_sylv_eqn_sd_t != SYLV_SPARSE_DENSE_GENERALIZED){
                onebyoneblock = (j == 0 || data->Hhat->values[j+(j-1)*data->Hhat->ld] == 0);
            }else{
                onebyoneblock = (j == 0 || (data->Fhat->values[j+(j-1)*data->Fhat->ld] == 0 && data->Hhat->values[j+(j-1)*data->Hhat->ld] == 0));
            }
        }else{
            //complex data 1-by-1 block on diagonal because of complex schur.
            onebyoneblock = 1;
        }

        if(onebyoneblock){
            /*-----------------------------------------------------------------------------
             * Real eigenvalue and data->isreal:
             * Solve (F(j,j)A + S(j,j)*E) xt = rhs, for xt.
             * or
             * Full complex Case. No 2x2 Block in Diagonal of S.
             *-----------------------------------------------------------------------------*/
            if(data->_sylv_eqn_sd_t != SYLV_SPARSE_DENSE_GENERALIZED){
                if(data->isreal){
                    ret = mess_matrix_addc(1,data->A,data->Hhat->values[j+j*data->Hhat->ld],E);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
                }else{
                    ret = mess_matrix_addc(1,data->A,data->Hhat->values_cpx[j+j*data->Hhat->ld],E);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
                }
            }else{
                if(data->isreal){
                    ret = mess_matrix_addc(data->Fhat->values[j+j*data->Fhat->ld], data->A, data->Hhat->values[j+j*data->Hhat->ld], E);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
                }else{
                    ret = mess_matrix_addc(data->Fhat->values_cpx[j+j*data->Fhat->ld], data->A, data->Hhat->values_cpx[j+j*data->Hhat->ld], E);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
                }
            }
            ret = mess_direct_init(&LUsolver);                                                                                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
            ret = mess_direct_create_sparse_lu(E,LUsolver);                                                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sparse_lu);
            ret = mess_direct_solve(MESS_OP_HERMITIAN,LUsolver,rhs,xt);                                                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solve);
        } else{
            /*-----------------------------------------------------------------------------
             * Complex eigenvalues and data->isreal:
             * The Idea is to solve
             * A*XTilde(:,j:j+1)*CompEvFhat + E*XTilde(:,j:j+1)*CompEvHhat = MTilde(:,j:j+1)
             * for XTilde.
             * We have a 2x2 Block on Diagonal of S.
             *-----------------------------------------------------------------------------*/
            if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_GENERALIZED){
                CompEvFhat->values[0+0*CompEvFhat->ld] = data->Fhat->values[j-1+(j-1)*data->Fhat->ld];
                CompEvFhat->values[1+0*CompEvFhat->ld] = data->Fhat->values[j-1+j*data->Fhat->ld];
                CompEvFhat->values[0+1*CompEvFhat->ld] = data->Fhat->values[j+(j-1)*data->Fhat->ld];
                CompEvFhat->values[1+1*CompEvFhat->ld] = data->Fhat->values[j+j*data->Fhat->ld];
                CompEvHhat->values[0+0*CompEvHhat->ld] = data->Hhat->values[j-1+(j-1)*data->Hhat->ld];
                CompEvHhat->values[1+0*CompEvHhat->ld] = data->Hhat->values[j-1+j*data->Hhat->ld];
                CompEvHhat->values[0+1*CompEvHhat->ld] = data->Hhat->values[j+(j-1)*data->Hhat->ld];
                CompEvHhat->values[1+1*CompEvHhat->ld] = data->Hhat->values[j+j*data->Hhat->ld];
                ret = mess_eigen_gschur_complex(CompEvFhat, CompEvHhat, F2, H2, Q2, Z2, NULL, NULL, NULL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_gschur_complex);
            }else{
                CompEvHhat->values[0+0*CompEvHhat->ld] = data->Hhat->values[j-1+(j-1)*data->Hhat->ld];
                CompEvHhat->values[1+0*CompEvHhat->ld] = data->Hhat->values[j-1+j*data->Hhat->ld];
                CompEvHhat->values[0+1*CompEvHhat->ld] = data->Hhat->values[j+(j-1)*data->Hhat->ld];
                CompEvHhat->values[1+1*CompEvHhat->ld] = data->Hhat->values[j+j*data->Hhat->ld];
                ret = mess_eigen_schur_complex(CompEvHhat, H2, Z2, NULL);                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_schur_complex);
            }

            ret = mess_matrix_setcol(HELP, 1, rhs);                                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
            ret = mess_matrix_getcol(Mtilde, j-1, rhs);                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
            ret = mess_matrix_setcol(HELP, 0, rhs);                                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
            ret = mess_matrix_multiply(MESS_OP_NONE, HELP, MESS_OP_NONE, Z2, MZ2);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            ret = mess_matrix_getcol(MZ2, 0, rhs);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);

            /*-----------------------------------------------------------------------------
             * Complex eigenvalues and data->isreal:
             * Solve (F(0,0)*A^T + H2(0,0)*E^T) * xt = rhs2, for xt.
             *-----------------------------------------------------------------------------*/
            ret = mess_direct_init(&LUsolver);                                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
            if(data->_sylv_eqn_sd_t != SYLV_SPARSE_DENSE_GENERALIZED){
                ret = mess_matrix_addc(1,data->A,H2->values_cpx[0],E);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
            }else{
                ret = mess_matrix_addc(F2->values_cpx[0], data->A, H2->values_cpx[0], E);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
            }
            ret = mess_direct_create_sparse_lu(E,LUsolver);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sparse_lu);
            ret = mess_direct_solve(MESS_OP_TRANSPOSE,LUsolver,rhs,xt);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solve);
            ret = mess_matrix_setcol(HELP, 0, xt);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);

            /*-----------------------------------------------------------------------------
             *  update right hand side MZ2 with solution xt and get new right hand side
             *-----------------------------------------------------------------------------*/
            if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_STANDARD){
                ret = mess_matrix_colaxpy(-H2->values_cpx[0+1*H2->ld],xt,1,MZ2);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
            }else if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_HALF_GENERALIZED){
                ret = mess_matrix_mvp(MESS_OP_TRANSPOSE, data->E, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                ret = mess_matrix_colaxpy(-H2->values_cpx[0+1*H2->ld],rhs,1,MZ2);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
            }else{
                ret = mess_matrix_mvp(MESS_OP_TRANSPOSE, data->A, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                ret = mess_matrix_colaxpy(-(F2->values_cpx[0+1*F2->ld]), rhs, 1, MZ2);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                ret = mess_matrix_mvp(MESS_OP_TRANSPOSE, data->E, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                ret = mess_matrix_colaxpy(-(H2->values_cpx[0+1*H2->ld]), rhs, 1, MZ2);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
            }
            ret = mess_matrix_getcol(MZ2, 1, rhs);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);

            /*-----------------------------------------------------------------------------
             * Solve for second column
             * diagonal element is the conjugated EV for standard
             * and half generalized case.
             * F2(1,1)*A+H2(1,1)*E)Xhat_j+1= RHS - (F2(1,2)A+H2(1,2)E)Xhat_j
             *-----------------------------------------------------------------------------*/
            if(data->_sylv_eqn_sd_t != SYLV_SPARSE_DENSE_GENERALIZED){
                ret = mess_direct_solve(MESS_OP_HERMITIAN,LUsolver,rhs,xt);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solve);
            }else{
                ret = mess_direct_clear(&LUsolver);                                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_clear);
                ret = mess_matrix_copy(data->E,E);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
                ret = mess_matrix_addc(F2->values_cpx[1+1*F2->ld], data->A, H2->values_cpx[1+1*H2->ld], E);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_addc);
                ret = mess_direct_init(&LUsolver);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_init);
                ret = mess_direct_create_sparse_lu(E,LUsolver);                                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_create_sparse_lu);
                ret = mess_direct_solve(MESS_OP_NONE,LUsolver,rhs,xt);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_solve);
            }
            ret = mess_matrix_setcol(HELP, 1, xt);                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);

            /*-----------------------------------------------------------------------------
             * get the resulting columns of Xtilde in MZ2
             *-----------------------------------------------------------------------------*/
            if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_GENERALIZED){
                ret = mess_matrix_multiply(MESS_OP_NONE, HELP, MESS_OP_HERMITIAN, Q2, MZ2);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            }else{
                ret = mess_matrix_multiply(MESS_OP_NONE, HELP, MESS_OP_HERMITIAN, Z2, MZ2);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
            }

            ret = mess_matrix_getcol(MZ2, 1, xt);                                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);

            if(realM){
                ret = mess_vector_toreal_nowarn(xt);                                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
            }
            ret = mess_matrix_setcol(Xtilde, j, xt);                                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);

            /*-----------------------------------------------------------------------------
             * Here the columns of MTilde are updated with the new solution xt.
             *-----------------------------------------------------------------------------*/
            for (i=j-1; i >= 0 ; i--){
                if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_STANDARD){
                    if(data->isreal){
                        ret = mess_matrix_colaxpy(-data->Hhat->values[i+j*data->Hhat->ld],xt,i,Mtilde);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    }else{
                        ret = mess_matrix_colaxpy(-conj(data->Hhat->values_cpx[i+j*data->Hhat->ld]),xt,i,Mtilde);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    }
                }else if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_HALF_GENERALIZED){
                    if(data->isreal){
                        ret = mess_matrix_mvp(MESS_OP_HERMITIAN, data->E, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                        ret = mess_matrix_colaxpy(-data->Hhat->values[i+j*data->Hhat->ld],rhs,i,Mtilde);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    }else{
                        ret = mess_matrix_mvp(MESS_OP_HERMITIAN, data->E, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                        ret = mess_matrix_colaxpy(-conj(data->Hhat->values_cpx[i+j*data->Hhat->ld]),rhs,i,Mtilde);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    }
                }else{
                    if(data->isreal){
                        ret = mess_matrix_mvp(MESS_OP_HERMITIAN, data->A, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                        ret = mess_matrix_colaxpy(-data->Fhat->values[i+j*data->Fhat->ld],rhs,i,Mtilde);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                        ret = mess_matrix_mvp(MESS_OP_HERMITIAN, data->E, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                        ret = mess_matrix_colaxpy(-data->Hhat->values[i+j*data->Hhat->ld],rhs,i,Mtilde);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    }else{
                        ret = mess_matrix_mvp(MESS_OP_HERMITIAN, data->A, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                        ret = mess_matrix_colaxpy(-conj(data->Fhat->values_cpx[i+j*data->Fhat->ld]),rhs,i,Mtilde);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                        ret = mess_matrix_mvp(MESS_OP_HERMITIAN, data->E, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                        ret = mess_matrix_colaxpy(-conj(data->Hhat->values_cpx[i+j*data->Hhat->ld]),rhs,i,Mtilde);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    }
                }
            }

            ret = mess_matrix_getcol(MZ2, 0, xt);                                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
            if(realM){
                ret = mess_vector_toreal_nowarn(xt);                                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
            }
            j--;
        }
        ret = mess_matrix_setcol(Xtilde, j, xt);                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);


        /*-----------------------------------------------------------------------------
         * Here the columns of MTilde are updated with the new solution xt.
         *-----------------------------------------------------------------------------*/
        for (i=j-1; i >= 0 ; i--){
            if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_STANDARD){
                if(data->isreal){
                    ret = mess_matrix_colaxpy(-data->Hhat->values[i+j*data->Hhat->ld],xt,i,Mtilde);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                }else{
                    ret = mess_matrix_colaxpy(-conj(data->Hhat->values_cpx[i+j*data->Hhat->ld]),xt,i,Mtilde);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                }
            }else if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_HALF_GENERALIZED){
                if(data->isreal){
                    ret = mess_matrix_mvp(MESS_OP_HERMITIAN, data->E, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_colaxpy(-data->Hhat->values[i+j*data->Hhat->ld],rhs,i,Mtilde);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                }else{
                    ret = mess_matrix_mvp(MESS_OP_HERMITIAN, data->E, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_colaxpy(-conj(data->Hhat->values_cpx[i+j*data->Hhat->ld]),rhs,i,Mtilde);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                }
            }else{
                if(data->isreal){
                    ret = mess_matrix_mvp(MESS_OP_HERMITIAN, data->A, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_colaxpy(-data->Fhat->values[i+j*data->Fhat->ld],rhs,i,Mtilde);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    ret = mess_matrix_mvp(MESS_OP_HERMITIAN, data->E, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_colaxpy(-data->Hhat->values[i+j*data->Hhat->ld],rhs,i,Mtilde);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                }else{
                    ret = mess_matrix_mvp(MESS_OP_HERMITIAN, data->A, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_colaxpy(-conj(data->Fhat->values_cpx[i+j*data->Fhat->ld]),rhs,i,Mtilde);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                    ret = mess_matrix_mvp(MESS_OP_HERMITIAN, data->E, xt, rhs);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_colaxpy(-conj(data->Hhat->values_cpx[i+j*data->Hhat->ld]),rhs,i,Mtilde);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colaxpy);
                }
            }
        }
        mess_direct_clear(&LUsolver);
    }

    if(data->_sylv_eqn_sd_t != SYLV_SPARSE_DENSE_GENERALIZED){
        ret = mess_matrix_multiply(MESS_OP_NONE, Xtilde, MESS_OP_HERMITIAN, data->Q, X);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    }else{
        ret = mess_matrix_multiply(MESS_OP_NONE, Xtilde, MESS_OP_HERMITIAN, data->Z, X);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_multiply);
    }

    /*-----------------------------------------------------------------------------
     *  clear memory
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&Mtilde, &Xtilde, &E, &CompEvFhat, &CompEvHhat, &F2, &H2, &Q2, &Z2, &HELP, &MZ2);
    MESS_CLEAR_VECTORS(&rhs,&xt);
    return 0;
}



/**
 * @internal
 * @brief Solve \f$ A^TX+XH^T+M=0 \f$, \f$A^TX+E^TXH^T+M=0\f$, \f$A^TXF^T+E^TXH^T+M=0\f$.
 * @param[in] datain  input solver data
 * @param[in] M  input right hand sides dense matrix.
 * @param[out] X solution
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref sylvester_sd_solvemt function solves sparse-dense Sylvester equations:
 *  <ul>
 *  <li> \f$A^TX+XH^T+M=0\f$
 *  <li> \f$A^TX+E^TXH^T+M=0\f$
 *  <li> \f$A^TXF^T+E^TXH^T+M=0\f$
 *  </ul>
 *
 * @attention Internal use only.
 */
static int sylvester_sd_solvemt(void* datain, mess_matrix M, mess_matrix X) {
    MSG_FNAME(__func__);
    int ret = 0;
    ret = mess_matrix_conj(M);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_conj);
    ret = sylvester_sd_solvemh(datain,M,X);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), sylvester_sd_solvemh);
    ret = mess_matrix_conj(X);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_conj);
    ret = mess_matrix_conj(M);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_conj);
    return ret;
}


/**
 * @brief Generate a solver for a sparse-dense Sylvester equation with a small second matrix.
 * @param[in] A  input first matrix
 * @param[in] F  input generalized second matrix
 * @param[in] E  input generalized matrix for the right term
 * @param[in] H  input second matrix
 * @param[out] solver solver to create
 * @return zero on success or a non-zero error value otherwise

 * The @ref mess_direct_create_sylvester_sparsedense function creates a solver for
 *
 * <ul>
 * <li> Standard:  \f$AX+XH+M=0\f$, \f$A^TX+XH^T+M=0\f$,\f$A^HX+XH^H+M=0\f$
 * <li> Half Generalized Sylvester: \f$AX+EXH+M=0\f$, \f$A^TX+E^TXH^T+M=0\f$, \f$A^HX+E^HXH^H+M=0\f$
 * <li> Generalized Sylvester: \f$AXF+EXH+M=0\f$, \f$A^TXF^T+E^TXH^T+M=0\f$, \f$A^HXF^H+E^HXH^H+M=0\f$
 * </ul>
 * with \f$ A \f$ and \f$ E \f$ sparse and large and \f$ F \f$ and \f$ H \f$ are small and dense.
 *
 *
 * If it gets @c NULL pointers for \f$ F \f$ and \f$ E \f$, it creates the solver for \f$ AX+XH+M=0 \f$.
 *
 * If only \f$ F \f$ is a @c NULL pointer, it creates the solver for \f$ AX+EXH+M=0 \f$.
 *
 * If none of the arguments are @c NULL, it creates the solver for \f$ AXF+EXH+M=0 \f$.
 *
 * See @cite GarLAM92 and  @cite morBenKS11 for references.
 *
 * @author @koehlerm
 * @author @dykstra
 * @author @mbehr
 *
 */
int mess_direct_create_sylvester_sparsedense(mess_matrix A, mess_matrix F, mess_matrix E, mess_matrix H, mess_direct solver){
    MSG_FNAME(__func__);
    int ret =0;
    _sylv_solver_sd* data;
    _sylv_equation_sd_t type = SYLV_SPARSE_DENSE_STANDARD;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(H);
    mess_check_nullpointer(solver);
    mess_check_square(A);
    mess_check_sparse(A);
    mess_check_square(H);
    mess_check_dense(H);
    mess_check_real_or_complex(A);
    mess_check_real_or_complex(H);
    if(E){
        mess_check_square(E);
        mess_check_sparse(E);
        mess_check_real_or_complex(E);
        mess_check_same_size(A,E);
        type = SYLV_SPARSE_DENSE_HALF_GENERALIZED;
    }
    if(F){
        mess_check_square(F);
        mess_check_dense(F);
        mess_check_real_or_complex(F);
        mess_check_same_size(F,H);
        type = SYLV_SPARSE_DENSE_GENERALIZED;
    }

    if(!E && F){
        MSG_ERROR("E is given but not F, this kind of sparse-dense Sylvester equation is not supported.\n");
        return MESS_ERROR_NOSUPPORT;
    }

    /*-----------------------------------------------------------------------------
     *  init matrices, fill fields of data
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(data, _sylv_solver_sd*, sizeof(_sylv_solver_sd));
    data->_sylv_eqn_sd_t = type;
    data->rowsX = A->rows;
    data->colsX = H->rows;
    MESS_INIT_MATRICES(&(data->A), &(data->E), &(data->Fhat), &(data->Hhat), &(data->Q), &(data->Z));

    ret = mess_matrix_convert(A, data->A, MESS_CSC);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
    if(data->_sylv_eqn_sd_t != SYLV_SPARSE_DENSE_STANDARD){
        ret = mess_matrix_convert(E, data->E, MESS_CSC);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
    }

    if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_STANDARD){
        data->isreal = (MESS_IS_REAL(A) && MESS_IS_REAL(H)) ? 1:0;
    }else if (data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_HALF_GENERALIZED){
        data->isreal = (MESS_IS_REAL(A) && MESS_IS_REAL(H) && MESS_IS_REAL(E)) ? 1:0;
    }else{
        data->isreal = (MESS_IS_REAL(A) && MESS_IS_REAL(H) && MESS_IS_REAL(E) && MESS_IS_REAL(F)) ? 1:0;
    }


    /*-----------------------------------------------------------------------------
     *  compute the schur decompositon
     *-----------------------------------------------------------------------------*/
    if(data->isreal){
        if(data->_sylv_eqn_sd_t != SYLV_SPARSE_DENSE_GENERALIZED){
            ret = mess_eigen_schur(H, data->Hhat, data->Q, NULL);                                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_schur);
        }else{
            ret = mess_eigen_gschur(F,H,data->Fhat,data->Hhat,data->Q, data->Z,NULL,NULL,NULL);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_gschur);
        }
    } else {
        if(data->_sylv_eqn_sd_t != SYLV_SPARSE_DENSE_GENERALIZED){
            ret = mess_eigen_schur_complex(H, data->Hhat, data->Q, NULL);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_schur_complex);
        }else{
            ret = mess_eigen_gschur_complex(F,H,data->Fhat,data->Hhat,data->Q, data->Z,NULL,NULL,NULL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_gschur_complex);
        }
    }


    /*-----------------------------------------------------------------------------
     *  setup solver
     *-----------------------------------------------------------------------------*/
    solver->data_type = data->isreal?MESS_REAL:MESS_COMPLEX;
    solver->rows = data->rowsX;
    solver->cols = data->colsX;
    if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_STANDARD){
        SET_SOLVERNAME(solver->name, "SYLV_SPARSE_DENSE_STANDARD");
    }else if(data->_sylv_eqn_sd_t == SYLV_SPARSE_DENSE_HALF_GENERALIZED){
        SET_SOLVERNAME(solver->name, "SYLV_SPARSE_DENSE_HALF_GENERALIZED");
    }else{
        SET_SOLVERNAME(solver->name, "SYLV_SPARSE_DENSE_GENERALIZED");
    }
    solver->data = (void *)data;
    solver->solve  = NULL;
    solver->solvem = NULL;
    solver->solvem = NULL;
    solver->solvem  = sylvester_sd_solvem;
    solver->solvemt = sylvester_sd_solvemt;
    solver->solvemh = sylvester_sd_solvemh;
    solver->clear = sylvester_sd_clear;
    solver->getL = NULL;
    solver->getU = NULL;
    solver->getpermp = NULL;
    solver->getpermq = NULL;
    solver->getscalerow = NULL;
    solver->getscalecol = NULL;
    solver->det = NULL;
    solver->detc = NULL;
    solver->inverse = NULL;

    return 0;
}    /* -----  end of function mess_direct_create_sylvester_sparsedense  ----- */
