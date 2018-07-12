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
 * @file tutorials/lrcf_adi/tutorial_care.c
 * @brief Tutorial for solving a (generalized) Riccati Equation  using \ref mess_care solver.
 * @author  @mbehr
 *
 * # Tutorial Solving a (generalized) Riccati Equation
 * In this tutorial we solve the (generalized) Riccati Equation, i.e.
 * \f$ A^TX +      XA - XBB^{T}X        + C^{T}C = 0\f$ or
 * \f$ A^TXE    + E^{T}XA - E^{T}XBB^{T}XE  + C^{T}C = 0\f$.\n
 * We assume that \f$A,E\in \mathbb{R}^{n\times n} \f$ is large and sparse
 * and the spectrum
 * \f$ \Lambda(A)\subseteq \mathbb{C}^{-}\f$ or \f$ \Lambda(A,E)\subseteq \mathbb{C}^{-}\f$
 * is contained in the left open complex half-plane. Under stability and detectability assumption
 * theory admits  existence of a unique solution of the (generalized) algebraic Riccati Equation and ensures convergence of
 * the algorithm. \n
 * The right-hand side is defined by \f$C \in \mathbb{R}^{p\times n}\f$ and \f$ B \in\mathbb{R}^{n\times m}\f$.
 * Furthermore the number of outputs and inputs should be small compared to the size of the system  \f$ p,m \ll n\f$. \n
 * The solver  \ref mess_care returns the solution \f$ X\f$ as a low-rank factor approximation \f$Z\f$ such that \f$ X\approx ZZ^{T}\f$.\n
 * In the last step we compute the relative residual
 * \f$\|A^TZZ^T + ZZ^TA - ZZ^TBB^{T}ZZ^T + C^{T}C \|_2/\|C^TC\|_2 \f$ or
 * \f$\|\|A^TZZ^TE + E^TZZ^TA - E^TZZ^TBB^{T}ZZ^TE + C^{T}C \|_2/\|C^TC\|_2 \f$
 * and clear the memory.
 *
 * ## 1. Include Header Files
 * We include standard @c C header files and in addtion the @c mess.h header file
 * to make the @mess library available.
 * @snippet "tutorials/lrcf_adi/tutorial_care.c" HEADER
 * \n
 *
 *
 * ## 2. Init Matrices and Read Data
 * Before we are able to use we have to init the fields of the \ref mess_matrix structure variables using \ref mess_matrix_init
 * @snippet "tutorials/lrcf_adi/tutorial_care.c" INIT
 * \n
 *
 * Now we can use \ref mess_matrix_read_formated to read our matrices from a Matrix Market File Format and convert it to
 * the storage type \ref MESS_CSR and \ref MESS_DENSE. \n
 * \ref mess_matrix_read_formated dynamical allocates memory at this point.
 * @snippet "tutorials/lrcf_adi/tutorial_care.c" READ
 * \n
 *
 * ## 3. Solve the (Generalized) Riccati Equation
 * As already announced \ref mess_care solves the (generalized) Riccati Equation.
 * We can indicate the case of the nongeneralized Riccati Equation by providing  @c NULL as  the second argument to \ref mess_care.
 * @snippet "tutorials/lrcf_adi/tutorial_care.c" SOLVE
 * \n
 *
 * ## 4. Compute the Relative Residual
 * We declare two @c double variables for storing the right-hand side norm and the absolute residual.
 * \ref mess_matrix_dynorm2 computes the norm of the right-hand side and \ref mess_lrcfadi_res2nm or \ref mess_lrcfadi_res2nmg
 * are helpfull to compute the desired 2-norms.
 * @snippet "tutorials/lrcf_adi/tutorial_care.c" RES
 * \n
 *
 * ## 5. Don't Forget to Clear Memory
 * A rule of thumb is that every structure which needs an @c init call before usage must be cleared at the
 * at the end of your code. In our cases we have used \ref mess_matrix_init which means we have to use
 * \ref mess_matrix_clear to clear dynamical allocated memory.
 * @snippet "tutorials/lrcf_adi/tutorial_care.c" CLEAR
 * \n
 *
 *
 * @sa @ref tutorials
 *
 */

//cond -> code is ignored by doxygen documentation, but snippet references are working


///@cond
///[HEADER]
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "mess/mess.h"
///[HEADER]

int main ( int argc, char **argv){

    /// [DECLARE]
    mess_matrix A,B,C,E,Z;
    /// [DECLARE]

    /// [INPUT]
    if(argc!=4 && argc!=5){
        printf("usage: %s A.mtx B.mtx C.mtx [E.mtx]\n",argv[0]);
        return 1;
    }
    /// [INPUT]

    /// [INIT]
    mess_matrix_init(&A);
    mess_matrix_init(&B);
    mess_matrix_init(&C);
    if(argc==5){mess_matrix_init(&E);}
    mess_matrix_init(&Z);
    /// [INIT]

    /// [READ]
    mess_matrix_read_formated(argv[1],A,MESS_CSR);
    mess_matrix_read_formated(argv[2],B,MESS_DENSE);
    mess_matrix_read_formated(argv[3],C,MESS_DENSE);
    if(argc==5){
        mess_matrix_read_formated(argv[4],E,MESS_CSR);
    }
    /// [READ]

    /// [SOLVE]
    if(argc==4){
        mess_care(A,NULL,B,C,Z);
    }else{
        mess_care(A,E,B,C,Z);
    }
    /// [SOLVE]

    /// [RES]
    double nrmRHS=1, nrm=1;
    mess_matrix_dynorm2(C,&nrmRHS);
    if(argc==4){
        mess_lrcfadi_res2nm(A,B,C,Z,&nrm);
    }else{
        mess_lrcfadi_res2nmg(A,B,C,E,Z,&nrm);
    }
    printf("relative Residual= %10.15e\n",nrm/nrmRHS);
    /// [RES]

    /// [CLEAR]
    mess_matrix_clear(&A);
    mess_matrix_clear(&B);
    mess_matrix_clear(&C);
    mess_matrix_clear(&Z);
    if(argc==5)
        mess_matrix_clear(&E);

    return nrm/nrmRHS>1e-9;
    /// [CLEAR]
}
///@endcond












