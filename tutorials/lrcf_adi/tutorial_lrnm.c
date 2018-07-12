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
 * @file tutorials/lrcf_adi/tutorial_lrnm.c
 * @brief Tutorial for solving a (generalized) Riccati  Equation  using \ref mess_lrnm solver.
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
 * The solver  \ref mess_lrnm returns the solution \f$ X\f$ as a low-rank factor approximation \f$Z\f$ such that \f$ X\approx ZZ^{T}\f$.\n
 * In the last step we compute the relative residual
 * \f$\|A^TZZ^T + ZZ^TA - ZZ^TBB^{T}ZZ^T + C^{T}C \|_2/\|C^TC\|_2 \f$ or
 * \f$\|\|A^TZZ^TE + E^TZZ^TA - E^TZZ^TBB^{T}ZZ^TE + C^{T}C \|_2/\|C^TC\|_2 \f$
 * and clear the memory.
 *
 * ## 1. Include Header Files
 * We include standard @c C header files and in addtion the @c mess.h header file
 * to make the @mess library available.
 * @snippet "tutorials/lrcf_adi/tutorial_lrnm.c" HEADER
 * \n
 *
 * ## 2. Declare Matrices and Check the Number of Input Arguments.
 * We use the \ref mess_matrix structure to declare our matrices \f$ A,B,C,E,Z\f$.
 * Furthermore we need a @ref mess_equation, @ref mess_options and @ref mess_status structure.
 * @snippet "tutorials/lrcf_adi/tutorial_lrnm.c" DECLARE
 * By using the \c argc variable we can check the number of input arguments.
 * @snippet "tutorials/lrcf_adi/tutorial_lrnm.c" INPUT
 * \n
 *
 * ## 3. Init Matrices and Read Data
 * Before we are able to use the matrices we have to init the fields of the \ref mess_matrix structure variables using \ref mess_matrix_init and read the matrices from @mm file format.
 * @snippet "tutorials/lrcf_adi/tutorial_lrnm.c" INITMATRICES
 * \n
 *
 * ## 4. Init Options, Status and Equations
 * We use @ref mess_options_init to init the @ref mess_options structure.
 * The @ref mess_options.adi_shifts_paratype describes the shift parameter strategy.
 * @ref mess_options.adi_output makes @ref mess_lrnm verbose and with @ref mess_options.type
 * we can switch between the transposed and non-transposed Riccati Equation.
 *   |      Operation Type      | Riccati Equation                       | generalized Riccati Equation                  |
 *   |:------------------------:|:--------------------------------------:|:---------------------------------------------:|
 *   |  @ref MESS_OP_NONE       | \f$ AX + XA^T - XBB^TX + C^TC = 0\f$   | \f$ AXE^T + EXA^T - EXBB^TX^T + C^TC = 0\f$   |
 *   |  @ref MESS_OP_TRANSPOSE  | \f$ A^TX + XA - XC^TCX + BB^T = 0\f$   | \f$ A^TXE + E^TXA - E^TXC^TCXE + BB^T = 0\f$  |
 *
 * \n
 * @snippet "tutorials/lrcf_adi/tutorial_lrnm.c" INITOPTIONS
 * \n
 * Now we init the @ref mess_status structure.
 * @snippet "tutorials/lrcf_adi/tutorial_lrnm.c" INITSTATUS
 * \n
 * @ref mess_equation_riccati builds an @ref mess_equation structure for the (generalized) Riccati Equation.
 * @snippet "tutorials/lrcf_adi/tutorial_lrnm.c" INITEQUATION
 * \n
 *
 * ## 5. Solve the (Generalized) Riccati Equation
 * Now we use @ref mess_lrcfadi_parameter to compute a shift parameter set.
 * @ref mess_status structure @c stat is used to collect runtime and residual information.
 * @ref mess_lrnm solves the Riccati Equation and after convergence
 * the low-rank factor is stored in @c Z.
 * @snippet "tutorials/lrcf_adi/tutorial_lrnm.c" SOLVE
 * \n
 *
 * ## 6. Compute the Residual
* We declare two @c double variables for storing the right-hand side norm and the absolute residual.
* \ref mess_matrix_dynorm2 computes the norm of the right-hand side and \ref mess_lrcfadi_residual
* computes the desired 2-norms.
* @snippet "tutorials/lrcf_adi/tutorial_lrnm.c" RESIDUAL
* \n
*
* ## 7. Clear Memory
* Clearing the @ref mess_status, @ref mess_options, @ref mess_equation is done by the corresponding
* functions.
* @snippet "tutorials/lrcf_adi/tutorial_lrnm.c" CLEAR
* \n
*
*
* @sa @ref tutorials
*/

///@cond
///[HEADER]
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "mess/mess.h"
///[HEADER]

int main ( int argc, char **argv){

    ///[DECLARE]
    mess_matrix A,B,C,E,Z;
    mess_equation eqn;
    mess_options opt;
    mess_status stat;
    ///[DECLARE]

    ///[INPUT]
    if ( argc !=4 && argc != 5 ) {
        printf("usage: %s A.mtx B.mtx C.mtx [E.mtx] \n", argv[0]);
        return 1;
    }
    ///[INPUT]

    ///[INITMATRICES]
    mess_matrix_init(&A);
    mess_matrix_init(&B);
    mess_matrix_init(&C);
    mess_matrix_init(&Z);
    if(argc==5)
        mess_matrix_init(&E);

    mess_matrix_read_formated(argv[1],A,MESS_CSR);
    mess_matrix_read_formated(argv[2],B,MESS_DENSE);
    mess_matrix_read_formated(argv[3],C,MESS_DENSE);
    if(argc==5){
        mess_matrix_read_formated(argv[4],E,MESS_CSR);
    }
    ///[INITMATRICES]

    ///[INITOPTIONS]
    mess_options_init(&opt);
    opt->adi_shifts_paratype = MESS_LRCFADI_PARA_MINMAX;
    opt->adi_output = 1;
    opt->type = MESS_OP_TRANSPOSE;
    mess_options_print(opt);
    ///[INITOPTIONS]

    ///[INITSTATUS]
    mess_status_init(&stat);
    ///[INITSTATUS]

    ///[INITEQUATION]
    mess_equation_init(&eqn);
    if(argc==4){
        mess_equation_riccati(eqn,opt,A,NULL,B,C);
    }else{
        mess_equation_riccati(eqn,opt,A,E,B,C);
    }
    ///[INITEQUATION]

    ///[SOLVE]
    mess_parameter(eqn,opt,stat);
    mess_status_print(stat);
    mess_lrnm(eqn,opt,stat,Z);
    mess_status_print(stat);
    ///[SOLVE]

    ///[RESIDUAL]
    double nrm=1, nrmRHS=1;
    mess_matrix_dynorm2(C,&nrmRHS);
    mess_lrcfadi_residual(eqn,opt,Z,&nrm);
    printf("relative Residual= %10.15e\n",nrm/nrmRHS);
    ///[RESIDUAL]

    ///[CLEAR]
    mess_status_clear(&stat);
    mess_options_clear(&opt);
    mess_equation_clear(&eqn);
    if(argc==4){
        mess_matrix_clear(&A);
        mess_matrix_clear(&B);
        mess_matrix_clear(&C);
        mess_matrix_clear(&Z);
    }else{
        mess_matrix_clear(&A);
        mess_matrix_clear(&B);
        mess_matrix_clear(&C);
        mess_matrix_clear(&E);
        mess_matrix_clear(&Z);
    }
    return nrm/nrmRHS > 1e-9;
    ///[CLEAR]
}
///@endcond

