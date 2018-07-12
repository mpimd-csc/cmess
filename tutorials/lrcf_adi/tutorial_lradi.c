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
 * @file tutorials/lrcf_adi/tutorial_lradi.c
 * @brief Tutorial for solving a (generalized) Lyapunov  Equation  using \ref mess_lradi solver.
 * @author  @mbehr
 *
 * # Tutorial Solving a (generalized) Lyapunov Equation
 * In this tutorial we show solve the (generalized) Lyapunov Equation, i.e.
 * \f$ FX  + XF^T    = -GG^T\f$ or
 * \f$ FXE^T + EXF^T = -GG^T\f$.\n
 * We assume that \f$F,E\in \mathbb{R}^{n\times n} \f$ is large and sparse
 * and the spectrum
 * \f$ \Lambda(F)\subseteq \mathbb{C}^{-}\f$ or \f$ \Lambda(F,E)\subseteq \mathbb{C}^{-}\f$
 * is contained in the left open complex half-plane. \n
 * The right-hand side is defined by \f$G \in \mathbb{R}^{p\times n}\f$ and \f$p \ll n\f$. \n
 * The solver  \ref mess_lradi returns the solution \f$ X\f$ as a low-rank factor approximation \f$Z\f$ such that \f$X\approx ZZ^{T}\f$.\n
 * In the last step we compute the relative residual
 * \f$\|FZZ^T + ZZ^TF^T+GG^T\|_2/\|GG^T\|_2 \f$ or \f$\|FZZ^TE^T + EZZ^TF^T+GG^T\|_2/\|GG^T\|_2 \f$
 * and clear the memory.
 *
 * ## 1. Include Header Files
 * We include standard @c C header files and in addtion the @c mess.h header file
 * to make the @mess library available.
 * @snippet "tutorials/lrcf_adi/tutorial_lradi.c" HEADER
 * \n
 *
 * ## 2. Declare Matrices and Check the Number of Input Arguments.
 * We use the \ref mess_matrix structure to declare our matrices \f$ F,G,E,Z\f$.
 * Furthermore we need a @ref mess_equation, @ref mess_options and @ref mess_status structure.
 * @snippet "tutorials/lrcf_adi/tutorial_lradi.c" DECLARE
 * By using the \c argc variable we can check the number of input arguments.
 * @snippet "tutorials/lrcf_adi/tutorial_lradi.c" INPUT
 * \n
 *
 * ## 3. Init Matrices and Read Data
 * Before we are able to use the matrices we have to init the fields of the \ref mess_matrix structure variables using \ref mess_matrix_init and read the matrices from @mm file format.
 * @snippet "tutorials/lrcf_adi/tutorial_lradi.c" INITMATRICES
 * \n
 *
 * ## 4. Init Options, Status and Equations
 * We use @ref mess_options_init to init the @ref mess_options structure.
 * The @ref mess_options.adi_shifts_paratype describes the shift parameter strategy.
 * @ref mess_options.adi_output makes @ref mess_lradi verbose and with @ref mess_options.type
 * we can switch between the transposed and non-transposed Lyapunov Equation.
 *      Operation Type      | Lyapunov Equation         | generalized Lyapunov Equation
 *  :----------------------:|:-------------------------:|:----------------------------:
 *  @ref MESS_OP_NONE       | \f$ FX+XF^T=-GG^T\f$      | \f$ FXE^T+EXF^T=-GG^T\f$
 *  @ref MESS_OP_TRANSPOSE  | \f$ F^TX +XF = -G^TG \f$  | \f$ F^TXE + E^TXF = -G^T G \f$
 *
 * \n
 * @snippet "tutorials/lrcf_adi/tutorial_lradi.c" INITOPTIONS
 * \n
 * Now we init the @ref mess_status structure.
 * @snippet "tutorials/lrcf_adi/tutorial_lradi.c" INITSTATUS
 * \n
 * @ref mess_equation_lyap builds an @ref mess_equation structure for the (generalized) Lyapunov Equation.
 * @snippet "tutorials/lrcf_adi/tutorial_lradi.c" INITEQUATION
 * \n
 *
 * ## 5. Solve the (Generalized) Lyapunov Equation
 * Now we use @ref mess_lrcfadi_parameter to compute a shift parameter set.
 * @ref mess_status structure @c stat is used to collect runtime and residual information.
 * @ref mess_lradi solves the Lyapunov Equation and after convergence
 * the low-rank factor is stored in @c Z.
 * @snippet "tutorials/lrcf_adi/tutorial_lradi.c" SOLVE
 * \n
 *
 * ## 6. Compute the Residual
 * We declare two @c double variables for storing the right-hand side norm and the absolute residual.
 * \ref mess_matrix_dynorm2 computes the norm of the right-hand side and \ref mess_lrcfadi_residual
 * computes the desired 2-norms.
 * @snippet "tutorials/lrcf_adi/tutorial_lradi.c" RESIDUAL
* \n
*
* ## 7. Clear Memory
* Clearing the @ref mess_status, @ref mess_options, @ref mess_equation is done by the corresponding
* functions.
* @snippet "tutorials/lrcf_adi/tutorial_lradi.c" CLEAR
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
    mess_matrix F,G,E,Z;
    mess_equation eqn;
    mess_options opt;
    mess_status stat;
    ///[DECLARE]

    ///[INPUT]
    if ( argc !=3&& argc != 4 ) {
        printf("usage: %s F.mtx G.mtx [E.mtx] \n", argv[0]);
        return 1;
    }
    ///[INPUT]

    ///[INITMATRICES]
    mess_matrix_init(&F);
    mess_matrix_init(&G);
    mess_matrix_init(&Z);
    if(argc==4)
        mess_matrix_init(&E);

    mess_matrix_read(argv[1],F);
    mess_matrix_read(argv[2],G);
    printf("Matrix F:\n"); mess_matrix_printinfo(F);
    printf("Matrix G:\n"); mess_matrix_printinfo(G);
    if(argc==4){
        mess_matrix_read(argv[3],E);
        printf("Matrix E:\n");
        mess_matrix_printinfo(E);
    }
    ///[INITMATRICES]

    ///[INITOPTIONS]
    mess_options_init(&opt);
    opt->adi_shifts_paratype = MESS_LRCFADI_PARA_MINMAX;
    opt->adi_output = 1;
    opt->type = MESS_OP_NONE;
    mess_options_print(opt);
    ///[INITOPTIONS]

    ///[INITSTATUS]
    mess_status_init(&stat);
    ///[INITSTATUS]

    ///[INITEQUATION]
    mess_equation_init(&eqn);
    if(argc==3){
        mess_equation_lyap(eqn,opt,F,NULL,G);
    }else{
        mess_equation_lyap(eqn,opt,F,E,G);
    }
    ///[INITEQUATION]

    ///[SOLVE]
    mess_parameter(eqn,opt,stat);
    mess_status_print(stat);
    mess_lradi(eqn,opt,stat,Z);
    mess_status_print(stat);
    ///[SOLVE]

    ///[RESIDUAL]
    double nrm=1, nrmRHS=1;
    mess_matrix_dynorm2(G,&nrmRHS);
    mess_lrcfadi_residual(eqn,opt,Z,&nrm);
    printf("relative Residual= %10.15e\n",nrm/nrmRHS);
    ///[RESIDUAL]

    ///[CLEAR]
    mess_status_clear(&stat);
    mess_options_clear(&opt);
    mess_equation_clear(&eqn);
    if(argc==3){
        mess_matrix_clear(&F);
        mess_matrix_clear(&G);
        mess_matrix_clear(&Z);
    }else{
        mess_matrix_clear(&F);
        mess_matrix_clear(&G);
        mess_matrix_clear(&E);
        mess_matrix_clear(&Z);
    }
    return nrm/nrmRHS > 1e-9;
    ///[CLEAR]
}
///@endcond

