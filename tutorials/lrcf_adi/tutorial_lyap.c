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
 * @file tutorials/lrcf_adi/tutorial_lyap.c
 * @brief Tutorial for solving a (generalized) Lyapunov  Equation  using \ref mess_lyap solver.
 * @author  @mbehr
 *
 * # Tutorial Solving a (generalized) Lyapunov Equation
 * In this tutorial we solve the (generalized) Lyapunov Equation, i.e.
 * \f$ FX  + XF^T    = -GG^T\f$ or
 * \f$ FXE^T + EXF^T = -GG^T\f$.\n
 * We assume that \f$F,E\in \mathbb{R}^{n\times n} \f$ is large and sparse
 * and the spectrum
 * \f$ \Lambda(F)\subseteq \mathbb{C}^{-}\f$ or \f$ \Lambda(F,E)\subseteq \mathbb{C}^{-}\f$
 * is contained in the left open complex half-plane. \n
 * The right-hand side is defined by \f$G \in \mathbb{R}^{n\times p}\f$ and \f$p \ll n\f$. \n
 * The solver  \ref mess_lyap returns the solution \f$ X\f$ as a low-rank factor approximation \f$Z\f$ such that \f$X\approx ZZ^{T}\f$.\n
 * In the last step we compute the relative residual
 * \f$\|FZZ^T + ZZ^TF^T+GG^T\|_2/\|GG^T\|_2 \f$ or \f$\|FZZ^TE^T + EZZ^TF^T+GG^T\|_2/\|GG^T\|_2 \f$
 * and clear the memory.
 *
 * ## 1. Include Header Files
 * We include standard @c C header files and in addtion the @c mess.h header file
 * to make the @mess library available.
 * @snippet "tutorials/lrcf_adi/tutorial_lyap.c" HEADER
 * \n
 *
 * ## 2. Declare Matrices and Check the Number of Input Arguments.
 * We use the \ref mess_matrix structure to declare our matrices \f$ F,G,E,Z\f$.
 * @snippet "tutorials/lrcf_adi/tutorial_lyap.c" DECLARE
 * By using the \c argc variable we can check the number of input arguments.
 * @snippet "tutorials/lrcf_adi/tutorial_lyap.c" INPUT
 * \n
 *
 * ## 3. Init Matrices and Read Data
 * Before we are able to use we have to init the fields of the \ref mess_matrix structure variables using \ref mess_matrix_init
 * @snippet "tutorials/lrcf_adi/tutorial_lyap.c" INIT
 * \n
 *
 * Now we can use \ref mess_matrix_read_formated to read our matrices from a Matrix Market File Format and convert it to
 * the storage type \ref MESS_CSR and \ref MESS_DENSE. \n
 * \ref mess_matrix_read_formated dynamical allocates memory at this point.
 * @snippet "tutorials/lrcf_adi/tutorial_lyap.c" READ
 * \n
 *
 * ## 4. Solve the (Generalized) Lyapunov Equation
 * As already announced \ref mess_lyap will solve the (generalized) Lyapunov Equation.
 * We can indicate the case of the normal Lyapunov Equation by providing  @c NULL as  the second argument to \ref mess_lyap.
 * @snippet "tutorials/lrcf_adi/tutorial_lyap.c" SOLVE
 * \n
 *
 * ## 5. Compute the Relative Residual
 * We declare two @c double variables for storing the right-hand side norm and the absolute residual.
 * \ref mess_matrix_dynorm2 computes the norm of the right-hand side and \ref mess_lrcfadi_res2 or \ref mess_lrcfadi_res2g
 * are helpfull to compute the desired 2-norms.
 * @snippet "tutorials/lrcf_adi/tutorial_lyap.c" RES
 * \n
 *
 * ## 6. Don't Forget to Clear Memory
 * A rule of thumb is that every structure which needs an @c init call before usage must be cleared at the
 * at the end of your code. In our cases we have used \ref mess_matrix_init which means we have to use
 * \ref mess_matrix_clear to clear dynamical allocated memory.
 * @snippet "tutorials/lrcf_adi/tutorial_lyap.c" CLEAR
 * \n
 *
 *
 *
 * @sa @ref tutorials
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
    mess_matrix F,G,E,Z;
    /// [DECLARE]


    /// [INPUT]
    if(argc!=3 && argc!=4){
        printf("usage: %s F.mtx G.mtx [E.mtx]\n",argv[0]);
        return 1;
    }
    /// [INPUT]


    /// [INIT]
    mess_matrix_init(&F);
    mess_matrix_init(&G);
    if(argc==4){ mess_matrix_init(&E);}
    mess_matrix_init(&Z);
    /// [INIT]


    /// [READ]
    mess_matrix_read_formated(argv[1],F,MESS_CSR);
    mess_matrix_read_formated(argv[2],G,MESS_DENSE);
    if(argc==4){ mess_matrix_read_formated(argv[3],E,MESS_CSR);}
    /// [READ]


    /// [SOLVE]
    if(argc==3){
        mess_lyap(F,NULL,G,Z);
    }else{
        mess_lyap(F,E,G,Z);
    }
    /// [SOLVE]


    /// [RES]
    double nrmRHS=1, nrm=1;
    mess_matrix_dynorm2(G,&nrmRHS);
    if(argc==3){
        mess_lrcfadi_res2(F,G,Z,&nrm);
    }else{
        mess_lrcfadi_res2g(F,E,G,Z,&nrm);
    }
    printf("relative Residual= %10.15e\n",nrm/nrmRHS);
    /// [RES]


    /// [CLEAR]
    mess_matrix_clear(&F);
    mess_matrix_clear(&G);
    mess_matrix_clear(&Z);
    if(argc==4)
        mess_matrix_clear(&E);
    /// [CLEAR]

    return nrm/nrmRHS>1e-9;
}
///@endcond












