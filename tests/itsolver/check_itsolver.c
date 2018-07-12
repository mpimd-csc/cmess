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
 * @addtogroup test_itsolver
 * @{
 * @file tests/itsolver/check_itsolver.c
 * @brief  Check different iterative solvers.
 * @author @koehlerm
 * @test
 * This function checks the @ref mess_solver_bicgstab function defined in bicgstab.c, @ref mess_solver_cg function defined in cg.c,
 * @ref mess_solver_gmres and @ref mess_solver_gmres_restart functions defined in gmres.c, @ref mess_solver_jacobi fucntion defined in
 * jacobi.c and @ref mess_solver_gaussseidel, @ref mess_solver_sor and @ref mess_solver_ssor functions defined in sor.c. \n
 * That means it checks if a linear system of equations
 * \f[ Ax=b \f]
 * is solved correctly using a certain solver.\n
 * Possible solver names are:
 * * "BICGSTAB"
 * * "CG"
 * * "GMRES"
 * * "GS"
 * * "JACOBI"
 * * "RGMRES"
 * * "SOR"
 * * "SSOR"
 *
 *
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>

#include "mess/mess.h"
#include "mess/error_macro.h"

#include "../call_macro.h"



/**
 * @brief Check different iterative solvers.
 * @param[in] name   input name of solver
 * @param[in] A      input matrix
 * @param[in] maxit  input maximum number of iteration
 * @param[in] tol    input tolerance for the iteration
 * @return zero on success or a non zero error code.
 *
 * The @ref check_solver function selects a solver by the \f$ name \f$, solves the linear system of equations
 * \f[ Ax=b \f]
 * and checks if the solution is computed correctly. \n
 * The entries of the right hand side \f$ b \f$ are computed by the sum of each row of a \f$ A \f$. \n
 * Possible solver names are:
 * <ul>
 * <li> "BICGSTAB",
 * <li> "CG",
 * <li> "GMRES",
 * <li> "GS",
 * <li> "JACOBI",
 * <li> "RGMRES",
 * <li> "SOR",
 * <li> "SSOR".
 * </ul>
 *
 */
int check_solver(const char *name, mess_matrix A, mess_int_t maxit, double tol) {
    mess_init();
    mess_solver_options opt;
    mess_solver_status stat;
    mess_vector x,b;
    double normb, res, resid;
    mess_mvpcall mvpcall;
    int ok = 0 ;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  Setup right hand side
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&x,&b);
    CALL(mess_vector_alloc(x,A->rows,A->data_type));
    CALL(mess_vector_alloc(b,A->rows,A->data_type));

    CALL(mess_matrix_rowsums(A,b));
    CALL(mess_vector_zeros(x));
    CALL(mess_vector_norm2(b,&normb));

    /*-----------------------------------------------------------------------------
     *  Setup Solver
     *-----------------------------------------------------------------------------*/
    CALL(mess_solver_options_init(&opt));
    CALL(mess_solver_status_init(&stat));
    opt->maxit = maxit;
    opt->tol = tol ;

    CALL(mess_mvpcall_matrix(&mvpcall, MESS_OP_NONE, A));
    /*-----------------------------------------------------------------------------
     *  Select the solver
     *-----------------------------------------------------------------------------*/
    if (strcasecmp(name,"CG") == 0) {
        CALL(mess_solver_cg(mvpcall, NULL, b, x, opt, stat));
    } else if ( strcasecmp(name, "BICGSTAB") == 0 ) {
        CALL(mess_solver_bicgstab(mvpcall, NULL, b, x, opt, stat));
    } else if ( strcasecmp(name, "GMRES") == 0 ) {
        CALL(mess_solver_gmres(mvpcall, NULL, b, x, opt, stat));
    } else if ( strcasecmp(name, "RGMRES") == 0 ) {
        CALL(mess_solver_gmres_restart(mvpcall, NULL, b, x, opt, stat));
    } else if ( strcasecmp(name, "JACOBI") == 0 ) {
        CALL(mess_solver_jacobi(A, NULL, b, x, opt, stat));
    } else if ( strcasecmp(name, "GS") == 0 ) {
        CALL(mess_solver_gaussseidel(A, NULL, b, x, opt, stat));
    } else if ( strcasecmp(name,"SOR") == 0 ) {
        opt->omega = 1.1;
        CALL(mess_solver_sor(A, NULL, b, x, opt, stat));
    } else if ( strcasecmp(name,"SSOR") == 0 ) {
        opt->omega = 1.0;
        CALL(mess_solver_ssor(A, NULL, b, x, opt, stat));
    }

    mess_mvpcall_clear(&mvpcall);

    printf("Status\n");
    CALL(mess_solver_status_print(stat));

    CALL(mess_direct_res2(MESS_OP_NONE, A, x,b,&res, &resid));

    if ( resid < tol ) ok = 0;
    else ok = 1;
    printf("Achieved residual: %18.12e  - %d \n", resid, ok?0:1);

    mess_solver_options_clear(&opt);
    mess_solver_status_clear(&stat);
    mess_vector_clear(&x);
    mess_vector_clear(&b);
    return ok;
}

int main(int argc, const char **argv)
{
    mess_matrix A;
    mess_int_t maxit;
    double tol;
    int ret = 0;

    if ( argc != 5 ) {
        fprintf (stderr, "ITSolver Test\n");
        fprintf (stderr, "Usage: %s SOLVER file.mtx maxit tol\n", argv[0]) ;
        exit(1);
    }

    /*-----------------------------------------------------------------------------
     *  Initialize the Matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_read_formated(argv[2], A, MESS_CSR));
    maxit = atoi(argv[3]);
    tol = atof(argv[4]);
    CALL(check_solver(argv[1],A,maxit,tol));

    mess_matrix_clear(&A);
    return 0;
}

