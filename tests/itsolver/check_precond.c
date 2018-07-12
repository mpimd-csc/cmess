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
 * @file tests/itsolver/check_precond.c
 * @brief  Check different preconditioner.
 * @author @koehlerm
 * @test
 * This function checks the @ref mess_precond_diag function defined in diag.c,mess_precond_ilu0 function defined in ilu0.c,
 * @ref mess_precond_iluk  function defined in iluk.c, mess_precond_sor and mess_precond_ssor functions defined in sor.c. \n
 * That means it checks if a linear system of equations
 * \f[ Ax=b \f]
 * is solved correctly using a certain preconditioner.\n
 * Possible preconditioner names are:
 * <ul>
 * <li> "DIAG",
 * <li> "ILU0",
 * <li> "ILUK0",
 * <li> "ILUK1",
 * <li> "ILUK2",
 * <li> "SOR",
 * <li> "SSOR".
 * </ul>
 * that should be put in the second argument of the second input.
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
 * @brief Check different preconditioner.
 * @param[in] name   input name of preconditioner
 * @param[in] A      input matrix
 * @param[in] maxit  input maximum number of iteration
 * @param[in] tol    input tolerance for the iteration
 * @return zero on success or a non zero error code.
 *
 * The @ref check_precond function selects a preconditioner by the \f$ name \f$, solves the linear system of equations
 * \f[ Ax=b \f]
 * using the computed preconditioner and checks if the solution is computed correctly. \n
 * The entries of the right hand side \f$ b \f$ are computed by the sum of each row of a \f$ A \f$. \n
 * Possible preconditioner names are:
 * <ul>
 * <li> "DIAG",
 * <li> "ILU0",
 * <li> "ILUK0",
 * <li> "ILUK1",
 * <li> "ILUK2",
 * <li> "SOR",
 * <li> "SSOR".
 * </ul>
 *
 */
int check_precond(const char *name, mess_matrix A, mess_int_t maxit, double tol) {
    mess_init();
    mess_solver_options opt;
    mess_solver_status stat;
    mess_precond pre;
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
    CALL(mess_precond_init(&pre));

    /*-----------------------------------------------------------------------------
     *  Setup Solver
     *-----------------------------------------------------------------------------*/
    CALL(mess_solver_options_init(&opt));
    CALL(mess_solver_status_init(&stat));
    opt->maxit = maxit;
    opt->tol = tol ;

    CALL(mess_mvpcall_matrix(&mvpcall, MESS_OP_NONE, A));
    /*-----------------------------------------------------------------------------
     *  Select the preconditioner
     *-----------------------------------------------------------------------------*/
    if (strcasecmp(name,"DIAG") == 0 ) {
        CALL(mess_precond_diag(pre, A));
    } else if (strcasecmp(name,"ILU0") == 0 ) {
        CALL(mess_precond_ilu0(pre, A));
    } else if (strcasecmp(name,"ILUK0") == 0 ) {
        CALL(mess_precond_iluk(pre, A,0));
    } else if (strcasecmp(name,"ILUK1") == 0 ) {
        CALL(mess_precond_iluk(pre, A,1));
    } else if (strcasecmp(name,"ILUK2") == 0 ) {
        CALL(mess_precond_iluk(pre, A,2));
    } else if (strcasecmp(name,"SOR") == 0 ) {
        CALL(mess_precond_sor(pre, A,20,1.1));
    } else if (strcasecmp(name,"SSOR") == 0 ) {
        CALL(mess_precond_ssor(pre, A,20,1.1));
    }
    CALL(mess_solver_bicgstab(mvpcall, pre, b, x, opt, stat));

    mess_mvpcall_clear(&mvpcall);
    mess_precond_clear(&pre);

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
    CALL(check_precond(argv[1],A,maxit,tol));

    mess_matrix_clear(&A);
    return 0;
}

