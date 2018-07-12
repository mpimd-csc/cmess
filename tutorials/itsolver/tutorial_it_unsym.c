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
 *
 * @file tutorials/itsolver/tutorial_it_unsym.c
 * @brief Demonstrate how to use various iterative solvers for unsymmetric systems of linear equations.
 *
 * @sa mess_solver_bicgstab
 * @sa mess_solver_gmres
 * @sa mess_solver_gmres_restart
 * @sa mess_solver_sor
 * @sa mess_precond_ilu0
 * @sa mess_precond_iluk
 * @sa mess_precond_diag
 *
 * # Tutorial: Iterative Solvers
 *
 * This function demonstrates how to use different iterative solvers to solve
 * \f[Ax = b \f]
 * for an unsymmetric matrix \f$ A \f$. \n
 * Demonstrated iterative solvers are
 * <ul>
 * <li> BICGSTAB,
 * <li> generalized minimal residual (GMRES),
 * <li> restarted GMRES,
 * <li> successive overrelaxation (SOR).
 * </ul>
 * Restarted GMRES and BICGSTAB are demonstrated using different preconditioners:
 * <ul>
 * <li> diagonal preconditioner,
 * <li> \f$ ILU(0) \f$ preconditioner (only used for BICGSTAB),
 * <li> \f$ IlU(k) \f$ preconditioner with a level of fill \f$ = 3 \f$.
 * </ul>
 *
 * @snippet "tutorials/itsolver/tutorial_it_unsym.c" CODE
 *
 *
 */

///@cond
///[CODE]

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/mess.h"


int main ( int argc, char **argv )
{
    mess_matrix matrix = NULL;
    mess_matrix input = NULL;
    mess_vector x = NULL;
    mess_vector b = NULL;
    mess_solver_options opt;
    mess_solver_status  stat;
    mess_mvpcall mvpcall;
    mess_precond pre_diag = NULL;
    mess_precond pre_ilu0 = NULL;
    mess_precond pre_sor  = NULL;
    mess_precond pre_iluk3 = NULL;
    int use_ilu0 = 0;
    int use_diag = 1;
    int use_ilu3 = 1;
    // int ret = 0;

    printf("mess - demo - iterative solvers for general systems\n");
    printf("====================================================\n");
    mess_version();


    mess_set_errorlevel(3);

    if ( argc != 5 ) {
        printf (" Usage: %s file.mtx maxit restarts tol\n", argv[0]) ;
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  read matrix and prepare
     *-----------------------------------------------------------------------------*/
    mess_matrix_init(&input);
    mess_matrix_init(&matrix);
    if (mess_matrix_read(argv[1], input)) {
        printf("error reading matrix\n");
        return 1;
    }
    mess_matrix_convert(input, matrix, MESS_CSR);
    mess_matrix_clear(&input);

    if (matrix->cols != matrix->rows){
        printf("Only with square matrices possible.\n");
    }

    MESS_INIT_VECTORS(&x,&b);
    mess_vector_alloc(x, matrix->rows, MESS_REAL);
    mess_vector_alloc(b, matrix->rows, MESS_REAL);
    mess_matrix_rowsums(matrix,b);


    /*-----------------------------------------------------------------------------
     *  Setup the preconditioners
     *-----------------------------------------------------------------------------*/
    mess_precond_init(&pre_diag);
    mess_precond_init(&pre_ilu0);
    mess_precond_init(&pre_sor);
    mess_precond_init(&pre_iluk3);

    if ( mess_precond_diag(pre_diag, matrix) !=0 ){
        printf("-> diagonal preconditioner can not be created. Not used.\n");
        use_diag = 0;
    }
    if ( mess_precond_ilu0(pre_ilu0, matrix) !=0) {
        printf("-> ilu 0 preconditioner can not be created. Not used.\n");
        use_ilu0 = 0;
    }
    if ( mess_precond_iluk(pre_iluk3, matrix, 3) !=0) {
        printf("-> ilul preconditioner can not be created. Not used.\n");
        use_ilu3 = 0;
    }

    /*-----------------------------------------------------------------------------
     *  fill options
     *-----------------------------------------------------------------------------*/
    mess_solver_options_init(&opt);
    opt->tol = atof(argv[4]);
    opt->maxit = atoi(argv[2]);
    opt->restarts = atoi(argv[3]);

    printf("\noptions:\n");
    mess_solver_options_print(opt);

    mess_solver_status_init(&stat);

    mess_mvpcall_matrix(&mvpcall, MESS_OP_NONE, matrix);

    /*-----------------------------------------------------------------------------
     *  GMRES
     *-----------------------------------------------------------------------------*/
    printf("\n\n-> GMRES\n");
    mess_vector_zeros(x);
    mess_solver_gmres(mvpcall , NULL, b, x, opt, stat);
    printf("Status:\n");
    mess_solver_status_print(stat);
    printf("solution:\n");
    mess_vector_printshort(x);
    mess_solver_status_clean(stat);

    /*-----------------------------------------------------------------------------
     *  restarted GMRES
     *-----------------------------------------------------------------------------*/
    printf("\n\n-> restarted GMRES\n");
    mess_vector_zeros(x);
    mess_solver_gmres_restart(mvpcall, NULL, b, x, opt, stat);
    printf("Status:\n");
    mess_solver_status_print(stat);
    printf("solution:\n");
    mess_vector_printshort(x);
    if(stat->relres>opt->tol){printf("FAILED\n"); return 1;}
    mess_solver_status_clean(stat);

    /*-----------------------------------------------------------------------------
     *  restarted GMRES with Diag
     *-----------------------------------------------------------------------------*/
    if ( use_diag) {
        printf("\n\n-> restarted GMRES with diagonal preconditioner\n");
        mess_vector_zeros(x);
        mess_solver_gmres_restart(mvpcall , pre_diag, b, x, opt, stat);
        printf("Status:\n");
        mess_solver_status_print(stat);
        printf("solution:\n");
        mess_vector_printshort(x);
        if(stat->relres>opt->tol){printf("FAILED\n"); return 1;}
        mess_solver_status_clean(stat);
    }

    /*-----------------------------------------------------------------------------
     *  restarted GMRES with ilu(3)
     *-----------------------------------------------------------------------------*/
    if ( use_ilu3) {
        printf("\n\n-> restarted GMRES with ilu(3) preconditioner\n");
        mess_vector_zeros(x);
        mess_solver_gmres_restart(mvpcall, pre_iluk3, b, x, opt, stat);
        printf("Status:\n");
        mess_solver_status_print(stat);
        printf("solution:\n");
        mess_vector_printshort(x);
        if(stat->relres>opt->tol){printf("FAILED\n"); return 1;}
        mess_solver_status_clean(stat);
    }

    /*-----------------------------------------------------------------------------
     * BICGSTAB
     *-----------------------------------------------------------------------------*/
    opt->maxit = opt->maxit * opt->restarts;
    printf("\n\n-> BiCGStab\n");
    mess_vector_zeros(x);
    mess_solver_bicgstab(mvpcall, NULL, b, x, opt, stat);
    printf("Status:\n");
    mess_solver_status_print(stat);
    printf("solution:\n");
    mess_vector_printshort(x);
    if(stat->relres>opt->tol){printf("FAILED\n"); return 1;}
    mess_solver_status_clean(stat);

    /*-----------------------------------------------------------------------------
     * BiCGStab with Diag
     *-----------------------------------------------------------------------------*/
    if ( use_diag) {
        printf("\n\n-> BiCGStab with diagonal preconditioner\n");
        mess_vector_zeros(x);
        mess_solver_bicgstab(mvpcall, pre_diag, b, x, opt, stat);
        printf("Status:\n");
        mess_solver_status_print(stat);
        printf("solution:\n");
        mess_vector_printshort(x);
        if(stat->relres>opt->tol){printf("FAILED\n"); return 1;}
        mess_solver_status_clean(stat);
    }

    /*-----------------------------------------------------------------------------
     *  BiCGStab with ilu(0)
     *-----------------------------------------------------------------------------*/
    if ( use_ilu0) {
        printf("\n\n-> BiCGStab with ilu(0) preconditioner\n");
        mess_vector_zeros(x);
        mess_solver_bicgstab(mvpcall, pre_ilu0, b, x, opt, stat);
        printf("Status:\n");
        mess_solver_status_print(stat);
        printf("solution:\n");
        mess_vector_printshort(x);
        if(stat->relres>opt->tol){printf("FAILED\n"); return 1;}
        mess_solver_status_clean(stat);
    }
    /*-----------------------------------------------------------------------------
     *  BiCGStab with iluk(3)
     *-----------------------------------------------------------------------------*/
    if ( use_ilu3) {
        printf("\n\n-> BiCGStab with ilu(3) preconditioner\n");
        mess_vector_zeros(x);
        mess_solver_bicgstab(mvpcall, pre_iluk3, b, x, opt, stat);
        printf("Status:\n");
        mess_solver_status_print(stat);
        printf("solution:\n");
        mess_vector_printshort(x);
        if(stat->relres>opt->tol){printf("FAILED\n"); return 1;}
        mess_solver_status_clean(stat);
    }

    /*-----------------------------------------------------------------------------
     *  SOR omega = 1
     *-----------------------------------------------------------------------------*/
    printf("\n\n-> SOR omega = 1\n");
    mess_vector_zeros(x);
    opt->omega = 1;
    mess_solver_sor(matrix, pre_ilu0, b, x, opt, stat);
    printf("Status:\n");
    mess_solver_status_print(stat);
    printf("solution:\n");
    mess_vector_printshort(x);
    if(stat->relres>opt->tol){printf("FAILED\n"); return 1;}
    mess_solver_status_clean(stat);

    /*-----------------------------------------------------------------------------
     *  SOR omega = 1.5
     *-----------------------------------------------------------------------------*/
    printf("\n\n-> SOR omega = 1.1\n");
    mess_vector_zeros(x);
    opt->omega = 1.1;
    mess_solver_sor(matrix, pre_ilu0, b, x, opt, stat);
    printf("Status:\n");
    mess_solver_status_print(stat);
    printf("solution:\n");
    mess_vector_printshort(x);
    if(stat->relres>opt->tol){printf("FAILED\n"); return 1;}
    mess_solver_status_clean(stat);

    /*-----------------------------------------------------------------------------
     *  SOR omega = 0.5
     *-----------------------------------------------------------------------------*/
    printf("\n\n-> SOR omega = 0.5\n");
    mess_vector_zeros(x);
    opt->omega = 0.5;
    mess_solver_sor(matrix, pre_ilu0, b, x, opt, stat);
    printf("Status:\n");
    mess_solver_status_print(stat);
    printf("solution:\n");
    mess_vector_printshort(x);
    if(stat->relres>opt->tol){printf("FAILED\n"); return 1;}
    mess_solver_status_clean(stat);



    /*-----------------------------------------------------------------------------
     *  cleanup
     *-----------------------------------------------------------------------------*/

    mess_mvpcall_clear(&mvpcall);
    mess_solver_options_clear(&opt);
    mess_solver_status_clear(&stat);
    mess_matrix_clear(&matrix);
    mess_precond_clear(&pre_ilu0);
    mess_precond_clear(&pre_diag);
    mess_precond_clear(&pre_sor);
    mess_precond_clear(&pre_iluk3);
    mess_vector_clear(&b);
    mess_vector_clear(&x);
    return 0;
}

///[CODE]
///@endcond
