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
 * @file tutorials/lrcf_adi/tutorial_lrnm-g.c
 * @brief Demonstrate how to solve a generalized algebraic Riccati Equation using LRCF-ADI-NM.
 *
 * # Tutorial: Solve generalized algebraic Riccati Equatio, Plot and Stepdebug Function
 *
 * This function demonstrates computing the solution factor \f$ Z \f$ of the solution \f$ X \approx ZZ^T \f$ of
 * a generalized algebraic Riccati Equation
 * \f[ A^T X E + E^T X A - C^T C + E^T X B B^T X E = 0\f]
 * using the low-rank Cholesky factor alternate direction implicit Newton method (LRCF-ADI-NM).
 *
 *
 * @snippet "tutorials/lrcf_adi/tutorial_lrnm-g.c" CODE
 *
 *
 */


///@cond
///[CODE]

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"

//#define PLOT
#undef PLOT

#ifdef PLOT
void lrcf_stepfunction(void *data,mess_equation eqn, mess_options opt, mess_lrcfadi_step * step){
    mess_plot p = (mess_plot) data;
    if ( step->it == 0) {
        mess_plot_clearData(p,0);
        mess_plot_clearData(p,1);
        mess_plot_clearData(p,2);
        mess_plot_update(p);
    }
    mess_plot_addData(p, 0, step->it, step->res2);
    mess_plot_addData(p, 1, step->it, step->rel_change);
    mess_plot_addData(p, 2, step->it, (step->res2_change));
    mess_plot_update(p);
}
void lrnm_stepfunction(void *data, mess_equation eqn, mess_equation lyap, mess_options opt, mess_lrcfadi_step * step){
    mess_plot p = (mess_plot) data;
    if ( step->it == 0) {
        mess_plot_clearData(p,0);
        mess_plot_clearData(p,1);
        mess_plot_update(p);
    }
    mess_plot_addData(p, 0, step->it, step->res2);
    mess_plot_addData(p, 1, step->it, (step->res2_change));
    mess_plot_update(p);
}
#endif


int main ( int argc, char **argv){
    mess_matrix A,B,C,E,Z;
    mess_options opt;
    mess_equation eqn;
    mess_status stat;
#ifdef PLOT
    mess_plot plot;
    mess_plot plot_nm;
#endif

    mess_version();
    printf("mess - low-rank-cholesky-factor-nm-g  demo plotting\n");
    printf("====================================================\n");


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc!=5) {
        printf("usage: %s A.mtx B.mtx C.mtx E.mtx\n", argv[0]);
        return 1;
    }

    mess_init();

    /*-----------------------------------------------------------------------------
     *  init and read matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&B,&C,&E,&Z);

    mess_matrix_read_formated(argv[1],A, MESS_CSR);
    mess_matrix_read_formated(argv[2],B, MESS_DENSE);
    mess_matrix_read_formated(argv[3],C, MESS_DENSE);
    mess_matrix_read_formated(argv[4],E, MESS_CSR);

    //scale to force more lrnm steps
    mess_matrix_scale(1e+2,C);
    /*-----------------------------------------------------------------------------
     *  Plot for LRCFADI
     *-----------------------------------------------------------------------------*/
#ifdef PLOT
    mess_plot_create(&plot,800,600,"LRCF-ADI residual", MESS_PLOT_LIN, MESS_PLOT_LOG, 3);
    mess_plot_setLabel(plot,0, "res2-norm");
    mess_plot_setLabel(plot,1, "rel-change");
    mess_plot_setLabel(plot,2, "res2-change");
    mess_plot_setColor(plot,0, "red");
    mess_plot_setColor(plot,1, "blue");
    mess_plot_setColor(plot,2, "green");
    mess_plot_update(plot);

    /*-----------------------------------------------------------------------------
     *  PLOT for LRCF-NM
     *-----------------------------------------------------------------------------*/
    mess_plot_create(&plot_nm,800,600,"LRCF-NM residual", MESS_PLOT_LIN, MESS_PLOT_LOG, 3);
    mess_plot_setLabel(plot_nm,0, "res2-norm");
    mess_plot_setLabel(plot_nm,1, "rel-change");
    mess_plot_setLabel(plot_nm,2, "res2-change");
    mess_plot_setColor(plot_nm,0, "red");
    mess_plot_setColor(plot_nm,1, "blue");
    mess_plot_setColor(plot_nm,2, "green");
    mess_plot_update(plot_nm);
#endif

    /*-----------------------------------------------------------------------------
     *  init status and options
     *-----------------------------------------------------------------------------*/
    mess_status_init(&stat);
    mess_options_init(&opt);
    opt->nm_gpStep=0;
    mess_options_print(opt);

#ifdef PLOT
    opt->stepfunction           = lrcf_stepfunction;
    opt->stepfunction_aux       = plot;
    opt->nm_stepfunction        = lrnm_stepfunction;
    opt->nm_stepfunction_aux    = plot_nm;
#endif

    /*-----------------------------------------------------------------------------
     *  create equation
     *-----------------------------------------------------------------------------*/
    mess_equation_init(&eqn);
    mess_equation_riccati(eqn,opt,A,E,B,C);

    /*-----------------------------------------------------------------------------
     *  solve riccati equation
     *-----------------------------------------------------------------------------*/
    mess_lrcfadi_nm(eqn,opt,stat,Z);

    /*-----------------------------------------------------------------------------
     *  print status information
     *-----------------------------------------------------------------------------*/
    printf("\n state structure :\n");
    mess_status_print(stat);

    /*-----------------------------------------------------------------------------
     *  write matrix
     *-----------------------------------------------------------------------------*/
    printf("-> writing Z factor to %s\n","Z.mtx" );
    mess_matrix_write( "Z.mtx" , Z);

    /*-----------------------------------------------------------------------------
     *  plot
     *-----------------------------------------------------------------------------*/
#ifdef PLOT
    sleep(10);
    mess_plot_save(plot,"plot.xpm");
    mess_plot_close(&plot);
    mess_plot_close(&plot_nm);
#endif

    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&B,&C,&E,&Z);
    mess_equation_clear(&eqn);
    mess_status_clear(&stat);
    mess_options_clear(&opt);

    mess_exit();
    return 0;
}

///[CODE]
///@endcond
