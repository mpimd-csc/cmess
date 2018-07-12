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
 * @file tutorials/h2/tutorial_tsia.c
 * @brief Demonstrate how to reduce a LTI system using the TSIA algorithm.
 *
 * @sa mess_h2_tsia
 *
 * # Tutorial: TSIA Algorithm
 * This function demonstrates computing a \f$ \mathcal{H}_2 \f$ reduced linear time invariant (LTI)
 * system using the TSIA algorithm.
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
#include <complex.h>


//#define PLOT
#undef PLOT

void irka_stepdebug(void *data, mess_int_t it, double ds, double h2err){
    mess_plot p = (mess_plot) data;
    mess_plot_addData(p, 0, it, ds);
    if ( !isinf(h2err) ) mess_plot_addData(p, 1, it, h2err);
    mess_plot_update(p);
}

int main ( int argc, char **argv){
    mess_matrix A,B,C,V,W,E;
    mess_dynsys orig, reduced;
    mess_h2_options opt;
    mess_h2_status status;
    mess_int_t useE=0;
#ifdef PLOT
    mess_plot plot;
#endif
    mess_int_t r0, kp,km;
    mess_version();

    printf("mess - TSIA\n");
    printf("============\n");


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc!=8 && argc!=9 ) {
        printf("usage: %s A.mtx B.mtx C.mtx [E.mtx] r0 maxit tol h2\n", argv[0]);
        return 0;
    }
    mess_init();

    //mess_error_level = 2;

    /*-----------------------------------------------------------------------------
     *  init and read matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A,&B,&C,&V,&W);
    mess_dynsys_init(&orig);
    mess_h2_options_init(&opt);
    mess_h2_status_init(&status);

    // read data
    mess_matrix_read_formated(argv[1], A, MESS_CSR);
    mess_matrix_read_formated(argv[2], B, MESS_DENSE);
    mess_matrix_read_formated(argv[3], C, MESS_DENSE);
    if(argc==9){
        mess_matrix_init(&E);
        mess_matrix_read_formated(argv[4], E, MESS_CSR);
        useE = 1;
    }else{
        E=NULL;
    }


    /*-----------------------------------------------------------------------------
     *  create plot
     *-----------------------------------------------------------------------------*/

#ifdef PLOT
    mess_plot_create(&plot,800,600,"IRKA BiOrth", MESS_PLOT_LIN, MESS_PLOT_LOG, 2);
    mess_plot_setLabel(plot,0, "rel |sigma-sigmaold|");
    mess_plot_setLabel(plot,1, "h2 error");

    mess_plot_setColor(plot,0, "red");
    mess_plot_setColor(plot,1, "blue");

    mess_plot_update(plot);

    opt->stepdebug = irka_stepdebug;
    opt->stepdebug_data = (void*) plot;
#endif


    /*-----------------------------------------------------------------------------
     *  H2 options
     *-----------------------------------------------------------------------------*/
    opt->calc_h2err     = atoi(argv[7+useE]);
    opt->calc_finalh2   = 1;
    opt->output         = 1;
    r0                  = atoi(argv[4+useE]);
    opt->maxit          = atoi(argv[5+useE]);
    opt->tol            = atof(argv[6+useE]);
    opt->rdim           = r0;
    kp = 3*r0;
    km = 2*r0;

    printf("Dimension: " MESS_PRINTF_INT "\n", A->rows);
    printf("wanted number of parameters: " MESS_PRINTF_INT "\n", r0);
    printf("Arnoldi steps w.r.t. A: " MESS_PRINTF_INT "\n", kp);
    printf("Arnoldi steps w.r.t. inv(A): " MESS_PRINTF_INT "\n", km);
    printf("\nTSIA-options:\n");
    mess_h2_options_print(opt);
    printf("\n");


    /*-----------------------------------------------------------------------------
     *  create dynsys system and start h2 tsia
     *-----------------------------------------------------------------------------*/
    if(useE){
        mess_dynsys_glti(orig,A,B,C,E);
    }else{
        mess_dynsys_lti(orig,A,B,C);
    }

    mess_dynsys_init(&reduced);
    if(useE){
        mess_h2_tsiag(orig,NULL,opt,reduced,V,W, status);
    }else{
        mess_h2_tsia(orig,NULL,opt,reduced,V,W, status);
    }
    mess_dynsys_clear(&reduced);
    printf("status: \n");
    mess_h2_status_print(status);
    mess_h2_status_write_mfile(status,"tsia_output.m");


    /*-----------------------------------------------------------------------------
     *  save plot
     *-----------------------------------------------------------------------------*/
#ifdef PLOT
    sleep(10);
    mess_plot_save(plot,"plot.xpm");
    mess_plot_close(&plot);
#endif


    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&V,&W);
    mess_dynsys_clear(&orig); //A, B, C, E cleared here
    mess_dynsys_clear(&reduced);
    mess_h2_options_clear(&opt);
    mess_h2_status_clear(&status);
    mess_exit();
    return 0;
}

///[CODE]
///@endcond
