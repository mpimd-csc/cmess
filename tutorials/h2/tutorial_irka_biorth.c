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
 * @file tutorials/h2/tutorial_irka_biorth.c
 * @brief Demonstrate the IRKA algorithm.
 *
 * @sa mess_h2_irka_biorth
 *
 * # Tutorial: Irka Biorth algorithm
 *
 * This function demonstrates computing a \f$ \mathcal{H}_2 \f$ reduced (generalized) linear time invariant ((G)LTI)
 * system using the IRKA algorithm.
 *
 * @snippet "tutorials/h2/tutorial_irka_biorth.c" CODE
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
#include "mess/mess.h"
#include <complex.h>


void irka_stepdebug(void *data, mess_int_t it, double ds, double h2err){
    mess_plot p = (mess_plot) data;
    mess_plot_addData(p, 0, it, ds);
    if ( !isinf(h2err) ) mess_plot_addData(p, 1, it, h2err);
    mess_plot_update(p);
}

//#define PLOT
#undef PLOT

int main ( int argc, char **argv){
    mess_int_t r0, kp,km, ret = 0,useE=0;
    mess_vector sigma;
    mess_matrix A,B,C,E=NULL,V,W;
    mess_direct Asolver;
    mess_dynsys orig, reduced;
    mess_h2_options opt;
    mess_h2_status status;

#ifdef PLOT
    mess_plot plot;
#endif

    mess_version();

    printf("mess - IRKA BiOrth\n");
    printf("===================\n");


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if ( argc!=8 && argc!=9) {
        printf("usage: %s A.mtx B.mtx C.mtx (E.mtx) r0 maxit tol h2\n", argv[0]);      return 1;
    }
    if ( argc==9){
        printf(" -> for generalized LTI systems.\n");
        mess_matrix_init(&E);
    }
    mess_init();
    MESS_INIT_MATRICES(&A,&B,&C,&V,&W);

    /*-----------------------------------------------------------------------------
     *  init dynamical system strucutre and H2 options
     *-----------------------------------------------------------------------------*/
    mess_dynsys_init(&orig);
    mess_h2_options_init(&opt);
    mess_h2_status_init(&status);


    /*-----------------------------------------------------------------------------
     *  read matrices
     *-----------------------------------------------------------------------------*/
    mess_matrix_read_formated(argv[1], A, MESS_CSR);
    mess_matrix_read_formated(argv[2], B, MESS_DENSE);
    mess_matrix_read_formated(argv[3], C, MESS_DENSE);
    if (E){
        mess_matrix_read_formated(argv[4],E,MESS_CSR);
        useE=1;
    }
    if (B->cols!=1){
        //only siso systems
        mess_matrix tmp; MESS_INIT_MATRICES(&tmp); mess_matrix_colsub(B,0,0,tmp);
        mess_matrix_copy(tmp,B); mess_matrix_clear(&tmp);
    }

    if (C->rows!=1){
        mess_matrix tmp; MESS_INIT_MATRICES(&tmp); mess_matrix_rowsub(C,0,0,tmp);
        mess_matrix_copy(tmp,C); mess_matrix_clear(&tmp);
    }

#ifdef PLOT
    /*-----------------------------------------------------------------------------
     *  plot
     *-----------------------------------------------------------------------------*/
    mess_plot_create(&plot,800,600,"IRKA BiOrth", MESS_PLOT_LIN, MESS_PLOT_LOG, 2);
    mess_plot_setLabel(plot,0, "rel |sigma-sigmaold|");
    mess_plot_setLabel(plot,1, "h2 error");

    mess_plot_setColor(plot,0, "red");
    mess_plot_setColor(plot,1, "blue");

    mess_plot_update(plot);

#endif

    /*-----------------------------------------------------------------------------
     *  H2 options
     *-----------------------------------------------------------------------------*/
#ifdef PLOT
    opt->stepdebug = irka_stepdebug;
    opt->stepdebug_data = (void*) plot;
#endif
    r0=atoi(argv[4+useE]);
    opt->maxit = atoi(argv[5+useE]);
    opt->tol=atof(argv[6+useE]);
    opt->calc_h2err = atoi(argv[7+useE]);
    opt->calc_finalh2  = 0;
    opt->output = 1;

    kp=3*r0;
    km=2*r0;


    printf("Dimension: " MESS_PRINTF_INT "\n", A->rows);
    printf("wanted number of parameters: " MESS_PRINTF_INT "\n", r0);
    printf("Arnoldi steps w.r.t. A: " MESS_PRINTF_INT "\n", kp);
    printf("Arnoldi steps w.r.t. inv(A): " MESS_PRINTF_INT "\n", km);
    printf("\nIRKA-options:\n");
    mess_h2_options_print(opt);
    printf("\n");


    /*-----------------------------------------------------------------------------
     *  compute irka parameters
     *-----------------------------------------------------------------------------*/
    mess_direct_init(&Asolver);
    mess_direct_lu(A,Asolver);

    mess_vector_init(&sigma);
    mess_vector_alloc(sigma, r0, MESS_COMPLEX);
    ret = mess_h2_irka_init(A,Asolver, E, &r0, kp,km, sigma);
    mess_direct_clear(&Asolver);

    if ( ret !=0) {
        printf("error computing IRKA parameters: %s\n", mess_get_error(ret));
    } else {
        printf("number of parameters: " MESS_PRINTF_INT "\n", r0);
        mess_vector_print(sigma);
    }

    if  (useE)  mess_dynsys_glti(orig, A,B,C,E);
    else mess_dynsys_lti(orig, A,B,C);

    /*-----------------------------------------------------------------------------
     *  reduce system
     *-----------------------------------------------------------------------------*/
    mess_dynsys_init(&reduced);
    mess_h2_irka_biorth(orig,sigma,opt,reduced,V,W, status);
    mess_dynsys_clear(&reduced);

    /*-----------------------------------------------------------------------------
     *  write status
     *-----------------------------------------------------------------------------*/
    printf("status: \n");
    mess_h2_status_print(status);
    mess_h2_status_write_mfile(status,"irka_output.m");


    /*-----------------------------------------------------------------------------
     *  evaluate transfer function
     *-----------------------------------------------------------------------------*/
    {
        mess_vector omega, G;
        MESS_INIT_VECTORS(&omega,&G);
        mess_vector_alloc(omega,200,MESS_REAL);
        mess_vector_alloc(G,200, MESS_REAL);
        mess_dynsys_evaltransfer(orig,-1,6,200,omega,G,NULL);

        mess_vector_printshort(omega);
        mess_vector_printshort(G);

        mess_vector_clear(&G);
        mess_vector_clear(&omega);
    }

#ifdef PLOT
    /*-----------------------------------------------------------------------------
     *  save plot
     *-----------------------------------------------------------------------------*/
    sleep(5);
    mess_plot_save(plot,"plot.xpm");
#endif

    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
#ifdef PLOT
    mess_plot_close(&plot);
#endif
    MESS_CLEAR_MATRICES(&V,&W);
    //MESS_CLEAR_MATRICES(&A,&B,&C,&E); cleared by dynsys_clear
    mess_dynsys_clear(&orig);
    mess_h2_options_clear(&opt);
    mess_h2_status_clear(&status);
    mess_vector_clear(&sigma);
    mess_exit();

    return 0;
}
///[CODE]
///@endcond
