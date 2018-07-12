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
 * @file include/mess/plot_x11.h
 * @author @koehlerm
 * @brief Interface definition of the plotting module.
 */

/** @addtogroup plot_x11
  @{
 * This module provides a small and easy to use values plotter based on the X11 library. It
 * can be used to create linear and logarithmic plots. \n
 * An example to plot sine and cosine will be:
 *
 * \code{.c}
 * @ref mess_plot plot;
 * @ref mess_int_t i, N=1000;
 * double x,y1,y2;
 * double h = 2.0*3.14159/N;
 * @ref mess_plot_create(&plot, 400, 300, "Test-Plot", MESS_PLOT_LIN, MESS_PLOT_LIN, 2);
 * @ref mess_plot_setLabel(plot, 0, "sin(x)");
 * @ref mess_plot_setColor(plot, 0, "Red");
 * @ref mess_plot_setLabel(plot, 1, "cosine(x)");
 * @ref mess_plot_setColor(plot, 1, "Blue");
 * for ( i = 0; i <= N; i++){
 *  x = h*i;
 *  y1 = sin(x);
 *  y2 = cos(x);
 *  @ref mess_plot_addData(plot,0, x,y1 );
 *  @ref mess_plot_addData(plot,1, x,y2 );
 * }
 * @ref mess_plot_update(plot);
 * sleep(10);
 * @ref mess_plot_close(&plot);
 * \endcode
 */

#ifndef PLOT_X11_C_
#define PLOT_X11_C_

#ifdef  __cplusplus
extern "C" {
#endif

#ifdef MESS_HAVE_X11
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <pthread.h>
#endif

    /**
     * @brief Structure for one data element in a plot.
     *
     * The @ref mess_plotdata structure holds all information of a data element in a plot.
     * */
    typedef struct __mess_plotdata {
        double *x;  /**< \f$  x \f$  values of the data element */
        double *y;  /**< \f$  y \f$  values of the data element */
        char color[40];     /**< Name of the plotting color */
        char label[40];     /**< Label of the plot */
        mess_int_t len; /**< Number of values */
    } mess_plotdata;


    /**
     * @internal
     * @brief Structure to hold a complete plotter.
     *
     * The @ref __mess_plotter is the basic structure to represent a complete plot.
     * @attention Interal use only.
     * */
    struct __mess_plotter {
#ifdef MESS_HAVE_X11
        Display *display;       /**< X11 display */
        Window win;             /**< X11 window */
        GC gc;                  /**< X11 graphic context */
        int width;              /**< Width of the plot */
        int height;             /**< Height of the plot */
        int screen_num;         /**< X11 screen number */
        Atom delete_atom;       /**< X11 destroy window atom */
        int done;               /**< Done flag */
        char *title;            /**< Title of the window */
        pthread_t thread;       /**< Worker thread for the plot */
        pthread_mutex_t lock;   /**< Lock for the thread */
        pthread_cond_t  sig;    /**< Signal for the lock */
        int editing;            /**< State variable for the lock */
        mess_plotdata *plot;    /**< Data elements  */
        int plot_len;           /**< Number of data elements */
        int xscale;             /**< Scaling of the \f$ x \f$ axis */
        int yscale;             /**< Scaling of the \f$ y \f$ axis */
#else
        int dummy;
        /* dummy, to ensure cplusplus compatibility
         * empty structs have size 0 in C and size 1 in C++
         */
#endif
    };
    /**
     * @brief Type definition for the plotter structure.
     *
     * The @ref mess_plot type definition is for the plotter structure to use in programms. */
    typedef struct __mess_plotter *mess_plot;

    int mess_plot_create(mess_plot *p, mess_int_t width, mess_int_t height, char *title, int xscale, int yscale, int data);
    int mess_plot_addData(mess_plot p, int data, double x, double y);
    int mess_plot_clearData (mess_plot p, int data) ;
    int mess_plot_setLabel(mess_plot p, int data, char *label);
    int mess_plot_setColor(mess_plot p, int data, char *color);
    int mess_plot_save(mess_plot p, char *filename) ;
    int mess_plot_update(mess_plot p);
    int mess_plot_close(mess_plot *p);



#ifdef  __cplusplus
}
#endif


#endif /* PLOT_X11_C_ */
/* @} */
