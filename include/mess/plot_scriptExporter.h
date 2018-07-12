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
 * @file include/mess/plot_scriptExporter.h
 * @author @weiss
 * @brief Interface definition of the Gnuplot module.
 */

/** @addtogroup plot_export
  @{
 * This module provides a small and easy to use values plotter based on the X11 library. It
 * can be used to create linear and logarithmic plots. \n
 * An example to plot sine and cosine will be:
 *
 * \code{.c}
 * @ref mess_plotScript plot;
 * @ref mess_int_t i, N=1000;
 * double x,y1,y2;
 * double h = 2.0*3.14159/N;
 * @ref mess_gnuplot_create(&plot, 400, 300, "Test-Plot", MESS_PLOT_LIN, MESS_PLOT_LIN, 2);
 * @ref mess_gnuplot_setlabel(plot, 0, "sin(x)");
 * @ref mess_gnuplot_setcolor(plot, 0, "Red");
 * @ref mess_gnuplot_setlabel(plot, 1, "cosine(x)");
 * @ref mess_gnuplot_setcolor(plot, 1, "Blue");
 * for ( i = 0; i <= N; i++){
 *  x = h*i;
 *  y1 = sin(x);
 *  y2 = cos(x);
 *  @ref mess_gnuplot_adddata(plot,0, x,y1 );
 *  @ref mess_gnuplot_adddata(plot,1, x,y2 );
 * }
 * @ref mess_gnuplot_createScript();
 * \endcode
 */

#ifdef  __cplusplus
extern "C" {
#endif

    /**
     * @brief Enumeration type for setting axes scaling.
     *
     */
    typedef enum {
        MESS_PLOT_LIN = 0,        /**< Linear axis scaling for plot.*/
        MESS_PLOT_LOG = 1         /**< Logarithmic axis scaling for plot.*/
    } mess_plot_axis_scale_t;


    /**
     * @internal
     * @brief Structure for one data element in a plot.
     *
     * The @ref mess_plotExportData structure holds all information of a data element in a plot.
     *
     * @attention Internal use only.
     * */
    typedef struct __mess_plotExportData {
        double *x;              /**< \f$ x \f$ values of the data element */
        double *y;              /**< \f$ y \f$ values of the data element */
        char color[40];         /**< Colorname to be used for the plot*/
        char label[40];         /**< Label of the plot */
        char type[40];          /**< Type of the plot */
        mess_int_t len;         /**< Number of values */
    } mess_plotExportData;

    /**
     * @internal
     * @brief Structure to hold a complete plotter.
     *
     * The @ref __mess_plotterScript is the basic structure to represent a complete plot.
     *
     * @attention Internal use only.
     * */
    struct __mess_plotterScript {
        char *title;                    /**< Title of the window */
        char *legendGnuPos;             /**< Orientation of key-box for Gnu */
        char *legendTikzPos;            /**< Orientation of key-box for Tikz */
        mess_plotExportData *plot;      /**< Data elements  */
        int plot_len;                   /**< Number of data elements */
        mess_plot_axis_scale_t xscale;  /**< Scaling of the \f$ x \f$ axis */
        mess_plot_axis_scale_t yscale;  /**< Scaling of the \f$ y \f$ axis */
        char *xLabel;                   /**< Label of the \f$ x \f$ axis */
        char *yLabel;                   /**< Label of the \f$ y \f$ axis */

    };

    /**
     * @internal
     * @brief Type definition for the plotter structure.
     *
     * The @ref mess_plotExport type definition is for the plotter structure to use in programms.
     *
     * @attention Internal use only.
     * */
    typedef struct __mess_plotterScript *mess_plotExport;

    int mess_plotExport_init(mess_plotExport *p, char *title, mess_plot_axis_scale_t xscale, mess_plot_axis_scale_t yscale, char *xLabel, char *yLabel, int data);
    int mess_plotExport_addData(mess_plotExport p, int data, double x, double y);
    int mess_plotExport_clearData (mess_plotExport p, int data);
    int mess_plotExport_setLegendPos(mess_plotExport p, char *orientation);
    int mess_plotExport_setLabel(mess_plotExport p, int data, char *label);
    int mess_plotExport_setColor(mess_plotExport p, int data, char *color);
    int mess_plotExport_setType(mess_plotExport p, int data, char *type);
    int mess_plotExport_save(mess_plotExport p, const char *path) ;
    int mess_plotExport_createGnuScript(mess_plotExport p);
    int mess_plotExport_createGnuScript_filename(mess_plotExport p, const char *path);
    int mess_plotExport_createTikzScript(mess_plotExport p);
    int mess_plotExport_createTikzScript_filename(mess_plotExport p, const char *path);
    int mess_plotExport_clear(mess_plotExport *plot);
#ifdef  __cplusplus
}
#endif
/** @}  */
